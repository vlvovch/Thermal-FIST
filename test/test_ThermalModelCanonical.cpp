/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include <utility>
#include <algorithm>
#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelParameters.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGBase/ThermalModelCanonical.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  // ---- Helpers ---------------------------------------------------------------

  /// Build a default parameter set (Boltzmann, zero-width, T=155 MeV)
  ThermalModelParameters MakeParams(double V, int B = 0, int Q = 0,
                                    int S = 0, int C = 0, double muB = 0.)
  {
    ThermalModelParameters p;
    p.T = 0.155;
    p.muB = muB;  p.muS = 0.;  p.muQ = 0.;  p.muC = 0.;
    p.gammaS = 1.;  p.gammaq = 1.;
    p.V = V;  p.SVc = V;
    p.B = B;  p.Q = Q;  p.S = S;  p.C = C;
    return p;
  }

  /// Create, configure, and run a canonical model
  ThermalModelCanonical* RunCanonical(
      ThermalParticleSystem* TPS,
      const ThermalModelParameters& params,
      int method = GaussLegendre,
      int multiplier = 1,
      bool ceB = true, bool ceQ = true, bool ceS = true, bool ceC = true)
  {
    auto* m = new ThermalModelCanonical(TPS, params);
    m->SetStatistics(false);
    m->SetUseWidth(ThermalParticle::ZeroWidth);
    m->SetMethod(static_cast<CanonicalMethod>(method));
    if (method == GaussLegendre)
      m->SetIntegrationIterationsMultiplier(multiplier);
    m->ConserveBaryonCharge(ceB);
    m->ConserveElectricCharge(ceQ);
    m->ConserveStrangeness(ceS);
    m->ConserveCharm(ceC);
    m->CalculateDensities();
    return m;
  }

  /// Count species whose densities differ by more than relTol.
  /// Also counts NaN/Inf.  Returns {mismatches, nanCount}.
  std::pair<int,int> CompareDensities(
      const ThermalModelCanonical& a,
      const std::vector<double>& ref,
      const ThermalParticleSystem* TPS,
      double relTol, int maxPrint = 3,
      const std::string& tag = "")
  {
    int mis = 0, nan = 0;
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      double d  = a.Densities()[i];
      double dr = ref[i];
      if (std::isnan(d) || std::isinf(d)) { nan++; continue; }
      if (std::abs(dr) > 1.e-15 && std::abs(d) > 1.e-15) {
        double r = std::abs(d - dr) / std::max(std::abs(d), std::abs(dr));
        if (r > relTol) {
          mis++;
          if (mis <= maxPrint)
            std::cout << "  " << tag << TPS->Particle(i).Name()
                      << ": ref=" << dr << " got=" << d << " rel=" << r << "\n";
        }
      }
    }
    return {mis, nan};
  }

  // ---- Fixtures --------------------------------------------------------------

  class ThermalModelCanonicalTest : public ::testing::Test {
  protected:
    void SetUp() override {
      std::string dir = std::string(ThermalFIST_INPUT_FOLDER);
      TPS = new ThermalParticleSystem(dir + "/list/PDG2020/list.dat",
                                      dir + "/list/PDG2020/decays.dat");
    }
    void TearDown() override { delete TPS; }
    ThermalParticleSystem* TPS;
  };

  class ThermalModelCanonicalCharmTest : public ::testing::Test {
  protected:
    void SetUp() override {
      std::string dir = std::string(ThermalFIST_INPUT_FOLDER);
      std::vector<std::string> lists = {
        dir + "/list/PDG2020_modular/list-hadrons.dat",
        dir + "/list/PDG2020_modular/list-charm.dat"
      };
      std::vector<std::string> decays;
      TPS = new ThermalParticleSystem(lists, decays);
    }
    void TearDown() override { delete TPS; }
    ThermalParticleSystem* TPS;
  };

  // ---- GL tests --------------------------------------------------------------

  // Test 1: GL charge conservation at N=0 and N!=0
  TEST_F(ThermalModelCanonicalTest, GLChargeConservation) {
    // (a) B=Q=S=0, V=5000
    {
      auto p = MakeParams(5000.);
      std::unique_ptr<ThermalModelCanonical> m(RunCanonical(TPS, p));
      EXPECT_NEAR(m->CalculateBaryonDensity()      * p.SVc, 0., 1.e-4);
      EXPECT_NEAR(m->CalculateChargeDensity()       * p.SVc, 0., 1.e-4);
      EXPECT_NEAR(m->CalculateStrangenessDensity()  * p.SVc, 0., 1.e-4);
    }
    // (b) B=2, Q=1, S=0, V=100
    {
      auto p = MakeParams(100., 2, 1, 0);
      std::unique_ptr<ThermalModelCanonical> m(RunCanonical(TPS, p, GaussLegendre, 2));
      EXPECT_NEAR(m->CalculateBaryonDensity()      * p.SVc, 2., 1.e-3);
      EXPECT_NEAR(m->CalculateChargeDensity()       * p.SVc, 1., 1.e-3);
      EXPECT_NEAR(m->CalculateStrangenessDensity()  * p.SVc, 0., 1.e-3);
    }
  }

  // Test 2: GL GCE limit and convergence
  TEST_F(ThermalModelCanonicalTest, GLGCELimitAndConvergence) {
    // (a) GCE limit: CE at V=100000 vs GCE
    {
      auto p = MakeParams(100000.);
      std::unique_ptr<ThermalModelCanonical> mCE(RunCanonical(TPS, p));

      ThermalModelIdeal mGCE(TPS, p);
      mGCE.SetStatistics(false);
      mGCE.SetUseWidth(ThermalParticle::ZeroWidth);
      mGCE.CalculateDensities();

      auto [mis, nan] = CompareDensities(*mCE, mGCE.Densities(), TPS, 1.e-2, 5, "GCE ");
      EXPECT_EQ(mis, 0) << mis << " species exceed GCE tolerance";
      EXPECT_EQ(nan, 0) << nan << " species with NaN/Inf density (GCE limit)";
    }
    // (b) Convergence: multiplier=1 vs 2 at V=5000
    {
      auto p = MakeParams(5000.);
      std::unique_ptr<ThermalModelCanonical> m1(RunCanonical(TPS, p, GaussLegendre, 1));
      std::unique_ptr<ThermalModelCanonical> m2(RunCanonical(TPS, p, GaussLegendre, 2));
      auto [mis, nan] = CompareDensities(*m1, m2->Densities(), TPS, 2.e-4, 5, "conv ");
      EXPECT_EQ(mis, 0) << mis << " species exceed convergence tolerance";
      EXPECT_EQ(nan, 0) << nan << " species with NaN/Inf density (convergence)";
    }
  }

  // ---- Saddle-point tests (all three methods) --------------------------------

  // Test 3: SP and SP2 vs GL at large volume
  TEST_F(ThermalModelCanonicalTest, SaddlePointVsGL) {
    auto p = MakeParams(50000.);
    std::unique_ptr<ThermalModelCanonical> mGL(RunCanonical(TPS, p, GaussLegendre, 2));
    std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint));
    std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO));

    auto [mis1, nan1] = CompareDensities(*mSP,  mGL->Densities(), TPS, 1.e-3, 5, "SP ");
    auto [mis2, nan2] = CompareDensities(*mSP2, mGL->Densities(), TPS, 1.e-2, 5, "SP2 ");
    EXPECT_EQ(nan1, 0);
    EXPECT_EQ(nan2, 0);
    EXPECT_EQ(mis1, 0) << mis1 << " species: SP vs GL mismatch";
    EXPECT_EQ(mis2, 0) << mis2 << " species: SP2 vs GL mismatch";
  }

  // Test 4: SP GCE limit
  TEST_F(ThermalModelCanonicalTest, SaddlePointGCELimit) {
    auto p = MakeParams(1000000.);
    std::unique_ptr<ThermalModelCanonical> mCE(RunCanonical(TPS, p, SaddlePoint));

    ThermalModelIdeal mGCE(TPS, p);
    mGCE.SetStatistics(false);
    mGCE.SetUseWidth(ThermalParticle::ZeroWidth);
    mGCE.CalculateDensities();

    auto [mis, nan] = CompareDensities(*mCE, mGCE.Densities(), TPS, 1.e-3, 5, "SP GCE ");
    EXPECT_EQ(nan, 0);
    EXPECT_EQ(mis, 0) << mis << " species exceed GCE tolerance (SP)";
  }

  // Test 5: SP and SP2 charge conservation at B=2, Q=1
  TEST_F(ThermalModelCanonicalTest, SaddlePointChargeConservation) {
    auto p = MakeParams(1000., 2, 1, 0);

    // SP — conservation to O(1/V)
    std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint));
    EXPECT_NEAR(mSP->CalculateBaryonDensity()      * p.SVc, 2., 1.e-2);
    EXPECT_NEAR(mSP->CalculateChargeDensity()       * p.SVc, 1., 1.e-2);
    EXPECT_NEAR(mSP->CalculateStrangenessDensity()  * p.SVc, 0., 1.e-2);

    // SP2 — exact by construction
    std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO));
    double tB = 0., tQ = 0., tS = 0.;
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      tB += TPS->Particle(i).BaryonCharge()    * mSP2->Densities()[i] * p.V;
      tQ += TPS->Particle(i).ElectricCharge()   * mSP2->Densities()[i] * p.V;
      tS += TPS->Particle(i).Strangeness()      * mSP2->Densities()[i] * p.V;
    }
    EXPECT_NEAR(tB, 2., 1.e-6) << "SP2 baryon not conserved";
    EXPECT_NEAR(tQ, 1., 1.e-6) << "SP2 charge not conserved";
    EXPECT_NEAR(tS, 0., 1.e-6) << "SP2 strangeness not conserved";
  }

  // ---- Mixed-canonical tests -------------------------------------------------

  // Test 6: Mixed-canonical ensembles (BS-only, S-only, BS with B=2)
  TEST_F(ThermalModelCanonicalTest, MixedCanonical) {
    // (a) Canonical B,S only at V=50000, N=0: GL vs SP vs SP2
    {
      auto p = MakeParams(50000.);
      std::unique_ptr<ThermalModelCanonical> mGL(RunCanonical(TPS, p, GaussLegendre, 2, true, false, true, false));
      std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint, 1, true, false, true, false));
      std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO, 1, true, false, true, false));

      // Conservation
      EXPECT_NEAR(mSP2->CalculateBaryonDensity()     * p.SVc, 0., 1.e-6);
      EXPECT_NEAR(mSP2->CalculateStrangenessDensity() * p.SVc, 0., 1.e-6);

      // Density agreement
      auto [mis1, nan1] = CompareDensities(*mSP,  mGL->Densities(), TPS, 1.e-2, 3, "BS SP ");
      auto [mis2, nan2] = CompareDensities(*mSP2, mGL->Densities(), TPS, 1.e-2, 3, "BS SP2 ");
      EXPECT_EQ(mis1, 0) << mis1 << " species: SP vs GL (BS-only)";
      EXPECT_EQ(mis2, 0) << mis2 << " species: SP2 vs GL (BS-only)";
    }

    // (b) Canonical S only at muB=0: GL vs SP vs SP2
    {
      auto p = MakeParams(50000.);
      std::unique_ptr<ThermalModelCanonical> mGL(RunCanonical(TPS, p, GaussLegendre, 2, false, false, true, false));
      std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint, 1, false, false, true, false));
      std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO, 1, false, false, true, false));

      EXPECT_NEAR(mSP2->CalculateStrangenessDensity() * p.SVc, 0., 1.e-6);
      auto [mis1, nan1] = CompareDensities(*mSP,  mGL->Densities(), TPS, 1.e-2, 3, "S SP ");
      auto [mis2, nan2] = CompareDensities(*mSP2, mGL->Densities(), TPS, 1.e-2, 3, "S SP2 ");
      EXPECT_EQ(mis1, 0);
      EXPECT_EQ(mis2, 0);
    }

    // (c) Canonical S only at muB=0.3: SP vs SP2
    //     GL at V=50000 with muB=0.3 is marginal even with the contour shift
    //     (integrand peak is extremely narrow at this volume). See the
    //     GLContourShiftFiniteMuB test for GL vs SP comparison at moderate V.
    {
      auto p = MakeParams(50000., 0, 0, 0, 0, 0.3);
      std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint, 1, false, false, true, false));
      std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO, 1, false, false, true, false));

      EXPECT_NEAR(mSP->CalculateStrangenessDensity()  * p.SVc, 0., 1.e-2);
      EXPECT_NEAR(mSP2->CalculateStrangenessDensity() * p.SVc, 0., 1.e-6);
      auto [mis, nan] = CompareDensities(*mSP, mSP2->Densities(), TPS, 1.e-2, 3, "S muB=0.3 ");
      EXPECT_EQ(mis, 0) << mis << " species: SP vs SP2 (S-only, muB=0.3)";
    }

    // (d) Canonical B,S with B=2 at V=5000: all three methods
    {
      auto p = MakeParams(5000., 2, 0, 0);
      std::unique_ptr<ThermalModelCanonical> mGL(RunCanonical(TPS, p, GaussLegendre, 2, true, false, true, false));
      std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint, 1, true, false, true, false));
      std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO, 1, true, false, true, false));

      EXPECT_NEAR(mSP2->CalculateBaryonDensity()     * p.SVc, 2., 1.e-6);
      EXPECT_NEAR(mSP2->CalculateStrangenessDensity() * p.SVc, 0., 1.e-6);
      auto [mis1, nan1] = CompareDensities(*mSP,  mGL->Densities(), TPS, 1.e-2, 3, "BS B=2 SP ");
      auto [mis2, nan2] = CompareDensities(*mSP2, mGL->Densities(), TPS, 1.e-2, 3, "BS B=2 SP2 ");
      EXPECT_EQ(mis1, 0);
      EXPECT_EQ(mis2, 0);
    }
  }

  // ---- Charm tests -----------------------------------------------------------

  // Test 7: Charm conservation (GL), SP vs GL, charm-only canonical, SP consistency
  TEST_F(ThermalModelCanonicalCharmTest, CharmCanonical) {
    // (a) GL conservation: all four charges, V=5000
    {
      auto p = MakeParams(5000.);
      std::unique_ptr<ThermalModelCanonical> m(RunCanonical(TPS, p, GaussLegendre, 2));

      EXPECT_NEAR(m->CalculateBaryonDensity()      * p.SVc, 0., 1.e-3);
      EXPECT_NEAR(m->CalculateChargeDensity()       * p.SVc, 0., 1.e-3);
      EXPECT_NEAR(m->CalculateStrangenessDensity()  * p.SVc, 0., 1.e-3);
      EXPECT_NEAR(m->CalculateCharmDensity()        * p.SVc, 0., 1.e-3);

      // At least one charmed particle has nonzero density
      bool found = false;
      for (int i = 0; i < TPS->ComponentsNumber(); ++i)
        if (TPS->Particle(i).Charm() != 0 && m->Densities()[i] > 1.e-20)
          { found = true; break; }
      EXPECT_TRUE(found) << "No charmed particle with nonzero density";
    }

    // (b) SP vs GL with charm at small volume (R=7 fm, V≈1437 fm³, Nccbar≈0.1).
    //     Charm is exponentially rare here; this tests that the saddle-point solver
    //     handles ill-conditioned BQS-only sectors without diverging.
    {
      auto p = MakeParams(1437.);
      std::unique_ptr<ThermalModelCanonical> mGL(RunCanonical(TPS, p, GaussLegendre, 2));
      std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint));

      auto [mis, nan] = CompareDensities(*mSP, mGL->Densities(), TPS, 5.e-2, 3, "small-V charm SP ");
      EXPECT_EQ(nan, 0) << "NaN/Inf in SP densities at small charm volume";
      EXPECT_EQ(mis, 0) << mis << " species: SP vs GL (small-V charm)";
    }

    // (c) SP and SP2 vs GL with charm, large volume V=50000
    {
      auto p = MakeParams(50000.);
      std::unique_ptr<ThermalModelCanonical> mGL(RunCanonical(TPS, p, GaussLegendre, 2));
      std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint));
      std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO));

      EXPECT_NEAR(mSP->CalculateCharmDensity()  * p.SVc, 0., 1.e-2);
      EXPECT_NEAR(mSP2->CalculateCharmDensity() * p.SVc, 0., 1.e-6);

      auto [mis1, nan1] = CompareDensities(*mSP,  mGL->Densities(), TPS, 1.e-2, 3, "charm SP ");
      auto [mis2, nan2] = CompareDensities(*mSP2, mGL->Densities(), TPS, 1.e-2, 3, "charm SP2 ");
      EXPECT_EQ(nan1 + nan2, 0);
      EXPECT_EQ(mis1, 0) << mis1 << " species: SP vs GL (charm)";
      EXPECT_EQ(mis2, 0) << mis2 << " species: SP2 vs GL (charm)";
    }

    // (d) Canonical C only (GCE for B,Q,S), V=50000
    {
      auto p = MakeParams(50000.);
      std::unique_ptr<ThermalModelCanonical> mGL(RunCanonical(TPS, p, GaussLegendre, 2, false, false, false, true));
      std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint, 1, false, false, false, true));
      std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO, 1, false, false, false, true));

      EXPECT_NEAR(mGL->CalculateCharmDensity()  * p.SVc, 0., 1.e-3);
      EXPECT_NEAR(mSP2->CalculateCharmDensity() * p.SVc, 0., 1.e-6);

      auto [mis1, nan1] = CompareDensities(*mSP,  mGL->Densities(), TPS, 1.e-2, 3, "C-only SP ");
      auto [mis2, nan2] = CompareDensities(*mSP2, mGL->Densities(), TPS, 1.e-2, 3, "C-only SP2 ");
      EXPECT_EQ(mis1, 0);
      EXPECT_EQ(mis2, 0);
    }

    // (e) Full BQSC: SP vs SP2 agreement + SP2 conservation, V=10^6
    //     (charmed particles are too rare for GCE convergence at this V)
    {
      auto p = MakeParams(1000000.);
      std::unique_ptr<ThermalModelCanonical> mSP(RunCanonical(TPS, p, SaddlePoint));
      std::unique_ptr<ThermalModelCanonical> mSP2(RunCanonical(TPS, p, SaddlePointNLO));

      EXPECT_NEAR(mSP2->CalculateBaryonDensity()      * p.SVc, 0., 1.e-6);
      EXPECT_NEAR(mSP2->CalculateChargeDensity()       * p.SVc, 0., 1.e-6);
      EXPECT_NEAR(mSP2->CalculateStrangenessDensity()  * p.SVc, 0., 1.e-6);
      EXPECT_NEAR(mSP2->CalculateCharmDensity()        * p.SVc, 0., 1.e-6);

      auto [mis, nan] = CompareDensities(*mSP, mSP2->Densities(), TPS, 1.e-2, 3, "BQSC ");
      EXPECT_EQ(nan, 0);
      EXPECT_EQ(mis, 0) << mis << " species: SP vs SP2 (BQSC)";
    }
  }

  // Test 8: FixParameters for mixed canonical at finite muB using SaddlePoint.
  //
  // This tests the scenario where muQ needs to be constrained in a
  // strangeness-canonical setup at finite muB. The FixParameters method
  // uses a GCE solve to find stable chemical potentials (avoiding GL/SP
  // inside the Broyden loop). Then CalculatePrimordialDensities uses SP
  // for the final result.
  //
  // Note: GL is unreliable for S-canonical at finite muB and large volumes
  // because the integrand becomes highly oscillatory. SP is used instead.
  TEST_F(ThermalModelCanonicalTest, FixParametersMixedCanonical) {
    // Strangeness-canonical at finite muB, using SaddlePoint
    ThermalModelParameters p;
    p.T = 0.155;
    p.muB = 0.500;  // 500 MeV
    p.muS = 0.;  p.muQ = 0.;  p.muC = 0.;
    p.gammaS = 1.;  p.gammaq = 1.;
    p.V = 5000.;  p.SVc = 5000.;
    p.B = 0;  p.Q = 0;  p.S = 0;  p.C = 0;

    ThermalModelCanonical model(TPS, p);
    model.SetStatistics(false);
    model.SetUseWidth(ThermalParticle::ZeroWidth);
    model.SetMethod(SaddlePoint);
    // S canonical, B/Q/C grand-canonical
    model.ConserveBaryonCharge(false);
    model.ConserveElectricCharge(false);
    model.ConserveStrangeness(true);
    model.ConserveCharm(false);

    // Constrain muQ to satisfy Q/B = 0.4 (typical for heavy-ion collisions)
    model.ConstrainMuB(false);
    model.ConstrainMuQ(true);
    model.ConstrainMuS(false);  // S is canonical, muS is irrelevant
    model.ConstrainMuC(false);
    model.SetQoverB(0.4);

    // This triggers FixParameters → GCE solve for muQ, then SP calculation
    model.ConstrainChemicalPotentials(true);

    // Verify the calculation produced finite results
    double nb = model.CalculateBaryonDensity();
    double nq = model.CalculateChargeDensity();
    double ns = model.CalculateStrangenessDensity();
    double entropy = model.CalculateEntropyDensity();

    EXPECT_TRUE(std::isfinite(nb)) << "Baryon density is not finite";
    EXPECT_TRUE(std::isfinite(nq)) << "Charge density is not finite";
    EXPECT_TRUE(std::isfinite(ns)) << "Strangeness density is not finite";
    EXPECT_TRUE(std::isfinite(entropy)) << "Entropy density is not finite";

    // Strangeness should be conserved (net S = 0)
    EXPECT_NEAR(ns * p.SVc, 0., 1.e-2)
      << "Net strangeness not conserved: S = " << ns * p.SVc;

    // Q/B should be close to 0.4 (GCE-solved muQ should give ~0.4)
    EXPECT_TRUE(nb > 0.) << "Baryon density should be positive at muB=500 MeV";
    if (nb > 0.) {
      double qOverB = nq / nb;
      EXPECT_NEAR(qOverB, 0.4, 0.02)
        << "Q/B = " << qOverB << " (expected ~0.4)";
    }
  }

  // Test 9: GL convergence at finite muB after contour shift.
  //
  // Before the contour shift fix, GL integration failed for S-canonical
  // at finite muB and moderate volumes (R >= 3 fm) due to rapid phase
  // oscillations from strange baryon-antibaryon asymmetry.
  // The contour shift (using saddle-point mu_S in the GL integrand)
  // eliminates linear oscillations and makes GL converge.
  TEST_F(ThermalModelCanonicalTest, GLContourShiftFiniteMuB) {
    // S-canonical at T=155 MeV, muB=500 MeV, R=3 fm (V ~ 113 fm^3)
    // This previously failed completely with GL.
    double R = 3.0;
    double V = 4. / 3. * xMath::Pi() * R * R * R;
    auto p = MakeParams(V, 0, 0, 0, 0, 0.5);

    std::unique_ptr<ThermalModelCanonical> mGL(
      RunCanonical(TPS, p, GaussLegendre, 2, false, false, true, false));
    std::unique_ptr<ThermalModelCanonical> mSP(
      RunCanonical(TPS, p, SaddlePoint, 1, false, false, true, false));

    // GL should now produce finite results (previously gave NaN/garbage)
    double ns_GL = mGL->CalculateStrangenessDensity() * p.SVc;
    double ns_SP = mSP->CalculateStrangenessDensity() * p.SVc;
    double ed_GL = mGL->CalculateEnergyDensity();
    double ed_SP = mSP->CalculateEnergyDensity();

    EXPECT_TRUE(std::isfinite(ed_GL)) << "GL energy density not finite at muB=500";
    EXPECT_NEAR(ns_GL, 0., 1.e-2)
      << "GL strangeness not conserved: S = " << ns_GL;
    EXPECT_NEAR(ns_SP, 0., 1.e-2)
      << "SP strangeness not conserved: S = " << ns_SP;

    // GL and SP should agree on energy density (within ~1%)
    if (std::isfinite(ed_GL) && ed_SP > 0.)
      EXPECT_NEAR(ed_GL / ed_SP, 1., 5.e-2)
        << "GL/SP energy density mismatch: " << ed_GL << " vs " << ed_SP;

    // Compare particle densities: GL vs SP
    auto [mis, nan] = CompareDensities(*mSP, mGL->Densities(), TPS, 5.e-2, 3,
                                        "GL-shift muB=0.5 ");
    EXPECT_EQ(nan, 0) << "NaN in GL densities at muB=500";
    EXPECT_EQ(mis, 0) << mis << " species: GL vs SP (S-only, muB=0.5, R=3)";
  }

} // namespace
