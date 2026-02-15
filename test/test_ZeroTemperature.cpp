/*
 * Thermal-FIST package
 *
 * Copyright (c) 2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

 /**
  * @file test_ZeroTemperature.cpp
  *
  * Verifies that thermodynamic quantities computed at T = 0 are
  * consistent with the T -> 0 limit obtained from small but finite
  * temperatures (T = 1e-4, 1e-3, 1e-2 GeV).
  *
  * For each model (Ideal, EV Diagonal, VDW, RealGas) at fixed
  * muB = 1.2 GeV, muQ = muS = muC = 0, the test checks:
  *   - Pressure, energy density, baryon density approach T=0 values
  *   - Entropy density vanishes as T -> 0
  *   - Adiabatic speed of sound squared (B-conserved) approaches T=0 value
  *   - ds/dT|_mu (analytic Sommerfeld vs C_mu/T) converges monotonically
  *
  * The HRG contains light pions (~140 MeV), so even T = 10 MeV gives
  * O(1%) corrections.  At T = 0.1 MeV, numerical precision of the
  * integration routines starts to dominate.  Therefore we test at
  * T = {1e-2, 1e-3} GeV with generous tolerances and verify that
  * corrections shrink as T decreases.
  */
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "ThermalFISTConfig.h"
#include "HRGBase.h"
#include "HRGEV.h"
#include "HRGRealGas.h"
#include "HRGVDW.h"
#include "HRGEV/ExcludedVolumeHelper.h"

#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  // ----------------------------------------------------------------
  //  Test fixture
  // ----------------------------------------------------------------
  class ZeroTemperatureTest : public ::testing::Test {
  protected:
    void SetUp() override {
      std::string inputDir = std::string(ThermalFIST_INPUT_FOLDER);
      TPS = new ThermalParticleSystem(
        inputDir + "/list/PDG2025/list.dat",
        inputDir + "/list/PDG2025/decays.dat");
    }

    void TearDown() override {
      delete TPS;
    }

    // Quantities collected at a single (T, muB) point.
    struct ThermoPoint {
      double P;       // pressure
      double e;       // energy density
      double s;       // entropy density
      double rhoB;    // baryon density
      double cs2_B;   // adiabatic speed of sound (B conserved)
      double cT2_B;   // isothermal speed of sound (B conserved)
      double dsdT_mu; // ds/dT at fixed mu  (= C_mu / T  for T > 0)
    };

    // Evaluate all quantities for a given model at its current parameters.
    ThermoPoint Evaluate(ThermalModelBase* model) {
      model->CalculatePrimordialDensities();
      model->CalculateFluctuations();
      model->CalculateTemperatureDerivatives();

      ThermoPoint pt;
      pt.P    = model->Pressure();
      pt.e    = model->EnergyDensity();
      pt.s    = model->EntropyDensity();
      pt.rhoB = model->BaryonDensity();

      double cs2 = model->cs2(true, false, false, false);
      pt.cs2_B = std::isfinite(cs2) ? cs2 : std::numeric_limits<double>::quiet_NaN();

      double cT2 = model->cT2(true, false, false, false);
      pt.cT2_B = std::isfinite(cT2) ? cT2 : std::numeric_limits<double>::quiet_NaN();

      pt.dsdT_mu = model->CalculateEntropyDensityDerivativeT();
      return pt;
    }

    // Helper: relative difference, returns NaN if denominator is tiny.
    static double RelDiff(double val, double ref) {
      if (std::abs(ref) < 1.e-30)
        return std::numeric_limits<double>::quiet_NaN();
      return std::abs(val - ref) / std::abs(ref);
    }

    // Run the convergence check for one model.
    void RunConvergenceCheck(ThermalModelBase* model, const std::string& tag) {
      const double muB = 1.200;  // GeV, well above nucleon mass
      model->SetBaryonChemicalPotential(muB);
      model->SetElectricChemicalPotential(0.);
      model->SetStrangenessChemicalPotential(0.);
      model->SetCharmChemicalPotential(0.);

      // ----- T = 0 reference -----
      model->SetTemperature(0.);
      ThermoPoint ref = Evaluate(model);

      // Basic T = 0 sanity checks
      EXPECT_TRUE(std::isfinite(ref.P))    << tag << ": P(T=0) is not finite";
      EXPECT_TRUE(std::isfinite(ref.e))    << tag << ": e(T=0) is not finite";
      EXPECT_TRUE(std::isfinite(ref.rhoB)) << tag << ": rhoB(T=0) is not finite";
      EXPECT_GT(ref.P, 0.)    << tag << ": P(T=0) should be > 0 at large muB";
      EXPECT_GT(ref.e, 0.)    << tag << ": e(T=0) should be > 0 at large muB";
      EXPECT_GT(ref.rhoB, 0.) << tag << ": rhoB(T=0) should be > 0";
      EXPECT_NEAR(ref.s, 0., 1.e-12) << tag << ": s(T=0) should be 0";
      EXPECT_GT(ref.cs2_B, 0.)  << tag << ": cs2(T=0) should be > 0";
      EXPECT_GT(ref.cT2_B, 0.)  << tag << ": cT2(T=0) should be > 0";
      EXPECT_TRUE(std::isfinite(ref.dsdT_mu)) << tag << ": ds/dT(T=0) should be finite";
      EXPECT_GT(ref.dsdT_mu, 0.) << tag << ": ds/dT(T=0) should be > 0";

      // ----- Evaluate at small finite T -----
      // T = 10 MeV and T = 1 MeV.
      // The HRG has many species including light pions, so corrections
      // are non-negligible at T = 10 MeV but should shrink at T = 1 MeV.
      const double T1 = 1.e-2;  // 10 MeV
      const double T2 = 1.e-3;  // 1 MeV

      model->SetTemperature(T1);
      ThermoPoint pt1 = Evaluate(model);

      model->SetTemperature(T2);
      ThermoPoint pt2 = Evaluate(model);

      // --- Check 1: Finite-T results are valid ---
      EXPECT_TRUE(std::isfinite(pt1.P))    << tag << " T=" << T1 << ": P not finite";
      EXPECT_TRUE(std::isfinite(pt1.e))    << tag << " T=" << T1 << ": e not finite";
      EXPECT_TRUE(std::isfinite(pt1.rhoB)) << tag << " T=" << T1 << ": rhoB not finite";
      EXPECT_TRUE(std::isfinite(pt2.P))    << tag << " T=" << T2 << ": P not finite";
      EXPECT_TRUE(std::isfinite(pt2.e))    << tag << " T=" << T2 << ": e not finite";
      EXPECT_TRUE(std::isfinite(pt2.rhoB)) << tag << " T=" << T2 << ": rhoB not finite";

      // --- Check 2: Quantities at T = 10 MeV are within ~5% of T = 0 ---
      // (generous bound; pion thermal contributions are non-negligible)
      EXPECT_LT(RelDiff(pt1.P, ref.P), 0.05)
        << tag << " T=" << T1 << ": P too far from T=0 value";
      EXPECT_LT(RelDiff(pt1.e, ref.e), 0.05)
        << tag << " T=" << T1 << ": e too far from T=0 value";
      EXPECT_LT(RelDiff(pt1.rhoB, ref.rhoB), 0.05)
        << tag << " T=" << T1 << ": rhoB too far from T=0 value";

      // --- Check 3: Quantities at T = 1 MeV are within ~0.1% of T = 0 ---
      EXPECT_LT(RelDiff(pt2.P, ref.P), 1.e-3)
        << tag << " T=" << T2 << ": P too far from T=0 value";
      EXPECT_LT(RelDiff(pt2.e, ref.e), 1.e-3)
        << tag << " T=" << T2 << ": e too far from T=0 value";
      EXPECT_LT(RelDiff(pt2.rhoB, ref.rhoB), 1.e-3)
        << tag << " T=" << T2 << ": rhoB too far from T=0 value";

      // --- Check 4: Deviations shrink from T1 to T2 (convergence) ---
      // |Q(T2) - Q(0)| < |Q(T1) - Q(0)|  for each quantity Q.
      EXPECT_LT(std::abs(pt2.P - ref.P), std::abs(pt1.P - ref.P))
        << tag << ": P should converge monotonically toward T=0";
      EXPECT_LT(std::abs(pt2.e - ref.e), std::abs(pt1.e - ref.e))
        << tag << ": e should converge monotonically toward T=0";
      EXPECT_LT(std::abs(pt2.rhoB - ref.rhoB), std::abs(pt1.rhoB - ref.rhoB))
        << tag << ": rhoB should converge monotonically toward T=0";

      // --- Check 5: Entropy density is O(T) and shrinks ---
      EXPECT_GT(pt1.s, 0.) << tag << " T=" << T1 << ": s should be > 0";
      EXPECT_GT(pt2.s, 0.) << tag << " T=" << T2 << ": s should be > 0";
      EXPECT_LT(pt2.s, pt1.s)
        << tag << ": s should decrease with decreasing T";

      // --- Check 6: Speed of sound at T = 1 MeV within 1% of T = 0 ---
      if (std::isfinite(pt2.cs2_B) && ref.cs2_B > 0.) {
        EXPECT_LT(RelDiff(pt2.cs2_B, ref.cs2_B), 0.01)
          << tag << " T=" << T2 << ": cs2_B too far from T=0 value";
      }

      // --- Check 7: ds/dT convergence (Sommerfeld vs C_mu/T) ---
      // The leading Sommerfeld correction to C_mu/T is O(T^2), but
      // in the full HRG with many species including light pions the
      // coefficient can be large.  Use generous bounds.
      if (std::isfinite(ref.dsdT_mu) && ref.dsdT_mu > 0.) {
        if (std::isfinite(pt1.dsdT_mu)) {
          // At T = 10 MeV: within 20%
          EXPECT_LT(RelDiff(pt1.dsdT_mu, ref.dsdT_mu), 0.20)
            << tag << " T=" << T1 << ": ds/dT too far from Sommerfeld value";
        }
        if (std::isfinite(pt2.dsdT_mu)) {
          // At T = 1 MeV: within 1%
          EXPECT_LT(RelDiff(pt2.dsdT_mu, ref.dsdT_mu), 0.01)
            << tag << " T=" << T2 << ": ds/dT too far from Sommerfeld value";
        }
        // Convergence: T2 closer than T1
        if (std::isfinite(pt1.dsdT_mu) && std::isfinite(pt2.dsdT_mu)) {
          EXPECT_LT(std::abs(pt2.dsdT_mu - ref.dsdT_mu),
                     std::abs(pt1.dsdT_mu - ref.dsdT_mu) + 1.e-10)
            << tag << ": ds/dT should converge toward Sommerfeld value";
        }
      }

      // Restore T = 0 state
      model->SetTemperature(0.);
      model->CalculatePrimordialDensities();
    }

    ThermalParticleSystem* TPS;
  };

  // ----------------------------------------------------------------
  //  Ideal HRG
  // ----------------------------------------------------------------
  TEST_F(ZeroTemperatureTest, IdealHRG) {
    ThermalModelIdeal model(TPS);
    model.SetStatistics(true);
    model.SetUseWidth(false);
    model.SetCalculationType(IdealGasFunctions::Quadratures);
    RunConvergenceCheck(&model, "IdealHRG");
  }

  // ----------------------------------------------------------------
  //  EV Diagonal HRG
  // ----------------------------------------------------------------
  TEST_F(ZeroTemperatureTest, EVDiagonalHRG) {
    ThermalModelEVDiagonal model(TPS);
    model.SetStatistics(true);
    model.SetUseWidth(false);
    model.SetCalculationType(IdealGasFunctions::Quadratures);

    // Set baryon eigenvolume b = 1 fm^3 for all baryon-baryon pairs
    double b = 1.0; // fm^3
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS->ComponentsNumber(); ++j) {
        if (TPS->Particle(i).BaryonCharge() != 0 && TPS->Particle(j).BaryonCharge() != 0)
          model.SetVirial(i, j, b);
      }
    }

    RunConvergenceCheck(&model, "EVDiagHRG");
  }

  // ----------------------------------------------------------------
  //  QvdW HRG
  // ----------------------------------------------------------------
  TEST_F(ZeroTemperatureTest, QvdWHRG) {
    ThermalModelVDW model(TPS);
    model.SetStatistics(true);
    model.SetUseWidth(false);
    model.SetCalculationType(IdealGasFunctions::Quadratures);
    SetVDWHRGInteractionParameters(&model, 0.329, 3.42);

    RunConvergenceCheck(&model, "QvdWHRG");
  }

  // ----------------------------------------------------------------
  //  RealGas HRG  (Carnahan-Starling EV + VDW mean field)
  // ----------------------------------------------------------------
  TEST_F(ZeroTemperatureTest, RealGasHRG) {
    ThermalModelRealGas model(TPS);
    model.SetStatistics(true);
    model.SetUseWidth(false);
    model.SetCalculationType(IdealGasFunctions::Quadratures);

    double a = 0.329;  // GeV fm^3
    double b = 3.42;   // fm^3

    // Carnahan-Starling excluded volume
    {
      ExcludedVolumeModelBase* evmod = new ExcludedVolumeModelCS();
      auto bij = GetBaryonBaryonInteractionMatrix(model.TPS(), b);
      model.SetExcludedVolumeModel(
        new ExcludedVolumeModelCrosstermsGeneralized(evmod, bij));
    }
    // VDW mean field
    {
      auto aij = GetBaryonBaryonInteractionMatrix(model.TPS(), a);
      model.SetMeanFieldModel(new MeanFieldModelMultiVDW(aij));
    }

    RunConvergenceCheck(&model, "RealGasHRG");
  }

}  // namespace
