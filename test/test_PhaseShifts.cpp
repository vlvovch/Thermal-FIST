/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <set>
#include "ThermalFISTConfig.h"
#include "HRGPhaseShifts/PhaseShiftDensity.h"
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"
#include "HRGPhaseShifts/MesonBaryonPhaseShifts.h"
#include "HRGPhaseShifts/PhaseShiftModel.h"
#include "HRGPhaseShifts/PhaseShiftChannel.h"
#include "HRGBase/IdealGasFunctions.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGBase/ThermalModelCanonical.h"
#include "gtest/gtest.h"

// M_PI is a POSIX extension, not part of the C++ standard. MSVC does not
// define it in <cmath> unless _USE_MATH_DEFINES is set before the first
// inclusion of <cmath>, which cannot be guaranteed across transitive
// includes. Provide a local fallback that is safe on every toolchain.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace thermalfist;

namespace {

  // Build the pi-pi I=2 channel S-matrix d.o.f. (repulsive, S+D waves).
  // Uses the catalog waves, which carry the exact analytic d(delta)/dM.
  PhaseShiftDensity makePiPiI2(int nodes = 64) {
    return PhaseShiftDensity(PhaseShifts::PiPi_I2_Waves(),
                             PhaseShifts::PionMass(), PhaseShifts::PionMass(),
                             PhaseShifts::PiPiI2_Mmax(), /*Bose*/ -1, nodes);
  }

  // Same channel but WITHOUT analytic derivatives -> engine uses the
  // central-finite-difference fallback for d(delta)/dq.
  PhaseShiftDensity makePiPiI2FiniteDiff(int nodes = 64) {
    std::vector<PhaseShiftPartialWave> waves;
    waves.push_back(PhaseShiftPartialWave(1, &PhaseShifts::PiPi_delta_I2_S)); // no ddelta_dM
    waves.push_back(PhaseShiftPartialWave(5, &PhaseShifts::PiPi_delta_I2_D)); // no ddelta_dM
    return PhaseShiftDensity(waves,
                             PhaseShifts::PionMass(), PhaseShifts::PionMass(),
                             PhaseShifts::PiPiI2_Mmax(), /*Bose*/ -1, nodes);
  }

  // chi_2^Q contribution of the pi-pi I=2 multiplet at mu_Q = 0:
  // Iz = +-2 (Q=+-2), +-1 (Q=+-1), 0 (Q=0) -> sum Q^2 = 2*(4+1) = 10 times chi2(channel).
  double chi2QPiPiI2(PhaseShiftDensity& psd, double T) {
    return 10.0 * psd.Quantity(IdealGasFunctions::chi2, T, 0.0);
  }

  TEST(PhaseShifts, PiPiI2IsRepulsive) {
    PhaseShiftDensity psd = makePiPiI2();
    // The I=2 phase shift is negative and decreasing -> spectral weight < 0,
    // and it must be finite everywhere on (0, qmax) including near threshold.
    for (double q : {0.02, 0.05, 0.1, 0.2, 0.3}) {
      double w = psd.SpectralWeight(q);
      EXPECT_TRUE(std::isfinite(w)) << "q=" << q;
      EXPECT_LT(w, 0.0) << "q=" << q;   // repulsive
    }
  }

  TEST(PhaseShifts, PhaseShiftSignAndThreshold) {
    // delta_I2 < 0 above threshold; exactly 0 at/below threshold.
    double Mthr = 2.0 * PhaseShifts::PionMass();
    EXPECT_DOUBLE_EQ(PhaseShifts::PiPi_delta_I2_S(Mthr), 0.0);
    EXPECT_LT(PhaseShifts::PiPi_delta_I2_S(0.5), 0.0);
    EXPECT_LT(PhaseShifts::PiPi_delta_I2_S(0.8), 0.0);
    EXPECT_LT(PhaseShifts::PiPi_delta_I2_D(0.8), 0.0);
  }

  TEST(PhaseShifts, Chi2QContributionIsNegative) {
    // A repulsive channel suppresses electric-charge fluctuations.
    PhaseShiftDensity psd = makePiPiI2();
    double c2 = chi2QPiPiI2(psd, 0.150);
    EXPECT_TRUE(std::isfinite(c2));
    EXPECT_LT(c2, 0.0);
  }

  TEST(PhaseShifts, AnalyticMatchesFiniteDiff) {
    // The analytic d(delta)/dM overload and the finite-difference fallback must
    // agree to ~FD accuracy across the whole q range (away from the q=0 kink).
    PhaseShiftDensity psA = makePiPiI2();           // analytic
    PhaseShiftDensity psF = makePiPiI2FiniteDiff(); // finite diff
    for (double q : {0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7}) {
      double wA = psA.SpectralWeight(q);
      double wF = psF.SpectralWeight(q);
      EXPECT_NEAR(wA, wF, 1e-6 * std::fabs(wA) + 1e-9) << "q=" << q;
    }
    // and the integrated thermodynamics must match to quadrature precision.
    double cA = chi2QPiPiI2(psA, 0.150);
    double cF = chi2QPiPiI2(psF, 0.150);
    EXPECT_NEAR(cA, cF, 1e-6 * std::fabs(cA));
  }

  TEST(PhaseShifts, SpectralWeightReferenceValues) {
    // Analytic spectral weight (1/pi) sum_l (2J_l+1) ddelta_l/dq vs an
    // independent mpmath high-precision reference (S+D, charged-pion mass).
    PhaseShiftDensity psd = makePiPiI2();
    EXPECT_NEAR(psd.SpectralWeight(0.05), -1.430522040942e-01, 1e-9);
    EXPECT_NEAR(psd.SpectralWeight(0.20), -3.366924650159e-01, 1e-9);
    EXPECT_NEAR(psd.SpectralWeight(0.50), -4.885358082049e-01, 1e-9);
  }

  TEST(PhaseShifts, QuadratureConvergence) {
    // The q-integral must be converged: 32 -> 64 -> 128 nodes should agree well.
    double T = 0.150;
    PhaseShiftDensity p32 = makePiPiI2(32);
    PhaseShiftDensity p64 = makePiPiI2(64);
    PhaseShiftDensity p128 = makePiPiI2(128);
    double c32  = chi2QPiPiI2(p32,  T);
    double c64  = chi2QPiPiI2(p64,  T);
    double c128 = chi2QPiPiI2(p128, T);
    ASSERT_LT(std::fabs(c128), 1.0);          // sanity: O(<1) contribution
    EXPECT_NEAR(c64,  c128, 1e-4 * std::fabs(c128));
    EXPECT_NEAR(c32,  c128, 1e-2 * std::fabs(c128));
  }

  // ---- channel loader: PDG-id scheme + isospin expansion + model routing ----

  TEST(PhaseShifts, PdgIdSchemeFollowsPdgConventions) {
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I2_Channel();
    // PDG-aligned scheme: 99 F (2I) (2|Iz|) nn (2J+1), last digit = 2J+1.
    // pi-pi I=2: family 1, 2I=4. S-wave 2J+1=1, D-wave 2J+1=5.
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 4, 1), 99144001LL); // Iz=+2, S-wave
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 4, 5), 99144005LL); // Iz=+2, D-wave
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 2, 1), 99142001LL); // Iz=+1, S-wave
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 0, 5), 99140005LL); // Iz=0,  D-wave
    // Magnitude uses |2Iz| (antiparticle uses the negative of the whole id).
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, -4, 5), 99144005LL);
    // Last digit is a genuine 2J+1, per PDG.
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 4, 1) % 10, 1LL);
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 4, 5) % 10, 5LL);
    // 8-digit 99-prefix, clear of real PDG codes (<=7-digit hadrons, 10-digit nuclei).
    long long id = PhaseShifts::PhaseShiftPdgId(ch, 4, 5);
    EXPECT_EQ(id / 1000000LL,        99LL);   // reserved S-matrix block
    EXPECT_EQ((id / 100000LL) % 10,  1LL);    // family = pi pi
    EXPECT_EQ((id / 10000LL)  % 10,  4LL);    // 2I
    EXPECT_EQ((id / 1000LL)   % 10,  4LL);    // 2|Iz|
    EXPECT_EQ((id / 10LL)     % 100, 0LL);    // reserved excitation slot
    EXPECT_GE(id, 10000000LL);                // 8 digits
    EXPECT_LT(id, 100000000LL);
    // Charge/baryon/strangeness are NOT in the code (1 wouldn't matter); they
    // are carried by the particle fields (checked in the expansion test).
  }

  // Write a minimal base list with pi+ and pi0 (pi- auto-generated) for the
  // module's decays to resolve against. New (whitespace) list format:
  //   pdgid name stable mass deg stat B Q S C |S| |C| width threshold
  void writePionList(const std::string& path) {
    std::ofstream f(path.c_str());
    f << "# base pions\n";
    f << "211 pi+ 1 0.13957  1 -1 0 1 0 0 0 0 0 0\n";
    f << "111 pi0 1 0.134977 1 -1 0 0 0 0 0 0 0 0\n";
  }

  // Pions + kaons (anti-K auto-generated) for the pi-K module's decays.
  void writePiKList(const std::string& path) {
    std::ofstream f(path.c_str());
    f << "# base pions + kaons\n";
    f << "211 pi+ 1 0.13957  1 -1 0  1 0 0 0 0 0 0\n";
    f << "111 pi0 1 0.134977 1 -1 0  0 0 0 0 0 0 0\n";
    f << "321 K+  1 0.493677 1 -1 0  1 1 0 1 0 0 0\n";
    f << "311 K0  1 0.497611 1 -1 0  0 1 0 1 0 0 0\n";
  }

  // Pions + nucleons + the Delta(1232) multiplet (antiparticles auto-generated)
  // for the pi-N module. Baryons carry stat = +1 (Fermi), B = +1.
  void writePiNList(const std::string& path) {
    std::ofstream f(path.c_str());
    f << "# base pions + nucleons + Delta(1232)\n";
    f << "211  pi+          1 0.13957  1 -1 0  1 0 0 0 0 0 0\n";
    f << "111  pi0          1 0.134977 1 -1 0  0 0 0 0 0 0 0\n";
    f << "2212 p            1 0.938272 2  1 1  1 0 0 0 0 0 0\n";
    f << "2112 n            1 0.939565 2  1 1  0 0 0 0 0 0 0\n";
    f << "2224 Delta(1232)++ 0 1.232   4  1 1  2 0 0 0 0 0.117 1.078\n";
    f << "2214 Delta(1232)+  0 1.232   4  1 1  1 0 0 0 0 0.117 1.075\n";
    f << "2114 Delta(1232)0  0 1.232   4  1 1  0 0 0 0 0 0.117 1.076\n";
    f << "1114 Delta(1232)-  0 1.232   4  1 1 -1 0 0 0 0 0.117 1.079\n";
  }

  TEST(PhaseShifts, ClebschGordanDecays) {
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I2_Channel();
    // Iz=+2 -> pi+ pi+
    auto d2 = PhaseShifts::ChannelDecays(ch, 4);
    ASSERT_EQ(d2.size(), 1u);
    EXPECT_NEAR(d2[0].first, 1.0, 1e-12);
    EXPECT_EQ(d2[0].second.first, 211); EXPECT_EQ(d2[0].second.second, 211);
    // Iz=+1 -> pi+ pi0 (sorted daughters 111, 211)
    auto d1 = PhaseShifts::ChannelDecays(ch, 2);
    ASSERT_EQ(d1.size(), 1u);
    EXPECT_NEAR(d1[0].first, 1.0, 1e-12);
    EXPECT_EQ(d1[0].second.first, 111); EXPECT_EQ(d1[0].second.second, 211);
    // Iz=0 -> 1/3 pi+ pi- + 2/3 pi0 pi0
    auto d0 = PhaseShifts::ChannelDecays(ch, 0);
    ASSERT_EQ(d0.size(), 2u);
    double sum = 0., brPM = 0., br00 = 0.;
    for (size_t i = 0; i < d0.size(); ++i) {
      sum += d0[i].first;
      if (d0[i].second.first == -211 && d0[i].second.second == 211) brPM = d0[i].first;
      if (d0[i].second.first ==  111 && d0[i].second.second == 111) br00 = d0[i].first;
    }
    EXPECT_NEAR(sum, 1.0, 1e-12);
    EXPECT_NEAR(brPM, 1.0 / 3.0, 1e-9);
    EXPECT_NEAR(br00, 2.0 / 3.0, 1e-9);
  }

  TEST(PhaseShifts, ListBuildAttachChi2Q) {
    // Build via the list path (LoadList generates antiparticles + decays), attach
    // the chosen model, and verify the channel's contribution to Susc(Q,Q)
    // equals the direct 10 * chi2 and is negative.
    const double T = 0.150;
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_a.dat";
    writePionList(pionF);

    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I2_Channel();
    ch.name = "pipi_I2_a";    // unique file names (PDG ids unchanged: family+2I)
    auto model = PhaseShifts::PiPi_I2_Waves();
    const std::string listF = dir + "list-" + ch.name + ".dat";
    const std::string decF  = dir + "decays-" + ch.name + ".dat";
    PhaseShifts::WritePhaseShiftListFile(ch, model, listF);
    PhaseShifts::WritePhaseShiftDecaysFile(ch, model, decF);

    auto suscQ = [&](ThermalParticleSystem& TPS) {
      ThermalModelIdeal m(&TPS);
      m.SetTemperature(T);
      m.SetBaryonChemicalPotential(0.0);
      m.SetElectricChemicalPotential(0.0);
      m.SetStrangenessChemicalPotential(0.0);
      m.CalculatePrimordialDensities();
      m.CalculateFluctuations();
      return m.Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    };

    ThermalParticleSystem TPSbase(std::vector<std::string>(1, pionF));
    double base = suscQ(TPSbase);

    std::vector<std::string> lists; lists.push_back(pionF); lists.push_back(listF);
    ThermalParticleSystem TPS(lists, std::vector<std::string>(1, decF));
    PhaseShifts::AttachDensities(TPS, ch, model);
    PhaseShifts::SubsumeResonances(TPS, ch);   // no-op for I=2
    double full = suscQ(TPS);

    PhaseShiftDensity psd = makePiPiI2(64);
    double direct = 10.0 * psd.Quantity(IdealGasFunctions::chi2, T, 0.0);

    EXPECT_LT(full - base, 0.0);
    ASSERT_NE(direct, 0.0);
    EXPECT_NEAR(full - base, direct, 1e-6 * std::fabs(direct));

    std::remove(pionF.c_str()); std::remove(listF.c_str()); std::remove(decF.c_str());
  }

  TEST(PhaseShifts, FeeddownReducesPionDensity) {
    // The repulsive clusters decay to pions, so their (negative) feeddown must
    // reduce the total pion density relative to the free pion gas.
    const double T = 0.150;
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_b.dat";
    writePionList(pionF);

    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I2_Channel();
    ch.name = "pipi_I2_b";
    auto model = PhaseShifts::PiPi_I2_Waves();
    const std::string listF = dir + "list-" + ch.name + ".dat";
    const std::string decF  = dir + "decays-" + ch.name + ".dat";
    PhaseShifts::WritePhaseShiftListFile(ch, model, listF);
    PhaseShifts::WritePhaseShiftDecaysFile(ch, model, decF);

    auto totalPiPlus = [&](ThermalParticleSystem& TPS, bool attach) {
      ThermalModelIdeal m(&TPS);
      m.SetTemperature(T);
      m.SetBaryonChemicalPotential(0.0);
      m.SetElectricChemicalPotential(0.0);
      m.SetStrangenessChemicalPotential(0.0);
      m.CalculateDensities();
      return m.GetDensity(211, Feeddown::StabilityFlag);
    };

    ThermalParticleSystem TPSbase(std::vector<std::string>(1, pionF));
    double piBase = totalPiPlus(TPSbase, false);

    std::vector<std::string> lists; lists.push_back(pionF); lists.push_back(listF);
    ThermalParticleSystem TPS(lists, std::vector<std::string>(1, decF));
    PhaseShifts::AttachDensities(TPS, ch, model);
    double piFull = totalPiPlus(TPS, true);

    EXPECT_GT(piBase, 0.0);
    EXPECT_LT(piFull, piBase);    // repulsion -> negative feeddown -> fewer pions

    std::remove(pionF.c_str()); std::remove(listF.c_str()); std::remove(decF.c_str());
  }

  TEST(PhaseShifts, AnalyticVsTabulated) {
    // Tabulate the Garcia-Martin delta(M) for the S and D waves, build a tabulated
    // model (cubic-spline), and check it reproduces the analytic spectral weight
    // and chi2 closely.
    const std::string dir = ::testing::TempDir();
    const std::string sF = dir + "ps_tab_S.dat", dF = dir + "ps_tab_D.dat";
    const double mpi = PhaseShifts::PionMass();
    const double thr = 2. * mpi, Mmax = PhaseShifts::PiPiI2_Mmax();
    const int N = 2000;
    {
      std::ofstream fs(sF.c_str()), fd(dF.c_str());
      fs << "# M delta_S\n"; fd << "# M delta_D\n";
      for (int i = 0; i <= N; ++i) {
        double M = thr + (Mmax - thr) * i / N;
        fs << M << " " << PhaseShifts::PiPi_delta_I2_S(M) << "\n";
        fd << M << " " << PhaseShifts::PiPi_delta_I2_D(M) << "\n";
      }
    }
    std::vector<PhaseShiftPartialWave> tab;
    tab.push_back(PhaseShifts::TabulatedWave(1, sF));   // S-wave (2J+1=1)
    tab.push_back(PhaseShifts::TabulatedWave(5, dF));   // D-wave (2J+1=5)

    PhaseShiftDensity ana(PhaseShifts::PiPi_I2_Waves(), mpi, mpi, Mmax, -1, 64);
    PhaseShiftDensity tabd(tab,                  mpi, mpi, Mmax, -1, 64);

    // Sample only within the parametrization range (qmax ~ 0.696 for Mmax=1.42);
    // the small residual is the cubic-spline derivative error near the threshold cusp.
    for (double q : {0.1, 0.2, 0.4, 0.6}) {
      double wa = ana.SpectralWeight(q), wt = tabd.SpectralWeight(q);
      EXPECT_NEAR(wt, wa, 5e-3 * std::fabs(wa) + 1e-4) << "q=" << q;
    }
    double ca = ana.Quantity(IdealGasFunctions::chi2, 0.150, 0.0);
    double ct = tabd.Quantity(IdealGasFunctions::chi2, 0.150, 0.0);
    EXPECT_NEAR(ct, ca, 2e-3 * std::fabs(ca));

    std::remove(sF.c_str()); std::remove(dF.c_str());
  }

  // ---- high-level convenience API ----

  TEST(PhaseShifts, AddPhaseShiftChannelOneCall) {
    // The one-call in-memory helper must reproduce the file-based result: the
    // channel's Susc(Q,Q) contribution equals 10*chi2, and its (repulsive)
    // feeddown reduces the total pion density.
    const double T = 0.150;
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_c.dat";
    writePionList(pionF);

    ThermalParticleSystem TPSbase(std::vector<std::string>(1, pionF));
    ThermalModelIdeal mb(&TPSbase);
    mb.SetTemperature(T); mb.SetBaryonChemicalPotential(0.0);
    mb.SetElectricChemicalPotential(0.0); mb.SetStrangenessChemicalPotential(0.0);
    mb.CalculateDensities(); mb.CalculateFluctuations();
    double base = mb.Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    double piBase = mb.GetDensity(211, Feeddown::StabilityFlag);

    ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I2_Channel(),
                                      PhaseShifts::PiPi_I2_Waves());
    ThermalModelIdeal m(&TPS);
    m.SetTemperature(T); m.SetBaryonChemicalPotential(0.0);
    m.SetElectricChemicalPotential(0.0); m.SetStrangenessChemicalPotential(0.0);
    m.CalculateDensities(); m.CalculateFluctuations();
    double full = m.Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    double piFull = m.GetDensity(211, Feeddown::StabilityFlag);

    PhaseShiftDensity psd = makePiPiI2(64);
    double direct = 10.0 * psd.Quantity(IdealGasFunctions::chi2, T, 0.0);
    EXPECT_NEAR(full - base, direct, 1e-6 * std::fabs(direct));
    EXPECT_GT(piBase, 0.0);
    EXPECT_LT(piFull, piBase);    // repulsion -> negative feeddown -> fewer pions

    std::remove(pionF.c_str());
  }

  TEST(PhaseShifts, AddChannelAttachesToExistingClusters) {
    // If the clusters are already in the list as plain particles (e.g. loaded
    // from a list-<name>.dat file with no densities), AddPhaseShiftChannel must
    // attach densities to them rather than duplicate-add - so a phase-shift list
    // loaded directly still ends up with the negative Beth-Uhlenbeck densities.
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_f.dat";
    writePionList(pionF);

    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I2_Channel();
    ch.name = "pipi_I2_f";
    auto model = PhaseShifts::PiPi_I2_Waves();
    const std::string listF = dir + "list-" + ch.name + ".dat";
    const std::string decF  = dir + "decays-" + ch.name + ".dat";
    PhaseShifts::WritePhaseShiftListFile(ch, model, listF);
    PhaseShifts::WritePhaseShiftDecaysFile(ch, model, decF);

    // load base + raw module files (clusters present, but WITHOUT densities)
    std::vector<std::string> lists; lists.push_back(pionF); lists.push_back(listF);
    ThermalParticleSystem TPS(lists, std::vector<std::string>(1, decF));
    ASSERT_GE(TPS.PdgToId(99144005), 0);                        // cluster is present
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 0);   // but carries no density
    double idealdens = TPS.ParticleByPDG(99144005).GetGeneralizedDensity() == nullptr ? 1. : -1.;
    EXPECT_GT(idealdens, 0.0);                                  // no GeneralizedDensity yet
    const int nbefore = TPS.ComponentsNumber();

    // now attach via the channel: must not duplicate, must bind densities
    PhaseShifts::AddPhaseShiftChannel(TPS, ch, model);
    EXPECT_EQ(TPS.ComponentsNumber(), nbefore);                 // no duplicates added
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 10);  // densities now attached
    ASSERT_NE(TPS.ParticleByPDG(99144005).GetGeneralizedDensity(), nullptr);
    double c2 = TPS.ParticleByPDG(99144005).GetGeneralizedDensity()
                  ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0);
    EXPECT_LT(c2, 0.0);                                         // repulsive -> negative

    std::remove(pionF.c_str()); std::remove(listF.c_str()); std::remove(decF.c_str());
  }

  TEST(PhaseShifts, SurvivesClearDensityModels) {
    // ClearDensityModels() resets the EMM pion-BEC models and is called by the
    // GUI on every calculation; it must NOT delete phase-shift densities, or the
    // channels silently revert to ideal mesons (positive density, dead toggle).
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_g.dat";
    writePionList(pionF);

    ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I2_Channel(),
                                      PhaseShifts::PiPi_I2_Waves());
    ThermalModelIdeal model(&TPS);
    ASSERT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 10);

    model.ClearDensityModels();   // what SetThermalModelConfiguration does

    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 10);   // must survive
    ASSERT_NE(TPS.ParticleByPDG(99144005).GetGeneralizedDensity(), nullptr);
    double c2 = TPS.ParticleByPDG(99144005).GetGeneralizedDensity()
                  ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0);
    EXPECT_LT(c2, 0.0);

    std::remove(pionF.c_str());
  }

  TEST(PhaseShifts, EnableDisableToggle) {
    // Disabling the channels (without removing them) must zero their contribution
    // exactly, and re-enabling must restore it bit-for-bit.
    const double T = 0.150;
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_e.dat";
    writePionList(pionF);

    ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I2_Channel(),
                                      PhaseShifts::PiPi_I2_Waves());
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 10);

    auto suscQ = [&](ThermalParticleSystem& T_) {
      ThermalModelIdeal m(&T_);
      m.SetTemperature(T); m.SetBaryonChemicalPotential(0.0);
      m.SetElectricChemicalPotential(0.0); m.SetStrangenessChemicalPotential(0.0);
      m.CalculatePrimordialDensities(); m.CalculateFluctuations();
      return m.Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    };

    double on = suscQ(TPS);
    PhaseShifts::SetPhaseShiftsEnabled(TPS, false);
    double off = suscQ(TPS);
    PhaseShifts::SetPhaseShiftsEnabled(TPS, true);
    double on2 = suscQ(TPS);

    ThermalParticleSystem TPSbase(std::vector<std::string>(1, pionF));
    double base = suscQ(TPSbase);

    EXPECT_LT(on, off);                       // enabled repulsive channel lowers chi2Q
    EXPECT_DOUBLE_EQ(on, on2);                // toggling back restores exactly
    EXPECT_NEAR(off, base, 1e-12 * std::fabs(base) + 1e-15);  // disabled == pion-only

    std::remove(pionF.c_str());
  }

  TEST(PhaseShifts, AddPhaseShiftChannelsFromConfig) {
    // A single config file ("<channel>:<wave> <list> <decays> <model>") wires
    // everything in one call by loading the per-wave .dat files it references.
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_d.dat";
    writePionList(pionF);

    // The per-wave list/decay module files the conf will reference (one per wave).
    PhaseShifts::WritePhaseShiftFiles(PhaseShifts::PiPi_I2_Channel(),
                                      PhaseShifts::PiPi_I2_Waves(), dir);

    const std::string confF = dir + "ps_config.conf";
    {
      std::ofstream f(confF.c_str());
      f << "# wave      list                 decays                model\n";
      f << "pipi_I2:S  list-pipi_I2_S.dat  decays-pipi_I2_S.dat  GarciaMartin2011_S\n";
      f << "pipi_I2:D  list-pipi_I2_D.dat  decays-pipi_I2_D.dat  GarciaMartin2011_D\n";
    }

    ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
    std::vector<long long> added = PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, confF);
    EXPECT_EQ(added.size(), 10u);                    // 5 Iz x 2 waves (signed)
    // both the S (last digit 1) and D (last digit 5) clusters present, each with
    // its density attached - loaded from the referenced per-wave .dat files
    ASSERT_GE(TPS.PdgToId(99144001), 0);             // Iz=+2 S-wave
    ASSERT_GE(TPS.PdgToId(99144005), 0);             // Iz=+2 D-wave
    EXPECT_NE(TPS.ParticleByPDG(99144001).GetGeneralizedDensity(), nullptr);
    EXPECT_NE(TPS.ParticleByPDG(99144005).GetGeneralizedDensity(), nullptr);
    EXPECT_NE(TPS.ParticleByPDG(-99144005).GetGeneralizedDensity(), nullptr);
    // repulsive -> negative chi2 on a loaded cluster
    EXPECT_LT(TPS.ParticleByPDG(99144005).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);
    // decays loaded from file: Iz=+2 D-wave -> pi+ pi+ (211 211)
    ASSERT_EQ(TPS.ParticleByPDG(99144005).Decays().size(), 1u);
    // spectroscopic naming: name ends with the wave letter, not "2J+1=..."
    EXPECT_NE(TPS.ParticleByPDG(99144001).Name().find("S]"), std::string::npos);
    EXPECT_NE(TPS.ParticleByPDG(99144005).Name().find("D]"), std::string::npos);

    std::remove(pionF.c_str()); std::remove(confF.c_str());
    std::remove((dir + "list-pipi_I2_S.dat").c_str());
    std::remove((dir + "list-pipi_I2_D.dat").c_str());
    std::remove((dir + "decays-pipi_I2_S.dat").c_str());
    std::remove((dir + "decays-pipi_I2_D.dat").c_str());
  }

  TEST(PhaseShifts, ConfigWaveModelMismatchThrows) {
    // The analytic model name is per-wave; if the wave key (col 1) and the model's
    // wave disagree (e.g. ":S" but "GarciaMartin2011_D"), loading must throw.
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_mm.dat";
    writePionList(pionF);
    PhaseShifts::WritePhaseShiftFiles(PhaseShifts::PiPi_I2_Channel(),
                                      PhaseShifts::PiPi_I2_Waves(), dir);
    const std::string confF = dir + "ps_config_mm.conf";
    {
      std::ofstream f(confF.c_str());
      f << "pipi_I2:S  list-pipi_I2_S.dat  decays-pipi_I2_S.dat  GarciaMartin2011_D\n";
    }

    ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
    EXPECT_THROW(PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, confF), std::invalid_argument);

    std::remove(pionF.c_str()); std::remove(confF.c_str());
    std::remove((dir + "list-pipi_I2_S.dat").c_str());
    std::remove((dir + "list-pipi_I2_D.dat").c_str());
    std::remove((dir + "decays-pipi_I2_S.dat").c_str());
    std::remove((dir + "decays-pipi_I2_D.dat").c_str());
  }

  // ---- pi-K (strange, non-self-conjugate multiplet) ----

  TEST(PhaseShifts, PiKDeltaSigns) {
    // I=3/2 repulsive (delta < 0), I=1/2 attractive (delta > 0, the kappa).
    EXPECT_DOUBLE_EQ(PhaseShifts::PiK_delta_I32_S(PhaseShifts::PionMass()
                                                + PhaseShifts::KaonMass()), 0.0);  // threshold
    EXPECT_LT(PhaseShifts::PiK_delta_I32_S(0.9), 0.0);
    EXPECT_GT(PhaseShifts::PiK_delta_I12_S(0.9), 0.0);
    // I=1/2 is elastic only below the K-eta threshold -> 0 above it
    EXPECT_DOUBLE_EQ(PhaseShifts::PiK_delta_I12_S(1.2), 0.0);
  }

  TEST(PhaseShifts, PiKPdgIdStrangeMultiplet) {
    // Strange (S=+1) -> not self-conjugate: every Iz is a distinct member,
    // the Iz field encodes I+Iz in 0..2I (antiparticle = the negative code).
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiK_I32_Channel();
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, -3, 1), 99230001LL);  // Iz=-3/2
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, -1, 1), 99231001LL);  // Iz=-1/2
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, +1, 1), 99232001LL);  // Iz=+1/2
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, +3, 1), 99233001LL);  // Iz=+3/2
    EXPECT_EQ((99233001LL / 100000LL) % 10, 2LL);                    // family = pi-K
    EXPECT_EQ((99233001LL / 10000LL)  % 10, 3LL);                    // 2I = 3
    // I=1/2 reuses the real kappa/K0*(700) codes (subsumption by coincidence)
    PhaseShifts::PhaseShiftChannel c12 = PhaseShifts::PiK_I12_Channel();
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(c12, -1, 1), 9000311LL);  // K0*(700)0
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(c12, +1, 1), 9000321LL);  // K0*(700)+
  }

  TEST(PhaseShifts, PiKClebschGordanDecays) {
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiK_I32_Channel();
    // Iz=3/2 -> pi+ K+ (single channel)
    auto d3 = PhaseShifts::ChannelDecays(ch, 3);
    ASSERT_EQ(d3.size(), 1u);
    EXPECT_NEAR(d3[0].first, 1.0, 1e-9);
    EXPECT_EQ(d3[0].second.first, 211); EXPECT_EQ(d3[0].second.second, 321);
    // Iz=1/2 -> 2/3 pi0 K+ (111,321) + 1/3 pi+ K0 (211,311); branchings sum to 1
    auto d1 = PhaseShifts::ChannelDecays(ch, 1);
    ASSERT_EQ(d1.size(), 2u);
    double tot = 0.; for (size_t i = 0; i < d1.size(); ++i) tot += d1[i].first;
    EXPECT_NEAR(tot, 1.0, 1e-9);
    for (size_t i = 0; i < d1.size(); ++i) {
      if (d1[i].second.first == 111 && d1[i].second.second == 321)
        EXPECT_NEAR(d1[i].first, 2.0 / 3.0, 1e-9);   // pi0 K+
      if (d1[i].second.first == 211 && d1[i].second.second == 311)
        EXPECT_NEAR(d1[i].first, 1.0 / 3.0, 1e-9);   // pi+ K0
    }
  }

  TEST(PhaseShifts, PiKStrangeMultipletBuild) {
    // A strange channel builds the FULL multiplet (all Iz, all S=+1) plus the
    // S=-1 antimultiplet - not just Iz>=0 + antiparticles.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_pik_list.dat";
    writePiKList(lst);
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    const int nbefore = TPS.ComponentsNumber();

    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiK_I32_Channel(),
                                      PhaseShifts::PiK_I32_Waves());
    // (4 members + 4 antiparticles) per wave, S+P+D -> 24
    EXPECT_EQ(TPS.ComponentsNumber() - nbefore, 24);

    // Iz=3/2: Q=2, S=+1; Iz=-3/2: Q=-1, S=+1 (distinct member, NOT an antiparticle)
    ASSERT_GE(TPS.PdgToId(99233001), 0);
    EXPECT_EQ(TPS.ParticleByPDG(99233001).ElectricCharge(), 2);
    EXPECT_EQ(TPS.ParticleByPDG(99233001).Strangeness(),   1);
    ASSERT_GE(TPS.PdgToId(99230001), 0);
    EXPECT_EQ(TPS.ParticleByPDG(99230001).ElectricCharge(), -1);
    EXPECT_EQ(TPS.ParticleByPDG(99230001).Strangeness(),    1);
    // antiparticle carries S=-1
    ASSERT_GE(TPS.PdgToId(-99233001), 0);
    EXPECT_EQ(TPS.ParticleByPDG(-99233001).Strangeness(), -1);
    // density attached, repulsive -> negative chi2
    ASSERT_NE(TPS.ParticleByPDG(99233001).GetGeneralizedDensity(), nullptr);
    EXPECT_LT(TPS.ParticleByPDG(99233001).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);
    std::remove(lst.c_str());
  }

  TEST(PhaseShifts, PiKChi2QFromConfig) {
    // Load both pi-K channels from a config and check the net electric-charge
    // susceptibility contribution is negative (repulsive I=3/2, sum Q^2=12,
    // dominates the attractive I=1/2, sum Q^2=2). I=3/2 uses synthetic-cluster
    // files; I=1/2 reuses the real kappa codes ("-", created since absent).
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_pik_cfg_list.dat";
    writePiKList(lst);
    PhaseShifts::WritePhaseShiftFiles(PhaseShifts::PiK_I32_Channel(),
                                      PhaseShifts::PiK_I32_Waves(), dir);
    const std::string confF = dir + "ps_piK.conf";
    {
      std::ofstream f(confF.c_str());
      f << "piK_I32:S  list-piK_I32_S.dat  decays-piK_I32_S.dat  PelaezRodas2016_S\n";
      f << "piK_I12:S  -  -  PelaezRodas2016_S\n";   // kappa reuses 9000321/9000311
    }
    auto suscQ = [&](ThermalParticleSystem& T_) {
      ThermalModelIdeal m(&T_);
      m.SetTemperature(0.150); m.SetBaryonChemicalPotential(0.0);
      m.SetElectricChemicalPotential(0.0); m.SetStrangenessChemicalPotential(0.0);
      m.CalculatePrimordialDensities(); m.CalculateFluctuations();
      return m.Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    };
    ThermalParticleSystem TPSbase(std::vector<std::string>(1, lst));
    double base = suscQ(TPSbase);
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    std::vector<long long> added = PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, confF);
    EXPECT_EQ(added.size(), 12u);                 // (4+4) I=3/2 + (2+2) I=1/2 (kappa)
    EXPECT_GE(TPS.PdgToId(9000321), 0);           // kappa created with its real code
    EXPECT_GE(TPS.PdgToId(9000311), 0);
    double full = suscQ(TPS);
    EXPECT_LT(full - base, 0.0);                  // net repulsive

    std::remove(lst.c_str()); std::remove(confF.c_str());
    std::remove((dir + "list-piK_I32_S.dat").c_str());
    std::remove((dir + "decays-piK_I32_S.dat").c_str());
  }

  // ---- pi-pi I=0 (the sigma; resonant, branch-tracked) ----

  TEST(PhaseShifts, PiPiI0DeltaBranchTracking) {
    // delta_0^0 is attractive (>0), 0 at threshold, rises monotonically THROUGH
    // 90 deg (the sigma) without a discontinuity, and equals d0=226.5 deg exactly
    // at the K-Kbar threshold (the matched intermediate parametrization).
    const double r2d = 180.0 / M_PI;
    const double thr = 2.0 * PhaseShifts::PionMass();
    EXPECT_DOUBLE_EQ(PhaseShifts::PiPi_delta_I0_S(thr), 0.0);
    EXPECT_GT(PhaseShifts::PiPi_delta_I0_S(0.5), 0.0);            // attractive
    // monotonic through 90 deg
    double prev = -1.0; bool crossed90 = false, monotonic = true;
    for (double M = thr + 1e-3; M <= 2.0 * 0.496; M += 0.01) {
      double d = PhaseShifts::PiPi_delta_I0_S(M) * r2d;
      if (d < prev - 1e-6) monotonic = false;
      if (prev < 90.0 && d >= 90.0) crossed90 = true;
      prev = d;
    }
    EXPECT_TRUE(monotonic);
    EXPECT_TRUE(crossed90);                                       // passes the 90 deg point
    // d0 at the K-Kbar threshold (Table V CFD): 226.5 deg
    EXPECT_NEAR(PhaseShifts::PiPi_delta_I0_S(2.0 * 0.496) * r2d, 226.5, 1e-3);
    // continuity (no jump) at the low/intermediate matching point sqrt(sM)=0.85:
    // delta rises steeply here (~3.5 rad/GeV), so use a tiny step - a real jump
    // would be O(0.01+ rad) regardless of step, continuity gives ~slope*2eps.
    EXPECT_NEAR(PhaseShifts::PiPi_delta_I0_S(0.85 - 1e-5),
                PhaseShifts::PiPi_delta_I0_S(0.85 + 1e-5), 1e-3);
  }

  TEST(PhaseShifts, PiPiI0DecaysAndMember) {
    // I=0 isoscalar: single neutral self-conjugate member, decays 2/3 pi+pi- + 1/3 pi0pi0.
    // It reuses the real sigma/f0(500) code (9000221).
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I0_Channel();
    EXPECT_EQ(ch.twoI, 0);
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 0, 1), 9000221LL);
    auto d = PhaseShifts::ChannelDecays(ch, 0);
    ASSERT_EQ(d.size(), 2u);
    double tot = 0.; for (size_t i = 0; i < d.size(); ++i) tot += d[i].first;
    EXPECT_NEAR(tot, 1.0, 1e-9);
    for (size_t i = 0; i < d.size(); ++i) {
      if (d[i].second.first == -211 && d[i].second.second == 211)
        EXPECT_NEAR(d[i].first, 2.0 / 3.0, 1e-9);   // pi+ pi-
      if (d[i].second.first == 111 && d[i].second.second == 111)
        EXPECT_NEAR(d[i].first, 1.0 / 3.0, 1e-9);   // pi0 pi0
    }
  }

  TEST(PhaseShifts, PiPiI0AttractiveRaisesPions) {
    // The attractive I=0 cluster has POSITIVE density and decays to pions, so its
    // feeddown INCREASES the total pion density (opposite to the repulsive I=2).
    const double T = 0.150;
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_i0.dat";
    writePionList(pionF);

    auto totalPiPlus = [&](ThermalParticleSystem& TPS) {
      ThermalModelIdeal m(&TPS);
      m.SetTemperature(T); m.SetBaryonChemicalPotential(0.0);
      m.SetElectricChemicalPotential(0.0); m.SetStrangenessChemicalPotential(0.0);
      m.CalculateDensities();
      return m.GetDensity(211, Feeddown::StabilityFlag);
    };
    ThermalParticleSystem TPSbase(std::vector<std::string>(1, pionF));
    double piBase = totalPiPlus(TPSbase);

    ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I0_Channel(),
                                      PhaseShifts::PiPi_I0_Waves());
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 1);     // single neutral member
    // sigma absent here -> created with its real code 9000221 (a real-resonance
    // cluster, so it counts as "overridden" -> the toggle rebuilds)
    ASSERT_GE(TPS.PdgToId(9000221), 0);
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 1);
    // attractive -> positive chi2 for the cluster
    EXPECT_GT(TPS.ParticleByPDG(9000221).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, T, 0.0), 0.0);
    double piFull = totalPiPlus(TPS);
    EXPECT_GT(piBase, 0.0);
    EXPECT_GT(piFull, piBase);   // attraction -> positive feeddown -> more pions

    std::remove(pionF.c_str());
  }

  TEST(PhaseShifts, PiPiI0f0980ReusesRealResonance) {
    // The f0(980) channel reuses the real f0(980) PDG code (subsumption by
    // coincidence): it does NOT add a particle, overrides the existing f0(980)'s
    // thermal contribution with the phase-shift density, and keeps it (with its
    // decays) in the list.
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I0_f0980_Channel();
    EXPECT_EQ(ch.memberPdg[0], 9010221LL);
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 0, 1), 9010221LL);   // reuses real code
    EXPECT_TRUE(ch.subsumedPdg.empty());                            // no removal

    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_f0_list.dat";
    {
      std::ofstream f(lst.c_str());
      f << "211 pi+      1 0.13957  1 -1 0 1 0 0 0 0 0 0\n";
      f << "111 pi0      1 0.134977 1 -1 0 0 0 0 0 0 0 0\n";
      f << "9010221 f0(980) 0 0.990  1 -1 0 0 0 0 0 0 0 0\n";   // a plain BW resonance
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    const int nbefore = TPS.ComponentsNumber();
    ASSERT_GE(TPS.PdgToId(9010221), 0);
    EXPECT_EQ(TPS.ParticleByPDG(9010221).GetGeneralizedDensity(), nullptr);  // plain BW

    PhaseShifts::AddPhaseShiftChannel(TPS, ch, PhaseShifts::PiPi_I0_f0980_Waves());
    EXPECT_EQ(TPS.ComponentsNumber(), nbefore);                    // NOT a new particle
    ASSERT_GE(TPS.PdgToId(9010221), 0);                            // still present
    ASSERT_NE(TPS.ParticleByPDG(9010221).GetGeneralizedDensity(), nullptr);  // now overridden
    // it is an overridden real resonance (not a synthetic cluster) -> the cheap
    // toggle is not exact, so the GUI rebuilds instead.
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 1);
    // the f0(980) region of delta_0^0 rises -> positive (attractive) contribution
    EXPECT_GT(TPS.ParticleByPDG(9010221).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);
    std::remove(lst.c_str());
  }

  // ---- pi-pi I=1 P-wave (the rho; resonant, branch-tracked) ----

  TEST(PhaseShifts, PiPiI1RhoDelta) {
    // delta_1^1 = 0 at threshold, exactly 90 deg at the rho mass (the (Mrho^2-s)
    // factor zeroes cot delta there), and rises monotonically through it.
    const double r2d = 180.0 / M_PI;
    EXPECT_DOUBLE_EQ(PhaseShifts::PiPi_delta_I1_P(2.0 * PhaseShifts::PionMass()), 0.0);
    EXPECT_NEAR(PhaseShifts::PiPi_delta_I1_P(0.7736) * r2d, 90.0, 1e-6);   // at Mrho
    EXPECT_GT(PhaseShifts::PiPi_delta_I1_P(0.6), 0.0);                     // attractive
    double prev = -1.0; bool crossed90 = false, monotonic = true;
    for (double M = 2.0 * PhaseShifts::PionMass() + 1e-3; M <= 2.0 * 0.496; M += 0.01) {
      double d = PhaseShifts::PiPi_delta_I1_P(M) * r2d;
      if (d < prev - 1e-6) monotonic = false;
      if (prev < 90.0 && d >= 90.0) crossed90 = true;
      prev = d;
    }
    EXPECT_TRUE(monotonic);
    EXPECT_TRUE(crossed90);
  }

  TEST(PhaseShifts, PiPiI1RhoReuseAndDecays) {
    // P-wave I=1 = rho(770): reuses the real rho codes (rho0=113, rho+=213),
    // self-conjugate multiplet. Decays: rho0 -> pi+ pi- ; rho+ -> pi+ pi0.
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I1_Channel();
    EXPECT_EQ(ch.twoI, 2);
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 0, 3), 113LL);   // rho0
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, +2, 3), 213LL);  // rho+
    auto d0 = PhaseShifts::ChannelDecays(ch, 0);                // rho0 -> pi+ pi- (pi0pi0 forbidden)
    ASSERT_EQ(d0.size(), 1u);
    EXPECT_NEAR(d0[0].first, 1.0, 1e-9);
    EXPECT_EQ(d0[0].second.first, -211); EXPECT_EQ(d0[0].second.second, 211);
    auto d1 = PhaseShifts::ChannelDecays(ch, 2);                // rho+ -> pi0 pi+
    ASSERT_EQ(d1.size(), 1u);
    EXPECT_NEAR(d1[0].first, 1.0, 1e-9);
    EXPECT_EQ(d1[0].second.first, 111); EXPECT_EQ(d1[0].second.second, 211);

    // build: a list with the rho -> the density overrides rho0/rho+/rho-
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_rho_list.dat";
    {
      std::ofstream f(lst.c_str());
      f << "211 pi+  1 0.13957  1 -1 0 1 0 0 0 0 0 0\n";
      f << "111 pi0  1 0.134977 1 -1 0 0 0 0 0 0 0 0\n";
      f << "113 rho0 0 0.77526  3 -1 0 0 0 0 0 0 0.147 0.279\n";
      f << "213 rho+ 0 0.77511  3 -1 0 1 0 0 0 0 0.149 0.279\n";
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    const int nbefore = TPS.ComponentsNumber();
    PhaseShifts::AddPhaseShiftChannel(TPS, ch, PhaseShifts::PiPi_I1_Waves());
    EXPECT_EQ(TPS.ComponentsNumber(), nbefore);                 // rho not re-added
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 3);  // rho0, rho+, rho-
    ASSERT_NE(TPS.ParticleByPDG(113).GetGeneralizedDensity(), nullptr);
    ASSERT_NE(TPS.ParticleByPDG(213).GetGeneralizedDensity(), nullptr);
    ASSERT_NE(TPS.ParticleByPDG(-213).GetGeneralizedDensity(), nullptr);  // rho- too
    // attractive resonance -> positive chi2
    EXPECT_GT(TPS.ParticleByPDG(113).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);
    std::remove(lst.c_str());
  }

  TEST(PhaseShifts, PiPiI1FWaveNonResonant) {
    // I=1 F-wave: small, attractive, non-resonant (no 90 deg crossing below 1.42).
    const double r2d = 180.0 / M_PI;
    EXPECT_DOUBLE_EQ(PhaseShifts::PiPi_delta_I1_F(2.0 * PhaseShifts::PionMass()), 0.0);
    EXPECT_GT(PhaseShifts::PiPi_delta_I1_F(1.0), 0.0);            // attractive
    EXPECT_LT(PhaseShifts::PiPi_delta_I1_F(1.42) * r2d, 30.0);    // small (a few deg)

    // synthetic cluster (no resonance to reuse): 2J+1=7, F field, family pi-pi
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I1_F_Channel();
    EXPECT_TRUE(ch.memberPdg.empty());
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, 0, 7), 99120007LL);  // Iz=0
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, +2, 7), 99122007LL); // Iz=+1

    // build (catalog create): synthetic F clusters with pi-pi I=1 decays
    const std::string dir = ::testing::TempDir();
    const std::string pionF = dir + "ps_pions_fw.dat";
    writePionList(pionF);
    ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
    PhaseShifts::AddPhaseShiftChannel(TPS, ch, PhaseShifts::PiPi_I1_F_Waves());
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 3);     // F0, F+, F- (synthetic)
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 0);    // synthetic, not a reuse
    ASSERT_GE(TPS.PdgToId(99120007), 0);
    EXPECT_GT(TPS.ParticleByPDG(99120007).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);  // attractive
    std::remove(pionF.c_str());
  }

  // ---- audit fixes: constituent fugacities, GCE gammas, canonical, ownership --

  TEST(PhaseShifts, ConstituentFugacityMetadata) {
    // P2a: a cluster must carry the SUMMED constituent quark content, not the
    // single-meson powers GetAbsQ() = 2 - |S| - |C| would give. pi-pi -> gammaq^4;
    // pi-K -> gammaq^3 gammaS^1.
    const std::string dir = ::testing::TempDir();
    {
      const std::string pionF = dir + "ps_meta_pipi.dat";
      writePionList(pionF);
      ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
      PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I2_Channel(),
                                        PhaseShifts::PiPi_I2_Waves());
      const ThermalParticle& cl = TPS.ParticleByPDG(99144005);
      EXPECT_DOUBLE_EQ(cl.AbsoluteQuark(), 4.0);        // 2 pions = 4 light quarks
      EXPECT_DOUBLE_EQ(cl.AbsoluteStrangeness(), 0.0);
      std::remove(pionF.c_str());
    }
    {
      const std::string lst = dir + "ps_meta_pik.dat";
      writePiKList(lst);
      ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
      PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiK_I32_Channel(),
                                        PhaseShifts::PiK_I32_Waves());
      const ThermalParticle& cl = TPS.ParticleByPDG(99233001);
      EXPECT_DOUBLE_EQ(cl.AbsoluteQuark(), 3.0);        // pi (2) + K (1 light)
      EXPECT_DOUBLE_EQ(cl.AbsoluteStrangeness(), 1.0);  // K carries one strange quark
      std::remove(lst.c_str());
    }
  }

  TEST(PhaseShifts, GammaFactorsReachGeneralizedDensity) {
    // P1b: gammaq/gammaS must be applied to a generalized density. With a
    // BOLTZMANN cluster the particle density scales EXACTLY as gammaq^Nq gammaS^Ns
    // (before the fix the gammas were ignored: ratio would be 1).
    const std::string dir = ::testing::TempDir();
    // pi-pi I=2 D-wave: gammaq^4
    {
      const std::string pionF = dir + "ps_g_pipi.dat";
      writePionList(pionF);
      ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
      PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I2_Channel(),
                                        PhaseShifts::PiPi_I2_Waves());
      ThermalParticle& cl = TPS.ParticleByPDG(99144005);
      ASSERT_DOUBLE_EQ(cl.AbsoluteQuark(), 4.0);
      // override with a Boltzmann density so the scaling is exact
      std::vector<PhaseShiftPartialWave> dwave(1,
        PhaseShiftPartialWave(5, &PhaseShifts::PiPi_delta_I2_D));
      cl.SetGeneralizedDensity(new PhaseShiftDensity(dwave, PhaseShifts::PionMass(),
        PhaseShifts::PionMass(), PhaseShifts::PiPiI2_Mmax(), /*Boltzmann*/ 0, 64));

      ThermalModelParameters p1; p1.T = 0.150; p1.gammaq = 1.0;
      ThermalModelParameters pg = p1; pg.gammaq = 1.7;
      double n1 = cl.Density(p1, IdealGasFunctions::ParticleDensity, false, 0.0);
      double ng = cl.Density(pg, IdealGasFunctions::ParticleDensity, false, 0.0);
      ASSERT_LT(n1, 0.0);                                   // repulsive
      EXPECT_NEAR(ng / n1, std::pow(1.7, 4), 1e-9);         // gammaq^4 exactly
      std::remove(pionF.c_str());
    }
    // pi-K I=3/2 S-wave: gammaq^3 gammaS^1
    {
      const std::string lst = dir + "ps_g_pik.dat";
      writePiKList(lst);
      ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
      PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiK_I32_Channel(),
                                        PhaseShifts::PiK_I32_Waves());
      ThermalParticle& cl = TPS.ParticleByPDG(99233001);
      ASSERT_DOUBLE_EQ(cl.AbsoluteQuark(), 3.0);
      ASSERT_DOUBLE_EQ(cl.AbsoluteStrangeness(), 1.0);
      std::vector<PhaseShiftPartialWave> swave(1,
        PhaseShiftPartialWave(1, &PhaseShifts::PiK_delta_I32_S));
      cl.SetGeneralizedDensity(new PhaseShiftDensity(swave, PhaseShifts::PionMass(),
        PhaseShifts::KaonMass(), PhaseShifts::PiK_I32_Mmax(), /*Boltzmann*/ 0, 64));

      ThermalModelParameters p1; p1.T = 0.150; p1.gammaq = 1.0; p1.gammaS = 1.0;
      ThermalModelParameters pg = p1; pg.gammaq = 1.3; pg.gammaS = 1.5;
      double n1 = cl.Density(p1, IdealGasFunctions::ParticleDensity, false, 0.0);
      double ng = cl.Density(pg, IdealGasFunctions::ParticleDensity, false, 0.0);
      ASSERT_LT(n1, 0.0);                                   // repulsive
      EXPECT_NEAR(ng / n1, std::pow(1.3, 3) * std::pow(1.5, 1), 1e-9);
      std::remove(lst.c_str());
    }
  }

  TEST(PhaseShifts, CanonicalDensityClusterUsesBethUhlenbeck) {
    // P1a: DensityCluster() (the canonical-ensemble path) must route through the
    // generalized density, not fall back to the ideal pole particle. The repulsive
    // pi-K I=3/2 cluster therefore has NEGATIVE DensityCluster (the BU value);
    // an ideal pole particle would be positive.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_can_pik.dat";
    writePiKList(lst);
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiK_I32_Channel(),
                                      PhaseShifts::PiK_I32_Waves());
    ThermalParticle& cl = TPS.ParticleByPDG(99233001);     // S-wave, S=1, repulsive
    ThermalModelParameters p; p.T = 0.150;
    double nc1 = cl.DensityCluster(1, p, IdealGasFunctions::ParticleDensity, false, 0.0);
    EXPECT_LT(nc1, 0.0);                                    // BU sign, not ideal pole
    // DensityCluster(1) routes exactly to QuantityCluster(1) (meson sign +1)
    double bu = cl.GetGeneralizedDensity()
                  ->QuantityCluster(1, IdealGasFunctions::ParticleDensity, p.T, 0.0);
    EXPECT_DOUBLE_EQ(nc1, bu);
    // The n-th Boltzmann cluster terms (Boltzmann gas at T/n) sum to the full
    // QUANTUM (here Bose) GCE density: Sum_n QuantityCluster(n) == Quantity().
    // (The n=1 term alone, used by the canonical single-term path, is the
    // Boltzmann approximation - it differs by the small Bose enhancement.)
    double sum = 0.;
    for (int n = 1; n <= 12; ++n)
      sum += cl.GetGeneralizedDensity()
               ->QuantityCluster(n, IdealGasFunctions::ParticleDensity, p.T, 0.0);
    double bose = cl.GetGeneralizedDensity()
                    ->Quantity(IdealGasFunctions::ParticleDensity, p.T, 0.0);
    EXPECT_NEAR(sum, bose, 1e-6 * std::fabs(bose) + 1e-12);
    EXPECT_GT(std::fabs(bose), std::fabs(nc1));   // Bose > Boltzmann n=1 in magnitude

    // Contract guard: the T/n reconstruction is valid only for the fugacity-linear
    // quantities (number/pressure/energy) the canonical ensemble sums. A supported
    // quantity has nonzero higher cluster terms; an unsupported one (a
    // susceptibility) keeps ONLY the n==1 Boltzmann term (higher terms -> 0).
    GeneralizedDensity* gd = cl.GetGeneralizedDensity();
    EXPECT_NE(gd->QuantityCluster(2, IdealGasFunctions::ParticleDensity, p.T, 0.0), 0.0);
    EXPECT_NE(gd->QuantityCluster(2, IdealGasFunctions::EnergyDensity, p.T, 0.0), 0.0);
    EXPECT_EQ(gd->QuantityCluster(2, IdealGasFunctions::chi2, p.T, 0.0), 0.0);  // guarded
    EXPECT_EQ(gd->QuantityCluster(3, IdealGasFunctions::chi2, p.T, 0.0), 0.0);
    // n==1 is the Boltzmann value for ANY quantity (even unsupported ones)
    EXPECT_NE(gd->QuantityCluster(1, IdealGasFunctions::chi2, p.T, 0.0), 0.0);
    std::remove(lst.c_str());
  }

  TEST(PhaseShifts, CanonicalModelKeepsBethUhlenbeckSign) {
    // End-to-end audit scenario: in a strangeness-canonical calculation the strange
    // pi-K cluster is treated canonically (DensityCluster path). Its primordial
    // density must stay negative (repulsive BU), i.e. not revert to an ideal pole.
    std::string dir = std::string(ThermalFIST_INPUT_FOLDER);
    ThermalParticleSystem TPS(dir + "/list/PDG2020/list.dat",
                              dir + "/list/PDG2020/decays.dat");
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiK_I32_Channel(),
                                      PhaseShifts::PiK_I32_Waves());
    ASSERT_GE(TPS.PdgToId(99233001), 0);
    EXPECT_EQ(TPS.ParticleByPDG(99233001).Strangeness(), 1);   // canonical sector

    ThermalModelParameters p;
    p.T = 0.150; p.muB = 0.; p.muS = 0.; p.muQ = 0.; p.muC = 0.;
    p.gammaq = 1.; p.gammaS = 1.;
    p.V = 5000.; p.SVc = 5000.;
    p.B = 0; p.Q = 0; p.S = 0; p.C = 0;

    ThermalModelCanonical m(&TPS, p);
    m.SetStatistics(false);
    m.SetUseWidth(ThermalParticle::ZeroWidth);
    m.ConserveBaryonCharge(false);
    m.ConserveElectricCharge(false);
    m.ConserveStrangeness(true);    // strangeness canonical -> pi-K cluster canonical
    m.ConserveCharm(false);
    m.CalculateDensities();

    int id = TPS.PdgToId(99233001);
    double n = m.Densities()[id];
    EXPECT_TRUE(std::isfinite(n));
    EXPECT_LT(n, 0.0);   // BU repulsive survived the canonical (DensityCluster) path
  }

  TEST(PhaseShifts, GeneralizedDensityOwnershipSurvivesCopyAndRebuild) {
    // P2b: the density is shared_ptr-owned, so it survives value-copies of the
    // particle (the list vector copies), the antiparticle does NOT share it, and
    // reassigning the owning system frees the old densities without dangling.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_own_pik.dat";
    writePiKList(lst);

    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiK_I32_Channel(),
                                      PhaseShifts::PiK_I32_Waves());
    // particle and antiparticle each have their OWN density (not shared)
    GeneralizedDensity* dp = TPS.ParticleByPDG(99233001).GetGeneralizedDensity();
    GeneralizedDensity* da = TPS.ParticleByPDG(-99233001).GetGeneralizedDensity();
    ASSERT_NE(dp, nullptr);
    ASSERT_NE(da, nullptr);
    EXPECT_NE(dp, da);

    // a value-copy of the particle shares the same model and stays valid
    ThermalParticle copy = TPS.ParticleByPDG(99233001);
    ASSERT_EQ(copy.GetGeneralizedDensity(), dp);
    EXPECT_LT(copy.GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);

    // reassigning the whole system (the GUI rebuild path) frees the old densities
    // and rebuilds cleanly - no crash, fresh densities attached
    TPS = ThermalParticleSystem(std::vector<std::string>(1, lst));
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiK_I32_Channel(),
                                      PhaseShifts::PiK_I32_Waves());
    ASSERT_NE(TPS.ParticleByPDG(99233001).GetGeneralizedDensity(), nullptr);
    EXPECT_LT(TPS.ParticleByPDG(99233001).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);
    std::remove(lst.c_str());
  }

  // ---- pi-K elastic set: I=3/2 P & D (non-resonant) + K*(892) (resonant) ----

  TEST(PhaseShifts, PiKI32PDWavesSmall) {
    // The I=3/2 P and D waves are repulsive (delta < 0), zero at threshold, and
    // small (high-l, momentum-suppressed) across the elastic range.
    const double r2d = 180.0 / M_PI;
    const double thr = PhaseShifts::PionMass() + PhaseShifts::KaonMass();
    EXPECT_DOUBLE_EQ(PhaseShifts::PiK_delta_I32_P(thr), 0.0);
    EXPECT_DOUBLE_EQ(PhaseShifts::PiK_delta_I32_D(thr), 0.0);
    EXPECT_LT(PhaseShifts::PiK_delta_I32_P(1.0), 0.0);             // repulsive
    EXPECT_LT(PhaseShifts::PiK_delta_I32_D(1.0), 0.0);
    // both are at most a couple of degrees over the whole elastic range
    for (double M = thr + 0.05; M <= PhaseShifts::PiK_I32_Mmax(); M += 0.1) {
      EXPECT_LT(std::fabs(PhaseShifts::PiK_delta_I32_P(M)) * r2d, 5.0) << "M=" << M;
      EXPECT_LT(std::fabs(PhaseShifts::PiK_delta_I32_D(M)) * r2d, 5.0) << "M=" << M;
    }
  }

  TEST(PhaseShifts, PiKI32HasThreeWaves) {
    // The I=3/2 catalog now carries S (2J+1=1), P (3) and D (5).
    auto w = PhaseShifts::PiK_I32_Waves();
    ASSERT_EQ(w.size(), 3u);
    EXPECT_EQ(w[0].twoJplus1, 1);
    EXPECT_EQ(w[1].twoJplus1, 3);
    EXPECT_EQ(w[2].twoJplus1, 5);
    // the analytic registry resolves each wave by its per-wave model name
    EXPECT_EQ(PhaseShifts::AnalyticWave("piK_I32", "PelaezRodas2016_S").twoJplus1, 1);
    EXPECT_EQ(PhaseShifts::AnalyticWave("piK_I32", "PelaezRodas2016_P").twoJplus1, 3);
    EXPECT_EQ(PhaseShifts::AnalyticWave("piK_I32", "PelaezRodas2016_D").twoJplus1, 5);
  }

  TEST(PhaseShifts, PiKKstar892ResonantBranchTracked) {
    // delta_1^{1/2} = the K*(892): 0 at threshold, branch-tracked through exactly
    // 90 deg at the resonance mass (m_r = 0.8957, where (m_r^2 - s) zeroes cot
    // delta), rising monotonically; elastic only below the K-eta threshold.
    const double r2d = 180.0 / M_PI;
    const double thr = PhaseShifts::PionMass() + PhaseShifts::KaonMass();
    EXPECT_DOUBLE_EQ(PhaseShifts::PiK_delta_I12_P(thr), 0.0);
    EXPECT_NEAR(PhaseShifts::PiK_delta_I12_P(0.8957) * r2d, 90.0, 1e-6);  // at m_r
    EXPECT_GT(PhaseShifts::PiK_delta_I12_P(0.8), 0.0);                    // attractive
    // elastic region only -> 0 at/above the K-eta threshold
    EXPECT_DOUBLE_EQ(PhaseShifts::PiK_delta_I12_P(PhaseShifts::PiK_I12_Mmax()), 0.0);
    EXPECT_DOUBLE_EQ(PhaseShifts::PiK_delta_I12_P(1.2), 0.0);
    // monotonic through 90 deg up to (just below) the K-eta threshold
    double prev = -1.0; bool crossed90 = false, monotonic = true;
    for (double M = thr + 1e-3; M < PhaseShifts::PiK_I12_Mmax() - 1e-3; M += 0.005) {
      double d = PhaseShifts::PiK_delta_I12_P(M) * r2d;
      if (d < prev - 1e-6) monotonic = false;
      if (prev < 90.0 && d >= 90.0) crossed90 = true;
      prev = d;
    }
    EXPECT_TRUE(monotonic);
    EXPECT_TRUE(crossed90);
  }

  TEST(PhaseShifts, PiKKstar892ReusesRealResonance) {
    // The K*(892) channel reuses the real K*(892) codes (subsumption by PDG
    // coincidence): K*(892)+ = 323, K*(892)0 = 313, the P-wave 2J+1=3 carrying the
    // K* spin. Strange multiplet -> all Iz are distinct members.
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiK_K892_Channel();
    EXPECT_EQ(ch.twoI, 1);
    EXPECT_EQ(ch.memberPdg[+1], 323LL);
    EXPECT_EQ(ch.memberPdg[-1], 313LL);
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, +1, 3), 323LL);   // K*(892)+
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, -1, 3), 313LL);   // K*(892)0
    EXPECT_TRUE(ch.subsumedPdg.empty());                          // override, no removal

    // build a list with the K* present -> the density overrides K*+/K*0 and their
    // antiparticles (K*- = -323, K*bar0 = -313); they stay as decay products.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_kstar_list.dat";
    {
      std::ofstream f(lst.c_str());
      f << "211 pi+      1 0.13957  1 -1 0 1 0 0 0 0 0 0\n";
      f << "111 pi0      1 0.134977 1 -1 0 0 0 0 0 0 0 0\n";
      f << "321 K+       1 0.493677 1 -1 0 1 1 0 1 0 0 0\n";
      f << "311 K0       1 0.497611 1 -1 0 0 1 0 1 0 0 0\n";
      f << "323 K*(892)+ 0 0.89167  3 -1 0 1 1 0 1 0 0.0514 0.633\n";
      f << "313 K*(892)0 0 0.89555  3 -1 0 0 1 0 1 0 0.0473 0.633\n";
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    const int nbefore = TPS.ComponentsNumber();
    EXPECT_EQ(TPS.ParticleByPDG(323).GetGeneralizedDensity(), nullptr);  // plain pole/BW
    PhaseShifts::AddPhaseShiftChannel(TPS, ch, PhaseShifts::PiK_K892_Waves());
    EXPECT_EQ(TPS.ComponentsNumber(), nbefore);                  // K* not re-added
    // K*+, K*0 and their two antiparticles -> 4 overridden real resonances
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 4);
    ASSERT_NE(TPS.ParticleByPDG(323).GetGeneralizedDensity(), nullptr);
    ASSERT_NE(TPS.ParticleByPDG(313).GetGeneralizedDensity(), nullptr);
    ASSERT_NE(TPS.ParticleByPDG(-323).GetGeneralizedDensity(), nullptr);
    ASSERT_NE(TPS.ParticleByPDG(-313).GetGeneralizedDensity(), nullptr);
    // attractive resonance -> positive chi2
    EXPECT_GT(TPS.ParticleByPDG(323).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);
    std::remove(lst.c_str());
  }

  TEST(PhaseShifts, PiKElasticSetFromConfig) {
    // Load the full elastic pi-K set from a config: I=3/2 S/P/D (synthetic files),
    // I=1/2 S = kappa (reuse 9000321/9000311) and the K*(892) P-wave (reuse
    // 323/313). Check the cluster counts and that the densities are attached.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_pik_elastic_list.dat";
    {
      std::ofstream f(lst.c_str());
      f << "211 pi+      1 0.13957  1 -1 0 1 0 0 0 0 0 0\n";
      f << "111 pi0      1 0.134977 1 -1 0 0 0 0 0 0 0 0\n";
      f << "321 K+       1 0.493677 1 -1 0 1 1 0 1 0 0 0\n";
      f << "311 K0       1 0.497611 1 -1 0 0 1 0 1 0 0 0\n";
      f << "323 K*(892)+ 0 0.89167  3 -1 0 1 1 0 1 0 0.0514 0.633\n";
      f << "313 K*(892)0 0 0.89555  3 -1 0 0 1 0 1 0 0.0473 0.633\n";
    }
    PhaseShifts::WritePhaseShiftFiles(PhaseShifts::PiK_I32_Channel(),
                                      PhaseShifts::PiK_I32_Waves(), dir);
    const std::string confF = dir + "ps_piK_elastic.conf";
    {
      std::ofstream f(confF.c_str());
      f << "piK_I32:S   list-piK_I32_S.dat  decays-piK_I32_S.dat  PelaezRodas2016_S\n";
      f << "piK_I32:P   list-piK_I32_P.dat  decays-piK_I32_P.dat  PelaezRodas2016_P\n";
      f << "piK_I32:D   list-piK_I32_D.dat  decays-piK_I32_D.dat  PelaezRodas2016_D\n";
      f << "piK_I12:S   -  -  PelaezRodas2016_S\n";   // kappa reuses 9000321/9000311
      f << "piK_K892:P  -  -  PelaezRodas2016_P\n";   // K*(892) reuses 323/313
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    std::vector<long long> added = PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, confF);
    // I=3/2: 8 per wave x 3 waves = 24; kappa created (4); K*(892) overrides (4).
    EXPECT_EQ(added.size(), 32u);
    // K*(892) overridden (real codes), kappa created, I=3/2 synthetic clusters.
    EXPECT_NE(TPS.ParticleByPDG(323).GetGeneralizedDensity(), nullptr);   // K* overridden
    EXPECT_NE(TPS.ParticleByPDG(313).GetGeneralizedDensity(), nullptr);
    EXPECT_GE(TPS.PdgToId(9000321), 0);                                   // kappa created
    EXPECT_GE(TPS.PdgToId(9000311), 0);
    EXPECT_GE(TPS.PdgToId(99233003), 0);                                  // I=3/2 P cluster
    EXPECT_GE(TPS.PdgToId(99233005), 0);                                  // I=3/2 D cluster
    // Both the K*(892) (323/313 + antis) and the kappa (9000321/9000311 + antis)
    // sit on real (non-synthetic) codes, so both count as overridden resonances:
    // 4 (K*) + 4 (kappa) = 8. (Created vs pre-existing does not matter - it is the
    // real code that makes the cheap toggle inexact, so the GUI rebuilds.)
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 8);

    std::remove(lst.c_str()); std::remove(confF.c_str());
    std::remove((dir + "list-piK_I32_S.dat").c_str());
    std::remove((dir + "list-piK_I32_P.dat").c_str());
    std::remove((dir + "list-piK_I32_D.dat").c_str());
    std::remove((dir + "decays-piK_I32_S.dat").c_str());
    std::remove((dir + "decays-piK_I32_P.dat").c_str());
    std::remove((dir + "decays-piK_I32_D.dat").c_str());
  }

  // ---- pi-N: the Delta(1232) (meson-baryon, B=1, Fermi, resonant) ----

  TEST(PhaseShifts, PiNDeltaP33Resonant) {
    // delta_1+^{3/2} = the Delta(1232): 0 at threshold, branch-tracked through
    // 90 deg at the resonance (~1.232 GeV), rising to ~150 deg by the matching
    // point 1.38 GeV; elastic across the resonance.
    const double r2d = 180.0 / M_PI;
    const double thr = PhaseShifts::PionMass() + PhaseShifts::NucleonMass();
    EXPECT_DOUBLE_EQ(PhaseShifts::PiN_delta_P33(thr), 0.0);
    EXPECT_GT(PhaseShifts::PiN_delta_P33(1.15), 0.0);                  // attractive
    EXPECT_NEAR(PhaseShifts::PiN_delta_P33(1.232) * r2d, 90.0, 5.0);  // ~90 at the Delta
    EXPECT_GT(PhaseShifts::PiN_delta_P33(1.38) * r2d, 140.0);         // well past 90
    // monotonic through 90 deg up to the matching point
    double prev = -1.0; bool crossed90 = false, monotonic = true;
    for (double M = thr + 1e-3; M <= PhaseShifts::PiN_Delta_Mmax(); M += 0.005) {
      double d = PhaseShifts::PiN_delta_P33(M) * r2d;
      if (d < prev - 1e-6) monotonic = false;
      if (prev < 90.0 && d >= 90.0) crossed90 = true;
      prev = d;
    }
    EXPECT_TRUE(monotonic);
    EXPECT_TRUE(crossed90);
  }

  TEST(PhaseShifts, PiNDeltaChannelAndDecays) {
    // Baryon channel: I=3/2, B=+1, fermionic (stat=+1), reuses the real Delta
    // codes; not self-conjugate so every Iz is a distinct member.
    PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiN_Delta_Channel();
    EXPECT_EQ(ch.twoI, 3);
    EXPECT_EQ(ch.B, 1);
    EXPECT_EQ(ch.statistics, 1);                  // Fermi
    EXPECT_EQ(ch.memberPdg[+3], 2224LL);          // Delta++
    EXPECT_EQ(ch.memberPdg[+1], 2214LL);          // Delta+
    EXPECT_EQ(ch.memberPdg[-1], 2114LL);          // Delta0
    EXPECT_EQ(ch.memberPdg[-3], 1114LL);          // Delta-
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(ch, +3, 4), 2224LL);   // reuses real code
    EXPECT_TRUE(ch.subsumedPdg.empty());
    // electric charge Q = Iz + (B+S+C)/2: Delta++ Q=2, Delta- Q=-1
    EXPECT_EQ((3 + ch.B) / 2, 2);
    EXPECT_EQ((-3 + ch.B) / 2, -1);
    // isospin-CG decays: Delta++ -> pi+ p; Delta+ -> 2/3 pi0 p + 1/3 pi+ n
    auto d3 = PhaseShifts::ChannelDecays(ch, 3);
    ASSERT_EQ(d3.size(), 1u);
    EXPECT_NEAR(d3[0].first, 1.0, 1e-9);
    EXPECT_EQ(d3[0].second.first, 211); EXPECT_EQ(d3[0].second.second, 2212);
    auto d1 = PhaseShifts::ChannelDecays(ch, 1);
    ASSERT_EQ(d1.size(), 2u);
    double tot = 0.; for (size_t i = 0; i < d1.size(); ++i) tot += d1[i].first;
    EXPECT_NEAR(tot, 1.0, 1e-9);
    for (size_t i = 0; i < d1.size(); ++i) {
      if (d1[i].second.first == 111 && d1[i].second.second == 2212)
        EXPECT_NEAR(d1[i].first, 2.0 / 3.0, 1e-9);   // pi0 p
      if (d1[i].second.first == 211 && d1[i].second.second == 2112)
        EXPECT_NEAR(d1[i].first, 1.0 / 3.0, 1e-9);   // pi+ n
    }
  }

  TEST(PhaseShifts, PiNDeltaBuildOverride) {
    // Build against a list with the Delta present: the density overrides all four
    // Delta charge states AND their four antiparticles (8 overridden resonances),
    // the cluster is a B=1 fermion with pi+N constituent fugacity (gammaq^5), and
    // the attractive resonance gives a positive chi2.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_piN_list.dat";
    writePiNList(lst);
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    const int nbefore = TPS.ComponentsNumber();
    EXPECT_EQ(TPS.ParticleByPDG(2224).GetGeneralizedDensity(), nullptr);  // plain resonance

    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiN_Delta_Channel(),
                                      PhaseShifts::PiN_Delta_Waves());
    EXPECT_EQ(TPS.ComponentsNumber(), nbefore);          // nothing added (all overridden)
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 8);  // 4 Delta + 4 anti-Delta
    for (long long pdg : {2224LL, 2214LL, 2114LL, 1114LL, -2224LL, -2214LL, -2114LL, -1114LL})
      ASSERT_NE(TPS.ParticleByPDG(pdg).GetGeneralizedDensity(), nullptr) << "pdg=" << pdg;

    const ThermalParticle& dpp = TPS.ParticleByPDG(2224);
    EXPECT_EQ(dpp.BaryonCharge(), 1);
    EXPECT_EQ(dpp.Statistics(), 1);                      // Fermi cluster
    EXPECT_DOUBLE_EQ(dpp.AbsoluteQuark(), 5.0);          // pi (2) + N (3) light quarks
    EXPECT_DOUBLE_EQ(dpp.AbsoluteStrangeness(), 0.0);
    // attractive resonance -> positive chi2 (Beth-Uhlenbeck)
    EXPECT_GT(dpp.GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);
    std::remove(lst.c_str());
  }

  TEST(PhaseShifts, PiNDeltaRaisesNucleons) {
    // Against a Delta-less pi+N base, the attractive Delta cluster (created here,
    // its codes being absent) decays to pi N, so its POSITIVE Beth-Uhlenbeck
    // feeddown increases the proton density. (Note: vs a base that already has a
    // zero-width pole Delta the BU effective density is smaller - the pole-mass
    // treatment overcounts the broad resonance; that is the point of the S-matrix
    // term - so the right additive test uses a Delta-less base.)
    const double T = 0.150;
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_piN_fd.dat";   // pi + N only, NO Delta
    {
      std::ofstream f(lst.c_str());
      f << "211  pi+ 1 0.13957  1 -1 0  1 0 0 0 0 0 0\n";
      f << "111  pi0 1 0.134977 1 -1 0  0 0 0 0 0 0 0\n";
      f << "2212 p   1 0.938272 2  1 1  1 0 0 0 0 0 0\n";
      f << "2112 n   1 0.939565 2  1 1  0 0 0 0 0 0 0\n";
    }
    auto protons = [&](ThermalParticleSystem& TPS) {
      ThermalModelIdeal m(&TPS);
      m.SetTemperature(T); m.SetBaryonChemicalPotential(0.0);
      m.SetElectricChemicalPotential(0.0); m.SetStrangenessChemicalPotential(0.0);
      m.CalculateDensities();
      return m.GetDensity(2212, Feeddown::StabilityFlag);
    };
    ThermalParticleSystem TPSbase(std::vector<std::string>(1, lst));
    double pBase = protons(TPSbase);
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiN_Delta_Channel(),
                                      PhaseShifts::PiN_Delta_Waves());
    ASSERT_GE(TPS.PdgToId(2224), 0);          // Delta created (was absent)
    EXPECT_GT(PhaseShifts::CountOverriddenResonances(TPS), 0);  // real codes -> rebuild toggle
    double pFull = protons(TPS);
    EXPECT_GT(pBase, 0.0);
    EXPECT_GT(pFull, pBase);   // attraction -> positive feeddown -> more protons
    std::remove(lst.c_str());
  }

  TEST(PhaseShifts, PiNDeltaFromConfig) {
    // Config path: the Delta reuses its real codes ("-" files); the wave key is
    // the numeric 2J+1 = 4 (J=3/2; the S/P/D letters encode 2J+1 = 2l+1 for
    // mesons, which does not apply to baryon waves).
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_piN_cfg.dat";
    writePiNList(lst);
    const std::string confF = dir + "ps_piN.conf";
    {
      std::ofstream f(confF.c_str());
      f << "piN_Delta:4   -   -   RoySteiner2016_P33\n";   // Delta, reuses 2224/2214/2114/1114
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    std::vector<long long> added = PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, confF);
    EXPECT_EQ(added.size(), 8u);                          // 4 Delta + 4 anti-Delta
    EXPECT_NE(TPS.ParticleByPDG(2224).GetGeneralizedDensity(), nullptr);
    EXPECT_NE(TPS.ParticleByPDG(-2224).GetGeneralizedDensity(), nullptr);
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 8);
    std::remove(lst.c_str()); std::remove(confF.c_str());
  }

  // ---- pi-N non-resonant Roy-Steiner background waves (S31, S11, P31, P11, P13) ----

  TEST(PhaseShifts, PiNBackgroundDeltaValues) {
    // Lock the transcription against the verified Roy-Steiner values (degrees) at
    // sqrt(s) = 1.15 and 1.20 GeV (see piN_RoySteiner_reference.md). Signs: S31
    // repulsive, S11 attractive, P-waves small.
    const double r2d = 180.0 / M_PI;
    EXPECT_NEAR(PhaseShifts::PiN_delta_S31(1.15) * r2d, -7.3, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_S31(1.20) * r2d, -11.2, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_S11(1.15) * r2d,  8.4, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_S11(1.20) * r2d, 10.1, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_P13(1.15) * r2d, -1.0, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_P13(1.20) * r2d, -2.0, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_P31(1.15) * r2d, -1.8, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_P31(1.20) * r2d, -3.6, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_P11(1.15) * r2d, -1.0, 0.15);
    EXPECT_NEAR(PhaseShifts::PiN_delta_P11(1.20) * r2d,  0.3, 0.15);
    // zero at threshold; non-resonant (no 90 deg crossing in the elastic range)
    const double thr = PhaseShifts::PionMass() + PhaseShifts::NucleonMass();
    EXPECT_DOUBLE_EQ(PhaseShifts::PiN_delta_S31(thr), 0.0);
    EXPECT_LT(PhaseShifts::PiN_delta_S31(1.30) * r2d, 0.0);
    EXPECT_GT(PhaseShifts::PiN_delta_S11(1.30) * r2d, 0.0);
    EXPECT_LT(std::fabs(PhaseShifts::PiN_delta_P13(1.30) * r2d), 30.0);  // small
  }

  TEST(PhaseShifts, PiNBackgroundPdgDisambiguation) {
    // Baryon waves with the same (I, J) but different orbital l (S31/P31 and
    // S11/P11, all 2J+1=2) must get DISTINCT synthetic codes via the excitation
    // (n n) slot = orbital l. Meson codes (excitation 0) are unaffected.
    PhaseShifts::PhaseShiftChannel s31 = PhaseShifts::PiN_S31_Channel();
    PhaseShifts::PhaseShiftChannel p31 = PhaseShifts::PiN_P31_Channel();
    EXPECT_EQ(s31.excitation, 0);   // S-wave (l=0)
    EXPECT_EQ(p31.excitation, 1);   // P-wave (l=1)
    long long cs = PhaseShifts::PhaseShiftPdgId(s31, 3, 2);   // S31 Iz=+3/2
    long long cp = PhaseShifts::PhaseShiftPdgId(p31, 3, 2);   // P31 Iz=+3/2
    EXPECT_NE(cs, cp);
    EXPECT_EQ((cs / 10) % 100, 0);  // n n = 0 for S
    EXPECT_EQ((cp / 10) % 100, 1);  // n n = 1 for P
    EXPECT_EQ(cs % 10, 2);          // 2J+1 = 2 (last digit) for both
    EXPECT_EQ(cp % 10, 2);
    // family/2I unchanged
    EXPECT_EQ((cs / 100000) % 10, 4);   // FamilyPiN
    EXPECT_EQ((cs / 10000) % 10, 3);    // 2I = 3
    // S11 vs P11 (I=1/2) likewise distinct
    EXPECT_NE(PhaseShifts::PhaseShiftPdgId(PhaseShifts::PiN_S11_Channel(), 1, 2),
              PhaseShifts::PhaseShiftPdgId(PhaseShifts::PiN_P11_Channel(), 1, 2));
    // regression: a meson channel (excitation 0) is unchanged (pi-pi I=2 D-wave)
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(PhaseShifts::PiPi_I2_Channel(), 4, 5), 99144005LL);
  }

  TEST(PhaseShifts, PiNBackgroundSyntheticBuild) {
    // The background waves are synthetic clusters (no resonance to reuse): build
    // creates the full I=3/2 multiplet + antiparticles (8), none overriding a real
    // code; B=1 Fermi clusters with pi+N fugacity (gammaq^5); S31 repulsive -> chi2<0,
    // S11 attractive -> chi2>0.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_piN_bg.dat";
    {
      std::ofstream f(lst.c_str());        // pi + N only (no resonances)
      f << "211  pi+ 1 0.13957  1 -1 0  1 0 0 0 0 0 0\n";
      f << "111  pi0 1 0.134977 1 -1 0  0 0 0 0 0 0 0\n";
      f << "2212 p   1 0.938272 2  1 1  1 0 0 0 0 0 0\n";
      f << "2112 n   1 0.939565 2  1 1  0 0 0 0 0 0 0\n";
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    const int nbefore = TPS.ComponentsNumber();
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiN_S31_Channel(),
                                      PhaseShifts::PiN_S31_Waves());
    EXPECT_EQ(TPS.ComponentsNumber() - nbefore, 8);            // 4 members + 4 anti (synthetic)
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 8);
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 0); // synthetic, no real-code reuse
    const long long s31pp = PhaseShifts::PhaseShiftPdgId(PhaseShifts::PiN_S31_Channel(), 3, 2);
    ASSERT_GE(TPS.PdgToId(s31pp), 0);
    const ThermalParticle& cl = TPS.ParticleByPDG(s31pp);
    EXPECT_EQ(cl.BaryonCharge(), 1);
    EXPECT_EQ(cl.ElectricCharge(), 2);          // Iz=+3/2, B=1 -> Q=2
    EXPECT_EQ(cl.Statistics(), 1);              // Fermi
    EXPECT_DOUBLE_EQ(cl.AbsoluteQuark(), 5.0);  // pi (2) + N (3)
    EXPECT_LT(cl.GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);  // repulsive

    // S11 (I=1/2, attractive): 2 members + 2 anti = 4; chi2 > 0
    ThermalParticleSystem TPS2(std::vector<std::string>(1, lst));
    PhaseShifts::AddPhaseShiftChannel(TPS2, PhaseShifts::PiN_S11_Channel(),
                                      PhaseShifts::PiN_S11_Waves());
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS2), 4);
    const long long s11p = PhaseShifts::PhaseShiftPdgId(PhaseShifts::PiN_S11_Channel(), 1, 2);
    ASSERT_GE(TPS2.PdgToId(s11p), 0);
    EXPECT_GT(TPS2.ParticleByPDG(s11p).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);   // attractive
    std::remove(lst.c_str());
  }

  TEST(PhaseShifts, PiNFullSetFromConfig) {
    // Load the full pi-N set (Delta + 5 backgrounds) from a config against a list
    // with the Delta present: the Delta overrides its 8 real codes; the 5
    // background waves create synthetic clusters - I=3/2 S31+P31 = 2x8 = 16, I=1/2
    // S11+P11+P13 = 3x4 = 12. Total added = 8 + 16 + 12 = 36.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_piN_full.dat";
    writePiNList(lst);
    const std::string confF = dir + "ps_piN_full.conf";
    {
      std::ofstream f(confF.c_str());
      f << "piN_Delta:4  -  -  RoySteiner2016_P33\n";
      f << "piN_S31:2    -  -  RoySteiner2016_S31\n";
      f << "piN_P13:4    -  -  RoySteiner2016_P13\n";
      f << "piN_S11:2    -  -  RoySteiner2016_S11\n";
      f << "piN_P31:2    -  -  RoySteiner2016_P31\n";
      f << "piN_P11:2    -  -  RoySteiner2016_P11\n";
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    std::vector<long long> added = PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, confF);
    EXPECT_EQ(added.size(), 36u);                  // 8 (Delta) + 16 (I=3/2 x2) + 12 (I=1/2 x3)
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 8);   // only the Delta reuses real codes
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 36);
    // a thermal calculation with the full pi-N set runs
    ThermalModelIdeal m(&TPS); m.ClearDensityModels();
    m.SetTemperature(0.150); m.SetBaryonChemicalPotential(0.0);
    m.SetElectricChemicalPotential(0.0); m.SetStrangenessChemicalPotential(0.0);
    m.CalculateDensities();
    EXPECT_TRUE(std::isfinite(m.GetDensity(2212, Feeddown::StabilityFlag)));
    std::remove(lst.c_str()); std::remove(confF.c_str());
  }

  // ---- K-N (kaon-nucleon, exotic S=+1; Gibbs-Arceo) ----

  TEST(PhaseShifts, KNDeltaValues) {
    // Lock the transcription against the verified Gibbs-Arceo values (degrees) at
    // sqrt(s) = 1.494 (q~0.2) and 1.565 (q~0.3); see KN_GibbsArceo_reference.md.
    // I=1 S-wave strongly repulsive; I=0 P-waves spin-orbit split (P01 +, P03 -).
    const double r2d = 180.0 / M_PI;
    EXPECT_NEAR(PhaseShifts::KN_delta_S11(1.494) * r2d, -18.1, 0.3);
    EXPECT_NEAR(PhaseShifts::KN_delta_S11(1.565) * r2d, -27.5, 0.3);
    EXPECT_NEAR(PhaseShifts::KN_delta_S01(1.565) * r2d, -10.0, 0.3);
    EXPECT_NEAR(PhaseShifts::KN_delta_P01(1.565) * r2d,  25.7, 0.4);   // I=0 P1/2 attractive
    EXPECT_NEAR(PhaseShifts::KN_delta_P03(1.565) * r2d,  -4.4, 0.3);   // I=0 P3/2 repulsive
    // zero at threshold; signs
    const double thr = PhaseShifts::KaonMass() + PhaseShifts::NucleonMass();
    EXPECT_DOUBLE_EQ(PhaseShifts::KN_delta_S11(thr), 0.0);
    EXPECT_LT(PhaseShifts::KN_delta_S11(1.50), 0.0);    // I=1 S-wave repulsive
    EXPECT_LT(PhaseShifts::KN_delta_S01(1.50), 0.0);    // I=0 S-wave weakly repulsive
    EXPECT_GT(PhaseShifts::KN_delta_P01(1.50), 0.0);    // spin-orbit: P01 > 0
    EXPECT_LT(PhaseShifts::KN_delta_P03(1.50), 0.0);    //              P03 < 0
  }

  TEST(PhaseShifts, KNChannelStructure) {
    // Exotic S=+1 baryon channels: B=+1, S=+1, fermionic; not self-conjugate.
    PhaseShifts::PhaseShiftChannel s11 = PhaseShifts::KN_S11_Channel();
    EXPECT_EQ(s11.twoI, 2);            // I=1
    EXPECT_EQ(s11.B, 1);
    EXPECT_EQ(s11.S, 1);
    EXPECT_EQ(s11.statistics, 1);      // Fermi
    EXPECT_TRUE(s11.memberPdg.empty()); // synthetic (no S=+1 resonance to reuse)
    EXPECT_EQ(PhaseShifts::KN_S01_Channel().twoI, 0);   // I=0
    // S11 vs P11 (same I=1, J=1/2, 2J+1=2) -> distinct codes via excitation = l
    long long cS = PhaseShifts::PhaseShiftPdgId(s11, 2, 2);                       // Iz=+1
    long long cP = PhaseShifts::PhaseShiftPdgId(PhaseShifts::KN_P11_Channel(), 2, 2);
    EXPECT_NE(cS, cP);
    EXPECT_EQ((cS / 100000) % 10, 6);  // FamilyKN
    EXPECT_EQ((cS / 10) % 100, 0);     // S-wave n n = 0
    EXPECT_EQ((cP / 10) % 100, 1);     // P-wave n n = 1
  }

  TEST(PhaseShifts, KNSyntheticBuild) {
    // Build the I=1 S-wave against a K+N list: synthetic clusters (3 members +
    // 3 antibaryons), B=1 S=1 Fermi with gammaq^4 gammaS, repulsive -> chi2 < 0.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_KN.dat";
    {
      std::ofstream f(lst.c_str());     // K + N only (no S=+1 resonances exist)
      f << "321  K+ 1 0.493677 1 -1 0  1 1 0 1 0 0 0\n";
      f << "311  K0 1 0.497611 1 -1 0  0 1 0 1 0 0 0\n";
      f << "2212 p  1 0.938272 2  1 1  1 0 0 0 0 0 0\n";
      f << "2112 n  1 0.939565 2  1 1  0 0 0 0 0 0 0\n";
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    const int nbefore = TPS.ComponentsNumber();
    PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::KN_S11_Channel(),
                                      PhaseShifts::KN_S11_Waves());
    EXPECT_EQ(TPS.ComponentsNumber() - nbefore, 6);            // 3 members + 3 anti
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 6);
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 0); // synthetic
    const long long kpp = PhaseShifts::PhaseShiftPdgId(PhaseShifts::KN_S11_Channel(), 2, 2);
    ASSERT_GE(TPS.PdgToId(kpp), 0);
    const ThermalParticle& cl = TPS.ParticleByPDG(kpp);
    EXPECT_EQ(cl.BaryonCharge(), 1);
    EXPECT_EQ(cl.Strangeness(), 1);            // exotic S=+1
    EXPECT_EQ(cl.ElectricCharge(), 2);         // Iz=+1, (2Iz+B+S)/2 = 2 (K+ p)
    EXPECT_EQ(cl.Statistics(), 1);             // Fermi
    EXPECT_DOUBLE_EQ(cl.AbsoluteQuark(), 4.0); // K (1 light) + N (3 light)
    EXPECT_DOUBLE_EQ(cl.AbsoluteStrangeness(), 1.0);
    EXPECT_LT(cl.GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);  // repulsive
    std::remove(lst.c_str());
  }

  TEST(PhaseShifts, KNFromConfig) {
    // Load the full K-N set (6 waves) from a config: all synthetic (no overrides).
    // I=1 (S11/P11/P13): 3x6 = 18; I=0 (S01/P01/P03): 3x2 = 6. Total 24.
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_KN_cfg.dat";
    {
      std::ofstream f(lst.c_str());
      f << "321  K+ 1 0.493677 1 -1 0  1 1 0 1 0 0 0\n";
      f << "311  K0 1 0.497611 1 -1 0  0 1 0 1 0 0 0\n";
      f << "2212 p  1 0.938272 2  1 1  1 0 0 0 0 0 0\n";
      f << "2112 n  1 0.939565 2  1 1  0 0 0 0 0 0 0\n";
    }
    const std::string confF = dir + "ps_KN.conf";
    {
      std::ofstream f(confF.c_str());
      f << "KN_S11:2 - - GibbsArceo2007_S11\n";
      f << "KN_P11:2 - - GibbsArceo2007_P11\n";
      f << "KN_P13:4 - - GibbsArceo2007_P13\n";
      f << "KN_S01:2 - - GibbsArceo2007_S01\n";
      f << "KN_P01:2 - - GibbsArceo2007_P01\n";
      f << "KN_P03:4 - - GibbsArceo2007_P03\n";
    }
    ThermalParticleSystem TPS(std::vector<std::string>(1, lst));
    std::vector<long long> added = PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, confF);
    EXPECT_EQ(added.size(), 24u);                          // 18 (I=1 x3) + 6 (I=0 x3)
    EXPECT_EQ(PhaseShifts::CountPhaseShiftDensities(TPS), 24);
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS), 0);  // all synthetic
    ThermalModelIdeal m(&TPS); m.ClearDensityModels();
    m.SetTemperature(0.150); m.SetBaryonChemicalPotential(0.0);
    m.SetElectricChemicalPotential(0.0); m.SetStrangenessChemicalPotential(0.0);
    m.CalculateDensities(); m.CalculateFluctuations();
    // exotic S=+1 baryons are present and the calculation runs
    EXPECT_TRUE(std::isfinite(m.Pressure()));
    std::remove(lst.c_str()); std::remove(confF.c_str());
  }

  // ---- subsumed-resonance per-resonance files (create-if-absent / keep-if-present) ----

  TEST(PhaseShifts, SkipDuplicatesLoaderFlag) {
    // The skip_duplicates load flag must KEEP an already-present code (with its QN)
    // instead of throwing on the duplicate, while still adding genuinely new codes.
    const std::string dir = ::testing::TempDir();
    const std::string base = dir + "ps_dup_base.dat";
    {
      std::ofstream f(base.c_str());
      f << "113 rho0 0 0.77526 3 -1 0 0 0 0 0 0 0.149 0.279\n";
    }
    const std::string more = dir + "ps_dup_more.dat";
    {
      std::ofstream f(more.c_str());
      f << "113 rho0_alt 0 0.990  1 -1 0 0 0 0 0 0 0.5    0.279\n";  // duplicate pdg, other QN
      f << "223 omega    0 0.7827 3 -1 0 0 0 0 0 0 0.0085 0.41\n";   // genuinely new
    }
    // without the flag: a duplicate pdg throws
    ThermalParticleSystem TPS1(std::vector<std::string>(1, base), std::vector<std::string>());
    EXPECT_THROW(TPS1.AddParticlesToListFromFile(more), std::invalid_argument);
    // with the flag: existing 113 kept (deg=3, not 1), omega added
    ThermalParticleSystem TPS2(std::vector<std::string>(1, base), std::vector<std::string>());
    std::set<std::string> flags;
    flags.insert(ThermalParticleSystem::flag_skip_duplicates);
    EXPECT_NO_THROW(TPS2.AddParticlesToListFromFile(more, flags));
    TPS2.FinalizeList();
    ASSERT_GE(TPS2.PdgToId(113), 0);
    EXPECT_DOUBLE_EQ(TPS2.ParticleByPDG(113).Degeneracy(), 3.0);  // kept the original, not 1
    EXPECT_NEAR(TPS2.ParticleByPDG(113).Mass(), 0.77526, 1e-9);   // kept the original mass
    EXPECT_GE(TPS2.PdgToId(223), 0);                               // new code added
    std::remove(base.c_str()); std::remove(more.c_str());
  }

  TEST(PhaseShifts, SubsumedResonanceFilesCarryRealQN) {
    // The generated per-resonance list file writes the REAL resonance QN (mass,
    // deg=2J+1, width), not the cluster nominal (m1+m2, deg 1). Check the rho file.
    const std::string dir = ::testing::TempDir();
    PhaseShifts::WritePhaseShiftFiles(PhaseShifts::PiPi_I1_Channel(),
                                      PhaseShifts::PiPi_I1_Waves(), dir);
    const std::string base = dir + "ps_sub_pions.dat";
    writePionList(base);
    // Load the rho list file directly (no densities) and check the QN on 113.
    std::vector<std::string> lists; lists.push_back(base);
    lists.push_back(dir + "list-pipi_I1_P.dat");
    ThermalParticleSystem TPS(lists);
    ASSERT_GE(TPS.PdgToId(113), 0);
    EXPECT_DOUBLE_EQ(TPS.ParticleByPDG(113).Degeneracy(), 3.0);    // 2J+1, not 1
    EXPECT_NEAR(TPS.ParticleByPDG(113).Mass(), 0.77526, 1e-6);     // real rho mass, not 2*mpi
    EXPECT_GE(TPS.PdgToId(213), 0);
    EXPECT_GE(TPS.PdgToId(-213), 0);                               // rho- auto-generated
    std::remove(base.c_str());
    std::remove((dir + "list-pipi_I1_P.dat").c_str());
    std::remove((dir + "decays-pipi_I1_P.dat").c_str());
  }

  TEST(PhaseShifts, SubsumedResonanceFromFileCreateOrKeep) {
    // The headline feature: a subsumed-resonance channel referenced by per-resonance
    // files (a) CREATES the resonance with its real QN + CG decays where the base
    // list lacks it, and (b) KEEPS the existing entry (QN + its own decays) where the
    // list has it, attaching the density either way (no duplicates).
    const std::string dir = ::testing::TempDir();
    PhaseShifts::WritePhaseShiftFiles(PhaseShifts::PiPi_I1_Channel(),
                                      PhaseShifts::PiPi_I1_Waves(), dir);
    const std::string confF = dir + "ps_rho_file.conf";
    {
      std::ofstream f(confF.c_str());
      f << "pipi_I1:P  list-pipi_I1_P.dat  decays-pipi_I1_P.dat  GarciaMartin2011_P\n";
    }

    // (a) rho ABSENT -> created with real QN + the isospin-CG decay
    const std::string pionF = dir + "ps_rho_pions.dat";
    writePionList(pionF);
    ThermalParticleSystem TPS(std::vector<std::string>(1, pionF));
    ASSERT_LT(TPS.PdgToId(113), 0);                               // absent before
    PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, confF);
    ASSERT_GE(TPS.PdgToId(113), 0);                               // created
    EXPECT_DOUBLE_EQ(TPS.ParticleByPDG(113).Degeneracy(), 3.0);   // real deg, not 1
    EXPECT_NEAR(TPS.ParticleByPDG(113).Mass(), 0.77526, 1e-6);    // real mass
    ASSERT_NE(TPS.ParticleByPDG(113).GetGeneralizedDensity(), nullptr);   // density attached
    EXPECT_GT(TPS.ParticleByPDG(113).GetGeneralizedDensity()
                ->Quantity(IdealGasFunctions::chi2, 0.150, 0.0), 0.0);    // attractive rho
    ASSERT_EQ(TPS.ParticleByPDG(113).Decays().size(), 1u);        // CG rho0 -> pi+ pi-

    // (b) rho PRESENT with a distinctive mass and TWO of its own decays -> kept,
    //     density attached, not duplicated, decays NOT overwritten by the CG single.
    const std::string rhoL = dir + "ps_rho_with.dat";
    const std::string rhoD = dir + "ps_rho_with_dec.dat";
    {
      std::ofstream f(rhoL.c_str());
      f << "211 pi+  1 0.13957  1 -1 0 1 0 0 0 0 0 0\n";
      f << "111 pi0  1 0.134977 1 -1 0 0 0 0 0 0 0 0\n";
      f << "113 rho0 0 0.77000  3 -1 0 0 0 0 0 0 0.150 0.279\n";  // distinctive mass 0.770
      f << "213 rho+ 0 0.77000  3 -1 0 1 0 0 0 0 0.150 0.279\n";
    }
    {   // give rho0 two (placeholder) decay channels to prove they survive
        // (leading '#' header selects the new free decay format)
      std::ofstream f(rhoD.c_str());
      f << "# rho0 placeholder decays\n113\n2\n0.5 -211 211\n0.5 111 111\n";
    }
    ThermalParticleSystem TPS2(std::vector<std::string>(1, rhoL),
                               std::vector<std::string>(1, rhoD));
    const int nbefore = TPS2.ComponentsNumber();
    ASSERT_EQ(TPS2.ParticleByPDG(113).Decays().size(), 2u);       // its own decays
    PhaseShifts::AddPhaseShiftChannelsFromFile(TPS2, confF);
    EXPECT_EQ(TPS2.ComponentsNumber(), nbefore);                  // no duplicate added
    EXPECT_NEAR(TPS2.ParticleByPDG(113).Mass(), 0.77000, 1e-6);   // KEPT the list's QN
    EXPECT_DOUBLE_EQ(TPS2.ParticleByPDG(113).Degeneracy(), 3.0);
    EXPECT_EQ(TPS2.ParticleByPDG(113).Decays().size(), 2u);       // KEPT its own decays (not CG)
    ASSERT_NE(TPS2.ParticleByPDG(113).GetGeneralizedDensity(), nullptr);  // density attached
    EXPECT_EQ(PhaseShifts::CountOverriddenResonances(TPS2), 3);   // rho0, rho+, rho-

    std::remove(pionF.c_str()); std::remove(confF.c_str());
    std::remove(rhoL.c_str()); std::remove(rhoD.c_str());
    std::remove((dir + "list-pipi_I1_P.dat").c_str());
    std::remove((dir + "decays-pipi_I1_P.dat").c_str());
  }

}
