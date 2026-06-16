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
#include "HRGPhaseShifts/PhaseShiftDensity.h"
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"
#include "HRGPhaseShifts/PhaseShiftModel.h"
#include "HRGPhaseShifts/PhaseShiftChannel.h"
#include "HRGBase/IdealGasFunctions.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "gtest/gtest.h"

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
    PhaseShifts::PhaseShiftChannel c12 = PhaseShifts::PiK_I12_Channel();
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(c12, -1, 1), 99210001LL);
    EXPECT_EQ(PhaseShifts::PhaseShiftPdgId(c12, +1, 1), 99211001LL);
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
    EXPECT_EQ(TPS.ComponentsNumber() - nbefore, 8);   // 4 members + 4 antiparticles

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
    // Load both pi-K channels from a config (per-wave files) and check the net
    // electric-charge susceptibility contribution is negative (repulsive I=3/2,
    // sum Q^2=12, dominates the attractive I=1/2, sum Q^2=2).
    const std::string dir = ::testing::TempDir();
    const std::string lst = dir + "ps_pik_cfg_list.dat";
    writePiKList(lst);
    PhaseShifts::WritePhaseShiftFiles(PhaseShifts::PiK_I32_Channel(),
                                      PhaseShifts::PiK_I32_Waves(), dir);
    PhaseShifts::WritePhaseShiftFiles(PhaseShifts::PiK_I12_Channel(),
                                      PhaseShifts::PiK_I12_Waves(), dir);
    const std::string confF = dir + "ps_piK.conf";
    {
      std::ofstream f(confF.c_str());
      f << "piK_I32:S  list-piK_I32_S.dat  decays-piK_I32_S.dat  PelaezRodas2016_S\n";
      f << "piK_I12:S  list-piK_I12_S.dat  decays-piK_I12_S.dat  PelaezRodas2016_S\n";
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
    EXPECT_EQ(added.size(), 12u);                 // (4+4) I=3/2 + (2+2) I=1/2
    double full = suscQ(TPS);
    EXPECT_LT(full - base, 0.0);                  // net repulsive

    std::remove(lst.c_str()); std::remove(confF.c_str());
    std::remove((dir + "list-piK_I32_S.dat").c_str());
    std::remove((dir + "decays-piK_I32_S.dat").c_str());
    std::remove((dir + "list-piK_I12_S.dat").c_str());
    std::remove((dir + "decays-piK_I12_S.dat").c_str());
  }

}
