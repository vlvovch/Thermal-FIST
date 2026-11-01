/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

/**
 * Unit tests for the HRG in partial chemical equilibrium (PCE):
 * ThermalModelPCE, ThermalModelPCEAnnihilation, and the code paths that use them.
 *
 * Tests:
 *  - Regression: the PCE chemical potentials survive ThermalModelBase::SetParameters()
 *    (after a2ec20a every PCE solve silently returned full chemical equilibrium, mu = 0)
 *  - Reference PCE values for T_ch = 155 MeV -> T = 90 MeV (PDG2014 list)
 *  - Conservation laws for several PCE configurations: effective numbers of stable hadrons,
 *    entropy, net baryon number, electric charge, strangeness, final yields after strong decays
 *  - Law of mass action for resonances, Saha equation for light nuclei,
 *    frozen long-lived resonances
 *  - Return to T_ch, fixed-volume mode as the inverse of the fixed-temperature mode
 *  - Baryon annihilation: conservation laws, mu_N + mu_Nbar = <npi>/3 sum_pi mu_pi,
 *    no annihilating species reproduces ThermalModelPCE
 *  - The PCE state survives the follow-up calculations (feeddown, fluctuations,
 *    temperature derivatives) and is used by the event generator and the thermal fit
 *  - The ThermalModelBase chemical potential behavior that PCE relies on
 */

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelParameters.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGPCE/ThermalModelPCE.h"
#include "HRGPCE/ThermalModelPCEAnnihilation.h"
#include "HRGEventGenerator/SphericalBlastWaveEventGenerator.h"
#include "HRGFit/ThermalModelFit.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  // Chemical freeze-out in central Pb-Pb collisions at the LHC, as in arXiv:1903.10024
  const double kTch = 0.155;  // GeV
  const double kVch = 4700.;  // fm^3

  // Temperature in the hadronic phase
  const double kT = 0.090;    // GeV

  // Relative accuracy of the PCE solution (the Broyden solver converges to 1e-10)
  const double kAccuracy = 1.e-7;

  std::string ListFile(const std::string& name) {
    return std::string(ThermalFIST_INPUT_FOLDER) + "/list/" + name;
  }

  // Chemical freeze-out state of the HRG model; at finite muB, muQ and muS are fixed by Q/B = 0.4 and S = 0
  void SetChemicalFreezeoutState(ThermalModelBase& model, double muB = 0.) {
    model.SetParameters(ThermalModelParameters(kTch, muB, 0., 0., 1., kVch));
    if (muB != 0.) {
      model.SetQoverB(0.4);
      model.ConstrainMuQ(true);
      model.ConstrainMuS(true);
      model.ConstrainChemicalPotentials();
    }
    model.FillChemicalPotentials();
  }

  double Mu(ThermalModelBase& model, long long pdg) {
    return model.ChemicalPotential(model.TPS()->PdgToId(pdg));
  }

  double Entropy(ThermalModelBase& model) {
    return model.EntropyDensity() * model.Volume();
  }

  // Index of the PCE stable component of a given species (-1 if the species is not PCE-stable)
  int StableIndex(const ThermalModelPCE& pce, long long pdg) {
    const int id = pce.ThermalModel()->TPS()->PdgToId(pdg);
    for (size_t is = 0; is < pce.StableMapTo().size(); ++is)
      if (pce.StableMapTo()[is] == id)
        return static_cast<int>(is);
    return -1;
  }

  // Effective numbers of the PCE-stable hadrons in the current state of the HRG model,
  // \bar{N}_s = V \sum_i q_i^s n_i, where q_i^s is the mean number of hadrons s from species i
  std::vector<double> EffectiveNumbers(const ThermalModelPCE& pce) {
    const ThermalModelBase* model = pce.ThermalModel();
    const std::vector< std::vector<double> >& q = pce.EffectiveCharges();
    std::vector<double> ret(pce.StableMapTo().size(), 0.);
    for (size_t i = 0; i < q.size(); ++i)
      for (size_t is = 0; is < ret.size(); ++is)
        ret[is] += q[i][is] * model->Densities()[i] * model->Volume();
    return ret;
  }

  // Net conserved charges and the total absolute charge content (used as the scale)
  struct ConservedCharges {
    double B, Q, S;
    double absB, absQ, absS;
  };

  ConservedCharges CalculateConservedCharges(ThermalModelBase& model) {
    ConservedCharges ret = { 0., 0., 0., 0., 0., 0. };
    for (int i = 0; i < model.TPS()->ComponentsNumber(); ++i) {
      const ThermalParticle& part = model.TPS()->Particles()[i];
      const double N = model.Densities()[i] * model.Volume();
      ret.B += part.BaryonCharge() * N;
      ret.Q += part.ElectricCharge() * N;
      ret.S += part.Strangeness() * N;
      ret.absB += std::abs(part.BaryonCharge()) * N;
      ret.absQ += std::abs(part.ElectricCharge()) * N;
      ret.absS += std::abs(part.Strangeness()) * N;
    }
    return ret;
  }

  struct PCEConfiguration {
    std::string name;
    std::string list;
    double muB;
    bool sahaForNuclei;
    bool freezeLonglived;
    bool normalizeBranchingRatios;
  };

  // Net charges are conserved exactly only if the decays carry the full charges of the resonances.
  // In the PDG2014 decay table the branching ratios of K4*(2045) add up to 0.408 only,
  // so at finite chemical potentials it is used with normalized branching ratios.
  std::vector<PCEConfiguration> PCEConfigurations() {
    std::vector<PCEConfiguration> ret;
    ret.push_back({ "PDG2014",                                   ListFile("PDG2014/list.dat"),            0.,  true,  false, false });
    ret.push_back({ "PDG2014, muB = 200 MeV, normalized BRs",    ListFile("PDG2014/list.dat"),            0.2, true,  false, true  });
    ret.push_back({ "PDG2014, frozen long-lived resonances",     ListFile("PDG2014/list.dat"),            0.,  true,  true,  false });
    ret.push_back({ "PDG2014 with nuclei, Saha equation",        ListFile("PDG2014/list-withnuclei.dat"), 0.,  true,  false, false });
    ret.push_back({ "PDG2014 with nuclei, frozen nuclei",        ListFile("PDG2014/list-withnuclei.dat"), 0.,  false, false, false });
    ret.push_back({ "Default list",                              ThermalFIST_DEFAULT_LIST_FILE,           0.,  true,  false, false });
    ret.push_back({ "Default list, muB = 200 MeV",               ThermalFIST_DEFAULT_LIST_FILE,           0.2, true,  false, false });
    return ret;
  }

} // namespace

// Regression test. Since a2ec20a ThermalModelBase::SetParameters() refills the chemical potentials
// from (muB, muQ, muS, muC). The PCE equations set the PCE chemical potentials before calling it,
// so CalculatePCE() silently returned full chemical equilibrium: mu = 0 and V/V_ch = (T_ch/T)^3.
TEST(PCE, RegressionChemicalPotentialsSurviveSetParameters) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  SetChemicalFreezeoutState(model);

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pce.CalculatePCE(kT);

  EXPECT_GT(Mu(model, 211), 0.05) << "PCE returned full chemical equilibrium";
  EXPECT_GT(std::abs(model.Volume() / kVch / std::pow(kTch / kT, 3) - 1.), 0.03)
    << "PCE volume equals the initial guess (T_ch/T)^3 V_ch";

  // The HRG model is left in the PCE state
  ASSERT_EQ(model.ChemicalPotentials().size(), pce.ChemicalPotentials().size());
  for (size_t i = 0; i < pce.ChemicalPotentials().size(); ++i)
    EXPECT_DOUBLE_EQ(model.ChemicalPotentials()[i], pce.ChemicalPotentials()[i]) << TPS.Particles()[i].Name();
  EXPECT_DOUBLE_EQ(model.Volume(), pce.Volume());
  EXPECT_DOUBLE_EQ(model.Parameters().T, kT);
}

// Reference values for T_ch = 155 MeV -> T = 90 MeV with the PDG2014 list used in arXiv:1903.10024
TEST(PCE, ReferenceValuesPDG2014) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  SetChemicalFreezeoutState(model);

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());

  // The chemical potentials grow and the volume increases as the system cools down
  double mupi = 0., V = kVch;
  const double temperatures[] = { 0.150, 0.130, 0.110 };
  for (double T : temperatures) {
    pce.CalculatePCE(T);
    EXPECT_GT(Mu(model, 211), mupi) << "T = " << T;
    EXPECT_GT(model.Volume(), V) << "T = " << T;
    mupi = Mu(model, 211);
    V = model.Volume();
  }

  pce.CalculatePCE(kT);
  EXPECT_NEAR(model.Volume() / kVch, 4.767595532, 1.e-6);
  EXPECT_NEAR(Mu(model, 211), 0.084618476, 1.e-6);
  EXPECT_NEAR(Mu(model, 321), 0.189264723, 1.e-6);
  EXPECT_NEAR(Mu(model, 2212), 0.376928025, 1.e-6);
  EXPECT_NEAR(Mu(model, 3122), 0.435258871, 1.e-6);
  EXPECT_NEAR(Mu(model, 3334), 0.640271769, 1.e-6);
}

// Effective numbers of all PCE-stable hadrons, entropy, and the conserved charges stay at their T_ch values
TEST(PCE, ConservationLaws) {
  const std::vector<PCEConfiguration> configurations = PCEConfigurations();
  for (const PCEConfiguration& config : configurations) {
    SCOPED_TRACE(config.name);

    ThermalParticleSystem TPS(config.list);
    ThermalModelIdeal model(&TPS);
    model.SetNormBratio(config.normalizeBranchingRatios);
    SetChemicalFreezeoutState(model, config.muB);

    ThermalModelPCE pce(&model);
    pce.UseSahaForNuclei(config.sahaForNuclei);
    pce.FreezeLonglivedResonances(config.freezeLonglived);
    pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());

    const std::vector<double> Nch = EffectiveNumbers(pce);
    const double Sch = Entropy(model);
    const ConservedCharges Cch = CalculateConservedCharges(model);

    pce.CalculatePCE(kT);

    EXPECT_GT(Mu(model, 211), 0.05) << "PCE returned full chemical equilibrium";

    const std::vector<double> N = EffectiveNumbers(pce);
    ASSERT_EQ(N.size(), Nch.size());
    for (size_t is = 0; is < N.size(); ++is)
      EXPECT_NEAR(N[is], Nch[is], kAccuracy * Nch[is]) << TPS.Particles()[pce.StableMapTo()[is]].Name();

    EXPECT_NEAR(Entropy(model) / Sch, 1., kAccuracy);

    const ConservedCharges C = CalculateConservedCharges(model);
    EXPECT_NEAR(C.B, Cch.B, kAccuracy * Cch.absB);
    EXPECT_NEAR(C.Q, Cch.Q, kAccuracy * Cch.absQ);
    EXPECT_NEAR(C.S, Cch.S, kAccuracy * Cch.absS);
  }
}

// Independent of the PCE effective charges: the final yields after strong decays are frozen
TEST(PCE, FinalYieldsAfterStrongDecaysAreFrozen) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  SetChemicalFreezeoutState(model);
  model.CalculateDensities();

  const long long pdgs[] = { 211, -211, 111, 321, -321, 2212, -2212, 3122, -3122, 3312, -3312, 3334, -3334, 221 };
  std::vector<double> Ych;
  for (long long pdg : pdgs)
    Ych.push_back(model.GetYield(pdg, Feeddown::Strong));

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pce.CalculatePCE(kT);

  for (size_t i = 0; i < Ych.size(); ++i)
    EXPECT_NEAR(model.GetYield(pdgs[i], Feeddown::Strong) / Ych[i], 1., kAccuracy) << pdgs[i];
}

// Law of mass action: the chemical potential of a resonance is the average one of its decay products
TEST(PCE, ResonanceChemicalPotentials) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  SetChemicalFreezeoutState(model);

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pce.CalculatePCE(kT);

  EXPECT_GT(Mu(model, 211), 0.05) << "PCE returned full chemical equilibrium";
  // rho0 -> pi+ pi-
  EXPECT_NEAR(Mu(model, 113), Mu(model, 211) + Mu(model, -211), 1.e-10);
  // Delta++ -> p pi+
  EXPECT_NEAR(Mu(model, 2224), Mu(model, 2212) + Mu(model, 211), 1.e-10);
  // K*(892)+ -> K+ pi0 (1/3), K0 pi+ (2/3)
  EXPECT_NEAR(Mu(model, 323), 0.3333 * (Mu(model, 321) + Mu(model, 111)) + 0.6667 * (Mu(model, 311) + Mu(model, 211)), 1.e-10);
}

// Light nuclei: either the Saha equation, mu_A = Z mu_p + (A - Z) mu_n (+ mu_Lambda), or frozen yields
TEST(PCE, LightNuclei) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list-withnuclei.dat"));
  const long long pdgd = 1000010020;
  ASSERT_NE(TPS.PdgToId(pdgd), -1);

  const bool sahaFlags[] = { true, false };
  for (bool saha : sahaFlags) {
    SCOPED_TRACE(saha ? "Saha equation" : "Frozen nuclei");

    ThermalModelIdeal model(&TPS);
    SetChemicalFreezeoutState(model);
    model.CalculateDensities();
    const double dpch = model.GetDensity(pdgd, Feeddown::StabilityFlag) / model.GetDensity(2212, Feeddown::StabilityFlag);

    ThermalModelPCE pce(&model);
    pce.UseSahaForNuclei(saha);
    pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
    pce.CalculatePCE(kT);

    const double dp = model.GetDensity(pdgd, Feeddown::StabilityFlag) / model.GetDensity(2212, Feeddown::StabilityFlag);

    if (saha) {
      EXPECT_EQ(StableIndex(pce, pdgd), -1);
      const double mup = Mu(model, 2212), mun = Mu(model, 2112), muL = Mu(model, 3122);
      EXPECT_NEAR(Mu(model, pdgd), mup + mun, 1.e-10);                   // d
      EXPECT_NEAR(Mu(model, 1000010030), mup + 2. * mun, 1.e-10);        // t
      EXPECT_NEAR(Mu(model, 1000020030), 2. * mup + mun, 1.e-10);        // 3He
      EXPECT_NEAR(Mu(model, 1000020040), 2. * mup + 2. * mun, 1.e-10);   // 4He
      EXPECT_NEAR(Mu(model, 1010010030), mup + mun + muL, 1.e-10);       // hypertriton
      EXPECT_NEAR(Mu(model, -pdgd), Mu(model, -2212) + Mu(model, -2112), 1.e-10);

      // The d/p ratio stays close to its chemical freeze-out value (arXiv:1903.10024),
      // whereas full chemical equilibrium would reduce it by more than an order of magnitude
      EXPECT_NEAR(dp / dpch, 1.387, 0.01);
    }
    else {
      EXPECT_NE(StableIndex(pce, pdgd), -1);
      EXPECT_NEAR(dp / dpch, 1., kAccuracy);
    }
  }
}

// Stability flags: hadrons not decaying strongly are frozen, optionally also long-lived resonances
TEST(PCE, StabilityFlags) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list-withnuclei.dat"));

  std::vector<int> flags = ThermalModelPCE::ComputePCEStabilityFlags(&TPS, true, false, 0.015);
  EXPECT_EQ(flags[TPS.PdgToId(211)], 1);
  EXPECT_EQ(flags[TPS.PdgToId(2212)], 1);
  EXPECT_EQ(flags[TPS.PdgToId(3334)], 1);
  EXPECT_EQ(flags[TPS.PdgToId(311)], 1);         // (anti)K0 instead of K0S and K0L
  EXPECT_EQ(flags[TPS.PdgToId(310)], 0);
  EXPECT_EQ(flags[TPS.PdgToId(113)], 0);
  EXPECT_EQ(flags[TPS.PdgToId(223)], 0);
  EXPECT_EQ(flags[TPS.PdgToId(1000010020)], 0);  // Saha equation for nuclei

  flags = ThermalModelPCE::ComputePCEStabilityFlags(&TPS, false, true, 0.015);
  EXPECT_EQ(flags[TPS.PdgToId(223)], 1);         // omega(782), width 8.5 MeV
  EXPECT_EQ(flags[TPS.PdgToId(333)], 1);         // phi(1020), width 4.3 MeV
  EXPECT_EQ(flags[TPS.PdgToId(113)], 0);         // rho(770), width 149 MeV
  EXPECT_EQ(flags[TPS.PdgToId(2224)], 0);        // Delta(1232), width 117 MeV
  EXPECT_EQ(flags[TPS.PdgToId(1000010020)], 1);  // frozen nuclei
}

// Cooling down and heating back up to T_ch returns the chemical freeze-out state
TEST(PCE, ReturnToChemicalFreezeout) {
  ThermalParticleSystem TPS(ThermalFIST_DEFAULT_LIST_FILE);
  ThermalModelIdeal model(&TPS);
  SetChemicalFreezeoutState(model, 0.2);
  const std::vector<double> much = model.ChemicalPotentials();

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pce.CalculatePCE(0.100);
  pce.CalculatePCE(kTch);

  EXPECT_NEAR(model.Volume() / kVch, 1., kAccuracy);
  for (size_t i = 0; i < much.size(); ++i)
    EXPECT_NEAR(model.ChemicalPotentials()[i], much[i], 1.e-8) << TPS.Particles()[i].Name();
}

// PCE at fixed volume is the inverse of PCE at fixed temperature
TEST(PCE, FixedVolumeMode) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS), modelV(&TPS);
  SetChemicalFreezeoutState(model);
  SetChemicalFreezeoutState(modelV);

  ThermalModelPCE pce(&model), pceV(&modelV);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pceV.SetChemicalFreezeout(modelV.Parameters(), modelV.ChemicalPotentials());

  pce.CalculatePCE(0.100);
  pceV.CalculatePCE(model.Volume(), ThermalModelPCE::AtFixedVolume);

  EXPECT_NEAR(modelV.Parameters().T, 0.100, 1.e-8);
  EXPECT_NEAR(modelV.Volume() / model.Volume(), 1., 1.e-12);
  for (size_t i = 0; i < model.ChemicalPotentials().size(); ++i)
    EXPECT_NEAR(modelV.ChemicalPotentials()[i], model.ChemicalPotentials()[i], 1.e-8) << TPS.Particles()[i].Name();
}

TEST(PCE, CalculatePCEWithoutChemicalFreezeoutThrows) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  ThermalModelPCE pce(&model);
  EXPECT_THROW(pce.CalculatePCE(kT), std::invalid_argument);
}

// The calculation sequence of the GUI (QtThermalFIST model tab) after a PCE solve keeps the PCE state
TEST(PCE, StateSurvivesFollowUpCalculations) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  SetChemicalFreezeoutState(model);
  model.CalculateDensities();

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pce.CalculatePCE(kT);
  EXPECT_GT(Mu(model, 211), 0.05) << "PCE returned full chemical equilibrium";

  const std::vector<double> mu = model.ChemicalPotentials();
  const std::vector<double> dens = model.Densities();
  const double V = model.Volume();

  model.CalculateFeeddown();
  model.CalculateFluctuations();
  model.CalculateTemperatureDerivatives();

  EXPECT_DOUBLE_EQ(model.Parameters().T, kT);
  EXPECT_DOUBLE_EQ(model.Volume(), V);
  for (size_t i = 0; i < mu.size(); ++i) {
    EXPECT_DOUBLE_EQ(model.ChemicalPotentials()[i], mu[i]) << TPS.Particles()[i].Name();
    EXPECT_NEAR(model.Densities()[i], dens[i], 1.e-12 * dens[i]) << TPS.Particles()[i].Name();
  }
}

// Baryon annihilation N Nbar <-> <npi> pi in PCE (arXiv:2210.15641)
TEST(PCEAnnihilation, ConservationLaws) {
  struct AnnihilationCase {
    double muB;
    double npi;
  };
  const AnnihilationCase cases[] = { { 0., 5. }, { 0.2, 6. } };

  for (const AnnihilationCase& c : cases) {
    SCOPED_TRACE("muB = " + std::to_string(c.muB) + ", <npi> = " + std::to_string(c.npi));

    ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
    ThermalModelIdeal model(&TPS);
    SetChemicalFreezeoutState(model, c.muB);

    ThermalModelPCEAnnihilation pce(&model);
    pce.SetPionAnnihilationNumber(c.npi);
    pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());

    const std::vector<double> Nch = EffectiveNumbers(pce);
    const double Sch = Entropy(model);

    pce.CalculatePCE(kT);

    EXPECT_GT(Mu(model, 211), 0.05) << "PCE returned full chemical equilibrium";

    const std::vector<double> N = EffectiveNumbers(pce);
    ASSERT_EQ(N.size(), Nch.size());

    const int ip = StableIndex(pce, 2212), iap = StableIndex(pce, -2212);
    const int in = StableIndex(pce, 2112), ian = StableIndex(pce, -2112);
    const int ipip = StableIndex(pce, 211), ipim = StableIndex(pce, -211), ipiz = StableIndex(pce, 111);
    ASSERT_TRUE(ip >= 0 && iap >= 0 && in >= 0 && ian >= 0 && ipip >= 0 && ipim >= 0 && ipiz >= 0);

    // Hadrons not taking part in annihilation
    for (size_t is = 0; is < N.size(); ++is) {
      const int ind = static_cast<int>(is);
      if (ind == ip || ind == iap || ind == in || ind == ian || ind == ipip || ind == ipim || ind == ipiz)
        continue;
      EXPECT_NEAR(N[is], Nch[is], kAccuracy * Nch[is]) << TPS.Particles()[pce.StableMapTo()[is]].Name();
    }

    // Net nucleon numbers and the pion isospin composition
    EXPECT_NEAR(N[ip] - N[iap], Nch[ip] - Nch[iap], kAccuracy * (Nch[ip] + Nch[iap]));
    EXPECT_NEAR(N[in] - N[ian], Nch[in] - Nch[ian], kAccuracy * (Nch[in] + Nch[ian]));
    EXPECT_NEAR(N[ipip] - N[ipim], Nch[ipip] - Nch[ipim], kAccuracy * (Nch[ipip] + Nch[ipim]));
    EXPECT_NEAR(N[ipip] - N[ipiz], Nch[ipip] - Nch[ipiz], kAccuracy * (Nch[ipip] + Nch[ipiz]));

    // (N + Nbar) + 2/<npi> (pi+ + pi0 + pi-) is conserved by the annihilation reactions
    const double combch = Nch[ip] + Nch[iap] + Nch[in] + Nch[ian] + 2. / c.npi * (Nch[ipip] + Nch[ipiz] + Nch[ipim]);
    const double comb = N[ip] + N[iap] + N[in] + N[ian] + 2. / c.npi * (N[ipip] + N[ipiz] + N[ipim]);
    EXPECT_NEAR(comb / combch, 1., kAccuracy);

    EXPECT_NEAR(Entropy(model) / Sch, 1., kAccuracy);

    // Chemical equilibrium of the annihilation reactions
    const double sumpi = Mu(model, 211) + Mu(model, 111) + Mu(model, -211);
    EXPECT_NEAR(Mu(model, 2212) + Mu(model, -2212), c.npi / 3. * sumpi, 1.e-10);
    EXPECT_NEAR(Mu(model, 2112) + Mu(model, -2112), c.npi / 3. * sumpi, 1.e-10);

    // Annihilation reduces the (anti)nucleon yields
    EXPECT_LT(N[ip] + N[iap], Nch[ip] + Nch[iap]);
  }
}

// Without annihilating species ThermalModelPCEAnnihilation reduces to ThermalModelPCE
TEST(PCEAnnihilation, WithoutAnnihilationReproducesPCE) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal modelA(&TPS), model(&TPS);
  SetChemicalFreezeoutState(modelA);
  SetChemicalFreezeoutState(model);

  ThermalModelPCEAnnihilation pceA(&modelA);
  pceA.SetStabilityFlags(pceA.RecalculateStabilityFlags(std::vector<long long>()));
  pceA.SetChemicalFreezeout(modelA.Parameters(), modelA.ChemicalPotentials());
  pceA.CalculatePCE(kT);

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pce.CalculatePCE(kT);

  EXPECT_GT(Mu(modelA, 211), 0.05) << "PCE returned full chemical equilibrium";
  EXPECT_NEAR(modelA.Volume() / model.Volume(), 1., 1.e-9);
  for (size_t i = 0; i < model.ChemicalPotentials().size(); ++i)
    EXPECT_NEAR(modelA.ChemicalPotentials()[i], model.ChemicalPotentials()[i], 1.e-9) << TPS.Particles()[i].Name();
}

// The event generator in PCE mode, configured as in the QtThermalFIST event generator tab
TEST(PCEEventGenerator, UsesPCEChemicalPotentials) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  SetChemicalFreezeoutState(model);

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pce.CalculatePCE(kT);

  const EventGeneratorConfiguration::Ensemble ensembles[] = { EventGeneratorConfiguration::GCE, EventGeneratorConfiguration::CE };
  for (EventGeneratorConfiguration::Ensemble ensemble : ensembles) {
    SCOPED_TRACE(ensemble == EventGeneratorConfiguration::GCE ? "GCE" : "CE");

    EventGeneratorConfiguration config;
    config.fEnsemble = ensemble;
    config.fModelType = EventGeneratorConfiguration::PointParticle;
    config.CFOParameters = model.Parameters();
    config.fUsePCE = true;
    config.fPCEChems = pce.ChemicalPotentials();

    SphericalBlastWaveEventGenerator generator(&TPS, config, kT, 0.5);
    ThermalModelBase* genmodel = generator.ThermalModel();

    EXPECT_GT(Mu(*genmodel, 211), 0.05) << "PCE returned full chemical equilibrium";
    ASSERT_EQ(genmodel->ChemicalPotentials().size(), model.ChemicalPotentials().size());
    EXPECT_DOUBLE_EQ(genmodel->Volume(), model.Volume());
    for (size_t i = 0; i < model.ChemicalPotentials().size(); ++i) {
      EXPECT_DOUBLE_EQ(genmodel->ChemicalPotentials()[i], model.ChemicalPotentials()[i]) << TPS.Particles()[i].Name();
      EXPECT_NEAR(genmodel->Densities()[i], model.Densities()[i], 1.e-10 * model.Densities()[i]) << TPS.Particles()[i].Name();
    }
  }
}

#ifdef USE_MINUIT
// Thermal fit with the yields evaluated at T_kin in PCE (ThermalModelFit::UseTkin)
TEST(PCEFit, YieldsAtTkinMatchPCE) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));

  const long long pdgs[] = { 211, -211, 321, -321, 2212, -2212, 3122, 3312, 3334 };

  // Direct PCE calculation
  ThermalModelIdeal model(&TPS);
  SetChemicalFreezeoutState(model);
  model.CalculateDensities();
  std::vector<double> Ych;
  for (long long pdg : pdgs)
    Ych.push_back(model.GetYield(pdg, Feeddown::StabilityFlag));

  ThermalModelPCE pce(&model);
  pce.SetChemicalFreezeout(model.Parameters(), model.ChemicalPotentials());
  pce.CalculatePCE(kT);

  // Thermal fit with all parameters fixed: the model is evaluated once
  ThermalModelIdeal modelfit(&TPS);
  SetChemicalFreezeoutState(modelfit);
  ThermalModelFit fit(&modelfit);
  fit.UseTkin(true);
  fit.SetParameterValue("T", kTch);
  fit.SetParameterValue("R", std::pow(3. * kVch / (4. * std::acos(-1.)), 1. / 3.));
  fit.SetParameterValue("muB", 0.);
  fit.SetParameterValue("gammaq", 1.);
  fit.SetParameterValue("gammaS", 1.);
  fit.SetParameterValue("Tkin", kT);
  const std::string names[] = { "T", "R", "Rc", "muB", "muQ", "muS", "muC", "gammaq", "gammaS", "gammaC", "Tkin" };
  for (const std::string& name : names)
    fit.SetParameterFitFlag(name, false);

  for (long long pdg : pdgs)
    fit.AddMultiplicity(ExperimentMultiplicity(pdg, 1., 1.));

  fit.PerformFit(false, false);

  // The yields of stable hadrons at T_kin coincide with the direct PCE calculation and with the T_ch yields
  ASSERT_EQ(fit.ModelDataSize(), static_cast<int>(Ych.size()));
  for (int i = 0; i < fit.ModelDataSize(); ++i) {
    EXPECT_NEAR(fit.ModelData(i) / model.GetYield(pdgs[i], Feeddown::StabilityFlag), 1., 1.e-6) << pdgs[i];
    EXPECT_NEAR(fit.ModelData(i) / Ych[i], 1., 1.e-6) << pdgs[i];
  }
}
#endif

// SetParameters() sets the chemical potentials of all species from (muB, muQ, muS, muC)
TEST(ChemicalPotentials, SetParametersRefillsChemicalPotentials) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  model.SetChemicalPotentials(std::vector<double>(TPS.ComponentsNumber(), 0.05));

  // T, muB, muS, muQ
  const ThermalModelParameters params(0.150, 0.200, 0.040, -0.010, 1., kVch);
  model.SetParameters(params);

  for (int i = 0; i < TPS.ComponentsNumber(); ++i) {
    const ThermalParticle& part = TPS.Particles()[i];
    const double mu = part.BaryonCharge() * params.muB + part.ElectricCharge() * params.muQ
                    + part.Strangeness() * params.muS + part.Charm() * params.muC;
    EXPECT_NEAR(model.ChemicalPotentials()[i], mu, 1.e-14) << part.Name();
  }
}

// Chemical potentials set after SetParameters() (as done in PCE) are used and kept by all calculations
TEST(ChemicalPotentials, ExplicitChemicalPotentialsAreKept) {
  ThermalParticleSystem TPS(ListFile("PDG2014/list.dat"));
  ThermalModelIdeal model(&TPS);
  model.SetParameters(ThermalModelParameters(0.100, 0., 0., 0., 1., kVch));

  std::vector<double> mu(TPS.ComponentsNumber(), 0.);
  const int ipi = TPS.PdgToId(211), ip = TPS.PdgToId(2212);
  mu[ipi] = 0.080;
  mu[ip] = 0.400;
  model.SetChemicalPotentials(mu);

  model.CalculateDensities();
  model.CalculateFluctuations();
  model.CalculateTemperatureDerivatives();

  for (size_t i = 0; i < mu.size(); ++i)
    EXPECT_DOUBLE_EQ(model.ChemicalPotentials()[i], mu[i]) << TPS.Particles()[i].Name();

  const int ids[] = { ipi, ip };
  for (int id : ids) {
    const double n = TPS.Particles()[id].Density(model.Parameters(), IdealGasFunctions::ParticleDensity, model.UseWidth(), mu[id]);
    EXPECT_NEAR(model.Densities()[id] / n, 1., 1.e-12) << TPS.Particles()[id].Name();
  }
}
