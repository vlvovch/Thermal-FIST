/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ThermalModelBase.h"

#include <cstdio>
#include <algorithm>

#include <Eigen/Dense>

#include "HRGBase/Utility.h"
#include "HRGBase/ThermalParticleSystem.h"

using namespace Eigen;

using namespace std;

namespace thermalfist {

  ThermalModelBase::ThermalModelBase(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
    m_TPS(TPS_), 
    m_Parameters(params),
    m_UseWidth(false),
    m_PCE(false),
    m_Calculated(false),
    m_FeeddownCalculated(false),
    m_FluctuationsCalculated(false),
    m_GCECalculated(false),
    m_NormBratio(false),
    m_QuantumStats(true),
    m_MaxDiff(0.),
    m_useOpenMP(0)
  {
    if (!Disclaimer::DisclaimerPrinted) 
      Disclaimer::DisclaimerPrinted = Disclaimer::PrintDisclaimer();
    
    m_QBgoal = 0.4;
    m_SBgoal = 50.;
    m_Chem.resize(m_TPS->Particles().size());
    m_Volume = params.V;
    m_densities.resize(m_TPS->Particles().size());
    m_densitiestotal.resize(m_TPS->Particles().size());
    m_densitiesbyfeeddown = std::vector< std::vector<double> >(ParticleDecayType::NumberOfDecayTypes, m_densitiestotal);

    m_wprim.resize(m_TPS->Particles().size());
    m_wtot.resize(m_TPS->Particles().size());
    m_skewprim.resize(m_TPS->Particles().size());
    m_skewtot.resize(m_TPS->Particles().size());
    m_kurtprim.resize(m_TPS->Particles().size());
    m_kurttot.resize(m_TPS->Particles().size());

    m_ConstrainMuB = false;
    m_ConstrainMuC = m_ConstrainMuQ = m_ConstrainMuS = true;

    m_Susc.resize(4);
    for (int i = 0; i < 4; ++i) m_Susc[i].resize(4);

    m_NormBratio = false;
  
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      const ThermalParticle &tpart = m_TPS->Particles()[i];
      for (size_t j = 0; j < tpart.Decays().size(); ++j) {
        if (tpart.DecaysOriginal().size() == tpart.Decays().size() && tpart.Decays()[j].mBratio != tpart.DecaysOriginal()[j].mBratio)
          m_NormBratio = true;
      }
    }

    m_PrimCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(TPS()->ComponentsNumber(), 0.));
    m_TotalCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(TPS()->ComponentsNumber(), 0.));
    m_PrimChargesCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(4, 0.));
    m_FinalChargesCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(4, 0.));

    m_Ensemble = GCE;
    m_InteractionModel = Ideal;

    //SetStatistics(m_QuantumStats);
    //SetCalculationType(IdealGasFunctions::Quadratures);
    SetUseWidth(TPS()->ResonanceWidthIntegrationType());

    ResetCalculatedFlags();

    m_ValidityLog = "";
  }


  void ThermalModelBase::FillVirial(const std::vector<double>& /*ri*/)
  {
  }

  void ThermalModelBase::SetUseWidth(bool useWidth)
  {
    if (!useWidth && m_TPS->ResonanceWidthIntegrationType() != ThermalParticle::ZeroWidth) {
      m_TPS->SetResonanceWidthIntegrationType(ThermalParticle::ZeroWidth);
      //m_TPS->ProcessDecays();
    }
    if (useWidth && m_TPS->ResonanceWidthIntegrationType() == ThermalParticle::ZeroWidth) {
      m_TPS->SetResonanceWidthIntegrationType(ThermalParticle::BWTwoGamma);
    }
    m_UseWidth = useWidth;
  }

  void ThermalModelBase::SetUseWidth(ThermalParticle::ResonanceWidthIntegration type)
  {
    m_UseWidth = (type != ThermalParticle::ZeroWidth);
    m_TPS->SetResonanceWidthIntegrationType(type);
  }


  void ThermalModelBase::SetNormBratio(bool normBratio) {
    if (normBratio != m_NormBratio) {
      m_NormBratio = normBratio;
      if (m_NormBratio) {
        m_TPS->NormalizeBranchingRatios();
      }
      else {
        m_TPS->RestoreBranchingRatios();
      }
    }
  }


  void ThermalModelBase::ResetChemicalPotentials() {
    m_Parameters.muS = m_Parameters.muB / 5.;
    m_Parameters.muQ = -m_Parameters.muB / 50.;
    m_Parameters.muC = 0.;
  }


  void ThermalModelBase::SetParameters(const ThermalModelParameters& params) {
    m_Parameters = params;
    m_Volume = m_Parameters.V;
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetTemperature(double T)
  {
    m_Parameters.T = T;
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetBaryonChemicalPotential(double muB)
  {
    m_Parameters.muB = muB;
    FillChemicalPotentials();
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetElectricChemicalPotential(double muQ)
  {
    m_Parameters.muQ = muQ;
    FillChemicalPotentials();
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetStrangenessChemicalPotential(double muS)
  {
    m_Parameters.muS = muS;
    FillChemicalPotentials();
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetCharmChemicalPotential(double muC)
  {
    m_Parameters.muC = muC;
    FillChemicalPotentials();
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetGammaS(double gammaS)
  {
    m_Parameters.gammaS = gammaS;
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetGammaC(double gammaC)
  {
    m_Parameters.gammaC = gammaC;
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetBaryonCharge(int B)
  {
    m_Parameters.B = B;
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetElectricCharge(int Q)
  {
    m_Parameters.Q = Q;
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetStrangeness(int S)
  {
    m_Parameters.S = S;
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetCharm(int C)
  {
    m_Parameters.C = C;
    ResetCalculatedFlags();
  }

  void ThermalModelBase::DisableMesonMesonVirial()
  {
    for (int i = 0; i < TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS()->ComponentsNumber(); ++j) {
        const ThermalParticle &part1 = TPS()->Particles()[i];
        const ThermalParticle &part2 = TPS()->Particles()[j];
        if (part1.BaryonCharge() == 0 && part2.BaryonCharge() == 0)
          SetVirial(i, j, 0.);
      }
    }
  }

  void ThermalModelBase::DisableMesonMesonAttraction()
  {
    for (int i = 0; i < TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS()->ComponentsNumber(); ++j) {
        const ThermalParticle &part1 = TPS()->Particles()[i];
        const ThermalParticle &part2 = TPS()->Particles()[j];
        if (part1.BaryonCharge() == 0 && part2.BaryonCharge() == 0)
          SetAttraction(i, j, 0.);
      }
    }
  }

  void ThermalModelBase::DisableMesonBaryonVirial()
  {
    for (int i = 0; i < TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS()->ComponentsNumber(); ++j) {
        const ThermalParticle &part1 = TPS()->Particles()[i];
        const ThermalParticle &part2 = TPS()->Particles()[j];
        if ((part1.BaryonCharge() == 0 && part2.BaryonCharge() != 0)
          || (part1.BaryonCharge() != 0 && part2.BaryonCharge() == 0))
          SetVirial(i, j, 0.);
      }
    }
  }

  void ThermalModelBase::DisableMesonBaryonAttraction()
  {
    for (int i = 0; i < TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS()->ComponentsNumber(); ++j) {
        const ThermalParticle &part1 = TPS()->Particles()[i];
        const ThermalParticle &part2 = TPS()->Particles()[j];
        if ((part1.BaryonCharge() == 0 && part2.BaryonCharge() != 0)
          || (part1.BaryonCharge() != 0 && part2.BaryonCharge() == 0))
          SetAttraction(i, j, 0.);
      }
    }
  }

  void ThermalModelBase::DisableBaryonBaryonVirial()
  {
    for (int i = 0; i < TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS()->ComponentsNumber(); ++j) {
        const ThermalParticle &part1 = TPS()->Particles()[i];
        const ThermalParticle &part2 = TPS()->Particles()[j];
        if ((part1.BaryonCharge() > 0 && part2.BaryonCharge() > 0)
          || (part1.BaryonCharge() < 0 && part2.BaryonCharge() < 0))
          SetVirial(i, j, 0.);
      }
    }
  }

  void ThermalModelBase::DisableBaryonBaryonAttraction()
  {
    for (int i = 0; i < TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS()->ComponentsNumber(); ++j) {
        const ThermalParticle &part1 = TPS()->Particles()[i];
        const ThermalParticle &part2 = TPS()->Particles()[j];
        if ((part1.BaryonCharge() > 0 && part2.BaryonCharge() > 0)
          || (part1.BaryonCharge() < 0 && part2.BaryonCharge() < 0))
          SetAttraction(i, j, 0.);
      }
    }
  }

  void ThermalModelBase::DisableBaryonAntiBaryonVirial()
  {
    for (int i = 0; i < TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS()->ComponentsNumber(); ++j) {
        const ThermalParticle &part1 = TPS()->Particles()[i];
        const ThermalParticle &part2 = TPS()->Particles()[j];
        if ((part1.BaryonCharge() > 0 && part2.BaryonCharge() < 0)
          || (part1.BaryonCharge() < 0 && part2.BaryonCharge() > 0))
          SetVirial(i, j, 0.);
      }
    }
  }

  void ThermalModelBase::DisableBaryonAntiBaryonAttraction()
  {
    for (int i = 0; i < TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < TPS()->ComponentsNumber(); ++j) {
        const ThermalParticle &part1 = TPS()->Particles()[i];
        const ThermalParticle &part2 = TPS()->Particles()[j];
        if ((part1.BaryonCharge() > 0 && part2.BaryonCharge() < 0)
          || (part1.BaryonCharge() < 0 && part2.BaryonCharge() > 0))
          SetAttraction(i, j, 0.);
      }
    }
  }

  void ThermalModelBase::SetGammaq(double gammaq)
  {
    m_Parameters.gammaq = gammaq;
    ResetCalculatedFlags();
  }


  void ThermalModelBase::ChangeTPS(ThermalParticleSystem *TPS_) {
    m_TPS = TPS_;
    m_Chem.resize(m_TPS->Particles().size());
    m_densities.resize(m_TPS->Particles().size());
    m_densitiestotal.resize(m_TPS->Particles().size());
    m_densitiesbyfeeddown = std::vector< std::vector<double> >(ParticleDecayType::NumberOfDecayTypes, m_densitiestotal);
    m_wprim.resize(m_TPS->Particles().size());
    m_wtot.resize(m_TPS->Particles().size());
    m_skewprim.resize(m_TPS->Particles().size());
    m_skewtot.resize(m_TPS->Particles().size());
    m_kurtprim.resize(m_TPS->Particles().size());
    m_kurttot.resize(m_TPS->Particles().size());
    m_PrimCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(TPS()->ComponentsNumber(), 0.));
    m_TotalCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(TPS()->ComponentsNumber(), 0.));
    m_PrimChargesCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(4, 0.));
    m_FinalChargesCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(4, 0.));
    ResetCalculatedFlags();
  }

  void ThermalModelBase::SetStatistics(bool stats) {
    m_QuantumStats = stats;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      m_TPS->Particle(i).UseStatistics(stats);
  }

  void ThermalModelBase::SetResonanceWidthIntegrationType(ThermalParticle::ResonanceWidthIntegration type)
  {
    if (!m_UseWidth) {
      printf("**WARNING** ThermalModelBase::SetResonanceWidthIntegrationType: Using resonance widths is switched off!\n");
      m_TPS->SetResonanceWidthIntegrationType(ThermalParticle::BWTwoGamma);
    }
    else
      m_TPS->SetResonanceWidthIntegrationType(type);
  }

  void ThermalModelBase::FillChemicalPotentials() {
    m_Chem.resize(m_TPS->Particles().size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      m_Chem[i] = m_TPS->Particles()[i].BaryonCharge() * m_Parameters.muB + m_TPS->Particles()[i].Strangeness() * m_Parameters.muS + m_TPS->Particles()[i].ElectricCharge() * m_Parameters.muQ + m_TPS->Particles()[i].Charm() * m_Parameters.muC;
  }

  void ThermalModelBase::SetChemicalPotentials(const std::vector<double>& chem)
  {
    if (chem.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::SetChemicalPotentials(const std::vector<double> & chem): size of chem does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_Chem = chem;
  }

  double ThermalModelBase::ChemicalPotential(int i) const
  {
    if (i < 0 || i >= static_cast<int>(m_Chem.size())) {
      printf("**ERROR** ThermalModelBase::ChemicalPotential(int i): i is out of bounds!");
      exit(1);
    }
    return m_Chem[i];
  }

  void ThermalModelBase::SetChemicalPotential(int i, double chem)
  {
    if (i < 0 || i >= static_cast<int>(m_Chem.size())) {
      printf("**ERROR** ThermalModelBase::SetChemicalPotential(int i): i is out of bounds!");
      exit(1);
    }
    m_Chem[i] = chem;
  }

  double ThermalModelBase::FullIdealChemicalPotential(int i) const
  {
    if (i < 0 || i >= static_cast<int>(m_Chem.size())) {
      printf("**ERROR** ThermalModelBase::FullIdealChemicalPotential(int i): i is out of bounds!");
      exit(1);
    }
    
    double ret = ChemicalPotential(i);

    ret += MuShift(i);

    const ThermalParticle& part = m_TPS->Particles()[i];

    if (!(m_Parameters.gammaq == 1.))                  ret += log(m_Parameters.gammaq) * part.AbsoluteQuark() * m_Parameters.T;
    if (!(m_Parameters.gammaS == 1. || part.AbsoluteStrangeness() == 0.))  ret += log(m_Parameters.gammaS) * part.AbsoluteStrangeness() * m_Parameters.T;
    if (!(m_Parameters.gammaC == 1. || part.AbsoluteCharm() == 0.))  ret += log(m_Parameters.gammaC) * part.AbsoluteCharm() * m_Parameters.T;

    return ret;
  }

  void ThermalModelBase::CalculateFeeddown() {
    if (m_UseWidth && m_TPS->ResonanceWidthIntegrationType() == ThermalParticle::eBW) {
      for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
        m_TPS->Particle(i).CalculateThermalBranchingRatios(m_Parameters, m_UseWidth, m_Chem[i] + MuShift(i));
      }
      m_TPS->ProcessDecays();
    }

    // Primordial
    m_densitiesbyfeeddown[static_cast<int>(Feeddown::Primordial)] = m_densities;

    // According to stability flags
    int feed_index = static_cast<int>(Feeddown::StabilityFlag);
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_densitiestotal[i] = m_densities[i];
      //const std::vector< std::pair<double, int> >& decayContributions = m_TPS->Particles()[i].DecayContributionsByFeeddown()[feed_index];
      const ThermalParticleSystem::DecayContributionsToParticle& decayContributions = m_TPS->DecayContributionsByFeeddown()[feed_index][i];
      for (size_t j = 0; j < decayContributions.size(); ++j)
        if (i != decayContributions[j].second) 
          m_densitiestotal[i] += decayContributions[j].first * m_densities[decayContributions[j].second];
    }

    m_densitiesbyfeeddown[feed_index] = m_densitiestotal;

    // Weak, EM, strong
    for (feed_index = static_cast<int>(Feeddown::Weak); feed_index <= static_cast<int>(Feeddown::Strong); ++feed_index) {
      for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
        m_densitiesbyfeeddown[feed_index][i] = m_densities[i];
        //const std::vector< std::pair<double, int> >& decayContributions = m_TPS->Particles()[i].DecayContributionsByFeeddown()[feed_index];
        const ThermalParticleSystem::DecayContributionsToParticle& decayContributions = m_TPS->DecayContributionsByFeeddown()[feed_index][i];
        for (size_t j = 0; j < decayContributions.size(); ++j)
          if (i != decayContributions[j].second)
            m_densitiesbyfeeddown[feed_index][i] += decayContributions[j].first * m_densities[decayContributions[j].second];
      }
    }

    m_FeeddownCalculated = true;
  }


  void ThermalModelBase::ConstrainChemicalPotentials(bool resetInitialValues)
  {
    if (resetInitialValues)
      FixParameters();
    else
      FixParametersNoReset();
  }

  void ThermalModelBase::FixParameters() {
    if (fabs(m_Parameters.muB) < 1e-6 && !m_ConstrainMuB) {
      if (m_ConstrainMuS)
        m_Parameters.muS = 0.;
      if (m_ConstrainMuQ)
        m_Parameters.muQ = 0.;
      if (m_ConstrainMuC)
        m_Parameters.muC = 0.;
      FillChemicalPotentials();
      CalculatePrimordialDensities();
      return;
    }
    if (m_ConstrainMuB) {
      m_Parameters.muB = xMath::mnucleon() / 2.;
    }
    double suppr = 10;
    if (m_Parameters.muB > 0.150) suppr = 8.;
    if (m_Parameters.muB > 0.300) suppr = 7.;
    if (m_Parameters.muB > 0.450) suppr = 6.;
    if (m_Parameters.muB > 0.600) suppr = 6.;
    if (m_Parameters.muB > 0.750) suppr = 5.;
    if (m_Parameters.muB > 0.900) suppr = 4.;
    if (m_Parameters.muB > 1.000) suppr = 3.;
    if (m_ConstrainMuS)
      m_Parameters.muS = m_Parameters.muB / suppr;
    if (m_ConstrainMuQ)
      m_Parameters.muQ = -m_Parameters.muB / suppr / 10.;
    if (m_ConstrainMuC)
      m_Parameters.muC = -m_Parameters.muS;

    FixParametersNoReset();
  }

  void ThermalModelBase::FixParametersNoReset() {
    if (fabs(m_Parameters.muB) < 1e-6 && !m_ConstrainMuB) {
      if (m_ConstrainMuS)
        m_Parameters.muS = 0.;
      if (m_ConstrainMuQ)
        m_Parameters.muQ = 0.;
      if (m_ConstrainMuC)
        m_Parameters.muC = 0.;
      //m_Parameters.muS = m_Parameters.muQ = m_Parameters.muC = 0.;
      FillChemicalPotentials();
      CalculatePrimordialDensities();
      return;
    }

    m_ConstrainMuB &= m_TPS->hasBaryons();
    m_ConstrainMuQ &= (m_TPS->hasCharged() && m_TPS->hasBaryons());
    m_ConstrainMuS &= m_TPS->hasStrange();
    m_ConstrainMuC &= m_TPS->hasCharmed();

    vector<double> x22(4);
    x22[0] = m_Parameters.muB;
    x22[1] = m_Parameters.muQ;
    x22[2] = m_Parameters.muS;
    x22[3] = m_Parameters.muC;
    vector<double> x2(4), xinit(4);
    xinit[0] = x2[0] = m_Parameters.muB;
    xinit[1] = x2[1] = m_Parameters.muQ;
    xinit[2] = x2[2] = m_Parameters.muS;
    xinit[3] = x2[3] = m_Parameters.muC;
    int iter = 0, iterMAX = 2;
    while (iter < iterMAX) {
      BroydenEquationsChem eqs(this);
      BroydenJacobianChem jaco(this);
      BroydenChem broydn(this, &eqs, &jaco);
      Broyden::BroydenSolutionCriterium crit(1.0E-8);
      broydn.Solve(x22, &crit);
      break;
    }
  }

  bool ThermalModelBase::SolveChemicalPotentials(double totB, double totQ, double totS, double totC,
    double muBinit, double muQinit, double muSinit, double muCinit,
    bool ConstrMuB, bool ConstrMuQ, bool ConstrMuS, bool ConstrMuC) {
    if (UsePartialChemicalEquilibrium()) {
      printf("**WARNING** PCE enabled, cannot assume chemical equilibrium to do optimization...");
      return false;
    }

    m_Parameters.muB = muBinit;
    m_Parameters.muS = muSinit;
    m_Parameters.muQ = muQinit;
    m_Parameters.muC = muCinit;
    bool allzero = true;
    allzero &= (totB == 0.0 && ConstrMuB) || (muBinit == 0 && !ConstrMuB);
    allzero &= (totQ == 0.0 && ConstrMuQ) || (muQinit == 0 && !ConstrMuQ);
    allzero &= (totS == 0.0 && ConstrMuS) || (muSinit == 0 && !ConstrMuS);
    allzero &= (totC == 0.0 && ConstrMuC) || (muCinit == 0 && !ConstrMuC);
    if (allzero) {
      m_Parameters.muB = 0.;
      m_Parameters.muS = 0.;
      m_Parameters.muQ = 0.;
      m_Parameters.muC = 0.;
      FillChemicalPotentials();
      CalculatePrimordialDensities();
      return true;
    }
    vector<int> vConstr(4, 1);
    vector<int> vType(4, 0);

    vConstr[0] = m_TPS->hasBaryons() && ConstrMuB;
    vConstr[1] = m_TPS->hasCharged() && ConstrMuQ;
    vConstr[2] = m_TPS->hasStrange() && ConstrMuS;
    vConstr[3] = m_TPS->hasCharmed() && ConstrMuC;

    vType[0] = (int)(totB == 0.0);
    vType[1] = (int)(totQ == 0.0);
    vType[2] = (int)(totS == 0.0);
    vType[3] = (int)(totC == 0.0);

    vector<double> vTotals(4);
    vTotals[0] = totB;
    vTotals[1] = totQ;
    vTotals[2] = totS;
    vTotals[3] = totC;

    vector<double> xin(4, 0.);
    xin[0] = muBinit;
    xin[1] = muQinit;
    xin[2] = muSinit;
    xin[3] = muCinit;

    vector<double> xinactual;
    for (int i = 0; i < 4; ++i) {
      if (vConstr[i]) {
        xinactual.push_back(xin[i]);
      }
    }

    BroydenEquationsChemTotals eqs(vConstr, vType, vTotals, this);
    BroydenJacobianChemTotals jaco(vConstr, vType, vTotals, this);
    Broyden broydn(&eqs, &jaco);
    Broyden::BroydenSolutionCriterium crit(1.0E-8);
    broydn.Solve(xinactual, &crit);

    return (broydn.Iterations() < broydn.MaxIterations());
  }

  void ThermalModelBase::CalculateDensities()
  {
    CalculatePrimordialDensities();

    CalculateFeeddown();
  }

  void ThermalModelBase::ValidateCalculation()
  {
    m_ValidityLog = "";

    char cc[1000];

    m_LastCalculationSuccessFlag = true;
    for (size_t i = 0; i < m_densities.size(); ++i) {
      if (m_densities[i] != m_densities[i]) {
        m_LastCalculationSuccessFlag = false;
      
        sprintf(cc, "**WARNING** Density for particle %lld (%s) is NaN!\n\n", m_TPS->Particle(i).PdgId(), m_TPS->Particle(i).Name().c_str());
        printf("%s", cc);

        m_ValidityLog.append(cc);
      }
      //m_LastCalculationSuccessFlag &= (m_densities[i] == m_densities[i]);
    }
  }

  std::vector<double> ThermalModelBase::CalculateChargeFluctuations(const std::vector<double>& /*chgs*/, int /*order*/)
  {
    printf("**WARNING** %s::CalculateChargeFluctuations(const std::vector<double>& chgs, int order) not implemented!\n", m_TAG.c_str());
    return std::vector<double>();
  }

  double ThermalModelBase::CalculateHadronDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_densities[i];

    return ret;
  }

  double ThermalModelBase::CalculateBaryonDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_TPS->Particles()[i].BaryonCharge() * m_densities[i];

    return ret;
  }

  double ThermalModelBase::CalculateChargeDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_TPS->Particles()[i].ElectricCharge() * m_densities[i];

    return ret;
  }

  double ThermalModelBase::CalculateStrangenessDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_TPS->Particles()[i].Strangeness() * m_densities[i];

    return ret;
  }

  double ThermalModelBase::CalculateCharmDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_TPS->Particles()[i].Charm() * m_densities[i];
    return ret;
  }

  double ThermalModelBase::CalculateAbsoluteBaryonDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += fabs((double)m_TPS->Particles()[i].BaryonCharge()) * m_densities[i];
    return ret;
  }

  double ThermalModelBase::CalculateAbsoluteChargeDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += fabs((double)m_TPS->Particles()[i].ElectricCharge()) * m_densities[i];
    return ret;
  }

  double ThermalModelBase::CalculateAbsoluteStrangenessDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_TPS->Particles()[i].AbsoluteStrangeness() * m_densities[i];
    return ret;
  }

  double ThermalModelBase::CalculateAbsoluteStrangenessDensityModulo() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += fabs((double)m_TPS->Particles()[i].Strangeness()) * m_densities[i];
    return ret;
  }

  double ThermalModelBase::CalculateAbsoluteCharmDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_TPS->Particles()[i].AbsoluteCharm() * m_densities[i];
    return ret;
  }

  double ThermalModelBase::CalculateAbsoluteCharmDensityModulo() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += fabs((double)m_TPS->Particles()[i].Charm()) * m_densities[i];
    return ret;
  }

  double ThermalModelBase::CalculateArbitraryChargeDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_TPS->Particles()[i].ArbitraryCharge() * m_densities[i];
    return ret;
  }


  double ThermalModelBase::GetDensity(long long PDGID, const std::vector<double> *dens)
  {
    if (m_TPS->PdgToId(PDGID) != -1)
      return dens->operator[](m_TPS->PdgToId(PDGID));

    // 1 - Npart
    if (PDGID == 1) return CalculateBaryonDensity();

    // K0S or K0L
    if (PDGID == 310 || PDGID == 130)
      if (m_TPS->PdgToId(311) != -1 && m_TPS->PdgToId(-311) != -1)
        return (dens->operator[](m_TPS->PdgToId(311)) + dens->operator[](m_TPS->PdgToId(-311))) / 2.;

    // Id Pdg code has a trailing zero, try to construct a particle + anti-particle yield
    if (PDGID % 10 == 0) {
      long long tpdgid = PDGID / 10;
      if (m_TPS->PdgToId(tpdgid) != -1 && m_TPS->PdgToId(-tpdgid) != -1)
        return dens->operator[](m_TPS->PdgToId(tpdgid)) + dens->operator[](m_TPS->PdgToId(-tpdgid));
    }

    // 22122112 - nucleons
    if (PDGID == 22122112 && m_TPS->PdgToId(2212) != -1 && m_TPS->PdgToId(2112) != -1)
      return  dens->operator[](m_TPS->PdgToId(2212)) + dens->operator[](m_TPS->PdgToId(2112));

    printf("**WARNING** %s: Density with PDG ID %lld not found!\n", m_TAG.c_str(), PDGID);

    return 0.;
  }

  double ThermalModelBase::GetDensity(long long PDGID, Feeddown::Type feeddown)
  {
    std::vector<double> *dens = NULL;
    if (feeddown == Feeddown::Primordial) 
      dens = &m_densities;
    else if (feeddown == Feeddown::StabilityFlag) 
      dens = &m_densitiestotal;
    else if (static_cast<size_t>(feeddown) < m_densitiesbyfeeddown.size()) 
      dens = &m_densitiesbyfeeddown[static_cast<int>(feeddown)];
    else {
      printf("**WARNING** %s: GetDensity: Unknown feeddown: %d\n", m_TAG.c_str(), static_cast<int>(feeddown));
      return 0.;
    }

    if (!m_Calculated)
      CalculatePrimordialDensities();

    if (feeddown != Feeddown::Primordial && !m_FeeddownCalculated)
      CalculateFeeddown();

    double ret = GetDensity(PDGID, dens);

    // Weak decay contributions from K0S if this particle is not in the list
    if (feeddown == Feeddown::Weak && m_TPS->PdgToId(310) == -1) {
      // pi0
      if (PDGID == 111) {
        ret += 2. * 0.308 * GetDensity(310, dens);
      }
      // pi+,-
      if (PDGID == 211 || PDGID == -211) {
        ret += 0.692 * GetDensity(310, dens);
      }
    }

    return ret;
  }


  std::vector<double> ThermalModelBase::GetIdealGasDensities() const {
    std::vector<double> ret = m_densities;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ret[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
    }
    return ret;
  }

  void ThermalModelBase::ResetCalculatedFlags()
  {
    m_Calculated = false;
    m_FeeddownCalculated = false;
    m_FluctuationsCalculated = false;
    m_GCECalculated = false;
  }

  double ThermalModelBase::ConservedChargeDensity(ConservedCharge::Name chg)
  {
    if (chg == ConservedCharge::BaryonCharge)
      return BaryonDensity();
    if (chg == ConservedCharge::ElectricCharge)
      return ElectricChargeDensity();
    if (chg == ConservedCharge::StrangenessCharge)
      return StrangenessDensity();
    if (chg == ConservedCharge::CharmCharge)
      return CharmDensity();
    return 0.0;
  }

  double ThermalModelBase::ChargedMultiplicity(int type)
  {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.0;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      int tQ = m_TPS->Particles()[i].ElectricCharge();
      bool fl = false;
      if (type == 0  && tQ != 0)
        fl = true;
      if (type == 1  && tQ > 0)
        fl = true;
      if (type == -1 && tQ < 0)
        fl = true;
      if (fl)
        ret += m_densities[i];
    }
    return ret * Volume();
  }

  double ThermalModelBase::ChargedScaledVariance(int type)
  {
    if (!m_FluctuationsCalculated) {
      printf("**WARNING** %s: ChargedScaledVariance(int): Fluctuations were not calculated\n", m_TAG.c_str());
      return 1.;
    }
    double ret = 0.0;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      int tQ = m_TPS->Particles()[i].ElectricCharge();
      bool fl = false;
      if (type == 0 && tQ != 0)
        fl = true;
      if (type == 1 && tQ > 0)
        fl = true;
      if (type == -1 && tQ < 0)
        fl = true;
      if (fl) {
        for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) {
          int tQ2 = m_TPS->Particles()[j].ElectricCharge();
          bool fl2 = false;
          if (type == 0 && tQ2 != 0)
            fl2 = true;
          if (type == 1 && tQ2 > 0)
            fl2 = true;
          if (type == -1 && tQ2 < 0)
            fl2 = true;

          if (fl2) {
            ret += m_PrimCorrel[i][j];
          }
        }
      }
    }
    return ret * m_Parameters.T * Volume() / ChargedMultiplicity(type);
  }

  double ThermalModelBase::ChargedMultiplicityFinal(int type)
  {
    if (!m_Calculated) CalculateDensities();

    int op = type;
    if (type == -1)
      op = 2;

    double ret = 0.0;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ret += m_densities[i] * m_TPS->Particles()[i].Nch()[op];
    }
    return ret * Volume();
  }

  double ThermalModelBase::ChargedScaledVarianceFinal(int type)
  {
    if (!m_FluctuationsCalculated) {
      printf("**WARNING** %s: ChargedScaledVarianceFinal(int): Fluctuations were not calculated\n", m_TAG.c_str());
      return 1.;
    }
    int op = type;
    if (type == -1)
      op = 2;
    double ret = 0.0;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ret += m_densities[i] * Volume() * m_TPS->Particles()[i].DeltaNch()[op];
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) {
        ret += m_PrimCorrel[i][j] * m_Parameters.T * Volume() * m_TPS->Particles()[i].Nch()[op] * m_TPS->Particles()[j].Nch()[op];
      }
    }
    return ret / ChargedMultiplicityFinal(type);
  }

  void ThermalModelBase::CalculateTwoParticleCorrelations() {
    printf("**WARNING** %s: Calculation of two-particle correlations and fluctuations is not implemented\n", m_TAG.c_str());
  }


  void ThermalModelBase::CalculateTwoParticleFluctuationsDecays()
  {
    // Decay contributions here are done according to Eq. (47) in nucl-th/0606036
  
    int NN = m_densities.size();

    // Fluctuations for all
    for (int i = 0; i < NN; ++i)
      //for(int j=0;j<NN;++j) 
    {
      m_TotalCorrel[i][i] = m_PrimCorrel[i][i];
      //for (int r = 0; r < m_TPS->Particles()[i].DecayContributions().size(); ++r) {
      const ThermalParticleSystem::DecayContributionsToParticle& decayContributions = m_TPS->DecayContributionsByFeeddown()[Feeddown::StabilityFlag][i];
      for (size_t r = 0; r < decayContributions.size(); ++r) {
        int rr = decayContributions[r].second;
      
        m_TotalCorrel[i][i] += 2. * m_PrimCorrel[i][rr] * decayContributions[r].first;
      
        for (size_t r2 = 0; r2 < decayContributions.size(); ++r2) {
          int rr2 = decayContributions[r2].second;
          m_TotalCorrel[i][i] += m_PrimCorrel[rr][rr2] * decayContributions[r].first * decayContributions[r2].first;
        }

        // Probabilistic decays
        m_TotalCorrel[i][i] += m_densities[rr] / m_Parameters.T * m_TPS->DecayCumulants()[i][r].first[1];
        //m_TotalCorrel[i][i] += m_densities[rr] / m_Parameters.T * m_TPS->Particles()[i].DecayCumulants()[r].first[1];
        //if (rr != i) { // && !m_TPS->Particles()[r].IsStable()) {
        //  double nij = 0., ni = 0., nj = 0., dnij = 0.;
        //  const auto& decayDistributions = TPS()->ResonanceFinalStatesDistributions()[rr];
        //  for (size_t br = 0; br < decayDistributions.size(); ++br) {
        //    nij += decayDistributions[br].first * decayDistributions[br].second[i] * decayDistributions[br].second[i];
        //    ni += decayDistributions[br].first * decayDistributions[br].second[i];
        //    nj += decayDistributions[br].first * decayDistributions[br].second[i];
        //  }
        //  dnij = nij - ni * nj;
        //  m_TotalCorrel[i][i] += m_densities[rr] / Parameters().T * dnij;
        //}
      }
    }


    // Correlations only for stable
    for (int i = 0; i < NN; ++i) {
      if (m_TPS->Particles()[i].IsStable()) {
        for (int j = 0; j < NN; ++j) {
          if (j != i && m_TPS->Particles()[j].IsStable()) {
            m_TotalCorrel[i][j] = m_PrimCorrel[i][j];

            const ThermalParticleSystem::DecayContributionsToParticle& decayContributionsI = m_TPS->DecayContributionsByFeeddown()[Feeddown::StabilityFlag][i];
            const ThermalParticleSystem::DecayContributionsToParticle& decayContributionsJ = m_TPS->DecayContributionsByFeeddown()[Feeddown::StabilityFlag][j];
            
            for (size_t r = 0; r < decayContributionsJ.size(); ++r) {
              int rr = decayContributionsJ[r].second;
              m_TotalCorrel[i][j] += m_PrimCorrel[i][rr] * decayContributionsJ[r].first;
            }

            for (size_t r = 0; r < decayContributionsI.size(); ++r) {
              int rr = decayContributionsI[r].second;
              m_TotalCorrel[i][j] += m_PrimCorrel[j][rr] * decayContributionsI[r].first;
            }

            for (size_t r = 0; r < decayContributionsI.size(); ++r) {
              int rr = decayContributionsI[r].second;

              for (size_t r2 = 0; r2 < decayContributionsJ.size(); ++r2) {
                int rr2 = decayContributionsJ[r2].second;
                m_TotalCorrel[i][j] += m_PrimCorrel[rr][rr2] * decayContributionsI[r].first * decayContributionsJ[r2].first;
              }
            }

            // Contribution from probabilistic decays
            for (int r = 0; r < m_TPS->ComponentsNumber(); ++r) {
              if (r != i && r != j) { // && !m_TPS->Particles()[r].IsStable()) {
                double nij = 0., ni = 0., nj = 0., dnij = 0.;
                //const ThermalParticle &tpart = m_TPS->Particle(r);
                const ThermalParticleSystem::ResonanceFinalStatesDistribution &decayDistributions = m_TPS->ResonanceFinalStatesDistributions()[r];
                for (size_t br = 0; br < decayDistributions.size(); ++br) {
                  nij += decayDistributions[br].first * decayDistributions[br].second[i] * decayDistributions[br].second[j];
                  ni  += decayDistributions[br].first * decayDistributions[br].second[i];
                  nj  += decayDistributions[br].first * decayDistributions[br].second[j];
                }
                dnij = nij - ni * nj;
                m_TotalCorrel[i][j] += m_densities[r] / m_Parameters.T * dnij;
              }
            }
          }
        }
      }
    }

    for (int i = 0; i < NN; ++i) {
      m_wtot[i] = m_TotalCorrel[i][i];
      if (m_densitiestotal[i] > 0.) m_wtot[i] *= m_Parameters.T / m_densitiestotal[i];
      else m_wtot[i] = 1.;
    }
  }

  double ThermalModelBase::TwoParticleSusceptibilityPrimordial(int i, int j) const
  {
    if (!IsFluctuationsCalculated()) {
      printf("**ERROR** ThermalModelBase::TwoParticleSusceptibilityPrimordial: fluctuations were not computed beforehand! Quitting...\n");
      exit(1);
    }

    return m_PrimCorrel[i][j] / m_Parameters.T / m_Parameters.T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm();
  }

  double ThermalModelBase::TwoParticleSusceptibilityPrimordialByPdg(long long id1, long long id2)
  {
    int i = TPS()->PdgToId(id1);
    int j = TPS()->PdgToId(id2);

    if (i == -1) {
      printf("**WARNING** ThermalModelBase::TwoParticleSusceptibilityPrimordialByPdg: unknown pdg code %lld", id1);
      return 0.;
    }
    if (j == -1) {
      printf("**WARNING** ThermalModelBase::TwoParticleSusceptibilityPrimordialByPdg: unknown pdg code %lld", id2);
      return 0.;
    }

    return TwoParticleSusceptibilityPrimordial(i, j);
  }

  double ThermalModelBase::NetParticleSusceptibilityPrimordialByPdg(long long id1, long long id2)
  {
    int i1 = TPS()->PdgToId(id1);
    int j1 = TPS()->PdgToId(id2);

    if (i1 == -1) {
      printf("**WARNING** ThermalModelBase::NetParticleSusceptibilityPrimordialByPdg: unknown pdg code %lld", id1);
      return 0.;
    }
    if (j1 == -1) {
      printf("**WARNING** ThermalModelBase::NetParticleSusceptibilityPrimordialByPdg: unknown pdg code %lld", id2);
      return 0.;
    }

    int i2 = TPS()->PdgToId(-id1);
    int j2 = TPS()->PdgToId(-id2);

    // Both particles are neutral
    if (i2 == -1 && j2 == -1) {
      return TwoParticleSusceptibilityPrimordial(i1, j1);
    }

    // First particle species is neutral
    if (i2 == -1) {
      return TwoParticleSusceptibilityPrimordial(i1, j1) - TwoParticleSusceptibilityPrimordial(i1, j2);
    }

    // Second particle species is neutral
    if (j2 == -1) {
      return TwoParticleSusceptibilityPrimordial(i1, j1) - TwoParticleSusceptibilityPrimordial(i2, j1);
    }

    // Both particles have antiparticles
    return TwoParticleSusceptibilityPrimordial(i1, j1) + TwoParticleSusceptibilityPrimordial(i2, j2) - TwoParticleSusceptibilityPrimordial(i1, j2) - TwoParticleSusceptibilityPrimordial(i2, j1);
  }

  double ThermalModelBase::TwoParticleSusceptibilityFinal(int i, int j) const
  {
    if (!IsFluctuationsCalculated()) {
      printf("**ERROR** ThermalModelBase::TwoParticleSusceptibilityFinal: fluctuations were not computed beforehand! Quitting...\n");
      exit(1);
    }

    if (!m_TPS->Particle(i).IsStable() || !m_TPS->Particle(j).IsStable()) {
      int tid = i;
      if (!m_TPS->Particle(j).IsStable())
        tid = j;
      printf("**ERROR** ThermalModelBase::TwoParticleSusceptibilityFinal: Particle %d is not stable! Final correlations not computed for unstable particles. Quitting...\n", tid);
      exit(1);
    }

    return m_TotalCorrel[i][j] / m_Parameters.T / m_Parameters.T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm();
  }

  double ThermalModelBase::TwoParticleSusceptibilityFinalByPdg(long long id1, long long id2)
  {
    int i = TPS()->PdgToId(id1);
    int j = TPS()->PdgToId(id2);

    if (i == -1) {
      printf("**WARNING** ThermalModelBase::TwoParticleSusceptibilityFinalByPdg: unknown pdg code %lld", id1);
      return 0.;
    }
    if (j == -1) {
      printf("**WARNING** ThermalModelBase::TwoParticleSusceptibilityFinalByPdg: unknown pdg code %lld", id2);
      return 0.;
    }

    return TwoParticleSusceptibilityFinal(i, j);
  }

  double ThermalModelBase::NetParticleSusceptibilityFinalByPdg(long long id1, long long id2)
  {
    int i1 = TPS()->PdgToId(id1);
    int j1 = TPS()->PdgToId(id2);

    if (i1 == -1) {
      printf("**WARNING** ThermalModelBase::NetParticleSusceptibilityFinalByPdg: unknown pdg code %lld", id1);
      return 0.;
    }
    if (j1 == -1) {
      printf("**WARNING** ThermalModelBase::NetParticleSusceptibilityFinalByPdg: unknown pdg code %lld", id2);
      return 0.;
    }

    int i2 = TPS()->PdgToId(-id1);
    int j2 = TPS()->PdgToId(-id2);

    // Both particles are neutral
    if (i2 == -1 && j2 == -1) {
      return TwoParticleSusceptibilityFinal(i1, j1);
    }

    // First particle species is neutral
    if (i2 == -1) {
      return TwoParticleSusceptibilityFinal(i1, j1) - TwoParticleSusceptibilityFinal(i1, j2);
    }

    // Second particle species is neutral
    if (j2 == -1) {
      return TwoParticleSusceptibilityFinal(i1, j1) - TwoParticleSusceptibilityFinal(i2, j1);
    }

    // Both particles have antiparticles
    return TwoParticleSusceptibilityFinal(i1, j1) + TwoParticleSusceptibilityFinal(i2, j2) - TwoParticleSusceptibilityFinal(i1, j2) - TwoParticleSusceptibilityFinal(i2, j1);
  }

  double ThermalModelBase::PrimordialParticleChargeSusceptibility(int i, ConservedCharge::Name chg) const
  {
    if (!IsFluctuationsCalculated()) {
      printf("**ERROR** ThermalModelBase::PrimordialParticleChargeSusceptibility: fluctuations were not computed beforehand! Quitting...\n");
      exit(1);
    }

    return m_PrimChargesCorrel[i][static_cast<int>(chg)] / m_Parameters.T / m_Parameters.T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm();
  }

  double ThermalModelBase::PrimordialParticleChargeSusceptibilityByPdg(long long id1, ConservedCharge::Name chg)
  {
    int i = TPS()->PdgToId(id1);
    if (i == -1) {
      printf("**WARNING** ThermalModelBase::PrimordialParticleChargeSusceptibilityByPdg: unknown pdg code %lld", id1);
      return 0.;
    }

    return PrimordialParticleChargeSusceptibility(i, chg);
  }

  double ThermalModelBase::PrimordialNetParticleChargeSusceptibilityByPdg(long long id1, ConservedCharge::Name chg)
  {
    int i1 = TPS()->PdgToId(id1);
    if (i1 == -1) {
      printf("**WARNING** ThermalModelBase::PrimordialNetParticleChargeSusceptibilityByPdg: unknown pdg code %lld", id1);
      return 0.;
    }

    int i2 = TPS()->PdgToId(-id1);
    if (i2 == -1)
      return PrimordialParticleChargeSusceptibility(i1, chg);

    return PrimordialParticleChargeSusceptibility(i1, chg) - PrimordialParticleChargeSusceptibility(i2, chg);
  }

  double ThermalModelBase::FinalParticleChargeSusceptibility(int i, ConservedCharge::Name chg) const
  {
    if (!IsFluctuationsCalculated()) {
      printf("**ERROR** ThermalModelBase::FinalParticleChargeSusceptibility: fluctuations were not computed beforehand! Quitting...\n");
      exit(1);
    }

    return m_FinalChargesCorrel[i][static_cast<int>(chg)] / m_Parameters.T / m_Parameters.T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm();
  }

  double ThermalModelBase::FinalParticleChargeSusceptibilityByPdg(long long id1, ConservedCharge::Name chg)
  {
    int i = TPS()->PdgToId(id1);
    if (i == -1) {
      printf("**WARNING** ThermalModelBase::FinalParticleChargeSusceptibilityByPdg: unknown pdg code %lld", id1);
      return 0.;
    }

    return FinalParticleChargeSusceptibility(i, chg);
  }

  double ThermalModelBase::FinalNetParticleChargeSusceptibilityByPdg(long long id1, ConservedCharge::Name chg)
  {
    int i1 = TPS()->PdgToId(id1);
    if (i1 == -1) {
      printf("**WARNING** ThermalModelBase::FinalNetParticleChargeSusceptibilityByPdg: unknown pdg code %lld", id1);
      return 0.;
    }

    int i2 = TPS()->PdgToId(-id1);
    if (i2 == -1)
      return FinalParticleChargeSusceptibility(i1, chg);

    return FinalParticleChargeSusceptibility(i1, chg) - FinalParticleChargeSusceptibility(i2, chg);
  }

  void ThermalModelBase::CalculateSusceptibilityMatrix()
  {
    m_Susc.resize(4);
    for (int i = 0; i < 4; ++i) m_Susc[i].resize(4);

    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        m_Susc[i][j] = 0.;
        for (size_t k = 0; k < m_PrimCorrel.size(); ++k) {
          int c1 = 0;
          if (i == 0) c1 = m_TPS->Particles()[k].BaryonCharge();
          if (i == 1) c1 = m_TPS->Particles()[k].ElectricCharge();
          if (i == 2) c1 = m_TPS->Particles()[k].Strangeness();
          if (i == 3) c1 = m_TPS->Particles()[k].Charm();
          for (size_t kp = 0; kp < m_PrimCorrel.size(); ++kp) {
            int c2 = 0;
            if (j == 0) c2 = m_TPS->Particles()[kp].BaryonCharge();
            if (j == 1) c2 = m_TPS->Particles()[kp].ElectricCharge();
            if (j == 2) c2 = m_TPS->Particles()[kp].Strangeness();
            if (j == 3) c2 = m_TPS->Particles()[kp].Charm();
            m_Susc[i][j] += c1 * c2 * m_PrimCorrel[k][kp];
          }
        }
        m_Susc[i][j] = m_Susc[i][j] / m_Parameters.T / m_Parameters.T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm();
      }
    }
  }

  void ThermalModelBase::CalculateProxySusceptibilityMatrix()
  {
    m_ProxySusc.resize(4);
    for (int i = 0; i < 4; ++i) 
      m_ProxySusc[i].resize(4);

    // Up to 3, no charm here yet
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        m_ProxySusc[i][j] = 0.;
        for (size_t k = 0; k < m_TotalCorrel.size(); ++k) {
          if (m_TPS->Particles()[k].IsStable()) {
            int c1 = 0;
            //if (i == 0) c1 = m_TPS->Particles()[k].BaryonCharge();
            if (i == 0) c1 = 1 * (m_TPS->Particles()[k].PdgId() == 2212) - 1 * (m_TPS->Particles()[k].PdgId() == -2212);
            if (i == 1) c1 = m_TPS->Particles()[k].ElectricCharge();
            //if (i == 1) c1 = 1 * (m_TPS->Particles()[k].PdgId() == 211) - 1 * (m_TPS->Particles()[k].PdgId() == -211);
            if (i == 2) c1 = 1 * (m_TPS->Particles()[k].PdgId() == 321) - 1 * (m_TPS->Particles()[k].PdgId() == -321);
            if (i == 3) c1 = m_TPS->Particles()[k].Charm();
            for (size_t kp = 0; kp < m_TotalCorrel.size(); ++kp) {
              if (m_TPS->Particles()[kp].IsStable()) {
                int c2 = 0;
                //if (j == 0) c2 = m_TPS->Particles()[kp].BaryonCharge();
                if (j == 0) c2 = 1 * (m_TPS->Particles()[kp].PdgId() == 2212) - 1 * (m_TPS->Particles()[kp].PdgId() == -2212);
                if (j == 1) c2 = m_TPS->Particles()[kp].ElectricCharge();
                //if (j == 1) c2 = 1 * (m_TPS->Particles()[kp].PdgId() == 211) - 1 * (m_TPS->Particles()[kp].PdgId() == -211);
                if (j == 2) c2 = 1 * (m_TPS->Particles()[kp].PdgId() == 321) - 1 * (m_TPS->Particles()[kp].PdgId() == -321);
                if (j == 3) c2 = m_TPS->Particles()[kp].Charm();
                m_ProxySusc[i][j] += c1 * c2 * m_TotalCorrel[k][kp];
              }
            }
          }
        }
        m_ProxySusc[i][j] = m_ProxySusc[i][j] / m_Parameters.T / m_Parameters.T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm();
      }
    }

    //printf("chi2netp/chi2skellam = %lf\n", m_ProxySusc[0][0] / (m_densitiestotal[m_TPS->PdgToId(2212)] + m_densitiestotal[m_TPS->PdgToId(-2212)]) * pow(m_Parameters.T * xMath::GeVtoifm(), 3));
    //printf("chi2netpi/chi2skellam = %lf\n", m_ProxySusc[1][1] / (m_densitiestotal[m_TPS->PdgToId(211)] + m_densitiestotal[m_TPS->PdgToId(-211)]) * pow(m_Parameters.T * xMath::GeVtoifm(), 3));
  }

  void ThermalModelBase::CalculateParticleChargeCorrelationMatrix()
  {
    m_PrimChargesCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(4, 0.));
    m_FinalChargesCorrel = std::vector< std::vector<double> >(TPS()->ComponentsNumber(), std::vector<double>(4, 0.));

    for (size_t i = 0; i < TPS()->ComponentsNumber(); ++i) {
      int c1 = 1;
      for (int chg = 0; chg < 4; ++chg) {
        for (size_t j = 0; j < TPS()->ComponentsNumber(); ++j) {
          int c2 = TPS()->Particle(j).ConservedCharge((ConservedCharge::Name)chg);
          m_PrimChargesCorrel[i][chg] += c1 * c2 * m_PrimCorrel[i][j];
        }
      }
    }

    for (size_t i = 0; i < TPS()->ComponentsNumber(); ++i) {
      if (m_TPS->Particles()[i].IsStable()) {
        int c1 = 1;
        for (int chg = 0; chg < 4; ++chg) {
          for (size_t j = 0; j < TPS()->ComponentsNumber(); ++j) {
            if (m_TPS->Particles()[j].IsStable()) {
              int c2 = TPS()->Particle(j).ConservedCharge((ConservedCharge::Name)chg);
              m_FinalChargesCorrel[i][chg] += c1 * c2 * m_TotalCorrel[i][j];
            }
          }
        }
      }
    }
  }

  void ThermalModelBase::CalculateFluctuations() {
    printf("**WARNING** %s: Calculation of fluctuations is not implemented\n", m_TAG.c_str());
  }

  std::vector<double> ThermalModelBase::BroydenEquationsChem::Equations(const std::vector<double>& x)
  {
    std::vector<double> ret(m_N, 0.);

    int i1 = 0;
    if (m_THM->ConstrainMuB()) { m_THM->SetBaryonChemicalPotential(x[i1]); i1++; }
    if (m_THM->ConstrainMuQ()) { m_THM->SetElectricChemicalPotential(x[i1]); i1++; }
    if (m_THM->ConstrainMuS()) { m_THM->SetStrangenessChemicalPotential(x[i1]); i1++; }
    if (m_THM->ConstrainMuC()) { m_THM->SetCharmChemicalPotential(x[i1]); i1++; }
    m_THM->FillChemicalPotentials();
    m_THM->CalculatePrimordialDensities();

    i1 = 0;

    // Baryon charge
    if (m_THM->ConstrainMuB()) {
      double fBd = m_THM->CalculateBaryonDensity();
      double fSd = m_THM->CalculateEntropyDensity();

      ret[i1] = (fBd / fSd - 1. / m_THM->SoverB()) * m_THM->SoverB();

      i1++;
    }

    // Electric charge
    if (m_THM->ConstrainMuQ()) {
      double fBd = m_THM->CalculateBaryonDensity();
      double fQd = m_THM->CalculateChargeDensity();

      // Update: remove division by Q/B to allow for charge neutrality
      ret[i1] = (fQd / fBd - m_THM->QoverB());
      //ret[i1] = (fQd / fBd - m_THM->QoverB()) / m_THM->QoverB();

      i1++;
    }


    // Strangeness
    if (m_THM->ConstrainMuS()) {
      double fSd = m_THM->CalculateStrangenessDensity();
      double fASd = m_THM->CalculateAbsoluteStrangenessDensity();
      if (fASd < 1.e-25)
        fASd = m_THM->CalculateAbsoluteStrangenessDensityModulo();

      ret[i1] = fSd / fASd;

      i1++;
    }


    // Charm
    if (m_THM->ConstrainMuC()) {
      double fCd = m_THM->CalculateCharmDensity();
      double fACd = m_THM->CalculateAbsoluteCharmDensity();
      if (fACd < 1.e-25)
        fACd = m_THM->CalculateAbsoluteCharmDensityModulo();

      ret[i1] = fCd / fACd;

      i1++;
    }

    return ret;
  }

  std::vector<double> ThermalModelBase::BroydenJacobianChem::Jacobian(const std::vector<double>& x)
  {
    int i1 = 0;
    // Analytic calculations of Jacobian not yet supported if entropy per baryon is involved
    if (m_THM->ConstrainMuB()) { 
      printf("**ERROR** Constraining chemical potentials: analytic calculation of the Jacobian not supported if muB is constrained\n");
      exit(1); 
    }

    if (m_THM->ConstrainMuQ()) { m_THM->SetElectricChemicalPotential(x[i1]); i1++; }
    if (m_THM->ConstrainMuS()) { m_THM->SetStrangenessChemicalPotential(x[i1]); i1++; }
    if (m_THM->ConstrainMuC()) { m_THM->SetCharmChemicalPotential(x[i1]); i1++; }
    m_THM->FillChemicalPotentials();
    m_THM->CalculatePrimordialDensities();

    double fBd  = m_THM->CalculateBaryonDensity();
    double fQd  = m_THM->CalculateChargeDensity();
    double fSd  = m_THM->CalculateStrangenessDensity();
    double fASd = m_THM->CalculateAbsoluteStrangenessDensity();
    if (fASd < 1.e-25) {
      fASd = m_THM->CalculateAbsoluteStrangenessDensityModulo();
    }
    double fCd  = m_THM->CalculateCharmDensity();
    double fACd = m_THM->CalculateAbsoluteCharmDensity();
    if (fACd < 1.e-25) {
      fACd = m_THM->CalculateAbsoluteCharmDensityModulo();
    }

    vector<double> chi2s;
    chi2s.resize(m_THM->Densities().size());
    for (size_t i = 0; i < chi2s.size(); ++i) {
      chi2s[i] = m_THM->TPS()->Particle(i).chiDimensionfull(2, m_THM->Parameters(), m_THM->UseWidth(), m_THM->ChemicalPotential(i) + m_THM->MuShift(i))
        * xMath::GeVtoifm3();
    }

    int NNN = 0;
    if (m_THM->ConstrainMuQ()) NNN++;
    if (m_THM->ConstrainMuS()) NNN++;
    if (m_THM->ConstrainMuC()) NNN++;
    MatrixXd ret(NNN, NNN);

    i1 = 0;
    // Electric charge derivatives
    if (m_THM->ConstrainMuQ()) {
      int i2 = 0;

      double d1 = 0., d2 = 0.;

      if (m_THM->ConstrainMuQ()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).ElectricCharge() * m_THM->TPS()->Particle(i).ElectricCharge() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).BaryonCharge() * m_THM->TPS()->Particle(i).ElectricCharge() * chi2s[i];

        // Update: remove division by Q/B to allow for charge neutrality
        ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2);
        //ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2) / m_THM->QoverB();

        i2++;
      }


      if (m_THM->ConstrainMuS()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).ElectricCharge() * m_THM->TPS()->Particle(i).Strangeness() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).BaryonCharge() * m_THM->TPS()->Particle(i).Strangeness() * chi2s[i];

        // Update: remove division by Q/B to allow for charge neutrality
        ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2);
        //ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2) / m_THM->QoverB();

        i2++;
      }


      if (m_THM->ConstrainMuC()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).ElectricCharge() * m_THM->TPS()->Particle(i).Charm() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).BaryonCharge() * m_THM->TPS()->Particle(i).Charm() * chi2s[i];

        // Update: remove division by Q/B to allow for charge neutrality
        ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2);
        //ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2) / m_THM->QoverB();

        i2++;
      }

      i1++;
    }


    // Strangeness derivatives
    if (m_THM->ConstrainMuS()) {
      int i2 = 0;

      double d1 = 0., d2 = 0.;

      if (m_THM->ConstrainMuQ()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).Strangeness() * m_THM->TPS()->Particle(i).ElectricCharge() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).AbsoluteStrangeness() * m_THM->TPS()->Particle(i).ElectricCharge() * chi2s[i];

        ret(i1, i2) = d1 / fASd - fSd / fASd / fASd * d2;

        i2++;
      }


      if (m_THM->ConstrainMuS()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).Strangeness() * m_THM->TPS()->Particle(i).Strangeness() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).AbsoluteStrangeness() * m_THM->TPS()->Particle(i).Strangeness() * chi2s[i];

        ret(i1, i2) = d1 / fASd - fSd / fASd / fASd * d2;

        i2++;
      }


      if (m_THM->ConstrainMuC()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).Strangeness() * m_THM->TPS()->Particle(i).Charm() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).AbsoluteStrangeness() * m_THM->TPS()->Particle(i).Charm() * chi2s[i];

        ret(i1, i2) = d1 / fASd - fSd / fASd / fASd * d2;

        i2++;
      }

      i1++;
    }


    // Charm derivatives
    if (m_THM->ConstrainMuC()) {
      int i2 = 0;

      double d1 = 0., d2 = 0.;

      if (m_THM->ConstrainMuQ()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).Charm() * m_THM->TPS()->Particle(i).ElectricCharge() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).AbsoluteCharm() * m_THM->TPS()->Particle(i).ElectricCharge() *chi2s[i];

        ret(i1, i2) = d1 / fACd - fCd / fACd / fACd * d2;

        i2++;
      }


      if (m_THM->ConstrainMuS()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).Charm() * m_THM->TPS()->Particle(i).Strangeness() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).AbsoluteCharm() * m_THM->TPS()->Particle(i).Strangeness() * chi2s[i];

        ret(i1, i2) = d1 / fACd - fCd / fACd / fACd * d2;

        i2++;
      }


      if (m_THM->ConstrainMuC()) {
        d1 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d1 += m_THM->TPS()->Particle(i).Charm() * m_THM->TPS()->Particle(i).Charm() * chi2s[i];

        d2 = 0.;
        for (size_t i = 0; i < chi2s.size(); ++i)
          d2 += m_THM->TPS()->Particle(i).AbsoluteCharm() * m_THM->TPS()->Particle(i).Charm() * chi2s[i];

        ret(i1, i2) = d1 / fACd - fCd / fACd / fACd * d2;

        i2++;
      }

      i1++;
    }

    std::vector<double> retVec(NNN*NNN, 0);
    for (int i = 0; i < NNN; ++i)
      for (int j = 0; j < NNN; ++j)
        retVec[i*NNN + j] = ret(i, j);


    return retVec;
  }

  std::vector<double> ThermalModelBase::BroydenChem::Solve(const std::vector<double>& x0, BroydenSolutionCriterium * solcrit, int max_iterations)
  {
    if (m_Equations == NULL) {
      printf("**ERROR** Broyden::Solve: Equations to solve not specified!\n");
      exit(1);
    }

    int NNN = 0;
    std::vector<double> xcur;
    if (m_THM->ConstrainMuB()) { xcur.push_back(x0[0]); NNN++; }
    if (m_THM->ConstrainMuQ()) { xcur.push_back(x0[1]); NNN++; }
    if (m_THM->ConstrainMuS()) { xcur.push_back(x0[2]); NNN++; }
    if (m_THM->ConstrainMuC()) { xcur.push_back(x0[3]); NNN++; }
    if (NNN == 0) {
      m_THM->FillChemicalPotentials();
      return xcur;
    }

    m_Equations->SetDimension(NNN);

    m_MaxIterations = max_iterations;

    BroydenSolutionCriterium *SolutionCriterium = solcrit;
    bool UseDefaultSolutionCriterium = false;
    if (SolutionCriterium == NULL) {
      SolutionCriterium = new BroydenSolutionCriterium(TOL);
      UseDefaultSolutionCriterium = true;
    }

    BroydenJacobian *JacobianInUse = m_Jacobian;
    bool UseDefaultJacobian = false;
    if (JacobianInUse == NULL || m_THM->ConstrainMuB()) {
      JacobianInUse = new BroydenJacobian(m_Equations);
      UseDefaultJacobian = true;
    }
    m_Iterations = 0;
    double &maxdiff = m_MaxDifference;
    int N = m_Equations->Dimension();

    

    std::vector<double> tmpvec, xdeltavec = xcur;
    VectorXd xold(N), xnew(N), xdelta(N);
    VectorXd fold(N), fnew(N), fdelta(N);

    xold = VectorXd::Map(&xcur[0], xcur.size());

    MatrixXd Jac = Eigen::Map< Matrix<double, Dynamic, Dynamic, RowMajor> >(&JacobianInUse->Jacobian(xcur)[0], N, N);

    bool constrmuB = m_THM->ConstrainMuB();
    bool constrmuQ = m_THM->ConstrainMuQ();
    bool constrmuS = m_THM->ConstrainMuS();
    bool constrmuC = m_THM->ConstrainMuC();
    bool repeat = false;
    NNN = 0;
    if (m_THM->ConstrainMuB()) {
      for (int j = 0; j < Jac.rows(); ++j)
        if (Jac(NNN, j) > 1.e8) { repeat = true; m_THM->ConstrainMuB(false); }
      double  S = m_THM->CalculateEntropyDensity();
      double nB = m_THM->CalculateBaryonDensity();
      if (abs(S) < 1.e-25 || abs(nB) < 1.e-25) { repeat = true; m_THM->ConstrainMuB(false); }
      NNN++;
    }
    if (m_THM->ConstrainMuQ()) {
      for (int j = 0; j < Jac.rows(); ++j)
        if (Jac(NNN, j) > 1.e8) { repeat = true; m_THM->ConstrainMuQ(false); }
      double nQ = m_THM->CalculateChargeDensity();
      double nB = m_THM->CalculateBaryonDensity();
      if (abs(nQ) < 1.e-25 || abs(nB) < 1.e-25) { repeat = true; m_THM->ConstrainMuQ(false); }
      NNN++;
    }
    if (m_THM->ConstrainMuS()) {
      for (int j = 0; j < Jac.rows(); ++j)
        if (Jac(NNN, j) > 1.e8) { repeat = true; m_THM->ConstrainMuS(false); }

      double nS = m_THM->CalculateAbsoluteStrangenessDensity();
      if (abs(nS) < 1.e-25) { nS = m_THM->CalculateAbsoluteStrangenessDensityModulo();  }
      if (abs(nS) < 1.e-25) { repeat = true; m_THM->ConstrainMuS(false); }
      NNN++;
    }
    if (m_THM->ConstrainMuC()) {
      for (int j = 0; j < Jac.rows(); ++j)
        if (Jac(NNN, j) > 1.e8) { repeat = true; m_THM->ConstrainMuC(false); }
      double nC = m_THM->CalculateAbsoluteCharmDensity();
      if (abs(nC) < 1.e-25) { nC = m_THM->CalculateAbsoluteCharmDensityModulo(); }
      if (abs(nC) < 1.e-25) { repeat = true; m_THM->ConstrainMuC(false); }
      NNN++;
    }
    if (repeat) {
      std::vector<double> ret = Solve(x0, solcrit, max_iterations);
      m_THM->ConstrainMuQ(constrmuB);
      m_THM->ConstrainMuQ(constrmuQ);
      m_THM->ConstrainMuS(constrmuS);
      m_THM->ConstrainMuC(constrmuC);
      return ret;
    }

    if (Jac.determinant() == 0.0)
    {
      printf("**WARNING** Singular Jacobian in Broyden::Solve\n");
      return xcur;
    }

    MatrixXd Jinv = Jac.inverse();
    tmpvec = m_Equations->Equations(xcur);
    fold = VectorXd::Map(&tmpvec[0], tmpvec.size());

    for (m_Iterations = 1; m_Iterations < max_iterations; ++m_Iterations) {
      xnew = xold - Jinv * fold;

      VectorXd::Map(&xcur[0], xcur.size()) = xnew;

      tmpvec = m_Equations->Equations(xcur);
      fnew = VectorXd::Map(&tmpvec[0], tmpvec.size());


      maxdiff = 0.;
      for (size_t i = 0; i < xcur.size(); ++i) {
        maxdiff = std::max(maxdiff, fabs(fnew[i]));
      }

      xdelta = xnew - xold;
      fdelta = fnew - fold;

      VectorXd::Map(&xdeltavec[0], xdeltavec.size()) = xdelta;

      if (SolutionCriterium->IsSolved(xcur, tmpvec, xdeltavec))
        break;

      if (!m_UseNewton) // Use Broyden's method
      {

        double norm = 0.;
        for (int i = 0; i < N; ++i)
          for (int j = 0; j < N; ++j)
            norm += xdelta[i] * Jinv(i, j) * fdelta[j];
        VectorXd p1(N);
        p1 = (xdelta - Jinv * fdelta);
        for (int i = 0; i < N; ++i) p1[i] *= 1. / norm;
        Jinv = Jinv + (p1 * xdelta.transpose()) * Jinv;
      }
      else // Use Newton's method
      {
        Jac = Eigen::Map< Matrix<double, Dynamic, Dynamic, RowMajor> >(&JacobianInUse->Jacobian(xcur)[0], N, N);
        Jinv = Jac.inverse();
      }

      xold = xnew;
      fold = fnew;
    }

    if (m_Iterations == max_iterations) {
      printf("**WARNING** Reached maximum number of iterations in Broyden procedure\n");
    }

    if (UseDefaultSolutionCriterium) {
      delete SolutionCriterium;
      SolutionCriterium = NULL;
    }
    if (UseDefaultJacobian) {
      delete JacobianInUse;
      JacobianInUse = NULL;
    }
    return xcur;
  }


  ThermalModelBase::BroydenEquationsChemTotals::BroydenEquationsChemTotals(const std::vector<int>& vConstr, const std::vector<int>& vType, const std::vector<double>& vTotals, ThermalModelBase * model) :
    BroydenEquations(), m_Constr(vConstr), m_Type(vType), m_Totals(vTotals), m_THM(model)
  {
    m_N = 0;
    for (size_t i = 0; i < m_Constr.size(); ++i)
      m_N += m_Constr[i];
  }

  std::vector<double> ThermalModelBase::BroydenEquationsChemTotals::Equations(const std::vector<double>& x)
  {
    std::vector<double> ret(m_N, 0.);

    int i1 = 0;
    for (int i = 0; i < 4; ++i) {
      if (m_Constr[i]) {
        if (i == 0) m_THM->SetBaryonChemicalPotential(x[i1]);
        if (i == 1) m_THM->SetElectricChemicalPotential(x[i1]);
        if (i == 2) m_THM->SetStrangenessChemicalPotential(x[i1]);
        if (i == 3) m_THM->SetCharmChemicalPotential(x[i1]);
        i1++;
      }
    }
    m_THM->FillChemicalPotentials();
    m_THM->CalculatePrimordialDensities();

    vector<double> dens(4, 0.), absdens(4, 0.);
    if (m_Constr[0]) {
      dens[0] = m_THM->CalculateBaryonDensity();
      absdens[0] = m_THM->CalculateAbsoluteBaryonDensity();
    }
    if (m_Constr[1]) {
      dens[1] = m_THM->CalculateChargeDensity();
      absdens[1] = m_THM->CalculateAbsoluteChargeDensity();
    }
    if (m_Constr[2]) {
      dens[2] = m_THM->CalculateStrangenessDensity();
      absdens[2] = m_THM->CalculateAbsoluteStrangenessDensity();
      if (absdens[2] < 1.e-25) {
        absdens[2] = m_THM->CalculateAbsoluteStrangenessDensityModulo();
      }
    }
    if (m_Constr[3]) {
      dens[3] = m_THM->CalculateCharmDensity();
      absdens[3] = m_THM->CalculateAbsoluteCharmDensity();
      if (absdens[3] < 1.e-25) {
        absdens[3] = m_THM->CalculateAbsoluteCharmDensityModulo();
      }
    }

    i1 = 0;

    for (int i = 0; i < 4; ++i) {
      if (m_Constr[i]) {
        if (m_Type[i] == 0)
          ret[i1] = (dens[i] * m_THM->Parameters().V - m_Totals[i]) / m_Totals[i];
        else
          ret[i1] = dens[i] / absdens[i];
        i1++;
      }
    }

    return ret;
  }

  std::vector<double> ThermalModelBase::BroydenJacobianChemTotals::Jacobian(const std::vector<double>& x)
  {
    int NNN = 0;
    for (int i = 0; i < 4; ++i) NNN += m_Constr[i];

    int i1 = 0;

    for (int i = 0; i < 4; ++i) {
      if (m_Constr[i]) {
        if (i == 0) m_THM->SetBaryonChemicalPotential(x[i1]);
        if (i == 1) m_THM->SetElectricChemicalPotential(x[i1]);
        if (i == 2) m_THM->SetStrangenessChemicalPotential(x[i1]);
        if (i == 3) m_THM->SetCharmChemicalPotential(x[i1]);
        i1++;
      }
    }

    m_THM->FillChemicalPotentials();
    m_THM->CalculatePrimordialDensities();

    vector<double> dens(4, 0.), absdens(4, 0.);
    if (m_Constr[0]) {
      dens[0] = m_THM->CalculateBaryonDensity();
      absdens[0] = m_THM->CalculateAbsoluteBaryonDensity();
    }
    if (m_Constr[1]) {
      dens[1] = m_THM->CalculateChargeDensity();
      absdens[1] = m_THM->CalculateAbsoluteChargeDensity();
    }
    if (m_Constr[2]) {
      dens[2] = m_THM->CalculateStrangenessDensity();
      absdens[2] = m_THM->CalculateAbsoluteStrangenessDensity();
      if (absdens[2] < 1.e-25) {
        absdens[2] = m_THM->CalculateAbsoluteStrangenessDensityModulo();
      }
    }
    if (m_Constr[3]) {
      dens[3] = m_THM->CalculateCharmDensity();
      absdens[3] = m_THM->CalculateAbsoluteCharmDensity();
      if (absdens[3] < 1.e-25) {
        absdens[3] = m_THM->CalculateAbsoluteCharmDensityModulo();
      }
    }

    vector<double> chi2s;
    chi2s.resize(m_THM->Densities().size());
    for (size_t i = 0; i < chi2s.size(); ++i) {
      chi2s[i] = m_THM->TPS()->Particle(i).chiDimensionfull(2, m_THM->Parameters(), m_THM->UseWidth(), m_THM->ChemicalPotential(i) + m_THM->MuShift(i))
        * xMath::GeVtoifm3();
    }

    vector< vector<double> > deriv(4, vector<double>(4)), derivabs(4, vector<double>(4));
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        deriv[i][j] = 0.;
        for (size_t part = 0; part < chi2s.size(); ++part)
          deriv[i][j] += m_THM->TPS()->Particles()[part].GetCharge(i) * m_THM->TPS()->Particles()[part].GetCharge(j) * chi2s[part];

        derivabs[i][j] = 0.;
        for (size_t part = 0; part < chi2s.size(); ++part)
          derivabs[i][j] += m_THM->TPS()->Particles()[part].GetAbsCharge(i) * m_THM->TPS()->Particles()[part].GetCharge(j) * chi2s[part];
      }


    MatrixXd ret(NNN, NNN);

    i1 = 0;

    for (int i = 0; i < 4; ++i) {
      if (m_Constr[i]) {
        int i2 = 0;
        for (int j = 0; j < 4; ++j)
          if (m_Constr[j]) {
            ret(i1, i2) = 0.;
            if (m_Type[i] == 0)
              ret(i1, i2) = deriv[i][j] * m_THM->Parameters().V / m_Totals[i];
            else
              ret(i1, i2) = deriv[i][j] / absdens[i] - dens[i] / absdens[i] / absdens[i] * derivabs[i][j];
            i2++;
          }
        i1++;
      }
    }

    std::vector<double> retVec(NNN*NNN, 0);
    for (int i = 0; i < NNN; ++i)
      for (int j = 0; j < NNN; ++j)
        retVec[i*NNN + j] = ret(i, j);


    return retVec;
  }

} // namespace thermalfist