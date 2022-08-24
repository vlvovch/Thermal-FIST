/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/EventGeneratorBase.h"

#include <functional>
#include <algorithm>


#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonical.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGEV/ThermalModelEVCrossterms.h"
#include "HRGEV/ExcludedVolumeHelper.h"
#include "HRGVDW/ThermalModelVDW.h"
#include "HRGEventGenerator/ParticleDecaysMC.h"

namespace thermalfist {

  int EventGeneratorBase::fCEAccepted, EventGeneratorBase::fCETotal;
  double EventGeneratorBase::m_LastWeight = 1., EventGeneratorBase::m_LastLogWeight = 0., EventGeneratorBase::m_LastNormWeight = 1.;

  std::vector<double> LorentzBoost(const std::vector<double>& fourvector, double vx, double vy, double vz)
  {
    std::vector<double> ret(4, 0);
    double v2 = vx * vx + vy * vy + vz * vz;
    if (v2 == 0.0)
      return fourvector;
    double gamma = 1. / sqrt(1. - v2);

    const double& r0 = fourvector[0];
    const double& rx = fourvector[1];
    const double& ry = fourvector[2];
    const double& rz = fourvector[3];

    ret[0] = gamma * r0 - gamma * (vx * rx + vy * ry + vz * rz);
    ret[1] = -gamma * vx * r0 + (1. + (gamma - 1.) * vx * vx / v2) * rx +
      (gamma - 1.) * vx * vy / v2 * ry + (gamma - 1.) * vx * vz / v2 * rz;
    ret[2] = -gamma * vy * r0 + (1. + (gamma - 1.) * vy * vy / v2) * ry +
      (gamma - 1.) * vy * vx / v2 * rx + (gamma - 1.) * vy * vz / v2 * rz;
    ret[3] = -gamma * vz * r0 + (1. + (gamma - 1.) * vz * vz / v2) * rz +
      (gamma - 1.) * vz * vx / v2 * rx + (gamma - 1.) * vz * vy / v2 * ry;

    return ret;
  }

  EventGeneratorBase::EventGeneratorBase() : 
    m_THM(NULL), m_ParametersSet(false),
    m_MeanB(0.), m_MeanAB(0.),
    m_MeanSM(0.), m_MeanASM(0.),
    m_MeanCM(0.), m_MeanACM(0.),
    m_MeanCHRM(0.), m_MeanACHRM(0.),
    m_MeanCHRMM(0.), m_MeanACHRMM(0.)
  {
    fCEAccepted = fCETotal = 0; 
  }

  EventGeneratorBase::~EventGeneratorBase()
  {
    ClearMomentumGenerators();
  }

  void EventGeneratorBase::ClearMomentumGenerators()
  {
    for (size_t i = 0; i < m_MomentumGens.size(); ++i) {
      if (m_MomentumGens[i] != NULL)
        delete m_MomentumGens[i];
    }
    m_MomentumGens.resize(0);

    for (size_t i = 0; i < m_BWGens.size(); ++i) {
      if (m_BWGens[i] != NULL)
        delete m_BWGens[i];
    }
    m_BWGens.resize(0);
  }

  //void EventGeneratorBase::SetConfiguration(const ThermalModelParameters& params, EventGeneratorConfiguration::Ensemble ensemble, EventGeneratorConfiguration::ModelType modeltype, ThermalParticleSystem *TPS, ThermalModelBase *THMEVVDW)
  void EventGeneratorBase::SetConfiguration(ThermalParticleSystem *TPS,
      const EventGeneratorConfiguration& config)
  {
    m_Config = config;
    SetEVUseSPR(config.fUseEVUseSPRApproximation);

    if (TPS == NULL)
      return;


    if (m_Config.fEnsemble == EventGeneratorConfiguration::CCE) {
      m_Config.CFOParameters.muC = 0.;
    }
    else if (m_Config.fEnsemble == EventGeneratorConfiguration::SCE) {
      m_Config.CFOParameters.muS = 0.;
      m_Config.CFOParameters.muC = 0.;
    }
    else if (m_Config.fEnsemble == EventGeneratorConfiguration::CE) {
      if (m_Config.CanonicalB) 
        m_Config.CFOParameters.muB = 0.;
      if (m_Config.CanonicalS)
        m_Config.CFOParameters.muS = 0.;
      if (m_Config.CanonicalQ)
        m_Config.CFOParameters.muQ = 0.;
      if (m_Config.CanonicalC)
        m_Config.CFOParameters.muC = 0.;
    }

    if (m_Config.fModelType == EventGeneratorConfiguration::PointParticle) {
      m_THM = new ThermalModelIdeal(TPS, m_Config.CFOParameters);
    }
    else if (m_Config.fModelType == EventGeneratorConfiguration::DiagonalEV) {
      m_THM = new ThermalModelEVDiagonal(TPS, m_Config.CFOParameters);
    }
    else if (m_Config.fModelType == EventGeneratorConfiguration::CrosstermsEV) {
      m_THM = new ThermalModelEVCrossterms(TPS, m_Config.CFOParameters);
    }
    else if (m_Config.fModelType == EventGeneratorConfiguration::QvdW) {
      m_THM = new ThermalModelVDWFull(TPS, m_Config.CFOParameters);
    }

    m_THM->SetUseWidth(TPS->ResonanceWidthIntegrationType());

    m_THM->ConstrainMuB(false);
    m_THM->ConstrainMuQ(false);
    m_THM->ConstrainMuS(false);
    m_THM->ConstrainMuC(false);

    if (!m_Config.fUsePCE)
      m_THM->FillChemicalPotentials();
    else
      m_THM->SetChemicalPotentials(m_Config.fPCEChems);

    if (m_Config.fModelType != EventGeneratorConfiguration::PointParticle) {
      for (size_t i = 0; i < m_THM->Densities().size(); ++i) {
        for (size_t j = 0; j < m_THM->Densities().size(); ++j) {
          if (m_Config.bij.size() == m_THM->Densities().size()
            && m_Config.bij[i].size() == m_THM->Densities().size())
            m_THM->SetVirial(i, j, m_Config.bij[i][j]);

          if (m_Config.aij.size() == m_THM->Densities().size()
            && m_Config.aij[i].size() == m_THM->Densities().size())
            m_THM->SetAttraction(i, j, m_Config.aij[i][j]);
        }
      }
    }

    if (m_Config.fEnsemble == EventGeneratorConfiguration::CE) {
      // The following procedure currently not used, 
      // but can be considered if SolveChemicalPotentials routine fails
      // Note to take into account PCE possiblity if this option is considered
      // TODO: Properly for mixed-canonical ensembles and charm
      if (0 && !(m_Config.B == 0 && m_Config.Q == 0 && m_Config.S == 0) && !m_Config.fUsePCE) {
        if (m_Config.S == 0 && !(m_Config.Q == 0 || m_Config.B == 0)) {
          double QBrat = (double)(m_Config.Q) / m_Config.B;

          double left = 0.000, right = 0.900, center;

          m_THM->SetBaryonChemicalPotential(left);
          m_THM->SetQoverB(QBrat);
          m_THM->FixParameters();
          m_THM->CalculatePrimordialDensities();
          double valleft = m_THM->CalculateBaryonDensity() * m_THM->Volume() - m_Config.B;

          m_THM->SetBaryonChemicalPotential(right);
          m_THM->SetQoverB(QBrat);
          m_THM->FixParameters();
          m_THM->CalculatePrimordialDensities();
          double valright = m_THM->CalculateBaryonDensity() * m_THM->Volume() - m_Config.B;

          double valcenter;

          while ((right - left) > 0.00001) {
            center = (left + right) / 2.;
            m_THM->SetBaryonChemicalPotential(center);
            m_THM->SetQoverB(QBrat);
            m_THM->FixParameters();
            m_THM->CalculatePrimordialDensities();
            valcenter = m_THM->CalculateBaryonDensity() * m_THM->Volume() - m_Config.B;

            if (valleft*valcenter > 0.) {
              left = center;
              valleft = valcenter;
            }
            else {
              right = center;
              valright = valcenter;
            }
          }

          m_THM->SetBaryonChemicalPotential(center);
          m_THM->SetQoverB(QBrat);
          m_THM->FixParameters();

          m_Config.CFOParameters.muB = m_THM->Parameters().muB;
          m_Config.CFOParameters.muS = m_THM->Parameters().muS;
          m_Config.CFOParameters.muQ = m_THM->Parameters().muQ;

          m_THM->FillChemicalPotentials();
        }
      }

      if (!m_Config.fUsePCE) 
      {
        bool solve_chems = m_THM->SolveChemicalPotentials(m_Config.B, m_Config.Q, m_Config.S, m_Config.C,
          m_THM->Parameters().muB, m_THM->Parameters().muQ, m_THM->Parameters().muS, m_THM->Parameters().muC,
          m_Config.CanonicalB, m_Config.CanonicalQ, m_Config.CanonicalS, m_Config.CanonicalC);
        if (m_THM->Parameters().muB != m_THM->Parameters().muB ||
          m_THM->Parameters().muQ != m_THM->Parameters().muQ ||
          m_THM->Parameters().muS != m_THM->Parameters().muS ||
          m_THM->Parameters().muC != m_THM->Parameters().muC ||
          !solve_chems)
        {
          printf("**WARNING** Could not constrain chemical potentials. Setting all for exactly conserved charges to zero...\n");
          if (m_Config.CanonicalB)
            m_THM->SetBaryonChemicalPotential(0.);
          else
            m_THM->SetBaryonChemicalPotential(m_Config.CFOParameters.muB);

          if (m_Config.CanonicalQ)
            m_THM->SetElectricChemicalPotential(0.);
          else
            m_THM->SetElectricChemicalPotential(m_Config.CFOParameters.muQ);

          if (m_Config.CanonicalS)
            m_THM->SetStrangenessChemicalPotential(0.);
          else
            m_THM->SetStrangenessChemicalPotential(m_Config.CFOParameters.muS);

          if (m_Config.CanonicalC)
            m_THM->SetCharmChemicalPotential(0.);
          else
            m_THM->SetCharmChemicalPotential(m_Config.CFOParameters.muC);

          m_THM->FillChemicalPotentials();
        }
        m_Config.CFOParameters.muB = m_THM->Parameters().muB;
        m_Config.CFOParameters.muS = m_THM->Parameters().muS;
        m_Config.CFOParameters.muQ = m_THM->Parameters().muQ;
        m_Config.CFOParameters.muC = m_THM->Parameters().muC;

        //std::cout << "Chemical potentials constrained!\n" 
        //  <<  "muB = " << m_Config.CFOParameters.muB 
        //  << " muQ = " << m_Config.CFOParameters.muQ 
        //  << " muS = " << m_Config.CFOParameters.muS 
        //  << " muC = " << m_Config.CFOParameters.muC << "\n";

        //std::cout << "B = " << m_THM->CalculateBaryonDensity() * m_THM->Parameters().V
        //  << " Q = " << m_THM->CalculateChargeDensity() * m_THM->Parameters().V
        //  << " S = " << m_THM->CalculateStrangenessDensity() * m_THM->Parameters().V
        //  << " C = " << m_THM->CalculateCharmDensity() * m_THM->Parameters().V
        //  << "\n";
      }
    }

    if (m_Config.fEnsemble == EventGeneratorConfiguration::SCE && !m_Config.fUsePCE) {
      m_THM->ConstrainMuS(true);
      m_THM->ConstrainMuC(true);
      m_THM->ConstrainChemicalPotentials();
      m_THM->FillChemicalPotentials();
      m_Config.CFOParameters.muS = m_THM->Parameters().muS;
      m_Config.CFOParameters.muC = m_THM->Parameters().muC;
    }

    if (m_Config.fEnsemble == EventGeneratorConfiguration::CCE && !m_Config.fUsePCE) {
      m_THM->ConstrainMuC(true);
      m_THM->ConstrainChemicalPotentials();
      m_THM->FillChemicalPotentials();
      m_Config.CFOParameters.muC = m_THM->Parameters().muC;
    }

    //m_THM->CalculateDensitiesGCE();
    m_THM->CalculatePrimordialDensities();
    m_DensitiesIdeal = m_THM->GetIdealGasDensities();

    if (m_Config.fEnsemble != EventGeneratorConfiguration::GCE)
      PrepareMultinomials();

    m_Radii = ComputeEVRadii();
  }


  void EventGeneratorBase::PrepareMultinomials() {
    if (!m_THM->IsCalculated()) {
      m_THM->CalculatePrimordialDensities();
    }

    m_Baryons.resize(0);
    m_AntiBaryons.resize(0);
    m_StrangeMesons.resize(0);
    m_AntiStrangeMesons.resize(0);
    m_ChargeMesons.resize(0);
    m_AntiChargeMesons.resize(0);
    m_CharmMesons.resize(0);
    m_AntiCharmMesons.resize(0);
    m_CharmAll.resize(0);
    m_AntiCharmAll.resize(0);
    m_MeanB = 0.;
    m_MeanAB = 0.;
    m_MeanSM = 0.;
    m_MeanASM = 0.;
    m_MeanCM = 0.;
    m_MeanACM = 0.;
    m_MeanCHRMM = 0.;
    m_MeanACHRMM = 0.;
    m_MeanCHRM = 0.;
    m_MeanACHRM = 0.;

    std::vector<double> yields = m_THM->Densities();
    for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i)
      yields[i] *= m_THM->Volume();

    for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
      if (m_THM->TPS()->Particles()[i].BaryonCharge() == 1) {
        m_Baryons.push_back(std::make_pair(yields[i], i));
        m_MeanB += yields[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == -1) {
        m_AntiBaryons.push_back(std::make_pair(yields[i], i));
        m_MeanAB += yields[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 1) {
        m_StrangeMesons.push_back(std::make_pair(yields[i], i));
        m_MeanSM += yields[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == -1) {
        m_AntiStrangeMesons.push_back(std::make_pair(yields[i], i));
        m_MeanASM += yields[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == 1) {
        m_ChargeMesons.push_back(std::make_pair(yields[i], i));
        m_MeanCM += yields[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == -1) {
        m_AntiChargeMesons.push_back(std::make_pair(yields[i], i));
        m_MeanACM += yields[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == 0 && m_THM->TPS()->Particles()[i].Charm() == 1) {
        m_CharmMesons.push_back(std::make_pair(yields[i], i));
        m_MeanCHRMM += yields[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == 0 && m_THM->TPS()->Particles()[i].Charm() == -1) {
        m_AntiCharmMesons.push_back(std::make_pair(yields[i], i));
        m_MeanACHRMM += yields[i];
      }

      if (m_THM->TPS()->Particles()[i].Charm() == 1) {
        m_CharmAll.push_back(std::make_pair(yields[i], i));
        m_MeanCHRM += yields[i];
      }
      else if (m_THM->TPS()->Particles()[i].Charm() == -1) {
        m_AntiCharmAll.push_back(std::make_pair(yields[i], i));
        m_MeanACHRM += yields[i];
      }
    }

    // sort in descending order and convert to prefix sums
    std::sort(m_Baryons.begin(), m_Baryons.end(), std::greater< std::pair<double, int> >());
    std::sort(m_AntiBaryons.begin(), m_AntiBaryons.end(), std::greater< std::pair<double, int> >());
    std::sort(m_StrangeMesons.begin(), m_StrangeMesons.end(), std::greater< std::pair<double, int> >());
    std::sort(m_AntiStrangeMesons.begin(), m_AntiStrangeMesons.end(), std::greater< std::pair<double, int> >());
    std::sort(m_ChargeMesons.begin(), m_ChargeMesons.end(), std::greater< std::pair<double, int> >());
    std::sort(m_AntiChargeMesons.begin(), m_AntiChargeMesons.end(), std::greater< std::pair<double, int> >());
    std::sort(m_CharmMesons.begin(), m_CharmMesons.end(), std::greater< std::pair<double, int> >());
    std::sort(m_AntiCharmMesons.begin(), m_AntiCharmMesons.end(), std::greater< std::pair<double, int> >());
    std::sort(m_CharmAll.begin(), m_CharmAll.end(), std::greater< std::pair<double, int> >());
    std::sort(m_AntiCharmAll.begin(), m_AntiCharmAll.end(), std::greater< std::pair<double, int> >());

    for (int i = 1; i < static_cast<int>(m_Baryons.size()); ++i)            m_Baryons[i].first += m_Baryons[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_AntiBaryons.size()); ++i)        m_AntiBaryons[i].first += m_AntiBaryons[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_StrangeMesons.size()); ++i)      m_StrangeMesons[i].first += m_StrangeMesons[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_AntiStrangeMesons.size()); ++i)  m_AntiStrangeMesons[i].first += m_AntiStrangeMesons[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_ChargeMesons.size()); ++i)       m_ChargeMesons[i].first += m_ChargeMesons[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_AntiChargeMesons.size()); ++i)   m_AntiChargeMesons[i].first += m_AntiChargeMesons[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_CharmMesons.size()); ++i)        m_CharmMesons[i].first += m_CharmMesons[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_AntiCharmMesons.size()); ++i)    m_AntiCharmMesons[i].first += m_AntiCharmMesons[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_CharmAll.size()); ++i)           m_CharmAll[i].first += m_CharmAll[i - 1].first;
    for (int i = 1; i < static_cast<int>(m_AntiCharmAll.size()); ++i)       m_AntiCharmAll[i].first += m_AntiCharmAll[i - 1].first;
  }

  std::vector<int> EventGeneratorBase::GenerateTotals() const {
    const_cast<EventGeneratorBase*>(this)->CheckSetParameters();

    if (!m_THM->IsCalculated())
      m_THM->CalculatePrimordialDensities();

    std::vector<int> totals(m_THM->TPS()->Particles().size());

    m_LastWeight = 1.;
    m_LastNormWeight = 1.;

    while (true) {
      // First generate a configuration which satisfies conservation laws
      if (m_Config.fEnsemble == EventGeneratorConfiguration::CCE)
        totals = GenerateTotalsCCE();
      else if (m_Config.fEnsemble == EventGeneratorConfiguration::SCE)
        totals = GenerateTotalsSCE();
      else if (m_Config.fEnsemble == EventGeneratorConfiguration::CE)
        totals = GenerateTotalsCE();
      else
        totals = GenerateTotalsGCE();

      double weight = ComputeWeightNew(totals);
      //std::cout << ComputeWeight(totals) << " " << ComputeWeightNew(totals) << "\n";
      if (weight < 0.)
        continue;

      break;
    }

    fCEAccepted++;

    return totals;
  }

  std::vector<int> EventGeneratorBase::GenerateTotalsGCE() const
  {
    fCETotal++;

    //if (!m_THM->IsGCECalculated()) 
    //  m_THM->CalculateDensitiesGCE();
    if (!m_THM->IsCalculated())
      m_THM->CalculatePrimordialDensities();
    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    const std::vector<double>& densities = m_THM->Densities();
    for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
      double mean = densities[i] * m_THM->Volume();
      int total = RandomGenerators::RandomPoisson(mean);
      totals[i] = total;
    }

    return totals;
  }

  std::vector<int> EventGeneratorBase::GenerateTotalsSCE() const
  {
    if (!m_THM->IsCalculated())
      m_THM->CalculatePrimordialDensities();

    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    // Generate SCE configuration depending on whether Vc > V, Vc = V, or Vc < V
    while (true) {
      const std::vector<double>& densities = m_THM->Densities();

      // If Vc = V then just generate yields from a single ensemble
      if (m_THM->Volume() == m_THM->CanonicalVolume()) {
        totals = GenerateTotalsSCESubVolume(m_THM->Volume());
      }
      // If Vc > V then generate yields from a single ensemble with volume Vc
      // particle from this ensemble is within a smaller volume V
      // with probability V/Vc
      else if (m_THM->Volume() < m_THM->CanonicalVolume()) {
        std::vector<int> totalsaux = GenerateTotalsSCESubVolume(m_THM->CanonicalVolume());
        double prob = m_THM->Volume() / m_THM->CanonicalVolume();
        for (size_t i = 0; i < totalsaux.size(); ++i) {
          for (int j = 0; j < totalsaux[i]; ++j) {
            if (RandomGenerators::randgenMT.rand() < prob)
              totals[i]++;
          }
        }
      }
      // If V > Vc then generate yields from (int)(V/Vc) SCE ensembles
      // plus special treatment of one more ensemble if (V mod Vc) != 0
      else {
        // Number of normal SCE sub-ensembles
        int multiples = static_cast<int>(m_THM->Volume() / m_THM->CanonicalVolume());

        for (int iter = 0; iter < multiples; ++iter) {
          std::vector<int> totalsaux = GenerateTotalsSCESubVolume(m_THM->CanonicalVolume());
          for (size_t i = 0; i < totalsaux.size(); ++i)
            totals[i] += totalsaux[i];
        }

        // For (V mod Vc) != 0 there is one more subvolume Vsub < Vc
        double fraction = (m_THM->Volume() - multiples * m_THM->CanonicalVolume()) / m_THM->CanonicalVolume();

        if (fraction > 0.0) {

          bool successflag = false;

          while (!successflag) {

            std::vector<int> totalsaux = GenerateTotalsSCESubVolume(m_THM->CanonicalVolume());
            std::vector<int> totalsaux2(m_THM->TPS()->Particles().size(), 0);
            int netS = 0;
            for (size_t i = 0; i < totalsaux.size(); ++i) {
              if (m_THM->TPS()->Particles()[i].Strangeness() > 0) {
                for (int j = 0; j < totalsaux[i]; ++j) {
                  if (RandomGenerators::randgenMT.rand() < fraction) {
                    totalsaux2[i]++;
                    netS += m_THM->TPS()->Particles()[i].Strangeness();
                  }
                }
              }
            }

            // Check if possible to match S+ with S-
            // If S = 1 then must be at least one S = -1 hadron
            if (netS == 1) {
              bool fl = false;
              for (size_t i = 0; i < totalsaux.size(); ++i) {
                if (m_THM->TPS()->Particles()[i].Strangeness() < 0 && totalsaux[i] > 0) {
                  fl = true;
                  break;
                }
              }
              if (!fl) {
                printf("**WARNING** SCE Event generator: Cannot match S- with S+ = 1, discarding subconfiguration...\n");
                continue;
              }
            }

            int repeatmax = 10000;
            int repeat = 0;
            while (true) {
              int netS2 = netS;
              for (size_t i = 0; i < totalsaux.size(); ++i) {
                if (m_THM->TPS()->Particles()[i].Strangeness() < 0) {
                  totalsaux2[i] = 0;
                  for (int j = 0; j < totalsaux[i]; ++j) {
                    if (RandomGenerators::randgenMT.rand() < fraction) {
                      totalsaux2[i]++;
                      netS2 += m_THM->TPS()->Particles()[i].Strangeness();
                    }
                  }
                }
              }
              if (netS2 == 0) {
                for (size_t i = 0; i < totalsaux2.size(); ++i)
                  totals[i] += totalsaux2[i];
                successflag = true;
                break;
              }
              repeat++;
              if (repeat >= repeatmax) {
                printf("**WARNING** SCE event generator: Cannot match S- with S+, too many tries, discarding configuration...\n");
                break;
              }
            }
          }
        }
      }

      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        if (m_THM->TPS()->Particles()[i].Strangeness() == 0) {
          double mean = densities[i] * m_THM->Volume();
          int total = RandomGenerators::RandomPoisson(mean);
          totals[i] = total;
        }
      }

      return totals;
    }

    return totals;
  }

  std::vector<int> EventGeneratorBase::GenerateTotalsSCESubVolume(double VolumeSC) const
  {
    if (!m_THM->IsCalculated())
      m_THM->CalculatePrimordialDensities();
    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    std::vector< std::pair<double, int> > fStrangeMesonsc = m_StrangeMesons;
    std::vector< std::pair<double, int> > fAntiStrangeMesonsc = m_AntiStrangeMesons;

    for (size_t i = 0; i < fStrangeMesonsc.size(); ++i)
      fStrangeMesonsc[i].first *= VolumeSC / m_THM->Volume();

    for (size_t i = 0; i < fAntiStrangeMesonsc.size(); ++i)
      fAntiStrangeMesonsc[i].first *= VolumeSC / m_THM->Volume();

    double fMeanSMc = m_MeanSM * VolumeSC / m_THM->Volume();
    double fMeanASMc = m_MeanASM * VolumeSC / m_THM->Volume();

    while (1) {
      fCETotal++;
      const std::vector<double>& densities = m_THM->Densities();

      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) totals[i] = 0;
      int netS = 0;
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        if (m_THM->TPS()->Particles()[i].BaryonCharge() != 0 && m_THM->TPS()->Particles()[i].Strangeness() != 0) {
          double mean = densities[i] * VolumeSC;
          int total = RandomGenerators::RandomPoisson(mean);
          totals[i] = total;
          netS += totals[i] * m_THM->TPS()->Particles()[i].Strangeness();
        }
      }
      int tSM = RandomGenerators::RandomPoisson(fMeanSMc);
      int tASM = RandomGenerators::RandomPoisson(fMeanASMc);

      if (netS != tASM - tSM) continue;

      for (int i = 0; i < tSM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fStrangeMesonsc.begin(), fStrangeMesonsc.end(), std::make_pair(fMeanSMc*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fStrangeMesonsc.size())) tind = fStrangeMesonsc.size() - 1;
        totals[fStrangeMesonsc[tind].second]++;
      }
      for (int i = 0; i < tASM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiStrangeMesonsc.begin(), fAntiStrangeMesonsc.end(), std::make_pair(fMeanASMc*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fAntiStrangeMesonsc.size())) tind = fAntiStrangeMesonsc.size() - 1;
        totals[fAntiStrangeMesonsc[tind].second]++;
      }

      // Cross-check that all resulting strangeness is zero
      int finS = 0;
      for (size_t i = 0; i < totals.size(); ++i) {
        finS += totals[i] * m_THM->TPS()->Particles()[i].Strangeness();
      }

      if (finS != 0) {
        printf("**ERROR** EventGeneratorBase::GenerateTotalsSCESubVolume(): Generated strangeness is non-zero!");
        exit(1);
      }

      return totals;
    }
    return totals;
  }

  std::vector<int> EventGeneratorBase::GenerateTotalsCCE() const
  {
    if (!m_THM->IsCalculated())
      m_THM->CalculatePrimordialDensities();

    // Check there are no multi-charmed particles, otherwise error
    for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
      if (abs(m_THM->TPS()->Particles()[i].Charm()) > 1) {
        printf("**ERROR** CCE Event generator does not support lists with multi-charmed particles\n");
        exit(1);
      }
    }

    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    // Generate CCE configuration depending on whether Vc > V, Vc = V, or Vc < V
    const std::vector<double>& densities = m_THM->Densities();

    // If Vc = V then just generate yields from a single ensemble
    if (m_THM->Volume() == m_THM->CanonicalVolume()) {
      totals = GenerateTotalsCCESubVolume(m_THM->Volume());
    }
    // If Vc > V then generate yields from a single ensemble with volume Vc
    // particle from this ensemble is within a smaller volume V
    // with probability V/Vc
    else if (m_THM->Volume() < m_THM->CanonicalVolume()) {
      std::vector<int> totalsaux = GenerateTotalsCCESubVolume(m_THM->CanonicalVolume());
      double prob = m_THM->Volume() / m_THM->CanonicalVolume();
      for (size_t i = 0; i < totalsaux.size(); ++i) {
        for (int j = 0; j < totalsaux[i]; ++j) {
          if (RandomGenerators::randgenMT.rand() < prob)
            totals[i]++;
        }
      }
    }
    // If V > Vc then generate yields from (int)(V/Vc) SCE ensembles
    // plus special treatment of one more ensemble if (V mod Vc) != 0
    else {
      // Number of normal CCE sub-ensembles
      int multiples = static_cast<int>(m_THM->Volume() / m_THM->CanonicalVolume());

      for (int iter = 0; iter < multiples; ++iter) {
        std::vector<int> totalsaux = GenerateTotalsCCESubVolume(m_THM->CanonicalVolume());
        for (size_t i = 0; i < totalsaux.size(); ++i)
          totals[i] += totalsaux[i];
      }

      // For (V mod Vc) != 0 there is one more subvolume Vsub < Vc
      double fraction = (m_THM->Volume() - multiples * m_THM->CanonicalVolume()) / m_THM->CanonicalVolume();

      if (fraction > 0.0) {

        bool successflag = false;

        while (!successflag) {

          std::vector<int> totalsaux = GenerateTotalsCCESubVolume(m_THM->CanonicalVolume());
          std::vector<int> totalsaux2(m_THM->TPS()->Particles().size(), 0);
          int netC = 0;
          for (size_t i = 0; i < totalsaux.size(); ++i) {
            if (m_THM->TPS()->Particles()[i].Charm() > 0) {
              for (int j = 0; j < totalsaux[i]; ++j) {
                if (RandomGenerators::randgenMT.rand() < fraction) {
                  totalsaux2[i]++;
                  netC += m_THM->TPS()->Particles()[i].Charm();
                }
              }
            }
          }

          // Check if possible to match C+ with C-
          // If C = 1 then must be at least one C = -1 hadron
          if (netC == 1) {
            bool fl = false;
            for (size_t i = 0; i < totalsaux.size(); ++i) {
              if (m_THM->TPS()->Particles()[i].Charm() < 0 && totalsaux[i] > 0) {
                fl = true;
                break;
              }
            }
            if (!fl) {
              printf("**WARNING** CCE Event generator: Cannot match C- with C+ = 1, discarding...\n");
              continue;
            }
          }

          int repeatmax = 10000;
          int repeat = 0;
          while (true) {
            int netC2 = netC;
            for (size_t i = 0; i < totalsaux.size(); ++i) {
              if (m_THM->TPS()->Particles()[i].Charm() < 0) {
                totalsaux2[i] = 0;
                for (int j = 0; j < totalsaux[i]; ++j) {
                  if (RandomGenerators::randgenMT.rand() < fraction) {
                    totalsaux2[i]++;
                    netC2 += m_THM->TPS()->Particles()[i].Charm();
                  }
                }
              }
            }
            if (netC2 == 0) {
              for (size_t i = 0; i < totalsaux2.size(); ++i)
                totals[i] += totalsaux2[i];
              successflag = true;
              break;
            }
            repeat++;
            if (repeat >= repeatmax) {
              printf("**WARNING** CCE event generator: Cannot match S- with S+, too many tries, discarding...\n");
              break;
            }
          }
        }
      }
    }

    for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
      if (m_THM->TPS()->Particles()[i].Charm() == 0) {
        double mean = densities[i] * m_THM->Volume();
        int total = RandomGenerators::RandomPoisson(mean);
        totals[i] = total;
      }
    }

    return totals;
  }

  std::vector<int> EventGeneratorBase::GenerateTotalsCCESubVolume(double VolumeSC) const
  {
    if (!m_THM->IsCalculated())
      m_THM->CalculatePrimordialDensities();

    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    std::vector< std::pair<double, int> > fCharmAllc = m_CharmAll;
    std::vector< std::pair<double, int> > fAntiCharmAllc = m_AntiCharmAll;

    for (size_t i = 0; i < fCharmAllc.size(); ++i)
      fCharmAllc[i].first *= VolumeSC / m_THM->Volume();

    for (size_t i = 0; i < fAntiCharmAllc.size(); ++i)
      fAntiCharmAllc[i].first *= VolumeSC / m_THM->Volume();

    // Assuming no multi-charmed particles
    double fMeanCharmc = m_MeanCHRM * VolumeSC / m_THM->Volume();
    double fMeanAntiCharmc = m_MeanACHRM * VolumeSC / m_THM->Volume();

    fCETotal++;

    int netC = 0;
    int tC = RandomGenerators::RandomPoisson(fMeanCharmc);
    int tAC = RandomGenerators::RandomPoisson(fMeanAntiCharmc);
    while (tC - tAC != m_Config.C - netC) {
      fCETotal++;
      tC = RandomGenerators::RandomPoisson(fMeanCharmc);
      tAC = RandomGenerators::RandomPoisson(fMeanAntiCharmc);
    }

    for (int i = 0; i < tC; ++i) {
      std::vector< std::pair<double, int> >::iterator it = lower_bound(fCharmAllc.begin(), fCharmAllc.end(), std::make_pair(fMeanCharmc*RandomGenerators::randgenMT.rand(), 0));
      int tind = std::distance(fCharmAllc.begin(), it);
      if (tind < 0) tind = 0;
      if (tind >= static_cast<int>(fCharmAllc.size())) tind = fCharmAllc.size() - 1;
      totals[fCharmAllc[tind].second]++;
    }
    for (int i = 0; i < tAC; ++i) {
      std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiCharmAllc.begin(), fAntiCharmAllc.end(), std::make_pair(fMeanAntiCharmc*RandomGenerators::randgenMT.rand(), 0));
      int tind = std::distance(fAntiCharmAllc.begin(), it);
      if (tind < 0) tind = 0;
      if (tind >= static_cast<int>(fAntiCharmAllc.size())) tind = fAntiCharmAllc.size() - 1;
      totals[fAntiCharmAllc[tind].second]++;
    }

    // Cross-check that total resulting net charm is zero
    int finC = 0;
    for (size_t i = 0; i < totals.size(); ++i) {
      finC += totals[i] * m_THM->TPS()->Particles()[i].Charm();
    }

    if (finC != m_Config.C) {
      printf("**ERROR** EventGeneratorBase::GenerateTotalsCCESubVolume(): Generated charm is non-zero!");
      exit(1);
    }

    return totals;
  }

  void EventGeneratorBase::SetParameters()
  {
    SetMomentumGenerators();
    m_ParametersSet = true;
  }

  std::vector<std::vector<double>> EventGeneratorBase::ComputeEVRadii() const
  {
    int Nspecies = m_THM->TPS()->ComponentsNumber();
    std::vector<std::vector<double>> radii = std::vector<std::vector<double>>(Nspecies,
      std::vector<double>(Nspecies, 0.0));

    for (int i = 0; i < Nspecies; ++i) {
      for (int j = 0; j < Nspecies; ++j) {
        double b = 0.0;
        if (i < m_Config.bij.size() && j < m_Config.bij[i].size()) 
          b = 0.5 * (m_Config.bij[i][j] + m_Config.bij[j][i]);
        radii[i][j] = CuteHRGHelper::rv(b);
      }
    }

    return radii;
  }

  bool EventGeneratorBase::CheckEVOverlap(
    const std::vector<SimpleParticle>& particles,
    const SimpleParticle& cand,
    const std::vector<int>& ids,
    const std::vector<std::vector<double>>& radii)
   const
  {
    if (m_Config.fModelType != EventGeneratorConfiguration::DiagonalEV 
      && m_Config.fModelType != EventGeneratorConfiguration::CrosstermsEV
      && m_Config.fModelType != EventGeneratorConfiguration::QvdW)
      return false;
    int idcand = m_THM->TPS()->PdgToId(cand.PDGID);
    for (int ipart = 0; ipart < particles.size(); ++ipart) {
      const SimpleParticle& part = particles[ipart];
      int idpart = 0;
      if (ids.size())
        idpart = ids[ipart];
      else
        idpart = m_THM->TPS()->PdgToId(part.PDGID);

      double r = 0.0;
      if (radii.size()) {
        r = radii[idpart][idcand];
      }
      else {
        double b = 0.5 * (m_Config.bij[idpart][idcand] + m_Config.bij[idcand][idpart]);
        r = CuteHRGHelper::rv(b);
      }

      if (r == 0.0)
        continue;

      double dist2 = ParticleDecaysMC::ParticleDistanceSquared(part, cand);
      if (dist2 <= 4. * r * r)
        return true;
    }
    return false;
  }


  std::vector<int> EventGeneratorBase::GenerateTotalsCE() const {
    if (!m_THM->IsCalculated())
      m_THM->CalculatePrimordialDensities();
    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    const std::vector< std::pair<double, int> >& fBaryonsc = m_Baryons;
    const std::vector< std::pair<double, int> >& fAntiBaryonsc = m_AntiBaryons;
    const std::vector< std::pair<double, int> >& fStrangeMesonsc = m_StrangeMesons;
    const std::vector< std::pair<double, int> >& fAntiStrangeMesonsc = m_AntiStrangeMesons;
    const std::vector< std::pair<double, int> >& fChargeMesonsc = m_ChargeMesons;
    const std::vector< std::pair<double, int> >& fAntiChargeMesonsc = m_AntiChargeMesons;
    const std::vector< std::pair<double, int> >& fCharmMesonsc = m_CharmMesons;
    const std::vector< std::pair<double, int> >& fAntiCharmMesonsc = m_AntiCharmMesons;

    // Primitive rejection sampling (not used, but can be explored for comparisons)
    while (0) {
      fCETotal++;
      int netB = 0, netS = 0, netQ = 0, netC = 0;
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        double mean = m_THM->Densities()[i] * m_THM->Volume();
        int total = RandomGenerators::RandomPoisson(mean);
        totals[i] = total;
        netB += totals[i] * m_THM->TPS()->Particles()[i].BaryonCharge();
        netS += totals[i] * m_THM->TPS()->Particles()[i].Strangeness();
        netQ += totals[i] * m_THM->TPS()->Particles()[i].ElectricCharge();
        netC += totals[i] * m_THM->TPS()->Particles()[i].Charm();
      }
      //if ((!m_Config.CanonicalB || netB == m_THM->Parameters().B) 
      //  && (!m_Config.CanonicalS || netS == m_THM->Parameters().S)
      //  && (!m_Config.CanonicalQ || netQ == m_THM->Parameters().Q)
      //  && (!m_Config.CanonicalC || netC == m_THM->Parameters().C)) {
      //  fCEAccepted++;
      //  return totals;
      //}
      if ((!m_Config.CanonicalB || netB == m_Config.B)
        && (!m_Config.CanonicalS || netS == m_Config.S)
        && (!m_Config.CanonicalQ || netQ == m_Config.Q)
        && (!m_Config.CanonicalC || netC == m_Config.C)) {
        fCEAccepted++;
        return totals;
      }
    }

    // Multi-step procedure as described in F. Becattini, L. Ferroni, hep-ph/0307061
    while (1) {
      fCETotal++;

      const std::vector<double>& densities = m_THM->Densities();

      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) totals[i] = 0;
      int netB = 0, netS = 0, netQ = 0, netC = 0;


      // Light nuclei first
      bool flNuclei = false; // Whether light nuclei appear at all

      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        if (abs(m_THM->TPS()->Particles()[i].BaryonCharge()) > 1) {
          flNuclei = true;
          double mean = densities[i] * m_THM->Volume();
          int total = RandomGenerators::RandomPoisson(mean);
          totals[i] = total;
          netB += totals[i] * m_THM->TPS()->Particles()[i].BaryonCharge();
          netS += totals[i] * m_THM->TPS()->Particles()[i].Strangeness();
          netQ += totals[i] * m_THM->TPS()->Particles()[i].ElectricCharge();
          netC += totals[i] * m_THM->TPS()->Particles()[i].Charm();
        }
      }

      // Then all hadrons

      int tB = 0, tAB = 0;
      // First total baryons and antibaryons from the Poisson distribution
      if (flNuclei || !m_Config.CanonicalB) {
        tB = RandomGenerators::RandomPoisson(m_MeanB);
        tAB = RandomGenerators::RandomPoisson(m_MeanAB);
        if (m_Config.CanonicalB && tB - tAB != m_Config.B - netB) continue;
        //if (RandomGenerators::randgenMT.rand() > RandomGenerators::SkellamProbability(m_THM->Parameters().B - netB, m_MeanB, m_MeanAB))
        //  continue;
      }
      else
      // Generate from the Bessel distribution, using Devroye's method, if no light nuclei
      {
        int nu = m_Config.B - netB;
        if (nu < 0) nu = -nu;
        double a = 2. * sqrt(m_MeanB * m_MeanAB);
        //int BessN = RandomGenerators::BesselDistributionGenerator::RandomBesselDevroye3(a, nu);
        //int BessN = RandomGenerators::BesselDistributionGenerator::RandomBesselPoisson(a, nu);
        //int BessN = RandomGenerators::BesselDistributionGenerator::RandomBesselCombined(a, nu);
        int BessN = RandomGenerators::BesselDistributionGenerator::RandomBesselDevroye1(a, nu);
        if (m_Config.B - netB < 0) {
          tB = BessN;
          tAB = nu + tB;
        }
        else {
          tAB = BessN;
          tB = nu + tAB;
        }
      }

      // Then individual baryons and antibaryons from the multinomial distribution
      for (int i = 0; i < tB; ++i) {
        std::vector< std::pair<double, int> >::const_iterator it = lower_bound(fBaryonsc.begin(), fBaryonsc.end(), std::make_pair(m_MeanB*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fBaryonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fBaryonsc.size())) tind = fBaryonsc.size() - 1;
        totals[fBaryonsc[tind].second]++;
        netS += m_THM->TPS()->Particles()[fBaryonsc[tind].second].Strangeness();
        netQ += m_THM->TPS()->Particles()[fBaryonsc[tind].second].ElectricCharge();
        netC += m_THM->TPS()->Particles()[fBaryonsc[tind].second].Charm();
      }
      for (int i = 0; i < tAB; ++i) {
        std::vector< std::pair<double, int> >::const_iterator it = lower_bound(fAntiBaryonsc.begin(), fAntiBaryonsc.end(), std::make_pair(m_MeanAB*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiBaryonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fAntiBaryonsc.size())) tind = fAntiBaryonsc.size() - 1;
        totals[fAntiBaryonsc[tind].second]++;
        netS += m_THM->TPS()->Particles()[fAntiBaryonsc[tind].second].Strangeness();
        netQ += m_THM->TPS()->Particles()[fAntiBaryonsc[tind].second].ElectricCharge();
        netC += m_THM->TPS()->Particles()[fAntiBaryonsc[tind].second].Charm();
      }

      // Total numbers of (anti)strange mesons
      
      int tSM = RandomGenerators::RandomPoisson(m_MeanSM);
      int tASM = RandomGenerators::RandomPoisson(m_MeanASM);
      if (m_Config.CanonicalS && netS != tASM - tSM + m_Config.S) continue;


      // Multinomial distribution for individual numbers of (anti)strange mesons
      for (int i = 0; i < tSM; ++i) {
        std::vector< std::pair<double, int> >::const_iterator it = lower_bound(fStrangeMesonsc.begin(), fStrangeMesonsc.end(), std::make_pair(m_MeanSM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fStrangeMesonsc.size())) tind = fStrangeMesonsc.size() - 1;
        totals[fStrangeMesonsc[tind].second]++;
        netQ += m_THM->TPS()->Particles()[fStrangeMesonsc[tind].second].ElectricCharge();
        netC += m_THM->TPS()->Particles()[fStrangeMesonsc[tind].second].Charm();
      }
      for (int i = 0; i < tASM; ++i) {
        std::vector< std::pair<double, int> >::const_iterator it = lower_bound(fAntiStrangeMesonsc.begin(), fAntiStrangeMesonsc.end(), std::make_pair(m_MeanASM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fAntiStrangeMesonsc.size())) tind = fAntiStrangeMesonsc.size() - 1;
        totals[fAntiStrangeMesonsc[tind].second]++;
        netQ += m_THM->TPS()->Particles()[fAntiStrangeMesonsc[tind].second].ElectricCharge();
        netC += m_THM->TPS()->Particles()[fAntiStrangeMesonsc[tind].second].Charm();
      }

      // Total numbers of remaining electrically charged mesons
      int tCM = RandomGenerators::RandomPoisson(m_MeanCM);
      int tACM = RandomGenerators::RandomPoisson(m_MeanACM);
      if (m_Config.CanonicalQ && netQ != tACM - tCM + m_Config.Q) continue;

      // Multinomial distribution for individual numbers of remaining electrically charged mesons
      for (int i = 0; i < tCM; ++i) {
        std::vector< std::pair<double, int> >::const_iterator it = lower_bound(fChargeMesonsc.begin(), fChargeMesonsc.end(), std::make_pair(m_MeanCM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fChargeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fChargeMesonsc.size())) tind = fChargeMesonsc.size() - 1;
        totals[fChargeMesonsc[tind].second]++;
        netC += m_THM->TPS()->Particles()[fChargeMesonsc[tind].second].Charm();
      }
      for (int i = 0; i < tACM; ++i) {
        std::vector< std::pair<double, int> >::const_iterator it = lower_bound(fAntiChargeMesonsc.begin(), fAntiChargeMesonsc.end(), std::make_pair(m_MeanACM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiChargeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fAntiChargeMesonsc.size())) tind = fAntiChargeMesonsc.size() - 1;
        totals[fAntiChargeMesonsc[tind].second]++;
        netC += m_THM->TPS()->Particles()[fAntiChargeMesonsc[tind].second].Charm();
      }

      // Total numbers of remaining charmed mesons
      int tCHRMM = RandomGenerators::RandomPoisson(m_MeanCHRMM);
      int tACHRNMM = RandomGenerators::RandomPoisson(m_MeanACHRMM);

      if (m_Config.CanonicalC && netC != tACHRNMM - tCHRMM + m_Config.C) continue;

      // Multinomial distribution for individual numbers of the remaining charmed mesons
      for (int i = 0; i < tCHRMM; ++i) {
        std::vector< std::pair<double, int> >::const_iterator it = lower_bound(fCharmMesonsc.begin(), fCharmMesonsc.end(), std::make_pair(m_MeanCHRMM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fCharmMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fCharmMesonsc.size())) tind = fCharmMesonsc.size() - 1;
        totals[fCharmMesonsc[tind].second]++;
      }
      for (int i = 0; i < tACHRNMM; ++i) {
        std::vector< std::pair<double, int> >::const_iterator it = lower_bound(fAntiCharmMesonsc.begin(), fAntiCharmMesonsc.end(), std::make_pair(m_MeanACHRMM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiCharmMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= static_cast<int>(fAntiCharmMesonsc.size())) tind = fAntiCharmMesonsc.size() - 1;
        totals[fAntiCharmMesonsc[tind].second]++;
      }

      // Poisson distribution for all neutral particles
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 
          && m_THM->TPS()->Particles()[i].Strangeness() == 0 
          && m_THM->TPS()->Particles()[i].ElectricCharge() == 0
          && m_THM->TPS()->Particles()[i].Charm() == 0) {
          double mean = densities[i] * m_THM->Volume();
          int total = RandomGenerators::RandomPoisson(mean);
          totals[i] = total;
        }
      }

      // Cross-check that all resulting charges are OK
      int finB = 0, finQ = 0, finS = 0, finC = 0;
      for (size_t i = 0; i < totals.size(); ++i) {
        finB += totals[i] * m_THM->TPS()->Particles()[i].BaryonCharge();
        finQ += totals[i] * m_THM->TPS()->Particles()[i].ElectricCharge();
        finS += totals[i] * m_THM->TPS()->Particles()[i].Strangeness();
        finC += totals[i] * m_THM->TPS()->Particles()[i].Charm();
      }

      if (m_Config.CanonicalB && finB != m_Config.B) {
        printf("**ERROR** EventGeneratorBase::GenerateTotalsCE(): Generated baryon number does not match the input!");
        exit(1);
      }

      if (m_Config.CanonicalQ && finQ != m_Config.Q) {
        printf("**ERROR** EventGeneratorBase::GenerateTotalsCE(): Generated electric charge does not match the input!");
        exit(1);
      }

      if (m_Config.CanonicalS && finS != m_Config.S) {
        printf("**ERROR** EventGeneratorBase::GenerateTotalsCE(): Generated strangeness does not match the input!");
        exit(1);
      }

      if (m_Config.CanonicalC && finC != m_Config.C) {
        printf("**ERROR** EventGeneratorBase::GenerateTotalsCE(): Generated charm does not match the input!");
        exit(1);
      }

      return totals;
    }
    return totals;
  }

  std::pair<std::vector<int>, double> EventGeneratorBase::SampleYields() const
  {
    std::vector<int> totals = GenerateTotals();
    return make_pair(totals, m_LastNormWeight);
    //std::vector<int> totals = GenerateTotalsGCE();
    //return make_pair(totals, 1.);
  }

  SimpleParticle EventGeneratorBase::SampleParticle(int id) const
  {
    const_cast<EventGeneratorBase*>(this)->CheckSetParameters();
    const ThermalParticle& species = m_THM->TPS()->Particles()[id];
    double tmass = species.Mass();
    if (m_THM->UseWidth() && !species.ZeroWidthEnforced() && !(species.GetResonanceWidthIntegrationType() == ThermalParticle::ZeroWidth))
      tmass = m_BWGens[id]->GetRandom();

    // Check for Bose-Einstein condensation
    // Force m = mu if the sampled mass is too small
    double tmu = m_THM->FullIdealChemicalPotential(id);
    if (species.Statistics() == -1 && tmu > tmass) {
      tmass = tmu;
    }

    std::vector<double> momentum = m_MomentumGens[id]->GetMomentum(tmass);

    return SimpleParticle(momentum[0], momentum[1], momentum[2], tmass, species.PdgId(), 0,
      momentum[3], momentum[4], momentum[5], momentum[6]);
  }

  SimpleParticle EventGeneratorBase::SampleParticleByPdg(long long pdgid) const
  {
    int id = m_THM->TPS()->PdgToId(pdgid);
    if (id == -1) {
      printf("**ERROR** EventGeneratorBase::SampleParticleByPdg(): The input pdg code does not exist in the particle list!");
      exit(1);
    }
    return SampleParticle(id);
  }

  SimpleEvent EventGeneratorBase::SampleMomenta(const std::vector<int>& yields) const
  {
    return SampleMomentaWithShuffle(yields);
    //const_cast<EventGeneratorBase*>(this)->CheckSetParameters();

    //SimpleEvent ret;

    //std::vector<int> ids;

    //for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
    //  const ThermalParticle& species = m_THM->TPS()->Particles()[i];
    //  //primParticles[i].resize(0);
    //  int total = yields[i];
    //  for (int part = 0; part < total; ++part) {
    //    SimpleParticle cand = SampleParticle(i);
    //    if (!CheckEVOverlap(ret.Particles, cand, ids, m_Radii)) {
    //      ret.Particles.push_back(cand);
    //      ids.push_back(i);
    //    }
    //    else {
    //      part--;
    //    }
    //  }
    //}

    //ret.AllParticles = ret.Particles;

    //ret.DecayMap.resize(ret.Particles.size());
    //fill(ret.DecayMap.begin(), ret.DecayMap.end(), -1);

    //ret.DecayMapFinal.resize(ret.Particles.size());
    //for (int i = 0; i < ret.DecayMapFinal.size(); ++i)
    //  ret.DecayMapFinal[i] = i;

    //return ret;
  }

  SimpleEvent EventGeneratorBase::SampleMomentaWithShuffle(const std::vector<int>& yields) const
  {
    const_cast<EventGeneratorBase*>(this)->CheckSetParameters();

    SimpleEvent ret;

    std::vector<int> ids;
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
      for (int part = 0; part < yields[i]; ++part)
        ids.push_back(i);
    std::random_shuffle(ids.begin(), ids.end());

    ret.Particles.resize(ids.size());

    bool flOverlap = true;
    while (flOverlap) {
      int sampled = 0;
      while (sampled < ids.size()) {
        flOverlap = false;
        int i = ids[sampled];
        const ThermalParticle& species = m_THM->TPS()->Particles()[i];
        SimpleParticle cand = SampleParticle(i);

        if (m_Config.fUseEVRejectionCoordinates && 
          (m_Config.fModelType == EventGeneratorConfiguration::DiagonalEV
          || m_Config.fModelType == EventGeneratorConfiguration::CrosstermsEV
          || m_Config.fModelType == EventGeneratorConfiguration::QvdW)) {
          for (int i = 0; i < sampled; ++i) {
            double r = m_Radii[ids[i]][ids[sampled]];
            if (r != 0.0) {
              flOverlap |= (ParticleDecaysMC::ParticleDistanceSquared(ret.Particles[i], cand) <= 4. * r * r);
            }
          }

          if (flOverlap) {
            if (EVUseSPR()) {
              continue;
            }
            else {
              break;
            }
          }
        }

        ret.Particles[sampled] = cand;
        sampled++;
      }
    }

    ret.AllParticles = ret.Particles;

    ret.DecayMap.resize(ret.Particles.size());
    fill(ret.DecayMap.begin(), ret.DecayMap.end(), -1);

    ret.DecayMapFinal.resize(ret.Particles.size());
    for (int i = 0; i < ret.DecayMapFinal.size(); ++i)
      ret.DecayMapFinal[i] = i;

    return ret;
  }

  SimpleEvent EventGeneratorBase::GetEvent(bool DoDecays) const
  {
    const_cast<EventGeneratorBase*>(this)->CheckSetParameters();
    
    if (!m_THM->IsCalculated())
      m_THM->CalculatePrimordialDensities();

    std::vector<int> totals = GenerateTotals();
    if ((m_Config.fModelType == EventGeneratorConfiguration::DiagonalEV
      || m_Config.fModelType == EventGeneratorConfiguration::CrosstermsEV)
      && m_Config.fUseEVRejectionMultiplicity) {
      while (RandomGenerators::randgenMT.rand() > m_LastNormWeight) {
        if (m_LastNormWeight > 1.) {
          printf("**WARNING** Event weight %lf > 1 in Monte Carlo rejection sampling!", m_LastNormWeight);
        }
        totals = GenerateTotals();
      }
      m_LastWeight = 1.0;
      m_LastLogWeight = 0.;
      m_LastNormWeight = 1.0;
    }

    SimpleEvent ret = SampleMomenta(totals);
    ret.weight = m_LastWeight;
    ret.logweight = m_LastLogWeight;
    ret.weight = m_LastNormWeight;

    if (DoDecays)
      return PerformDecays(ret, m_THM->TPS());
    else
      return ret;
  }

  // SimpleEvent EventGeneratorBase::PerformDecaysAlternativeWay(const SimpleEvent& evtin, ThermalParticleSystem* TPS)
  // {
  //   SimpleEvent ret;
  //   ret.weight = evtin.weight;
  //   ret.logweight = evtin.logweight;
  //   ret.AllParticles = evtin.Particles;
  //   reverse(ret.AllParticles.begin(), ret.AllParticles.end());
    
  //   for (auto&& part : ret.AllParticles)
  //     part.processed = false;


  //   bool flag_repeat = true;
  //   while (flag_repeat) {
  //     flag_repeat = false;

  //     bool current_stable_flag = false;
  //     long long current_pdgcode = 0;
  //     long long current_tid = -1;
  //     for (int i = ret.AllParticles.size() - 1; i >= 0; --i) {
  //       SimpleParticle& particle = ret.AllParticles[i];

  //       if (particle.processed)
  //         continue;

  //       long long tpdgcode = particle.PDGID;
  //       if (!(tpdgcode == current_pdgcode))
  //       {
  //         current_tid = TPS->PdgToId(tpdgcode);
  //         if (current_tid != -1)
  //           current_stable_flag = TPS->Particle(current_tid).IsStable();
  //         else
  //           current_stable_flag = true;
  //         current_pdgcode = tpdgcode;
  //       }


  //       if (current_stable_flag) {
  //         //SimpleParticle prt = primParticles[i][j];
  //         //double tpt = prt.GetPt();
  //         //double ty = prt.GetY();
  //         //if (static_cast<int>(m_acc.size()) < i || !m_acc[i].init || m_acc[i].getAcceptance(ty + m_ycm, tpt) > RandomGenerators::randgenMT.rand())
  //         //  ret.Particles.push_back(prt);
  //         //primParticles[i][j].processed = true;
  //         ret.Particles.push_back(particle);
  //         particle.processed = true;
  //       }
  //       else {
  //         flag_repeat = true;
  //         double DecParam = RandomGenerators::randgenMT.rand(), tsum = 0.;

  //         std::vector<double> Bratios;
  //         if (particle.MotherPDGID != 0 ||
  //           TPS->ResonanceWidthIntegrationType() != ThermalParticle::eBW) {
  //           Bratios = TPS->Particles()[current_tid].BranchingRatiosM(particle.m, false);
  //         }
  //         else {
  //           Bratios = TPS->Particles()[current_tid].BranchingRatiosM(particle.m, true);
  //         }

  //         int DecayIndex = 0;
  //         for (DecayIndex = 0; DecayIndex < static_cast<int>(Bratios.size()); ++DecayIndex) {
  //           tsum += Bratios[DecayIndex];
  //           if (tsum > DecParam) break;
  //         }
  //         if (DecayIndex < static_cast<int>(TPS->Particles()[current_tid].Decays().size())) {
  //           std::vector<double> masses(0);
  //           std::vector<long long> pdgids(0);
  //           for (size_t di = 0; di < TPS->Particles()[current_tid].Decays()[DecayIndex].mDaughters.size(); di++) {
  //             long long dpdg = TPS->Particles()[current_tid].Decays()[DecayIndex].mDaughters[di];
  //             if (TPS->PdgToId(dpdg) == -1) {
  //               continue;
  //             }
  //             masses.push_back(TPS->ParticleByPDG(dpdg).Mass());
  //             pdgids.push_back(dpdg);
  //           }
  //           std::vector<SimpleParticle> decres = ParticleDecaysMC::ManyBodyDecay(particle, masses, pdgids);
  //           for (size_t ind = 0; ind < decres.size(); ind++) {
  //             decres[ind].processed = false;
  //             if (TPS->PdgToId(decres[ind].PDGID) != -1)
  //               ret.AllParticles.push_back(decres[ind]);
  //           }
  //           ret.AllParticles[i].processed = true;
  //         }
  //         else {
  //           // Decay through unknown branching ratio, presumably radiative, no hadrons, just ignore the decay products
  //           ret.AllParticles[i].processed = true;
  //         }
  //       }
  //     }
  //   }

  //   return ret;
  // }

  SimpleEvent EventGeneratorBase::PerformDecays(const SimpleEvent& evtin, const ThermalParticleSystem* TPS)
  {
    SimpleEvent ret;
    ret.weight = evtin.weight;
    ret.logweight = evtin.logweight;

    ret.AllParticles = evtin.AllParticles;
    ret.DecayMap = evtin.DecayMap;

    std::vector< std::vector<SimpleParticle> > primParticles(TPS->Particles().size(), std::vector<SimpleParticle>(0));
    std::vector< std::vector<int> > AllParticlesMap(TPS->Particles().size(), std::vector<int>(0));
    for(int i = 0; i < evtin.Particles.size(); ++i) {
      const SimpleParticle& particle = evtin.Particles[i];
      long long tid = TPS->PdgToId(particle.PDGID);
      if (tid != -1) {
        primParticles[tid].push_back(particle);
        AllParticlesMap[tid].push_back(i);
      }
    }

    bool flag_repeat = true;
    while (flag_repeat) {
      flag_repeat = false;
      for (int i = static_cast<int>(primParticles.size()) - 1; i >= 0; --i) {
        if (TPS->Particles()[i].IsStable()) {
          for (size_t j = 0; j < primParticles[i].size(); ++j) {
            if (!primParticles[i][j].processed) {
              SimpleParticle prt = primParticles[i][j];
              ret.Particles.push_back(prt);
              primParticles[i][j].processed = true;
              int tid = AllParticlesMap[i][j];
              while (tid >= 0 && tid < ret.DecayMap.size() && ret.DecayMap[tid] != -1)
                tid = ret.DecayMap[tid];
              ret.DecayMapFinal.push_back(tid);
            }
          }
        }
        else {
          for (size_t j = 0; j < primParticles[i].size(); ++j) {
            if (!primParticles[i][j].processed) {
              flag_repeat = true;
              double DecParam = RandomGenerators::randgenMT.rand(), tsum = 0.;

              std::vector<double> Bratios;
              if (primParticles[i][j].MotherPDGID != 0 ||
                TPS->ResonanceWidthIntegrationType() != ThermalParticle::eBW) {
                Bratios = TPS->Particles()[i].BranchingRatiosM(primParticles[i][j].m, false);
              }
              else {
                Bratios = TPS->Particles()[i].BranchingRatiosM(primParticles[i][j].m, true);
              }

              int DecayIndex = 0;
              for (DecayIndex = 0; DecayIndex < static_cast<int>(Bratios.size()); ++DecayIndex) {
                tsum += Bratios[DecayIndex];
                if (tsum > DecParam) break;
              }
              if (DecayIndex < static_cast<int>(TPS->Particles()[i].Decays().size())) {
                std::vector<double> masses(0);
                std::vector<long long> pdgids(0);
                for (size_t di = 0; di < TPS->Particles()[i].Decays()[DecayIndex].mDaughters.size(); di++) {
                  long long dpdg = TPS->Particles()[i].Decays()[DecayIndex].mDaughters[di];
                  if (TPS->PdgToId(dpdg) == -1) {
                    // Try to see if the daughter particle is a photon/lepton
                    if (ExtraParticles::PdgToId(dpdg) == -1)
                      continue;
                    else
                      masses.push_back(ExtraParticles::ParticleByPdg(dpdg).Mass());
                  }
                  else {
                    masses.push_back(TPS->ParticleByPDG(dpdg).Mass());
                  }
                  pdgids.push_back(dpdg);
                }
                std::vector<SimpleParticle> decres = ParticleDecaysMC::ManyBodyDecay(primParticles[i][j], masses, pdgids);
                for (size_t ind = 0; ind < decres.size(); ind++) {
                  decres[ind].processed = false;
                  if (TPS->PdgToId(decres[ind].PDGID) != -1) {
                    int tid = TPS->PdgToId(decres[ind].PDGID);
                    SimpleParticle& dprt = decres[ind];
                    primParticles[tid].push_back(dprt);
                    ret.AllParticles.push_back(dprt);
                    AllParticlesMap[tid].push_back(static_cast<int>(ret.AllParticles.size()) - 1);
                    ret.DecayMap.push_back(AllParticlesMap[i][j]);
                  }
                  else if (ExtraParticles::PdgToId(decres[ind].PDGID) != -1) {
                    SimpleParticle& dprt = decres[ind];
                    ret.AllParticles.push_back(dprt);
                    ret.DecayMap.push_back(AllParticlesMap[i][j]);
                    ret.PhotonsLeptons.push_back(dprt);
                  }
                }
                primParticles[i][j].processed = true;
              }
              else {
                // Decay through unknown branching ratio, presumably radiative, no hadrons, just ignore decay products
                primParticles[i][j].processed = true;
              }
            }
          }
        }
      }
    }

    return ret;
  }

  std::vector<double> EventGeneratorBase::GCEMeanYields() const
  {
    std::vector<double> ret = m_THM->Densities();
    for(size_t i = 0; i < ret.size(); ++i) {
      ret[i] *= m_THM->Volume();
    }
    return ret;
  }

  void EventGeneratorBase::SetVolume(double V)
  {
    if (m_Config.fEnsemble != EventGeneratorConfiguration::GCE)
      RescaleCEMeans(V / m_THM->Volume());
    m_THM->SetVolume(V); 
    m_Config.CFOParameters.V = V; 
    m_THM->SetCanonicalVolume(V); 
    m_Config.CFOParameters.SVc = V; 
  }

  void EventGeneratorBase::RescaleCEMeans(double Vmod)
  {
    m_MeanB      *= Vmod;
    m_MeanAB     *= Vmod;
    m_MeanSM     *= Vmod;
    m_MeanASM    *= Vmod;
    m_MeanCM     *= Vmod;
    m_MeanACM    *= Vmod;
    m_MeanCHRMM  *= Vmod;
    m_MeanACHRMM *= Vmod;
    m_MeanCHRM   *= Vmod;
    m_MeanACHRM  *= Vmod;

    for (int i = 0; i < static_cast<int>(m_Baryons.size()); ++i)            m_Baryons[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_AntiBaryons.size()); ++i)        m_AntiBaryons[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_StrangeMesons.size()); ++i)      m_StrangeMesons[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_AntiStrangeMesons.size()); ++i)  m_AntiStrangeMesons[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_ChargeMesons.size()); ++i)       m_ChargeMesons[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_AntiChargeMesons.size()); ++i)   m_AntiChargeMesons[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_CharmMesons.size()); ++i)        m_CharmMesons[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_AntiCharmMesons.size()); ++i)    m_AntiCharmMesons[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_CharmAll.size()); ++i)           m_CharmAll[i].first *= Vmod;
    for (int i = 0; i < static_cast<int>(m_AntiCharmAll.size()); ++i)       m_AntiCharmAll[i].first *= Vmod;
  }

  double EventGeneratorBase::ComputeWeight(const std::vector<int>& totals) const
  {
    // Compute the normalized weight factor due to EV/vdW interactions
    // If V - bN < 0, returns -1.
    std::vector<double>* densitiesid = NULL;
    std::vector<double> tmpdens;
    const std::vector<double>& densities = m_THM->Densities();
    if (m_THM->InteractionModel() != ThermalModelBase::Ideal) {
      tmpdens = m_DensitiesIdeal;
      densitiesid = &tmpdens;
    }

    double ret = 1.;

    if (m_THM->InteractionModel() == ThermalModelBase::DiagonalEV) {
      ThermalModelEVDiagonal* model = static_cast<ThermalModelEVDiagonal*>(m_THM);
      double V = m_THM->Volume();
      double VVN = m_THM->Volume();

      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i)
        VVN -= model->ExcludedVolume(i) * totals[i];

      if (VVN < 0.)
        return -1.;

      double weight = 1.;
      double logweight = 0.;

      double normweight = 1.;
      double weightev = 1.;
      double VVNev = m_THM->Volume();
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i)
        VVNev -= model->ExcludedVolume(i) * densities[i] * m_THM->Volume();

      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        weight *= pow(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i], totals[i]);
        if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
          logweight += totals[i] * log(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i]);

        weightev *= pow(VVNev / V * densitiesid->operator[](i) / densities[i], densities[i] * V);

        if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
          //normweight *= pow(VVN / V, totals[i]) / pow(VVNev / V, densities[i] * V) * pow(densitiesid->operator[](i) / densities[i], totals[i] - (densities[i] * V));
          normweight *= pow(VVN / VVNev, totals[i]) * pow(VVNev / V, totals[i] - densities[i] * V) * pow(densitiesid->operator[](i) / densities[i], totals[i] - (densities[i] * V));
      }

      m_LastWeight = weight;
      m_LastLogWeight = logweight;
      m_LastNormWeight = normweight;

      ret = normweight;
    }


    if (m_THM->InteractionModel() == ThermalModelBase::CrosstermsEV) {
      ThermalModelEVCrossterms* model = static_cast<ThermalModelEVCrossterms*>(m_THM);
      double V = m_THM->Volume();

      double weight = 1.;
      double logweight = 0.;
      double normweight = 1.;
      double weightev = 1.;
      bool fl = 1;
      int Nspecies = m_THM->TPS()->Particles().size();
      for (size_t i = 0; i < Nspecies; ++i) {
        double VVN = m_THM->Volume();

        for (size_t j = 0; j < Nspecies; ++j)
          VVN -= model->VirialCoefficient(j, i) * totals[j];

        if (VVN < 0.) { fl = false; break; }

        double VVNev = m_THM->Volume();
        for (size_t j = 0; j < Nspecies; ++j)
          VVNev -= model->VirialCoefficient(j, i) * densities[j] * V;

        weight *= pow(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i], totals[i]);
        if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
          logweight += totals[i] * log(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i]);

        weightev *= pow(VVNev / m_THM->Volume() * densitiesid->operator[](i) / densities[i], densities[i] * V);
        if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
          normweight *= pow(VVN / VVNev, totals[i]) * pow(VVNev / V, totals[i] - densities[i] * V) * pow(densitiesid->operator[](i) / densities[i], totals[i] - (densities[i] * V));
      }

      if (!fl)
        return -1.;

      m_LastWeight = weight;
      m_LastLogWeight = logweight;
      m_LastNormWeight = normweight;

      ret = normweight;
    }

    if (m_THM->InteractionModel() == ThermalModelBase::QvdW) {
      ThermalModelVDW* model = static_cast<ThermalModelVDW*>(m_THM);
      double V = m_THM->Volume();

      double weight = 1.;
      double logweight = 0.;
      double normweight = 1.;
      double weightvdw = 1.;
      bool fl = 1;
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        double VVN = m_THM->Volume();

        for (size_t j = 0; j < m_THM->TPS()->Particles().size(); ++j)
          VVN -= model->VirialCoefficient(j, i) * totals[j];

        if (VVN < 0.) { fl = false; break; }

        double VVNev = m_THM->Volume();
        for (size_t j = 0; j < m_THM->TPS()->Particles().size(); ++j)
          VVNev -= model->VirialCoefficient(j, i) * densities[j] * V;

        weight *= pow(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i], totals[i]);
        if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
          logweight += totals[i] * log(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i]);

        for (size_t j = 0; j < m_THM->TPS()->Particles().size(); ++j) {
          double aij = model->AttractionCoefficient(i, j);
          weight *= exp(aij * totals[j] / m_THM->Parameters().T / m_THM->Volume() * totals[i]);
          logweight += totals[i] * aij * totals[j] / m_THM->Parameters().T / m_THM->Volume();
        }

        weightvdw *= pow(VVNev / m_THM->Volume() * densitiesid->operator[](i) / densities[i], densities[i] * V);
        if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
          normweight *= pow(VVN / VVNev, totals[i]) * pow(VVNev / V, totals[i] - densities[i] * V) * pow(densitiesid->operator[](i) / densities[i], totals[i] - (densities[i] * V));

        for (size_t j = 0; j < m_THM->TPS()->Particles().size(); ++j) {
          double aij = model->AttractionCoefficient(i, j);
          weightvdw *= exp(aij * densities[j] / m_THM->Parameters().T * densities[i] * V);
          normweight *= exp(aij * totals[j] / m_THM->Parameters().T / m_THM->Volume() * totals[i] - aij * densities[j] / m_THM->Parameters().T * densities[i] * V);
        }

      }
      if (!fl)
        return -1.;

      m_LastWeight = weight;
      m_LastLogWeight = logweight;
      m_LastNormWeight = normweight;

      ret = normweight;
    }

    return ret;
  }

  double EventGeneratorBase::ComputeWeightNew(const std::vector<int>& totals) const
  {
    // Compute the normlaized weight factor due to EV/vdW interactions
    // If V - bN < 0, returns -1.
    std::vector<double>* densitiesid = NULL;
    std::vector<double> tmpdens;
    const std::vector<double>& densities = m_THM->Densities();
    if (m_THM->InteractionModel() != ThermalModelBase::Ideal) {
      tmpdens = m_DensitiesIdeal;
      densitiesid = &tmpdens;
    }

    double ret = 1.;

    if (m_THM->InteractionModel() == ThermalModelBase::DiagonalEV) {
      ThermalModelEVDiagonal* model = static_cast<ThermalModelEVDiagonal*>(m_THM);
      double V = m_THM->Volume();
      double VVN = m_THM->Volume();

      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i)
        VVN -= model->ExcludedVolume(i) * totals[i];

      if (VVN < 0.)
        return -1.;

      double weight = 1.;
      double logweight = 0.;

      double normweight = 1.;
      double weightev = 1.;
      double VVNev = m_THM->Volume();
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i)
        VVNev -= model->ExcludedVolume(i) * densities[i] * m_THM->Volume();

      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        weight *= pow(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i], totals[i]);
        if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
          logweight += totals[i] * log(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i]);

        weightev *= pow(VVNev / V * densitiesid->operator[](i) / densities[i], densities[i] * V);

        if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
          //normweight *= pow(VVN / V, totals[i]) / pow(VVNev / V, densities[i] * V) * pow(densitiesid->operator[](i) / densities[i], totals[i] - (densities[i] * V));
          normweight *= pow(VVN / VVNev, totals[i]) * pow(VVNev / V, totals[i] - densities[i] * V) * pow(densitiesid->operator[](i) / densities[i], totals[i] - (densities[i] * V));
      }

      m_LastWeight = weight;
      m_LastLogWeight = logweight;
      m_LastNormWeight = normweight;

      ret = normweight;
    }


    if (m_THM->InteractionModel() == ThermalModelBase::CrosstermsEV) {
      ThermalModelEVCrossterms* model = static_cast<ThermalModelEVCrossterms*>(m_THM);
      double V = m_THM->Volume();

      double weight = 1.;
      double logweight = 0.;
      double normweight = 1.;
      double weightev = 1.;
      bool fl = true;
      int Nspecies = m_THM->TPS()->Particles().size();

      int NEVcomp = model->EVComponentIndices().size();
      std::vector<int> Nscomp(NEVcomp, 0);
      std::vector<double> Nevscomp(NEVcomp, 0.);
      std::vector<double> bns(NEVcomp, 0.), bnevs(NEVcomp, 0.), dmuTs(NEVcomp, 0.);
      const std::vector< std::vector<double> >& virial = model->VirialMatrix();

      for (size_t icomp = 0; icomp < NEVcomp; ++icomp) {
        const std::vector<int>& indis = model->EVComponentIndices()[icomp];
        int Nlocal = indis.size();
        for (size_t ilocal = 0; ilocal < Nlocal; ++ilocal) {
          int ip = indis[ilocal];
          Nscomp[icomp] += totals[ip];
          Nevscomp[icomp] += densities[ip] * V;
        }

        if (indis.size()) {
          int i1 = indis[0];

          for (size_t j = 0; j < Nspecies; ++j) {
            //bns[icomp] += model->VirialCoefficient(j, i1) * totals[j] / V;
            //bnevs[icomp] += model->VirialCoefficient(j, i1) * densities[j];
            bns[icomp] += virial[j][i1] * totals[j];// / V;
            bnevs[icomp] += virial[j][i1] * densities[j];
          }
          bns[icomp] /= V;

          if (bns[icomp] > 1.)
            fl = false;

          //dmuTs[icomp] = model->DeltaMu(i1) / model->Parameters().T;
          dmuTs[icomp] = log(densities[i1] / densitiesid->operator[](i1) / (1. - bnevs[icomp]));
        }

        normweight *= pow((1. - bns[icomp]) / (1. - bnevs[icomp]), Nscomp[icomp]) * exp(-dmuTs[icomp] * (Nscomp[icomp] - Nevscomp[icomp]));

      }

      if (!fl)
        return -1.;

      m_LastWeight = normweight;
      m_LastLogWeight = log(normweight);
      m_LastNormWeight = normweight;

      ret = normweight;
    }

    if (m_THM->InteractionModel() == ThermalModelBase::QvdW) {
      ThermalModelVDW* model = static_cast<ThermalModelVDW*>(m_THM);
      double V = m_THM->Volume();

      double weight = 1.;
      double logweight = 0.;
      double normweight = 1.;
      double weightvdw = 1.;
      bool fl = true;

      int Nspecies = m_THM->TPS()->Particles().size();

      int NVDWcomp = model->VDWComponentIndices().size();
      std::vector<int> Nscomp(NVDWcomp, 0);
      std::vector<double> Nvdwscomp(NVDWcomp, 0.);
      std::vector<double> bns(NVDWcomp, 0.), bnvdws(NVDWcomp, 0.), dmuTs(NVDWcomp, 0.), aijs(NVDWcomp, 0.), aijvdws(NVDWcomp, 0.);
      const std::vector< std::vector<double> >& virial = model->VirialMatrix();
      const std::vector< std::vector<double> >& attr   = model->AttractionMatrix();

      for (size_t icomp = 0; icomp < NVDWcomp; ++icomp) {
        const std::vector<int>& indis = model->VDWComponentIndices()[icomp];
        int Nlocal = indis.size();
        for (size_t ilocal = 0; ilocal < Nlocal; ++ilocal) {
          int ip = indis[ilocal];
          Nscomp[icomp] += totals[ip];
          Nvdwscomp[icomp] += densities[ip] * V;
        }

        if (indis.size()) {
          int i1 = indis[0];

          for (size_t j = 0; j < Nspecies; ++j) {
            //bns[icomp] += model->VirialCoefficient(j, i1) * totals[j] / V;
            //bnevs[icomp] += model->VirialCoefficient(j, i1) * densities[j];
            bns[icomp] += virial[j][i1] * totals[j];// / V;
            bnvdws[icomp] += virial[j][i1] * densities[j];

            aijs[icomp]     += attr[i1][j] * totals[j];
            aijvdws[icomp]  += attr[i1][j] * densities[j];
          }
          bns[icomp] /= V;
          aijs[icomp] /= V * model->Parameters().T;
          aijvdws[icomp] /= model->Parameters().T;


          if (bns[icomp] > 1.)
            fl = false;

          //dmuTs[icomp] = model->DeltaMu(i1) / model->Parameters().T;          
          dmuTs[icomp] = log(densities[i1] / densitiesid->operator[](i1) / (1. - bnvdws[icomp]));
        }

        normweight *= pow((1. - bns[icomp]) / (1. - bnvdws[icomp]), Nscomp[icomp]) 
          * exp(-dmuTs[icomp] * (Nscomp[icomp] - Nvdwscomp[icomp]))
          * exp(aijs[icomp] * Nscomp[icomp] - aijvdws[icomp] * Nvdwscomp[icomp]);
      }

      if (!fl)
        return -1.;

      m_LastWeight = normweight;
      m_LastLogWeight = log(normweight);
      m_LastNormWeight = normweight;

      ret = normweight;
    }

    return ret;
  }
  EventGeneratorConfiguration::EventGeneratorConfiguration()
  {
    fEnsemble = GCE;
    fModelType = PointParticle;
    CFOParameters = ThermalModelParameters();
    B = Q = S = C = 0;
    CanonicalB = CanonicalQ = CanonicalS = CanonicalC = true;
    fUsePCE = false;
    fUseEVRejectionMultiplicity = true;
    fUseEVRejectionCoordinates = true;
    fUseEVUseSPRApproximation = true;
  }

} // namespace thermalfist