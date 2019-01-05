/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
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
#include "HRGVDW/ThermalModelVDW.h"
#include "HRGEventGenerator/ParticleDecays.h"

namespace thermalfist {

  int EventGeneratorBase::fCEAccepted, EventGeneratorBase::fCETotal;
  double EventGeneratorBase::m_LastWeight, EventGeneratorBase::m_LastLogWeight;


  EventGeneratorBase::~EventGeneratorBase()
  {
    ClearMomentumGenerators();
  }

  void EventGeneratorBase::ClearMomentumGenerators()
  {
    for (int i = 0; i < m_MomentumGens.size(); ++i) {
      if (m_MomentumGens[i] != NULL)
        delete m_MomentumGens[i];
    }
    m_MomentumGens.resize(0);

    for (int i = 0; i < m_BWGens.size(); ++i) {
      if (m_BWGens[i] != NULL)
        delete m_BWGens[i];
    }
    m_BWGens.resize(0);
  }

  void EventGeneratorBase::SetConfiguration(const ThermalModelParameters& params, EventGeneratorConfiguration::Ensemble ensemble, EventGeneratorConfiguration::ModelType modeltype, ThermalParticleSystem *TPS, ThermalModelBase *original, ThermalModelBase *THMEVVDW)
  {
    ThermalModelParameters params2 = params;

    m_Config.fEnsemble = ensemble;
    m_Config.fModelType = modeltype;
    m_Config.muB = params2.muB;
    m_Config.muS = params2.muS;
    m_Config.muQ = params2.muQ;
    m_Config.muC = params2.muC;
    m_Config.gammaS = params2.gammaS;
    m_Config.gammaq = params2.gammaq;
    m_Config.gammaC = params2.gammaC;
    m_Config.R = 3. * pow(params2.V, 1. / 3.) / 4. / xMath::Pi();
    m_Config.B = params2.B;
    m_Config.S = params2.S;
    m_Config.Q = params2.Q;
    m_Config.C = params2.C;



    if (ensemble == EventGeneratorConfiguration::CCE) {
      m_Config.muC = 0.;
      params2.muC = 0.;
    }
    else if (ensemble == EventGeneratorConfiguration::SCE) {
      m_Config.muS = 0.;
      m_Config.muC = 0.;
      params2.muS = 0.;
      params2.muC = 0.;
    }
    else if (ensemble == EventGeneratorConfiguration::CE) {
      m_Config.muB = 0.;
      m_Config.muS = 0.;
      m_Config.muQ = 0.;
      m_Config.muC = 0.;
      params2.muB = 0.;
      params2.muS = 0.;
      params2.muQ = 0.;
      params2.muC = 0.;
    }

    if (modeltype == EventGeneratorConfiguration::PointParticle) {
      m_THM = new ThermalModelIdeal(TPS, params2);
    }
    else if (modeltype == EventGeneratorConfiguration::DiagonalEV) {
      m_THM = new ThermalModelEVDiagonal(TPS, params2);
    }
    else if (modeltype == EventGeneratorConfiguration::CrosstermsEV) {
      m_THM = new ThermalModelEVCrossterms(TPS, params2);
    }
    else if (modeltype == EventGeneratorConfiguration::QvdW) {
      m_THM = new ThermalModelVDWFull(TPS, params2);
    }

    m_THM->SetUseWidth(original->TPS()->ResonanceWidthIntegrationType());
    m_THM->SetStatistics(original->QuantumStatistics());
    m_THM->SetNormBratio(original->NormBratio());

    m_THM->ConstrainMuB(original->ConstrainMuB());
    m_THM->ConstrainMuQ(original->ConstrainMuQ());
    m_THM->ConstrainMuS(original->ConstrainMuS());
    m_THM->ConstrainMuC(original->ConstrainMuC());
    m_THM->FillChemicalPotentials();

    if (modeltype != EventGeneratorConfiguration::PointParticle) {
      for (int i = 0; i < m_THM->Densities().size(); ++i) {
        for (int j = 0; j < m_THM->Densities().size(); ++j) {
          m_THM->SetVirial(i, j, THMEVVDW->VirialCoefficient(i, j));
          m_THM->SetAttraction(i, j, THMEVVDW->AttractionCoefficient(i, j));
        }
      }
    }

    if (ensemble == EventGeneratorConfiguration::CE) {
      // The following procedure currently not used, but can be considered if SolveChemicalPotentials routine fails
      if (0 && !(m_Config.B == 0 && m_Config.Q == 0 && m_Config.S == 0)) {
        if (m_Config.S == 0 && !(m_Config.Q == 0 || m_Config.B == 0)) {
          double tmuB = 0., tmuQ = 0., tmuS = 0.;

          double QBrat = (double)(m_Config.Q) / m_Config.B;

          double left = 0.000, right = 0.900, center;

          m_THM->SetBaryonChemicalPotential(left);
          m_THM->SetQoverB(QBrat);
          m_THM->FixParameters();
          m_THM->CalculateDensities();
          double valleft = m_THM->CalculateBaryonDensity() * m_THM->Volume() - m_Config.B;

          m_THM->SetBaryonChemicalPotential(right);
          m_THM->SetQoverB(QBrat);
          m_THM->FixParameters();
          m_THM->CalculateDensities();
          double valright = m_THM->CalculateBaryonDensity() * m_THM->Volume() - m_Config.B;

          double valcenter;

          while ((right - left) > 0.00001) {
            center = (left + right) / 2.;
            m_THM->SetBaryonChemicalPotential(center);
            m_THM->SetQoverB(QBrat);
            m_THM->FixParameters();
            m_THM->CalculateDensities();
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

          m_Config.muB = m_THM->Parameters().muB;
          m_Config.muS = m_THM->Parameters().muS;
          m_Config.muQ = m_THM->Parameters().muQ;
          params2.muB = m_THM->Parameters().muB;
          params2.muS = m_THM->Parameters().muS;
          params2.muQ = m_THM->Parameters().muQ;

          m_THM->FillChemicalPotentials();
        }
      }

      m_THM->SolveChemicalPotentials(m_Config.B, m_Config.Q, m_Config.S, m_Config.C,
        m_THM->Parameters().muB, m_THM->Parameters().muQ, m_THM->Parameters().muS, m_THM->Parameters().muC);
      if (m_THM->Parameters().muB != m_THM->Parameters().muB ||
        m_THM->Parameters().muQ != m_THM->Parameters().muQ ||
        m_THM->Parameters().muS != m_THM->Parameters().muS ||
        m_THM->Parameters().muC != m_THM->Parameters().muC) {
        printf("**WARNING** Could not constrain chemical potentials. Setting all to zero...\n");
        m_THM->SetBaryonChemicalPotential(0.);
        m_THM->SetElectricChemicalPotential(0.);
        m_THM->SetStrangenessChemicalPotential(0.);
        m_THM->SetCharmChemicalPotential(0.);
        m_THM->FillChemicalPotentials();
      }
      m_Config.muB = m_THM->Parameters().muB;
      m_Config.muS = m_THM->Parameters().muS;
      m_Config.muQ = m_THM->Parameters().muQ;
      m_Config.muC = m_THM->Parameters().muC;

      std::cout << "Chemical potentials constrained!\n" << "muB = " 
        << m_Config.muB << " muQ = " 
        << m_Config.muQ << " muS = " 
        << m_Config.muS << " muC = " 
        << m_Config.muC << "\n";

      std::cout << "B = " << m_THM->CalculateBaryonDensity() * m_THM->Parameters().V 
        << " Q = " << m_THM->CalculateChargeDensity() * m_THM->Parameters().V 
        << " S = " << m_THM->CalculateStrangenessDensity() * m_THM->Parameters().V 
        << " C = " << m_THM->CalculateCharmDensity() * m_THM->Parameters().V
        << "\n";
    }

    if (ensemble == EventGeneratorConfiguration::SCE) {
      m_THM->SetQoverB(original->QoverB());
      m_THM->ConstrainMuS(true);
      m_THM->ConstrainMuC(true);
      m_THM->ConstrainChemicalPotentials();
      m_THM->FillChemicalPotentials();
      m_Config.muS = m_THM->Parameters().muS;
      m_Config.muC = m_THM->Parameters().muC;
    }

    if (ensemble == EventGeneratorConfiguration::CCE) {
      m_THM->SetQoverB(original->QoverB());
      m_THM->ConstrainMuC(true);
      m_THM->ConstrainChemicalPotentials();
      m_THM->FillChemicalPotentials();
      m_Config.muC = m_THM->Parameters().muC;
    }

    m_THM->CalculateDensitiesGCE();
    m_DensitiesIdeal = m_THM->GetIdealGasDensities();

    if (ensemble != EventGeneratorConfiguration::GCE)
      PrepareMultinomials();
  }

  void EventGeneratorBase::ReadAcceptance(std::string accfolder)
  {
    m_acc.resize(m_THM->TPS()->Particles().size());
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
      std::string filename = accfolder + "pty_acc_" + to_string_fix(m_THM->TPS()->Particles()[i].PdgId()) + ".txt";
      Acceptance::ReadAcceptanceFunction(m_acc[i], filename);
    }
  }

  void EventGeneratorBase::PrepareMultinomials() {
    if (!m_THM->IsGCECalculated())
      m_THM->CalculateDensitiesGCE();

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

    std::vector<double> densities = m_THM->Densities();

    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
      densities[i] *= m_THM->Volume();

    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
      if (m_THM->TPS()->Particles()[i].BaryonCharge() == 1) {
        m_Baryons.push_back(std::make_pair(densities[i], i));
        m_MeanB += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == -1) {
        m_AntiBaryons.push_back(std::make_pair(densities[i], i));
        m_MeanAB += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 1) {
        m_StrangeMesons.push_back(std::make_pair(densities[i], i));
        m_MeanSM += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == -1) {
        m_AntiStrangeMesons.push_back(std::make_pair(densities[i], i));
        m_MeanASM += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == 1) {
        m_ChargeMesons.push_back(std::make_pair(densities[i], i));
        m_MeanCM += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == -1) {
        m_AntiChargeMesons.push_back(std::make_pair(densities[i], i));
        m_MeanACM += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == 1) {
        m_ChargeMesons.push_back(std::make_pair(densities[i], i));
        m_MeanCM += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == 0 && m_THM->TPS()->Particles()[i].Charm() == 1) {
        m_CharmMesons.push_back(std::make_pair(densities[i], i));
        m_MeanCHRMM += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 && m_THM->TPS()->Particles()[i].Strangeness() == 0 && m_THM->TPS()->Particles()[i].ElectricCharge() == 0 && m_THM->TPS()->Particles()[i].Charm() == -1) {
        m_AntiCharmMesons.push_back(std::make_pair(densities[i], i));
        m_MeanACHRMM += densities[i];
      }

      if (m_THM->TPS()->Particles()[i].Charm() == 1) {
        m_CharmAll.push_back(std::make_pair(densities[i], i));
        m_MeanCHRM += densities[i];
      }
      else if (m_THM->TPS()->Particles()[i].Charm() == -1) {
        m_AntiCharmAll.push_back(std::make_pair(densities[i], i));
        m_MeanACHRM += densities[i];
      }
    }

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

    for (int i = 1; i < m_Baryons.size(); ++i)            m_Baryons[i].first += m_Baryons[i - 1].first;
    for (int i = 1; i < m_AntiBaryons.size(); ++i)        m_AntiBaryons[i].first += m_AntiBaryons[i - 1].first;
    for (int i = 1; i < m_StrangeMesons.size(); ++i)      m_StrangeMesons[i].first += m_StrangeMesons[i - 1].first;
    for (int i = 1; i < m_AntiStrangeMesons.size(); ++i)  m_AntiStrangeMesons[i].first += m_AntiStrangeMesons[i - 1].first;
    for (int i = 1; i < m_ChargeMesons.size(); ++i)       m_ChargeMesons[i].first += m_ChargeMesons[i - 1].first;
    for (int i = 1; i < m_AntiChargeMesons.size(); ++i)   m_AntiChargeMesons[i].first += m_AntiChargeMesons[i - 1].first;
    for (int i = 1; i < m_CharmMesons.size(); ++i)        m_CharmMesons[i].first += m_CharmMesons[i - 1].first;
    for (int i = 1; i < m_AntiCharmMesons.size(); ++i)    m_AntiCharmMesons[i].first += m_AntiCharmMesons[i - 1].first;
    for (int i = 1; i < m_CharmAll.size(); ++i)           m_CharmAll[i].first += m_CharmAll[i - 1].first;
    for (int i = 1; i < m_AntiCharmAll.size(); ++i)       m_AntiCharmAll[i].first += m_AntiCharmAll[i - 1].first;
  }

  std::vector<int> EventGeneratorBase::GenerateTotals() const {
    if (!m_THM->IsGCECalculated())
      m_THM->CalculateDensitiesGCE();
    std::vector<int> totals(m_THM->TPS()->Particles().size());

    m_LastWeight = 1.;

    while (true) {
      // First generate a configuration which satisfies conservation laws
      if (m_Config.fEnsemble == EventGeneratorConfiguration::CCE)
        totals = GenerateTotalsCCE();
      else if (m_Config.fEnsemble == EventGeneratorConfiguration::SCE)
        totals = GenerateTotalsSCEnew();
      else if (m_Config.fEnsemble == EventGeneratorConfiguration::CE)
        totals = GenerateTotalsCE();
      else
        totals = GenerateTotalsGCE();

      fCETotal++;

      // Then take into account EV/vdW interactions using importance sampling
      std::vector<double> *densitiesid = NULL;
      std::vector<double> tmpdens;
      const std::vector<double>& densities = m_THM->Densities();
      if (m_THM->InteractionModel() != ThermalModelBase::Ideal) {
        tmpdens = m_DensitiesIdeal;
        densitiesid = &tmpdens;
      }

      if (m_THM->InteractionModel() == ThermalModelBase::DiagonalEV) {
        ThermalModelEVDiagonal *model = static_cast<ThermalModelEVDiagonal*>(m_THM);
        double VVN = m_THM->Volume();

        for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
          VVN -= model->ExcludedVolume(i) * totals[i];

        if (VVN < 0.) continue;
        double weight = 1.;
        double logweight = 0.;
        for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
          weight *= pow(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i], totals[i]);
          if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
            logweight += totals[i] * log(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i]);
        }

        m_LastWeight = weight;
        m_LastLogWeight = logweight;
      }


      if (m_THM->InteractionModel() == ThermalModelBase::CrosstermsEV) {
        ThermalModelEVCrossterms *model = static_cast<ThermalModelEVCrossterms*>(m_THM);
        double weight = 1.;
        double logweight = 0.;
        bool fl = 1;
        for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
          double VVN = m_THM->Volume();

          for (int j = 0; j < m_THM->TPS()->Particles().size(); ++j)
            VVN -= model->VirialCoefficient(j, i) * totals[j];

          if (VVN < 0.) { fl = false; break; }

          weight *= pow(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i], totals[i]);
          if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
            logweight += totals[i] * log(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i]);
        }
        if (!fl) continue;
        m_LastWeight = weight;
        m_LastLogWeight = logweight;
      }

      if (m_THM->InteractionModel() == ThermalModelBase::QvdW) {
        ThermalModelVDW *model = static_cast<ThermalModelVDW*>(m_THM);
        double weight = 1.;
        double logweight = 0.;
        bool fl = 1;
        for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
          double VVN = m_THM->Volume();

          for (int j = 0; j < m_THM->TPS()->Particles().size(); ++j)
            VVN -= model->VirialCoefficient(j, i) * totals[j];

          if (VVN < 0.) { fl = false; break; }

          weight *= pow(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i], totals[i]);
          if (densitiesid->operator[](i) > 0. && densities[i] > 0.)
            logweight += totals[i] * log(VVN / m_THM->Volume() * densitiesid->operator[](i) / densities[i]);

          for (int j = 0; j < m_THM->TPS()->Particles().size(); ++j) {
            double aij = model->AttractionCoefficient(i, j);
            weight *= exp(aij * totals[j] / m_THM->Parameters().T / m_THM->Volume() * totals[i]);
            logweight += totals[i] * aij * totals[j] / m_THM->Parameters().T / m_THM->Volume();
          }

        }
        if (!fl) continue;
        m_LastWeight = weight;
        m_LastLogWeight = logweight;
      }

      break;
    }

    fCEAccepted++;

    return totals;
  }

  std::vector<int> EventGeneratorBase::GenerateTotalsGCE() const
  {
    if (!m_THM->IsGCECalculated()) m_THM->CalculateDensitiesGCE();
    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    const std::vector<double>& densities = m_THM->Densities();
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
      double mean = densities[i] * m_THM->Volume();
      int total = RandomGenerators::RandomPoisson(mean);
      totals[i] = total;
    }

    return totals;
  }

  std::vector<int> EventGeneratorBase::GenerateTotalsSCE() const {
    if (!m_THM->IsGCECalculated())
      m_THM->CalculateDensitiesGCE();

    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    std::vector< std::pair<double, int> > fStrangeMesonsc = m_StrangeMesons;
    std::vector< std::pair<double, int> > fAntiStrangeMesonsc = m_AntiStrangeMesons;

    while (true) {
      const std::vector<double>& densities = m_THM->Densities();

      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) totals[i] = 0;
      int netS = 0;
      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        if (m_THM->TPS()->Particles()[i].BaryonCharge() != 0 && m_THM->TPS()->Particles()[i].Strangeness() != 0) {
          double mean = densities[i] * m_THM->Volume();
          int total = RandomGenerators::RandomPoisson(mean);
          totals[i] = total;
          netS += totals[i] * m_THM->TPS()->Particles()[i].Strangeness();
        }
      }
      int tSM = RandomGenerators::RandomPoisson(m_MeanSM);
      int tASM = RandomGenerators::RandomPoisson(m_MeanASM);

      if (netS != tASM - tSM) continue;

      for (int i = 0; i < tSM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fStrangeMesonsc.begin(), fStrangeMesonsc.end(), std::make_pair(m_MeanSM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fStrangeMesonsc.size()) tind = fStrangeMesonsc.size() - 1;
        totals[fStrangeMesonsc[tind].second]++;
      }
      for (int i = 0; i < tASM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiStrangeMesonsc.begin(), fAntiStrangeMesonsc.end(), std::make_pair(m_MeanASM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fAntiStrangeMesonsc.size()) tind = fAntiStrangeMesonsc.size() - 1;
        totals[fAntiStrangeMesonsc[tind].second]++;
      }

      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
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

  std::vector<int> EventGeneratorBase::GenerateTotalsSCEnew() const
  {
    if (!m_THM->IsGCECalculated())
      m_THM->CalculateDensitiesGCE();

    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    // Generate SCE configuration depending on whether Vc > V, Vc = V, or Vc < V
    while (true) {
      const std::vector<double>& densities = m_THM->Densities();

      // If Vc = V then just generate yields from a single ensemble
      if (m_THM->Volume() == m_THM->StrangenessCanonicalVolume()) {
        totals = GenerateTotalsSCESubVolume(m_THM->Volume());
      }
      // If Vc > V then generate yields from a single ensemble with volume Vc
      // particle from this ensemble is within a smaller volume V
      // with probability V/Vc
      else if (m_THM->Volume() < m_THM->StrangenessCanonicalVolume()) {
        std::vector<int> totalsaux = GenerateTotalsSCESubVolume(m_THM->StrangenessCanonicalVolume());
        double prob = m_THM->Volume() / m_THM->StrangenessCanonicalVolume();
        for (int i = 0; i < totalsaux.size(); ++i) {
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
        int multiples = static_cast<int>(m_THM->Volume() / m_THM->StrangenessCanonicalVolume());

        for (int iter = 0; iter < multiples; ++iter) {
          std::vector<int> totalsaux = GenerateTotalsSCESubVolume(m_THM->StrangenessCanonicalVolume());
          for (int i = 0; i < totalsaux.size(); ++i)
            totals[i] += totalsaux[i];
        }

        // For (V mod Vc) != 0 there is one more subvolume Vsub < Vc
        double fraction = (m_THM->Volume() - multiples * m_THM->StrangenessCanonicalVolume()) / m_THM->StrangenessCanonicalVolume();

        if (fraction > 0.0) {

          bool successflag = false;

          while (!successflag) {

            std::vector<int> totalsaux = GenerateTotalsSCESubVolume(m_THM->StrangenessCanonicalVolume());
            std::vector<int> totalsaux2(m_THM->TPS()->Particles().size(), 0);
            int netS = 0;
            for (int i = 0; i < totalsaux.size(); ++i) {
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
              for (int i = 0; i < totalsaux.size(); ++i) {
                if (m_THM->TPS()->Particles()[i].Strangeness() < 0 && totalsaux[i] > 0) {
                  fl = true;
                  break;
                }
              }
              if (!fl) {
                printf("**WARNING** SCE Event generator: Cannot match S- with S+ = 1, discarding...\n");
                continue;
              }
            }

            int repeatmax = 10000;
            int repeat = 0;
            while (true) {
              int netS2 = netS;
              for (int i = 0; i < totalsaux.size(); ++i) {
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
                for (int i = 0; i < totalsaux2.size(); ++i)
                  totals[i] += totalsaux2[i];
                successflag = true;
                break;
              }
              repeat++;
              if (repeat >= repeatmax) {
                printf("**WARNING** SCE event generator: Cannot match S- with S+, too many tries, discarding...\n");
                break;
              }
            }
          }
        }
      }

      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
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
    if (!m_THM->IsGCECalculated()) m_THM->CalculateDensitiesGCE();
    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    std::vector< std::pair<double, int> > fStrangeMesonsc = m_StrangeMesons;
    std::vector< std::pair<double, int> > fAntiStrangeMesonsc = m_AntiStrangeMesons;

    for (int i = 0; i < fStrangeMesonsc.size(); ++i)
      fStrangeMesonsc[i].first *= VolumeSC / m_THM->Volume();

    for (int i = 0; i < fAntiStrangeMesonsc.size(); ++i)
      fAntiStrangeMesonsc[i].first *= VolumeSC / m_THM->Volume();

    double fMeanSMc = m_MeanSM * VolumeSC / m_THM->Volume();
    double fMeanASMc = m_MeanASM * VolumeSC / m_THM->Volume();

    while (1) {
      const std::vector<double>& densities = m_THM->Densities();

      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) totals[i] = 0;
      int netS = 0;
      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
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
        if (tind >= fStrangeMesonsc.size()) tind = fStrangeMesonsc.size() - 1;
        totals[fStrangeMesonsc[tind].second]++;
      }
      for (int i = 0; i < tASM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiStrangeMesonsc.begin(), fAntiStrangeMesonsc.end(), std::make_pair(fMeanASMc*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fAntiStrangeMesonsc.size()) tind = fAntiStrangeMesonsc.size() - 1;
        totals[fAntiStrangeMesonsc[tind].second]++;
      }
      return totals;
    }
    return totals;
  }

  std::vector<int> EventGeneratorBase::GenerateTotalsCCE() const
  {
    if (!m_THM->IsGCECalculated())
      m_THM->CalculateDensitiesGCE();

    // Check there are no multi-charmed particles, otherwise error
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
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
      for (int i = 0; i < totalsaux.size(); ++i) {
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
        for (int i = 0; i < totalsaux.size(); ++i)
          totals[i] += totalsaux[i];
      }

      // For (V mod Vc) != 0 there is one more subvolume Vsub < Vc
      double fraction = (m_THM->Volume() - multiples * m_THM->CanonicalVolume()) / m_THM->CanonicalVolume();

      if (fraction > 0.0) {

        bool successflag = false;

        while (!successflag) {

          std::vector<int> totalsaux = GenerateTotalsCCESubVolume(m_THM->StrangenessCanonicalVolume());
          std::vector<int> totalsaux2(m_THM->TPS()->Particles().size(), 0);
          int netC = 0;
          for (int i = 0; i < totalsaux.size(); ++i) {
            if (m_THM->TPS()->Particles()[i].Strangeness() > 0) {
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
            for (int i = 0; i < totalsaux.size(); ++i) {
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
            for (int i = 0; i < totalsaux.size(); ++i) {
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
              for (int i = 0; i < totalsaux2.size(); ++i)
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

    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
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
    if (!m_THM->IsGCECalculated()) 
      m_THM->CalculateDensitiesGCE();

    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    std::vector< std::pair<double, int> > fCharmAllc = m_CharmAll;
    std::vector< std::pair<double, int> > fAntiCharmAllc = m_AntiCharmAll;

    for (int i = 0; i < fCharmAllc.size(); ++i)
      fCharmAllc[i].first *= VolumeSC / m_THM->Volume();

    for (int i = 0; i < fAntiCharmAllc.size(); ++i)
      fAntiCharmAllc[i].first *= VolumeSC / m_THM->Volume();

    // Assuming no multi-charmed particles
    double fMeanCharmc = m_MeanCHRM * VolumeSC / m_THM->Volume();
    double fMeanAntiCharmc = m_MeanACHRM * VolumeSC / m_THM->Volume();

    int netC = 0;
    int tC = RandomGenerators::RandomPoisson(fMeanCharmc);
    int tAC = RandomGenerators::RandomPoisson(fMeanAntiCharmc);
    while (tC - tAC != m_THM->Parameters().C - netC) {
      tC = RandomGenerators::RandomPoisson(fMeanCharmc);
      tAC = RandomGenerators::RandomPoisson(fMeanAntiCharmc);
    }

    for (int i = 0; i < tC; ++i) {
      std::vector< std::pair<double, int> >::iterator it = lower_bound(fCharmAllc.begin(), fCharmAllc.end(), std::make_pair(fMeanCharmc*RandomGenerators::randgenMT.rand(), 0));
      int tind = std::distance(fCharmAllc.begin(), it);
      if (tind < 0) tind = 0;
      if (tind >= fCharmAllc.size()) tind = fCharmAllc.size() - 1;
      totals[fCharmAllc[tind].second]++;
    }
    for (int i = 0; i < tAC; ++i) {
      std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiCharmAllc.begin(), fAntiCharmAllc.end(), std::make_pair(fMeanAntiCharmc*RandomGenerators::randgenMT.rand(), 0));
      int tind = std::distance(fAntiCharmAllc.begin(), it);
      if (tind < 0) tind = 0;
      if (tind >= fAntiCharmAllc.size()) tind = fAntiCharmAllc.size() - 1;
      totals[fAntiCharmAllc[tind].second]++;
    }
    return totals;
  }


  std::vector<int> EventGeneratorBase::GenerateTotalsCE() const {
    if (!m_THM->IsGCECalculated()) m_THM->CalculateDensitiesGCE();
    std::vector<int> totals(m_THM->TPS()->Particles().size(), 0);

    std::vector< std::pair<double, int> > fBaryonsc = m_Baryons;
    std::vector< std::pair<double, int> > fAntiBaryonsc = m_AntiBaryons;
    std::vector< std::pair<double, int> > fStrangeMesonsc = m_StrangeMesons;
    std::vector< std::pair<double, int> > fAntiStrangeMesonsc = m_AntiStrangeMesons;
    std::vector< std::pair<double, int> > fChargeMesonsc = m_ChargeMesons;
    std::vector< std::pair<double, int> > fAntiChargeMesonsc = m_AntiChargeMesons;
    std::vector< std::pair<double, int> > fCharmMesonsc = m_CharmMesons;
    std::vector< std::pair<double, int> > fAntiCharmMesonsc = m_AntiCharmMesons;

    // Primitive rejection sampling (not used, but can be explored for comparisons)
    while (0) {
      fCETotal++;
      int netB = 0, netS = 0, netQ = 0, netC = 0;
      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        double mean = m_THM->Densities()[i] * m_THM->Volume();
        int total = RandomGenerators::RandomPoisson(mean);
        totals[i] = total;
        netB += totals[i] * m_THM->TPS()->Particles()[i].BaryonCharge();
        netS += totals[i] * m_THM->TPS()->Particles()[i].Strangeness();
        netQ += totals[i] * m_THM->TPS()->Particles()[i].ElectricCharge();
        netC += totals[i] * m_THM->TPS()->Particles()[i].Charm();
      }
      if (netB == m_THM->Parameters().B 
        && netS == m_THM->Parameters().S 
        && netQ == m_THM->Parameters().Q
        && netC == m_THM->Parameters().C) {
        fCEAccepted++;
        return totals;
      }
    }

    // Multi-step procedure as described in F. Becattini, L. Ferroni, hep-ph/0307061
    while (1) {
      const std::vector<double>& densities = m_THM->Densities();

      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) totals[i] = 0;
      int netB = 0, netS = 0, netQ = 0, netC = 0;


      // Light nuclei first

      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        if (abs(m_THM->TPS()->Particles()[i].BaryonCharge()) > 1) {
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

      // First total baryons and antibaryons from the Poisson distribution
      int tB = RandomGenerators::RandomPoisson(m_MeanB);
      int tAB = RandomGenerators::RandomPoisson(m_MeanAB);
      if (tB - tAB != m_THM->Parameters().B - netB) continue;

      // Then individual baryons and antibaryons from the multinomial distribution
      for (int i = 0; i < tB; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fBaryonsc.begin(), fBaryonsc.end(), std::make_pair(m_MeanB*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fBaryonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fBaryonsc.size()) tind = fBaryonsc.size() - 1;
        totals[fBaryonsc[tind].second]++;
        netS += m_THM->TPS()->Particles()[fBaryonsc[tind].second].Strangeness();
        netQ += m_THM->TPS()->Particles()[fBaryonsc[tind].second].ElectricCharge();
        netC += m_THM->TPS()->Particles()[fBaryonsc[tind].second].Charm();
      }
      for (int i = 0; i < tAB; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiBaryonsc.begin(), fAntiBaryonsc.end(), std::make_pair(m_MeanAB*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiBaryonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fAntiBaryonsc.size()) tind = fAntiBaryonsc.size() - 1;
        totals[fAntiBaryonsc[tind].second]++;
        netS += m_THM->TPS()->Particles()[fAntiBaryonsc[tind].second].Strangeness();
        netQ += m_THM->TPS()->Particles()[fAntiBaryonsc[tind].second].Charm();
      }

      // Total numbers of (anti)strange mesons
      int tSM = RandomGenerators::RandomPoisson(m_MeanSM);
      int tASM = RandomGenerators::RandomPoisson(m_MeanASM);

      if (netS != tASM - tSM + m_THM->Parameters().S) continue;

      // Multinomial distribution for individual numbers of (anti)strange mesons
      for (int i = 0; i < tSM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fStrangeMesonsc.begin(), fStrangeMesonsc.end(), std::make_pair(m_MeanSM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fStrangeMesonsc.size()) tind = fStrangeMesonsc.size() - 1;
        totals[fStrangeMesonsc[tind].second]++;
        netQ += m_THM->TPS()->Particles()[fStrangeMesonsc[tind].second].ElectricCharge();
        netC += m_THM->TPS()->Particles()[fStrangeMesonsc[tind].second].Charm();
      }
      for (int i = 0; i < tASM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiStrangeMesonsc.begin(), fAntiStrangeMesonsc.end(), std::make_pair(m_MeanASM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiStrangeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fAntiStrangeMesonsc.size()) tind = fAntiStrangeMesonsc.size() - 1;
        totals[fAntiStrangeMesonsc[tind].second]++;
        netQ += m_THM->TPS()->Particles()[fAntiStrangeMesonsc[tind].second].ElectricCharge();
        netC += m_THM->TPS()->Particles()[fAntiStrangeMesonsc[tind].second].Charm();
      }

      // Total numbers of remaining electrically charged mesons
      int tCM = RandomGenerators::RandomPoisson(m_MeanCM);
      int tACM = RandomGenerators::RandomPoisson(m_MeanACM);

      if (netQ != tACM - tCM + m_THM->Parameters().Q) continue;

      // Multinomial distribution for individual numbers of remaining electrically charged mesons
      for (int i = 0; i < tCM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fChargeMesonsc.begin(), fChargeMesonsc.end(), std::make_pair(m_MeanCM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fChargeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fChargeMesonsc.size()) tind = fChargeMesonsc.size() - 1;
        totals[fChargeMesonsc[tind].second]++;
        netC += m_THM->TPS()->Particles()[fChargeMesonsc[tind].second].Charm();
      }
      for (int i = 0; i < tACM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiChargeMesonsc.begin(), fAntiChargeMesonsc.end(), std::make_pair(m_MeanACM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiChargeMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fAntiChargeMesonsc.size()) tind = fAntiChargeMesonsc.size() - 1;
        totals[fAntiChargeMesonsc[tind].second]++;
        netC += m_THM->TPS()->Particles()[fAntiChargeMesonsc[tind].second].Charm();
      }

      // Total numbers of remaining charmed mesons
      int tCHRMM = RandomGenerators::RandomPoisson(m_MeanCHRMM);
      int tACHRNMM = RandomGenerators::RandomPoisson(m_MeanACHRMM);

      if (netC != tACHRNMM - tCHRMM + m_THM->Parameters().C) continue;

      // Multinomial distribution for individual numbers of remaining charmed mesons
      for (int i = 0; i < tCHRMM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fCharmMesonsc.begin(), fCharmMesonsc.end(), std::make_pair(m_MeanCHRMM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fCharmMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fCharmMesonsc.size()) tind = fCharmMesonsc.size() - 1;
        totals[fCharmMesonsc[tind].second]++;
      }
      for (int i = 0; i < tACHRNMM; ++i) {
        std::vector< std::pair<double, int> >::iterator it = lower_bound(fAntiCharmMesonsc.begin(), fAntiCharmMesonsc.end(), std::make_pair(m_MeanACHRMM*RandomGenerators::randgenMT.rand(), 0));
        int tind = std::distance(fAntiCharmMesonsc.begin(), it);
        if (tind < 0) tind = 0;
        if (tind >= fAntiCharmMesonsc.size()) tind = fAntiCharmMesonsc.size() - 1;
        totals[fAntiCharmMesonsc[tind].second]++;
      }

      // Poisson distribution for all neutral particles
      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        if (m_THM->TPS()->Particles()[i].BaryonCharge() == 0 
          && m_THM->TPS()->Particles()[i].Strangeness() == 0 
          && m_THM->TPS()->Particles()[i].ElectricCharge() == 0
          && m_THM->TPS()->Particles()[i].Charm() == 0) {
          double mean = densities[i] * m_THM->Volume();
          int total = RandomGenerators::RandomPoisson(mean);
          totals[i] = total;
        }
      }

      return totals;
    }
    return totals;
  }

  SimpleEvent EventGeneratorBase::GetEvent(bool PerformDecays) const
  {
    SimpleEvent ret;
    if (!m_THM->IsGCECalculated()) m_THM->CalculateDensitiesGCE();

    std::vector<int> totals = GenerateTotals();
    ret.weight = m_LastWeight;
    ret.logweight = m_LastLogWeight;

    // Do not perform decay, to be deprecated
    if (m_OnlyStable) 
    {
      for (int i = m_THM->TPS()->Particles().size() - 1; i >= 0; --i) {
        if (m_THM->TPS()->Particles()[i].IsStable()) {
          int total = totals[i];
          for (int part = 0; part < total; ++part) {
            std::vector<double> momentum = m_MomentumGens[i]->GetMomentum();
            SimpleParticle prt = SimpleParticle(momentum[0], momentum[1], momentum[2], m_THM->TPS()->Particles()[i].Mass(), m_THM->TPS()->Particles()[i].PdgId());

            double tpt = prt.GetPt();
            double ty = prt.GetY();
            if (m_acc.size() < i || !m_acc[i].init || m_acc[i].getAcceptance(ty + m_ycm, tpt) > RandomGenerators::randgenMT.rand())
              ret.Particles.push_back(prt);
          }
        }
      }
    }
    // Generate all primordial hadrons and optionally perform decays of resonances
    else 
    {
      std::vector< std::vector<SimpleParticle> > primParticles(m_THM->TPS()->Particles().size());
      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        {
          primParticles[i].resize(0);
          int total = totals[i];
          for (int part = 0; part < total; ++part) {
            std::vector<double> momentum = m_MomentumGens[i]->GetMomentum();

            double tmass = m_THM->TPS()->Particles()[i].Mass();
            if (m_THM->UseWidth() && !(m_THM->TPS()->Particles()[i].ResonanceWidth() / m_THM->TPS()->Particles()[i].Mass() < 0.01))
              tmass = m_BWGens[i]->GetRandom();

            primParticles[i].push_back(SimpleParticle(momentum[0], momentum[1], momentum[2], tmass, m_THM->TPS()->Particles()[i].PdgId()));
          }
        }
      }

      bool flag_repeat = true;
      while (flag_repeat) {
        flag_repeat = false;
        for (int i = primParticles.size() - 1; i >= 0; --i) {
          if (!PerformDecays || m_THM->TPS()->Particles()[i].IsStable()) {
            for (int j = 0; j < primParticles[i].size(); ++j) {
              if (!primParticles[i][j].processed) {
                SimpleParticle prt = primParticles[i][j];
                double tpt = prt.GetPt();
                double ty = prt.GetY();
                if (m_acc.size() < i || !m_acc[i].init || m_acc[i].getAcceptance(ty + m_ycm, tpt) > RandomGenerators::randgenMT.rand()) ret.Particles.push_back(prt);
                primParticles[i][j].processed = true;
              }
            }
          }
          else {
            for (int j = 0; j < primParticles[i].size(); ++j) {
              if (!primParticles[i][j].processed) {
                flag_repeat = true;
                double DecParam = RandomGenerators::randgenMT.rand(), tsum = 0.;

                std::vector<double> Bratios;
                if (primParticles[i][j].MotherPDGID != 0 ||
                  m_THM->TPS()->ResonanceWidthIntegrationType() != ThermalParticle::eBW) {
                  Bratios = m_THM->TPS()->Particles()[i].BranchingRatiosM(primParticles[i][j].m, false);
                }
                else {
                  Bratios = m_THM->TPS()->Particles()[i].BranchingRatiosM(primParticles[i][j].m, true);
                }

                int DecayIndex = 0;
                for (DecayIndex = 0; DecayIndex < Bratios.size(); ++DecayIndex) {
                  tsum += Bratios[DecayIndex];
                  if (tsum > DecParam) break;
                }
                if (DecayIndex < m_THM->TPS()->Particles()[i].Decays().size()) {
                  std::vector<double> masses(0);
                  std::vector<int> pdgids(0);
                  for (int di = 0; di < m_THM->TPS()->Particles()[i].Decays()[DecayIndex].mDaughters.size(); di++) {
                    int dpdg = m_THM->TPS()->Particles()[i].Decays()[DecayIndex].mDaughters[di];
                    if (m_THM->TPS()->PdgToId(dpdg) == -1) {
                      continue;
                    }
                    masses.push_back(m_THM->TPS()->ParticleByPDG(dpdg).Mass());
                    pdgids.push_back(dpdg);
                  }
                  std::vector<SimpleParticle> decres = ParticleDecays::ManyBodyDecay(primParticles[i][j], masses, pdgids);
                  for (int ind = 0; ind < decres.size(); ind++) {
                    decres[ind].processed = false;
                    if (m_THM->TPS()->PdgToId(decres[ind].PDGID) != -1)
                      primParticles[m_THM->TPS()->PdgToId(decres[ind].PDGID)].push_back(decres[ind]);
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
    }
    return ret;
  }

} // namespace thermalfist