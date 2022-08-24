/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/CylindricalBlastWaveEventGenerator.h"

#include <iostream>
#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"
#include "HRGBase/ThermalModelBase.h"
#include "HRGEV/ExcludedVolumeModel.h"
#include "HRGEventGenerator/ParticleDecaysMC.h"

using namespace std;

namespace thermalfist {

  CylindricalBlastWaveEventGenerator::CylindricalBlastWaveEventGenerator(ThermalParticleSystem * TPS, const EventGeneratorConfiguration & config, double T, double betas, double etamax, double npow, double Rperp) : 
    EventGeneratorBase(),
    m_T(T), m_BetaS(betas), m_EtaMax(etamax), m_n(npow), m_Rperp(Rperp)
  {
    SetConfiguration(TPS, config);

    //SetMomentumGenerators();
  }

  CylindricalBlastWaveEventGenerator::CylindricalBlastWaveEventGenerator(ThermalModelBase *THM, double T, double betas, double etamax, double npow, bool /*onlyStable*/, EventGeneratorConfiguration::ModelType EV, ThermalModelBase *THMEVVDW) :
    EventGeneratorBase(),
    m_T(T), m_BetaS(betas), m_EtaMax(etamax), m_n(npow), m_Rperp(6.5) {
    EventGeneratorConfiguration::ModelType modeltype = EV;
    EventGeneratorConfiguration::Ensemble ensemble = EventGeneratorConfiguration::GCE;
    if (THM->Ensemble() == ThermalModelBase::CE)
      ensemble = EventGeneratorConfiguration::CE;
    if (THM->Ensemble() == ThermalModelBase::SCE)
      ensemble = EventGeneratorConfiguration::SCE;
    if (THM->Ensemble() == ThermalModelBase::CCE)
      ensemble = EventGeneratorConfiguration::CCE;

    EventGeneratorConfiguration config;

    config.fEnsemble = ensemble;
    config.fModelType = modeltype;
    config.CFOParameters = THM->Parameters();
    config.B = THM->Parameters().B;
    config.Q = THM->Parameters().Q;
    config.S = THM->Parameters().S;
    config.C = THM->Parameters().C;

    config.bij.resize(THMEVVDW->ComponentsNumber());
    for (size_t i = 0; i < config.bij.size(); ++i) {
      config.bij[i].resize(THMEVVDW->ComponentsNumber());
      for (size_t j = 0; j < config.bij.size(); ++j) {
        config.bij[i][j] = THMEVVDW->VirialCoefficient(i, j);
      }
    }

    config.aij.resize(THMEVVDW->ComponentsNumber());
    for (size_t i = 0; i < config.bij.size(); ++i) {
      config.aij[i].resize(THMEVVDW->ComponentsNumber());
      for (size_t j = 0; j < config.bij.size(); ++j) {
        config.aij[i][j] = THMEVVDW->AttractionCoefficient(i, j);
      }
    }

    SetConfiguration(THMEVVDW->TPS(), config);

    //SetMomentumGenerators();
  }

  void CylindricalBlastWaveEventGenerator::SetParameters(double T, double betas, double etamax, double npow) {
    m_T = T;
    m_BetaS = betas;
    m_EtaMax = etamax;
    m_n = npow;
    m_ParametersSet = false;

    //SetMomentumGenerators();
  }

  void CylindricalBlastWaveEventGenerator::SetMeanBetaT(double betaT)
  {
    m_BetaS = (2. + m_n) / 2. * betaT;
    m_ParametersSet = false;

    //SetMomentumGenerators();
  }

  void CylindricalBlastWaveEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);

    // Find \tau_H from Veff / (\delta \eta) = \pi \tau R^2 * I where I is an integral over transverse velocity profile computed numerically
    double tau = m_THM->Volume() / (2. * GetEtaMax()) / (2. * xMath::Pi()) / GetRperp() / GetRperp() / GetVeffIntegral();

    if (m_THM != NULL) {
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        const ThermalParticle& part = m_THM->TPS()->Particles()[i];
        m_MomentumGens.push_back(new RandomGenerators::BoostInvariantMomentumGenerator(new CylindricalBlastWaveParametrization(GetBetaSurface(), GetNPow(), tau, GetRperp()), GetTkin(), GetEtaMax(), part.Mass(), part.Statistics(), m_THM->FullIdealChemicalPotential(i)));

        double T = m_THM->Parameters().T;
        double Mu = m_THM->FullIdealChemicalPotential(i);
        if (m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW || m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBWconstBR)
          m_BWGens.push_back(new RandomGenerators::ThermalEnergyBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        else
          m_BWGens.push_back(new RandomGenerators::ThermalBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
      }
    }
  }

  void CylindricalBlastWaveEventGenerator::RecalculateTotalConservedNumbers()
  {
    m_THM->SetTemperature(m_Config.CFOParameters.T);
    m_THM->SetBaryonChemicalPotential(m_Config.CFOParameters.muB);
    m_THM->SetElectricChemicalPotential(m_Config.CFOParameters.muQ);
    m_THM->SetStrangenessChemicalPotential(m_Config.CFOParameters.muS);

    if (m_Config.CFOParameters.gammaq != 1.0)
      m_THM->SetGammaq(m_Config.CFOParameters.gammaq);

    if (m_Config.CFOParameters.gammaS != 1.0)
      m_THM->SetGammaS(m_Config.CFOParameters.gammaS);

    if (m_Config.CFOParameters.gammaC != 1.0)
      m_THM->SetGammaC(m_Config.CFOParameters.gammaC);

    double Veff = m_Config.CFOParameters.V;

    m_THM->CalculatePrimordialDensities();
    double totB = m_THM->BaryonDensity() * Veff;
    double totQ = m_THM->ElectricChargeDensity() * Veff;
    double totS = m_THM->StrangenessDensity() * Veff;
    double totC = m_THM->CharmDensity() * Veff;

    if (m_Config.fEnsemble == EventGeneratorConfiguration::CE) {
      cout << endl;
      cout << "Integrated values of conserved charges:" << endl;
      cout << "B = " << totB << endl;
      cout << "Q = " << totQ << endl;
      cout << "S = " << totS << endl;
      cout << "C = " << totC << endl;

      double Vcorr = 1.;
      if (m_Config.CanonicalB) {

        if (abs(totB) > 1.e-6) {
          Vcorr = round(totB) / totB;
          Veff *= Vcorr;
          cout << "Volume rescaling factor Vcorr = " << Vcorr << endl;
          cout << "B: " << totB << " -> " << round(totB) << endl;
        }

        m_Config.B = round(totB);
      }

      if (m_Config.CanonicalQ) {
        cout << "Q: " << totQ * Vcorr << " -> " << round(totQ * Vcorr) << endl;
        m_Config.Q = round(totQ * Vcorr);
      }

      if (m_Config.CanonicalS && m_Config.S != 0) {
        cout << "S: " << totS * Vcorr << " -> " << round(totS * Vcorr) << endl;
        m_Config.S = round(totS * Vcorr);
      }

      if (m_Config.CanonicalC && m_Config.C != 0) {
        cout << "C: " << totC * Vcorr << " -> " << round(totC * Vcorr) << endl;
        m_Config.C = round(totC * Vcorr);
      }
    }

    m_THM->SetVolume(Veff);
    m_THM->SetCanonicalVolume(Veff);
    if (m_Config.fEnsemble != EventGeneratorConfiguration::GCE)
      PrepareMultinomials();
  }

  double CylindricalBlastWaveEventGenerator::GetVeffIntegral() const
  {
    double ret = 0.0;

    std::vector<double> xleg, wleg;
    NumericalIntegration::GetCoefsIntegrateLegendre32(0., 1., &xleg, &wleg);

    for (int iint = 0; iint < xleg.size(); ++iint) {
      double zeta  = xleg[iint];
      double w     = wleg[iint];
      double betar = pow(zeta, GetNPow()) * GetBetaSurface();
      ret += w * zeta / sqrt(1. - betar * betar);
    }

    return ret;
  }

} // namespace thermalfist