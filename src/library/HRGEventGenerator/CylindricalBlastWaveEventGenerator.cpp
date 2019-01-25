/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/CylindricalBlastWaveEventGenerator.h"

#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"
#include "HRGEV/ExcludedVolumeModel.h"
#include "HRGEventGenerator/ParticleDecaysMC.h"

namespace thermalfist {

  CylindricalBlastWaveEventGenerator::CylindricalBlastWaveEventGenerator() {
    m_THM = NULL;
  }

  CylindricalBlastWaveEventGenerator::CylindricalBlastWaveEventGenerator(ThermalParticleSystem * TPS, const EventGeneratorConfiguration & config, double T, double beta, double etamax, double npow) : 
    m_T(T), m_Beta(beta), m_EtaMax(etamax), m_n(npow)
  {
    SetConfiguration(TPS, config);

    SetMomentumGenerators();
  }

  CylindricalBlastWaveEventGenerator::CylindricalBlastWaveEventGenerator(ThermalModelBase *THM, double T, double beta, double etamax, double npow, bool /*onlyStable*/, EventGeneratorConfiguration::ModelType EV, ThermalModelBase *THMEVVDW) :m_T(T), m_Beta(beta), m_EtaMax(etamax), m_n(npow) {
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

    SetMomentumGenerators();
  }

  void CylindricalBlastWaveEventGenerator::SetParameters(double T, double beta, double etamax, double npow) {
    m_T = T;
    m_Beta = beta;
    m_EtaMax = etamax;
    m_n = npow;

    SetMomentumGenerators();
  }

  void CylindricalBlastWaveEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);
    if (m_THM != NULL) {
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        m_MomentumGens.push_back(new RandomGenerators::SSHMomentumGenerator(m_T, m_Beta, m_EtaMax, m_n, m_THM->TPS()->Particles()[i].Mass()));

        double T = m_THM->Parameters().T;
        double Mu = m_THM->ChemicalPotential(i);
        if (m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW || m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBWconstBR)
          m_BWGens.push_back(new RandomGenerators::ThermalEnergyBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        else
          m_BWGens.push_back(new RandomGenerators::ThermalBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
      }
    }
  }

} // namespace thermalfist