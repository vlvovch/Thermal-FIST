/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/SphericalBlastWaveEventGenerator.h"

#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"
#include "HRGEV/ExcludedVolumeModel.h"
#include "HRGEventGenerator/ParticleDecaysMC.h"

namespace thermalfist {

  SphericalBlastWaveEventGenerator::SphericalBlastWaveEventGenerator() {
    m_THM = NULL;
  }

  SphericalBlastWaveEventGenerator::SphericalBlastWaveEventGenerator(ThermalParticleSystem * TPS, const EventGeneratorConfiguration & config, double T, double beta) : m_T(T), m_Beta(beta)
  {
    SetConfiguration(TPS, config);

    SetMomentumGenerators();
  }


  SphericalBlastWaveEventGenerator::SphericalBlastWaveEventGenerator(ThermalModelBase *THM, double T, double beta, bool onlyStable, EventGeneratorConfiguration::ModelType EV, ThermalModelBase *THMEVVDW) : m_T(T), m_Beta(beta) {
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

    SetConfiguration(THM->TPS(), config);

    SetMomentumGenerators();
  }

  void SphericalBlastWaveEventGenerator::SetBlastWaveParameters(double T, double beta) {
    m_T = T;
    m_Beta = beta;

    SetMomentumGenerators();
  }

  void SphericalBlastWaveEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);
    if (m_THM != NULL) {
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        m_MomentumGens.push_back(new RandomGenerators::SiemensRasmussenMomentumGenerator(m_T, m_Beta, m_THM->TPS()->Particles()[i].Mass()));

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