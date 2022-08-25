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

  SphericalBlastWaveEventGenerator::SphericalBlastWaveEventGenerator() : 
    EventGeneratorBase() 
  {
  }

  SphericalBlastWaveEventGenerator::SphericalBlastWaveEventGenerator(ThermalParticleSystem * TPS, const EventGeneratorConfiguration & config, double T, double beta) : 
    EventGeneratorBase(),
    m_T(T), m_Beta(beta)
  {
    SetConfiguration(TPS, config);

    //SetMomentumGenerators();
  }


  SphericalBlastWaveEventGenerator::SphericalBlastWaveEventGenerator(ThermalModelBase *THM, double T, double beta, bool /*onlyStable*/, EventGeneratorConfiguration::ModelType EV, ThermalModelBase *THMEVVDW) : 
    EventGeneratorBase(),
    m_T(T), m_Beta(beta) {
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
    m_ParametersSet = false;
    //SetMomentumGenerators();
  }

  void SphericalBlastWaveEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);

    double gamma = 1. / sqrt(1 - GetBeta() * GetBeta());
    double R = pow(3. * m_THM->Volume() / (4. * xMath::Pi()) / gamma, 1. / 3.);

    if (m_THM != NULL) {
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        const ThermalParticle& part = m_THM->TPS()->Particles()[i];
        //m_MomentumGens.push_back(new RandomGenerators::SiemensRasmussenMomentumGeneratorGeneralized(m_T, m_Beta, m_THM->TPS()->Particles()[i].Mass()));
        m_MomentumGens.push_back(new RandomGenerators::SiemensRasmussenMomentumGeneratorGeneralized(m_T, m_Beta, R, part.Mass(), part.Statistics(), m_THM->FullIdealChemicalPotential(i)));

        double T = m_THM->Parameters().T;
        double Mu = m_THM->FullIdealChemicalPotential(i);
        if (m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW || m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBWconstBR)
          m_BWGens.push_back(new RandomGenerators::ThermalEnergyBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        else
          m_BWGens.push_back(new RandomGenerators::ThermalBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
      }
    }
  }

  //void SphericalBlastWaveEventGenerator::SetParameters()
  //{
  //  SetMomentumGenerators();
  //  m_ParametersSet = true;
  //}

} // namespace thermalfist