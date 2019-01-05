/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/SphericalBlastWaveEventGenerator.h"

#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"
#include "HRGEV/ExcludedVolumeModel.h"
#include "HRGEventGenerator/ParticleDecays.h"

namespace thermalfist {

  SphericalBlastWaveEventGenerator::SphericalBlastWaveEventGenerator() {
    m_THM = NULL;
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

    SetConfiguration(THM->Parameters(), ensemble, modeltype, THM->TPS(), THM, THMEVVDW);

    m_OnlyStable = onlyStable;

    SetMomentumGenerators();

    if (m_THM != NULL)
      m_acc.resize(m_THM->TPS()->Particles().size());
  }

  void SphericalBlastWaveEventGenerator::SetParameters(double T, double beta) {
    m_T = T;
    m_Beta = beta;

    SetMomentumGenerators();
  }

  void SphericalBlastWaveEventGenerator::SetThermalModel(ThermalModelBase *THM, bool regen) {
    m_THM = THM;
    if (regen) {
      SetMomentumGenerators();
    }
  }

  void SphericalBlastWaveEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);
    if (m_THM != NULL) {
      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        m_MomentumGens.push_back(new RandomGenerators::SiemensRasmussenGenerator(m_T, m_Beta, m_THM->TPS()->Particles()[i].Mass()));

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