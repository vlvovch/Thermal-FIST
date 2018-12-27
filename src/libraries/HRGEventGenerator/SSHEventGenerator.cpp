/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/SSHEventGenerator.h"

#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"
#include "HRGEV/ExcludedVolumeModel.h"
#include "HRGEventGenerator/ParticleDecays.h"

namespace thermalfist {

  SSHEventGenerator::SSHEventGenerator() {
    m_THM = NULL;
  }

  //SSHEventGenerator::SSHEventGenerator(ThermalModelBase *THM, double T, double beta, double etamax, bool onlyStable, int EV, ExcludedVolumeModel *exmod) :m_T(T), m_Beta(beta), m_EtaMax(etamax) {
  SSHEventGenerator::SSHEventGenerator(ThermalModelBase *THM, double T, double beta, double etamax, double npow, bool onlyStable, int EV, ThermalModelBase *THMEVVDW) :m_T(T), m_Beta(beta), m_EtaMax(etamax), m_n(npow) {
    EventGeneratorConfiguration::ModelType modeltype = EventGeneratorConfiguration::PointParticle;
    if (EV == 1) modeltype = EventGeneratorConfiguration::DiagonalEV;
    else if (EV == 2) modeltype = EventGeneratorConfiguration::CrosstermsEV;
    else if (EV == 3) modeltype = EventGeneratorConfiguration::MeanFieldEV;
    else if (EV == 4) modeltype = EventGeneratorConfiguration::QvdW;
    EventGeneratorConfiguration::Ensemble ensemble = EventGeneratorConfiguration::GCE;
    if (THM->Ensemble() == ThermalModelBase::CE)
      ensemble = EventGeneratorConfiguration::CE;
    if (THM->Ensemble() == ThermalModelBase::SCE)
      ensemble = EventGeneratorConfiguration::SCE;
    //SetConfiguration(THM->Parameters(), ensemble, modeltype, THM->TPS(), THM, exmod);
    SetConfiguration(THM->Parameters(), ensemble, modeltype, THM->TPS(), THM, THMEVVDW);

    SetMomentumGenerators();

    if (m_THM != NULL)
      m_acc.resize(m_THM->TPS()->Particles().size());
  }

  void SSHEventGenerator::SetParameters(double T, double beta, double etamax) {
    m_T = T;
    m_Beta = beta;
    m_EtaMax = etamax;

    SetMomentumGenerators();
  }

  void SSHEventGenerator::SetThermalModel(ThermalModelBase *THM, bool regen) {
    m_THM = THM;
    if (regen) {
      SetMomentumGenerators();
    }
  }

  void SSHEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);
    if (m_THM != NULL) {
      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        m_MomentumGens.push_back(new RandomGenerators::SSHGenerator(m_T, m_Beta, m_EtaMax, m_n, m_THM->TPS()->Particles()[i].Mass()));

        double T = m_THM->Parameters().T;
        double Mu = m_THM->ChemicalPotential(i);
        if (m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW || m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBWconstBR)
          m_BWGens.push_back(new RandomGenerators::ThermalEnergyBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        else
          m_BWGens.push_back(new RandomGenerators::ThermalBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        //if (!m_THM->TPS()->Particles()[i].IsStable()) m_BWGens.push_back(RandomGenerators::BreitWignerGenerator(m_THM->TPS()->Particles()[i].Mass(), m_THM->TPS()->Particles()[i].ResonanceWidth(), std::max(m_THM->TPS()->Particles()[i].Mass() - 2.*m_THM->TPS()->Particles()[i].ResonanceWidth(), m_THM->TPS()->Particles()[i].DecayThresholdMass())));
        //else m_BWGens.push_back(RandomGenerators::BreitWignerGenerator(m_THM->TPS()->Particles()[i].Mass(), 0.1, 0.2));
      }
    }
  }

} // namespace thermalfist