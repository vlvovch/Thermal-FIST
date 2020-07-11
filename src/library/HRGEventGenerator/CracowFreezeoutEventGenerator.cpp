/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/CracowFreezeoutEventGenerator.h"

#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"
#include "HRGEV/ExcludedVolumeModel.h"
#include "HRGEventGenerator/ParticleDecaysMC.h"

namespace thermalfist {

  CracowFreezeoutEventGenerator::CracowFreezeoutEventGenerator()
  {
    m_THM = NULL;
  }

  CracowFreezeoutEventGenerator::CracowFreezeoutEventGenerator(ThermalParticleSystem* TPS, const EventGeneratorConfiguration& config, double T, double RoverTauH, double etamax) :
    m_T(T), m_RoverTauH(RoverTauH), m_EtaMax(etamax)
  {
    SetConfiguration(TPS, config);

    SetMomentumGenerators();
  }

  void CracowFreezeoutEventGenerator::SetParameters(double T, double RoverTauH, double etamax)
  {
    m_T = T;
    m_RoverTauH = RoverTauH;
    m_EtaMax = etamax;

    SetMomentumGenerators();
  }

  void CracowFreezeoutEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);
    if (m_THM != NULL) {
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        const ThermalParticle& part = m_THM->TPS()->Particles()[i];
        //m_MomentumGens.push_back(new RandomGenerators::CracowFreezeoutMomentumGenerator(m_T, m_RoverTauH, m_EtaMax, part.Mass(), part.Statistics(), m_THM->FullIdealChemicalPotential(i)));
        m_MomentumGens.push_back(new RandomGenerators::BoostInvariantMomentumGenerator(new CracowFreezeoutParametrization(m_RoverTauH), GetTkin(), GetEtaMax(), part.Mass(), part.Statistics(), m_THM->FullIdealChemicalPotential(i)));

        double T = m_THM->Parameters().T;
        double Mu = m_THM->FullIdealChemicalPotential(i);
        if (m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW || m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBWconstBR)
          m_BWGens.push_back(new RandomGenerators::ThermalEnergyBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        else
          m_BWGens.push_back(new RandomGenerators::ThermalBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
      }
    }
  }

} // namespace thermalfist