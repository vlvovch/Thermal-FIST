#include "HRGEventGenerator/SREventGenerator.h"

#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"
#include "HRGEV/ExcludedVolumeModel.h"
#include "HRGEventGenerator/ParticleDecays.h"

namespace thermalfist {

  SiemensRasmussenEventGenerator::SiemensRasmussenEventGenerator() {
    m_THM = NULL;
  }


  SiemensRasmussenEventGenerator::SiemensRasmussenEventGenerator(ThermalModelBase *THM, double T, double beta, bool onlyStable, int EV, ThermalModelBase *THMEVVDW) : m_T(T), m_Beta(beta) {
    //SiemensRasmussenEventGenerator::SiemensRasmussenEventGenerator(ThermalModelBase *THM, double T, double beta, bool onlyStable, int EV, ExcludedVolumeModel *exmod) : m_T(T), m_Beta(beta) { 
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

    m_OnlyStable = onlyStable;

    SetMomentumGenerators();

    if (m_THM != NULL)
      m_acc.resize(m_THM->TPS()->Particles().size());
  }

  void SiemensRasmussenEventGenerator::SetParameters(double T, double beta) {
    m_T = T;
    m_Beta = beta;

    SetMomentumGenerators();
  }

  void SiemensRasmussenEventGenerator::SetThermalModel(ThermalModelBase *THM, bool regen) {
    m_THM = THM;
    if (regen) {
      SetMomentumGenerators();
    }
  }

  void SiemensRasmussenEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);
    if (m_THM != NULL) {
      for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        m_MomentumGens.push_back(new RandomGenerators::SiemensRasmussenGenerator(m_T, m_Beta, m_THM->TPS()->Particles()[i].Mass()));

        double T = m_THM->Parameters().T;
        double Mu = m_THM->ChemicalPotential(i);
        if (m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW)
          m_BWGens.push_back(new RandomGenerators::ThermalEnergyBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        else
          m_BWGens.push_back(new RandomGenerators::ThermalBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        //if (!m_THM->TPS()->Particles()[i].IsStable()) m_BWGens.push_back(RandomGenerators::BreitWignerGenerator(m_THM->TPS()->Particles()[i].Mass(), m_THM->TPS()->Particles()[i].ResonanceWidth(), std::max(m_THM->TPS()->Particles()[i].Mass() - 2.*m_THM->TPS()->Particles()[i].ResonanceWidth(), m_THM->TPS()->Particles()[i].DecayThresholdMass())));
        //else m_BWGens.push_back(RandomGenerators::BreitWignerGenerator(m_THM->TPS()->Particles()[i].Mass(), 0.1, 0.2));
      }
    }
  }

} // namespace thermalfist