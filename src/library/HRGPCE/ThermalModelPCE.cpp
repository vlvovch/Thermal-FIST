#include "HRGPCE/ThermalModelPCE.h"

#include <iostream>

using namespace std;

#include <Eigen/Dense>

using namespace Eigen;

namespace thermalfist {

  

  ThermalModelPCE::ThermalModelPCE(ThermalModelBase * THMbase, bool FreezeLongLived) : m_model(THMbase), m_ChemicalFreezeoutSet(false)
  {
    PrepareNucleiForPCE(THMbase->TPS());
    
    for (int i = 0; i < m_model->TPS()->Particles().size(); ++i) {
      int stability = static_cast<int> (m_model->TPS()->Particles()[i].IsStable());
      m_StabilityFlags.push_back(stability);
      //m_StabilityFlags.push_back( static_cast<int> (m_model->TPS()->Particles()[i].DecayType() != ParticleDecay::Strong ) );
    }

    // Automatic, all strongly decaying resonances and all nuclei are not frozen
    m_StabilityFlags.resize(0);
    for (int i = 0; i < m_model->TPS()->Particles().size(); ++i) {
      const ThermalParticle& part = m_model->TPS()->Particles()[i];

      int frozen = 0;

      frozen = static_cast<int> (part.DecayType() != ParticleDecayType::Strong && abs(part.BaryonCharge()) <= 1);

      // Yields of long-lived resonances are also frozen
      if (FreezeLongLived && part.DecayType() == ParticleDecayType::Strong && part.ResonanceWidth() < 0.015 && abs(part.BaryonCharge()) <= 1)
        frozen = 1;

      // The special case of K0
      if (part.PdgId() == 310 || part.PdgId() == 130)
        frozen = 0;
      if (part.PdgId() == 311 || part.PdgId() == -311)
        frozen = 1;

      m_StabilityFlags.push_back(frozen);
    }

    SetStabilityFlags(m_StabilityFlags);
  }

  void ThermalModelPCE::SetStabilityFlags(const std::vector<int>& StabilityFlags)
  {
    m_StabilityFlags = StabilityFlags;

    m_StableComponentsNumber = 0;
    for (int i = 0; i < m_model->TPS()->Particles().size(); ++i) {
      m_model->TPS()->Particle(i).SetStable(m_StabilityFlags[i] != 0);
      if (m_StabilityFlags[i] != 0)
        m_StableComponentsNumber++;
    }

    m_model->TPS()->FillResonanceDecays();

    m_StableMapTo = std::vector<int>(0);
    m_EffectiveCharges = std::vector< std::vector<double> >(m_StabilityFlags.size(), std::vector<double>(m_StableComponentsNumber, 0.));
    int stab_index = 0;
    for (int i = 0; i < m_model->TPS()->Particles().size(); ++i) {
      const ThermalParticle& part = m_model->TPS()->Particles()[i];
      if (part.IsStable()) {
        if (stab_index >= m_EffectiveCharges[0].size()) {
          printf("**ERROR** ThermalModelPCE::SetStabilityFlags: Wrong number of stable components!\n");
          exit(1);
        }
        m_EffectiveCharges[i][stab_index] = 1.;
        const ThermalParticleSystem::DecayContributionsToParticle& decayContributions = m_model->TPS()->DecayContributionsByFeeddown()[Feeddown::StabilityFlag][i];
        for (int j = 0; j < decayContributions.size(); ++j) {
          const pair<double,int> Contrs = decayContributions[j];
          m_EffectiveCharges[Contrs.second][stab_index] = Contrs.first;
        }
        m_StableMapTo.push_back(i);
        stab_index++;
      }
    }


    m_ChemicalFreezeoutSet = false;
    m_IsCalculated = false;
  }

  void ThermalModelPCE::SetChemicalFreezeout(const ThermalModelParameters & params, const std::vector<double>& ChemInit)
  {
    m_ParametersInit = params;
    m_ChemInit = ChemInit;
    m_model->SetParameters(m_ParametersInit);
    m_model->SetChemicalPotentials(m_ChemInit);
    m_model->CalculateDensities();
    m_DensitiesInit = m_model->TotalDensities();
    m_EntropyDensityInit = m_model->CalculateEntropyDensity();

    m_ParametersCurrent = m_ParametersInit;
    m_ChemCurrent = m_ChemInit;

    m_ChemicalFreezeoutSet = true;
    m_IsCalculated = false;
  }

  void ThermalModelPCE::CalculatePCE(double T)
  {
    //m_ParametersCurrent   = m_ParametersInit;
    m_ParametersCurrent.V = m_ParametersCurrent.V * pow(m_ParametersCurrent.T / T, 3.);
    //m_ParametersCurrent.T = T;
    //m_ChemCurrent = m_ChemInit;

    BroydenEquationsPCE eqs(this);
    Broyden broydn(&eqs);
    //BroydenJacobianPCE jac(this);
    //Broyden broydn(&eqs, &jac);

    std::vector<double> PCEParams(m_StableComponentsNumber, 0.);
    int stab_index = 0;
    for (int i = 0; i < m_StabilityFlags.size(); ++i) {
      if (m_StabilityFlags[i]) {
        if (stab_index >= m_EffectiveCharges[0].size()) {
          printf("**ERROR** ThermalModelPCE::CalculatePCE: Wrong number of stable components!\n");
          exit(1);
        }

        //PCEParams[stab_index] = m_ChemCurrent[i];
        // Improved initial condition
        PCEParams[stab_index] = m_ChemCurrent[i] * T / m_ParametersCurrent.T + ThermalModel()->TPS()->Particle(i).Mass() * (1. - T / m_ParametersCurrent.T);
        stab_index++;
      }
    }
    
    m_ParametersCurrent.T = T;

    PCEParams.push_back(m_ParametersCurrent.V);

    PCEParams = broydn.Solve(PCEParams);

    m_ChemCurrent = m_model->ChemicalPotentials();
    //m_ChemCurrent = std::vector<double>(PCEParams.begin(), PCEParams.begin() + m_StableComponentsNumber);
    m_ParametersCurrent.V = PCEParams[PCEParams.size() - 1];
    
    m_IsCalculated = true;
  }

  void ThermalModelPCE::PrepareNucleiForPCE(ThermalParticleSystem * TPS)
  {
    ThermalParticleSystem &parts = *TPS;
    vector<int> nuclpdgs;
    vector<string> nuclnames;
    vector< vector<int> > nuclcontent;

    nuclpdgs.push_back(1000010020);
    nuclnames.push_back("d");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);

    nuclpdgs.push_back(1000020030);
    nuclnames.push_back("He3");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);

    nuclpdgs.push_back(1010010030);
    nuclnames.push_back("H3La");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1000020040);
    nuclnames.push_back("He4");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);

    nuclpdgs.push_back(1010000020);
    nuclnames.push_back("LambdaNeutron");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1010010020);
    nuclnames.push_back("LambdaProton");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1020000020);
    nuclnames.push_back("DiLambda");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(3122);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1000010030);
    nuclnames.push_back("Triton");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);

    nuclpdgs.push_back(1010020040);
    nuclnames.push_back("He4La");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1010010040);
    nuclnames.push_back("H4La");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1010020050);
    nuclnames.push_back("He5La");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1020010020);
    nuclnames.push_back("XiProton");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(3322);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);

    nuclpdgs.push_back(1030000020);
    nuclnames.push_back("OmegaProton");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(3334);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);

    nuclpdgs.push_back(1040000020);
    nuclnames.push_back("DiXi0");
    nuclcontent.push_back(vector<int>());
    nuclcontent[nuclcontent.size() - 1].push_back(3322);
    nuclcontent[nuclcontent.size() - 1].push_back(3322);

    // Fill "decays" of light nuclei
    for (int i = 0; i < parts.Particles().size(); ++i) {
      ThermalParticle &part = parts.Particle(i);

      for (int j = 0; j < nuclpdgs.size(); ++j) {
        if (part.PdgId() == nuclpdgs[j] || part.PdgId() == -nuclpdgs[j]) {
          part.Decays().resize(0);
          int sign = 1;
          if (part.PdgId() < 0)
            sign = -1;
          std::vector<long long> daughters;
          for (int k = 0; k < nuclcontent[j].size(); ++k) {
            daughters.push_back(sign * nuclcontent[j][k]);
          }
          part.Decays().push_back(ParticleDecayChannel(1., daughters));
        }
      }
    }
  }

  std::vector<double> ThermalModelPCE::BroydenEquationsPCE::Equations(const std::vector<double>& x)
  {
    std::vector<double> ret(x.size(), 0.);

    std::vector<double> Chem(m_THM->ThermalModel()->Densities().size(), 0.);
    for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
      for (int is = 0; is < m_THM->m_EffectiveCharges[i].size(); ++is) {
        Chem[i] += m_THM->m_EffectiveCharges[i][is] * x[is];
      }
    }

    ThermalModelBase *model = m_THM->ThermalModel();

    model->SetChemicalPotentials(Chem);
    const double& V = x[x.size() - 1];
    m_THM->m_ParametersCurrent.V = V;
    model->SetParameters(m_THM->m_ParametersCurrent);
    model->CalculateDensities();

    for (int is = 0; is < m_THM->m_StableComponentsNumber; ++is) {
      ret[is] = model->TotalDensities()[ m_THM->m_StableMapTo[is] ] * V - m_THM->m_DensitiesInit[m_THM->m_StableMapTo[is]] * m_THM->m_ParametersInit.V;
    }

    ret[ret.size() - 1] = model->CalculateEntropyDensity() * V - m_THM->m_EntropyDensityInit * m_THM->m_ParametersInit.V;

    return ret;
  }

} // namespace thermalfist

