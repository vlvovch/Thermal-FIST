#include "HRGPCE/ThermalModelPCEAnnihilation.h"

#include <iostream>

using namespace std;

namespace thermalfist {

  ThermalModelPCEAnnihilation::ThermalModelPCEAnnihilation(ThermalModelBase* THMbase,
                                                           bool FreezeLonglived,
                                                           double LonglivedResoWidthCut) :
                                                           ThermalModelPCE(THMbase, FreezeLonglived, LonglivedResoWidthCut)
  {
    m_PionAnnihilationNumber = 5.;
    auto stability_flags = RecalculateStabilityFlags();

    SetStabilityFlags(stability_flags);
  }

  void ThermalModelPCEAnnihilation::SetStabilityFlags(const std::vector<int>& StabilityFlags)
  {
    ThermalModelPCE::SetStabilityFlags(StabilityFlags);

    // Indices of conserved hadron numbers
    // In case of baryon-antibaryon annihilation, the index corresponds to the baryon
    m_StableNormal.clear();
    for (int i = 0; i < m_StabilityFlags.size(); ++i) {
      if (m_StabilityFlags[i] == 1) {
        long long tpdgid = ThermalModel()->TPS()->Particle(i).PdgId();
        int idanti = ThermalModel()->TPS()->PdgToId(-tpdgid);
        if (tpdgid == 211 || tpdgid == 111 || tpdgid == -211
          || (idanti != -1 && m_StabilityFlags[idanti] == 2))
          continue;

        m_StableNormal.push_back(i);
      }
    }

    // Indices of annihilating baryons and antibaryons
    m_StableAnnihilate.clear();
    m_StableAnnihilateAnti.clear();
    for (int i = 0; i < m_StabilityFlags.size(); ++i) {
      if (m_StabilityFlags[i] == 2 && ThermalModel()->TPS()->Particle(i).PdgId() > 0) {
        m_StableAnnihilate.push_back(i);
        long long tpdgid = ThermalModel()->TPS()->Particle(i).PdgId();
        int idanti = ThermalModel()->TPS()->PdgToId(-tpdgid);
        if (idanti == -1) {
          throw std::runtime_error("ThermalModelPCEAnnihilation: No antibaryon with pdg code " + std::to_string(-tpdgid) + " exists in the list!");
        }
        m_StableAnnihilateAnti.push_back(idanti);
      }
    }

    // Indices of the pions
    m_Pions.clear();
    for (int i = 0; i < m_StabilityFlags.size(); ++i) {
      long long tpdgid = ThermalModel()->TPS()->Particle(i).PdgId();
      if (tpdgid == 211 || tpdgid == -211 || tpdgid == 111) {
        m_Pions.push_back(i);
      }
    }

  }

  void ThermalModelPCEAnnihilation::CalculatePCE(double param, PCEMode mode)
  {
    if (!m_ChemicalFreezeoutSet) {
      throw std::runtime_error("ThermalModelPCE::CalculatePCE: Tried to make a PCE calculation without setting the chemical freze-out! Call ThermalModelPCE::SetChemicalFreezeout() first.");
    }

    if (!m_StabilityFlagsSet) {
      SetChemicalFreezeout(m_ParametersInit, m_ChemInit);
    }

    double T = param;
    if (mode == PCEMode::AtFixedVolume) {
      // Initial guess for the new temperature
      T = m_ParametersCurrent.T * pow(m_ParametersCurrent.V / param, 1. / 3.);
      m_ParametersCurrent.V = param;
    }
    else {
      // Initial guess for the new volume
      m_ParametersCurrent.V = m_ParametersCurrent.V * pow(m_ParametersCurrent.T / T, 3.);
    }


    BroydenEquationsPCEAnnihilation eqs(this, mode);
    Broyden broydn(&eqs);

    // PCE parameters: volume/temperature + chemical potentials
    std::vector<double> PCEParams(m_StableNormal.size() + m_StableAnnihilate.size() + m_Pions.size(), 0.);
    int stab_index = 0;

    // First go through all the stable species that do not take part in annihilations
    // This excluded the chemical potentials of annihilating antibaryons since they are not independent
    for (int ind = 0; ind < m_StableNormal.size(); ++ind) {
      int i = m_StableNormal[ind];

      PCEParams[stab_index] = m_ChemCurrent[i] * T / m_ParametersCurrent.T + ThermalModel()->TPS()->Particle(i).Mass() * (1. - T / m_ParametersCurrent.T);
      stab_index++;
    }

    // Then the chemical potentials of annihilating baryons
    for (int ind = 0; ind < m_StableAnnihilate.size(); ++ind) {
      int i = m_StableAnnihilate[ind];

      PCEParams[stab_index] = m_ChemCurrent[i] * T / m_ParametersCurrent.T + ThermalModel()->TPS()->Particle(i).Mass() * (1. - T / m_ParametersCurrent.T);
      stab_index++;
    }

    // Finally, the three pions
    for (int ind = 0; ind < m_Pions.size(); ++ind) {
      int i = m_Pions[ind];

      PCEParams[stab_index] = m_ChemCurrent[i] * T / m_ParametersCurrent.T + ThermalModel()->TPS()->Particle(i).Mass() * (1. - T / m_ParametersCurrent.T);
      stab_index++;
    }

    m_ParametersCurrent.T = T;

    // The final parameter is the volume/temperature
    //PCEParams.push_back(m_ParametersCurrent.V);
    if (mode == PCEMode::AtFixedTemperature)
      PCEParams.push_back(m_ParametersCurrent.V);
    else
      PCEParams.push_back(m_ParametersCurrent.T);

    Broyden::BroydenSolutionCriterium crit(1.0E-10);
    //broydn.UseNewton(true);
    //PCEParams = broydn.Solve(PCEParams);
    PCEParams = broydn.Solve(PCEParams, &crit);

    m_ChemCurrent = m_model->ChemicalPotentials();
    if (mode == 0)
      m_ParametersCurrent.V = PCEParams[PCEParams.size() - 1];
    else
      m_ParametersCurrent.T = PCEParams[PCEParams.size() - 1];

    m_model->CalculateFeeddown();

    m_IsCalculated = true;
  }

  int ThermalModelPCEAnnihilation::StableHadronIndexByGlobalId(int globalid)
  {
    for (int i = 0; i < m_StableMapTo.size(); ++i)
      if (m_StableMapTo[i] == globalid)
        return i;
    return -1;
  }

  std::vector<double> ThermalModelPCEAnnihilation::StableChemsFromBroydenInput(const std::vector<double>& x)
  {
    std::vector<double> ret(m_StableComponentsNumber);

    int stab_index = 0;

    // First go through all the stable species that do not take part in annihilations
    for (int ind = 0; ind < m_StableNormal.size(); ++ind) {
      int iglob = m_StableNormal[ind];

      int istab = StableHadronIndexByGlobalId(iglob);

      ret[istab] = x[stab_index];

      stab_index++;
    }

    // Then the chemical potentials of annihilating baryons
    for (int ind = 0; ind < m_StableAnnihilate.size(); ++ind) {
      int iglob = m_StableAnnihilate[ind];

      int istab = StableHadronIndexByGlobalId(iglob);

      ret[istab] = x[stab_index];

      stab_index++;
    }

    // The three pions
    for (int ind = 0; ind < m_Pions.size(); ++ind) {
      int iglob = m_Pions[ind];

      int istab = StableHadronIndexByGlobalId(iglob);

      ret[istab] = x[stab_index];

      stab_index++;
    }

    // Finally, the chemical potentials of the annihilating antibaryons
    for (int ind = 0; ind < m_StableAnnihilate.size(); ++ind) {
      int iglob     = m_StableAnnihilate[ind];
      int iglobanti = m_StableAnnihilateAnti[ind];

      int istab     = StableHadronIndexByGlobalId(iglob);
      int istabanti = StableHadronIndexByGlobalId(iglobanti);

      int istabpip = StableHadronIndexByGlobalId(ThermalModel()->TPS()->PdgToId( 211));
      int istabpiz = StableHadronIndexByGlobalId(ThermalModel()->TPS()->PdgToId( 111));
      int istabpim = StableHadronIndexByGlobalId(ThermalModel()->TPS()->PdgToId(-211));

      // mu_N + mu_{\bar{N}} = <npi> * (mu_pi+ + mu_pi- + mu_pi0) / 3
      ret[istabanti] = -ret[istab] + m_PionAnnihilationNumber / 3. * (ret[istabpip] + ret[istabpiz] + ret[istabpim]);

      stab_index++;
    }

    return ret;
  }

  std::vector<int> ThermalModelPCEAnnihilation::RecalculateStabilityFlags(const vector<long long int> &annihilationpdgs) {
    auto stability_flags = ComputePCEStabilityFlags(ThermalModel()->TPS(), UseSahaForNuclei(), FreezeLonglivedResonances(), LonglivedResonanceWidthCut());
    for(auto& tpdg: annihilationpdgs) {
      // Baryon
      int id = ThermalModel()->TPS()->PdgToId(tpdg);
      if (id != -1)
        stability_flags[id] = 2;

      // Antibaryon
      id = ThermalModel()->TPS()->PdgToId(-tpdg);
      if (id != -1)
        stability_flags[id] = 2;
    }
    return stability_flags;
  }

  /// Conserved quantities in PCE with annihilations
  /// Entropy and (generalizes) stable hadron numbers
  std::vector<double> ThermalModelPCEAnnihilation::BroydenEquationsPCEAnnihilation::Equations(const std::vector<double>& x)
  {
    std::vector<double> ret(x.size(), 0.);

    /// Chemical potentials of PCE hadrons
    std::vector<double> ChemStable = m_THM->StableChemsFromBroydenInput(x);

    /// Chemical potentials of all the particle species
    std::vector<double> Chem(m_THM->ThermalModel()->Densities().size(), 0.);
    for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
      for (int is = 0; is < m_THM->m_EffectiveCharges[i].size(); ++is) {
        Chem[i] += m_THM->m_EffectiveCharges[i][is] * ChemStable[is];
      }
    }

    ThermalModelBase *model = m_THM->ThermalModel();
    model->SetChemicalPotentials(Chem);


    if (m_Mode == PCEMode::AtFixedTemperature) {
      const double& Vtmp = x[x.size() - 1];
      m_THM->m_ParametersCurrent.V = Vtmp;
    }
    else {
      const double& Ttmp = x[x.size() - 1];
      m_THM->m_ParametersCurrent.T = Ttmp;
    }
    double V = m_THM->m_ParametersCurrent.V;
    model->SetParameters(m_THM->m_ParametersCurrent);

    /// Calculate densities of all hadrons
    model->CalculateDensities();

    // Calculate the conserved quantities
    int stab_index = 0;

    // Total yields of stable species that do not take part in annihilations
    for (int ind = 0; ind < m_THM->m_StableNormal.size(); ++ind) {
      int iglob = m_THM->m_StableNormal[ind];

      int istab = m_THM->StableHadronIndexByGlobalId(iglob);

      double totdens = 0.;
      for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
        totdens += m_THM->m_EffectiveCharges[i][istab] * model->Densities()[i];
      }

      ret[stab_index] = (totdens * V) / (m_THM->m_StableDensitiesInit[istab] * m_THM->m_ParametersInit.V) - 1.;

      stab_index++;
    }

    // Net numbers of annihilating baryons
    for (int ind = 0; ind < m_THM->m_StableAnnihilate.size(); ++ind) {
      int iglob = m_THM->m_StableAnnihilate[ind];
      int iglobanti = m_THM->m_StableAnnihilateAnti[ind];

      int istab = m_THM->StableHadronIndexByGlobalId(iglob);
      int istabanti = m_THM->StableHadronIndexByGlobalId(iglobanti);

      double totdens = 0.;
      for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
        totdens += m_THM->m_EffectiveCharges[i][istab] * model->Densities()[i];
        totdens -= m_THM->m_EffectiveCharges[i][istabanti] * model->Densities()[i];
      }

      double netInit = (m_THM->m_StableDensitiesInit[istab] - m_THM->m_StableDensitiesInit[istabanti]) * m_THM->m_ParametersInit.V;
      if (abs(netInit) > 1.e-12) {
        ret[stab_index] = (totdens * V) / netInit - 1.;
      }
      else {
        double sumInit = (m_THM->m_StableDensitiesInit[istab] + m_THM->m_StableDensitiesInit[istabanti]) * m_THM->m_ParametersInit.V;
        ret[stab_index] = totdens * V / sumInit;
      }

      stab_index++;
    }

    // Net-pion numbers
    int istabpip = m_THM->StableHadronIndexByGlobalId(m_THM->ThermalModel()->TPS()->PdgToId(211));
    int istabpim = m_THM->StableHadronIndexByGlobalId(m_THM->ThermalModel()->TPS()->PdgToId(-211));
    int istabpiz = m_THM->StableHadronIndexByGlobalId(m_THM->ThermalModel()->TPS()->PdgToId(111));
    {
      // pi+ - pi-
      double netpicurrent = 0.;
      for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
        netpicurrent += m_THM->m_EffectiveCharges[i][istabpip] * model->Densities()[i];
        netpicurrent -= m_THM->m_EffectiveCharges[i][istabpim] * model->Densities()[i];
      }
      netpicurrent *= V;

      double netpiinit = (m_THM->m_StableDensitiesInit[istabpip] - m_THM->m_StableDensitiesInit[istabpim]) * m_THM->m_ParametersInit.V;

      if (abs(netpiinit) > 1.e-12) {
        ret[stab_index] = netpicurrent / netpiinit - 1.;
      }
      else {
        double sumpiinit = (m_THM->m_StableDensitiesInit[istabpip] + m_THM->m_StableDensitiesInit[istabpim]) * m_THM->m_ParametersInit.V;
        ret[stab_index] = netpicurrent / sumpiinit;
      }

      stab_index++;

      // pi+ - pi0
      double pippi0current = 0.;
      for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
        pippi0current += m_THM->m_EffectiveCharges[i][istabpip] * model->Densities()[i];
        pippi0current -= m_THM->m_EffectiveCharges[i][istabpiz] * model->Densities()[i];
      }
      pippi0current *= V;

      double pippi0init = (m_THM->m_StableDensitiesInit[istabpip] - m_THM->m_StableDensitiesInit[istabpiz]) * m_THM->m_ParametersInit.V;

      if (abs(pippi0init) > 1.e-12) {
        ret[stab_index] = pippi0current / pippi0init - 1.;
      }
      else {
        double sumpippi0init = (m_THM->m_StableDensitiesInit[istabpip] + m_THM->m_StableDensitiesInit[istabpiz]) * m_THM->m_ParametersInit.V;
        ret[stab_index] = pippi0current / sumpippi0init;
      }

      stab_index++;
    }

    // Conservation of (N + \bar{N}) / 2 + (pi+ + pi- + pi0)/<npi>
    {
      double tret = 0.;
      double tretinit = 0.;

      // Pions
      for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
        tret += m_THM->m_EffectiveCharges[i][istabpip] * model->Densities()[i] * 2. / m_THM->m_PionAnnihilationNumber;
        tret += m_THM->m_EffectiveCharges[i][istabpiz] * model->Densities()[i] * 2. / m_THM->m_PionAnnihilationNumber;
        tret += m_THM->m_EffectiveCharges[i][istabpim] * model->Densities()[i] * 2. / m_THM->m_PionAnnihilationNumber;
      }

      // Initial pions
      tretinit += (m_THM->m_StableDensitiesInit[istabpip] + m_THM->m_StableDensitiesInit[istabpiz] + 
        m_THM->m_StableDensitiesInit[istabpim]) * 2. / m_THM->m_PionAnnihilationNumber;

      // Baryon + antibaryon
      for (int ind = 0; ind < m_THM->m_StableAnnihilate.size(); ++ind) {
        int iglob = m_THM->m_StableAnnihilate[ind];
        int iglobanti = m_THM->m_StableAnnihilateAnti[ind];

        int istab = m_THM->StableHadronIndexByGlobalId(iglob);
        int istabanti = m_THM->StableHadronIndexByGlobalId(iglobanti);


        double totdens = 0.;
        for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
          tret += m_THM->m_EffectiveCharges[i][istab] * model->Densities()[i];
          tret += m_THM->m_EffectiveCharges[i][istabanti] * model->Densities()[i];
        }

        tretinit += m_THM->m_StableDensitiesInit[istab] + m_THM->m_StableDensitiesInit[istabanti];
      }

      ret[stab_index] = (tret * V) / (tretinit * m_THM->m_ParametersInit.V) - 1.;

      stab_index++;
    }

    // Entropy
    ret[ret.size() - 1] = model->CalculateEntropyDensity() * V - m_THM->m_EntropyDensityInit * m_THM->m_ParametersInit.V;
    ret[ret.size() - 1] /= m_THM->m_EntropyDensityInit * m_THM->m_ParametersInit.V;

    return ret;
  }

} // namespace thermalfist
