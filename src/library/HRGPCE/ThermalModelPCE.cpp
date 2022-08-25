#include "HRGPCE/ThermalModelPCE.h"

#include <iostream>

using namespace std;

#include <Eigen/Dense>

using namespace Eigen;

namespace thermalfist {

  

  ThermalModelPCE::ThermalModelPCE(ThermalModelBase * THMbase, bool FreezeLonglived, double LonglivedResoWidthCut) :
    m_model(THMbase), 
    m_UseSahaForNuclei(true),
    m_FreezeLonglivedResonances(FreezeLonglived),
    m_ResoWidthCut(LonglivedResoWidthCut),
    m_ChemicalFreezeoutSet(false), 
    m_StabilityFlagsSet(false),
    m_IsCalculated(false)
  {
    m_model->UsePartialChemicalEquilibrium(true);
  }

  void ThermalModelPCE::SetStabilityFlags(const std::vector<int>& StabilityFlags)
  {
    // A helper ThermalParticleSystem instance to compute the PCE effective charges for all particles
    ThermalParticleSystem TPShelper = ThermalParticleSystem(*m_model->TPS());

    // Set the nucleon content of the known stable light nuclei as "decay" products
    PrepareNucleiForPCE(&TPShelper);

    m_StabilityFlags = StabilityFlags;

    m_StableComponentsNumber = 0;
    for (int i = 0; i < TPShelper.Particles().size(); ++i) {
      TPShelper.Particle(i).SetStable(m_StabilityFlags[i] != 0);
      if (m_StabilityFlags[i] != 0)
        m_StableComponentsNumber++;
    }

    TPShelper.FillResonanceDecays();

    m_StableMapTo = std::vector<int>(0);
    m_EffectiveCharges = std::vector< std::vector<double> >(m_StabilityFlags.size(), std::vector<double>(m_StableComponentsNumber, 0.));
    int stab_index = 0;
    for (int i = 0; i < TPShelper.Particles().size(); ++i) {
      const ThermalParticle& part = TPShelper.Particles()[i];
      if (part.IsStable()) {
        if (stab_index >= m_EffectiveCharges[0].size()) {
          printf("**ERROR** ThermalModelPCE::SetStabilityFlags: Wrong number of stable components!\n");
          exit(1);
        }
        m_EffectiveCharges[i][stab_index] = 1.;
        const ThermalParticleSystem::DecayContributionsToParticle& decayContributions = TPShelper.DecayContributionsByFeeddown()[Feeddown::StabilityFlag][i];
        for (int j = 0; j < decayContributions.size(); ++j) {
          const pair<double,int> Contrs = decayContributions[j];
          m_EffectiveCharges[Contrs.second][stab_index] = Contrs.first;
        }
        m_StableMapTo.push_back(i);
        stab_index++;
      }
    }

    ApplyFixForBoseCondensation();

    m_StabilityFlagsSet = true;
    m_ChemicalFreezeoutSet = false;
    m_IsCalculated = false;
  }

  void ThermalModelPCE::SetChemicalFreezeout(const ThermalModelParameters & params, const std::vector<double>& ChemInit)
  {
    if (!m_StabilityFlagsSet) {
      SetStabilityFlags(ComputePCEStabilityFlags(m_model->TPS(), UseSahaForNuclei(), FreezeLonglivedResonances(), LonglivedResonanceWidthCut()));
    }
    
    m_ParametersInit = params;
    
    m_model->SetParameters(m_ParametersInit);
    if (ChemInit.size() == m_model->ComponentsNumber()) {
      m_model->SetChemicalPotentials(ChemInit);
      m_ChemInit = ChemInit;
    }
    else {
      m_model->FillChemicalPotentials();
      m_ChemInit = m_model->ChemicalPotentials();
    }
    //m_model->CalculateDensities();
    m_model->CalculatePrimordialDensities();

    //m_DensitiesInit = m_model->TotalDensities();

    m_StableDensitiesInit = std::vector<double>(m_StableComponentsNumber, 0.);
    for (int is = 0; is < m_StableDensitiesInit.size(); ++is) {
      double totdens = 0.;
      for (int i = 0; i < m_EffectiveCharges.size(); ++i) {
        totdens += m_EffectiveCharges[i][is] * m_model->Densities()[i];
      }
      m_StableDensitiesInit[is] = totdens;
    }


    m_EntropyDensityInit = m_model->EntropyDensity();

    m_ParticleDensityInit = m_model->HadronDensity();

    m_ParametersCurrent = m_ParametersInit;
    m_ChemCurrent = m_ChemInit;

    m_ChemicalFreezeoutSet = true;
    m_IsCalculated = false;
  }

  void ThermalModelPCE::CalculatePCE(double param, PCEMode mode)
  {
    if (!m_ChemicalFreezeoutSet) {
      printf("**ERROR** ThermalModelPCE::CalculatePCE:"
             "Tried to make a PCE calculation without setting the chemical freze-out!"
             "Call ThermalModelPCE::SetChemicalFreezeout() first.\n");
      exit(1);
    }

    if (!m_StabilityFlagsSet) {
      SetChemicalFreezeout(m_ParametersInit, m_ChemInit);
    }
    
    double T = param;
    if (mode == 1) {
      // Initial guess for the new temperature
      T = m_ParametersCurrent.T * pow(m_ParametersCurrent.V / param, 1. / 3.);
      m_ParametersCurrent.V = param;
    }
    else {
      // Initial guess for the new volume
      m_ParametersCurrent.V = m_ParametersCurrent.V * pow(m_ParametersCurrent.T / T, 3.);
    }
    
    

    BroydenEquationsPCE eqs(this, mode);
    Broyden broydn(&eqs);

    std::vector<double> PCEParams(m_StableComponentsNumber, 0.);
    int stab_index = 0;
    for (int i = 0; i < m_StabilityFlags.size(); ++i) {
      if (m_StabilityFlags[i]) {
        if (stab_index >= m_EffectiveCharges[0].size()) {
          printf("**ERROR** ThermalModelPCE::CalculatePCE: Wrong number of stable components!\n");
          exit(1);
        }

        //PCEParams[stab_index] = m_ChemCurrent[i];
        // Improved initial guesses for the PCE chemical potentials
        PCEParams[stab_index] = m_ChemCurrent[i] * T / m_ParametersCurrent.T + ThermalModel()->TPS()->Particle(i).Mass() * (1. - T / m_ParametersCurrent.T);
        stab_index++;
      }
    }
    
    m_ParametersCurrent.T = T;

    if (mode == 0)
      PCEParams.push_back(m_ParametersCurrent.V);
    else
      PCEParams.push_back(m_ParametersCurrent.T);

    PCEParams = broydn.Solve(PCEParams);

    m_ChemCurrent = m_model->ChemicalPotentials();
    if (mode == 0)
      m_ParametersCurrent.V = PCEParams[PCEParams.size() - 1];
    else
      m_ParametersCurrent.T = PCEParams[PCEParams.size() - 1];

    m_model->CalculateFeeddown();
    
    m_IsCalculated = true;
  }

  void ThermalModelPCE::PrepareNucleiForPCE(ThermalParticleSystem * TPS)
  {
    ThermalParticleSystem &parts = *TPS;
    vector<long long> nuclpdgs;
    vector<string> nuclnames;
    vector< vector<long long> > nuclcontent;

    nuclpdgs.push_back(1000010020);
    nuclnames.push_back("d");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);

    nuclpdgs.push_back(1000020030);
    nuclnames.push_back("He3");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);

    nuclpdgs.push_back(1010010030);
    nuclnames.push_back("H3La");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1000020040);
    nuclnames.push_back("He4");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);

    nuclpdgs.push_back(1010000020);
    nuclnames.push_back("LambdaNeutron");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1010010020);
    nuclnames.push_back("LambdaProton");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1020000020);
    nuclnames.push_back("DiLambda");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(3122);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1000010030);
    nuclnames.push_back("Triton");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);

    nuclpdgs.push_back(1010020040);
    nuclnames.push_back("He4La");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1010010040);
    nuclnames.push_back("H4La");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1010020050);
    nuclnames.push_back("He5La");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);
    nuclcontent[nuclcontent.size() - 1].push_back(2112);
    nuclcontent[nuclcontent.size() - 1].push_back(3122);

    nuclpdgs.push_back(1020010020);
    nuclnames.push_back("XiProton");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(3322);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);

    nuclpdgs.push_back(1030000020);
    nuclnames.push_back("OmegaProton");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(3334);
    nuclcontent[nuclcontent.size() - 1].push_back(2212);

    nuclpdgs.push_back(1040000020);
    nuclnames.push_back("DiXi0");
    nuclcontent.push_back(vector<long long>());
    nuclcontent[nuclcontent.size() - 1].push_back(3322);
    nuclcontent[nuclcontent.size() - 1].push_back(3322);

    // Fill "decays" of light nuclei
    for (int i = 0; i < parts.Particles().size(); ++i) {
      ThermalParticle &part = parts.Particle(i);

      for (int j = 0; j < nuclpdgs.size(); ++j) {
        if (part.PdgId() == nuclpdgs[j] || part.PdgId() == -nuclpdgs[j]) {
          part.Decays().resize(0);
          long long sign = 1;
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

  vector<int> ThermalModelPCE::ComputePCEStabilityFlags(const ThermalParticleSystem* TPS, bool SahaEquationForNuclei, bool FreezeLongLived, double WidthCut)
  {

    // By default strongly decaying resonances and all nuclei are considered not to be frozen
    vector<int> stability_flags(0);
    for (int i = 0; i < TPS->Particles().size(); ++i) {
      const ThermalParticle& part = TPS->Particles()[i];

      int frozen = 0;

      // Yields of hadrons not decaying strongly are frozen, except light nuclei
      frozen = static_cast<int> (part.DecayType() != ParticleDecayType::Strong);
      if (SahaEquationForNuclei && abs(part.BaryonCharge()) > 1)
        frozen = 0;

      // Yields of long-lived resonances might also be frozen
      if (FreezeLongLived && part.DecayType() == ParticleDecayType::Strong && part.ResonanceWidth() < WidthCut && abs(part.BaryonCharge()) <= 1)
        frozen = 1;

      // The special case of K0. Instead of K0S and K0L work directly with (anti-)K0
      if (part.PdgId() == 310 || part.PdgId() == 130)
        frozen = 0;
      if (part.PdgId() == 311 || part.PdgId() == -311)
        frozen = 1;

      stability_flags.push_back(frozen);
    }
    return stability_flags;
  }

  void ThermalModelPCE::ApplyFixForBoseCondensation()
  {
    for (int ipart = 0; ipart < m_EffectiveCharges.size(); ++ipart) {
      ThermalParticle& part = m_model->TPS()->Particle(ipart);
      part.CalculateAndSetDynamicalThreshold();
      if (part.ZeroWidthEnforced() || part.Statistics() != -1)
        continue;

      double totmu = 0.;
      int nonzerocharges = 0;
      for (int ifeed = 0; ifeed < m_EffectiveCharges[ipart].size(); ++ifeed) {
        totmu += m_EffectiveCharges[ipart][ifeed] * m_model->TPS()->Particle(m_StableMapTo[ifeed]).Mass();
        if (m_EffectiveCharges[ipart][ifeed] != 0.0)
          nonzerocharges++;
      }

      if (totmu > part.DecayThresholdMassDynamical()) {
        cout << "Changing threshold mass for " << part.Name() << " from " << part.DecayThresholdMassDynamical() << " to " << totmu << "\n";
        part.SetDecayThresholdMassDynamical(totmu);
        part.FillCoefficientsDynamical();
      }
    }

    m_model->TPS()->ProcessDecays();
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
    if (m_Mode == 0) {
      const double& Vtmp = x[x.size() - 1];
      m_THM->m_ParametersCurrent.V = Vtmp;
    }
    else {
      const double& Ttmp = x[x.size() - 1];
      m_THM->m_ParametersCurrent.T = Ttmp;
    }
    double V = m_THM->m_ParametersCurrent.V;
    model->SetParameters(m_THM->m_ParametersCurrent);
    //model->CalculateDensities();
    model->CalculatePrimordialDensities();

    for (int is = 0; is < m_THM->m_StableComponentsNumber; ++is) {
      double totdens = 0.;
      for (int i = 0; i < m_THM->m_EffectiveCharges.size(); ++i) {
        totdens += m_THM->m_EffectiveCharges[i][is] * model->Densities()[i];
      }
      
      //ret[is] = totdens * V - m_THM->m_StableDensitiesInit[is] * m_THM->m_ParametersInit.V;
      ret[is] = (totdens * V) / (m_THM->m_StableDensitiesInit[is] * m_THM->m_ParametersInit.V) - 1.;
    }

    ret[ret.size() - 1] = (model->CalculateEntropyDensity() * V) / (m_THM->m_EntropyDensityInit * m_THM->m_ParametersInit.V) - 1.;

    return ret;
  }

} // namespace thermalfist

