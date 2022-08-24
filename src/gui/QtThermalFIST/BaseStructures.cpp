#include "BaseStructures.h"

#include "HRGEV/ExcludedVolumeHelper.h"
#include "HRGBase/ThermalModelCanonical.h"

using namespace thermalfist;

ThermalModelConfig ThermalModelConfig::fromThermalModel(ThermalModelBase * model)
{
  ThermalModelConfig ret;
  
  if (model->InteractionModel() == ThermalModelBase::Ideal) {
    if (model->Ensemble() == ThermalModelBase::CE)
      ret.ModelType = ThermalModelConfig::CE;
    if (model->Ensemble() == ThermalModelBase::SCE)
      ret.ModelType = ThermalModelConfig::SCE;
    else if (model->Ensemble() == ThermalModelBase::CCE)
      ret.ModelType = ThermalModelConfig::CCE;
    else
      ret.ModelType = ThermalModelConfig::Ideal;
  }
  else if (model->InteractionModel() == ThermalModelBase::DiagonalEV) {
    if (model->Ensemble() == ThermalModelBase::SCE)
      ret.ModelType = ThermalModelConfig::EVSCE;
    else
      ret.ModelType = ThermalModelConfig::DiagonalEV;
  }
  else if (model->InteractionModel() == ThermalModelBase::CrosstermsEV) {
    ret.ModelType = ThermalModelConfig::CrosstermsEV;
  }
  else if (model->InteractionModel() == ThermalModelBase::QvdW) {
    if (model->Ensemble() == ThermalModelBase::SCE)
      ret.ModelType = ThermalModelConfig::VDWSCE;
    else
      ret.ModelType = ThermalModelConfig::QvdW;
  }

  ret.Ensemble = model->Ensemble();

  ret.InteractionModel = model->InteractionModel();

  ret.QuantumStatistics = model->QuantumStatistics();

  ret.QuantumStatisticsType = model->TPS()->QStatsCalculationType();

  ret.QuantumStatisticsInclude = 0;

  ret.InteractionScaling = 0;

  ret.vdWA = 0.;

  ret.vdWB = 1.;

  ret.InteractionInput = "";

  ret.DisableMM = 1;

  ret.DisableMB = 1;

  ret.DisableBB = 0;

  ret.DisableBantiB = 1;

  ret.T = model->Parameters().T;

  ret.muB = model->Parameters().muB;

  ret.muQ = model->Parameters().muQ;

  ret.muS = model->Parameters().muS;

  ret.muC = model->Parameters().muC;

  ret.gq = model->Parameters().gammaq;

  ret.gS = model->Parameters().gammaS;

  ret.gC = model->Parameters().gammaC;

  ret.VolumeR = pow(3. / 16. / xMath::Pi() * model->Parameters().V, 1./3.);

  ret.VolumeRSC = pow(3. / 16. / xMath::Pi() * model->Parameters().SVc, 1. / 3.);

  ret.B = model->Parameters().B;

  ret.Q = model->Parameters().Q;

  ret.S = model->Parameters().S;

  ret.C = model->Parameters().C;

  ret.CanonicalB = true;// model->IsConservedChargeCanonical(ConservedCharge::BaryonCharge);

  ret.CanonicalQ = true;// model->IsConservedChargeCanonical(ConservedCharge::ElectricCharge);

  ret.CanonicalS = true;// model->IsConservedChargeCanonical(ConservedCharge::StrangenessCharge);

  ret.CanonicalC = true;// model->IsConservedChargeCanonical(ConservedCharge::CharmCharge);

  ret.SoverB = model->SoverB();

  ret.QoverB = model->QoverB();

  ret.ConstrainMuB = model->ConstrainMuB();

  ret.ConstrainMuQ = model->ConstrainMuQ();

  ret.ConstrainMuS = model->ConstrainMuS();

  ret.ConstrainMuC = model->ConstrainMuC();

  ret.FiniteWidth = static_cast<int>(model->TPS()->ResonanceWidthIntegrationType());

  ret.WidthShape = static_cast<int>(model->TPS()->ResonanceWidthShape());

  ret.RenormalizeBR = model->NormBratio();

  ret.ComputeFluctations = false;


  ret.UsePCE = model->UsePartialChemicalEquilibrium();
  ret.Tkin = 0.100;
  ret.PCEFreezeLongLived = false;
  ret.PCEWidthCut = 0.015;
  ret.PCESahaForNuclei = true;

  ret.fUseEVRejectionMultiplicity = true;
  ret.fUseEVRejectionCoordinates  = true;
  ret.fUseEVUseSPRApproximation   = true;

  return ret;
}

void SetThermalModelParameters(thermalfist::ThermalModelBase* model, const ThermalModelConfig& config)
{
  model->SetTemperature(config.T);
  model->SetBaryonChemicalPotential(config.muB);
  model->SetElectricChemicalPotential(config.muQ);
  model->SetStrangenessChemicalPotential(config.muS);
  model->SetCharmChemicalPotential(config.muC);
  model->SetGammaq(config.gq);
  model->SetGammaS(config.gS);
  model->SetGammaC(config.gC);
  model->SetVolumeRadius(config.VolumeR);
  model->SetCanonicalVolumeRadius(config.VolumeRSC);
}

void SetThermalModelConfiguration(thermalfist::ThermalModelBase * model, const ThermalModelConfig & config)
{
  model->SetBaryonCharge(config.B);
  model->SetElectricCharge(config.Q);
  model->SetStrangeness(config.S);
  model->SetCharm(config.C);

  if (model->Ensemble() == ThermalModelBase::CE) {
    ThermalModelCanonical *modcan = static_cast<ThermalModelCanonical*>(model);
    modcan->ConserveBaryonCharge(config.CanonicalB);
    modcan->ConserveElectricCharge(config.CanonicalQ);
    modcan->ConserveStrangeness(config.CanonicalS);
    modcan->ConserveCharm(config.CanonicalC);

    //modcan->SetIntegrationIterationsMultiplier(1);
  }

  if (config.WidthShape == 0)
    model->SetResonanceWidthShape(ThermalParticle::RelativisticBreitWigner);
  else
    model->SetResonanceWidthShape(ThermalParticle::NonRelativisticBreitWigner);

  if (config.FiniteWidth == 0)
    model->SetUseWidth(ThermalParticle::ZeroWidth);
  else if (config.FiniteWidth == 1)
    model->SetUseWidth(ThermalParticle::BWTwoGamma);
  else if (config.FiniteWidth == 2)
    model->SetUseWidth(ThermalParticle::eBW);
  else if (config.FiniteWidth == 3)
    model->SetUseWidth(ThermalParticle::eBWconstBR);
  else
    model->SetUseWidth(ThermalParticle::ZeroWidth);

  model->SetSoverB(config.SoverB);
  model->ConstrainMuB(config.ConstrainMuB);
  model->SetQoverB(config.QoverB);
  model->ConstrainMuQ(config.ConstrainMuQ);
  model->ConstrainMuS(config.ConstrainMuS);
  model->ConstrainMuC(config.ConstrainMuC);

  model->SetNormBratio(config.RenormalizeBR);

  model->SetStatistics(config.QuantumStatistics);
  //model->SetClusterExpansionOrder(10);

  if (config.QuantumStatisticsType)
    model->SetCalculationType(IdealGasFunctions::Quadratures);
  else
    model->SetCalculationType(IdealGasFunctions::ClusterExpansion);

  if (config.QuantumStatisticsInclude == 1 || config.QuantumStatisticsInclude == 2) {
    for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
      ThermalParticle &part = model->TPS()->Particle(i);
      if (config.QuantumStatisticsInclude == 2) {
        if (part.PdgId() != 211 && part.PdgId() != 111 && part.PdgId() != -211) {
          part.UseStatistics(false);
        }
      }
      else if (config.QuantumStatisticsInclude == 1) {
        if (part.BaryonCharge() != 0) {
          part.UseStatistics(false);
        }
      }
    }
  }
}

void SetThermalModelInteraction(ThermalModelBase * model, const ThermalModelConfig & config)
{
  if (config.InteractionScaling != 3) {
    // First repulsion
    double radius = CuteHRGHelper::rv(config.vdWB);
    std::vector<double> radii(model->TPS()->Particles().size(), 0.);
    // Uniform EV
    if (config.InteractionScaling == 0) {
      model->SetRadius(radius);
      std::fill(radii.begin(), radii.end(), radius);
    }

    // Bag model EV, nucleon mass here equals 0.938 GeV
    if (config.InteractionScaling == 1) {
      for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
        ThermalParticle &part = model->TPS()->Particle(i);
        radii[i] = radius * pow(part.Mass() / xMath::mnucleon(), 1. / 3.);
      }
      model->FillVirial(radii);
    }

    // Linear in baryon content
    if (config.InteractionScaling == 2) {
      for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
        ThermalParticle &part = model->TPS()->Particle(i);
        if (part.BaryonCharge() != 0)
          radii[i] = radius * pow(abs(part.BaryonCharge()), 1. / 3.);
        else
          radii[i] = 0.;
      }
      model->FillVirial(radii);
    }

    model->FillVirial(radii); // Just in case

    // Now attraction

    if (config.InteractionModel == ThermalModelConfig::InteractionQVDW) {
      std::vector<double> as(model->TPS()->Particles().size(), config.vdWA);

      // Mass-proportional, nucleon mass here equals 0.938 GeV
      if (config.InteractionScaling == 1) {
        for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
          ThermalParticle& part = model->TPS()->Particle(i);
          as[i] = config.vdWA * part.Mass() / xMath::mnucleon();
        }
      }

      // Linear in baryon content
      if (config.InteractionScaling == 2) {
        for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
          ThermalParticle& part = model->TPS()->Particle(i);
          if (part.BaryonCharge() != 0)
            as[i] = config.vdWA * abs(part.BaryonCharge());
          else
            as[i] = 0.;
        }
      }

      // Fill aij assuming "chemistry rule" aij = \sqrt{ai * aj}
      for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
        for (int j = 0; j < model->TPS()->Particles().size(); ++j) {
          model->SetAttraction(i, j, sqrt(as[i] * as[j]));
        }
      }
    }

    if (config.DisableMM) {
      model->DisableMesonMesonRepulsion();
      model->DisableMesonMesonAttraction();
    }

    if (config.DisableMB) {
      model->DisableMesonBaryonRepulsion();
      model->DisableMesonBaryonAttraction();
    }

    if (config.DisableBB) {
      model->DisableBaryonBaryonRepulsion();
      model->DisableBaryonBaryonAttraction();
    }

    if (config.DisableBantiB) {
      model->DisableBaryonAntiBaryonRepulsion();
      model->DisableBaryonAntiBaryonAttraction();
    }
  }

  // Read from file
  if (config.InteractionScaling == 3) {
    model->ReadInteractionParameters(config.InteractionInput);
  }
}
