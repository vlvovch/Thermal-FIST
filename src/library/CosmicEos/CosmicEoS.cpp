#include "CosmicEos/CosmicEoS.h"

#include <iostream>

using namespace std;

namespace thermalfist {

  /// Pion decay constant
  double CosmicEoS::fpi = 0.133;

  /// Electron mass
  double LeptonFlavor::m_e   = 0.000511;

  /// Muon mass
  double LeptonFlavor::m_mu  = 0.105658;

  /// Tauon mass
  double LeptonFlavor::m_tau = 1.77684;

  CosmicEoS::CosmicEoS(ThermalModelBase* THMbase, bool pionsinteract) :
    m_modelHRG(THMbase)
  {
    // No strangeness of charm conservation
    m_modelHRG->SetStrangenessChemicalPotential(0.);
    m_modelHRG->SetCharmChemicalPotential(0.);

    m_ChemCurrent = std::vector<double>(2 + LeptonFlavor::NumberOfFlavors, 0.);
    m_Asymmetries = std::vector<double>(2 + LeptonFlavor::NumberOfFlavors, 0.);

    // Default baryon asymmetry
    double b = 8.6e-11;

    // Default lepton asymmetry (uniform across the three flavors)
    double l = -(51. / 28.) * b;
    m_Asymmetries[0] = b;
    m_Asymmetries[1] = 0.;
    m_Asymmetries[2] = l / 3.;
    m_Asymmetries[3] = l / 3.;
    m_Asymmetries[4] = l / 3.;

    // Photon
    m_Photon = ThermalParticle(true, "photon", 22, 2., -1, 0.);
    
    // Charged leptons
    m_ChargedLeptons.push_back(ThermalParticle(true, "e", 11, 2., 1, LeptonFlavor::m_e, 0, 0, -1));
    m_ChargedLeptons.push_back(ThermalParticle(true, "mu", 13, 2., 1, LeptonFlavor::m_mu, 0, 0, -1));
    m_ChargedLeptons.push_back(ThermalParticle(true, "tau", 15, 2., 1, LeptonFlavor::m_tau, 0, 0, -1));

    // Neutrinos
    m_Neutrinos.push_back(ThermalParticle(true, "nu_e", 12, 1., 1, 0.));
    m_Neutrinos.push_back(ThermalParticle(true, "nu_m", 14, 1., 1, 0.));
    m_Neutrinos.push_back(ThermalParticle(true, "nu_t", 16, 1., 1, 0.));

    SetPionsInteracting(pionsinteract);
  }

  void CosmicEoS::SetTemperature(double T)
  {
    m_T = T; 
    m_modelHRG->SetTemperature(T);
  }

  void CosmicEoS::SetElectricChemicalPotential(double muQ)
  {
    m_ChemCurrent[1] = muQ; 
    m_modelHRG->SetElectricChemicalPotential(muQ);
  }

  void CosmicEoS::CalculatePrimordialDensities()
  {
    m_modelHRG->CalculatePrimordialDensities();
  }

  double CosmicEoS::EntropyDensity()
  {
    double ret = EntropyDensityHRG();

    ret += m_Photon.Density(m_modelHRG->Parameters(), IdealGasFunctions::EntropyDensity, false, 0.0);

    double muB = m_ChemCurrent[0];
    double muQ = m_ChemCurrent[1];
    for (int iL = 0; iL < 3; ++iL) {
      ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EntropyDensity, false, m_ChemCurrent[2 + iL] - muQ);
      ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EntropyDensity, false, -(m_ChemCurrent[2 + iL] - muQ));
      ret += m_Neutrinos[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EntropyDensity, false, m_ChemCurrent[2 + iL]);
      ret += m_Neutrinos[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EntropyDensity, false, -m_ChemCurrent[2 + iL]);
    }

    return ret;
  }

  double CosmicEoS::EntropyDensityHRG()
  {
    double ret = m_modelHRG->EntropyDensity();

    return ret;
  }

  double CosmicEoS::Pressure()
  {
    //double ret = m_modelHRG->Pressure();
    double ret = PressureHRG();

    ret += m_Photon.Density(m_modelHRG->Parameters(), IdealGasFunctions::Pressure, false, 0.0);

    double muB = m_ChemCurrent[0];
    double muQ = m_ChemCurrent[1];
    for (int iL = 0; iL < 3; ++iL) {
      ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::Pressure, false, m_ChemCurrent[2 + iL] - muQ);
      ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::Pressure, false, -(m_ChemCurrent[2 + iL] - muQ));
      ret += m_Neutrinos[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::Pressure, false, m_ChemCurrent[2 + iL]);
      ret += m_Neutrinos[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::Pressure, false, -m_ChemCurrent[2 + iL]);
    }

    return ret;
  }

  double CosmicEoS::PressureHRG()
  {
    double ret = m_modelHRG->Pressure();

    return ret;
  }

  double CosmicEoS::EnergyDensity()
  {
    double ret = EnergyDensityHRG();

    ret += m_Photon.Density(m_modelHRG->Parameters(), IdealGasFunctions::EnergyDensity, false, 0.0);

    double muB = m_ChemCurrent[0];
    double muQ = m_ChemCurrent[1];
    for (int iL = 0; iL < 3; ++iL) {
      ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EnergyDensity, false, m_ChemCurrent[2 + iL] - muQ);
      ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EnergyDensity, false, -(m_ChemCurrent[2 + iL] - muQ));
      ret += m_Neutrinos[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EnergyDensity, false, m_ChemCurrent[2 + iL]);
      ret += m_Neutrinos[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EnergyDensity, false, -m_ChemCurrent[2 + iL]);
    }

    return ret;
  }

  double CosmicEoS::EnergyDensityHRG()
  {
    double ret = m_modelHRG->EnergyDensity();

    return ret;
  }

  double CosmicEoS::NetDensityChargedLepton(int iL)
  {
    double ret = 0.0;

    double muQ = m_ChemCurrent[1];
    ret += (-1.) * m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, m_ChemCurrent[2 + iL] - muQ);
    ret += ( 1.) * m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, -(m_ChemCurrent[2 + iL] - muQ));

    return ret;
  }

  double CosmicEoS::PressureChargedLepton(int iL)
  {
    double ret = 0.0;
    double muQ = m_ChemCurrent[1];
    ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::Pressure, false, m_ChemCurrent[2 + iL] - muQ);
    ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::Pressure, false, -(m_ChemCurrent[2 + iL] - muQ));
    return ret;
  }

  double CosmicEoS::EnergyDensityChargedLepton(int iL)
  {
    double ret = 0.0;
    double muQ = m_ChemCurrent[1];
    ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EnergyDensity, false, m_ChemCurrent[2 + iL] - muQ);
    ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::EnergyDensity, false, -(m_ChemCurrent[2 + iL] - muQ));
    return ret;
  }

  double CosmicEoS::BaryonDensity(bool absolute)
  {
    if (!absolute) return m_modelHRG->BaryonDensity();
    else return m_modelHRG->AbsoluteBaryonDensity();
  }

  double CosmicEoS::ElectricChargeDensity(bool absolute)
  {
    double ret = ElectricChargeDensityHRG(absolute);

    double mn = -1.;
    if (absolute)
      mn = 1.;

    double muB = m_ChemCurrent[0];
    double muQ = m_ChemCurrent[1];
    for (int iL = 0; iL < 3; ++iL) {
      ret += mn * m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, m_ChemCurrent[2 + iL] - muQ);
      ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, -(m_ChemCurrent[2 + iL] - muQ));
    }

    return ret;
  }

  double CosmicEoS::ElectricChargeDensityHRG(bool absolute)
  {
    double ret = 0.;
    if (!absolute) ret = m_modelHRG->ElectricChargeDensity();
    else ret = m_modelHRG->AbsoluteElectricChargeDensity();

    return ret;
  }

  double CosmicEoS::IsospinChargeDensity(bool absolute)
  {
    double nI = 0.0;
    for (int ipart = 0; ipart < HRGModel()->TPS()->Particles().size(); ++ipart) {
      double iso_chg = IsospinCharge(HRGModel()->TPS()->Particle(ipart));
      if (absolute)
        iso_chg = abs(iso_chg);
      nI += iso_chg * HRGModel()->Densities()[ipart];
    }

    return nI;
  }

  double CosmicEoS::LeptonFlavorDensity(LeptonFlavor::Name flavor, bool absolute)
  {
    double ret = 0.0;

    double mn = -1.;
    if (absolute)
      mn = 1.;

    double muQ = m_ChemCurrent[1];
    {
      int iL = static_cast<int>(flavor);
      ret += m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, m_ChemCurrent[2 + iL] - muQ);
      ret += mn * m_ChargedLeptons[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, -(m_ChemCurrent[2 + iL] - muQ));
      ret += m_Neutrinos[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, m_ChemCurrent[2 + iL]);
      ret += mn * m_Neutrinos[iL].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, -m_ChemCurrent[2 + iL]);
    }

    return ret;
  }

  std::vector<double> CosmicEoS::SolveChemicalPotentials(double T, const std::vector<double>& muInit)
  {
    std::vector<double> ret = m_ChemCurrent;

    if (muInit.size() == m_Asymmetries.size())
      m_ChemCurrent = muInit;
    else
      m_ChemCurrent = std::vector<double>(m_Asymmetries.size(), 0);

    SetTemperature(T);

    BroydenEquationsCosmology eqs(this);
    Broyden broydn(&eqs);
    //broydn.UseNewton(true);

    Broyden::BroydenSolutionCriterium criterium(1.e-7);

    return broydn.Solve(m_ChemCurrent, &criterium);
  }

  double CosmicEoS::GetPionMass() const
  {
    int pionid = m_modelHRG->TPS()->PdgToId(211);
    if (pionid != -1)
      return m_modelHRG->TPS()->Particle(pionid).Mass();
    return 0.138;
  }

  void CosmicEoS::SetPionsInteracting(bool pionsinteract, double fpiChPT)
  {
    m_InteractingPions = pionsinteract;
    ClearEMMs();
    if (pionsinteract) {
      HRGModel()->TPS()->ParticleByPDG(211).SetGeneralizedDensity(
              new EffectiveMassModel(
                      HRGModel()->TPS()->ParticleByPDG(211),
                      new EMMFieldPressureChPT(HRGModel()->TPS()->ParticleByPDG(211).Mass(), fpiChPT)
                      ));
      HRGModel()->TPS()->ParticleByPDG(111).SetGeneralizedDensity(
              new EffectiveMassModel(
                      HRGModel()->TPS()->ParticleByPDG(111),
                      new EMMFieldPressureChPT(HRGModel()->TPS()->ParticleByPDG(111).Mass(), fpiChPT)
              ));
      HRGModel()->TPS()->ParticleByPDG(-211).SetGeneralizedDensity(
              new EffectiveMassModel(
                      HRGModel()->TPS()->ParticleByPDG(-211),
                      new EMMFieldPressureChPT(HRGModel()->TPS()->ParticleByPDG(-211).Mass(), fpiChPT)
              ));
    }
  }

  bool CosmicEoS::InPionCondensedPhase() const
  {
    if (!InteractingPions())
      return abs(ElectricChemicalPotential()) >= GetPionMass();
    else {
      return HRGModel()->TPS()->ParticleByPDG(211).GetGeneralizedDensity()->IsBECPhase()
      || HRGModel()->TPS()->ParticleByPDG(111).GetGeneralizedDensity()->IsBECPhase()
      || HRGModel()->TPS()->ParticleByPDG(-211).GetGeneralizedDensity()->IsBECPhase();
    }
    return false;
  }

  std::string CosmicEoS::GetSpeciesName(int id) const
  {
    if (id == 0)
      return m_Photon.Name();
    else {
      int iL = id - 1;
      if (iL/2 < m_ChargedLeptons.size())
        if (iL % 2 == 0)
          return m_ChargedLeptons[iL/2].Name();
        else
          return "anti-" + m_ChargedLeptons[iL/2].Name();
      else {
        iL -= 2 * m_ChargedLeptons.size();
        if (iL/2 < m_Neutrinos.size())
          if (iL % 2 == 0)
            return m_Neutrinos[iL/2].Name();
          else
            return "anti-" + m_Neutrinos[iL/2].Name();
        else {
          printf("**ERRORR** in CosmicEoS::GetSpeciesName(int id) const: id = %d is out of range!", id);
        }
      }
    }
  }

  double CosmicEoS::GetDensity(int id) const
  {
    if (id == 0)
      return m_Photon.Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, 0.);
    else {
      double muQ = m_ChemCurrent[1];
      int iL = id - 1;
      if (iL/2 < m_ChargedLeptons.size())
        if (iL % 2 == 0)
          return m_ChargedLeptons[iL/2].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, m_ChemCurrent[2 + iL/2] - muQ);
        else
          return m_ChargedLeptons[iL/2].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, -(m_ChemCurrent[2 + iL/2] - muQ));
      else {
        iL -= 2 * m_ChargedLeptons.size();
        if (iL/2 < m_Neutrinos.size())
          if (iL % 2 == 0)
            return m_Neutrinos[iL/2].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, m_ChemCurrent[2 + iL/2]);
          else
            return m_Neutrinos[iL/2].Density(m_modelHRG->Parameters(), IdealGasFunctions::ParticleDensity, false, -m_ChemCurrent[2 + iL/2]);
        else
          printf("**ERRORR** in CosmicEoS::GetDensity(int id) const: id = %d is out of range!", id);
      }
    }
  }

  void CosmicEoS::ClearEMMs()
  {
    for (auto& part :m_modelHRG->TPS()->Particles()) {
      part.ClearGeneralizedDensity();
    }
  }

  std::vector<double> CosmicEoS::BroydenEquationsCosmology::Equations(const std::vector<double>& x)
  {
    std::vector<double> ret(x.size(), 0.);

    double muB = x[0];
    double muQ = x[1];

    m_THM->SetBaryonChemicalPotential(muB);
    m_THM->SetElectricChemicalPotential(muQ);
    for (int iL = 0; iL < LeptonFlavor::NumberOfFlavors; iL++)
    {
      auto flavor = static_cast<LeptonFlavor::Name>(iL);
      m_THM->SetLeptonChemicalPotential(flavor, x[2 + iL]);
    }

    m_THM->CalculatePrimordialDensities();


    vector<double> constraints = m_THM->m_Asymmetries;

    double s = m_THM->EntropyDensity();
    double nB = m_THM->BaryonDensity();

    if (constraints[0] != 0.0)
      ret[0] = (nB / s - constraints[0]) / constraints[0];
    else
      ret[0] = nB / m_THM->BaryonDensity(true);

    double nQ = m_THM->ElectricChargeDensity();

    if (constraints[1] != 0.0)
      ret[1] = (nQ / s - constraints[1]) / constraints[1];
    // else if (constraints[0] != 0.0)
    //   ret[1] = nQ / nB;
    else  
      ret[1] = nQ / m_THM->ElectricChargeDensity(true);

    vector<double> nLs(LeptonFlavor::NumberOfFlavors, 0.);

    for (int iL = 0; iL < LeptonFlavor::NumberOfFlavors; iL++)
    {
      auto flavor = static_cast<LeptonFlavor::Name>(iL);
      nLs[iL] = m_THM->LeptonFlavorDensity(flavor);

      if (constraints[2 + iL] != 0.0)
        ret[2 + iL] = (nLs[iL] / s - constraints[2 + iL]) / constraints[2 + iL];
      else
        ret[2 + iL] = nLs[iL] / m_THM->LeptonFlavorDensity(flavor, true);
    }

    return ret;
  }

  double IsospinCharge(const ThermalParticle& part)
  {
    return 0.5 * (2. * part.ElectricCharge() - part.BaryonCharge() - part.Strangeness() - part.Charm());
  }
}