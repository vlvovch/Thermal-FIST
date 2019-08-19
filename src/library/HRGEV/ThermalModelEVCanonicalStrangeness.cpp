/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEV/ThermalModelEVCanonicalStrangeness.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"

using namespace std;

namespace thermalfist {

  ThermalModelEVCanonicalStrangeness::ThermalModelEVCanonicalStrangeness(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
    ThermalModelCanonicalStrangeness(TPS_, params)
  {
    m_modelEV = NULL;
    m_v.resize(m_TPS->Particles().size());
    m_TAG = "ThermalModelEVCanonicalStrangeness";

    m_Ensemble = SCE;
    m_InteractionModel = DiagonalEV;
  }

  ThermalModelEVCanonicalStrangeness::~ThermalModelEVCanonicalStrangeness(void)
  {
    ClearModelEV();
  }

  void ThermalModelEVCanonicalStrangeness::CalculateDensitiesGCE() {
    if (m_modelEV == NULL)
      PrepareModelEV();

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_densitiesGCE[i] = /*m_Suppression **/ m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i] - m_v[i] * m_PNS);
    }

    m_GCECalculated = true;
  }

  void ThermalModelEVCanonicalStrangeness::CalculateEnergyDensitiesGCE() {
    if (m_modelEV == NULL)
      PrepareModelEV();

    m_energydensitiesGCE.resize(m_densitiesGCE.size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_energydensitiesGCE[i] = /*m_Suppression **/ m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i] - m_v[i] * m_PNS);
    }
  }

  void ThermalModelEVCanonicalStrangeness::CalculatePressuresGCE()
  {
    if (m_modelEV == NULL)
      PrepareModelEV();

    m_pressuresGCE.resize(m_densitiesGCE.size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_pressuresGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i] - m_v[i] * m_PNS);
    }
  }


  void ThermalModelEVCanonicalStrangeness::CalculatePrimordialDensities() {
    m_FluctuationsCalculated = false;

    m_energydensitiesGCE.resize(0);
    m_pressuresGCE.resize(0);

    PrepareModelEV();

    CalculateDensitiesGCE();

    CalculateSums(m_Parameters.SVc - m_EVNS);

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
        m_densities[i] = m_Suppression * (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_densitiesGCE[i];
    }

    m_EVS = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      if (m_TPS->Particles()[i].Strangeness() != 0)
        m_EVS += m_v[i] * m_densities[i] * m_Parameters.SVc;
    }

    //printf("%15lf%15lf\n", CalculatePressure(), m_PNS);
    double tP = CalculatePressure();
    printf("EV-SCE calculation: The following two parameters must be much smaller than unity. Otherwise we are in trouble...\n");
    printf("%20s%lf\n%20s%lf\n", "PS/P = ", (tP - m_PNS) / tP, "EVS/(V-EVNS) = ", m_EVS / (m_Parameters.SVc - m_EVNS));

    m_Calculated = true;
    ValidateCalculation();
    //m_LastCalculationSuccessFlag = true;
  }



  double ThermalModelEVCanonicalStrangeness::CalculateEnergyDensity() {
    if (!m_Calculated)
      CalculateDensities();

    if (m_energydensitiesGCE.size() == 0)
      CalculateEnergyDensitiesGCE();

    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_Suppression * (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_energydensitiesGCE[i];
    return ret;
  }


  double ThermalModelEVCanonicalStrangeness::CalculateEntropyDensity() {
    double ret = log(m_Zsum[m_StrMap[0]]) / m_Parameters.SVc;

    if (m_energydensitiesGCE.size() == 0)
      CalculateEnergyDensitiesGCE();

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].Strangeness() != 0) {
        if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
          ret += m_Suppression * (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * ((m_energydensitiesGCE[i] - (m_Chem[i] - m_v[i] * m_PNS) * m_densitiesGCE[i]) / m_Parameters.T);
      }
      else {
        ret += m_Suppression * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i] - m_v[i] * m_PNS);
      }
    return ret;
  }


  double ThermalModelEVCanonicalStrangeness::CalculatePressure() {
    double ret = 0.;
    if (m_pressuresGCE.size() == 0)
      CalculatePressuresGCE();
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_pressuresGCE[i];
    return ret;
  }

  double ThermalModelEVCanonicalStrangeness::MuShift(int id) const
  {
    if (id >= 0. && id < static_cast<int>(m_v.size()))
      return -m_v[id] * m_PNS;
    else
      return 0.0;
  }


  void ThermalModelEVCanonicalStrangeness::PrepareModelEV()
  {
    ClearModelEV();

    ThermalParticleSystem *TPSnew = new ThermalParticleSystem(*m_TPS);

    // "Switch off" all strange particles
    for (size_t i = 0; i < TPSnew->Particles().size(); ++i) {
      ThermalParticle &part = TPSnew->Particle(i);
      if (part.Strangeness() != 0)
        part.SetDegeneracy(0.);
    }

    m_modelEV = new ThermalModelEVDiagonal(TPSnew);
    m_modelEV->FillVirialEV(m_v);
    m_modelEV->SetUseWidth(m_UseWidth);
    m_modelEV->SetParameters(m_Parameters);
    m_modelEV->SetVolume(m_Parameters.SVc);
    m_modelEV->SetChemicalPotentials(m_Chem);

    m_modelEV->CalculateDensities();

    m_PNS = m_modelEV->CalculatePressure();
    m_Suppression = m_modelEV->CommonSuppressionFactor();
    m_EVNS = 0.;
    for (size_t i = 0; i < m_modelEV->Densities().size(); ++i) {
      m_EVNS += m_v[i] * m_modelEV->Densities()[i] * m_Parameters.SVc;
    }
  }

  void ThermalModelEVCanonicalStrangeness::ClearModelEV()
  {
    if (m_modelEV != NULL) {
      ThermalParticleSystem *TPSold = m_modelEV->TPS();
      if (TPSold != NULL)
        delete TPSold;
      delete m_modelEV;
      m_modelEV = NULL;
    }
  }

  void ThermalModelEVCanonicalStrangeness::SetRadius(double rad)
  {
    if (static_cast<int>(m_v.size()) != m_TPS->ComponentsNumber())
      m_v.resize(m_TPS->Particles().size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_v[i] = CuteHRGHelper::vr(rad);
    }
  }

  void ThermalModelEVCanonicalStrangeness::FillVirial(const std::vector<double>& ri)
  {
    if (ri.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillVirial(const std::vector<double> & ri): size of ri does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_v.resize(m_TPS->Particles().size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      m_v[i] = CuteHRGHelper::vr(ri[i]);
  }

  void ThermalModelEVCanonicalStrangeness::FillVirialEV(const std::vector<double>& vi)
  {
    if (vi.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillVirialEV(const std::vector<double> & vi): size of vi does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_v = vi;
  }

  void ThermalModelEVCanonicalStrangeness::ReadInteractionParameters(const std::string & filename)
  {
    m_v = std::vector<double>(m_TPS->Particles().size(), 0.);

    ifstream fin(filename.c_str());
    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1)
        continue;
      istringstream iss(elems[0]);
      int pdgid;
      double b;
      if (iss >> pdgid >> b) {
        int ind = m_TPS->PdgToId(pdgid);
        if (ind != -1)
          m_v[ind] = b;
      }
    }
    fin.close();
  }

  void ThermalModelEVCanonicalStrangeness::WriteInteractionParameters(const std::string & filename)
  {
    ofstream fout(filename.c_str());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      fout << std::setw(15) << m_TPS->Particle(i).PdgId();
      fout << std::setw(15) << m_v[i];
      fout << std::endl;
    }
    fout.close();
  }

  double ThermalModelEVCanonicalStrangeness::ExcludedVolume(int i) const
  {
    if (i < 0 || i >= static_cast<int>(m_v.size()))
      return 0.;
    return m_v[i];
  }

  double ThermalModelEVCanonicalStrangeness::CalculateEigenvolumeFraction()
  {
    double tEV = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      tEV += m_v[i] * m_densities[i] / 4.;
    }
    return tEV;
  }

  void ThermalModelEVCanonicalStrangeness::SetRadius(int i, double rad)
  {
    if (i >= 0 && i < static_cast<int>(m_v.size()))
      m_v[i] = CuteHRGHelper::vr(rad);
  }

  bool ThermalModelEVCanonicalStrangeness::IsConservedChargeCanonical(ConservedCharge::Name charge) const {
    return (charge == ConservedCharge::StrangenessCharge);
  }

} // namespace thermalfist