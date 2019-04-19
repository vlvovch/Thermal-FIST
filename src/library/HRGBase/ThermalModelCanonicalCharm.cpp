/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ThermalModelCanonicalCharm.h"

#include "HRGBase/xMath.h"


using namespace std;

namespace thermalfist {

  ThermalModelCanonicalCharm::ThermalModelCanonicalCharm(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
    ThermalModelBase(TPS_, params)
  {
    m_densitiesGCE.resize(m_TPS->Particles().size());

    m_CharmValues.resize(0);
    m_CharmValues.push_back(0);
    m_CharmValues.push_back(1);
    m_CharmValues.push_back(-1);

    m_CharmMap.clear();
    for (unsigned int i = 0; i < m_CharmValues.size(); ++i) m_CharmMap[m_CharmValues[i]] = i;

    m_Parameters.muC = 0.;

    m_TAG = "ThermalModelCanonicalCharm";

    m_Ensemble = CCE;
    m_InteractionModel = Ideal;

  }

  ThermalModelCanonicalCharm::~ThermalModelCanonicalCharm(void)
  {
  }


  void ThermalModelCanonicalCharm::SetParameters(const ThermalModelParameters& params) {
    ThermalModelBase::SetParameters(params);
    m_Parameters.muC = 0.;
  }

  void ThermalModelCanonicalCharm::SetCharmChemicalPotential(double /*muC*/)
  {
    m_Parameters.muC = 0.0;
  }

  void ThermalModelCanonicalCharm::ChangeTPS(ThermalParticleSystem *TPS_) {
    ThermalModelBase::ChangeTPS(TPS_);
    m_densitiesGCE.resize(m_TPS->Particles().size());
  }

  void ThermalModelCanonicalCharm::SetStatistics(bool stats) {
    m_QuantumStats = stats;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].Charm() == 0) m_TPS->Particle(i).UseStatistics(stats);
      else m_TPS->Particle(i).UseStatistics(false);
  }

  void ThermalModelCanonicalCharm::CalculateDensitiesGCE() {
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_densitiesGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
    }
    m_GCECalculated = true;
  }

  void ThermalModelCanonicalCharm::CalculateEnergyDensitiesGCE() {
    m_energydensitiesGCE.resize(m_densitiesGCE.size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_energydensitiesGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);
    }
  }

  void ThermalModelCanonicalCharm::FixParameters()
  {
    m_ConstrainMuC = false;
    ThermalModelBase::FixParameters();
  }

  void ThermalModelCanonicalCharm::CalculatePrimordialDensities() {
    if (UsePartialChemicalEquilibrium()) {
      printf("**ERROR** ThermalModelCanonicalCharm::CalculatePrimordialDensities(): PCE not supported!\n");
      exit(1);
    }

    m_FluctuationsCalculated = false;
    m_energydensitiesGCE.resize(0);

    CalculateDensitiesGCE();

    m_Zsum.resize(m_CharmValues.size());

    m_partialZ.resize(m_CharmValues.size());
    vector<double> xi(1, 0.), yi(1, 0.);

    for (size_t i = 0; i < m_CharmValues.size(); ++i) {
      m_partialZ[i] = 0.;
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j)
        if (m_CharmValues[i] == m_TPS->Particles()[j].Charm()) m_partialZ[i] += m_densitiesGCE[j] * m_Volume;
      if (m_partialZ[i] < 1.e-10) m_partialZ[i] += 1e-10;
    }


    for (int i = 0; i < 1; ++i) {
      xi[i] = 2. * sqrt(m_partialZ[m_CharmMap[i + 1]] * m_partialZ[m_CharmMap[-(i + 1)]]);
      yi[i] = sqrt(m_partialZ[m_CharmMap[i + 1]] / m_partialZ[m_CharmMap[-(i + 1)]]);
    }

    for (size_t i = 0; i < m_CharmValues.size(); ++i) {
      double res = 0.;

      res = xMath::BesselI(abs(-m_CharmValues[i]), xi[0]) * pow(yi[0], m_CharmValues[i]);
      m_Zsum[i] = res;
    }

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      if (m_CharmMap.count(-m_TPS->Particles()[i].Charm())) m_densities[i] = (m_Zsum[m_CharmMap[-m_TPS->Particles()[i].Charm()]] / m_Zsum[m_CharmMap[0]]) * m_densitiesGCE[i];
    }

    m_Calculated = true;
    ValidateCalculation();
  }

  void ThermalModelCanonicalCharm::CalculateFluctuations() {
    m_wprim.resize(m_densities.size());
    m_wtot.resize(m_densities.size());
    for (size_t i = 0; i < m_wprim.size(); ++i) {
      m_wprim[i] = 1.;
      m_wtot[i] = 1.;
      m_skewprim[i] = 1.;
      m_skewtot[i] = 1.;
    }
  }

  double ThermalModelCanonicalCharm::CalculateEnergyDensity() {
    if (!m_Calculated) CalculateDensities();
    if (m_energydensitiesGCE.size() == 0) CalculateEnergyDensitiesGCE();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) ret += (m_Zsum[m_CharmMap[-m_TPS->Particles()[i].Charm()]] / m_Zsum[m_CharmMap[0]]) * m_energydensitiesGCE[i];
    return ret;
  }

  double ThermalModelCanonicalCharm::CalculateEntropyDensity() {
    double ret = m_partialZ[0] + log(m_Zsum[m_CharmMap[0]]);
    if (m_energydensitiesGCE.size() == 0) CalculateEnergyDensitiesGCE();
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) if (m_CharmMap.count(-m_TPS->Particles()[i].Charm())) ret += (m_Zsum[m_CharmMap[-m_TPS->Particles()[i].Charm()]] / m_Zsum[m_CharmMap[0]]) * ((m_energydensitiesGCE[i] - (m_Parameters.muB * m_TPS->Particles()[i].BaryonCharge() + m_Parameters.muQ * m_TPS->Particles()[i].ElectricCharge() + m_Parameters.muS * m_TPS->Particles()[i].Strangeness()) * m_densitiesGCE[i]) / m_Parameters.T);
    return ret;
  }


  double ThermalModelCanonicalCharm::CalculatePressure() {
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) ret += (m_Zsum[m_CharmMap[-m_TPS->Particles()[i].Charm()]] / m_Zsum[m_CharmMap[0]]) * m_Parameters.T * m_densitiesGCE[i];
    return ret;
  }

  bool ThermalModelCanonicalCharm::IsConservedChargeCanonical(ConservedCharge::Name charge) const {
    return (charge == ConservedCharge::CharmCharge);
  }

} // namespace thermalfist