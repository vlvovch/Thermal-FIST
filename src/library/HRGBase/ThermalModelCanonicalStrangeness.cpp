/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ThermalModelCanonicalStrangeness.h"

#include "HRGBase/xMath.h"

using namespace std;

namespace thermalfist {

  ThermalModelCanonicalStrangeness::ThermalModelCanonicalStrangeness(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
    ThermalModelBase(TPS_, params)
  {
    m_densitiesGCE.resize(m_TPS->Particles().size());

    m_StrVals.resize(0);
    m_StrVals.push_back(0);
    m_StrVals.push_back(1);
    m_StrVals.push_back(2);
    m_StrVals.push_back(3);
    m_StrVals.push_back(-1);
    m_StrVals.push_back(-2);
    m_StrVals.push_back(-3);

    m_StrMap.clear();
    for (unsigned int i = 0; i < m_StrVals.size(); ++i) m_StrMap[m_StrVals[i]] = i;

    m_Parameters.muS = 0.;

    m_TAG = "ThermalModelCanonicalStrangeness";

    m_Ensemble = SCE;
    m_InteractionModel = Ideal;
  }

  ThermalModelCanonicalStrangeness::~ThermalModelCanonicalStrangeness(void)
  {
  }


  void ThermalModelCanonicalStrangeness::SetParameters(const ThermalModelParameters& params) {
    ThermalModelBase::SetParameters(params);
    m_Parameters.muS = 0.;
  }

  void ThermalModelCanonicalStrangeness::SetStrangenessChemicalPotential(double /*muS*/)
  {
    m_Parameters.muS = 0.0;
  }

  void ThermalModelCanonicalStrangeness::ChangeTPS(ThermalParticleSystem *TPS_) {
    ThermalModelBase::ChangeTPS(TPS_);
    m_densitiesGCE.resize(m_TPS->Particles().size());
  }
  void ThermalModelCanonicalStrangeness::FixParameters() {
    m_ConstrainMuS = false;
    m_ConstrainMuC = false; // No sense doing strangeness CE but charm GCE

    ThermalModelBase::FixParameters();
  }


  void ThermalModelCanonicalStrangeness::SetStatistics(bool stats) {
    m_QuantumStats = stats;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].Strangeness() == 0) m_TPS->Particle(i).UseStatistics(stats);
      else m_TPS->Particle(i).UseStatistics(false);
  }

  void ThermalModelCanonicalStrangeness::CalculateDensitiesGCE() {
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_densitiesGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
    }
    m_GCECalculated = true;
  }

  void ThermalModelCanonicalStrangeness::CalculateEnergyDensitiesGCE() {
    m_energydensitiesGCE.resize(m_densitiesGCE.size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_energydensitiesGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);
    }
  }

  void ThermalModelCanonicalStrangeness::CalculatePressuresGCE()
  {
    m_pressuresGCE.resize(m_densitiesGCE.size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_pressuresGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i]);
    }
  }

  void ThermalModelCanonicalStrangeness::CalculateSums(double Vc)
  {
    if (!m_GCECalculated)
      CalculateDensitiesGCE();

    m_Zsum.resize(m_StrVals.size());

    m_partialS.resize(m_StrVals.size());
    vector<double> xi(3, 0.), yi(3, 0.);

    for (size_t i = 0; i < m_StrVals.size(); ++i) {
      m_partialS[i] = 0.;
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j)
        if (m_StrVals[i] == m_TPS->Particles()[j].Strangeness())
          m_partialS[i] += m_densitiesGCE[j] * Vc;
    }

    for (int i = 0; i < 3; ++i) {
      xi[i] = 2. * sqrt(m_partialS[m_StrMap[i + 1]] * m_partialS[m_StrMap[-(i + 1)]]);
      yi[i] = sqrt(m_partialS[m_StrMap[i + 1]] / m_partialS[m_StrMap[-(i + 1)]]);
    }

    // TODO: Choose iters dynamically based on total strangeness
    int iters = 20;

    for (unsigned int i = 0; i < m_StrVals.size(); ++i) {
      double res = 0.;

      for (int m = -iters; m <= iters; ++m)
        for (int n = -iters; n <= iters; ++n) {
          double tmp = xMath::BesselI(abs(3 * m + 2 * n - m_StrVals[i]), xi[0]) *
            xMath::BesselI(abs(n), xi[1]) *
            xMath::BesselI(abs(m), xi[2]) *
            pow(yi[0], m_StrVals[i] - 3 * m - 2 * n) *
            pow(yi[1], n) *
            pow(yi[2], m);
          if (tmp != tmp) continue;
          res += tmp;
        }
      m_Zsum[i] = res;
    }
  }

  void ThermalModelCanonicalStrangeness::CalculatePrimordialDensities() {
    if (UsePartialChemicalEquilibrium()) {
      printf("**ERROR** ThermalModelCanonicalStrangeness::CalculatePrimordialDensities(): PCE not supported!\n");
      exit(1);
    }

    m_FluctuationsCalculated = false;

    m_energydensitiesGCE.resize(0);
    m_pressuresGCE.resize(0);

    CalculateDensitiesGCE();

    CalculateSums(m_Parameters.SVc);

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
        m_densities[i] = (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_densitiesGCE[i];
    }

    m_Calculated = true;
    ValidateCalculation();
  }


  void ThermalModelCanonicalStrangeness::CalculateFluctuations() {
    m_wprim.resize(m_densities.size());
    m_wtot.resize(m_densities.size());
    for (size_t i = 0; i < m_wprim.size(); ++i) {
      m_wprim[i] = 1.;
      m_wtot[i] = 1.;
      m_skewprim[i] = 1.;
      m_skewtot[i] = 1.;
    }
  }


  double ThermalModelCanonicalStrangeness::CalculateEnergyDensity() {
    if (!m_Calculated)
      CalculateDensities();
    if (m_energydensitiesGCE.size() == 0)
      CalculateEnergyDensitiesGCE();
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
        ret += (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_energydensitiesGCE[i];
    return ret;
  }


  double ThermalModelCanonicalStrangeness::CalculateEntropyDensity() {
    double ret = log(m_Zsum[m_StrMap[0]]) / m_Parameters.SVc;
    if (m_energydensitiesGCE.size() == 0)
      CalculateEnergyDensitiesGCE();
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].Strangeness() != 0) {
        if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
          ret += (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * ((m_energydensitiesGCE[i] - m_Chem[i] * m_densitiesGCE[i]) / m_Parameters.T);
      }
      else {
        ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i]);
      }
    return ret;
  }


  double ThermalModelCanonicalStrangeness::CalculatePressure() {
    double ret = 0.;
    if (m_pressuresGCE.size() == 0)
      CalculatePressuresGCE();
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
        ret += (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_pressuresGCE[i];
    return ret;
  }

  bool ThermalModelCanonicalStrangeness::IsConservedChargeCanonical(ConservedCharge::Name charge) const {
    return (charge == ConservedCharge::StrangenessCharge);
  }

} // namespace thermalfist