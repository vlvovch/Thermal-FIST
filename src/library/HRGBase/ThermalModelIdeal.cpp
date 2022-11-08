/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ThermalModelIdeal.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <cmath>


using namespace std;

namespace thermalfist {

  ThermalModelIdeal::ThermalModelIdeal(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
    ThermalModelBase(TPS_, params)
  {
    m_TAG = "ThermalModelIdeal";

    m_Ensemble = GCE;
    m_InteractionModel = Ideal;
  }

  ThermalModelIdeal::~ThermalModelIdeal(void)
  {
  }

  void ThermalModelIdeal::CalculatePrimordialDensities() {
    m_FluctuationsCalculated = false;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_densities[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
    }

    m_Calculated = true;
    ValidateCalculation();
  }

  void ThermalModelIdeal::CalculateTwoParticleCorrelations() {
    int NN = m_densities.size();
    vector<double> tN(NN);
    for (int i = 0; i < NN; ++i) 
      tN[i] = m_densities[i];

    vector<double> chi2s(NN);
    for (int i = 0; i < NN; ++i) 
      chi2s[i] = m_TPS->Particles()[i].chiDimensionfull(2, m_Parameters, m_UseWidth, m_Chem[i]);

    m_PrimCorrel.resize(NN);
    for (int i = 0; i < NN; ++i) m_PrimCorrel[i].resize(NN);
    m_TotalCorrel = m_PrimCorrel;

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        m_PrimCorrel[i][j] = 0.;
        if (i == j) m_PrimCorrel[i][j] += chi2s[i] * pow(xMath::GeVtoifm(), 3);
      }

    for (int i = 0; i < NN; ++i) {
      m_wprim[i] = m_PrimCorrel[i][i];
      if (m_densities[i] > 0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
      else m_wprim[i] = 1.;
    }

  }


  void ThermalModelIdeal::CalculateFluctuations() {
    CalculateTwoParticleCorrelations();
    CalculateSusceptibilityMatrix();
    CalculateTwoParticleFluctuationsDecays();
    CalculateProxySusceptibilityMatrix();
    CalculateParticleChargeCorrelationMatrix();

    for (size_t i = 0; i < m_wprim.size(); ++i) {
      m_wprim[i] = ParticleScaledVariance(i);
      m_skewprim[i] = ParticleSkewness(i);
      m_kurtprim[i] = ParticleKurtosis(i);
    }
    for (size_t i = 0; i < m_wtot.size(); ++i) {
      double tmp1 = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
      tmp2 = m_densities[i] * m_wprim[i];
      tmp3 = m_densities[i] * m_wprim[i] * m_skewprim[i];
      tmp4 = m_densities[i] * m_wprim[i] * m_kurtprim[i];
      const ThermalParticleSystem::DecayContributionsToParticle& decayContributions = m_TPS->DecayContributionsByFeeddown()[Feeddown::StabilityFlag][i];
      for (size_t r = 0; r < decayContributions.size(); ++r) {
        const ThermalParticleSystem::SingleDecayContribution& decayContrib = decayContributions[r];
        const ThermalParticleSystem::SingleDecayCumulantsContribution& decayCumulantsSingle = m_TPS->DecayCumulants()[i][r];
        tmp2 += m_densities[decayContrib.second] *
          (m_wprim[decayContrib.second] * decayContrib.first * decayContrib.first
            + decayCumulantsSingle.first[1]);

        int rr = decayContrib.second;
        double ni = decayContrib.first;
        tmp3 += m_densities[rr] * m_wprim[rr] * (m_skewprim[rr] * ni * ni * ni + 3. * ni * decayCumulantsSingle.first[1]);
        tmp3 += m_densities[rr] * decayCumulantsSingle.first[2];

        tmp4 += m_densities[rr] * m_wprim[rr] * (m_kurtprim[rr] * ni * ni * ni * ni
          + 6. * m_skewprim[rr] * ni * ni * decayCumulantsSingle.first[1]
          + 3. * decayCumulantsSingle.first[1] * decayCumulantsSingle.first[1]
          + 4. * ni * decayCumulantsSingle.first[2]);

        tmp4 += m_densities[rr] * decayCumulantsSingle.first[3];
      }


      tmp1 = m_densitiestotal[i];

      m_wtot[i] = tmp2 / tmp1;
      m_skewtot[i] = tmp3 / tmp2;
      m_kurttot[i] = tmp4 / tmp2;
    }

    m_FluctuationsCalculated = true;
  }

  std::vector<double> ThermalModelIdeal::CalculateChargeFluctuations(const std::vector<double>& chgs, int order)
  {
    vector<double> ret(order + 1, 0.);

    // chi1
    for (size_t i = 0; i < m_densities.size(); ++i)
      ret[0] += chgs[i] * m_densities[i];

    ret[0] /= pow(m_Parameters.T * xMath::GeVtoifm(), 3);

    if (order < 2) return ret;

    for (size_t i = 0; i < m_densities.size(); ++i)
      ret[1] += chgs[i] * chgs[i] * m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_Chem[i]);

    if (order < 3) return ret;

    for (size_t i = 0; i < m_densities.size(); ++i)
      ret[2] += chgs[i] * chgs[i] * chgs[i] * m_TPS->Particles()[i].chi(3, m_Parameters, m_UseWidth, m_Chem[i]);

    if (order < 4) return ret;

    for (size_t i = 0; i < m_densities.size(); ++i)
      ret[3] += chgs[i] * chgs[i] * chgs[i] * chgs[i] * m_TPS->Particles()[i].chi(4, m_Parameters, m_UseWidth, m_Chem[i]);

    return ret;
  }

  double ThermalModelIdeal::CalculateEnergyDensity() {
    double ret = 0.;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);

    return ret;
  }

  double ThermalModelIdeal::CalculateEntropyDensity() {
    double ret = 0.;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i]);

    return ret;
  }

  double ThermalModelIdeal::CalculateBaryonMatterEntropyDensity() {
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].BaryonCharge() != 0)
        ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i]);
    return ret;
  }

  double ThermalModelIdeal::CalculateMesonMatterEntropyDensity() {
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].BaryonCharge() == 0)
        ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i]);
    return ret;
  }

  double ThermalModelIdeal::CalculatePressure() {
    double ret = 0.;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i]);

    return ret;
  }

  double ThermalModelIdeal::ParticleScaledVariance(int part) {
    return m_TPS->Particles()[part].ScaledVariance(m_Parameters, m_UseWidth, m_Chem[part]);
  }

  double ThermalModelIdeal::ParticleSkewness(int part) {
    return m_TPS->Particles()[part].Skewness(m_Parameters, m_UseWidth, m_Chem[part]);
  }

  double ThermalModelIdeal::ParticleKurtosis(int part) {
    return m_TPS->Particles()[part].Kurtosis(m_Parameters, m_UseWidth, m_Chem[part]);
  }

  double ThermalModelIdeal::ParticleScalarDensity(int part) {
    return m_TPS->Particles()[part].Density(m_Parameters, IdealGasFunctions::ScalarDensity, m_UseWidth, m_Chem[part]);
  }

} // namespace thermalfist