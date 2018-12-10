/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELCanonicalCharm_H
#define THERMALMODELCanonicalCharm_H

#include <map>
#include <iostream>

#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/SplineFunction.h"

namespace thermalfist {

  class ThermalModelCanonicalCharm : public ThermalModelBase
  {
  public:
    ThermalModelCanonicalCharm(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters());

    virtual ~ThermalModelCanonicalCharm(void);

    virtual void SetParameters(const ThermalModelParameters& params);
    virtual void SetCharmChemicalPotential(double muC);

    virtual void ChangeTPS(ThermalParticleSystem *TPS_);

    //For charm only Boltzmann supported
    virtual void SetStatistics(bool stats);

    void CalculateDensitiesGCE();

    void CalculateEnergyDensitiesGCE();

    void CalculatefPhi(int iters = 300);

    virtual void FixParameters();

    virtual void CalculateDensities();

    // TODO properly
    void CalculateFluctuations();

    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    // Dummy
    virtual double CalculateBaryonMatterEntropyDensity() { return 0.; }
    virtual double CalculateMesonMatterEntropyDensity() { return 0.; }

    virtual double CalculatePressure();

    virtual double CalculateShearViscosity() { return 0.; }

    virtual double CalculateParticleScaledVariance(int part) { return 1.; }

    virtual double CalculateHadronScaledVariance() { return 1.; }

    virtual double ParticleScalarDensity(int part) { return 0.; }

  private:
    std::vector<double> m_densitiesGCE;
    std::vector<double> m_energydensitiesGCE;
    SplineFunction        m_phiRe;
    SplineFunction        m_phiIm;
    std::vector<int>      m_StrVals;
    std::map<int, int>    m_StrMap;
    std::vector<double> m_Zsum;
    std::vector<double> m_partialS;
  };

} // namespace thermalfist

#endif // THERMALMODELCanonicalCharm_H
