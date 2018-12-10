/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELCANONICALSTRANGENESS_H
#define THERMALMODELCANONICALSTRANGENESS_H

#include <map>
#include <iostream>

#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/SplineFunction.h"

namespace thermalfist {

  class ThermalModelCanonicalStrangeness : public ThermalModelBase
  {
  public:
    ThermalModelCanonicalStrangeness(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters());

    virtual ~ThermalModelCanonicalStrangeness(void);

    virtual void SetParameters(const ThermalModelParameters& params);

    virtual void SetStrangenessChemicalPotential(double muS);

    virtual void ChangeTPS(ThermalParticleSystem *TPS_);

    virtual void FixParameters();

    //For strangeness only Boltzmann supported
    virtual void SetStatistics(bool stats);

    // FillChemicalPotential override! or maybe not...
    virtual void CalculateDensitiesGCE();

    virtual void CalculateEnergyDensitiesGCE();

    virtual void CalculatePressuresGCE();

    void CalculatefPhi(int iters = 300);

    virtual void CalculateSums(double Vc);
    virtual void CalculateDensities();

    // TODO properly all fluctuations
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

    const std::vector<double>& DensitiesGCE() const { return m_densitiesGCE; }

  protected:
    std::vector<double> m_densitiesGCE;
    std::vector<double> m_energydensitiesGCE;
    std::vector<double> m_pressuresGCE;
    SplineFunction        m_phiRe;
    SplineFunction        m_phiIm;
    std::vector<int>      m_StrVals;
    std::map<int, int>    m_StrMap;
    std::vector<double> m_Zsum;
    std::vector<double> m_partialS;
  };

} // namespace thermalfist

#endif // THERMALMODELCANONICALSTRANGENESS_H
