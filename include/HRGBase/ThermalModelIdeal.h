/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELIDEAL_H
#define THERMALMODELIDEAL_H

#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  class ThermalModelIdeal : public ThermalModelBase
  {
  public:
    ThermalModelIdeal(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters());

    virtual ~ThermalModelIdeal(void);

    virtual void CalculateDensities();

    virtual void CalculateTwoParticleCorrelations();

    virtual void CalculateFluctuations();

    virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);

    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    virtual double CalculateBaryonMatterEntropyDensity();

    virtual double CalculateMesonMatterEntropyDensity();

    virtual double CalculatePressure();

    virtual double CalculateShearViscosity();

    virtual double CalculateHadronScaledVariance();

    virtual double CalculateParticleScaledVariance(int part);

    virtual double CalculateParticleSkewness(int part);

    virtual double CalculateParticleKurtosis(int part);

    virtual double ParticleScalarDensity(int part);
  };

} // namespace thermalfist

#endif
