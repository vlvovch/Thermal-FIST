/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELIDEAL_H
#define THERMALMODELIDEAL_H

#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  /**
   * \brief Class implementing the Ideal HRG model.
   * 
   */
  class ThermalModelIdeal : public ThermalModelBase
  {
  public:
    /**
    * \brief Construct a new ThermalModelIdeal object.
    *
    * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
    * \param params ThermalModelParameters object with current thermal parameters
    */
    ThermalModelIdeal(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelIdeal object
     * 
     */
    virtual ~ThermalModelIdeal(void);

    // Override functions begin

    virtual void CalculatePrimordialDensities();

    virtual void CalculateTwoParticleCorrelations();

    virtual void CalculateFluctuations();

    virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);

    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    virtual double CalculateBaryonMatterEntropyDensity();

    virtual double CalculateMesonMatterEntropyDensity();

    virtual double CalculatePressure();

    virtual double ParticleScaledVariance(int part);

    virtual double ParticleSkewness(int part);

    virtual double ParticleKurtosis(int part);

    virtual double ParticleScalarDensity(int part);

    // Override functions end
  };

} // namespace thermalfist

#endif
