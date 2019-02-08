/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
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

  /**
   * \brief Class implementing the ideal HRG model
   * with exact conservation of charm (charm-canonical ensemble).
   * 
   * Assumes that the net charm should be equal to zero.
   * 
   * Assumes there are no multi-charmed particles.
   * 
   * Applies Maxwell-Boltzmann statistics to all
   * charmed particles.
   * 
   * Fluctuation observables are not supported.
   * 
   */
  class ThermalModelCanonicalCharm : public ThermalModelBase
  {
  public:
    /**
      * \brief Construct a new ThermalModelCanonicalCharm object.
      *
      * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
      * \param params ThermalModelParameters object with current thermal parameters
      */
    ThermalModelCanonicalCharm(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelCanonicalCharm object
     * 
     */
    virtual ~ThermalModelCanonicalCharm(void);

    /// Calculates the grand-canonical energy densities
    void CalculateEnergyDensitiesGCE();

    // Override functions begin

    /**
     * \copydoc thermalfist::ThermalModelBase::SetParameters(const thermalfist::ThermalModelParameters&)
     * Ensures that charm chemical potential is zero.
     */
    virtual void SetParameters(const ThermalModelParameters& params);

    /**
     * \brief Override the base class method to always
     *        set \f$ \mu_C \f$ to zero
     * 
     * \param muC Value irrelevant
     */
    virtual void SetCharmChemicalPotential(double muC);

    virtual void ChangeTPS(ThermalParticleSystem *TPS);

    //For charm only Boltzmann supported
    /**
     * \copydoc thermalfist::ThermalModelBase::SetStatistics()
     * Ensures that the Maxwell-Boltzmann statistics are applied to all
     * charmed particles.
     */
    virtual void SetStatistics(bool stats);

    void CalculateDensitiesGCE();

    /**
     * \copydoc thermalfist::ThermalModelBase::FixParameters()
     * Ensures that charm chemical potential should not be computed.
     */
    virtual void FixParameters();

    /**
     * \copydoc thermalfist::ThermalModelBase::CalculatePrimordialDensities()
     * Calculates the primordial densities by applying the canonical
     * corrections factors to the grand-canonical densities.
     */
    virtual void CalculatePrimordialDensities();

    /**
     * \brief Dummy function. Fluctuations not yet supported.
     * 
     */
    void CalculateFluctuations();
    
    virtual double CalculatePressure();

    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    // Dummy
    virtual double CalculateBaryonMatterEntropyDensity() { return 0.; }

    virtual double CalculateMesonMatterEntropyDensity() { return 0.; }

    virtual double ParticleScaledVariance(int /*part*/) { return 1.; }

    virtual double ParticleScalarDensity(int /*part*/) { return 0.; }

    virtual bool IsConservedChargeCanonical(ConservedCharge::Name charge) const;

    // Override functions end

  private:
    std::vector<double> m_densitiesGCE;
    std::vector<double> m_energydensitiesGCE;
    std::vector<int>    m_CharmValues;
    std::map<int, int>  m_CharmMap;
    std::vector<double> m_Zsum;
    std::vector<double> m_partialZ;
  };

} // namespace thermalfist

#endif // THERMALMODELCanonicalCharm_H
