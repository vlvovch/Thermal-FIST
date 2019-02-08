/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
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

  /**
   * \brief Class implementing the ideal HRG model
   * with exact conservation of strangeness (strangeness-canonical ensemble).
   * 
   * Assumes that the net strangeness should be equal to zero.
   * 
   * Assumes there are no particles with strangeness larger than 3. 
   * 
   * Applies Maxwell-Boltzmann statistics to all
   * strange particles.
   * 
   * Fluctuation observables are not supported.
   * 
   */
  class ThermalModelCanonicalStrangeness : public ThermalModelBase
  {
  public:
    /**
      * \brief Construct a new ThermalModelCanonicalStrangeness object.
      *
      * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
      * \param params ThermalModelParameters object with current thermal parameters
      */
    ThermalModelCanonicalStrangeness(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelCanonicalStrangeness object
     * 
     */
    virtual ~ThermalModelCanonicalStrangeness(void);

    /// Calculates the grand-canonical energy densities
    virtual void CalculateEnergyDensitiesGCE();

    /// Calculates the grand-canonical pressures
    virtual void CalculatePressuresGCE();

    /**
     * \brief Calculates the strangeness-canonical partition functions.
     * 
     * Uses the series over the Bessel functions 
     * as described in [https://arxiv.org/pdf/hep-ph/0106066.pdf](https://arxiv.org/pdf/hep-ph/0106066.pdf)
     * 
     * \param Vc The strangeness correlation volume (fm\f$^3\f$)
     */
    virtual void CalculateSums(double Vc);

    /// A vector of the grand-canonical particle number densities
    const std::vector<double>& DensitiesGCE() const { return m_densitiesGCE; }

    // Override functions begin

    /**
     * \copydoc thermalfist::ThermalModelBase::SetParameters(const thermalfist::ThermalModelParameters&)
     * Ensures that strangeness chemical potential is zero.
     */
    virtual void SetParameters(const ThermalModelParameters& params);

    /**
     * \brief Override the base class method to always
     *        set \f$ \mu_S \f$ to zero
     * 
     * \param muS Value irrelevant
     */
    virtual void SetStrangenessChemicalPotential(double muS);

    virtual void ChangeTPS(ThermalParticleSystem *TPS);

    /**
     * \copydoc thermalfist::ThermalModelBase::FixParameters()
     * Ensures that charm chemical potential should not be computed.
     */
    virtual void FixParameters();

    /**
     * \copydoc thermalfist::ThermalModelBase::SetStatistics()
     * Ensures that the Maxwell-Boltzmann statistics are applied to all
     * strange particles.
     */
    virtual void SetStatistics(bool stats);

    virtual void CalculateDensitiesGCE();

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

  protected:
    std::vector<double> m_densitiesGCE;
    std::vector<double> m_energydensitiesGCE;
    std::vector<double> m_pressuresGCE;
    std::vector<int>    m_StrVals;
    std::map<int, int>  m_StrMap;
    std::vector<double> m_Zsum;
    std::vector<double> m_partialS;
  };

} // namespace thermalfist

#endif // THERMALMODELCANONICALSTRANGENESS_H
