/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELCANONICAL_H
#define THERMALMODELCANONICAL_H

#include <map>


#include "HRGBase/ThermalModelIdeal.h"

namespace thermalfist {

  /**
   * \brief Struct containing a set of quantum numbers:
   *        Baryon number, electric charge, strangeness, and charm 
   * 
   */
  struct QuantumNumbers
  {
    int B; ///< Baryon number
    int Q; ///< Electric charge
    int S; ///< Strangeness
    int C; ///< Charm 

    /**
     * \brief Construct a new QuantumNumbers object
     * 
     * \param iB Baryon number
     * \param iQ Electric charge
     * \param iS Strangeness
     * \param iC Charm
     */
    QuantumNumbers(int iB = 0, int iQ = 0, int iS = 0, int iC = 0) :
      B(iB), Q(iQ), S(iS), C(iC) { }
    const bool operator < (const QuantumNumbers &r) const {
      if (B != r.B)
        return (B < r.B);
      else if (Q != r.Q)
        return (Q < r.Q);
      else if (S != r.S)
        return (S < r.S);
      else
        return (C < r.C);
    }
  };

  /**
   * \brief Class implementing the ideal HRG model
   *        in the canonical ensemble.
   * 
   * Calculation of particle densities proceeds as described in
   * F. Becattini, U. Heinz, Z. Phys. C76, 269 (1997), [https://arxiv.org/pdf/hep-ph/9702274.pdf](https://arxiv.org/pdf/hep-ph/9702274.pdf)
   * 
   * Particle number fluctuations are determined through the derivatives with respect to the
   * fictitous fugacities, as is usually done in HRG,
   * see e.g. [https://arxiv.org/pdf/nucl-th/0404056.pdf](https://arxiv.org/pdf/nucl-th/0404056.pdf) 
   * 
   */
  class ThermalModelCanonical :
    public ThermalModelBase
  {
  public:
    /**
      * \brief Construct a new ThermalModelCanonical object.
      *
      * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
      * \param params ThermalModelParameters object with current thermal parameters
      */
    ThermalModelCanonical(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelCanonical object
     * 
     */
    virtual ~ThermalModelCanonical(void);

    /**
     * \brief Calculates the range of quantum numbers values
     *        for which it is necessary to compute the 
     *        canonical partition functions.
     * 
     * \param computeFluctuations Whether it will be necessary to compute
     *                            fluctuations. If that is the case the range is doubled
     *                            for every quantum number. 
     */
    virtual void CalculateQuantumNumbersRange(bool computeFluctuations = false);
    
    /**
     * \brief Calculates all necessary canonical partition functions.
     * 
     * This corresponds to Eq. (8) in [https://arxiv.org/pdf/hep-ph/9702274.pdf](https://arxiv.org/pdf/hep-ph/9702274.pdf)
     * 
     * Integrals are performed numerically using quadratures.
     * 
     * If multi-baryon states (light nuclei) are not included in the list,
     * and quantum statistics for baryons is neglected,
     * the integral over the baryon fugacity is performed analytically as
     * described in (https://arxiv.org/pdf/nucl-th/0112021.pdf)[https://arxiv.org/pdf/nucl-th/0112021.pdf].
     * 
     * \param Vc 
     */
    virtual void CalculatePartitionFunctions(double Vc = -1.);

    /**
     * \brief Determines whether the specified ThermalParticle
     *        is treat canonically or grand-canonically in the present
     *        setup. 
     * 
     * This depends on whether the particle carries any of the
     * exactly conserved charges.
     * 
     * \param part The particle species.
     * \return true  Canonically.
     * \return false Grand-canonically.
     */
    virtual bool IsParticleCanonical(const ThermalParticle &part);

    /**
     * \brief Specifies whether the baryon number is treated canonically.
     * 
     * By default the baryon number is treated canonically.
     * 
     * \param conserve true -- canonically, false -- grand-canonically
     */
    virtual void ConserveBaryonCharge(bool conserve = true) { m_BCE = static_cast<int>(conserve); }

    /**
     * \brief Specifies whether the electric charge is treated canonically.
     * 
     * By default the electric charge is treated canonically.
     * 
     * \param conserve true -- canonically, false -- grand-canonically
     */
    virtual void ConserveElectricCharge(bool conserve = true) { m_QCE = static_cast<int>(conserve); }

    /**
     * \brief Specifies whether the strangeness charge is treated canonically.
     * 
     * By default the strangeness charge is treated canonically.
     * 
     * \param conserve true -- canonically, false -- grand-canonically
     */
    virtual void ConserveStrangeness(bool conserve = true) { m_SCE = static_cast<int>(conserve); }

    /**
     * \brief Specifies whether the charm charge is treated canonically.
     * 
     * By default the charm charge is treated canonically.
     * 
     * \param conserve true -- canonically, false -- grand-canonically
     */
    virtual void ConserveCharm(bool conserve = true) { m_CCE = static_cast<int>(conserve); }

    virtual bool IsConservedChargeCanonical(ConservedCharge::Name charge) const;

    /**
     * \brief Density of particle species i in the grand-canonical ensemble.
     * 
     * \param i       0-based index of particle species
     * \return double The grand-canonical density
     */
    virtual double GetGCEDensity(int i) const;

    /**
     * \brief The multiplier of the number of iterations in the numerical integration
     *
     * \return The multiplier
     */
    int IntegrationIterationsMultiplier() const { return m_IntegrationIterationsMultiplier; }

    /**
     * \brief Assigns the multiplier of the number of iterations in the numerical integration
     *
     * The minimum value of multiplier is 1. Increase to improve the numerical accuracy of the canonical ensemble calculations.
     *
     * \param The multiplier
     */
    void SetIntegrationIterationsMultiplier(int multiplier) { (multiplier > 0 ? m_IntegrationIterationsMultiplier = multiplier : m_IntegrationIterationsMultiplier = 1); }
    

    // Override functions begin

    void ChangeTPS(ThermalParticleSystem *TPS);

    virtual void SetStatistics(bool stats);

    virtual void FixParameters();

    virtual void FixParametersNoReset();

    virtual void CalculatePrimordialDensities();

    virtual void ValidateCalculation();

    virtual double ParticleScaledVariance(int part);

    virtual void CalculateTwoParticleCorrelations();

    /**
     * \copydoc thermalfist::ThermalModelBase::CalculateFluctuations()
     * Restricted to 2nd moments.
     * TODO: Higher moments
     * 
     */
    virtual void CalculateFluctuations();

    virtual double CalculateEnergyDensity();

    virtual double CalculatePressure();

    virtual double CalculateEntropyDensity();

    virtual double ParticleScalarDensity(int /*part*/) { return 0.; }

    // Override functions end

  private:
    //@{
    /**
     * Functions for a different evaluation method,
     * which uses finite conserved charges chemical potentials
     * instead of zero chemical potentials.
     * Currently not used. 
     * 
     */
    void PrepareModelGCE();  /**< Creates the ThermalModelIdeal copy */

    void CleanModelGCE();    /**< Cleares the ThermalModelIdeal copy */
    //@}

  protected:

    /**
     * \brief A set of QuantumNumbers combinations
     *        for which it is necessary to compute the
     *        canonical partition function.
     * 
     */
    std::vector<QuantumNumbers> m_QNvec;

    /**
     * \brief Maps QuantumNumbers combinations
     *        to a 1-dimensional index.
     * 
     */
    std::map<QuantumNumbers, int> m_QNMap;
    
    /**
     * \brief A vector of chemical factors.
     * 
     * Chemical factors define the canonical corrections
     * to the grand canonical thermodynamic functions.
     * 
     */
    std::vector<double> m_Corr;

    /**
     * \brief The computed canonical partition function.
     * 
     * The partition functions are computed up to
     * a factor which does not depend of the quantum
     * numbers since only ratios of the partition functions
     * -- the chemical factors -- are of relevance.
     * 
     */
    std::vector<double> m_PartialZ;

    /**
     * \brief A multiplier to increase the number of iterations during the numerical integration used to calculate the partition functions.
     *
     * Set with SetIntegrationIterationsMultiplier()
     *
     */
    int m_IntegrationIterationsMultiplier;

    int m_BMAX, m_QMAX, m_SMAX, m_CMAX;
    int m_BMAX_list, m_QMAX_list, m_SMAX_list, m_CMAX_list;

    double m_MultExp;
    double m_MultExpBanalyt;

    ThermalModelIdeal *m_modelgce;

    int m_BCE, m_QCE, m_SCE, m_CCE;

    bool m_Banalyt;
  };

} // namespace thermalfist

#endif
