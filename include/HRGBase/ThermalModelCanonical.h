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
   * \brief Method used for computing the canonical partition functions.
   *
   * The canonical ensemble constrains all conserved charges (B, Q, S, C) to
   * fixed integer values.  Computing canonical densities and fluctuations is
   * done by first computing the canonical partition functions for each
   * set of quantum numbers, and then constructing the chemical factors.
   *
   * Four algorithms are available, trading accuracy for speed:
   *
   * - **GaussLegendre** (default): Numerically integrates the Fourier
   *   representation of the canonical partition functions using composite
   *   Gauss-Legendre quadrature.  Most accurate for small systems but
   *   O(N_QN * N_int) per species.
   *
   * - **SaddlePoint**: Evaluates the canonical partition function for each
   *   quantum-number sector at a single saddle point.  Fast, but computes
   *   a separate saddle point for each QN combination.
   *
   * - **SaddlePointNLO**: A single saddle-point solve for the dominant QN
   *   sector, followed by a systematic 1/V expansion of means, covariances,
   *   and thermodynamic quantities up to next-to-leading order (NLO).
   *   Preserves exact conservation laws.  Ideal-gas models only.
   *
   * - **SaddlePointLO**: Same saddle-point solve as NLO but keeps only
   *   the leading-order terms.  Densities and thermodynamic quantities
   *   coincide with the GCE evaluated at the saddle-point chemical
   *   potentials μ*.  Fluctuations (covariances) are modified by canonical
   *   conservation laws.  Supports non-ideal interaction models via
   *   SetModelGCE().
   *
   * The SaddlePointLO covariance formula is equivalent to the subensemble
   * acceptance method (SAM) at acceptance fraction α = 1; see
   * V. Vovchenko, R.V. Poberezhnyuk, V. Koch,
   * JHEP 10 (2020) 089 [arXiv:2007.03850], Eq. (2.85).
   */
  enum CanonicalMethod {
    GaussLegendre   = 0,  ///< Composite 10-point Gauss-Legendre on subintervals (default)
    SaddlePoint     = 1,  ///< Saddle-point approximation (per-QN solve)
    SaddlePointNLO  = 2,  ///< Saddle-point with systematic NLO expansion (preserves exact conservation laws)
    SaddlePointLO   = 3   ///< Saddle-point LO only (GCE densities/thermo, canonical LO fluctuations)
  };

  /**
   * \brief Class implementing the HRG model in the canonical ensemble.
   *
   * Computes particle densities and fluctuations with exact conservation
   * of an arbitrary subset of {B, Q, S, C}.  Several computational methods
   * are available (see CanonicalMethod):
   *
   * - **GaussLegendre / SaddlePoint**: numerical integration or per-QN
   *   saddle-point evaluation; ideal gas only.
   * - **SaddlePointNLO**: systematic 1/V expansion including O(1/V_c)
   *   corrections to means, covariances, and thermodynamics; ideal gas only.
   * - **SaddlePointLO**: leading-order saddle-point; supports non-ideal
   *   equations of state via SetModelGCE().
   *
   * ## Non-ideal interactions (SaddlePointLO only)
   *
   * Provide a pre-configured GCE model via SetModelGCE():
   * \code
   * auto* rg = new ThermalModelRealGas(&TPS, params);
   * rg->SetExcludedVolumeModel(...);
   * model.SetModelGCE(rg);    // borrow; caller retains ownership
   * \endcode
   *
   * ## LO covariance formula
   *
   * The SaddlePointLO covariance is computed from the general formula
   * \f[
   *   \mathrm{cov}(N_l, N_m) = \chi^{\rm GCE}_{lm}
   *     - \sum_{a,b} v_{l,a}\,[\Sigma_\chi^{-1}]_{ab}\,v_{m,b}
   * \f]
   * where \f$\chi^{\rm GCE}_{jk}\f$ is the GCE species-level susceptibility
   * matrix (diagonal for ideal gas, full for interacting models),
   * \f$v_{l,a} = \sum_j q_{a,j}\,\chi^{\rm GCE}_{jl}\f$,
   * and \f$\Sigma_{ab} = \sum_{j,k} q_{a,j}\,q_{b,k}\,\chi^{\rm GCE}_{jk}\f$.
   *
   * This is equivalent to the subensemble acceptance method (SAM)
   * at full acceptance (\f$\alpha = 1\f$), specifically Eq. (2.86) of
   * V. Vovchenko, R.V. Poberezhnyuk, V. Koch,
   * JHEP 10 (2020) 089 [arXiv:2007.03850]:
   * \f$\chi_{pp}^{\rm ce} = \det\tilde\chi / \det\chi\f$,
   * via the Schur complement identity.
   *
   * ## References
   *
   * - F. Becattini, U. Heinz, Z. Phys. C76, 269 (1997)
   *   [hep-ph/9702274](https://arxiv.org/pdf/hep-ph/9702274.pdf)
   * - Fluctuations via fictitious fugacities:
   *   [nucl-th/0404056](https://arxiv.org/pdf/nucl-th/0404056.pdf)
   * - SAM formalism:
   *   V. Vovchenko, R.V. Poberezhnyuk, V. Koch,
   *   JHEP 10 (2020) 089 [arXiv:2007.03850](https://arxiv.org/abs/2007.03850)
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
    virtual void ConserveBaryonCharge(bool conserve = true) { m_BCE = static_cast<int>(conserve); m_PartialZCalculated = false;  }

    /**
     * \brief Specifies whether the electric charge is treated canonically.
     * 
     * By default the electric charge is treated canonically.
     * 
     * \param conserve true -- canonically, false -- grand-canonically
     */
    virtual void ConserveElectricCharge(bool conserve = true) { m_QCE = static_cast<int>(conserve); m_PartialZCalculated = false; }

    /**
     * \brief Specifies whether the strangeness charge is treated canonically.
     * 
     * By default the strangeness charge is treated canonically.
     * 
     * \param conserve true -- canonically, false -- grand-canonically
     */
    virtual void ConserveStrangeness(bool conserve = true) { m_SCE = static_cast<int>(conserve); m_PartialZCalculated = false; }

    /**
     * \brief Specifies whether the charm charge is treated canonically.
     * 
     * By default the charm charge is treated canonically.
     * 
     * \param conserve true -- canonically, false -- grand-canonically
     */
    virtual void ConserveCharm(bool conserve = true) { m_CCE = static_cast<int>(conserve); m_PartialZCalculated = false; }

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

    /**
     * \brief Returns the method used for computing the canonical partition functions.
     */
    CanonicalMethod Method() const { return m_Method; }

    /**
     * \brief Sets the method for computing the canonical partition functions.
     *
     * \param method GaussLegendre (default) or SaddlePoint
     */
    void SetMethod(CanonicalMethod method) { m_Method = method; m_PartialZCalculated = false; }

    /**
     * \brief Sets an external GCE model to use in the saddle-point approach.
     *
     * This is the only way to use non-ideal interactions with the canonical
     * ensemble (SaddlePointLO method only).  The provided model is borrowed
     * (not owned) by the canonical model.  The caller retains ownership and
     * must ensure the model remains valid during canonical calculations.
     *
     * The model's thermal parameters and chemical potentials will be
     * overwritten during the saddle-point solve.
     * The interaction parameters (virial, attraction, etc.) are preserved.
     *
     * Pass nullptr to revert to the default ideal gas.
     *
     * \param model  Pointer to a pre-configured GCE model, or nullptr to clear
     */
    void SetModelGCE(ThermalModelBase* model);

    /*
     * \brief Reset all flags which correspond to a calculation status
     */
    virtual void ResetCalculatedFlags() {
      ThermalModelBase::ResetCalculatedFlags();
      ResetPartialZCalculated();
    }
    /*
     * \brief Reset the flags indicating whether the partial partition functions are calculated
    */
    void ResetPartialZCalculated() { m_PartialZCalculated = false; m_PartialZCalculatedFlucts = false; }

    /*
     * \brief Whether the partial partition functions are calculated
     */
    bool IsPartialZCalculated() const { return m_PartialZCalculated; }
    
    /*
     * \brief Whether the partial partition functions are calculated with fluctuations
     */
    bool IsPartialZCalculatedFlucts() const { return m_PartialZCalculatedFlucts; }

    // Override functions begin

    void ChangeTPS(ThermalParticleSystem *TPS);

    virtual void SetStatistics(bool stats);

    virtual void FixParameters();

    virtual void FixParametersNoReset();

    /// Solve non-canonical chemical potentials using a temporary GCE model,
    /// avoiding the expensive canonical partition function computation.
    void SolveChemicalPotentialsGCE(bool resetInitialValues);

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

    virtual double CalculateEnergyDensityDerivativeT() { throw std::runtime_error("CalculateEnergyDensityDerivativeT not implemented"); return 0.; } // Not implemented

    virtual double CalculateEntropyDensityDerivativeT() { throw std::runtime_error("CalculateEntropyDensityDerivativeT not implemented"); return 0.; } // Not implemented

    /**
     * \brief Returns the scalar density of a particle species.
     * 
     * This method is not implemented for the canonical ensemble model,
     * hence it returns a constant value of 0.0.
     * 
     * \param part The particle species index.
     * \return double The scalar density, which is 0.0 in this implementation.
     */
    virtual double ParticleScalarDensity(int /*part*/) { return 0.; }

    // Override functions end

  private:
    //@{
    /**
     * \brief Internal GCE model management for the saddle-point approach.
     *
     * The saddle-point solve requires a GCE model to evaluate densities
     * and susceptibilities at the saddle-point chemical potentials μ*.
     *
     * **PrepareModelGCE()** selects the GCE model to use:
     *   - If an external model was set via SetModelGCE(), use it.
     *   - Otherwise, use the persistent default ThermalModelIdeal.
     *   Then solve for chemical potentials matching the target quantum numbers.
     *
     * **CleanModelGCE()** resets the active model pointer (no deletion).
     */
    virtual void PrepareModelGCE();
    void CleanModelGCE();

    /// Enforce that non-ideal interaction models only run with SaddlePointLO.
    /// Any other method falls back to SaddlePointLO (with a warning), since the
    /// remaining methods assume an ideal gas.
    void EnforceSaddlePointLOForNonIdeal();
    //@}

    //@{
    /**
     * Helper methods for the saddle-point approximation (SaddlePoint quadrature method).
     */
    void ComputeSaddlePointPartitionFunctions();  ///< Orchestrates the full saddle-point computation
    void SolveSaddlePointEquations();             ///< Finds saddle-point chemical potentials mu* /T
    void SetupSaddlePointChargeIndices();          ///< Determines which charges (B,Q,S,C) are canonical

    /// Build Sigma from m_modelgce->PrimCorrel() (full interacting chi^GCE).
    /// Returns log(det(Sigma)).  Requires m_modelgce with two-particle correlations computed.
    double ComputeSigmaLogDetFromModel();

    /**
     * \brief Compute analytic means and (optionally) covariances from the
     *        saddle-point cumulant generating function expansion.
     *
     * This implements both SaddlePointNLO and SaddlePointLO modes:
     *
     * **Algorithm outline:**
     * 1. Solve the saddle-point equations for mu* /T.
     * 2. Compute GCE two-particle correlations at the saddle point
     *    (via m_modelgce->CalculateTwoParticleCorrelations()).
     * 3. Build per-species cluster moments W_j^{(k)} at μ*.
     * 4. Build the d×d charge susceptibility matrix Σ and its inverse.
     * 5. Compute primordial densities (LO, or LO + NLO correction).
     * 6. If computeCovariance is true:
     *    - **Ideal gas** (fast path): cov(l,m) = δ_{lm} W2[l] − W2[l] W2[m] G_{lm},
     *      using precomputed Σ^{−1}. No N×N allocation needed.
     *    - **Non-ideal** (general path): use the full N×N PrimCorrel from
     *      m_modelgce, then cov(l,m) = χ^{GCE}_{lm} − v_l^T Σ^{−1} v_m.
     *    - (NLO) Add O(1/V_c) corrections B_lm, C_lm.
     *
     * The LO covariance formula is the Schur complement of the extended
     * susceptibility matrix χ̃ from SAM Eq. (2.74), giving the SAM result
     * Eq. (2.86): χ_{pp}^{ce} = det(χ̃) / det(χ).
     *
     * \param computeCovariance  If true, compute the full N×N covariance matrix.
     */
    void ComputeAnalyticCumulants(bool computeCovariance = false);

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

    CanonicalMethod m_Method; ///< Method for computing the canonical partition functions

    int m_BMAX; ///< Maximum baryon number
    int m_QMAX; ///< Maximum electric charge
    int m_SMAX; ///< Maximum strangeness
    int m_CMAX; ///< Maximum charm
    int m_BMAX_list, m_QMAX_list, m_SMAX_list, m_CMAX_list;

    double m_MultExp; ///< Exponential multiplier for canonical partition function calculations
    double m_MultExpBanalyt; ///< Exponential multiplier for analytical baryon fugacity calculations

    /// Active GCE model pointer (set by PrepareModelGCE, cleared by CleanModelGCE).
    /// Points to either m_modelgce_default or m_modelgce_ext; never owned.
    ThermalModelBase *m_modelgce;
    ThermalModelBase *m_modelgce_default; ///< Default ideal-gas GCE model (owned, created once)
    ThermalModelBase *m_modelgce_ext;     ///< External GCE model (not owned), set via SetModelGCE()

    int m_BCE; ///< Flag indicating if baryon charge is conserved canonically
    int m_QCE; ///< Flag indicating if electric charge is conserved canonically
    int m_SCE; ///< Flag indicating if strangeness is conserved canonically
    int m_CCE; ///< Flag indicating if charm is conserved canonically

    /**
     * \brief Flag indicating whether the analytical calculation of baryon fugacity is used.
     */
    bool m_Banalyt;

    bool m_PartialZCalculated;
    bool m_PartialZCalculatedFlucts;

    //@{
    /// Saddle-point approximation data (for SaddlePoint quadrature method)
    int m_SaddlePointDim;    ///< Number of canonically conserved charges (d <= 4)
    std::vector<int> m_SaddlePointChargeIndices;  ///< Which of {B=0,Q=1,S=2,C=3} are canonical
    std::vector<double> m_SaddlePointMu;          ///< Saddle-point mu*_a/T for each conserved charge (size 4)
    std::vector<double> m_SaddlePointMuStar;      ///< Per-particle full chemical potential at saddle point
    double m_SaddlePointLogDetSigma;              ///< log(det(Sigma))
    std::vector<double> m_W2;                      ///< Per-particle W^(2) = Vc * chi2 (for SaddlePointNLO)

    /// NLO thermodynamic quantities computed by ComputeAnalyticCumulants().
    /// These are exact at O(1/Vc) including quantum statistics (cluster expansion).
    double m_NLOEnergyDensity;   ///< ε = ε^GCE(μ*) - (1/2) Σ_j E_j G_jj (LO: ε^GCE(μ*))
    double m_NLOPressure;        ///< P = p^GCE(μ*) - d·T/(2·Vc) (LO: p^GCE(μ*))
    double m_NLOEntropyDensity;  ///< LO entropy: s^GCE(μ*) (used only by SaddlePointLO; NLO uses lnZ_SP)
    //@}
  };

} // namespace thermalfist

#endif
