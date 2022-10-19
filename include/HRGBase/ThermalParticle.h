/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALPARTICLE_H
#define THERMALPARTICLE_H


/**
 * \file ThermalParticle.h
 * Contains classes and methods which 
 * provide basic information about a
 * particle species
 * and perform calculations
 * of various thermodynamic function
 * in the grand canonical ensemble.
 * 
 */


#include <string>
#include <vector>
#include <cmath>

#include "HRGBase/ParticleDecay.h"
#include "HRGBase/ThermalModelParameters.h"
#include "HRGBase/IdealGasFunctions.h"
#include "HRGBase/xMath.h"

namespace thermalfist {

  

  /**
   * \brief An auxiliary struct containing the list of conserved charges.
   *
   */
  struct ConservedCharge {
    /**
     * \brief Set of all conserved charges considered.
     *
     */
    enum Name {
      BaryonCharge = 0,      ///< Baryon number
      ElectricCharge = 1,    ///< Electric charge
      StrangenessCharge = 2, ///< Strangeness
      CharmCharge = 3        ///< Charm
    };
    static const int NumberOfTypes = 5;
  };

  

  /**
   * \brief Class containing all information about a particle specie.
   * 
   * Also contains implementation of calculation of various thermodynamic quantities
   * in an ideal gas in the grand canonical ensemble.
   */
  class ThermalParticle
  {
  public:
    /// Vector of all decay channels of a particle.
    typedef std::vector<ParticleDecayChannel> ParticleDecaysVector;

    /**
     * \brief Relativistic vs non-relativistic Breit-Wigner shape.
     *
     */
    enum ResonanceWidthShape {
      RelativisticBreitWigner,
      NonRelativisticBreitWigner
    };

    /**
     * \brief Treatment of finite resonance widths.
     *
     */
    enum ResonanceWidthIntegration {
      ZeroWidth,             ///< Zero-width approximation
      BWTwoGamma,            ///< Energy-independent Breit-Wigner in +-2\Gamma interval
      FullInterval,          ///< Energy-independent Breit-Wigner in full energy interval
      FullIntervalWeighted,  ///< Energy-independent Breit-Wigner in full energy interval with weighted branching ratios
      eBW,                   ///< Energy-dependent Breit-Wigner scheme (eBW)
      eBWconstBR             ///< Energy-dependent Breit-Wigner scheme (eBW) with constant branching ratios when evaluating feeddown
    };

    /**
     * \brief Construct a new ThermalParticle object
     *
     * \param Stable    Particle's stability flag
     * \param Name      Particle's name
     * \param PDGID     Particle's PDG ID
     * \param Deg       Particle's internal degeneracy
     * \param Stat      Statistics: 1 -- Fermi-Dirac, -1 -- Bose-Einstein, 0 - Maxwell-Boltzmann
     * \param Mass      Particle's mass
     * \param Strange   Particle's strangeness
     * \param Baryon    Particle's baryon number
     * \param Charge    Particle's electric charge
     * \param AbsS      Particle's strange quark content
     * \param Width     Particle's width
     * \param Threshold Particle's decays threshold
     * \param Charm     Particle's charm
     * \param AbsC      Particle's charm quark content
     * \param Quark     Particle's light quark content
     */
    ThermalParticle(bool Stable = true, std::string Name = "hadron", long long PDGID = 0, double Deg = 1., int Stat = 0, double Mass = 0.,
      int Strange = 0, int Baryon = 0, int Charge = 0, double AbsS = 0., double Width = 0., double Threshold = 0., int Charm = 0, double AbsC = 0., int Quark = 0);
    ~ThermalParticle(void);

    /**
     * \brief Fills coefficients for mass integration in the energy independent BW scheme
     * 
     */
    void FillCoefficients();

    /// Fills coefficients for mass integration in the eBW scheme
    void FillCoefficientsDynamical();

    /// Total width (eBW scheme) at a given mass
    double TotalWidtheBW(double M) const;

    /**
     * \brief (Energy-dependent) branching ratios
     * 
     * \param M Energy [GeV]
     * \param eBW Whether branching ratios are energy-dependent or not
     * \return std::vector<double> A vector of branching ratios for all decay channels
     */
    std::vector<double> BranchingRatiosM(double M, bool eBW = true) const;

    /**
     * \brief Mass distribution of a resonance in a thermal environment
     * 
     * Mass distribution of a resonance in a thermal environment (not normalized!).
     * Width is specified manually.
     * 
     * \param M Mass [GeV]
     * \param T Temperature [GeV]
     * \param Mu Chemical potential [GeV]
     * \param width Resonance width [GeV]
     * \return Mass distribution function 
     */
    double ThermalMassDistribution(double M, double T, double Mu, double width);

    /**
     * \brief Mass distribution of a resonance in a thermal environment
     * 
     * Mass distribution of a resonance in a thermal environment (not normalized!).
     * Energy-dependent width is computed automatically.
     * 
     * \param M Mass [GeV]
     * \param T Temperature [GeV]
     * \param Mu Chemical potential [GeV]
     * \return Mass distribution function (not normalized!)
     */
    double ThermalMassDistribution(double M, double T, double Mu);

    /**
     * \brief Normalizes all branching ratios such that
     *        they sum up to 100%.
     * 
     */
    void NormalizeBranchingRatios();

    /**
     * \brief Restores all branching ratios to the original values.
     * 
     */
    void RestoreBranchingRatios();

    /**
     * \brief Computes a specified ideal gas thermodynamic function.
     * 
     * Computes a specified ideal gas thermodynamic function.
     * Takes into account chemical non-equilibrium fugacity factors
     * and finite resonance widths.
     * 
     * \param params   Structure containing the temperature value and the chemical factors.
     * \param type     The type of the thermodynamic function calculated.
     * \param useWidth Whether finite widths are taken into account.
     * \param mu       Chemical potential.
     * \return         Value of the computed thermodynamic function.
     */
    double Density(const ThermalModelParameters &params, IdealGasFunctions::Quantity type = IdealGasFunctions::ParticleDensity, bool useWidth = 0, double mu = 0.) const;

    /**
     * Computes contribution of a single term in the cluster expansion
     * to the quantity which is to be computed by the Density() method.
     * 
     * \param n        Number of the term.
     * \param params   Structure containing the temperature value and the chemical factors.
     * \param type     The type of the thermodynamic function calculated.
     * \param useWidth Whether finite widths are taken into account.
     * \param mu       Chemical potential.
     * \return         Value of the computed term.
     */
    double DensityCluster(int n, const ThermalModelParameters &params, IdealGasFunctions::Quantity type = IdealGasFunctions::ParticleDensity, bool useWidth = 0, double mu = 0.) const;

    /**
     * \brief Computes the ideal gas generalized susceptibility \f$ \chi_n \equiv \frac{\partial^n p/T^4}{\partial (mu/T)^n} \f$.
     * 
     * Computes the generalized susceptibility \f$ \chi_n \equiv \frac{\partial^n p/T^4}{\partial (mu/T)^n} \f$
     * of the corresponding ideal gas.
     * Takes into account chemical non-equilibrium fugacity factors
     * and finite resonance widths.
     * 
     * \param index    Order of the susceptibility.
     * \param params   Structure containing the temperature value and the chemical factors.
     * \param useWidth Whether finite widths are taken into account.
     * \param mu       Chemical potential.
     * \return         Value of the computed susceptility.
     */
    double chi(int index, const ThermalModelParameters &params, bool useWidth = 0, double mu = 0.) const;

    /**
     * \brief Computes the ideal gas dimensionfull susceptibility \f$ \chi_n \equiv \frac{\partial^n p}{\partial mu^n} \f$.
     *
     * Computes the dimensionfull susceptibility \f$ \chi_n \equiv \frac{\partial^n p}{\partial mu^n} \f$
     * of the corresponding ideal gas.
     * Takes into account chemical non-equilibrium fugacity factors
     * and finite resonance widths.
     *
     * \param index    Order of the susceptibility.
     * \param params   Structure containing the temperature value and the chemical factors.
     * \param useWidth Whether finite widths are taken into account.
     * \param mu       Chemical potential.
     * \return         Value of the computed susceptility [GeV^{4-n].
     */
    double chiDimensionfull(int index, const ThermalModelParameters& params, bool useWidth = 0, double mu = 0.) const;

    /**
     * \brief Computes the scaled variance of particle number fluctuations
     *        in the ideal gas.
     * Computes the scaled variance (\chi_2 / \chi_1) of particle number fluctuations
     * in the ideal gas.
     * Takes into account chemical non-equilibrium fugacity factors
     * and finite resonance widths.
     * 
     * \param params   Structure containing the temperature value and the chemical factors.
     * \param useWidth Whether finite widths are taken into account.
     * \param mu       Chemical potential.
     * \return         Value of the computed scaled variance.
     */
    double ScaledVariance(const ThermalModelParameters &params, bool useWidth = 0, double mu = 0.) const;

    /**
     * \brief Computes the normalized skewness of particle number fluctuations
     *        in the ideal gas.
     * 
     * Computes the normalized skewness (\chi_3 / \chi_2) of particle number fluctuations
     * in the ideal gas.
     * Takes into account chemical non-equilibrium fugacity factors
     * and finite resonance widths.
     * 
     * \param params   Structure containing the temperature value and the chemical factors.
     * \param useWidth Whether finite widths are taken into account.
     * \param mu       Chemical potential.
     * \return         Value of the computed normalized skewness.
     */
    double Skewness(const ThermalModelParameters &params, bool useWidth = 0, double mu = 0.) const;

    /**
     * \brief Computes the normalized excess kurtosis of particle number fluctuations
     *        in the ideal gas.
     * 
     * Computes the normalized excess kurtosis (\chi_4 / \chi_2) of particle number fluctuations
     * in the ideal gas.
     * Takes into account chemical non-equilibrium fugacity factors
     * and finite resonance widths.
     * 
     * \param params   Structure containing the temperature value and the chemical factors.
     * \param useWidth Whether finite widths are taken into account.
     * \param mu       Chemical potential.
     * \return         Value of the computed normalized excess kurtosis.
     */
    double Kurtosis(const ThermalModelParameters &params, bool useWidth = 0, double mu = 0.) const;

    /**
     * \brief Fermi-Dirac distribution function.
     * 
     * \param k  Momentum [GeV]
     * \param T  Temperature [GeV]
     * \param mu Chemical potential [GeV]
     * \param m  Mass [GeV]
     * \return   Computed Fermi-Dirac function.
     */
    double FD(double k, double T, double mu, double m) const;
    
    /**
     * Computes the light quark content as follows:
     * |u,d| = 3 * |B| - |s| - |c|
     * where |B| is the absolute baryon number
     * and |s| and |c| is the stange and charm quark contents, respectively.
     * 
     * \return Computed light quark content.
     */
    double GetAbsQ() const;

    /**
     * \brief Get the quantum number numbered by the index
     * 
     * \param index 0 -- baryon number, 1 -- electric charge, 
     *              2 -- strangeness, 3 -- charm
     * \return Particle's quantum number 
     */
    double GetCharge(int index) const;

    /**
     * \brief Get the absolute value of a quantum number
     * 
     * \param index 0 -- absolute baryon number, 1 -- absolute electric charge, 
     *              2 -- strange quark content, 3 -- charm quark content
     * \return Particle's absolute value of a quantum number
     */
    double GetAbsCharge(int index) const;

    /**
     * \brief Whether particle is neutral one.
     * 
     * \return true  Particle is neutral (all quantum numbers are zero).
     * \return false Particle is not neutral (anti-particle exists).
     */
    bool IsNeutral() const;

    /// Return particle stability flag
    bool IsStable() const { return m_Stable; }

    /// Sets particle stability flag
    void SetStable(bool stable = true) { m_Stable = stable; }

    /// Whether particle is an antiparticle, i.e. its PDG ID is < 0
    bool IsAntiParticle() const { return m_AntiParticle; }

    /// Set manually whether particle is an antiparticle
    void SetAntiParticle(bool antpar = true) { m_AntiParticle = antpar; }

    /// Particle's name
    const std::string& Name() const { return m_Name; }

    /// Set particle's name
    void SetName(const std::string &name) { m_Name = name; }

    /// Particle's Particle Data Group (PDG) ID number
    long long  PdgId() const { return m_PDGID; }

    /// Set particle's particle's Particle Data Group (PDG) ID number
    void SetPdgId(long long PdgId) { m_PDGID = PdgId; }

    /// Particle's internal degeneracy factor
    double Degeneracy() const { return m_Degeneracy; }

    /// Set particle's internal degeneracy factor
    void SetDegeneracy(double deg) { m_Degeneracy = deg; }

    /**
     * \brief Particle's statistics
     * 
     * 1 -- Fermi-Dirac,
     * -1 -- Bose-Einstein,
     * 0 - Maxwell-Boltzmann
     * 
     * \return Particle's statistics 
     */
    int  Statistics() const { return m_Statistics; }

    /**
     * \brief Set particle's statistics
     * 
     * 1 -- Fermi-Dirac
     * -1 -- Bose-Einstein
     * 0 - Maxwell-Boltzmann
     * 
     * \param stat Statistics
     */
    void SetStatistics(int stat) { m_Statistics = stat; }

    /**
     * \brief Use quantum statistics
     * 
     * \param enable true -- use quantum statistics
     *               false -- use Maxwell-Boltzmann statistics
     */
    void UseStatistics(bool enable);

    /// Particle's mass [GeV]
    double Mass() const { return m_Mass; }

    /// Set particle's mass [GeV]
    void SetMass(double mass);// { m_Mass = mass; }

    /// Particle's baryon number
    int BaryonCharge() const { return m_Baryon; }

    /// Set particle's baryon number
    void SetBaryonCharge(int chg) { m_Baryon = chg; SetAbsoluteQuark(GetAbsQ()); }

    /// Particle's electric charge
    int ElectricCharge() const { return m_ElectricCharge; }

    /// Set particle's electric charge
    void SetElectricCharge(int chg) { m_ElectricCharge = chg; }

    /// Particle's strangeness
    int Strangeness() const { return m_Strangeness; }
    /// Set particle's strangeness
    void SetStrangenessCharge(int chg) { m_Strangeness = chg; }

    /// Particle's charm
    int Charm() const { return m_Charm; }

    /// Set particle's charm
    void SetCharm(int chg) { m_Charm = chg; }

    /// One of the four QCD conserved charges
    int ConservedCharge(ConservedCharge::Name chg) const;

    /**
     * \brief Arbitrary (auxiliary) charge assigned to particle
     *
     * \return Arbitrary (auxiliary) charge
     */
    double ArbitraryCharge() const { return m_ArbitraryCharge; }

    /**
     * \brief Assigns arbitrary (auxiliary) charge to particle
     * 
     * \param Arbitrary (auxiliary) charge 
     */
    void SetArbitraryCharge(double arbchg) { m_ArbitraryCharge = arbchg; }

    /// Absolute light quark content |u,d|
    double AbsoluteQuark() const { return m_AbsQuark; }

    /// Set absolute light quark content |u,d|
    void SetAbsoluteQuark(double abschg) { m_AbsQuark = abschg; }

    /// Absolute strange quark content |s|
    double AbsoluteStrangeness() const { return m_AbsS; }

    /// Set absolute strange quark content |s|, light quark content then re-evaluted
    void SetAbsoluteStrangeness(double abschg) { m_AbsS = abschg; SetAbsoluteQuark(GetAbsQ()); }

    /// Absolute charm quark content |s|
    double AbsoluteCharm() const { return m_AbsC; }

    /// Set absolute charm quark content |s|, light quark content then re-evaluted
    void SetAbsoluteCharm(double abschg) { m_AbsC = abschg; SetAbsoluteQuark(GetAbsQ()); }

    /// Whether zero-width approximation is enforced for this particle species
    bool ZeroWidthEnforced() const;

    /// Particle's width at the pole mass (GeV)
    double ResonanceWidth() const { return m_Width; }

    /**
     * \brief Sets the particle's width at the pole mass
     * 
     * If width is non-zero,
     * the coefficients used for mass integration
     * are re-evaluated
     * 
     * \param width Width (GeV)
     */
    void SetResonanceWidth(double width);

    /**
     * \brief The decays threshold mass
     * 
     * The threshold mass for calculation in the
     * energy-independent Breit-Wigner scheme
     * 
     * \return Threshold mass 
     */
    double DecayThresholdMass() const { return m_Threshold; }
    /**
     * \brief Set the decays threshold mass
     * 
     * If width is non-zero,
     * the coefficients used for mass integration
     * in the energy independent scheme are re-evaluated
     * 
     * \param threshold Threshold mass (GeV)
     */
    void SetDecayThresholdMass(double threshold);

    

    /// Returns threshold mass as the minimum
    /// threshold among all the decay channels.
    /// Used in the eBW scheme
    double DecayThresholdMassDynamical() const { return m_ThresholdDynamical; }

    /// Set the threshold mass manually for use in the eBW scheme. 
    void SetDecayThresholdMassDynamical(double threshold);

    /// Evaluate the threshold mass as the minimum
    /// threshold among all the decay channels
    void CalculateAndSetDynamicalThreshold();

    /**
     * \brief Resonance width profile in use
     * 
     * Can be relativistic or non-relativistic Breit-Wigner
     * 
     * \return ResonanceWidthShape Width profile used
     */
    ResonanceWidthShape GetResonanceWidthShape() const { return m_ResonanceWidthShape; }

    /**
     * \brief Set the resonance width profile to use
     * 
     * \param shape Relativistic or non-relativistic Breit-Wigner
     */
    void SetResonanceWidthShape(ResonanceWidthShape shape);

    /**
     * \brief Resonance width integration scheme used to treat
     *        finite resonance widths
     * 
     * \return ResonanceWidthIntegration 
     */
    ResonanceWidthIntegration GetResonanceWidthIntegrationType() const { return m_ResonanceWidthIntegrationType; }
    
    /**
     * \brief Set the ResonanceWidthIntegration scheme used to treat
     *        finite resonance widths
     * 
     * \param type ResonanceWidthIntegration scheme
     */
    void SetResonanceWidthIntegrationType(ResonanceWidthIntegration type);

    /**
     * Resonance mass distribution: Relativistic or non-relativistic Breit-Wigner
     * evaluated at the given mass m (GeV) and pole mass's width
     * 
     * \param m Mass (GeV)
     * \return Mass distribution 
     */
    double MassDistribution(double m) const;

    // Resonance mass distribution with manually input width
    double MassDistribution(double m, double width) const;

    /**
     * \brief Particle's weight
     * 
     * Multiplies the degeneracy factor, equal to one by default.
     * Currently not used.
     * 
     * \return Weight. 
     */
    double Weight() const { return m_Weight; }

    /// Set particle's weight factor
    void SetWeight(double weight) { m_Weight = weight; }

    /**
     * \brief Decay type of the particle.
     * 
     * \return Decay type of the particle.
     */
    ParticleDecayType::DecayType DecayType() const { return m_DecayType; }

    /// Set particle's Decay Type
    void SetDecayType(ParticleDecayType::DecayType type) { m_DecayType = type; }

    /**
     * \brief A vector of particle's decays
     * 
     * A vector of ParticleDecay objects corresponding to
     * all decay channels of the particle.
     * 
     * \return const std::vector<ParticleDecay>& 
     */
    const ParticleDecaysVector& Decays() const { return m_Decays; }

    /// Returns a non-const reference to Decays()
    ParticleDecaysVector& Decays() { return m_Decays; }

    /**
     * \brief Set the Decays vector
     * 
     * Sets all decays of the particle
     * 
     * \param Decays ParticleDecay vector containing all particle decays 
     */
    void SetDecays(const ParticleDecaysVector &Decays) { m_Decays = Decays; }

    /// Remove all decays
    void ClearDecays() { m_Decays.resize(0); }

    //@{
    /// A backup copy of particle's decays
    const ParticleDecaysVector& DecaysOriginal() const { return m_DecaysOrig; }
    ParticleDecaysVector& DecaysOriginal() { return m_DecaysOrig; }
    void SetDecaysOriginal(const ParticleDecaysVector &DecaysOrig) { m_DecaysOrig = DecaysOrig; }
    //@}

    /// Read decays from a file and assign them to the particle
    void ReadDecays(std::string filename = "");

    
    /**
     * \brief Computes average decay branching ratios
     *        by integrating over the thermal mass distribution.
     * 
     * To be later used when evaluating feeddown contributions.
     * 
     * \param params   Structure containing the temperature value and the chemical factors.
     * \param useWidth Whether finite widths are taken into account.
     * \param mu       Chemical potential.
     */
    void CalculateThermalBranchingRatios(const ThermalModelParameters &params, bool useWidth = 0, double mu = 0.);

    /**
     * \brief Sets the CalculationType() method to evaluate quantum statistics.
     * 
     * \param type Method to evaluate quantum statistics.
     */
    void SetCalculationType(IdealGasFunctions::QStatsCalculationType type) { m_QuantumStatisticsCalculationType = type; }
    
    /**
     * \brief Method to evaluate quantum statistics.
     * 
     * Cluster expansion or numerical integration
     * using the quadratures.
     * 
     * \return IdealGasFunctions::QStatsCalculationType 
     */
    IdealGasFunctions::QStatsCalculationType CalculationType()       const { return m_QuantumStatisticsCalculationType; }

    /**
     * \brief Number of terms in the cluster expansion method.
     * 
     * \return Number of terms in the cluster expansion method.
     */
    int ClusterExpansionOrder() const { return m_ClusterExpansionOrder; }

    /// Set ClusterExpansionOrder()
    void SetClusterExpansionOrder(int order) { m_ClusterExpansionOrder = order; }

    std::vector<double> BranchingRatioWeights(const std::vector<double> & ms) const;

    const std::vector<double>& Nch() const { return m_Nch; }
    std::vector<double>&  Nch() { return m_Nch; }

    const std::vector<double>& DeltaNch() const { return m_DeltaNch; }
    std::vector<double>&  DeltaNch() { return m_DeltaNch; }

    /**
     * \brief Generates the anti-particle to the current particle specie
     * 
     * Note: Decay channels of anti-particle 
     * are NOT generated by this method and have to
     * be set elsewhere.
     * 
     * \return ThermalParticle Antiparticle
     */
    ThermalParticle GenerateAntiParticle(/*ThermalParticleSystem *TPS = NULL*/) const;

    bool operator==(const ThermalParticle &rhs) const; // TODO: improve
    bool operator!=(const ThermalParticle &rhs) const { return !(*this == rhs); }

  private:
    /**
    *  Auxiliary coefficients used for numerical integration using quadratures
    */
    std::vector<double> m_xlag32, m_wlag32;
    std::vector<double> m_xleg, m_wleg;
    std::vector<double> m_xleg32, m_wleg32;
    std::vector<double> m_brweight;


    /**
    *  For the eBW scheme
    */
    std::vector<double> m_xlegdyn, m_wlegdyn, m_vallegdyn;
    std::vector<double> m_xlegpdyn, m_wlegpdyn, m_vallegpdyn;
    std::vector<double> m_xlagdyn, m_wlagdyn, m_vallagdyn;

    std::vector<double> m_xalldyn, m_walldyn, m_densalldyn;


    bool m_Stable;                /**< Flag whether particle is marked stable. */
    ParticleDecayType::DecayType m_DecayType;        /**< Type wrt to decay: Stable, Default (placeholder), Weak, Electromagnetic, Strong */

    bool m_AntiParticle;          /**< Whether particle was created as an antiparticle to another one. */
    std::string m_Name;           /**< Particle name. */
    long long m_PDGID;            /**< PDG (HEP) ID of a particle. */
    double m_Degeneracy;          /**< (Spin) Degeneracy factor. */
    int m_Statistics;             /**< Statistics used (Bolzmann or Quantum). */
    int m_StatisticsOrig;         /**< Particle's original Fermi/Bose statistics. */
    double m_Mass;                /**< Mass (GeV) */

    /**
     *   0 - Cluster expansion (default), 1 - Numerical integration
     */
    IdealGasFunctions::QStatsCalculationType m_QuantumStatisticsCalculationType;

    /**
     *   Number of terms in cluster expansion.
     *   Default is 10 for m < 200 MeV (pions), 5 for m < 1000 MeV, and 3 otherwise
     *
     */
    int m_ClusterExpansionOrder;

    int m_Baryon;                 /**< Baryon charge */
    int m_ElectricCharge;         /**< Electric charge */
    int m_Strangeness;            /**< Strangeness charge */
    int m_Charm;                  /**< Charm charge */
    int m_Quark;                  /**< Absolute quarks content */

    double m_ArbitraryCharge;     /**< Arbitrary charge associated with particle (external) */
    double m_AbsQuark;            /**< Absolute light quark content */
    double m_AbsS;                /**< Absolute strangeness content */
    double m_AbsC;                /**< Absolute charm content */

    double m_Width;               /**< Resonance width (GeV) */
    double m_Threshold;           /**< Lower decay threshold (GeV) */
    double m_ThresholdDynamical;  /**< Lower decay threshold (GeV) calculated from daughter masses */
    ResonanceWidthShape m_ResonanceWidthShape;                  /**< Either relativistic or non-relativitic Breit-Wigner */
    ResonanceWidthIntegration m_ResonanceWidthIntegrationType;  /**< Plus-minus TwoGamma or from m0 to infty */
    double m_Radius;              /**< Hard-core radius (fm) */
    double m_Vo;                  /**< Eigenvolume parameter (fm^3). Obsolete. To be removed. */
    double m_Weight;              /**< Weight of a given particle. Default is 1 */

    ParticleDecaysVector m_Decays;      /**< All decay channels currently in use.  */

    /**
     *   All original decay channels. Needed to switch between renormalizing decay branching ratios and not.
     */
    ParticleDecaysVector m_DecaysOrig;


    /**
     *   For calculating final state charged particle multiplicities
     *   Arrays contain Nch, N+, N-
     */
    std::vector<double> m_Nch;
    std::vector<double> m_DeltaNch;

    /// Whether BEC was encountered
    bool m_LastDensityOk;
  };

} // namespace thermalfist

#endif
