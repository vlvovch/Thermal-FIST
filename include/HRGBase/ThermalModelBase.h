/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELBASE_H
#define THERMALMODELBASE_H

#include <string>

#include "HRGBase/ThermalParticleSystem.h"
#include "HRGBase/xMath.h"
#include "HRGBase/Broyden.h"


namespace thermalfist {

  /**
   * \brief Abstract base class for an HRG model implementation.
   * 
   * Contains the base implementation of an arbitrary HRG model.
   * The actual calculations of thermodynamic functions are implemented
   * in derived classes.
   * 
   */
  class ThermalModelBase
  {
  public:
    /**
     * \brief The list of statistical ensembles.
     * 
     */
    enum ThermalModelEnsemble { 
      GCE = 0, ///< Grand canonical ensemble
      CE = 1,  ///< Canonical ensemble
      SCE = 2, ///< Strangeness-canonical ensemble
      CCE = 3 ///< Charm-canonical ensemble
    };

    /**
     * \brief Type of interactions included in the HRG model.
     * 
     */
    enum ThermalModelInteraction { 
      Ideal = 0,        ///< Ideal HRG model
      DiagonalEV = 1,   ///< Diagonal excluded volume model
      CrosstermsEV = 2, ///< Crossterms excluded volume model
      QvdW = 3,         ///< Quantum van der Waals model
      RealGas = 4,      ///< Real gas model. Not yet fully implemented.
      MeanField = 5     ///< Mean field model. Not yet fully implemented.
    };

    /**
     * \brief Construct a new ThermalModelBase object.
     * 
     * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
     * \param params ThermalModelParameters object with current thermal parameters
     */
    ThermalModelBase(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    virtual ~ThermalModelBase(void) { }

    /// Number of different particle species in the list
    int ComponentsNumber() const { return static_cast<int>(m_densities.size()); }

    /**
     * \brief Fills the excluded volume coefficients \f$ \tilde{b}_{ij} \f$ based on the 
     *        provided radii parameters for all species.
     * 
     * Fills the coefficients in accordance with Eqs. (5) and (7)
     * here https://arxiv.org/pdf/1606.06218.pdf
     * 
     * \param ri A vector with radii parameters for all species.
     *           0-based indices of the vector must corresponds to the
     *           0-based indices of the particle list TPS()
     */
    virtual void FillVirial(const std::vector<double> & ri = std::vector<double>(0));

    /// Whether finite resonance widths are considered.
    bool UseWidth() const { return m_UseWidth; }

    /**
     * \brief Sets whether finite resonance widths are used. Deprecated.
     * 
     * If the widths are to be used, the ThermalParticle::ResonanceWidthIntegration::BWTwoGamma
     * scheme will be used.
     * 
     * \param useWidth Whether to use the finite resonance widths.
     */
    void SetUseWidth(bool useWidth);

    /**
     * \brief Sets the finite resonance widths scheme to use.
     * 
     * \param type ThermalParticle::ResonanceWidthIntegration scheme
     */
    void SetUseWidth(ThermalParticle::ResonanceWidthIntegration type);

    //@{
      /**
       * \brief Whether branching ratios are renormalized to 100%.
       * 
       * \param normBratio  Whether branching ratios shoul be renormalized to 100%.
       */
    void SetNormBratio(bool normBratio);
    bool NormBratio() const { return m_NormBratio; }
    //@}

    /// OpenMP support. Currently not used.
    void SetOMP(bool openMP) { m_useOpenMP = openMP; }
     
    //@{
      /**
       * \brief The thermal parameters
       * 
       * \param params The thermal parameters to be used 
       */
    virtual void SetParameters(const ThermalModelParameters& params);
    const ThermalModelParameters& Parameters() const { return m_Parameters; }
    //@}

    /**
     * \brief Calls SetParameters() with current m_Parameters.
     * 
     */
    void UpdateParameters() { SetParameters(m_Parameters); }

    /**
     * \brief Set the temperature
     * 
     * \param T Temperature (GeV)
     */
    virtual void SetTemperature(double T);

    /**
     * \brief Set the baryon chemical potential
     * 
     * \param muB Baryon chemical potential (GeV)
     */
    virtual void SetBaryonChemicalPotential(double muB);

    /**
     * \brief Set the electric chemical potential
     * 
     * \param muQ Electric chemical potential (GeV)
     */
    virtual void SetElectricChemicalPotential(double muQ);

    /**
     * \brief Set the strangeness chemical potential
     * 
     * \param muS Strangeness chemical potential (GeV)
     */
    virtual void SetStrangenessChemicalPotential(double muS);

    /**
     * \brief Set the charm chemical potential
     * 
     * \param muC Charm chemical potential (GeV)
     */
    virtual void SetCharmChemicalPotential(double muC);

    /**
     * \brief Set the light quark fugacity factor
     * 
     * \param gammaq Light quark fugacity factor
     */
    virtual void SetGammaq(double gammaq);

    /**
     * \brief Set the strange quark fugacity factor
     * 
     * \param gammaS Strange quark fugacity factor
     */
    virtual void SetGammaS(double gammaS);

    /**
     * \brief Set the charm quark fugacity factor
     * 
     * \param gammaC Charm quark fugacity factor
     */
    virtual void SetGammaC(double gammaC);

    /**
     * \brief Set the total baryon number (for canonical ensemble only)
     * 
     * \param B Total baryon number
     */
    virtual void SetBaryonCharge(int B);

    /**
     * \brief Set the total electric charge (for canonical ensemble only)
     * 
     * \param B Total electric charge
     */
    virtual void SetElectricCharge(int Q);

    /**
     * \brief Set the total strangeness (for canonical ensemble only)
     * 
     * \param B Total strangeness
     */
    virtual void SetStrangeness(int S);

    /**
     * \brief Set the total charm (for canonical ensemble only)
     * 
     * \param B Total charm
     */
    virtual void SetCharm(int C);

    /**
     * \brief Set the same excluded volume radius parameter for all species.
     * 
     * \param rad Radius parameter (fm)
     */
    virtual void SetRadius(double /*rad*/) { }

    /**
     * \brief Set the radius parameter for particle species i
     * 
     * \param i 0-based index of particle species
     * \param rad Radius parameter (fm)
     */
    virtual void SetRadius(int /*i*/, double /*rad*/) { }

    /**
     * \brief Set the excluded volume coefficient \f$ \tilde{b}_{ij} \f$
     * 
     * Excluded parameter for repulsive interaction between particle
     * species i and j. 
     * 
     * \param i 0-based index of the first particle species
     * \param j 0-based index of the second particle species
     * \param b Excluded volume parameter \f$ \tilde{b}_{ij} \f$ (fm\f$^3\f$)
     */
    virtual void SetVirial(int /*i*/, int /*j*/, double /*b*/) { }

    /**
     * \brief Set the vdW mean field attraction coefficient \f$ a_{ij} \f$
     * 
     * \param i 0-based index of the first particle species
     * \param j 0-based index of the second particle species
     * \param a vdW mean field attraction parameter \f$ a_{ij} \f$ (GeV fm\f$^3\f$) 
     */
    virtual void SetAttraction(int /*i*/, int /*j*/, double /*a*/) { }

    //@{
    /// Switches off excluded volume terms for all meson-meson pairs
    /// i.e. \f$ \tilde{b}_{ij} = 0 \f$ for meson-meson pairs
    virtual void DisableMesonMesonVirial();
    void DisableMesonMesonRepulsion() { return DisableMesonMesonVirial(); }
    //@}

    /// Switches off QvdW attraction terms for all meson-meson pairs
    /// i.e. \f$ a_{ij} = 0 \f$ for meson-meson pairs
    virtual void DisableMesonMesonAttraction();

    //@{
    /// Switches off excluded volume terms for all meson-baryon pairs
    /// i.e. \f$ \tilde{b}_{ij} = 0 \f$ for meson-baryon pairs
    virtual void DisableMesonBaryonVirial();
    void DisableMesonBaryonRepulsion() { return DisableMesonBaryonVirial(); }
    //@}

    /// Switches off QvdW attraction terms for all meson-baryon pairs
    /// i.e. \f$ a_{ij} = 0 \f$ for meson-baryon pairs
    virtual void DisableMesonBaryonAttraction();

    //@{
    /// Switches off excluded volume terms for all baryon-baryon pairs
    /// i.e. \f$ \tilde{b}_{ij} = 0 \f$ for baryon-baryon pairs
    virtual void DisableBaryonBaryonVirial();
    void DisableBaryonBaryonRepulsion() { return DisableBaryonBaryonVirial(); }
    //@}

    /// Switches off QvdW attraction terms for all baryon-baryon pairs
    /// i.e. \f$ a_{ij} = 0 \f$ for baryon-baryon pairs
    virtual void DisableBaryonBaryonAttraction();

    //@{
    /// Switches off eigenvolume terms for all baryon-antibaryon pairs
    /// i.e. \f$ \tilde{b}_{ij} = 0 \f$ for baryon-antibaryon pairs
    virtual void DisableBaryonAntiBaryonVirial();
    void DisableBaryonAntiBaryonRepulsion() { return DisableBaryonAntiBaryonVirial(); }
    //@}

    /// Switches off QvdW attraction terms for all baryon-antibaryon pairs
    /// i.e. \f$ a_{ij} = 0 \f$ for baryon-antibaryon pairs
    virtual void DisableBaryonAntiBaryonAttraction();

    /**
     * \brief Reads the QvdW interaction parameters from a file.
     * 
     * Actual implementation is in a derived class.
     * 
     * \param filename File with interaction parameters.
     */
    virtual void ReadInteractionParameters(const std::string & /*filename*/) { }

    /**
     * \brief Write the QvdW interaction parameters to a file.
     * 
     * Actual implementation is in a derived class.
     * 
     * \param filename Output file.
     */
    virtual void WriteInteractionParameters(const std::string & /*filename*/) { }

    /**
     * \brief Change the particle list.
     * 
     * \param TPS A pointer to new particle list.
     */
    virtual void ChangeTPS(ThermalParticleSystem *TPS);

    //@{
      /**
       * \brief Excluded volume coefficient \f$ \tilde{b}_{ij} = 0 \f$
       * 
       * Excluded parameter for repulsive interaction between particle
       * species i and j.
       * 
       * \param i 0-based index of the first particle species
       * \param j 0-based index of the second particle species
       * \return double Coefficient \f$ \tilde{b}_{ij} = 0 \f$
       */
    virtual double VirialCoefficient(int /*i*/, int /*j*/) const { return 0.; }
    double RepulsionCoefficient(int i, int j) const { return VirialCoefficient(i,j); }
    //@}

    /**
     * \brief QvdW mean field attraction coefficient \f$ a_{ij} \f$
     * 
     * \param i 0-based index of the first particle species
     * \param j 0-based index of the second particle species
     * \return double Coefficient \f$ a_{ij} = 0 \f$
     */
    virtual double AttractionCoefficient(int /*i*/, int /*j*/) const { return 0.; }

    /// Whether quantum statistics is used, 
    /// 0 - Boltzmann, 1 - Quantum
    bool QuantumStatistics() const { return m_QuantumStats; }

    /// Set whether quantum statistics is used, 
    /// 0 - Boltzmann, 1 - Quantum
    virtual void SetStatistics(bool stats);

    /**
     * \brief Sets the CalculationType() method to evaluate quantum statistics. 
     *        Calls the corresponding method in TPS().
     * 
     * \param type Method to evaluate quantum statistics.
     */
    virtual void SetCalculationType(IdealGasFunctions::QStatsCalculationType type) { m_TPS->SetCalculationType(type); }
    
    /**
     * \brief Set the number of terms in the cluster expansion method.
     *        Calls the corresponding method in TPS().
     * 
     * Sets the same value for all particles.
     * To set individually for each particle
     * use ThermalParticle::SetClusterExpansionOrder().
     * 
     * \param order Number of terms.
     */
    virtual void SetClusterExpansionOrder(int order) { m_TPS->SetClusterExpansionOrder(order); }
    
    /**
     * \brief Set the ThermalParticle::ResonanceWidthShape for all particles.
     *        Calls the corresponding method in TPS().
     * 
     * This global setting can be overriden for individual particles
     * by ThermalParticle::SetResonanceWidthShape().
     * 
     * \param shape ThermalParticle::ResonanceWidthShape
     */
    void SetResonanceWidthShape(ThermalParticle::ResonanceWidthShape shape) { m_TPS->SetResonanceWidthShape(shape); }
    
    /**
     * \brief Set the ThermalParticle::ResonanceWidthIntegration scheme for
     *        all particles.
     *        Calls the corresponding method in TPS().
     * 
     * \param type ThermalParticle::ResonanceWidthIntegration
     */
    void SetResonanceWidthIntegrationType(ThermalParticle::ResonanceWidthIntegration type);// { m_TPS->SetResonanceWidthIntegrationType(type); }

    /**
     * \brief Sets the chemical potentials of all particles.
     * 
     * Uses the current values of \f$ \mu_B,\,\mu_Q,\,\mu_S,\,\mu_Q \f$ 
     * to set \f$ \mu_i \f$.
     * 
     */
    virtual void FillChemicalPotentials();

    /**
     * \brief Sets the chemical potentials of all particles.
     *      
     * 
     * \param chem A vector with chemical potentials of all species.
     *             0-based indices of the vector must corresponds to the
     *             0-based indices of the particle list TPS()
     */
    virtual void SetChemicalPotentials(const std::vector<double> & chem = std::vector<double>(0));

    /**
     * \brief A vector of chemical potentials of all particles.
     * 
     * \return const std::vector<double>& A vector of chemical potentials of all particles
     */
    const std::vector<double>& ChemicalPotentials() const { return m_Chem; }

    /**
     * \brief Chemical potential of particle species i.
     * 
     * \param i 0-based index of particle species
     * \return double Chemical potential of particle species i
     */
    double ChemicalPotential(int i) const;

    /**
     * \brief Sets the chemical potential of  particle species i.
     *
     *
     * \param i    0-based index of particle species
     * \param chem value of the chemical potential
     */
    virtual void SetChemicalPotential(int i, double chem);

    /**
     * \brief Chemical potential entering the ideal gas expressions of particle species i.
     *
     * Includes chemical non-equilibrium and EV/vdW effects
     *
     * \param i 0-based index of particle species
     * \return double Chemical potential of particle species i
     */
    virtual double FullIdealChemicalPotential(int i) const;


    /// Whether the baryon chemical potential is to be constrained
    /// by a fixed entropy per baryon ratio SoverB().
    bool ConstrainMuB() const { return m_ConstrainMuB; }

    /// Sets whether the baryon chemical potential is to be constrained
    /// by a fixed entropy per baryon ratio SoverB().
    void ConstrainMuB(bool constrain) { m_ConstrainMuB = constrain; }

    /// Whether the electric chemical potential is to be constrained
    /// by a fixed electric-to-baryon ratio QoverB().
    bool ConstrainMuQ() const { return m_ConstrainMuQ; }

    /// Sets whether the electric chemical potential is to be constrained
    /// by a fixed electric-to-baryon ratio QoverB().
    void ConstrainMuQ(bool constrain) { m_ConstrainMuQ = constrain; }

    /// Whether the strangeness chemical potential is to be constrained
    /// by the condition of strangeness neutrality
    bool ConstrainMuS() const { return m_ConstrainMuS; }

    /// Sets whether the strangeness chemical potential is to be constrained
    /// by the condition of strangeness neutrality
    void ConstrainMuS(bool constrain) { m_ConstrainMuS = constrain; }

    /// Whether the charm chemical potential is to be constrained
    /// by the condition of charm neutrality
    bool ConstrainMuC() const { return m_ConstrainMuC; }

    /// Sets whether the charm chemical potential is to be constrained
    /// by the condition of charm neutrality
    void ConstrainMuC(bool constrain) { m_ConstrainMuC = constrain; }

    /// Sets whether partial chemical equilibrium with additional chemical potentials is used
    void UsePartialChemicalEquilibrium(bool usePCE) { m_PCE = usePCE; }

    /// Whether partial chemical equilibrium with additional chemical potentials is used
    bool UsePartialChemicalEquilibrium() { return m_PCE; }

    //@{
      /**
       * \brief The entropy per baryon ratio
       *        to be used to constrain the 
       *        baryon chemical potential.
       * 
       * \param SB The entropy per baryon
       */
    void SetSoverB(double SB) { m_SBgoal = SB; }
    double SoverB() const { return m_SBgoal; }
    //@}

    //@{
      /**
       * \brief The electric-to-baryon charge ratio
       *        to be used to constrain the 
       *        electric chemical potential.
       * 
       * \param QB The electric-to-baryon charge ratio
       */
    void SetQoverB(double QB) { m_QBgoal = QB; }
    double QoverB() const { return m_QBgoal; }
    //@}

    /**
     * \brief Sets the system volume
     * 
     * \param Volume System volume (fm\f$^3\f$)
     */
    void SetVolume(double Volume) { m_Volume = Volume; m_Parameters.V = Volume; }
    /// System volume (fm\f$^3\f$)
    double Volume() const { return m_Parameters.V; }
    
    /**
     * \brief Sets the system radius
     * 
     * The system volume is computed as
     * \f$V = \frac{4\pi}{3} \, R^3\f$
     * 
     * \param Volume System radius (fm)
     */
    void SetVolumeRadius(double R) { m_Volume = 4. / 3.*xMath::Pi() * R * R * R; m_Parameters.V = m_Volume; }

    /// The canonical correlation volume V\f$_c\f$ (fm\f$^3\f$)
    double CanonicalVolume() const { return m_Parameters.SVc; }
    /**
     * \brief Set the canonical correlation volume V\f$_c\f$ 
     * 
     * \param Volume Canonical correlation volume (fm\f$^3\f$)
     */

    void SetCanonicalVolume(double Volume) { m_Parameters.SVc = Volume; }

    /**
     * \brief Set the canonical correlation system radius
     * 
     * The canonical correlation volume is computed as
     * \f$V_c = \frac{4\pi}{3} \, R_c^3\f$
     * 
     * \param Radius 
     */
    void SetCanonicalVolumeRadius(double Rc) { m_Parameters.SVc = 4. / 3. * xMath::Pi() * Rc * Rc * Rc; }


    //double StrangenessCanonicalVolume() const { return CanonicalVolume(); }
    //void SetStrangenessCanonicalVolume(double Volume) { SetCanonicalVolume(Volume); }
    //void SetStrangenessCanonicalVolumeRadius(double Radius) { SetCanonicalVolumeRadius(Radius); }

    // Same as FixParameters but with a more clear name on what is actually does
    /**
     * \brief Constrains the chemical potentials \f$ \mu_B,\,\mu_Q,\,\mu_S,\,\mu_C \f$ 
     *        by the conservation laws imposed.
     * 
     * This procedure uses the Broyden's method to solve
     * the system of equations corresponding to the
     * conservation laws.
     * The actual implementation of the procedure
     * is in FixParameters() and FixParametersNoReset() methods.
     *
     * \param resetInitialValues Whether initial guess values for 
     *        \f$ \mu_B,\,\mu_Q,\,\mu_S,\,\mu_C \f$ 
     *        are reset or current values will be used
     * 
     */
    void ConstrainChemicalPotentials(bool resetInitialValues = true);

    
    /**
     * \brief Method which actually implements ConstrainChemicalPotentials()
     *        (for backward compatibility).
     * 
     */
    virtual void FixParameters();

    /**
     * \brief Method which actually implements ConstrainChemicalPotentialsNoReset()
     *        (for backward compatibility).
     * 
     */
    virtual void FixParametersNoReset();

    /**
     * \brief The procedure which calculates the chemical potentials
     *        \f$ \mu_B,\,\mu_Q,\,\mu_S,\,\mu_Q \f$ which reproduce
     *        the specified total baryon, electric, strangeness, and charm
     *        charges of the system.
     * 
     * Uses the Broyden's method to determine
     * the chemical potentials.
     * The resulting chemical potentials are
     * stored in Parameters().
     * Useful for canonical ensemble applications.
     * 
     * \param totB The total desired net baryon charge of the system.
     * \param totQ The total desired electric charge of the system.
     * \param totS The total desired net strangeness of the system.
     * \param totC The total desired net charm of the system.
     * \param muBinit The initial guess for the baryon chemical potential.
     * \param muQinit The initial guess for the electric chemical potential.
     * \param muSinit The initial guess for the strangeness chemical potential.
     * \param muCinit The initial guess for the charm chemical potential.
     * \param ConstrMuB Whether the baryon chemical potential should be constrained.
     * \param ConstrMuQ Whether the electric chemical potential should be constrained.
     * \param ConstrMuS Whether the strangeness chemical potential should be constrained.
     * \param ConstrMuC Whether the charm chemical potential should be constrained.
     *
     * \return true is chemical potentials were contrained successfully, false otherwise
     */
    virtual bool SolveChemicalPotentials(double totB = 0., double totQ = 0., double totS = 0., double totC = 0.,
      double muBinit = 0., double muQinit = 0., double muSinit = 0., double muCinit = 0.,
      bool ConstrMuB = true, bool ConstrMuQ = true, bool ConstrMuS = true, bool ConstrMuC = true);

    /**
    * \brief Calculates the primordial
    *        densities of all species.
    *
    */
    virtual void CalculatePrimordialDensities() = 0;

    /**
     * \brief Calculates the primordial and total (after decays)
     *        densities of all species.
     * 
     */
    virtual void CalculateDensities();

    /**
     * \brief Checks whether issues have occured
     *        during the calculation of particle densities
     *        in the CalculateDensities() method.
     * 
     */
    virtual void ValidateCalculation();

    /**
     * \brief All messaged which occured during the validation
     *        procedure in the ValidateCalculation() method.
     * 
     * \return std::string List of messages, one per line
     */
    std::string ValidityCheckLog() const { return m_ValidityLog; }

    /**
     * \brief Calculates the particle densities in a grand-canonical ensemble.
     * 
     * A non-GCE based derived class will override this method.
     * 
     */
    virtual void CalculateDensitiesGCE() { CalculateDensities(); m_GCECalculated = true; }

    /**
     * \brief Calculates the total densities which include feeddown
     *        contributions.
     * 
     * Calculation is based on the primordial densities computed by CalculateDensities().
     * 
     */
    virtual void CalculateFeeddown();

    /**
     * \brief Computes the fluctuation observables.
     * 
     * Includes the matrix of 2nd order susceptibilities of
     * conserved charges, as well as particle number correlations
     * and fluctuations, both for primordial yields and after decays.
     * 
     */
    virtual void CalculateFluctuations();

    /**
     * \brief Computes the fluctuations and
     *        correlations of the primordial particle numbers.
     * 
     * More specifically, computes the susceptibility matrix
     * \f$ \frac{1}{VT^3} \, \langle \Delta N_i^* \Delta N_j^* \rangle \f$,
     * where \f$ N_i^* \f$ is the primordial yield of species i.
     * 
     */
    virtual void CalculateTwoParticleCorrelations();

    /**
     * \brief Computes particle number correlations
     *        and fluctuations for all final-state particles which
     *        are marked stable.
     * 
     * More specifically, computes the susceptibility matrix
     * \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_j \rangle \f$,
     * where \f$ N_i \f$ is the final yield of species i,
     * and indices i and j corresponds to stable particles only.
     * 
     * To incorporate probabilistic decays uses the formalism from paper 
     * V.V. Begun, M.I. Gorenstein, M. Hauer, V. Konchakovski, O.S. Zozulya,
     * Phys.Rev. C 74, 044903 (2006), 
     * [https://arxiv.org/pdf/nucl-th/0606036.pdf](https://arxiv.org/pdf/nucl-th/0606036.pdf)
     * 
     * 
     */
    virtual void CalculateTwoParticleFluctuationsDecays();

    /**
     * \brief Returns the computed primordial particle number (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_j \rangle \f$
     *        for particles with ids i and j. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double TwoParticleSusceptibilityPrimordial(int i, int j) const;

    /**
     * \brief Returns the computed primordial particle number (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_j \rangle \f$
     *        for particles with pdg codes id1 and id2. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double TwoParticleSusceptibilityPrimordialByPdg(long long id1, long long id2);

    /**
     * \brief Returns the computed (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_j \rangle \f$
     *        between primordial net-particle numbers for pdg codes id1 and id2. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double NetParticleSusceptibilityPrimordialByPdg(long long id1, long long id2);

    /**
     * \brief Returns the computed final particle number (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_j \rangle \f$
     *        for particles with ids i and j. CalculateFluctuations() must be called beforehand. Both particle species must be those marked stable.
     *
     */
    virtual double TwoParticleSusceptibilityFinal(int i, int j) const;

    /**
     * \brief Returns the computed final particle number (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_j \rangle \f$
     *        for particles with pdg codes id1 and id2. CalculateFluctuations() must be called beforehand. Both particle species must be those marked stable.
     *
     */
    virtual double TwoParticleSusceptibilityFinalByPdg(long long id1, long long id2);

    /**
     * \brief Returns the computed (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_j \rangle \f$
     *        between final net-particle numbers for pdg codes id1 and id2. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double NetParticleSusceptibilityFinalByPdg(long long id1, long long id2);

    /**
     * \brief Returns the computed primordial particle vs conserved charge (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_chg \rangle \f$
     *        for particle with id i and conserved charge chg. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double PrimordialParticleChargeSusceptibility(int i, ConservedCharge::Name chg) const;

    /**
     * \brief Returns the computed primordial particle vs conserved charge (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_chg \rangle \f$
     *        for particle with pdg code id1 and conserved charge chg. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double PrimordialParticleChargeSusceptibilityByPdg(long long id1, ConservedCharge::Name chg);

    /**
     * \brief Returns the computed primordial net particle vs conserved charge (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_chg \rangle \f$
     *        for particle with pdg code id1 and conserved charge chg. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double PrimordialNetParticleChargeSusceptibilityByPdg(long long id1, ConservedCharge::Name chg);


    /**
     * \brief Returns the computed final (after decays) particle vs conserved charge (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_chg \rangle \f$
     *        for particle with id i and conserved charge chg. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double FinalParticleChargeSusceptibility(int i, ConservedCharge::Name chg) const;

    /**
     * \brief Returns the computed final (after decays) particle vs conserved charge (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_chg \rangle \f$
     *        for particle with pdg code id1 and conserved charge chg. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double FinalParticleChargeSusceptibilityByPdg(long long id1, ConservedCharge::Name chg);

    /**
     * \brief Returns the computed final (after decays) net particle vs conserved charge (cross-)susceptibility \f$ \frac{1}{VT^3} \, \langle \Delta N_i \Delta N_chg \rangle \f$
     *        for particle with pdg code id1 and conserved charge chg. CalculateFluctuations() must be called beforehand.
     *
     */
    virtual double FinalNetParticleChargeSusceptibilityByPdg(long long id1, ConservedCharge::Name chg);

    /**
     * \brief Calculates the conserved charges susceptibility matrix.
     *        
     * Computes \f$ \chi_{lmnk}^{BSQC}~=~\frac{\partial^{l+m+n+k}p/T^4}{\partial(\mu_B/T)^l 
     * \,\partial(\mu_S/T)^m \,\partial(\mu_Q/T)^n\,\partial(\mu_C/T)^k} \f$
     * where i+j+k+l = 2.
     * 
     * The calculation results are accessible through Susc()
     * 
     */
    virtual void CalculateSusceptibilityMatrix();

    /**
     * \brief Calculates the susceptibility matrix of conserved charges proxies.
     * 
     * The following proxies are used:
     * final net proton number for baryon number,
     * net charge as is,
     * and net kaon number for net strangeness.
     * Charm not yet considered.
     * 
     * Decay feeddown contributions are in accordance with the stability flags.
     * 
     * The calculation results are accessible through ProxySusc()
     * 
     */
    virtual void CalculateProxySusceptibilityMatrix();

    /**
     * \brief Calculates the matrix of correlators between primordial (and also final) particle numbers and conserved charges.
     *
     * Correlators are normalized by VT^3.
     *
     * The calculation results are accessible through PrimordialParticleChargeCorrelation() and FinalParticleChargeCorrelation().
     *
     */
    virtual void CalculateParticleChargeCorrelationMatrix();

    /**
     * \brief Calculates fluctuations (diagonal susceptibilities) of an arbitrary "conserved" charge.
     * 
     * Each particle specie is assumed to carry a conserved charge with
     * a value provided by an input vector.
     * 
     * Restricted to the grand canonical ensemble.
     * 
     * \param chgs A vector with conserved charge values for all species.
     *             0-based indices of the vector must correspond to the
     *             0-based indices of the particle list TPS()
     * \param order Up to which order the susceptibilities are computed
     * \return std::vector<double> A vector with computed values of diagonal susceptibilities
     */
    virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);

    //virtual double GetParticlePrimordialDensity(unsigned int);
    //virtual double GetParticleTotalDensity(unsigned int);

    // Equation of state etc.

    /// System pressure (GeV fm\f$^{-3}\f$)
    double Pressure() { return CalculatePressure(); }

    /// Energy density (GeV fm\f$^{-3}\f$)
    double EnergyDensity() { return CalculateEnergyDensity(); }

    /// Entropy density (fm\f$^{-3}\f$)
    double EntropyDensity() { return CalculateEntropyDensity(); }

    /// Total number density of all particle (fm\f$^{-3}\f$)
    double HadronDensity() { return CalculateHadronDensity(); }

    /// Net baryon density (fm\f$^{-3}\f$)
    double BaryonDensity() { return CalculateBaryonDensity(); }

    /// Electric charge density (fm\f$^{-3}\f$)
    double ElectricChargeDensity() { return CalculateChargeDensity(); }

    /// Net strangeness density (fm\f$^{-3}\f$)
    double StrangenessDensity() { return CalculateStrangenessDensity(); }

    /// Net charm density (fm\f$^{-3}\f$)
    double CharmDensity() { return CalculateCharmDensity(); }

    /// Absolute baryon number density (baryons + antibaryons) (fm\f$^{-3}\f$)
    double AbsoluteBaryonDensity() { return CalculateAbsoluteBaryonDensity(); }

    /// Absolute electric charge density (Q+ + Q-) (fm\f$^{-3}\f$)
    double AbsoluteElectricChargeDensity() { return CalculateAbsoluteChargeDensity(); }

    /// Absolute strange quark content density (fm\f$^{-3}\f$)
    double AbsoluteStrangenessDensity() { return CalculateAbsoluteStrangenessDensity(); }

    /// Absolute charm quark content density (fm\f$^{-3}\f$)
    double AbsoluteCharmDensity() { return CalculateAbsoluteCharmDensity(); }

    //@{
    /// Implementation of the equation of state functions
    virtual double CalculatePressure() = 0;
    virtual double CalculateEnergyDensity() = 0;
    virtual double CalculateEntropyDensity() = 0;
    virtual double CalculateHadronDensity();
    virtual double CalculateBaryonDensity();
    virtual double CalculateChargeDensity();
    virtual double CalculateStrangenessDensity();
    virtual double CalculateCharmDensity();
    virtual double CalculateAbsoluteBaryonDensity();
    virtual double CalculateAbsoluteChargeDensity();
    virtual double CalculateAbsoluteStrangenessDensity();
    virtual double CalculateAbsoluteCharmDensity();
    virtual double CalculateAbsoluteStrangenessDensityModulo();
    virtual double CalculateAbsoluteCharmDensityModulo();
    //@}

    /// Computes the density of the auxiliary ArbitraryCharge()
    virtual double CalculateArbitraryChargeDensity();
    
    /// The fraction of entropy carried by baryons (Ideal GCE only)
    virtual double CalculateBaryonMatterEntropyDensity() { return 0.; }

    /// The fraction of entropy carried by mesons (Ideal GCE only)
    virtual double CalculateMesonMatterEntropyDensity() { return 0.; }

    /// Scaled variance of primordial particle number fluctuations for species i
    virtual double ParticleScaledVariance(int) { return 1.; }

    /// Skewness of primordial particle number fluctuations for species i
    virtual double ParticleSkewness(int) { return 1.; }

    /// Kurtosis of primordial particle number fluctuations for species i
    virtual double ParticleKurtosis(int) { return 1.; }

    /// Fraction of the total volume occupied by the finite-sizes particles 
    /// (Diagonal excluded volume model only)
    virtual double CalculateEigenvolumeFraction() { return 0.; }

    /// The scalar density of the particle species i
    virtual double ParticleScalarDensity(int i) = 0;

    virtual double GetMaxDiff() const { return m_MaxDiff; }

    /// Whether particle densities calculation
    /// in the CalculateDensities() method
    /// were successfull 
    virtual bool   IsLastSolutionOK() const { return m_LastCalculationSuccessFlag; }

    /**
     * \brief Same as GetDensity(int,Feeddown::Type)
     * 
     * \deprecated
     */
    double GetDensity(long long PDGID, int feeddown) { return GetDensity(PDGID, static_cast<Feeddown::Type>(feeddown)); }
    
    /**
     * \brief Particle number density of species
     *        with a specified PDG ID and feeddown
     * 
     * \param PDGID    Particle Data Group ID of the needed specie
     * \param feeddown Which decay feeddown contributions to take into account
     * \return Particle number density
     */
    double GetDensity(long long PDGID, Feeddown::Type feeddown);

    /**
     * \brief Particle number yield of species
     *        with a specified PDG ID and feeddown
     *
     * \param PDGID    Particle Data Group ID of the needed specie
     * \param feeddown Which decay feeddown contributions to take into account
     * \return Particle number yield
     */
    double GetYield(long long PDGID, Feeddown::Type feeddown) { return GetDensity(PDGID, feeddown) * Volume(); }

    std::vector<double> GetIdealGasDensities() const;

    /// A pointer to the ThermalParticleSystem object
    /// with the particle list.
    ThermalParticleSystem* TPS() { return m_TPS; }

    /// A vector with primordial particle number densities.
    /// Each entry corresponds to a density of single species,
    /// in accordance with the 0-based index of each specie.
    /// PdgToId() maps PDG ID to the 0-based index.
    const std::vector<double>& Densities()      const { return m_densities; }
    std::vector<double>& Densities() { return m_densities; }

    /// A vector with total particle number densities,
    /// which include the feeddown contribution in accordance
    /// with the stability flags.
    /// Each entry corresponds to a density of single species,
    /// in accordance with the 0-based index of each specie.
    /// PdgToId() maps PDG ID to the 0-based index.
    const std::vector<double>& TotalDensities() const { return m_densitiestotal; }

    /// A vector of vectors of particle number densities
    /// with the different feeddown contributions listed
    /// in Feeddown::Type
    const std::vector< std::vector<double> >& AllDensities()  const { return m_densitiesbyfeeddown; }

    /**
     * \brief Set the tag for this ThermalModelBase object.
     * 
     * \param tag Tag.
     */
    void SetTAG(const std::string & tag) { m_TAG = tag; }

    /// The tag of this ThermalModelBase object.
    const std::string& TAG() const { return m_TAG; }

    /// Reset all flags which correspond to a calculation status
    void ResetCalculatedFlags();

    /// Whether densities had been calculated with
    /// the CalculateDensities() method
    bool IsCalculated() const { return m_Calculated; }

    /// Whether fluctuations were calculated with
    /// the CalculateFluctuations() method
    bool IsFluctuationsCalculated() const { return m_FluctuationsCalculated; }

    /// Whether the grand-canonical ensemble particle
    /// number densities were calculated
    bool IsGCECalculated() const { return m_GCECalculated; }

    /**
     * \brief Scaled variance of primordial particle
     *        number fluctuations for particle species id
     * 
     * \param id  0-based index of particle species
     * \return    Scaled variance
     */
    double ScaledVariancePrimordial(int id) const { return (id >= 0 && id < static_cast<int>(m_wprim.size())) ? m_wprim[id] : 1.; }
    
    /**
     * \brief Scaled variance of final particle
     *        number fluctuations for particle species id
     * 
     * Decay feeddown is in accordance with the
     * stability flags.
     * 
     * \param id  0-based index of particle species
     * \return    Scaled variance
     */
    double ScaledVarianceTotal(int id) const { return (id >= 0 && id < static_cast<int>(m_wtot.size())) ? m_wtot[id] : 1.; }
    
    /**
     * \brief Normalized skewness of primordial particle
     *        number fluctuations for particle species id
     * 
     * \param id  0-based index of particle species
     * \return    Normalized skewness
     */
    double SkewnessPrimordial(int id) const { return (id >= 0 && id < static_cast<int>(m_skewprim.size())) ? m_skewprim[id] : 1.; }
    
    /**
     * \brief Normalized skewness of final particle
     *        number fluctuations for particle species id
     * 
     * Decay feeddown is in accordance with the
     * stability flags.
     * 
     * \param id  0-based index of particle species
     * \return    Normalized skewness
     */
    double SkewnessTotal(int id) const { return (id >= 0 && id < static_cast<int>(m_skewtot.size())) ? m_skewtot[id] : 1.; }
    
    /**
     * \brief Normalized excess kurtosis of primordial particle
     *        number fluctuations for particle species id
     * 
     * \param id  0-based index of particle species
     * \return    Normalized excess kurtosis
     */
    double KurtosisPrimordial(int id) const { return (id >= 0 && id < static_cast<int>(m_kurtprim.size())) ? m_kurtprim[id] : 1.; }
    
    /**
     * \brief Normalized excess kurtosis of final particle
     *        number fluctuations for particle species id
     * 
     * Decay feeddown is in accordance with the
     * stability flags.
     * 
     * \param id  0-based index of particle species
     * \return    Normalized excess kurtosis
     */
    double KurtosisTotal(int id) const { return (id >= 0 && id < static_cast<int>(m_kurttot.size())) ? m_kurttot[id] : 1.; }

    /**
     * \brief A density of a conserved charge (in fm^-3)
     *
     * \f$ \rho_{c_i} \f$
     *
     * \param i Conserved charge
     * \return  Conserved charge density \f$ \rho_{c_i} \f$
     */
    double ConservedChargeDensity(ConservedCharge::Name chg);

    /**
     * \brief A 2nd order susceptibility of conserved charges
     * 
     * \f$ \chi_{11}^{c_i c_j} \f$
     * 
     * \param i First conserved charge
     * \param j Second conserved charge
     * \return  Susceptibility \f$ \chi_{11}^{c_i c_j} \f$
     */
    double Susc(ConservedCharge::Name i, ConservedCharge::Name j) const { return m_Susc[i][j]; }

    /**
     * \brief A 2nd order susceptibility of conserved charges proxies
     * 
     * The following proxies are used:
     * final net proton number for baryon number,
     * net charge as is,
     * and net kaon number for net strangeness.
     * Charm not yet considered.
     * 
     * Decay feeddown contributions are in accordance with the stability flags.
     * 
     * \param i First proxy charge
     * \param j Second proxy charge
     * \return Proxy susceptibility 
     */
    double ProxySusc(ConservedCharge::Name i, ConservedCharge::Name j) const { return m_ProxySusc[i][j]; }

    /**
     * \brief Multiplicity of charged particles.
     * 
     * 
     * \param type    Type of charged particle multiplicity:
     *                0 -- multiplicity of all charged particles,
     *                \f$ \pm 1 \f$ -- multiplicity of positively (negatively) charged particles,
     *                \f$ \pm 2 \f$ -- total charge of all positively (negatively) charged particles
     * \return Multiplicity 
     */
    double ChargedMultiplicity(int type = 0);

    /**
     * \brief Scaled variance of charged particles.
     * 
     * 
     * \param type    Type of the scaled variance:
     *                0 -- scaled variance of all charged particles,
     *                \f$ \pm 1 \f$ -- scaled variance of positively (negatively) charged particles,
     *                \f$ \pm 2 \f$ -- scaled variance of total charge of all positively (negatively) charged particles
     * \return Scaled variance 
     */
    double ChargedScaledVariance(int type = 0);

    /**
     * \brief Multiplicity of charged particles including the feeddown
     *        contributions in accordance with the stability flags.
     * 
     * The meaning is the same as for ChargedMultiplicity()
     * 
     */
    double ChargedMultiplicityFinal(int type = 0);

    /**
     * \brief Scaled variance of charged particles including the feeddown
     *        contributions in accordance with the stability flags.
     * 
     * The meaning is the same as for ChargedScaledVariance()
     * 
     */
    double ChargedScaledVarianceFinal(int type = 0);

    /**
     * \brief The statistical ensemble of the current HRG model.
     * 
     * \return ThermalModelEnsemble Ensemble
     */
    ThermalModelEnsemble Ensemble() { return m_Ensemble; }


    /**
     * \brief Whether the given conserved charge is treated canonically
     */
    virtual bool IsConservedChargeCanonical(ConservedCharge::Name charge) const { return 0; }

    /**
     * \brief The interactions present in the current HRG model.
     * 
     * \return ThermalModelInteraction 
     */
    ThermalModelInteraction InteractionModel() { return m_InteractionModel; }

  protected:
    ThermalModelParameters m_Parameters;
    ThermalParticleSystem* m_TPS;

    bool   m_LastCalculationSuccessFlag;
    double m_MaxDiff;

    bool m_Calculated;
    bool m_FeeddownCalculated;
    bool m_FluctuationsCalculated;
    bool m_GCECalculated;
    bool m_UseWidth;
    bool m_NormBratio;
    bool m_QuantumStats;
    double m_QBgoal;
    double m_SBgoal;
    double m_Volume;

    bool m_ConstrainMuB;
    bool m_ConstrainMuQ;
    bool m_ConstrainMuS;
    bool m_ConstrainMuC;

    bool m_PCE;

    bool m_useOpenMP;

    std::vector<double> m_densities;
    std::vector<double> m_densitiestotal;
    //std::vector<double> m_densitiestotalweak;
    std::vector< std::vector<double> > m_densitiesbyfeeddown;
    std::vector<double> m_Chem;

    // Scaled variance
    std::vector<double> m_wprim;
    std::vector<double> m_wtot;

    // Skewness
    std::vector<double> m_skewprim;
    std::vector<double> m_skewtot;

    // Kurtosis
    std::vector<double> m_kurtprim;
    std::vector<double> m_kurttot;

    // 2nd order correlations of primordial and total numbers
    std::vector< std::vector<double> > m_PrimCorrel;
    std::vector< std::vector<double> > m_TotalCorrel;

    // Particle number-conserved charge correlators
    std::vector< std::vector<double> > m_PrimChargesCorrel;
    std::vector< std::vector<double> > m_FinalChargesCorrel;

    // Conserved charges susceptibility matrix
    std::vector< std::vector<double> > m_Susc;

    // Susceptibility matrix of net-p, net-Q, and net-K
    std::vector< std::vector<double> > m_ProxySusc;

    // Cumulants of arbitrary charge calculation
    //std::vector< std::vector<double> > m_chi;

    // Contains log of possible errors when checking the calculation
    std::string m_ValidityLog;

    double m_wnSum;

    std::string m_TAG;

    ThermalModelEnsemble m_Ensemble;
    ThermalModelInteraction m_InteractionModel;

    

    /// Shift in chemical potential of particle species id due to interactions
    virtual double MuShift(int /*id*/) const { return 0.; }

  private:
    void ResetChemicalPotentials();

    double GetDensity(long long PDGID, const std::vector<double> *dens);

    class BroydenEquationsChem : public BroydenEquations
    {
    public:
      BroydenEquationsChem(ThermalModelBase *model) : BroydenEquations(), m_THM(model) { m_N = 2; }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelBase *m_THM;
    };

    class BroydenJacobianChem : public BroydenJacobian
    {
    public:
      BroydenJacobianChem(ThermalModelBase *model) : BroydenJacobian(), m_THM(model) { }
      std::vector<double> Jacobian(const std::vector<double> &x);
    private:
      ThermalModelBase *m_THM;
    };

    class BroydenChem : public Broyden
    {
    public:
      BroydenChem(ThermalModelBase *THM, BroydenEquations *eqs = NULL, BroydenJacobian *jaco = NULL) : Broyden(eqs, jaco) { m_THM = THM; }
      ~BroydenChem(void) { }
      std::vector<double> Solve(const std::vector<double> &x0, BroydenSolutionCriterium *solcrit = NULL, int max_iterations = MAX_ITERS);
    private:
      ThermalModelBase *m_THM;
    };


    class BroydenEquationsChemTotals : public BroydenEquations
    {
    public:
      BroydenEquationsChemTotals(const std::vector<int> & vConstr, const std::vector<int> & vType, const std::vector<double> & vTotals, ThermalModelBase *model);// : BroydenEquations(), m_THM(model) { m_N = 3; }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      std::vector<int> m_Constr;
      std::vector<int> m_Type;
      std::vector<double> m_Totals;
      ThermalModelBase *m_THM;
    };

    class BroydenJacobianChemTotals : public BroydenJacobian
    {
    public:
      BroydenJacobianChemTotals(const std::vector<int> & vConstr, const std::vector<int> & vType, const std::vector<double> & vTotals, ThermalModelBase *model) : BroydenJacobian(), m_Constr(vConstr), m_Type(vType), m_Totals(vTotals), m_THM(model) { }
      std::vector<double> Jacobian(const std::vector<double> &x);
    private:
      std::vector<int> m_Constr;
      std::vector<int> m_Type;
      std::vector<double> m_Totals;
      ThermalModelBase *m_THM;
    };
  };

} // namespace thermalfist

#endif // THERMALMODELBASE_H
