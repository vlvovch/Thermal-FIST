/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2016-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELVDW_H
#define THERMALMODELVDW_H

#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  /**
   * \brief Class implementing the quantum van der Waals HRG model.
   * 
   * The model formulation can be found in
   * 
   * V. Vovchenko, M.I. Gorenstein, H. Stoecker,
   * Phys. Rev. Lett. **118**, 182301 (2017),
   * [http://arxiv.org/pdf/1609.03975.pdf](http://arxiv.org/pdf/1609.03975.pdf)
   * 
   * and in, in more detail, in
   * 
   * V. Vovchenko, A. Motornenko, P. Alba, M.I. Gorenstein, L.M. Satarov, H. Stoecker,
   * Phys. Rev. C **96**, 045202 (2017),
   * [http://arxiv.org/pdf/1707.09215.pdf](http://arxiv.org/pdf/1707.09215.pdf) 
   * 
   * The system of transcendental equations for 
   * the "shifted" chemical potentials of hadrons,
   * Eq. (15) in the latter reference,
   * is solved using the Broyden's method.
   * 
   */
  class ThermalModelVDW : public ThermalModelBase
  {
  public:
    /**
     * \brief Construct a new ThermalModelVDW object
     * 
     * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
     * \param params ThermalModelParameters object with current thermal parameters
     */
    ThermalModelVDW(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelVDW object
     * 
     */
    virtual ~ThermalModelVDW(void);
    
    /**
     * \brief Same as FillVirial() but uses the matrix of excluded-volume
     *        coefficients \f$ v_i \equiv b_{ii} \f$ as input instead of radii.
     * 
     * \param bij A vector with crossterms excluded-volume coefficients for all pairs of particle species
     */
    void FillVirialEV(const std::vector< std::vector<double> > & bij = std::vector< std::vector<double> >(0));

    /**
     * \brief Set the temperature derivative
     *        of the eigenvolume parameter \f$ \tilde{b}_{ij} \f$
     * 
     * \param i    0-based index of the first particle species
     * \param j    0-based index of the second particle species
     * \param dbdT \f$ d \tilde{b}_{ij} / dT \f$ in the units of fm\f$^3\f$ GeV\f$^{-1}\f$
     */
    void SetVirialdT(int i, int j, double dbdT) { if (i >= 0 && i < static_cast<int>(m_VirialdT.size()) && j >= 0 && j < static_cast<int>(m_VirialdT[i].size())) m_VirialdT[i][j] = dbdT; }
    
    /**
     * \brief Set the temperature derivative
     *        of the QvdW attraction parameter \f$ a_{ij} \f$
     * 
     * \param i    0-based index of the first particle species
     * \param j    0-based index of the second particle species
     * \param dadT \f$ d a_{ij} / dT \f$ in the units of fm\f$^3\f$
     */
    void SetAttractiondT(int i, int j, double dadT) { if (i >= 0 && i < static_cast<int>(m_AttrdT.size()) && j >= 0 && j < static_cast<int>(m_AttrdT[i].size()))     m_AttrdT[i][j] = dadT; }

    /**
     * \brief The temperature derivative
     *        of the eigenvolume parameter \f$ \tilde{b}_{ij} \f$
     * 
     * \param i    0-based index of the first particle species
     * \param j    0-based index of the second particle species
     * \return     dbdT \f$ d \tilde{b}_{ij} / dT \f$ in the units of fm\f$^3\f$ GeV\f$^{-1}\f$
     */
    double VirialCoefficientdT(int i, int j) const;
    
    /**
     * \brief The temperature derivative
     *        of the QvdW attraction parameter \f$ a_{ij} \f$
     * 
     * \param i    0-based index of the first particle species
     * \param j    0-based index of the second particle species
     * \return     \f$ d a_{ij} / dT \f$ in the units of fm\f$^3\f$
     */
    double AttractionCoefficientdT(int i, int j) const;

    /**
     * \brief Sets whether temperature depedence of
     *        QvdW parameters should be considered.
     * 
     * \param Tdep true -- considered, false -- not considered
     */
    void SetTemperatureDependentAB(bool Tdep) { m_TemperatureDependentAB = Tdep; }
   
    /**
     * \brief Whether temperature depedence of
     *        QvdW parameters is considered.
     * 
     * \return true  considered
     * \return false not considered
     */
    bool TemperatureDependentAB() const { return m_TemperatureDependentAB; }
    

    /**
     * \brief Whether to search for multiple solutions of the QvdW equations
     * by considering different initial guesses in the Broyden's method.
     * 
     * Multiple solutions in the QvdW model appear e.g. below the
     * critical temperature of the liquid-gas phase transition.
     * 
     * \param search Whether multiple solutions of the QvdW equations
     *               should be considered. False by default.
     */
    virtual void SetMultipleSolutionsMode(bool search) { m_SearchMultipleSolutions = search; }

    /**
     * \brief Whether to search for multiple solutions of the QvdW equations
     * by considering different initial guesses in the Broyden's method.
     *  
     * \return true  Multiple solutions considered
     * \return false Multiple solutions not considered
     */
    bool UseMultipleSolutionsMode() const { return m_SearchMultipleSolutions; }

    /// The shifted chemical potential of particle species i.
    double MuStar(int i) const { return m_MuStar[i]; }

    /// Returns vector of shifted chemical potentials,
    /// one element per each species
    std::vector<double> GetMuStar() const { return m_MuStar; }
    
    /// Set the vector of shifted chemical potentials
    void SetMuStar(const std::vector<double> & MuStar) { m_MuStar = MuStar; }

    // Override functions begin

    void FillChemicalPotentials();

    virtual void SetChemicalPotentials(const std::vector<double> & chem = std::vector<double>(0));

    virtual void FillVirial(const std::vector<double> & ri = std::vector<double>(0));

    virtual void FillAttraction(const std::vector< std::vector<double> > & aij = std::vector< std::vector<double> >(0));

    virtual void ReadInteractionParameters(const std::string &filename);

    virtual void WriteInteractionParameters(const std::string &filename);

    virtual void SetVirial(int i, int j, double b) { if (i >= 0 && i < static_cast<int>(m_Virial.size()) && j >= 0 && j < static_cast<int>(m_Virial[i].size())) m_Virial[i][j] = b; m_VDWComponentMapCalculated = false; }
    
    virtual void SetAttraction(int i, int j, double a) { if (i >= 0 && i < static_cast<int>(m_Attr.size()) && j >= 0 && j < static_cast<int>(m_Attr[i].size()))     m_Attr[i][j] = a; m_VDWComponentMapCalculated = false; }

    double VirialCoefficient(int i, int j) const;

    double AttractionCoefficient(int i, int j) const;

    virtual void ChangeTPS(ThermalParticleSystem *TPS);

    virtual void CalculatePrimordialDensities();

    virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);

    virtual std::vector< std::vector<double> >  CalculateFluctuations(int order);

    void CalculateTwoParticleCorrelations();

    void CalculateFluctuations();
    
    virtual double CalculatePressure();

    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    // Dummy
    virtual double CalculateBaryonMatterEntropyDensity() { return 0.; }

    virtual double CalculateMesonMatterEntropyDensity() { return 0.; }

    virtual double ParticleScalarDensity(int part);

    bool   IsLastSolutionOK() const { return m_LastBroydenSuccessFlag && m_LastCalculationSuccessFlag; }

    double DensityId(int index) { return m_DensitiesId[index]; }

    // Override functions end

    const std::vector< std::vector<int> >& VDWComponentIndices() const { return m_dMuStarIndices; }
    virtual double DeltaMu(int i) const { return MuShift(i); }
    const std::vector< std::vector<double> >& VirialMatrix() const { return m_Virial; }
    const std::vector< std::vector<double> >& AttractionMatrix() const { return m_Attr; }

  protected:
    /// Returns vector of particle densities
    /// for given values of shifted chemical potentials
    /// \param dmustar A vector of shifted chemical potentials
    /// \return std::vector<double> Vector of particle densities
    std::vector<double> ComputeNp(const std::vector<double>& dmustar);
    
    /// Same as ComputeNp(const std::vector<double>&) but using
    /// the vector of ideal gas densities as input instead of
    /// calculating it
    std::vector<double> ComputeNp(const std::vector<double>& dmustar, const std::vector<double>& ns);

    /// Partitions particles species into sets that have identical VDW parameters
    void CalculateVDWComponentsMap();

    /**
     * \brief Uses the Broyden method with a provided initial guess
     *        to determine the shifted chemical potentials
     *        by solving the transcendental equations with the
     *        Broyden's method.
     * 
     * \param muStarInit Initial guess for the shifted chemical potentials
     * \return std::vector<double> The solved shifted chemical potentials
     */
    virtual std::vector<double> SearchSingleSolution(const std::vector<double> & muStarInit);

    /**
     * \brief Uses the Broyden method with different initial guesses
     *        to look for different possible solutions
     *        of the transcendental equations for
     *        shifted chemical potentials
     * 
     * Looks for the solution with the largest pressure.
     * 
     * \param iters Number of different initial guesses to try
     * \return std::vector<double> The solution with the largest pressure among those which were found
     */
    std::vector<double> SearchMultipleSolutions(int iters = 300);

    /// Solve the transcedental equations for the
    /// shifted chemical potentials
    void SolveEquations();

    /**
     * \brief The shift in the chemical potential
     *        of particle species i due to the
     *        QvdW interactions.
     * 
     * \param i 0-based particle specie index
     * \return  The shift in the chemical potential
     */
    virtual double MuShift(int id) const;

    /// Vector of ideal gas densities with shifted chemical potentials
    std::vector<double> m_DensitiesId;

    /// Vector of scalar densities. Not used.
    std::vector<double> m_scaldens;

    /// \copydoc thermalfist::ThermalModelEVCrossterms::m_Virial
    std::vector< std::vector<double> > m_Virial;

    /// Matrix of the attractive QvdW coefficients \f$ a_{ij} \f$
    std::vector< std::vector<double> > m_Attr;

    /// Matrix of the temperature derivatives
    /// of the virial (excluded-volume) coefficients \f$ d \tilde{b}_{ij} / dT \f$
    std::vector< std::vector<double> > m_VirialdT;

    /// Matrix of the temperature derivatives
    /// of the attractive QvdW coefficients \f$ d a_{ij} / dT \f$
    std::vector< std::vector<double> > m_AttrdT;

    /// Whether temperature depedence of
    /// QvdW parameters is considered.
    bool   m_TemperatureDependentAB;

    /// Whether multiple solutions are considered
    bool   m_SearchMultipleSolutions;

    /// Whether Broyden's method was successfull
    bool   m_LastBroydenSuccessFlag;

    /// Whether the mapping to components with the same VDW parameters has been calculated
    bool   m_VDWComponentMapCalculated;

    /// Vector of the shifted chemical potentials
    std::vector<double> m_MuStar;

    std::vector<int> m_MapTodMuStar;

    std::vector<int> m_MapFromdMuStar;

    std::vector< std::vector<int> > m_dMuStarIndices;

  private:
    std::vector< std::vector<double> > m_chi;

    std::vector<double> m_chiarb;
    
    
    virtual void CalculatePrimordialDensitiesOld();

    virtual void CalculatePrimordialDensitiesNew();

    class BroydenEquationsVDW : public BroydenEquations
    {
    public:
      BroydenEquationsVDW(ThermalModelVDW *model) : BroydenEquations(), m_THM(model) { m_N = model->m_MapFromdMuStar.size(); }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelVDW *m_THM;
    };

    class BroydenJacobianVDW : public BroydenJacobian
    {
    public:
      BroydenJacobianVDW(ThermalModelVDW *model) : BroydenJacobian(), m_THM(model) { }
      std::vector<double> Jacobian(const std::vector<double> &x);
    private:
      ThermalModelVDW *m_THM;
    };

    class BroydenSolutionCriteriumVDW : public Broyden::BroydenSolutionCriterium
    {
    public:
      BroydenSolutionCriteriumVDW(ThermalModelVDW *model, double relative_error = Broyden::TOL) : Broyden::BroydenSolutionCriterium(relative_error), m_THM(model) { }
      virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
    protected:
      ThermalModelVDW *m_THM;
    };
  };

  /// For backward compatibility
  typedef ThermalModelVDW ThermalModelVDWFull;

} // namespace thermalfist

#endif

