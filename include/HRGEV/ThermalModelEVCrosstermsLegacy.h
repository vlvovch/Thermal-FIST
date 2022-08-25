/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELEVCROSSTERMSLEGACY_H
#define THERMALMODELEVCROSSTERMSLEGACY_H

#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  /**
   * \brief Class implementing the crossterms
   *        excluded-volume model.
   * 
   * The model formulation can be found in
   * 
   * M.I. Gorenstein, A.P. Kostyuk, Y.D. Krivenko,
   * J.Phys. G **25**, L75 (1999),
   * [http://arxiv.org/pdf/nucl-th/9906068.pdf](http://arxiv.org/pdf/nucl-th/9906068.pdf)
   * 
   * and in
   * 
   * V. Vovchenko, H. Stoecker,
   * Phys. Rev. C **95**, 044904 (2017),
   * [http://arxiv.org/pdf/1606.06218.pdf](http://arxiv.org/pdf/1606.06218.pdf) 
   * 
   * The system of transcendental equations for 
   * "partial pressures" of hadrons
   * is solved using the Broyden's method.
   * 
   * \deprecated Since version 1.4 ThermalModelEVCrossterms() is implemented as partial case of ThermalModelVDW()
   *             The implementation in this file is not a legacy implementation
   */
  class ThermalModelEVCrosstermsLegacy : public ThermalModelBase
  {
  public:
    /**
     * \brief Construct a new ThermalModelEVCrossterms object
     * 
     * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
     * \param params ThermalModelParameters object with current thermal parameters
     */
    ThermalModelEVCrosstermsLegacy(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelEVCrossterms object
     * 
     */
    virtual ~ThermalModelEVCrosstermsLegacy(void);

    // Override functions begin

    virtual void FillVirial(const std::vector<double> & ri = std::vector<double>(0));

    virtual void ReadInteractionParameters(const std::string &filename);

    virtual void WriteInteractionParameters(const std::string &filename);

    void SetRadius(double rad);

    double VirialCoefficient(int i, int j) const;

    void SetVirial(int i, int j, double b);

    virtual void ChangeTPS(ThermalParticleSystem *TPS);
    
    virtual void CalculatePrimordialDensities();

    void CalculateTwoParticleCorrelations();

    void CalculateFluctuations();

    virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);

    virtual double CalculatePressure();

    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    // Dummy
    virtual double CalculateBaryonMatterEntropyDensity() { return 0.; }
    virtual double CalculateMesonMatterEntropyDensity() { return 0.; }

    // TODO properly with excluded volume
    virtual double ParticleScaledVariance(int /*part*/) { return 1.; }

    // TODO properly with excluded volume
    virtual double ParticleSkewness(int part) { return m_skewprim[part]; }

    // TODO properly with excluded volume
    virtual double ParticleKurtosis(int part) { return m_kurtprim[part]; }

    // TODO properly with excluded volume
    virtual double ParticleScalarDensity(int /*part*/) { return 0.; }

    // Override functions end


    const std::vector< std::vector<int> >& EVComponentIndices() const { return m_EVComponentIndices; }
    virtual double DeltaMu(int i) const { return MuShift(i); }
    const std::vector< std::vector<double> >& VirialMatrix() const { return m_Virial; }

  protected:
    /**
     * \brief Solves the system of transcdental equations 
     *        for the pressure using the Broyden's method.
     * 
     * Solves \f$ p_i(T,\mu) = p_i^{\rm id} (T, \mu_i - \sum_j \tilde{b_{ij}} p_j). \f$
     * 
     */
    virtual void SolvePressure(bool resetPartials = true);  // Using Broyden's method

    void SolvePressureIter();  // Using iteration method
    
    virtual void CalculatePrimordialDensitiesNoReset();

    virtual void CalculatePrimordialDensitiesIter();

    /**
     * \brief Solves the transcendental equation of 
     *        the corresponding diagonal EV model.
     * 
     * The diagonal EV model here has \f$v_i = \tilde{b}_{ii}\f$.
     * 
     * The partial pressures of the diagonal model
     * will be recorded into m_Ps.
     * 
     */
    virtual void SolveDiagonal();

    /**
     * \brief The "partial pressure" of hadron species i
     *        for the given total pressure in the diagonal model.
     * 
     * \param i 0-based particle species index.
     * \param P Input pressure (GeV fm\f$^{-3}\f$)
     * \return  Computed partial pressure (GeV fm\f$^{-3}\f$)
     */
    double PartialPressureDiagonal(int i, double P);

    /**
     * \brief The total pressure 
     *        for the given input pressure in the diagonal model.
     * 
     * \param   Input pressure (GeV fm\f$^{-3}\f$)
     * \return  Computed pressure (GeV fm\f$^{-3}\f$)
     */
    double PressureDiagonalTotal(double P);

    /**
     * \brief Calculate the ideal gas density of
     *        particle species i for the given values
     *        of partial pressures.
     * 
     * \param i        0-based particle specie index
     * \param Pressure Input vector of partial pressures (GeV fm\f$^{-3}\f$)
     * \return         Ideal gas density (fm\f$^{-3}\f$)
     */
    virtual double DensityId(int i, const std::vector<double>& pstars = std::vector<double>());

    /**
     * \brief Calculate the ideal gas pressure of
     *        particle species i for the given values
     *        of partial pressures.
     * 
     * \param i        0-based particle specie index
     * \param Pressure Input vector of partial pressures (GeV fm\f$^{-3}\f$)
     * \return         Ideal gas pressure (fm\f$^{-3}\f$)
     */
    virtual double Pressure(int i, const std::vector<double>& pstars = std::vector<double>());

    /**
     * \brief Calculate the ideal gas scaled variance of
     *        particle species i number fluctuations 
     *        for the given values
     *        of partial pressures.
     * 
     * \param i        0-based particle specie index
     * \param Pressure Input vector of partial pressures (GeV fm\f$^{-3}\f$)
     * \return         Ideal gas scaled variance
     */
    double ScaledVarianceId(int ind);
    
    /**
     * \brief The shift in the chemical potential
     *        of particle species i due to the
     *        excluded volume interactions.
     * 
     * Equal to \f$- \sum_j \tilde{b_{ij}} p_j\f$
     * 
     * \param i 0-based particle specie index
     * \return  The shift in the chemical potential
     */
    virtual double MuShift(int i) const;

    std::vector<double> m_densitiesid;            /**< Vector of ideal gas densities with shifted chemical potentials */
    std::vector<double> m_Ps;                     /**< Vector of (solved) partial pressures */
    std::vector< std::vector<double> > m_Virial;  /**< Matrix of virial (excluded-volume) coefficients \f$ \tilde{b}_{ij} \f$ */
    double m_Pressure;                            /**< The (solved) total pressure */
    double m_TotalEntropyDensity;                 /**< The (solved) entropy pressure */


    std::vector<int> m_MapToEVComponent;
    std::vector<int> m_MapFromEVComponent;
    std::vector< std::vector<int> > m_EVComponentIndices;

  private:
    class BroydenEquationsCRS : public BroydenEquations
    {
    public:
      BroydenEquationsCRS(ThermalModelEVCrosstermsLegacy*model) : BroydenEquations(), m_THM(model) { m_N = model->Densities().size(); }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelEVCrosstermsLegacy *m_THM;
    };

    class BroydenEquationsCRSDEV : public BroydenEquations
    {
    public:
      BroydenEquationsCRSDEV(ThermalModelEVCrosstermsLegacy*model) : BroydenEquations(), m_THM(model) { m_N = 1; }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelEVCrosstermsLegacy *m_THM;
    };

    class BroydenJacobianCRS : public BroydenJacobian
    {
    public:
      BroydenJacobianCRS(ThermalModelEVCrosstermsLegacy*model) : BroydenJacobian(), m_THM(model) { }
      std::vector<double> Jacobian(const std::vector<double> &x);
    private:
      ThermalModelEVCrosstermsLegacy *m_THM;
    };

    class BroydenSolutionCriteriumCRS : public Broyden::BroydenSolutionCriterium
    {
    public:
      BroydenSolutionCriteriumCRS(ThermalModelEVCrosstermsLegacy*model, double relative_error = Broyden::TOL) : Broyden::BroydenSolutionCriterium(relative_error), m_THM(model) { }
      virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
    protected:
      ThermalModelEVCrosstermsLegacy *m_THM;
    };
  };

} // namespace thermalfist

#endif

