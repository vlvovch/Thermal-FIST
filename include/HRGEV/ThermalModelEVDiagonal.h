/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELEVDIAGONAL_H
#define THERMALMODELEVDIAGONAL_H

#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  struct EVSolution {
    double T, mu;
    double n, P, en, s;
    double mutil;
    double nid, eid, sid, pid;
    double fav, f2av, f3av, wav, wfav, w2fav, w2av;
    double scv, scvid;
    double dnmuid, dndTid, demuid, deTid, dsmuid, dsTid, dmunT, dPnT, dmuTn, dPTn, denT, deTn, dsignT, dsigTn, dTnsig;
    double dwiddmu, dfdmu, df2dmu, d2widdmu2;
    double dndmu, dndT, dendmu, dendT;
    double skew, kurtosis;
    double Cs;
    double PT4, eT4, sT3;
    double deltaEN, sigmaEN;
  };

  /**
   * \brief Class implementing the diagonal
   *        excluded-volume model.
   * 
   * The model formulation can be found in
   * G.D. Yen, M.I. Gorenstein, W. Greiner, S.-N. Yang,
   * Phys. Rev. C **56** 2210, (1997),
   * [http://arxiv.org/pdf/nucl-th/9711062.pdf](http://arxiv.org/pdf/nucl-th/9711062.pdf)
   * 
   * The transcendental equation for the pressure
   * is solved
   * using the Broyden's method.
   * 
   */
  class ThermalModelEVDiagonal : public ThermalModelBase
  {
  public:
    /**
     * \brief Construct a new ThermalModelEVDiagonal object
     * 
     * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
     * \param params ThermalModelParameters object with current thermal parameters
     */
    ThermalModelEVDiagonal(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelEVDiagonal object
     * 
     */
    virtual ~ThermalModelEVDiagonal(void);

    /**
     * \brief Same as FillVirial() but uses the diagonal excluded-volume
     *        coefficients \f$ v_i \equiv b_{ii} \f$ as input instead of radii.
     * 
     * \param vi A vector with diagonal excluded-volume coefficients for all species.
     *           0-based indices of the vector must corresponds to the
     *           0-based indices of the particle list TPS()
     */
    void FillVirialEV(const std::vector<double> & vi = std::vector<double>(0));
    
    /**
     * \brief The excluded volume parameter of
     *        particle species i.
     * 
     * \param i 0-based particle species index.
     * \return  The excluded volume parameter.
     */
    double ExcludedVolume(int i) const;
    
    /**
     * \brief The density suppression factor
     *        \f$ (1 - \sum_i v_i n_i) \f$, 
     *        common for all species.
     * 
     * \return The density suppression factor 
     */
    double CommonSuppressionFactor();

    
    //@{
    /// Differential treatment of pair interactions 
    /// is not applicable for the diagonal
    /// excluded volume model
    virtual void DisableMesonMesonVirial() {}
    virtual void DisableMesonMesonAttraction() {}
    virtual void DisableMesonBaryonVirial() {}
    virtual void DisableMesonBaryonAttraction() {}
    virtual void DisableBaryonBaryonVirial() {}
    virtual void DisableBaryonBaryonAttraction() {}
    virtual void DisableBaryonAntiBaryonVirial() {}
    virtual void DisableBaryonAntiBaryonAttraction() {}
    //@}

    // Override functions begin

    void SetRadius(double rad);
    
    void SetRadius(int i, double rad);

    /**
     * \brief Fills the vector of particle
     *        eigenvolume parameters.
     * 
     * \param ri A vector of radii parameters for all particles (in fm)
     */
    void FillVirial(const std::vector<double> & ri = std::vector<double>(0));

    /**
     * \brief Read the set of eigenvolume parameters for
     *        all particles from an external file.
     * 
     * \param filename Input file name
     */
    virtual void ReadInteractionParameters(const std::string &filename);

    /**
     * \brief Write the set of eigenvolume parameters for
     *        all particles to an external file.
     * 
     * One particle specie per line.
     * 
     * \param filename Output file name
     */
    virtual void WriteInteractionParameters(const std::string &filename);

    virtual double CalculateEigenvolumeFraction();

    /**
     * \brief Return the eigenvolume parameter
     *        of particle species i.
     * 
     * \param i 0-based particle specie index
     * \param j Irrelevant
     * \return  The eigenvolume parameter (fm\f$^{-3}\f$)
     */
    double VirialCoefficient(int i, int j) const;

    /**
     * \brief Sets the eigenvolume parameter
     *        of particle species i.
     * 
     * \param i 0-based particle specie index
     * \param j Here must be equal to i
     * \param b The eigenvolume parameter (fm\f$^{-3}\f$)
     */
    void SetVirial(int i, int j, double b);

    virtual void ChangeTPS(ThermalParticleSystem *TPS);

    virtual void CalculatePrimordialDensities();
    
    virtual void CalculateFluctuations();

    virtual void CalculateTwoParticleCorrelations();

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
    virtual double ParticleSkewness(int /*part*/) { return 1.; }

    // TODO properly with excluded volume
    virtual double ParticleKurtosis(int /*part*/) { return 1.; }

    virtual double ParticleScalarDensity(int part);

    // Override functions end

  protected:
    /**
     * \brief Solves the transcdental equation for the pressure.
     * 
     * Solves \f$ p(T,\mu) = \sum_i p_i^{\rm id} (T, \mu_i - v_i p). \f$
     * 
     */
    virtual void SolvePressure();

    /**
     * \brief Computes the l.h.s. of the transcdental equation for the pressure.
     * 
     * \param P Input pressure (GeV fm\f$^{-3}\f$)
     * \return  Computed pressure (GeV fm\f$^{-3}\f$)
     */
    double Pressure(double P);

    /**
     * \brief Calculate the ideal gas density of
     *        particle species i for the given pressure value.
     * 
     * \param i        0-based particle specie index
     * \param Pressure Input pressure (GeV fm\f$^{-3}\f$)
     * \return         Ideal gas density (fm\f$^{-3}\f$)
     */
    double DensityId(int i, double Pressure);

    /**
     * \brief Calculate the ideal gas pressure of
     *        particle species i for the given pressure value.
     * 
     * \param i        0-based particle specie index
     * \param Pressure Input pressure (GeV fm\f$^{-3}\f$)
     * \return         Ideal gas pressure (fm\f$^{-3}\f$)
     */
    double PressureId(int i, double Pressure);

    /**
     * \brief Calculate the ideal gas scaled variance of
     *        particle species i number fluctuations 
     *        for the given pressure value.
     * 
     * \param i        0-based particle specie index
     * \param Pressure Input pressure (GeV fm\f$^{-3}\f$)
     * \return         Ideal gas scaled variance
     */
    double ScaledVarianceId(int i, double Pressure);

    /**
     * \brief The shift in the chemical potential
     *        of particle species i due to the
     *        excluded volume interactions.
     * 
     * Equal to \f$-v_i p\f$
     * 
     * \param i 0-based particle specie index
     * \return  The shift in the chemical potential
     */
    virtual double MuShift(int i) const;

    std::vector<double> m_densitiesid;             /**< Vector of ideal gas densities with shifted chemical potentials */
    std::vector<double> m_densitiesidnoshift;      /**< Vector of ideal gas densities without shifted chemical potentials */
    std::vector<double> m_v;                       /**< Vector of eigenvolumes of all hadrons */
    double m_Suppression;                          /**< The common density suppression factor */
    double m_Pressure;                             /**< The solved pressure */
    double m_Densityid;
    double m_TotalDensity;
    EVSolution m_sol;                              /**< Axuiliary */

  private:
    class BroydenEquationsDEV : public BroydenEquations
    {
    public:
      BroydenEquationsDEV(ThermalModelEVDiagonal *model) : BroydenEquations(), m_THM(model) { m_N = 1; m_mnc = 1.; }
      std::vector<double> Equations(const std::vector<double> &x);
      void SetMnc(double mnc) { m_mnc = mnc; }
    private:
      ThermalModelEVDiagonal *m_THM;
      double m_mnc;
    };

    class BroydenJacobianDEV : public BroydenJacobian
    {
    public:
      BroydenJacobianDEV(ThermalModelEVDiagonal *model) : BroydenJacobian(), m_THM(model) { m_mnc = 1.; }
      std::vector<double> Jacobian(const std::vector<double> &x);
      void SetMnc(double mnc) { m_mnc = mnc; }
    private:
      ThermalModelEVDiagonal *m_THM;
      double m_mnc;
    };

    class BroydenSolutionCriteriumDEV : public Broyden::BroydenSolutionCriterium
    {
    public:
      BroydenSolutionCriteriumDEV(ThermalModelEVDiagonal *model, double relative_error = Broyden::TOL) : Broyden::BroydenSolutionCriterium(relative_error), m_THM(model) { }
      virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
    protected:
      ThermalModelEVDiagonal *m_THM;
    };

    class BroydenEquationsDEVOrig : public BroydenEquations
    {
    public:
      BroydenEquationsDEVOrig(ThermalModelEVDiagonal *model) : BroydenEquations(), m_THM(model) { m_N = 1; }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelEVDiagonal *m_THM;
    };

    class BroydenJacobianDEVOrig : public BroydenJacobian
    {
    public:
      BroydenJacobianDEVOrig(ThermalModelEVDiagonal *model) : BroydenJacobian(), m_THM(model) {}
      std::vector<double> Jacobian(const std::vector<double> &x);
    private:
      ThermalModelEVDiagonal *m_THM;
    };
  };

} // namespace thermalfist

#endif

