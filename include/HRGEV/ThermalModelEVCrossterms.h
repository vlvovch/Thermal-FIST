/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELEVCROSSTERMS_H
#define THERMALMODELEVCROSSTERMS_H

#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  class ThermalModelEVCrossterms : public ThermalModelBase
  {
  public:
    ThermalModelEVCrossterms(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters(), double RHad_ = 0., int mode = 0);

    virtual ~ThermalModelEVCrossterms(void);

    virtual void FillVirial(const std::vector<double> & ri = std::vector<double>(0));

    virtual void SetParameters(double T, double muB, double muS, double muQ, double gammaS, double V, double R);

    virtual void ReadInteractionParameters(const std::string &filename);
    virtual void WriteInteractionParameters(const std::string &filename);
    void SetRadius(double rad);

    void DisableBBarRepulsion();

    double VirialCoefficient(int i, int j) const;
    void SetVirial(int i, int j, double b);

    virtual void SetParameters(const ThermalModelParameters& params);

    virtual void ChangeTPS(ThermalParticleSystem *TPS_);

    virtual void SolveDiagonal();
    virtual void SolvePressure(bool resetPartials = true);  // Using Broyden's method
    void SolvePressureIter();  // Using iteration method
    virtual void CalculateDensities();
    virtual void CalculateDensitiesNoReset();
    virtual void CalculateDensitiesIter();
    void CalculateTwoParticleCorrelations();
    void CalculateFluctuations();
    virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);
    virtual double DensityId(int ind);
    virtual double Pressure(int ind);
    virtual double DensityId(int ind, const std::vector<double>& pstars);
    virtual double Pressure(int ind, const std::vector<double>& pstars);
    double ScaledVarianceId(int ind);
    double PressureDiagonal(int ind, double P);
    double PressureDiagonalTotal(double P);


    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    // Dummy
    virtual double CalculateBaryonMatterEntropyDensity() { return 0.; }
    virtual double CalculateMesonMatterEntropyDensity() { return 0.; }

    virtual double CalculatePressure();

    virtual double CalculateHadronScaledVariance() { return 1.; } // TODO properly

    virtual double CalculateParticleScaledVariance(int part) { return 1.; }// { return m_wprim[part]; }

    // TODO properly with excluded volume
    virtual double CalculateParticleSkewness(int part) { return m_skewprim[part]; }

    // TODO properly with excluded volume
    virtual double CalculateParticleKurtosis(int part) { return m_kurtprim[part]; }

    virtual double CalculateBaryonScaledVariance(bool susc = false);
    virtual double CalculateChargeScaledVariance(bool susc = false);
    virtual double CalculateStrangenessScaledVariance(bool susc = false);
    virtual double ParticleScalarDensity(int part) { return 0.; }

  protected:
    // TODO: test
    virtual double MuShift(int id);

    //private:
    std::vector<double> m_densitiesid;
    std::vector<double> m_Ps;
    std::vector< std::vector<double> > m_Virial;
    double m_Suppression;
    double m_Pressure;
    double m_Densityid;
    double m_TotalDensity;
    double m_TotalEntropyDensity;
    double m_RHad;
    int m_Mode;

  private:
    class BroydenEquationsCRS : public BroydenEquations
    {
    public:
      BroydenEquationsCRS(ThermalModelEVCrossterms *model) : BroydenEquations(), m_THM(model) { m_N = model->Densities().size(); }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelEVCrossterms *m_THM;
    };

    class BroydenEquationsCRSDEV : public BroydenEquations
    {
    public:
      BroydenEquationsCRSDEV(ThermalModelEVCrossterms *model) : BroydenEquations(), m_THM(model) { m_N = 1; }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelEVCrossterms *m_THM;
    };

    class BroydenJacobianCRS : public BroydenJacobian
    {
    public:
      BroydenJacobianCRS(ThermalModelEVCrossterms *model) : BroydenJacobian(), m_THM(model) { }
      Eigen::MatrixXd Jacobian(const std::vector<double> &x);
    private:
      ThermalModelEVCrossterms *m_THM;
    };

    class BroydenSolutionCriteriumCRS : public Broyden::BroydenSolutionCriterium
    {
    public:
      BroydenSolutionCriteriumCRS(ThermalModelEVCrossterms *model, double relative_error = Broyden::TOL) : Broyden::BroydenSolutionCriterium(relative_error), m_THM(model) { }
      virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
    protected:
      ThermalModelEVCrossterms *m_THM;
    };
  };

} // namespace thermalfist

#endif

