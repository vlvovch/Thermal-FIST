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

    void SolveDiagonal();
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
  };

} // namespace thermalfist

#endif

