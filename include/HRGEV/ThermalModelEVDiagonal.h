/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
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

  class ThermalModelEVDiagonal : public ThermalModelBase
  {
  public:
    ThermalModelEVDiagonal(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters(), double RHad_ = 0., int mode = 0);

    virtual ~ThermalModelEVDiagonal(void);

    void SetRadius(double rad);
    void FillVirial(const std::vector<double> & ri = std::vector<double>(0));
    void FillVirialEV(const std::vector<double> & vi = std::vector<double>(0));

    virtual void ReadInteractionParameters(const std::string &filename);
    virtual void WriteInteractionParameters(const std::string &filename);
    double ExcludedVolume(int i) const;// { return m_v[i]; }
    virtual double CalculateEigenvolumeFraction();
    void SetRadius(int i, double rad);

    double VirialCoefficient(int i, int j) const;
    void SetVirial(int i, int j, double b);

    virtual void SetParameters(const ThermalModelParameters& params);

    virtual void ChangeTPS(ThermalParticleSystem *TPS_);

    void SolvePressure();

    virtual void CalculateDensities();

    virtual void CalculateTwoParticleCorrelations();

    virtual void CalculateFluctuations();

    virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);

    double DensityId(int ind);

    double PressureId(int ind);

    double ScaledVarianceId(int ind);

    double Pressure(double P);


    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    // Dummy
    virtual double CalculateBaryonMatterEntropyDensity();

    virtual double CalculateMesonMatterEntropyDensity();

    virtual double CalculatePressure();

    //virtual double CalculateShearViscosity();

    // TODO Properly for multi-component
    virtual double CalculateHadronScaledVariance();

    // TODO Properly for multi-component
    virtual double CalculateParticleScaledVariance(int part);

    // TODO properly with excluded volume
    virtual double CalculateParticleSkewness(int part);

    // TODO properly with excluded volume
    virtual double CalculateParticleKurtosis(int part);

    virtual double CalculateBaryonScaledVariance(bool susc = false);
    virtual double CalculateChargeScaledVariance(bool susc = false);
    virtual double CalculateStrangenessScaledVariance(bool susc = false);

    virtual double ParticleScalarDensity(int part);

    double CommonSuppressionFactor();

  protected:
    // TODO: test
    virtual double MuShift(int id);

  private:
    std::vector<double> m_densitiesid;
    std::vector<double> m_densitiesidnoshift;
    std::vector<double> m_v;                       /**< Vector of eigenvolumes of all hadrons */
    double m_Suppression;
    double m_Pressure;
    double m_Densityid;
    double m_TotalDensity;
    EVSolution m_sol;
  };

} // namespace thermalfist

#endif

