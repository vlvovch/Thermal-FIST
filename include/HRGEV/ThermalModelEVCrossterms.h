#ifndef ThermalModelEVMultiEigenEIGEN_H
#define ThermalModelEVMultiEigenEIGEN_H

#include "HRGBase/ThermalModelBase.h"


class ThermalModelEVCrossterms : public ThermalModelBase
{
public:
	ThermalModelEVCrossterms(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters(), double RHad_ = 0., int mode = 0);

	virtual ~ThermalModelEVCrossterms(void);

	void FillVirial(const std::vector<double> & ri = std::vector<double>(0));

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
	void SolvePressure(bool resetPartials = true);	// Using Broyden's method
	void SolvePressureIter();	// Using iteration method
	virtual void CalculateDensities();
	virtual void CalculateDensitiesNoReset();
	virtual void CalculateDensitiesIter();
	void CalculateTwoParticleCorrelations();
	void CalculateFluctuations();
	virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);
	double DensityId(int ind);
	double Pressure(int ind);
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

private:
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

#endif

