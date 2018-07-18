#ifndef THERMALMODELCANONICAL_H
#define THERMALMODELCANONICAL_H

#include <map>


#include "HRGBase/ThermalModelIdeal.h"

struct QuantumNumbers
{
	int B, Q, S, C;
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

class ThermalModelCanonical :
	public ThermalModelBase
{
public:
	ThermalModelCanonical(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters());

	virtual ~ThermalModelCanonical(void);

	void ChangeTPS(ThermalParticleSystem *TPS_);

	virtual void CalculateQuantumNumbersRange(bool doubleRange = false);

	void SetParameters(double T, double gammaS, double V, int B, int Q, int S = 0, int C = 0);

	void SetParameters(const ThermalModelParameters& params);


	virtual void SetStatistics(bool stats);

	virtual void FixParameters();
	virtual void FixParametersNoReset();
	//virtual void FixParameters(double) { };
	//virtual void FixParametersNoReset() { };

	virtual void CalculateDensities();
	virtual void ValidateCalculation();

	virtual void CalculatePartitionFunctions(double Vc = -1.);
	//virtual void CalculatePartitionFunctionsBoseOnly();

	// TODO check the validity
	virtual double CalculateParticleScaledVariance(int part);
	virtual void CalculateTwoParticleCorrelations();

	// TODO properly for higher moments
	virtual void CalculateFluctuations();


	virtual double CalculateEnergyDensity();

	virtual double CalculatePressure();

	virtual double CalculateEntropyDensity();

	virtual double GetGCEDensity(int i) const;

	virtual double ParticleScalarDensity(int part) { return 0.; }

	virtual bool IsParticleCanonical(const ThermalParticle &part);

	virtual void ConserveBaryonCharge(bool conserve = true) { m_BCE = static_cast<int>(conserve); }
	virtual void ConserveElectricCharge(bool conserve = true) { m_QCE = static_cast<int>(conserve); }
	virtual void ConserveStrangeness(bool conserve = true) { m_SCE = static_cast<int>(conserve); }
	virtual void ConserveCharm(bool conserve = true) { m_CCE = static_cast<int>(conserve); }

private:
	void PrepareModelGCE();  /**< Creates the ThermalModelIdeal copy */

	void CleanModelGCE();		/**< Cleares the ThermalModelIdeal copy */

protected:

	std::map<QuantumNumbers, int> m_QNMap;
	std::vector<QuantumNumbers> m_QNvec;

	std::vector<double> m_Corr;
	std::vector<double> m_PartialZ;

	int m_BMAX, m_QMAX, m_SMAX, m_CMAX;
	int m_BMAX_list, m_QMAX_list, m_SMAX_list, m_CMAX_list;

	//bool m_OnlyBose;

	double m_MultExp;

	ThermalModelIdeal *m_modelgce;

	int m_BCE, m_QCE, m_SCE, m_CCE;

	bool m_Banalyt;
};



#endif
