#ifndef THERMALMODELCANONICALBSQPLUS_H
#define THERMALMODELCANONICALBSQPLUS_H

#include <map>

#include "HRGBase/ThermalModelBase.h"

struct QuantumNumbers
{
	int B, Q, S, C;
	QuantumNumbers(int iB = 0, int iQ = 0, int iS = 0, int iC = 0) :
		B(iB), Q(iQ), S(iS), C(iC) { }
	const bool operator < (const QuantumNumbers &r) const {
		if (B != r.B)
			return (B < r.B);
		else if (Q!= r.Q)
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

		~ThermalModelCanonical(void);

		void ChangeTPS(ThermalParticleSystem *TPS_);

		void CalculateQuantumNumbersRange(bool doubleRange = false);

		void SetParameters(double T, double gammaS, double V, int B, int Q, int S = 0, int C = 0);

		void SetParameters(const ThermalModelParameters& params);


		virtual void SetStatistics(bool stats);

		virtual void FixParameters() { };
		virtual void FixParameters(double) { };
		virtual void FixParametersNoReset() { };

		void CalculateDensities();
		void CalculatePartitionFunctions();
		void CalculatePartitionFunctionsBoseOnly();

		// TODO check the validity
		virtual double CalculateParticleScaledVariance(int part);
		void CalculateTwoParticleCorrelations();

		// TODO properly for higher moments
		void CalculateFluctuations();


		virtual double CalculateEnergyDensity();

		virtual double CalculatePressure();

		virtual double CalculateEntropyDensity();

		double GetGCEDensity(int i) const;

		virtual double ParticleScalarDensity(int part) { return 0.; }

	private:

		std::map<QuantumNumbers, int> m_QNMap;
		std::vector<QuantumNumbers> m_QNvec;

		std::vector<double> m_Corr;
		std::vector<double> m_PartialZ;

		int m_BMAX, m_QMAX, m_SMAX, m_CMAX;

		bool m_OnlyBose;
};

#endif
