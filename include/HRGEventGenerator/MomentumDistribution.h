#ifndef MOMENTUMDISTRIBUTION_H
#define MOMENTUMDISTRIBUTION_H

#include <cmath>
#include <vector>

#include "HRGBase/SplineFunction.h"
#include "HRGEventGenerator/Acceptance.h"

class MomentumDistributionBase {
public:
	MomentumDistributionBase(int pdgid=0, double mass=0.):m_PDGID(pdgid), m_Mass(mass), m_Normalized(false) { }
	virtual ~MomentumDistributionBase() { }
	virtual void Normalize() = 0;
	virtual double dndp(double p) const = 0;
	virtual double dndy(double y) const = 0;
	virtual double dnmtdmt(double mt) const = 0;
	virtual double d2ndptdy(double pt, double y) const = 0;
	bool isNormalized() const { return m_Normalized; }
	void SetAcceptance(Acceptance::AcceptanceFunction *acc_, double ycm_=0.) {
		m_acc    = acc_;
		m_useacc = true;
		m_ycm    = ycm_;
	}

protected:
	int m_PDGID;
	double m_Mass;
	bool m_Normalized;
	Acceptance::AcceptanceFunction *m_acc;
	double m_ycm;
	bool m_useacc;
};


/**
* Spherically symmetric Blast-Wave Model
* P. Siemens, J. Rasmussen, Phys. Rev. Lett. 42, 880 (1979)
*/
class SiemensRasmussenDistribution : public MomentumDistributionBase {
	public:
		SiemensRasmussenDistribution(int pdgid=0, double mass=0., double T=0.100, double beta = 0.5):
			MomentumDistributionBase(pdgid,mass),
			m_T(T), m_Beta(beta)
			{
				m_Gamma = 1./sqrt(1.-m_Beta*m_Beta);
				Normalize();
				m_useacc = false;
			}
		virtual ~SiemensRasmussenDistribution() { }

		void SetParameters(double T, double beta, double mass, int pdgid = 0) {
			m_T = T;
			m_Beta = beta;
			m_Mass = mass;
			if (pdgid!=0) m_PDGID = pdgid;
			m_Gamma = 1./sqrt(1.-m_Beta*m_Beta);
			Normalize();
		}


		double w(double p) const {
			return sqrt(p*p + m_Mass*m_Mass);
		}

		double alpha(double p) const {
			return m_Gamma * m_Beta * p / m_T;
		}

		double PAv() const;

		void Normalize();

		virtual double dndp(double p) const;
		virtual double dndy(double y) const;
		virtual double dnmtdmt(double mt) const;
		virtual double d2ndptdy(double pt, double y) const;// { return 1.; }
	private:
		double m_T;
		double m_Beta;
		double m_Gamma;
		double m_Norm;
		std::vector<double> m_xlag, m_wlag;
};

/**
* Longitudinally symmetric Blast-Wave Model
* E. Schnedermann, J. Sollfrank, U. Heinz, Phys. Rev. C 48, 2462 (1993)
*/
class SSHDistribution : public MomentumDistributionBase {
public:
	SSHDistribution(int pdgid=0, double mass=0., double T=0.100, double beta = 0.5, double etamax = 0.5, double npow = 1., bool norm = false):
		MomentumDistributionBase(pdgid,mass),
		m_T(T), m_Beta(beta), m_EtaMax(etamax), m_n(npow)
		{
			m_NormY = m_NormPt = m_Norm = 1.;
			if (norm) Normalize();
			else Initialize();
			m_useacc = false;
		}
	virtual ~SSHDistribution() { }

	void SetParameters(double T, double beta, double etamax, double npow, double mass, int pdgid = 0, bool norm = true) {
		m_T = T;
		m_Beta = beta;
		m_EtaMax = etamax;
		m_n = npow;
		m_Mass = mass;
		if (pdgid!=0) m_PDGID = pdgid;
		m_NormY = m_NormPt = m_Norm = 1.;
		m_Normalized = false;
		if (norm) Normalize();
		else Initialize();
	}

	double w(double p) const {
		return sqrt(p*p + m_Mass*m_Mass);
	}

	double asinh(double x) const {
		return log(x + sqrt(1.+x*x));
	}

	double atanh(double x) const {
		return 0.5 * log((1.+x)/(1.-x));
	}

	double betar(double r) const {
		if (m_n == 1.)
			return m_Beta * r;
		else if (m_n == 2.)
			return m_Beta * r * m_Beta * r;
		else
			return pow(m_Beta * r, m_n);
	}

	double rho(double r) const {
		return atanh(betar(r));
	}

	void Initialize();
	void Normalize();

	virtual double dndp(double p) const {
		return 0.;
	}

	virtual double dndy(double y) const;
	virtual double dndysingle(double y) const;
	virtual double dnmtdmt(double mt) const;
	virtual double dndy(double y, double pt) const;
	virtual double dndysingle(double y, double pt) const;
	virtual double dndpt(double pt) const;
	virtual double dndpt(double pt, double y) const;
	virtual double dndptsingle(double pt, double y) const;
	virtual double d2ndptdy(double pt, double y) const;

	double MtAv() const;
	double y2Av() const;
private:
	double m_T;
	double m_Beta, m_EtaMax;
	double m_NormY, m_NormPt, m_Norm;
	double m_n;
	std::vector<double> m_xlag, m_wlag;
	std::vector<double> m_xlegT, m_wlegT;
	std::vector<double> m_xlegY, m_wlegY;
	std::vector<double> m_xlegeta, m_wlegeta;

	SplineFunction m_dndy, m_dndyint;
};

#endif
