#include "HRGBase/ThermalParticle.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"

using namespace std;

const double lowTlimit = 10.;


ThermalParticle::ThermalParticle(bool Stable_, std::string Name, int PDGID, double Deg, int Stat, double Mass,
								 int Strange, int Baryon, int Charge, double AbsS, double Width, double Threshold, int Charm, double AbsC, double radius, int Quark):
			m_Stable(Stable_), m_AntiParticle(false), m_Name(Name), m_PDGID(PDGID), m_Degeneracy(Deg), m_Statistics(Stat), m_StatisticsOrig(Stat), m_Mass(Mass),
			m_Strangeness(Strange), m_Baryon(Baryon), m_ElectricCharge(Charge), m_Charm(Charm), m_ArbitraryCharge(Baryon), m_AbsS(AbsS), m_AbsC(AbsC), m_Width(Width), m_Threshold(Threshold), m_Radius(radius), m_Quark(Quark), m_Weight(1.)
{
	SetCalculationType(IdealGasFunctions::ClusterExpansion);
	
	SetClusterExpansionOrder(3);
	if (m_Mass < 1.000) SetClusterExpansionOrder(5);
	if (m_Mass < 0.200) SetClusterExpansionOrder(10);
	
	SetResonanceWidthShape(RelativisticBreitWiger);
	SetResonanceWidthIntegrationType(TwoGamma);
	
	m_DecayContributions.resize(0);
	m_DecayContributionsSigmas.resize(0);
	m_DecayProbabilities.resize(0);

	FillCoefficients();

	SetAbsoluteQuark(GetAbsQ());

}
ThermalParticle::~ThermalParticle(void)
{
}

void ThermalParticle::SetResonanceWidth(double width)
{
	m_Width = width;
	if (m_Width != 0.0) FillCoefficients();
}

void ThermalParticle::SetDecayThresholdMass(double threshold)
{
	m_Threshold = threshold;
	if (m_Threshold < 0.0) {
		printf("**WARNING** Trying to set negative decay threshold for %s, setting to zero instead", m_Name.c_str());
	}
	if (m_Width != 0.0) FillCoefficients();
}

void ThermalParticle::SetResonanceWidthIntegrationType(ResonanceWidthIntegration type)
{
	m_ResonanceWidthIntegrationType = type;
	FillCoefficients();
}

double ThermalParticle::MassDistribution(double m) const
{
	if (m_ResonanceWidthShape==RelativisticBreitWiger) 
		return m_Mass * m_Width * m / ((m * m - m_Mass*m_Mass)*(m * m - m_Mass*m_Mass) + m_Mass*m_Mass*m_Width*m_Width);
	else 
		return 1. / ((m - m_Mass)*(m - m_Mass) + m_Width*m_Width / 4.);
}

void ThermalParticle::ReadDecays(string filename) {
	m_Decays.resize(0);
	ifstream fin(filename.c_str());
	if (fin.is_open()) {
		char cc[400];
		double tmpbr;
		while (fin >> tmpbr) {
			ParticleDecay decay;
			decay.mBratio = tmpbr / 100.;
			fin.getline(cc, 350);
			stringstream ss;
			ss << cc;
			int tmpid;
			while (ss >> tmpid) {
				decay.mDaughters.push_back(tmpid);
			}
			m_Decays.push_back(decay);
		}
	}
}

std::vector<double> ThermalParticle::BranchingRatioWeights(const std::vector<double>& ms) const
{
	std::vector<double> ret = ms;

	std::vector<double> mthr, br;
	double tsum = 0.;
	for (int i = 0; i < m_Decays.size(); ++i) {
		mthr.push_back(m_Decays[i].mM0);
		br.push_back(m_Decays[i].mBratio);
		tsum += m_Decays[i].mBratio;
	}
	if (tsum < 1.) {
		mthr.push_back(0.);
		br.push_back(1. - tsum);
	}
	else if (tsum > 1.) {
		for (int i = 0; i < br.size(); ++i)
			br[i] *= 1. / tsum;
	}

	for (int i = 0; i < ms.size(); ++i) {
		double tw = 0.;
		for (int j = 0.; j < br.size(); ++j) {
			if (mthr[j] <= ms[i])
				tw += br[j];
		}
		ret[i] = tw;
	}
	return ret;
}

void ThermalParticle::NormalizeBranchingRatios() {
	double sum = 0.;
	for(int i=0;i<m_Decays.size();++i) {
		sum += m_Decays[i].mBratio;
	}
	for(int i=0;i<m_Decays.size();++i) {
		m_Decays[i].mBratio *= 1./sum;
	}
}

void ThermalParticle::FillCoefficients() {
	double a, b;
	if (m_ResonanceWidthIntegrationType != TwoGamma && m_Threshold >= 0.) {
		a = m_Threshold;
		b = m_Mass + 2.*m_Width;
		NumericalIntegration::GetCoefsIntegrateLegendre32(a, b, &m_xleg, &m_wleg);
		m_brweight = BranchingRatioWeights(m_xleg);
	}
	else {
		a = max(m_Threshold, m_Mass - 2.*m_Width);
		b = m_Mass + 2.*m_Width;
		NumericalIntegration::GetCoefsIntegrateLegendre10(a, b, &m_xleg, &m_wleg);
		m_brweight = BranchingRatioWeights(m_xleg);
	}

	// Old version
	//if (m_Width / m_Mass<1e-1) { NumericalIntegration::GetCoefsIntegrateLegendre10(a, b, &m_xleg, &m_wleg); }
	//else { NumericalIntegration::GetCoefsIntegrateLegendre32(a, b, &m_xleg, &m_wleg); }
	
	// New version
	NumericalIntegration::GetCoefsIntegrateLegendre32(0., 1., &m_xleg32, &m_wleg32);
	NumericalIntegration::GetCoefsIntegrateLaguerre32(&m_xlag32, &m_wlag32);
}

void ThermalParticle::UseStatistics(bool enable) {
	if (!enable) m_Statistics = 0;
	else m_Statistics = m_StatisticsOrig;
}

void ThermalParticle::SetMass(double mass)
{
	m_Mass = mass;
	if (m_Width != 0.0) FillCoefficients();
}

bool ThermalParticle::IsNeutral() const {
	return (m_Baryon==0 && m_ElectricCharge==0 && m_Strangeness==0 && m_Charm==0);
}


double ThermalParticle::Density(const ThermalModelParameters &params, IdealGasFunctions::Quantity type, bool useWidth, double pMu, double dMu) const {
	double mu = pMu + dMu;
	if (!(params.gammaq == 1.))									mu += log(params.gammaq) * m_AbsQuark * params.T;
	if (!(params.gammaS == 1. || m_AbsS == 0.))	mu += log(params.gammaS) * m_AbsS     * params.T;
	if (!(params.gammaC == 1. || m_AbsC == 0.))	mu += log(params.gammaC) * m_AbsC     * params.T;
	
	if (!useWidth || m_Width/m_Mass < 1.e-2) {
		return IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, m_Statistics, params.T, mu, m_Mass, m_Degeneracy, m_ClusterExpansionOrder);
	}

	int ind = m_xleg.size();
	const vector<double> &x = m_xleg, &w = m_wleg;
	double ret1 = 0., ret2 = 0., tmp = 0.;

	// Integration from m0 or M-2*Gamma to M+2*Gamma
	for (int i = 0; i < ind; i++) {

		tmp = w[i] * MassDistribution(x[i]);

		if (m_ResonanceWidthIntegrationType == FullIntervalWeighted) 
			tmp *= m_brweight[i];

		double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, m_Statistics, params.T, mu, x[i], m_Degeneracy, m_ClusterExpansionOrder);

		ret1 += tmp * dens;
		ret2 += tmp;
	}
	
	// Integration from M+2*Gamma to infinity
	if (m_ResonanceWidthIntegrationType == FullInterval) {
		int ind2 = m_xlag32.size();
		for (int i = 0; i < ind2; ++i) {
			double tmass = m_Mass + 2.*m_Width + m_xlag32[i] * m_Width;
			tmp = m_wlag32[i] * m_Width * MassDistribution(tmass);
			double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, m_Statistics, params.T, mu, tmass, m_Degeneracy, m_ClusterExpansionOrder);

			ret1 += tmp * dens;
			ret2 += tmp;
		}
	}

	return ret1 / ret2;
}

double ThermalParticle::DensityCluster(int n, const ThermalModelParameters & params, IdealGasFunctions::Quantity type, bool useWidth, double pMu, double dMu) const
{
	double mn = 1.;
	if ((abs(BaryonCharge()) & 1) && !(n & 1))
		mn = -1.;

	double mu = pMu + dMu;
	if (!(params.gammaq == 1.))									mu += log(params.gammaq) * m_AbsQuark /*GetAbsQ()*/  * params.T;
	if (!(params.gammaS == 1. || m_AbsS == 0.))	mu += log(params.gammaS) * m_AbsS     * params.T;
	if (!(params.gammaC == 1. || m_AbsC == 0.))	mu += log(params.gammaC) * m_AbsC     * params.T;

	if (!useWidth || m_Width / m_Mass < 1.e-2) {
		return mn * IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, 0, params.T / static_cast<double>(n), mu, m_Mass, m_Degeneracy);
	}

	int ind = m_xleg.size();
	const vector<double> &x = m_xleg, &w = m_wleg;
	double ret1 = 0., ret2 = 0., tmp = 0.;

	// Integration from m0 or M-2*Gamma to M+2*Gamma
	for (int i = 0; i < ind; i++) {
		tmp = w[i] * MassDistribution(x[i]);

		if (m_ResonanceWidthIntegrationType == FullIntervalWeighted)
			tmp *= m_brweight[i];

		double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, 0, params.T / static_cast<double>(n), mu, x[i], m_Degeneracy);

		ret1 += tmp * dens;
		ret2 += tmp;
	}

	// Integration from M+2*Gamma to infinity
	if (m_ResonanceWidthIntegrationType == FullInterval) {
		int ind2 = m_xlag32.size();
		for (int i = 0; i < ind2; ++i) {
			double tmass = m_Mass + 2.*m_Width + m_xlag32[i] * m_Width;
			tmp = m_wlag32[i] * m_Width * MassDistribution(tmass);
			double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, 0, params.T / static_cast<double>(n), mu, tmass, m_Degeneracy);

			ret1 += tmp * dens;
			ret2 += tmp;
		}
	}

	return mn * ret1 / ret2;
}


double ThermalParticle::ScaledVariance(const ThermalModelParameters &params, bool useWidth, double pMu, double dMu) const {
	if (m_Degeneracy == 0) return 1.;
	if (m_Statistics == 0) return 1.;
	double dens = Density(params, IdealGasFunctions::ParticleDensity, useWidth, pMu, dMu);
	if (dens==0.) return 1.;
	double ret = chi(2, params, useWidth, pMu, dMu) / chi(1, params, useWidth, pMu, dMu);
	if (ret != ret) ret = 1.;
	return ret;
}

double ThermalParticle::Skewness(const ThermalModelParameters &params, bool useWidth, double pMu, double dMu) const
{
	if (m_Degeneracy == 0) return 1.;
	if (m_Statistics == 0) return 1.;
	double dens = Density(params, IdealGasFunctions::ParticleDensity, useWidth, pMu, dMu);
	if (dens == 0.) return 1.;
	double ret = chi(3, params, useWidth, pMu, dMu) / chi(2, params, useWidth, pMu, dMu);
	if (ret != ret) ret = 1.;
	return ret;
}

double ThermalParticle::Kurtosis(const ThermalModelParameters &params, bool useWidth, double pMu, double dMu) const
{
	if (m_Degeneracy == 0) return 1.;
	if (m_Statistics == 0) return 1.;
	double dens = Density(params, IdealGasFunctions::ParticleDensity, useWidth, pMu, dMu);
	if (dens == 0.) return 1.;
	double ret = chi(4, params, useWidth, pMu, dMu) / chi(2, params, useWidth, pMu, dMu);
	if (ret != ret) ret = 1.;
	return ret;
}

double ThermalParticle::FD(double k, double T, double mu, double m) const {
	double arg = (sqrt(k*k+m*m)-mu)/T;
	return 1./(exp(arg)+1.);
}

double ThermalParticle::GetAbsQ() const {
	if (m_Baryon==0) return 2. - m_AbsC - m_AbsS;
	else return abs(m_Baryon) * (3. - m_AbsC - m_AbsS);
}

double ThermalParticle::GetCharge(int index) const {
	if (index==0) return (double)m_Baryon;
	if (index==1) return (double)m_ElectricCharge;
	if (index==2) return (double)m_Strangeness;
	if (index==3) return (double)m_Charm;
	return 0.;
}

double ThermalParticle::GetAbsCharge(int index) const {
	if (index==0) return fabs((double)m_Baryon);
	if (index==1) return fabs((double)m_ElectricCharge);
	if (index==2) return m_AbsS;
	if (index==3) return m_AbsC;
	return 0.;
}

double ThermalParticle::chi(int index, const ThermalModelParameters &params, bool useWidth, double pMu, double dMu) const {
	if (index == 0) return Density(params, IdealGasFunctions::Pressure, useWidth, pMu, dMu) / pow(params.T, 4) / pow(xMath::GeVtoifm(), 3);
	if (index == 1) return Density(params, IdealGasFunctions::ParticleDensity, useWidth, pMu, dMu) / pow(params.T, 3) / pow(xMath::GeVtoifm(), 3);
	if (index == 2) return Density(params, IdealGasFunctions::chi2, useWidth, pMu, dMu);
	if (index == 3) return Density(params, IdealGasFunctions::chi3, useWidth, pMu, dMu);
	if (index == 4) return Density(params, IdealGasFunctions::chi4, useWidth, pMu, dMu);
	return 1.;
}
