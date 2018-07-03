#include "HRGEV/ThermalModelEVDiagonal.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"


using namespace std;
ThermalModelEVDiagonal::ThermalModelEVDiagonal(ThermalParticleSystem *TPS_, const ThermalModelParameters& params, double RHad_, int mode) :
	ThermalModelBase(TPS_, params)
{
	m_densitiesid.resize(m_TPS->Particles().size());
	m_v.resize(m_TPS->Particles().size());
	m_Volume = params.V;
	m_TAG = "ThermalModelEVDiagonal";

	m_Ensemble = GCE;
	m_InteractionModel = DiagonalEV;
}


ThermalModelEVDiagonal::~ThermalModelEVDiagonal(void)
{
}


void ThermalModelEVDiagonal::SetRadius(double rad) {
	if (m_v.size() != m_TPS->Particles().size())
		m_v.resize(m_TPS->Particles().size());
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		m_v[i] = CuteHRGHelper::vr(rad);
	}
}

void ThermalModelEVDiagonal::FillVirial(const std::vector<double>& ri)
{
	if (ri.size() != m_TPS->Particles().size()) {
		printf("**WARNING** %s::FillVirial(const std::vector<double> & ri): size of ri does not match number of hadrons in the list", m_TAG.c_str());
		return;
	}
	m_v.resize(m_TPS->Particles().size());
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		m_v[i] = CuteHRGHelper::vr(ri[i]);
}

void ThermalModelEVDiagonal::FillVirialEV(const std::vector<double>& vi)
{
	if (vi.size() != m_TPS->Particles().size()) {
		printf("**WARNING** %s::FillVirialEV(const std::vector<double> & vi): size of vi does not match number of hadrons in the list", m_TAG.c_str());
		return;
	}
	m_v = vi;
}

void ThermalModelEVDiagonal::ReadInteractionParameters(const std::string & filename)
{
	m_v = std::vector<double>(m_TPS->Particles().size(), 0.);
	
	ifstream fin(filename);
	char cc[2000];
	while (!fin.eof()) {
		fin.getline(cc, 2000);
		string tmp = string(cc);
		vector<string> elems = CuteHRGHelper::split(tmp, '#');
		if (elems.size() < 1)
			continue;
		istringstream iss(elems[0]);
		int pdgid;
		double b;
		if (iss >> pdgid >> b) {
			int ind = m_TPS->PdgToId(pdgid);
			if (ind != -1)
				m_v[ind] = b;
		}
	}
	fin.close();
}

void ThermalModelEVDiagonal::WriteInteractionParameters(const std::string & filename)
{
	ofstream fout(filename);
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		fout << std::setw(15) << m_TPS->Particle(i).PdgId();
		fout << std::setw(15) << m_v[i];
		fout << std::endl;
	}
	fout.close();
}

double ThermalModelEVDiagonal::ExcludedVolume(int i) const
{
	if (i<0 || i >= m_v.size())
		return 0.;
	return m_v[i];
}


void ThermalModelEVDiagonal::SetParameters(const ThermalModelParameters& params) {
	m_Parameters = params;
	m_Calculated = false;
}


void ThermalModelEVDiagonal::ChangeTPS(ThermalParticleSystem *TPS_) {
	ThermalModelBase::ChangeTPS(TPS_);
	m_densitiesid.resize(m_TPS->Particles().size());
	//m_Calculated = false;
}


double ThermalModelEVDiagonal::DensityId(int i) {
	double ret = 0.;

	double dMu = -m_v[i] * m_Pressure;

	return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i], dMu);
}

double ThermalModelEVDiagonal::PressureId(int i) {
	double ret = 0.;

	double dMu = -m_v[i] * m_Pressure;

	return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i], dMu);
}

double ThermalModelEVDiagonal::ScaledVarianceId(int i) {
	double ret = 0.;

	double dMu = -m_v[i] * m_Pressure;

	return m_TPS->Particles()[i].ScaledVariance(m_Parameters, m_UseWidth, m_Chem[i], dMu);
}

double ThermalModelEVDiagonal::Pressure(double P) {
	double ret = 0.;

#pragma omp parallel for reduction(+:ret) if(useOpenMP)
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		double dMu = -m_v[i] * P;
		ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i], dMu);
	}

	return ret;
}

void ThermalModelEVDiagonal::SolvePressure() {
	const double TOLF = 1.0e-8, EPS = 1.0e-8;
	const int MAXITS = 200;
	double x = 0.;
	if (0) {
		double h = x;
		h = EPS*abs(h);
		if (h == 0.0) h = EPS;
		double xh = x + h;
		double r1 = x - Pressure(x);
		double r2 = xh - Pressure(xh);
		double Jinv = h / (r1 - r2);
		double xold = x, rold = r1;
		for (int iter = 1; iter <= MAXITS; ++iter) {
			x = xold - Jinv*rold;
			r1 = x - Pressure(x);
			if (abs(r1 / x) < TOLF) break;
			Jinv = (x - xold) / (r1 - rold);

			xold = x;
			rold = r1;
		}
		m_Pressure = x;
	}
	else {
		double mnc = pow(m_Parameters.T, 4.) * pow(xMath::GeVtoifm(), 3.);
		m_Pressure = Pressure(0.);
		x = log(m_Pressure / mnc);

		double r1 = m_Pressure - Pressure(m_Pressure);
		if (abs(r1 / m_Pressure) < TOLF) return;
		double Jinv = 0.;
		for (int i = 0; i < m_densities.size(); ++i)
			Jinv += m_v[i] * DensityId(i);
		Jinv += 1.;
		Jinv *= m_Pressure;
		Jinv = 1. / Jinv;
		double xold = x, rold = r1;
		int iter;
		for (iter = 1; iter <= MAXITS; ++iter) {
			x = xold - Jinv*rold;
			m_Pressure = mnc * exp(x);
			r1 = m_Pressure - Pressure(m_Pressure);
			if (abs(r1 / m_Pressure) < TOLF) break;

			Jinv = 0.;
			for (int i = 0; i < m_densities.size(); ++i)
				Jinv += m_v[i] * DensityId(i);
			Jinv += 1.;
			Jinv *= m_Pressure;
			Jinv = 1. / Jinv;

			xold = x;
			rold = r1;
		}
		m_Pressure = mnc * exp(x);
		if (iter == MAXITS) m_LastCalculationSuccessFlag = false;
		else m_LastCalculationSuccessFlag = true;
		m_MaxDiff = abs(r1 / m_Pressure);
	}
}

void ThermalModelEVDiagonal::CalculateDensities() {
	m_FluctuationsCalculated = false;

	SolvePressure();

	m_wnSum = 0.;
	m_Densityid = 0.;
	m_TotalDensity = 0.;
	m_Suppression = 0.;
	double densityid = 0., suppression = 0.;

	m_densitiesidnoshift = m_densitiesid;

#pragma omp parallel for reduction(+:densityid) reduction(+:suppression) if(useOpenMP)
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		double dMu = -m_v[i] * m_Pressure;
		m_densitiesid[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i], dMu);
		densityid += m_densitiesid[i];
		suppression += m_v[i] * m_densitiesid[i];

		m_densitiesidnoshift[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i], 0.);
	}

	m_Densityid = densityid;
	m_Suppression = suppression;

	m_Suppression = 1. / (1. + m_Suppression);

	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		m_densities[i] = m_densitiesid[i] * m_Suppression;
		m_TotalDensity += m_densities[i];
		m_wnSum += m_densities[i] * m_TPS->Particles()[i].ScaledVariance(m_Parameters, m_UseWidth, m_Chem[i], -m_v[i] * m_Pressure);
	}
	
	CalculateFeeddown();

	m_Calculated = true;
}

void ThermalModelEVDiagonal::CalculateTwoParticleCorrelations() {
	int NN = m_densities.size();
	vector<double> tN(NN), tW(NN);
	for (int i = 0; i < NN; ++i) tN[i] = DensityId(i);
	for (int i = 0; i < NN; ++i) tW[i] = ScaledVarianceId(i);

	m_PrimCorrel.resize(NN);
	for (int i = 0; i < NN; ++i) m_PrimCorrel[i].resize(NN);
	m_TotalCorrel = m_PrimCorrel;

	for (int i = 0; i < NN; ++i)
		for (int j = 0; j < NN; ++j) {
			m_PrimCorrel[i][j] = 0.;
			if (i == j) m_PrimCorrel[i][j] += m_densities[i] * tW[i];
			m_PrimCorrel[i][j] += -m_v[i] * m_densities[i] * m_densities[j] * tW[i];
			m_PrimCorrel[i][j] += -m_v[j] * m_densities[i] * m_densities[j] * tW[j];
			double tmp = 0.;
			for (int k = 0; k < m_densities.size(); ++k)
				tmp += m_v[k] * m_v[k] * tW[k] * m_densities[k];
			m_PrimCorrel[i][j] += m_densities[i] * m_densities[j] * tmp;
			m_PrimCorrel[i][j] /= m_Parameters.T;
		}

	for (int i = 0; i < NN; ++i) {
		m_wprim[i] = m_PrimCorrel[i][i];
		if (m_densities[i] > 0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
		else m_wprim[i] = 1.;
	}

	for (int i = 0; i < NN; ++i)
		//for(int j=0;j<NN;++j) 
	{
		m_TotalCorrel[i][i] = m_PrimCorrel[i][i];
		for (int r = 0; r < m_TPS->Particles()[i].DecayContributions().size(); ++r) {
			int rr = m_TPS->Particles()[i].DecayContributions()[r].second;
			m_TotalCorrel[i][i] += m_densities[rr] / m_Parameters.T * m_TPS->Particles()[i].DecayCumulants()[r].first[1];
			m_TotalCorrel[i][i] += 2. * m_PrimCorrel[i][rr] * m_TPS->Particles()[i].DecayContributions()[r].first;
			for (int r2 = 0; r2 < m_TPS->Particles()[i].DecayContributions().size(); ++r2) {
				int rr2 = m_TPS->Particles()[i].DecayContributions()[r2].second;
				m_TotalCorrel[i][i] += m_PrimCorrel[rr][rr2] * m_TPS->Particles()[i].DecayContributions()[r].first * m_TPS->Particles()[i].DecayContributions()[r2].first;
			}
		}
	}

	for (int i = 0; i < NN; ++i) {
		m_wtot[i] = m_TotalCorrel[i][i];
		if (m_densitiestotal[i] > 0.) m_wtot[i] *= m_Parameters.T / m_densitiestotal[i];
		else m_wtot[i] = 1.;
	}

	m_Susc.resize(4);
	for (int i = 0; i < 4; ++i) m_Susc[i].resize(4);

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			m_Susc[i][j] = 0.;
			for (int k = 0; k < m_PrimCorrel.size(); ++k) {
				int c1 = 0;
				if (i == 0) c1 = m_TPS->Particles()[k].BaryonCharge();
				if (i == 1) c1 = m_TPS->Particles()[k].ElectricCharge();
				if (i == 2) c1 = m_TPS->Particles()[k].Strangeness();
				if (i == 3) c1 = m_TPS->Particles()[k].Charm();
				for (int kp = 0; kp < m_PrimCorrel.size(); ++kp) {
					int c2 = 0;
					if (j == 0) c2 = m_TPS->Particles()[kp].BaryonCharge();
					if (j == 1) c2 = m_TPS->Particles()[kp].ElectricCharge();
					if (j == 2) c2 = m_TPS->Particles()[kp].Strangeness();
					if (j == 3) c2 = m_TPS->Particles()[kp].Charm();
					m_Susc[i][j] += c1 * c2 * m_PrimCorrel[k][kp];
				}
			}
			m_Susc[i][j] = m_Susc[i][j] / m_Parameters.T / m_Parameters.T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm();
		}
	}
}

// TODO include correlations
void ThermalModelEVDiagonal::CalculateFluctuations() {
	for (int i = 0; i < m_wprim.size(); ++i) {
		//m_wprim[i] = CalculateParticleScaledVariance(i);
		m_skewprim[i] = CalculateParticleSkewness(i);
		m_kurtprim[i] = CalculateParticleKurtosis(i);
	}
	CalculateTwoParticleCorrelations();
	m_FluctuationsCalculated = true;

	for (int i = 0; i < m_wtot.size(); ++i) {
		double tmp1 = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
		tmp2 = m_densities[i] * m_wprim[i];
		tmp3 = m_densities[i] * m_wprim[i] * m_skewprim[i];
		tmp4 = m_densities[i] * m_wprim[i] * m_kurtprim[i];
		for (int r = 0; r < m_TPS->Particles()[i].DecayContributions().size(); ++r) {
			tmp2 += m_densities[m_TPS->Particles()[i].DecayContributions()[r].second] *
				(m_wprim[m_TPS->Particles()[i].DecayContributions()[r].second] * m_TPS->Particles()[i].DecayContributions()[r].first * m_TPS->Particles()[i].DecayContributions()[r].first
					+ m_TPS->Particles()[i].DecayContributionsSigmas()[r].first);

			int rr = m_TPS->Particles()[i].DecayContributions()[r].second;
			double ni = m_TPS->Particles()[i].DecayContributions()[r].first;
			tmp3 += m_densities[rr] * m_wprim[rr] * (m_skewprim[rr] * ni * ni * ni + 3. * ni * m_TPS->Particles()[i].DecayCumulants()[r].first[1]);
			tmp3 += m_densities[rr] * m_TPS->Particles()[i].DecayCumulants()[r].first[2];

			tmp4 += m_densities[rr] * m_wprim[rr] * (m_kurtprim[rr] * ni * ni * ni * ni
				+ 6. * m_skewprim[rr] * ni * ni * m_TPS->Particles()[i].DecayCumulants()[r].first[1]
				+ 3. * m_TPS->Particles()[i].DecayCumulants()[r].first[1] * m_TPS->Particles()[i].DecayCumulants()[r].first[1]
				+ 4. * ni * m_TPS->Particles()[i].DecayCumulants()[r].first[2]);

			tmp4 += m_densities[rr] * m_TPS->Particles()[i].DecayCumulants()[r].first[3];
		}

		tmp1 = m_densitiestotal[i];

		//m_wtot[i] = tmp2 / tmp1;
		m_skewtot[i] = tmp3 / tmp2;
		m_kurttot[i] = tmp4 / tmp2;
	}
}


double ThermalModelEVDiagonal::CalculateEnergyDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	double dMu = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		dMu = -m_v[i] * m_Pressure;
		ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i], dMu);
	}
	return ret * m_Suppression;
}

double ThermalModelEVDiagonal::CalculateEntropyDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	double dMu = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		dMu = -m_v[i] * m_Pressure;
		ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i], dMu);
	}
	return ret * m_Suppression;
}

// Dummy
double ThermalModelEVDiagonal::CalculateBaryonMatterEntropyDensity() {
	double ret = 0.;
	return ret;
}
double ThermalModelEVDiagonal::CalculateMesonMatterEntropyDensity() {
	double ret = 0.;
	return ret;
}

double ThermalModelEVDiagonal::CalculatePressure() {
	if (!m_Calculated) CalculateDensities();
	return m_Pressure;
}

// TODO: Properly for multi-component
double ThermalModelEVDiagonal::CalculateHadronScaledVariance() {
	return 1.;
}

// TODO: Properly for multi-component
double ThermalModelEVDiagonal::CalculateParticleScaledVariance(int part) {
	return 1.;
}

// TODO: properly with excluded volume
double ThermalModelEVDiagonal::CalculateParticleSkewness(int part) {
	double dMu = -m_v[part] * m_Pressure;
	return m_TPS->Particles()[part].Skewness(m_Parameters, m_UseWidth, m_Chem[part], dMu);
}

// TODO: properly with excluded volume
double ThermalModelEVDiagonal::CalculateParticleKurtosis(int part) {
	double dMu = -m_v[part] * m_Pressure;
	return m_TPS->Particles()[part].Kurtosis(m_Parameters, m_UseWidth, m_Chem[part], dMu);
}


double ThermalModelEVDiagonal::CalculateBaryonScaledVariance(bool susc) {
	return 1.;
}

double ThermalModelEVDiagonal::CalculateChargeScaledVariance(bool susc) {
	return 1.;
}

double ThermalModelEVDiagonal::CalculateStrangenessScaledVariance(bool susc) {
	return 1.;
}

double ThermalModelEVDiagonal::ParticleScalarDensity(int part) {
	if (!m_Calculated) CalculateDensities();

	double dMu = -m_v[part] * m_Pressure;
	double ret = m_TPS->Particles()[part].Density(m_Parameters, IdealGasFunctions::ScalarDensity, m_UseWidth, m_Chem[part], dMu);
	return ret * m_Suppression;
}

double ThermalModelEVDiagonal::CommonSuppressionFactor()
{
	if (!m_Calculated)
		CalculateDensities();
	return m_Suppression;
}

double ThermalModelEVDiagonal::MuShift(int id)
{
	if (id >= 0. && id < m_v.size())
		return -m_v[id] * m_Pressure;
	else
		return 0.0;
}

double ThermalModelEVDiagonal::CalculateEigenvolumeFraction() {
	double tEV = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		tEV += m_v[i] * m_densities[i] / 4.;
	}
	return tEV;
}

void ThermalModelEVDiagonal::SetRadius(int i, double rad)
{
	if (i >= 0 && i < m_v.size())
		m_v[i] = CuteHRGHelper::vr(rad);
}

double ThermalModelEVDiagonal::VirialCoefficient(int i, int j) const
{
	return m_v[i];
}

void ThermalModelEVDiagonal::SetVirial(int i, int j, double b)
{
	if (i == j)
		m_v[i] = b;
}
