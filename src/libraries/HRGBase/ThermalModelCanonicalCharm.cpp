#include "HRGBase/ThermalModelCanonicalCharm.h"

#include "HRGBase/xMath.h"


using namespace std;


ThermalModelCanonicalCharm::ThermalModelCanonicalCharm(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
	ThermalModelBase(TPS_, params)
{
	m_densitiesGCE.resize(m_TPS->Particles().size());

	m_StrVals.resize(0);
	m_StrVals.push_back(0);
	m_StrVals.push_back(1);
	m_StrVals.push_back(-1);

	m_StrMap.clear();
	for (unsigned int i = 0; i < m_StrVals.size(); ++i) m_StrMap[m_StrVals[i]] = i;

	m_Parameters.muC = 0.;

	m_TAG = "ThermalModelCanonicalCharm";

	m_Ensemble = CCE;
	m_InteractionModel = Ideal;

}

ThermalModelCanonicalCharm::~ThermalModelCanonicalCharm(void)
{
}


void ThermalModelCanonicalCharm::SetParameters(const ThermalModelParameters& params) {
	ThermalModelBase::SetParameters(params);
	m_Parameters.muC = 0.;
}

void ThermalModelCanonicalCharm::SetCharmChemicalPotential(double muC)
{
	m_Parameters.muC = 0.0;
}

void ThermalModelCanonicalCharm::ChangeTPS(ThermalParticleSystem *TPS_) {
	ThermalModelBase::ChangeTPS(TPS_);
	m_densitiesGCE.resize(m_TPS->Particles().size());
}

void ThermalModelCanonicalCharm::SetStatistics(bool stats) {
	m_QuantumStats = stats;
	for (unsigned int i = 0; i < m_TPS->Particles().size(); ++i)
		if (m_TPS->Particles()[i].Charm() == 0) m_TPS->Particle(i).UseStatistics(stats);
		else m_TPS->Particle(i).UseStatistics(false);
}

void ThermalModelCanonicalCharm::CalculateDensitiesGCE() {
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		m_densitiesGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
	}
	m_GCECalculated = true;
}

void ThermalModelCanonicalCharm::CalculateEnergyDensitiesGCE() {
	m_energydensitiesGCE.resize(m_densitiesGCE.size());
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		m_energydensitiesGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);
	}
}

void ThermalModelCanonicalCharm::CalculatefPhi(int iters) {
	CalculateDensitiesGCE();
	double dphi = 2 * xMath::Pi() / iters;
	m_phiRe.clearall();
	m_phiIm.clearall();
	for (unsigned int i = 0; i < iters; ++i) {
		double tphi = -xMath::Pi() + (i + 0.5) * dphi;
		double tRe = 0., tIm = 0.;
		for (unsigned int j = 0; j < m_TPS->Particles().size(); ++j) {
			if (m_TPS->Particles()[j].Strangeness() == 0) tRe += m_densitiesGCE[j];
			else {
				tRe += m_Volume * m_densitiesGCE[j] * cos(m_TPS->Particles()[j].Strangeness()*tphi);
				tIm += m_Volume * m_densitiesGCE[j] * sin(m_TPS->Particles()[j].Strangeness()*tphi);
			}
		}
		m_phiRe.add_val(tphi, exp(tRe) * cos(tIm));
		m_phiIm.add_val(tphi, exp(tRe) * sin(tIm));
	}
}

void ThermalModelCanonicalCharm::FixParameters()
{
	m_ConstrainMuC = false;
	ThermalModelBase::FixParameters();
}

void ThermalModelCanonicalCharm::CalculateDensities() {
	m_FluctuationsCalculated = false;
	m_energydensitiesGCE.resize(0);

	CalculateDensitiesGCE();

	m_Zsum.resize(m_StrVals.size());

	m_partialS.resize(m_StrVals.size());
	vector<double> xi(1, 0.), yi(1, 0.);

	for (unsigned int i = 0; i < m_StrVals.size(); ++i) {
		m_partialS[i] = 0.;
		for (unsigned int j = 0; j < m_TPS->Particles().size(); ++j)
			if (m_StrVals[i] == m_TPS->Particles()[j].Charm()) m_partialS[i] += m_densitiesGCE[j] * m_Volume;
		if (m_partialS[i] < 1.e-10) m_partialS[i] += 1e-10;
	}


	for (int i = 0; i < 1; ++i) {
		xi[i] = 2. * sqrt(m_partialS[m_StrMap[i + 1]] * m_partialS[m_StrMap[-(i + 1)]]);
		yi[i] = sqrt(m_partialS[m_StrMap[i + 1]] / m_partialS[m_StrMap[-(i + 1)]]);
	}

	for (unsigned int i = 0; i < m_StrVals.size(); ++i) {
		double res = 0.;

		res = xMath::BesselI(abs(-m_StrVals[i]), xi[0]) * pow(yi[0], m_StrVals[i]);
		m_Zsum[i] = res;
	}

	for (unsigned int i = 0; i < m_TPS->Particles().size(); ++i) {
		if (m_StrMap.count(-m_TPS->Particles()[i].Charm())) m_densities[i] = (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Charm()]] / m_Zsum[m_StrMap[0]]) * m_densitiesGCE[i];
	}

	CalculateFeeddown();

	m_Calculated = true;
	ValidateCalculation();
}

void ThermalModelCanonicalCharm::CalculateFluctuations() {
	m_wprim.resize(m_densities.size());
	m_wtot.resize(m_densities.size());
	for (int i = 0; i < m_wprim.size(); ++i) {
		m_wprim[i] = 1.;
		m_wtot[i] = 1.;
		m_skewprim[i] = 1.;
		m_skewtot[i] = 1.;
	}
}

double ThermalModelCanonicalCharm::CalculateEnergyDensity() {
	if (!m_Calculated) CalculateDensities();
	if (m_energydensitiesGCE.size() == 0) CalculateEnergyDensitiesGCE();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) ret += (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Charm()]] / m_Zsum[m_StrMap[0]]) * m_energydensitiesGCE[i];
	return ret;
}

double ThermalModelCanonicalCharm::CalculateEntropyDensity() {
	double ret = m_partialS[0] + log(m_Zsum[m_StrMap[0]]);
	if (m_energydensitiesGCE.size() == 0) CalculateEnergyDensitiesGCE();
	for (int i = 0; i < m_TPS->Particles().size(); ++i) if (m_StrMap.count(-m_TPS->Particles()[i].Charm())) ret += (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Charm()]] / m_Zsum[m_StrMap[0]]) * ((m_energydensitiesGCE[i] - (m_Parameters.muB * m_TPS->Particles()[i].BaryonCharge() + m_Parameters.muQ * m_TPS->Particles()[i].ElectricCharge() + m_Parameters.muS * m_TPS->Particles()[i].Strangeness()) * m_densitiesGCE[i]) / m_Parameters.T);
	return ret;
}


double ThermalModelCanonicalCharm::CalculatePressure() {
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) ret += (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Charm()]] / m_Zsum[m_StrMap[0]]) * m_Parameters.T * m_densitiesGCE[i];
	return ret;
}
