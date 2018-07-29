#include "HRGBase/ThermalModelIdeal.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <cmath>



using namespace std;


ThermalModelIdeal::ThermalModelIdeal(ThermalParticleSystem *TPS_, const ThermalModelParameters& params):
        ThermalModelBase(TPS_, params)
{
	m_TAG = "ThermalModelIdeal";

	m_Ensemble = GCE;
	m_InteractionModel = Ideal;
}

ThermalModelIdeal::~ThermalModelIdeal(void)
{
}

void ThermalModelIdeal::CalculateDensities() {
	m_FluctuationsCalculated = false;

  #pragma omp parallel for if(m_useOpenMP)
	for(int i=0;i<m_TPS->Particles().size();++i) {
		m_densities[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
	}

	CalculateFeeddown();

	m_Calculated = true;
	ValidateCalculation();
}

void ThermalModelIdeal::CalculateTwoParticleCorrelations() {
	int NN = m_densities.size();
	vector<double> tN(NN), tW(NN);
	for(int i=0;i<NN;++i) tN[i] = m_densities[i];

	//#pragma omp parallel for if(m_useOpenMP)
	for(int i=0;i<NN;++i) tW[i] = CalculateParticleScaledVariance(i);

	m_PrimCorrel.resize(NN);
	for(int i=0;i<NN;++i) m_PrimCorrel[i].resize(NN);
	m_TotalCorrel = m_PrimCorrel;

	for(int i=0;i<NN;++i)
		for(int j=0;j<NN;++j) {
			m_PrimCorrel[i][j] = 0.;
			if (i==j) m_PrimCorrel[i][j] += m_densities[i] * tW[i];
			m_PrimCorrel[i][j] /= m_Parameters.T;
		}

	for(int i=0;i<NN;++i) {
		m_wprim[i] = m_PrimCorrel[i][i];
		if (m_densities[i]>0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
		else m_wprim[i] = 1.;
	}

	CalculateSusceptibilityMatrix();
	CalculateTwoParticleFluctuationsDecays();
	CalculateProxySusceptibilityMatrix();
}


void ThermalModelIdeal::CalculateFluctuations() {
	CalculateTwoParticleCorrelations();
	for(int i=0;i<m_wprim.size();++i) {
		m_wprim[i] = CalculateParticleScaledVariance(i);
		m_skewprim[i] = CalculateParticleSkewness(i);
		m_kurtprim[i] = CalculateParticleKurtosis(i);
	}
	for(int i=0;i<m_wtot.size();++i) {
		double tmp1 = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
		tmp2 = m_densities[i] * m_wprim[i];
		tmp3 = m_densities[i] * m_wprim[i] * m_skewprim[i];
		tmp4 = m_densities[i] * m_wprim[i] * m_kurtprim[i];
		for(int r=0;r<m_TPS->Particles()[i].DecayContributions().size();++r) {
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

		m_wtot[i] = tmp2 / tmp1;
		m_skewtot[i] = tmp3 / tmp2;
		m_kurttot[i] = tmp4 / tmp2;
	}

	m_FluctuationsCalculated = true;
}

std::vector<double> ThermalModelIdeal::CalculateChargeFluctuations(const std::vector<double>& chgs, int order)
{
	vector<double> ret(order + 1, 0.);

	// chi1
	for (int i = 0; i<m_densities.size(); ++i)
		ret[0] += chgs[i] * m_densities[i];

	ret[0] /= pow(m_Parameters.T * xMath::GeVtoifm(), 3);

	if (order<2) return ret;

	for (int i = 0; i<m_densities.size(); ++i)
		ret[1] += chgs[i] * chgs[i] * m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_Chem[i], 0.);

	if (order<3) return ret;

	for (int i = 0; i<m_densities.size(); ++i)
		ret[2] += chgs[i] * chgs[i] * chgs[i] * m_TPS->Particles()[i].chi(3, m_Parameters, m_UseWidth, m_Chem[i], 0.);

	if (order<4) return ret;

	for (int i = 0; i<m_densities.size(); ++i)
		ret[3] += chgs[i] * chgs[i] * chgs[i] * chgs[i] * m_TPS->Particles()[i].chi(4, m_Parameters, m_UseWidth, m_Chem[i], 0.);

	return ret;
}

double ThermalModelIdeal::CalculateEnergyDensity() {
	double ret = 0.;

	//#pragma omp parallel for reduction(+:ret) if(m_useOpenMP)
	for(int i=0;i<m_TPS->Particles().size();++i) ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);

	return ret;
}

double ThermalModelIdeal::CalculateEntropyDensity() {
	double ret = 0.;

	//#pragma omp parallel for reduction(+:ret) if(m_useOpenMP)
	for(int i=0;i<m_TPS->Particles().size();++i) ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i]);

	return ret;
}

double ThermalModelIdeal::CalculateBaryonMatterEntropyDensity() {
    double ret = 0.;
    for(int i=0;i<m_TPS->Particles().size();++i)
        if (m_TPS->Particles()[i].BaryonCharge()!=0)
            ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i]);
    return ret;
}

double ThermalModelIdeal::CalculateMesonMatterEntropyDensity() {
    double ret = 0.;
    for(int i=0;i<m_TPS->Particles().size();++i)
        if (m_TPS->Particles()[i].BaryonCharge()==0)
            ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i]);
    return ret;
}

double ThermalModelIdeal::CalculatePressure() {
	double ret = 0.;

	//#pragma omp parallel for reduction(+:ret) if(m_useOpenMP)
	for(int i=0;i<m_TPS->Particles().size();++i) ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i]);
		
	return ret;
}

double ThermalModelIdeal::CalculateShearViscosity() {
    return 0.;
}

double ThermalModelIdeal::CalculateHadronScaledVariance() {
	return 1.;
}

double ThermalModelIdeal::CalculateParticleScaledVariance(int part) {
	return m_TPS->Particles()[part].ScaledVariance(m_Parameters, m_UseWidth, m_Chem[part]);
}

double ThermalModelIdeal::CalculateParticleSkewness(int part) {
	return m_TPS->Particles()[part].Skewness(m_Parameters, m_UseWidth, m_Chem[part]);
}

double ThermalModelIdeal::CalculateParticleKurtosis(int part) {
	return m_TPS->Particles()[part].Kurtosis(m_Parameters, m_UseWidth, m_Chem[part]);
}

double ThermalModelIdeal::ParticleScalarDensity(int part) {
	return m_TPS->Particles()[part].Density(m_Parameters, IdealGasFunctions::ScalarDensity, m_UseWidth, m_Chem[part]);
}
