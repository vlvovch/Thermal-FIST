#include "HRGEV/ThermalModelEVCrossterms.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"


using namespace std;


ThermalModelEVCrossterms::ThermalModelEVCrossterms(ThermalParticleSystem *TPS_, const ThermalModelParameters& params, double RHad_, int mode) :
	ThermalModelBase(TPS_, params),
	m_RHad(RHad_),
	m_Mode(mode)
{
	m_densitiesid.resize(m_TPS->Particles().size());
	m_Volume = params.V;
	m_Ps.resize(m_TPS->Particles().size());
	FillVirial(std::vector<double>(m_TPS->Particles().size(), 0.));
	m_TAG = "ThermalModelEVCrossterms";

	m_Ensemble = GCE;
	m_InteractionModel = CrosstermsEV;
}

ThermalModelEVCrossterms::~ThermalModelEVCrossterms(void)
{
}

void ThermalModelEVCrossterms::FillVirial(const std::vector<double> & ri) {
	if (ri.size() != m_TPS->Particles().size()) {
		printf("**WARNING** %s::FillVirial(const std::vector<double> & ri): size %d of ri does not match number of hadrons %d in the list", m_TAG.c_str(), ri.size(), m_TPS->Particles().size());
		return;
	}
	m_Virial.resize(m_TPS->Particles().size());
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		m_Virial[i].resize(m_TPS->Particles().size());
		for (int j = 0; j < m_TPS->Particles().size(); ++j)
			m_Virial[i][j] = CuteHRGHelper::brr(ri[i], ri[j]);
	}

	// Correction for non-diagonal terms R1=R2 and R2=0
	std::vector< std::vector<double> > fVirialTmp = m_Virial;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		for (int j = 0; j < m_TPS->Particles().size(); ++j) {
			if (i == j) m_Virial[i][j] = fVirialTmp[i][j];
			else if ((fVirialTmp[i][i] + fVirialTmp[j][j]) > 0.0) m_Virial[i][j] = 2. * fVirialTmp[i][j] * fVirialTmp[i][i] / (fVirialTmp[i][i] + fVirialTmp[j][j]);
		}
}

void ThermalModelEVCrossterms::SetParameters(double T, double muB, double muS, double muQ, double gammaS, double V, double R) {
	m_Parameters.T = T;
	m_Parameters.muB = muB;
	m_Parameters.muS = muS;
	m_Parameters.muQ = muQ;
	m_Parameters.gammaS = gammaS;
	m_Parameters.V = V;
	m_Calculated = false;
}

void ThermalModelEVCrossterms::ReadInteractionParameters(const std::string & filename)
{
	m_Virial = std::vector< std::vector<double> >(m_TPS->Particles().size(),  std::vector<double>(m_TPS->Particles().size(), 0.));

	ifstream fin(filename);
	char cc[2000];
	while (!fin.eof()) {
		fin.getline(cc, 2000);
		string tmp = string(cc);
		vector<string> elems = CuteHRGHelper::split(tmp, '#');
		if (elems.size() < 1)
			continue;
		istringstream iss(elems[0]);
		int pdgid1, pdgid2;
		double b;
		if (iss >> pdgid1 >> pdgid2 >> b) {
			int ind1 = m_TPS->PdgToId(pdgid1);
			int ind2 = m_TPS->PdgToId(pdgid2);
			if (ind1 != -1 && ind2 != -1)
				m_Virial[ind1][ind2] = b;
		}
	}
	fin.close();
}

void ThermalModelEVCrossterms::WriteInteractionParameters(const std::string & filename)
{
	ofstream fout(filename);
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		for (int j = 0; j < m_TPS->Particles().size(); ++j) {
			fout << std::setw(15) << m_TPS->Particle(i).PdgId();
			fout << std::setw(15) << m_TPS->Particle(j).PdgId();
			fout << std::setw(15) << m_Virial[i][j];
			fout << std::endl;
		}
	}
	fout.close();
}

void ThermalModelEVCrossterms::SetRadius(double rad) {
	m_RHad = rad;
	FillVirial(vector<double>(m_TPS->Particles().size(), rad));
}

void ThermalModelEVCrossterms::DisableBBarRepulsion() {
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		for (int j = 0; j < m_TPS->Particles().size(); ++j) {
			if (m_TPS->Particles()[i].BaryonCharge() != 0 && m_TPS->Particles()[j].BaryonCharge() != 0 && m_TPS->Particles()[i].BaryonCharge()*m_TPS->Particles()[j].BaryonCharge() < 0) {
				m_Virial[i][j] = m_Virial[j][i] = 0.;
			}
		}
}

double ThermalModelEVCrossterms::VirialCoefficient(int i, int j) const {
	if (i<0 || i>=m_Virial.size() || j<0 || j>m_Virial.size())
		return 0.;
	return m_Virial[i][j];
}

void ThermalModelEVCrossterms::SetVirial(int i, int j, double b) {
	if (i >= 0 && i < m_Virial.size() && j >= 0 && j < m_Virial.size()) m_Virial[i][j] = b;
	else printf("**WARNING** Index overflow in ThermalModelEVCrossterms::SetVirial\n");
}

void ThermalModelEVCrossterms::SetParameters(const ThermalModelParameters& params) {
	m_Parameters = params;
	m_Calculated = false;
}

void ThermalModelEVCrossterms::ChangeTPS(ThermalParticleSystem *TPS_) {
	ThermalModelBase::ChangeTPS(TPS_);
	m_densitiesid.resize(m_TPS->Particles().size());
	FillVirial();
}

double ThermalModelEVCrossterms::DensityId(int i) {
	double ret = 0.;

	double dMu = 0.;
	for (int j = 0; j < m_TPS->Particles().size(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];

	return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i], dMu);
}

double ThermalModelEVCrossterms::Pressure(int i) {
	double ret = 0.;

	double dMu = 0.;
	for (int j = 0; j < m_TPS->Particles().size(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];

	return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i], dMu);
}

double ThermalModelEVCrossterms::ScaledVarianceId(int i) {
	double ret = 0.;

	double dMu = 0.;
	for (int j = 0; j < m_TPS->Particles().size(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];

	return m_TPS->Particles()[i].ScaledVariance(m_Parameters, m_UseWidth, m_Chem[i], dMu);
}

double ThermalModelEVCrossterms::PressureDiagonal(int i, double P) {
	double ret = 0.;

	double dMu = 0.;
	dMu += -m_Virial[i][i] * P;

	return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i], dMu);
}


double ThermalModelEVCrossterms::PressureDiagonalTotal(double P) {
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += PressureDiagonal(i, P);
	return ret;
}

void ThermalModelEVCrossterms::SolveDiagonal() {
	m_Pressure = 0.;
	if (1) {
		const double TOLF = 1.0e-8, EPS = 1.0e-8;
		const int MAXITS = 200;
		double x = m_Pressure;
		double h = x;
		h = EPS*abs(h);
		if (h == 0.0) h = EPS;
		double xh = x + h;
		double r1 = x - PressureDiagonalTotal(x);
		double r2 = xh - PressureDiagonalTotal(xh);
		double Jinv = h / (r2 - r1);
		double xold = x, rold = r1;
		for (int iter = 1; iter <= MAXITS; ++iter) {
			x = xold - Jinv*rold;
			r1 = x - PressureDiagonalTotal(x);
			Jinv = (x - xold) / (r1 - rold);
			if (abs(r1) < TOLF) break;
			xold = x;
			rold = r1;
		}
		m_Pressure = x;
	}
	else {
		double vo = 4. * 4. * xMath::Pi() / 3. * m_RHad * m_RHad * m_RHad;
		const double TOLF = 1.0e-8, EPS = 1.0e-8;
		const int MAXITS = 200;
		double x = 0.;
		double mnc = pow(m_Parameters.T, 4.) * pow(xMath::GeVtoifm(), 3.);
		m_Pressure = PressureDiagonalTotal(0.);
		x = log(m_Pressure / mnc);

		double r1 = m_Pressure - PressureDiagonalTotal(m_Pressure);
		if (abs(r1 / m_Pressure) < TOLF) return;
		double Jinv = 0.;
		for (int i = 0; i < m_densities.size(); ++i)
			Jinv += m_Virial[i][i] * DensityId(i);
		Jinv += 1.;
		Jinv *= m_Pressure;
		Jinv = 1. / Jinv;
		double xold = x, rold = r1;
		for (int iter = 1; iter <= MAXITS; ++iter) {
			x = xold - Jinv*rold;
			m_Pressure = mnc * exp(x);
			r1 = m_Pressure - PressureDiagonalTotal(m_Pressure);
			if (abs(r1 / m_Pressure) < TOLF) break;

			Jinv = 0.;
			for (int i = 0; i < m_densities.size(); ++i)
				Jinv += m_Virial[i][i] * DensityId(i);
			Jinv += 1.;
			Jinv *= m_Pressure;
			Jinv = 1. / Jinv;

			xold = x;
			rold = r1;
		}
		m_Pressure = mnc * exp(x);
	}
	for (int i = 0; i < m_Ps.size(); ++i)
		m_Ps[i] = PressureDiagonal(i, m_Pressure);
}

#include <Eigen/Dense>

using namespace Eigen;

void ThermalModelEVCrossterms::SolvePressure(bool resetPartials) {
	if (resetPartials) {
		m_Ps.resize(m_TPS->Particles().size());
		for (int i = 0; i < m_Ps.size(); ++i) m_Ps[i] = 0.;
		SolveDiagonal();
	}
	vector<double> Pstmp = m_Ps;
	int iter = 0;
	double maxdiff = 0.;
	const double TOLF = 1.0e-10, EPS = 1.0e-8;
	const int MAXITS = 200;
	int NNN = m_densities.size();
	VectorXd Pold(NNN), Pnew(NNN), Pdelta(NNN);
	VectorXd fold(NNN), fnew(NNN), fdelta(NNN);
	for (int i = 0; i < m_Ps.size(); ++i) {
		Pold[i] = m_Ps[i];
		Pnew[i] = m_Ps[i] + EPS*m_Ps[i];
		if (Pnew[i] == 0.0) Pnew[i] = EPS;
	}
	for (int i = 0; i < m_Ps.size(); ++i) m_Ps[i] = Pold[i];
	vector<double> tP(NNN), tN(NNN);
	for (int i = 0; i < NNN; ++i) {
		tP[i] = Pressure(i);
		tN[i] = DensityId(i);
	}
	MatrixXd Jac(NNN, NNN), Jinv(NNN, NNN);
	for (int i = 0; i < NNN; ++i) {
		for (int j = 0; j < NNN; ++j) {
			Jac(i, j) = 0.;
			if (i == j) Jac(i, j) += 1.;
			Jac(i, j) += m_Virial[i][j] * tN[i];
		}
	}
	{
		for (int i = 0; i < m_Ps.size(); ++i) m_Ps[i] = Pold[i];
		for (int i = 0; i < m_Ps.size(); ++i) fold[i] = m_Ps[i] - Pressure(i);
	}
	Pnew = Pold;
	fnew = fold;
	Jinv = Jac.inverse();
	for (iter = 1; iter < MAXITS; ++iter) {
		Pnew = Pold - Jinv * fold;
		for (int i = 0; i < m_Ps.size(); ++i) m_Ps[i] = Pnew[i];
		for (int i = 0; i < m_Ps.size(); ++i) fnew[i] = m_Ps[i] - Pressure(i);
		Pdelta = Pnew - Pold;
		fdelta = fnew - fold;
		double norm = 0.;
		for (int i = 0; i < NNN; ++i)
			for (int j = 0; j < NNN; ++j)
				norm += Pdelta[i] * Jinv(i, j) * fdelta[j];
		VectorXd p1(NNN);
		p1 = (Pdelta - Jinv*fdelta);
		for (int i = 0; i < NNN; ++i) p1[i] *= 1. / norm;
		Jinv = Jinv + (p1 * Pdelta.transpose()) * Jinv;
		maxdiff = 0.;
		for (int i = 0; i < m_Ps.size(); ++i) {
			if (m_TPS->Particles()[i].Degeneracy() > 0.0) maxdiff = max(maxdiff, fabs(fnew[i] / Pnew[i]));
		}
		if (maxdiff < TOLF) break;
		Pold = Pnew;
		fold = fnew;
	}

	m_Pressure = 0.;
	for (int i = 0; i < m_Ps.size(); ++i) m_Pressure += m_Ps[i];

	if (iter == MAXITS) m_LastCalculationSuccessFlag = false;
	else m_LastCalculationSuccessFlag = true;
	m_MaxDiff = maxdiff;
}

void ThermalModelEVCrossterms::CalculateDensities() {
	m_FluctuationsCalculated = false;

	// Pressure
	SolvePressure();
	vector<double> tN(m_densities.size());
	for (int i = 0; i < m_Ps.size(); ++i) tN[i] = DensityId(i);

	// Densities

	int NN = m_densities.size();

	MatrixXd densMatrix(NN, NN);
	VectorXd solVector(NN), xVector(NN);

	for (int i = 0; i < NN; ++i)
		for (int j = 0; j < NN; ++j) {
			densMatrix(i, j) = m_Virial[i][j] * tN[i];
			if (i == j) densMatrix(i, j) += 1.;
		}

	PartialPivLU<MatrixXd> decomp(densMatrix);

	for (int i = 0; i < NN; ++i) xVector[i] = 0.;
	for (int i = 0; i < NN; ++i) {
		xVector[i] = tN[i];
		solVector = decomp.solve(xVector);
		if (1) {
			m_densities[i] = 0.;
			for (int j = 0; j < NN; ++j)
				m_densities[i] += solVector[j];
		}
		else {
			cout << "Could not recover m_densities from partial pressures!\n";
		}
		xVector[i] = 0.;
	}

	std::vector<double> fEntropyP(m_densities.size());
	for (int i = 0; i < NN; ++i) {
		double dMu = 0.;
		for (int j = 0; j < m_TPS->Particles().size(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];
		xVector[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i], dMu);
	}

	solVector = decomp.solve(xVector);

	if (1) {
		m_TotalEntropyDensity = 0.;
		for (int i = 0; i < NN; ++i)
			m_TotalEntropyDensity += solVector[i];
	}
	else {
		cout << "**ERROR** Could not recover m_densities from partial pressures!\n";
		return;
	}

	// Decays

	CalculateFeeddown();

	m_Calculated = true;
}

void ThermalModelEVCrossterms::CalculateDensitiesNoReset() {
	// Pressure
	SolvePressure(false);
	vector<double> tN(m_densities.size());
	for (int i = 0; i < m_Ps.size(); ++i) tN[i] = DensityId(i);

	// Densities

	int NN = m_densities.size();

	MatrixXd densMatrix(NN, NN);
	VectorXd solVector(NN), xVector(NN);

	for (int i = 0; i < NN; ++i)
		for (int j = 0; j < NN; ++j) {
			densMatrix(i, j) = m_Virial[i][j] * tN[i];
			if (i == j) densMatrix(i, j) += 1.;
		}

	PartialPivLU<MatrixXd> decomp(densMatrix);

	for (int i = 0; i < NN; ++i) xVector[i] = 0.;
	for (int i = 0; i < NN; ++i) {
		xVector[i] = tN[i];//m_Ps[i] / m_Parameters.T;
		//solVector = lu.Solve(xVector, ok);
		solVector = decomp.solve(xVector);
		//if (ok) {
		if (1) {
			//if (decomp.info()==Eigen::Success) {
			m_densities[i] = 0.;
			for (int j = 0; j < NN; ++j)
				m_densities[i] += solVector[j];
		}
		else {
			cout << "**ERROR** Could not recover m_densities from partial pressures!\n";
			return;
		}
		xVector[i] = 0.;
	}

	std::vector<double> fEntropyP(m_densities.size());
	for (int i = 0; i < NN; ++i) {
		double dMu = 0.;
		for (int j = 0; j < m_TPS->Particles().size(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];
		xVector[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i], dMu);
	}

	solVector = decomp.solve(xVector);

	//if (ok) {
	if (1) {
		//if (decomp.info()==Eigen::Success) {
		m_TotalEntropyDensity = 0.;
		for (int i = 0; i < NN; ++i)
			m_TotalEntropyDensity += solVector[i];
	}
	else {
		cout << "**ERROR** Could not recover m_densities from partial pressures!\n";
		return;
	}

	// Decays

	CalculateFeeddown();

	m_Calculated = true;
}

void ThermalModelEVCrossterms::SolvePressureIter() {
	m_Ps.resize(m_TPS->Particles().size());
	for (int i = 0; i < m_Ps.size(); ++i) m_Ps[i] = 0.;
	SolveDiagonal();
	vector<double> Pstmp = m_Ps;
	int iter = 0;
	double maxdiff = 0.;
	for (iter = 0; iter < 1000; ++iter)
	{
		maxdiff = 0.;
		for (int i = 0; i < m_Ps.size(); ++i) {
			Pstmp[i] = Pressure(i);
			maxdiff = max(maxdiff, fabs((Pstmp[i] - m_Ps[i]) / Pstmp[i]));
		}
		m_Ps = Pstmp;
		//cout << iter << "\t" << maxdiff << "\n";
		if (maxdiff < 1.e-10) break;
	}
	if (iter == 1000) cout << iter << "\t" << maxdiff << "\n";
	m_Pressure = 0.;
	for (int i = 0; i < m_Ps.size(); ++i) m_Pressure += m_Ps[i];
}

void ThermalModelEVCrossterms::CalculateDensitiesIter() {
	// Pressure
	SolvePressureIter();

	int NN = m_densities.size();
	vector<double> tN(NN);
	for (int i = 0; i < NN; ++i) tN[i] = DensityId(i);

	MatrixXd densMatrix(NN, NN);
	for (int i = 0; i < NN; ++i)
		for (int j = 0; j < NN; ++j) {
			densMatrix(i, j) = m_Virial[j][i] * tN[i];
			if (i == j) densMatrix(i, j) += 1.;
		}
	//densMatrix.SetMatrixArray(matr.GetArray());

	VectorXd solVector(NN), xVector(NN);
	for (int i = 0; i < NN; ++i) xVector[i] = tN[i];

	PartialPivLU<MatrixXd> decomp(densMatrix);

	solVector = decomp.solve(xVector);

	if (1) {
		//if (decomp.info()==Eigen::Success) {
		for (int i = 0; i < NN; ++i)
			m_densities[i] = solVector[i];
	}
	else {
		cout << "**ERROR** Could not recover m_densities from partial pressures!\n";
		return;
	}

	// Entropy
	for (int i = 0; i < NN; ++i)
		for (int j = 0; j < NN; ++j) {
			densMatrix(i, j) = m_Virial[i][j] * tN[i];
			if (i == j) densMatrix(i, j) += 1.;
		}

	std::vector<double> fEntropyP(m_densities.size());
	for (int i = 0; i < NN; ++i) {
		double dMu = 0.;
		for (int j = 0; j < m_TPS->Particles().size(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];
		xVector[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i], dMu);
	}

	decomp = PartialPivLU<MatrixXd>(densMatrix);
	solVector = decomp.solve(xVector);

	if (1) {
		//if (decomp.info()==Eigen::Success) {
		m_TotalEntropyDensity = 0.;
		for (int i = 0; i < NN; ++i)
			m_TotalEntropyDensity += solVector[i];
	}
	else {
		cout << "Could not recover entropy m_densities from partial pressures!\n";
	}

	// Decays

	CalculateFeeddown();

	m_Calculated = true;
}

void ThermalModelEVCrossterms::CalculateTwoParticleCorrelations() {
	int NN = m_densities.size();
	vector<double> tN(NN), tW(NN);
	for (int i = 0; i < NN; ++i) tN[i] = DensityId(i);
	for (int i = 0; i < NN; ++i) tW[i] = ScaledVarianceId(i);
	MatrixXd densMatrix(NN, NN);
	VectorXd solVector(NN), xVector(NN), xVector2(NN);

	for (int i = 0; i < NN; ++i)
		for (int j = 0; j < NN; ++j) {
			densMatrix(i, j) = m_Virial[i][j] * tN[i];
			if (i == j) densMatrix(i, j) += 1.;
		}

	PartialPivLU<MatrixXd> decomp(densMatrix);

	vector< vector<double> > ders, coefs;

	ders.resize(NN);
	coefs.resize(NN);

	for (int i = 0; i < NN; ++i) {
		ders[i].resize(NN);
		coefs[i].resize(NN);
	}

	for (int i = 0; i < NN; ++i) xVector[i] = 0.;
	for (int i = 0; i < NN; ++i) {
		xVector[i] = tN[i];
		solVector = decomp.solve(xVector);
		if (1) {
			//if (decomp.info()==Eigen::Success) {
			for (int j = 0; j < NN; ++j) {
				ders[j][i] = solVector[j];
			}

			for (int l = 0; l < NN; ++l) {
				coefs[l][i] = 0.;
				for (int k = 0; k < NN; ++k) {
					coefs[l][i] += -m_Virial[l][k] * ders[k][i];
				}
				if (l == i) coefs[l][i] += 1.;
			}
		}
		else {
			cout << "**WARNING** Could not recover fluctuations!\n";
		}
		xVector[i] = 0.;
	}


	m_PrimCorrel.resize(NN);
	for (int i = 0; i < NN; ++i) m_PrimCorrel[i].resize(NN);
	m_TotalCorrel = m_PrimCorrel;

	for (int i = 0; i < NN; ++i)
		for (int j = i; j < NN; ++j) {
			for (int l = 0; l < NN; ++l)
				xVector[l] = tN[l] / m_Parameters.T * tW[l] * coefs[l][i] * coefs[l][j];
			solVector = decomp.solve(xVector);
			if (1) {
				//if (decomp.info()==Eigen::Success) {
				m_PrimCorrel[i][j] = 0.;
				for (int k = 0; k < NN; ++k) {
					m_PrimCorrel[i][j] += solVector[k];
				}
				m_PrimCorrel[j][i] = m_PrimCorrel[i][j];
			}
			else {
				cout << "**WARNING** Could not recover fluctuations!\n";
			}
		}

	//cout << "Primaries solved!\n";

	for (int i = 0; i < NN; ++i) {
		m_wprim[i] = m_PrimCorrel[i][i];
		if (m_densities[i] > 0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
		else m_wprim[i] = 1.;
	}

	CalculateSusceptibilityMatrix();
	CalculateTwoParticleFluctuationsDecays();
}

// TODO include correlations
void ThermalModelEVCrossterms::CalculateFluctuations() {

	CalculateTwoParticleCorrelations();

	m_FluctuationsCalculated = true;

	for (int i = 0; i < m_wprim.size(); ++i) {
		//m_wprim[i] = CalculateParticleScaledVariance(i);
		m_skewprim[i] = 1.;//CalculateParticleSkewness(i);
		m_kurtprim[i] = 1.;//CalculateParticleKurtosis(i);
	}
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

		m_skewtot[i] = tmp3 / tmp2;
		m_kurttot[i] = tmp4 / tmp2;
	}
}

double ThermalModelEVCrossterms::CalculateEnergyDensity() {
	double ret = 0.;
	ret += m_Parameters.T * CalculateEntropyDensity();
	ret += -CalculatePressure();
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_Chem[i] * m_densities[i];
	return ret;
}

double ThermalModelEVCrossterms::CalculateEntropyDensity() {
	if (!m_Calculated) CalculateDensities();
	return m_TotalEntropyDensity;
}

double ThermalModelEVCrossterms::CalculatePressure() {
	if (!m_Calculated) CalculateDensities();
	return m_Pressure;
}

double ThermalModelEVCrossterms::CalculateBaryonScaledVariance(bool susc) {
	return 1.;
}

double ThermalModelEVCrossterms::CalculateChargeScaledVariance(bool susc) {
	return 1.;
}

double ThermalModelEVCrossterms::CalculateStrangenessScaledVariance(bool susc) {
	return 1.;
}

double ThermalModelEVCrossterms::MuShift(int id)
{
	if (id >= 0. && id < m_Virial.size()) {
		double dMu = 0.;
		for (int j = 0; j < m_TPS->Particles().size(); ++j) 
			dMu += -m_Virial[id][j] * m_Ps[j];
		return dMu;
	}
	else
		return 0.0;
}
