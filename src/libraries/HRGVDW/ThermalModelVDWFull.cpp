#include "HRGVDW/ThermalModelVDWFull.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;

ThermalModelVDWFull::ThermalModelVDWFull(ThermalParticleSystem *TPS_, const ThermalModelParameters& params, int type):
    ThermalModelBase(TPS_, params), m_SearchMultipleSolutions(false), m_TemperatureDependentAB(false)
{
	m_chi.resize(6);
	for(int i=0;i<6;++i) m_chi[i].resize(3);
	m_chiarb.resize(6);
	m_DensitiesId.resize(m_densities.size());
	m_MuStar.resize(m_densities.size());
	m_Virial.resize(m_densities.size(), vector<double>(m_densities.size(),0.));
	m_Attr = m_Virial;
	m_VirialdT = m_Virial;
	m_AttrdT   = m_AttrdT;
	m_Volume = params.V;
	m_TAG = "ThermalModelVDWFull";

	m_Ensemble = GCE;
	m_InteractionModel = QvdW;
}


ThermalModelVDWFull::~ThermalModelVDWFull(void)
{
}

void ThermalModelVDWFull::UpdateParameters() {
	SetParameters(m_Parameters);
}

void ThermalModelVDWFull::FillChemicalPotentials() {
	ThermalModelBase::FillChemicalPotentials();
	for(int i=0; i<m_MuStar.size(); ++i) 
		m_MuStar[i] = m_Chem[i];
}

void ThermalModelVDWFull::SetChemicalPotentials(const std::vector<double>& chem)
{
	ThermalModelBase::SetChemicalPotentials(chem);
	for (int i = 0; i<m_MuStar.size(); ++i)
		m_MuStar[i] = m_Chem[i];
}

void ThermalModelVDWFull::FillVirial(const std::vector<double> & ri) {
	if (ri.size() != m_TPS->Particles().size()) {
		printf("**WARNING** %s::FillVirial(const std::vector<double> & ri): size of ri does not match number of hadrons in the list", m_TAG.c_str());
		return;
	}
	m_Virial.resize(m_TPS->Particles().size());
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		m_Virial[i].resize(m_TPS->Particles().size());
		for (int j = 0; j < m_TPS->Particles().size(); ++j)
			m_Virial[i][j] = CuteHRGHelper::brr(ri[i], ri[j]);
	}

	// Correct R1=R2 and R2=0
	std::vector< std::vector<double> > fVirialTmp = m_Virial;
	for(int i=0;i<m_TPS->Particles().size();++i)
		for(int j=0;j<m_TPS->Particles().size();++j) {
			if (i==j) m_Virial[i][j] = fVirialTmp[i][j];
			else if ((fVirialTmp[i][i] + fVirialTmp[j][j]) > 0.0) m_Virial[i][j] = 2. * fVirialTmp[i][j] * fVirialTmp[i][i] / (fVirialTmp[i][i] + fVirialTmp[j][j]);
		}

	m_VirialdT = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(),0.));
	m_AttrdT   = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));
}

void ThermalModelVDWFull::FillVirialEV(const std::vector< std::vector<double> >& bij)
{
	if (bij.size() != m_TPS->Particles().size()) {
		printf("**WARNING** %s::FillVirialEV(const std::vector<double> & bij): size of bij does not match number of hadrons in the list", m_TAG.c_str());
		return;
	}
	m_Virial = bij;
}

void ThermalModelVDWFull::FillAttraction(const std::vector<std::vector<double> >& aij)
{
	if (aij.size() != m_TPS->Particles().size()) {
		printf("**WARNING** %s::FillAttraction(const std::vector<double> & aij): size of aij does not match number of hadrons in the list", m_TAG.c_str());
		return;
	}
	m_Attr = aij;
}

void ThermalModelVDWFull::ReadInteractionParameters(const std::string & filename)
{
	m_Virial = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));
	m_Attr   = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));

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
		double b, a;
		if (iss >> pdgid1 >> pdgid2 >> b) {
			if (!(iss >> a))
				a = 0.;
			int ind1 = m_TPS->PdgToId(pdgid1);
			int ind2 = m_TPS->PdgToId(pdgid2);
			if (ind1 != -1 && ind2 != -1) {
				m_Virial[ind1][ind2] = b;
				m_Attr[ind1][ind2]   = a;
			}
		}
	}
	fin.close();
}

void ThermalModelVDWFull::WriteInteractionParameters(const std::string & filename)
{
	ofstream fout(filename);
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		for (int j = 0; j < m_TPS->Particles().size(); ++j) {
			fout << std::setw(15) << m_TPS->Particle(i).PdgId();
			fout << std::setw(15) << m_TPS->Particle(j).PdgId();
			fout << std::setw(15) << m_Virial[i][j];
			fout << std::setw(15) << m_Attr[i][j];
			fout << std::endl;
		}
	}
	fout.close();
}

void ThermalModelVDWFull::SetParameters(double T, double muB, double muS, double muQ, double gammaS, double V, double R) {
	m_Parameters.T = T;
	m_Parameters.muB = muB;
	m_Parameters.muS = muS;
	m_Parameters.muQ = muQ;
	m_Parameters.gammaS = gammaS;
	//m_Parameters.R = R;
    m_Parameters.V  = V;
	//m_Parameters.Vc = V;
	m_Calculated = false;
}

void ThermalModelVDWFull::SetParameters(const ThermalModelParameters& params) {
	m_Parameters = params;
	m_Calculated = false;
}

void ThermalModelVDWFull::ChangeTPS(ThermalParticleSystem *TPS_) {
    ThermalModelBase::ChangeTPS(TPS_);
		m_Virial = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));
		m_Attr = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));
		m_VirialdT = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));
		m_AttrdT = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));
}

#include <Eigen/Dense>

using namespace Eigen;

vector<double> ThermalModelVDWFull::SearchSolutionsSingle(const vector<double> &muStarInit) {
	const double EPS = 1.0e-10;

	int NN = m_densities.size();

	MatrixXd densMatrix(NN, NN);
	VectorXd solVector(NN), xVector(NN);

	vector<double> muscur = muStarInit;
	vector<double> musold;
	int iters = 0;
	double maxdif = 0.;
	const int MAXITS = 200;
	while (iters<MAXITS) {
		musold = muscur;
		vector<double> Ps(m_densities.size(), 0.);
		for(int i=0;i<NN;++i) Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, muscur[i], 0.);
		vector<double> ns(m_densities.size(), 0.);
		for(int i=0;i<NN;++i) ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, muscur[i], 0.);

		for(int i=0;i<NN;++i)
			for(int j=0;j<NN;++j) {
				densMatrix(i,j) = m_Virial[j][i] * ns[i];
				if (i==j) densMatrix(i,j) += 1.;
			}

		PartialPivLU<MatrixXd> decomp(densMatrix);

		for(int i=0;i<NN;++i) xVector[i] = ns[i];
		solVector = decomp.solve(xVector);

		for(int i=0;i<NN;++i) {
			muscur[i] = m_Chem[i];
			for(int j=0;j<NN;++j) muscur[i] += -m_Virial[i][j] * Ps[j];
			for(int j=0;j<NN;++j) muscur[i] += (m_Attr[i][j] + m_Attr[j][i]) * solVector[j];
		}

		
		for(int i=0;i<NN;++i) maxdif = max(maxdif, abs(muscur[i] - musold[i]));

		iters++;

		if (maxdif < EPS) break;
	}
	if (iters==MAXITS) m_LastBroydenSuccessFlag = false;
	else m_LastBroydenSuccessFlag = true;
	m_MaxDiff = maxdif;
	//printf("Iteration %d: %lf\n", iters, maxdif);
	return muscur;
}

vector<double> ThermalModelVDWFull::SearchSolutionsSingleBroyden(const vector<double> &muStarInit) {
	

	const double EPS = 1.0e-10;

	int NN = m_densities.size();

	bool attrfl = false;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			if (m_Attr[i][j] != 0.0) {
				attrfl = true;
				break;
			}
		}
		if (attrfl) break;
	}

	MatrixXd densMatrix(NN, NN);
	VectorXd solVector(NN), xVector(NN);

	vector<double> muscur = muStarInit;

	VectorXd musold(NN), musnew(NN), musdelta(NN);
	VectorXd fold(NN), fnew(NN), fdelta(NN);
	for(int i=0;i<NN;++i) { 
		musold[i] = muscur[i];
	}

	MatrixXd Jac(NN, NN), Jinv(NN, NN);
	{
		vector<double> Ps(m_densities.size(), 0.);
		for(int i=0;i<NN;++i) Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, musold[i], 0.);
		vector<double> ns(m_densities.size(), 0.);
		for(int i=0;i<NN;++i) ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, musold[i], 0.);
		vector<double> chi2s(m_densities.size(), 0.);
		for(int i=0;i<NN;++i) chi2s[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, musold[i], 0.);
		

		for(int i=0;i<NN;++i)
			for(int j=0;j<NN;++j) {
				densMatrix(i,j) = m_Virial[j][i] * ns[i];
				if (i==j) densMatrix(i,j) += 1.;
			}

		//printf("Before n-matrix LU-decomposition\n");
		clock_t tbeg = clock();
		PartialPivLU<MatrixXd> decomp(densMatrix);
		//printf("After n-matrix LU-decomposition. Time: %lf s\n", (clock() - tbeg) / static_cast<double>(CLOCKS_PER_SEC));


		for(int i=0;i<NN;++i) xVector[i] = ns[i];
		vector<double> np(m_densities.size(), 0.);
		solVector = decomp.solve(xVector);
		for(int i=0;i<NN;++i)
			np[i] = solVector[i];


		
		for(int i=0;i<NN;++i) xVector[i] = 0.;

		for(int k=0;k<NN;++k) {
			double tmps = 0.;
			for(int j=0;j<NN;++j) 
				tmps += m_Virial[j][k] * np[j];
			xVector[k] = chi2s[k] * m_Parameters.T * m_Parameters.T * pow(xMath::GeVtoifm(), 3) * (1. - tmps);
			solVector = decomp.solve(xVector);
			for(int i=0;i<NN;++i)
				if (solVector[i]>1.) solVector[i] = 1.;	// Stabilizer
			xVector[k] = 0.;
			for(int i=0;i<NN;++i) {
				Jac(i,k) = 0.;
				if (i==k) Jac(i,k) += 1.;
				Jac(i,k) += m_Virial[i][k] * ns[k];
				if (attrfl) 
					for(int j=0;j<NN;++j)
						Jac(i,k) += -(m_Attr[i][j] + m_Attr[j][i]) * solVector[j];
			}
		}

		{
			for(int i=0;i<fold.size();++i) {
				fold[i] = musold[i];
				for(int j=0;j<NN;++j) 
					fold[i] += m_Virial[i][j] * Ps[j] - (m_Attr[i][j] + m_Attr[j][i]) * np[j];
				fold[i] += -m_Chem[i];
			}
		}
	}
	musnew = musold;
	fnew   = fold;
	//printf("Before matrix inversion\n");
	clock_t tbeg = clock();
	Jinv = Jac.inverse();
	//printf("After matrix inversion. Time: %lf s\n", (clock() - tbeg) / static_cast<double>(CLOCKS_PER_SEC));

	int iters = 1;
	double maxdif = 0.;
	const int MAXITS = 200;
	while (iters<MAXITS) {
		musnew = musold - Jinv * fold;

		vector<double> Ps(m_densities.size(), 0.);
		for(int i=0;i<NN;++i) Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, musnew[i], 0.);
		vector<double> ns(m_densities.size(), 0.);
		for(int i=0;i<NN;++i) ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, musnew[i], 0.);

		for(int i=0;i<NN;++i)
			for(int j=0;j<NN;++j) {
				densMatrix(i,j) = m_Virial[j][i] * ns[i];
				if (i==j) densMatrix(i,j) += 1.;
			}

		PartialPivLU<MatrixXd> decomp(densMatrix);

		for(int i=0;i<NN;++i) xVector[i] = ns[i];
		vector<double> np(m_densities.size(), 0.);
		solVector = decomp.solve(xVector);
		for(int i=0;i<NN;++i)
			np[i] = solVector[i];

		for(int i=0;i<fnew.size();++i) {
			fnew[i] = musnew[i];
			for(int j=0;j<NN;++j) 
				fnew[i] += m_Virial[i][j] * Ps[j] - (m_Attr[i][j] + m_Attr[j][i]) * np[j];
			fnew[i] += -m_Chem[i];
		}

		musdelta = musnew - musold;
		fdelta   = fnew - fold;
		double norm = 0.;
		for(int i=0;i<NN;++i)
			for(int j=0;j<NN;++j)
				norm += musdelta[i] * Jinv(i,j) * fdelta[j];
		VectorXd p1(NN);
		p1 = (musdelta - Jinv*fdelta);
		for(int i=0;i<NN;++i) p1[i] *= 1. / norm;
		Jinv = Jinv + (p1 * musdelta.transpose()) * Jinv;

		if (0) {
			vector<double> Ps(m_densities.size(), 0.);
			for(int i=0;i<NN;++i) Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, musnew[i], 0.);
			vector<double> ns(m_densities.size(), 0.);
			for(int i=0;i<NN;++i) ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, musnew[i], 0.);
			vector<double> chi2s(m_densities.size(), 0.);
			for(int i=0;i<NN;++i) chi2s[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, musnew[i], 0.);
		

			for(int i=0;i<NN;++i)
				for(int j=0;j<NN;++j) {
					densMatrix(i,j) = m_Virial[j][i] * ns[i];
					if (i==j) densMatrix(i,j) += 1.;
				}

			PartialPivLU<MatrixXd> decomp(densMatrix);


			for(int i=0;i<NN;++i) xVector[i] = ns[i];
			vector<double> np(m_densities.size(), 0.);
			solVector = decomp.solve(xVector);
			for(int i=0;i<NN;++i)
				np[i] = solVector[i];


		
			for(int i=0;i<NN;++i) xVector[i] = 0.;

			for(int k=0;k<NN;++k) {
				double tmps = 0.;
				for(int j=0;j<NN;++j) 
					tmps += m_Virial[j][k] * np[j];
				xVector[k] = chi2s[k] * m_Parameters.T * m_Parameters.T * pow(xMath::GeVtoifm(), 3) * (1. - tmps);
				solVector = decomp.solve(xVector);
				xVector[k] = 0.;
				for(int i=0;i<NN;++i) {
					Jac(i,k) = 0.;
					if (i==k) Jac(i,k) += 1.;
					Jac(i,k) += m_Virial[i][k] * ns[k];
					for(int j=0;j<NN;++j) 
						Jac(i,k) += -(m_Attr[i][j] + m_Attr[j][i]) * solVector[j];
				}
			}

			Jinv = Jac.inverse();
		}
		
		maxdif = 0.;
		for(int i=0;i<musnew.size();++i) {
			if (m_TPS->Particles()[i].Degeneracy()>0.0) maxdif = max(maxdif, abs(musnew[i] - musold[i]));
		}
		iters++;
		//printf("Iteration %d: %E\n", iters, maxdif);

		if (maxdif<EPS) break;
		musold = musnew;
		fold   = fnew;
	}
	//printf("Iteration %d: %E\n", iters, maxdif);
	for(int i=0; i<muscur.size(); ++i)
		muscur[i] = musnew[i];
	if (iters==MAXITS) m_LastBroydenSuccessFlag = false;
	else m_LastBroydenSuccessFlag = true;
	m_MaxDiff = maxdif;
	return muscur;
}

std::vector<double> ThermalModelVDWFull::SearchSolutionsSingleBroydenOptimized(const std::vector<double>& muStarInit)
{
	const double EPS = 1.0e-10;

	int NN = m_densities.size();

	bool attrfl = false;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			if (m_Attr[i][j] != 0.0) {
				attrfl = true;
				break;
			}
		}
		if (attrfl) break;
	}

	

	int NNdmu = m_MapFromdMuStar.size();

	vector<double> dmuscur(NNdmu, 0.);
	for (int i = 0; i < NNdmu; ++i)
		dmuscur[i] = muStarInit[m_MapFromdMuStar[i]] - m_Chem[m_MapFromdMuStar[i]];

	VectorXd dmusold(NNdmu), dmusnew(NNdmu), dmusdelta(NNdmu);
	VectorXd fold(NNdmu), fnew(NNdmu), fdelta(NNdmu);
	for (int i = 0; i<NNdmu; ++i) {
		dmusold[i] = dmuscur[i];
	}


	MatrixXd densMatrix(NN, NN);
	VectorXd solVector(NN), xVector(NN);

	MatrixXd Jac(NNdmu, NNdmu), Jinv(NNdmu, NNdmu);
	{
		vector<double> Ps(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i) 
			Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i] + dmusold[m_MapTodMuStar[i]], 0.);
		
		vector<double> ns(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i) 
			ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i] + dmusold[m_MapTodMuStar[i]], 0.);
		
		vector<double> chi2s(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i) 
			chi2s[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_Chem[i] + dmusold[m_MapTodMuStar[i]], 0.);


		for (int i = 0; i<NN; ++i)
			for (int j = 0; j<NN; ++j) {
				densMatrix(i, j) = m_Virial[j][i] * ns[i];
				if (i == j) 
					densMatrix(i, j) += 1.;
			}

		//printf("Before n-matrix LU-decomposition\n");
		clock_t tbeg = clock();
		PartialPivLU<MatrixXd> decomp(densMatrix);
		//printf("After n-matrix LU-decomposition. Time: %lf s\n", (clock() - tbeg) / static_cast<double>(CLOCKS_PER_SEC));


		for (int i = 0; i<NN; ++i) 
			xVector[i] = ns[i];
		vector<double> np(m_densities.size(), 0.);
		solVector = decomp.solve(xVector);
		for (int i = 0; i<NN; ++i)
			np[i] = solVector[i];


		for (int k = 0; k < NNdmu; ++k) {
			double tmps = 0.;

			for (int j = 0; j < NN; ++j)
				tmps += m_Virial[j][m_MapFromdMuStar[k]] * np[j]; // Check for generality?!

			for (int i = 0; i < NN; ++i) {
				if (m_MapTodMuStar[i] == k)
					xVector[i] = chi2s[i] * m_Parameters.T * m_Parameters.T * pow(xMath::GeVtoifm(), 3) * (1. - tmps);
				else
					xVector[i] = 0.;
			}

			solVector = decomp.solve(xVector);
			for (int i = 0; i<NN; ++i)
				if (solVector[i]>1.) solVector[i] = 1.;	// Stabilizer

			for (int i = 0; i<NNdmu; ++i) {
				Jac(i, k) = 0.;
				if (i == k) 
					Jac(i, k) += 1.;
				for (int j = 0; j < NN; ++j)
					if (m_MapTodMuStar[j] == k)
						Jac(i, k) += m_Virial[m_MapFromdMuStar[i]][j] * ns[j];

				if (attrfl)
					for (int j = 0; j < NN; ++j)
						Jac(i, k) += -(m_Attr[m_MapFromdMuStar[i]][j] + m_Attr[j][m_MapFromdMuStar[i]]) * solVector[j];
			}
		}

		{
			for (int i = 0; i < fold.size(); ++i) {
				fold[i] = dmusold[i];
				for (int j = 0; j < NN; ++j)
					fold[i] += m_Virial[m_MapFromdMuStar[i]][j] * Ps[j] - (m_Attr[m_MapFromdMuStar[i]][j] + m_Attr[j][m_MapFromdMuStar[i]]) * np[j];
			}
		}
	}
	dmusnew = dmusold;
	fnew = fold;

	//printf("Before matrix inversion\n");
	clock_t tbeg = clock();
	Jinv = Jac.inverse();
	//printf("After matrix inversion. Time: %lf s\n", (clock() - tbeg) / static_cast<double>(CLOCKS_PER_SEC));

	int iters = 1;
	double maxdif = 0.;
	const int MAXITS = 200;
	while (iters<MAXITS) {
		dmusnew = dmusold - Jinv * fold;

		vector<double> Ps(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i) 
			Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i] + dmusnew[m_MapTodMuStar[i]], 0.);
		vector<double> ns(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i) 
			ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i] + dmusnew[m_MapTodMuStar[i]], 0.);

		for (int i = 0; i<NN; ++i)
			for (int j = 0; j<NN; ++j) {
				densMatrix(i, j) = m_Virial[j][i] * ns[i];
				if (i == j) 
					densMatrix(i, j) += 1.;
			}

		PartialPivLU<MatrixXd> decomp(densMatrix);

		for (int i = 0; i<NN; ++i) 
			xVector[i] = ns[i];
		vector<double> np(m_densities.size(), 0.);
		solVector = decomp.solve(xVector);
		for (int i = 0; i<NN; ++i)
			np[i] = solVector[i];

		for (int i = 0; i<fnew.size(); ++i) {
			fnew[i] = dmusnew[i];
			for (int j = 0; j<NN; ++j)
				fnew[i] += m_Virial[m_MapFromdMuStar[i]][j] * Ps[j] - (m_Attr[m_MapFromdMuStar[i]][j] + m_Attr[j][m_MapFromdMuStar[i]]) * np[j];
		}

		dmusdelta = dmusnew - dmusold;
		fdelta = fnew - fold;
		double norm = 0.;
		for (int i = 0; i<NNdmu; ++i)
			for (int j = 0; j<NNdmu; ++j)
				norm += dmusdelta[i] * Jinv(i, j) * fdelta[j];
		//TVectorD p1(NNN);
		VectorXd p1(NNdmu);
		p1 = (dmusdelta - Jinv*fdelta);
		for (int i = 0; i<NNdmu; ++i) 
			p1[i] *= 1. / norm;
		Jinv = Jinv + (p1 * dmusdelta.transpose()) * Jinv;

		maxdif = 0.;
		for (int i = 0; i<dmusnew.size(); ++i) {
				maxdif = max(maxdif, abs(dmusnew[i] - dmusold[i]));
		}
		iters++;
		//printf("Iteration %d: %E\n", iters, maxdif);

		if (maxdif<EPS) break;
		dmusold = dmusnew;
		fold = fnew;
	}
	//printf("Iteration %d: %E\n", iters, maxdif);
	for (int i = 0; i<dmuscur.size(); ++i)
		dmuscur[i] = dmusnew[i];
	if (iters == MAXITS) m_LastBroydenSuccessFlag = false;
	else m_LastBroydenSuccessFlag = true;
	m_MaxDiff = maxdif;

	vector<double> ret(NN);
	for (int i = 0; i < NN; ++i)
		ret[i] = m_Chem[i] + dmuscur[m_MapTodMuStar[i]];

	return ret;
}

std::vector<double> ThermalModelVDWFull::SearchSolutionsSingleBroydenOptimized2(const std::vector<double>& muStarInit)
{
	const double EPS = 1.0e-10;

	int NN = m_densities.size();

	bool attrfl = false;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			if (m_Attr[i][j] != 0.0) {
				attrfl = true;
				break;
			}
		}
		if (attrfl) break;
	}



	int NNdmu = m_MapFromdMuStar.size();

	vector<double> dmuscur(NNdmu, 0.);
	for (int i = 0; i < NNdmu; ++i)
		dmuscur[i] = muStarInit[m_MapFromdMuStar[i]] - m_Chem[m_MapFromdMuStar[i]];

	VectorXd dmusold(NNdmu), dmusnew(NNdmu), dmusdelta(NNdmu);
	VectorXd fold(NNdmu), fnew(NNdmu), fdelta(NNdmu);
	for (int i = 0; i<NNdmu; ++i) {
		dmusold[i] = dmuscur[i];
	}


	MatrixXd densMatrix(NNdmu, NNdmu);
	VectorXd solVector(NNdmu), xVector(NNdmu);

	MatrixXd Jac(NNdmu, NNdmu), Jinv(NNdmu, NNdmu);
	{
		vector<double> Ps(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i)
			Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i] + dmusold[m_MapTodMuStar[i]], 0.);

		vector<double> ns(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i)
			ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i] + dmusold[m_MapTodMuStar[i]], 0.);

		vector<double> chi2s(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i)
			chi2s[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_Chem[i] + dmusold[m_MapTodMuStar[i]], 0.);


		for (int i = 0; i < NNdmu; ++i) {
			for (int j = 0; j < NNdmu; ++j) {
				densMatrix(i, j) = 0.;
				if (i == j)
					densMatrix(i, j) += 1.;

				for (int m = 0; m < m_dMuStarIndices[i].size(); ++m) {
					densMatrix(i, j) += m_Virial[m_MapFromdMuStar[j]][m_dMuStarIndices[i][m]] * ns[m_dMuStarIndices[i][m]];
				}
			}
		}

		clock_t tbeg = clock();
		PartialPivLU<MatrixXd> decomp(densMatrix);
		

		for (int kp = 0; kp < NNdmu; ++kp) {
			xVector[kp] = 0.;
			for (int m = 0; m < m_dMuStarIndices[kp].size(); ++m) {
				xVector[kp] += ns[m_dMuStarIndices[kp][m]];
			}
		}

		
		solVector = decomp.solve(xVector);

		vector<double> ntildek(NNdmu, 0.);
		for (int i = 0; i < NNdmu; ++i)
			ntildek[i] = solVector[i];

		vector<double> np(m_densities.size(), 0.);
		for (int i = 0; i < NN; ++i) {
			np[i] = 0.;
			for (int k = 0; k < NNdmu; ++k) {
				np[i] += m_Virial[m_MapFromdMuStar[k]][i] * solVector[k];
			}
			np[i] = (1. - np[i]) * ns[i];
		}

		for (int kp = 0; kp < NNdmu; ++kp) {

			if (attrfl) {
				for (int l = 0; l < NNdmu; ++l) {
					xVector[l] = 0.;
					for (int m = 0; m < m_dMuStarIndices[l].size(); ++m) {
						int ti = m_dMuStarIndices[l][m];
						if (m_MapTodMuStar[ti] != kp)
							continue;

						double tmps = 0.;
						for (int k = 0; k < NNdmu; ++k) {
							tmps += m_Virial[m_MapFromdMuStar[k]][ti] * ntildek[k];
						}
						xVector[l] += chi2s[ti] * m_Parameters.T * m_Parameters.T * pow(xMath::GeVtoifm(), 3) * (1. - tmps);
					}
				}

				solVector = decomp.solve(xVector);
				for (int i = 0; i < NNdmu; ++i)
					if (solVector[i] > 1.) solVector[i] = 1.;	// Stabilizer
			}

			std::vector<double> dnjdmukp(NN, 0.);
			if (attrfl) {
				for (int j = 0; j < NN; ++j) {
					for (int kk = 0; kk < NNdmu; ++kk) {
						dnjdmukp[j] += -m_Virial[m_MapFromdMuStar[kk]][j] * solVector[kk] * ns[j];
					}

					if (m_MapTodMuStar[j] == kp) {
						double tmps = 0.;
						for (int kk = 0; kk < NNdmu; ++kk) {
							tmps += m_Virial[m_MapFromdMuStar[kk]][j] * ntildek[kk];
						}
						dnjdmukp[j] += (1. - tmps) * chi2s[j] * m_Parameters.T * m_Parameters.T * pow(xMath::GeVtoifm(), 3);
					}
				}
			}


			for (int k = 0; k < NNdmu; ++k) {
				Jac(k, kp) = 0.;
				if (k == kp)
					Jac(k, kp) += 1.;
				for (int m = 0; m < m_dMuStarIndices[kp].size(); ++m) { 
					int tj = m_dMuStarIndices[kp][m];
					Jac(k, kp) += m_Virial[m_MapFromdMuStar[k]][tj] * ns[tj];
				}

				if (attrfl) {
					for (int j = 0; j < NN; ++j) {
						Jac(k, kp) += -(m_Attr[m_MapFromdMuStar[k]][j] + m_Attr[j][m_MapFromdMuStar[k]]) * dnjdmukp[j];
					}
				}
			}

		}

		{
			for (int i = 0; i < fold.size(); ++i) {
				fold[i] = dmusold[i];
				for (int j = 0; j < NN; ++j)
					fold[i] += m_Virial[m_MapFromdMuStar[i]][j] * Ps[j] - (m_Attr[m_MapFromdMuStar[i]][j] + m_Attr[j][m_MapFromdMuStar[i]]) * np[j];
			}
		}
	}
	dmusnew = dmusold;
	fnew = fold;

	//printf("Before matrix inversion\n");
	clock_t tbeg = clock();
	Jinv = Jac.inverse();
	//printf("After matrix inversion. Time: %lf s\n", (clock() - tbeg) / static_cast<double>(CLOCKS_PER_SEC));

	int iters = 1;
	double maxdif = 0.;
	const int MAXITS = 200;
	while (iters<MAXITS) {
		dmusnew = dmusold - Jinv * fold;

		vector<double> Ps(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i)
			Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i] + dmusnew[m_MapTodMuStar[i]], 0.);
		vector<double> ns(m_densities.size(), 0.);
		for (int i = 0; i<NN; ++i)
			ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i] + dmusnew[m_MapTodMuStar[i]], 0.);

		for (int i = 0; i < NNdmu; ++i) {
			for (int j = 0; j < NNdmu; ++j) {
				densMatrix(i, j) = 0.;
				if (i == j)
					densMatrix(i, j) += 1.;

				for (int m = 0; m < m_dMuStarIndices[i].size(); ++m) {
					densMatrix(i, j) += m_Virial[m_MapFromdMuStar[j]][m_dMuStarIndices[i][m]] * ns[m_dMuStarIndices[i][m]];
				}
			}
		}

		PartialPivLU<MatrixXd> decomp(densMatrix);

		for (int kp = 0; kp < NNdmu; ++kp) {
			xVector[kp] = 0.;
			for (int m = 0; m < m_dMuStarIndices[kp].size(); ++m) {
				xVector[kp] += ns[m_dMuStarIndices[kp][m]];
			}
		}

		
		solVector = decomp.solve(xVector);

		vector<double> ntildek(NNdmu, 0.);
		for (int i = 0; i < NNdmu; ++i)
			ntildek[i] = solVector[i];

		vector<double> np(m_densities.size(), 0.);
		for (int i = 0; i < NN; ++i) {
			np[i] = 0.;
			for (int k = 0; k < NNdmu; ++k) {
				np[i] += m_Virial[m_MapFromdMuStar[k]][i] * solVector[k];
			}
			np[i] = (1. - np[i]) * ns[i];
		}

		for (int i = 0; i<fnew.size(); ++i) {
			fnew[i] = dmusnew[i];
			for (int j = 0; j<NN; ++j)
				fnew[i] += m_Virial[m_MapFromdMuStar[i]][j] * Ps[j] - (m_Attr[m_MapFromdMuStar[i]][j] + m_Attr[j][m_MapFromdMuStar[i]]) * np[j];
		}

		dmusdelta = dmusnew - dmusold;
		fdelta = fnew - fold;
		double norm = 0.;
		for (int i = 0; i<NNdmu; ++i)
			for (int j = 0; j<NNdmu; ++j)
				norm += dmusdelta[i] * Jinv(i, j) * fdelta[j];
		VectorXd p1(NNdmu);
		p1 = (dmusdelta - Jinv*fdelta);
		for (int i = 0; i<NNdmu; ++i)
			p1[i] *= 1. / norm;
		Jinv = Jinv + (p1 * dmusdelta.transpose()) * Jinv;

		maxdif = 0.;
		for (int i = 0; i<dmusnew.size(); ++i) {
			maxdif = max(maxdif, abs(dmusnew[i] - dmusold[i]));
		}
		iters++;
		if (maxdif > 1.e1)
			printf("Iteration %d: %E\n", iters, maxdif);

		if (maxdif<EPS) break;
		dmusold = dmusnew;
		fold = fnew;
	}
	//printf("Iteration %d: %E\n", iters, maxdif);
	for (int i = 0; i<dmuscur.size(); ++i)
		dmuscur[i] = dmusnew[i];
	if (iters == MAXITS) m_LastBroydenSuccessFlag = false;
	else m_LastBroydenSuccessFlag = true;
	m_MaxDiff = maxdif;

	vector<double> ret(NN);
	for (int i = 0; i < NN; ++i)
		ret[i] = m_Chem[i] + dmuscur[m_MapTodMuStar[i]];

	return ret;
}

vector<double> ThermalModelVDWFull::SearchSolutionsMulti(int iters) {
	vector<double> csol(m_densities.size(), 0.);
	double Psol = 0.;
	bool solved = false;
	double muBmin = m_Parameters.muB - 0.5 * xMath::mnucleon();
	double muBmax = m_Parameters.muB + 0.5 * xMath::mnucleon();
	double dmu = (muBmax - muBmin) / iters;
	vector<double> curmust(m_densities.size(), 0.);
	double maxdif = 0.;
	for(int i=0; i<iters; ++i) {
		double tmu = muBmin + (0.5 + i) * dmu;
		for(int j=0; j<curmust.size(); ++j) {
			curmust[j] = m_Chem[j] + (tmu - m_Parameters.muB) * m_Chem[j] / m_Parameters.muB;
			if (m_TPS->Particles()[j].Statistics()==-1 && curmust[j] > m_TPS->Particles()[j].Mass()) 
				curmust[j] = 0.98 * m_TPS->Particles()[j].Mass();
		}

		vector<double> sol = SearchSolutionsSingleBroyden(curmust);
		bool fl = true;
		for(int i=0; i<sol.size(); ++i)
			if (sol[i]!=sol[i]) fl = false;
		fl &= m_LastBroydenSuccessFlag;
		if (!fl) continue;

		for(int i=0;i<m_TPS->Particles().size();++i) 
			m_DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, sol[i], 0.);

		int NN = m_densities.size();

		MatrixXd densMatrix(NN, NN);
		VectorXd solVector(NN), xVector(NN);

		for(int i=0;i<NN;++i)
			for(int j=0;j<NN;++j) {
				densMatrix(i,j) = m_Virial[j][i] * m_DensitiesId[i];
				if (i==j) densMatrix(i,j) += 1.;
			}

		PartialPivLU<MatrixXd> decomp(densMatrix);

		for(int i=0;i<NN;++i) xVector[i] = m_DensitiesId[i];
		solVector = decomp.solve(xVector);
		for(int i=0;i<NN;++i) m_densities[i] = solVector[i];

		double tP = 0.;
		for(int i=0;i<m_densities.size();++i) tP += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, sol[i], 0.);
		for(int i=0;i<m_densities.size();++i)
			for(int j=0;j<m_densities.size();++j)
				tP += -m_Attr[i][j] * m_densities[i] * m_densities[j];

		if (!solved || tP > Psol) {
			solved = true;
			Psol = tP;
			csol = sol;
			maxdif = m_MaxDiff;
		}
	}
	m_LastBroydenSuccessFlag = solved;
	m_MaxDiff = maxdif;
	return csol;
}

void ThermalModelVDWFull::SearchSolutions() {
	if (!m_SearchMultipleSolutions) {

		vector<double> muStarInit = m_MuStar;

		for(int i=0;i<muStarInit.size();++i) {
			if (m_TPS->Particles()[i].Statistics()==-1 && muStarInit[i] > m_TPS->Particles()[i].Mass()) 
				muStarInit[i] = 0.98 * m_TPS->Particles()[i].Mass();
		}

		m_MuStar = SearchSolutionsSingleBroydenOptimized2(muStarInit);

	}
	else {
		m_MuStar = SearchSolutionsMulti(100);
	}
}

void ThermalModelVDWFull::CalculateDensities() {
	CalculateDensitiesNew();
}

void ThermalModelVDWFull::CalculateDensitiesOld() {
	m_FluctuationsCalculated = false;

	std::map< std::vector<double> , int> m_MapVDW;

	int NN = m_densities.size();

	{
		m_MapTodMuStar.resize(NN);
		m_MapFromdMuStar.clear();
		m_MapVDW.clear();
		m_dMuStarIndices.clear();

		int tind = 0;
		for (int i = 0; i < NN; ++i) {
			std::vector<double> VDWParam(0);
			for (int j = 0; j < NN; ++j) {
				VDWParam.push_back(m_Virial[i][j]);
			}
			for (int j = 0; j < NN; ++j) {
				VDWParam.push_back(m_Attr[i][j] + m_Attr[j][i]);
			}

			if (m_MapVDW.count(VDWParam) == 0) {
				m_MapVDW[VDWParam] = tind;
				m_MapTodMuStar[i]  = tind;
				m_MapFromdMuStar.push_back(i);
				m_dMuStarIndices.push_back(std::vector<int>(1, i));
				tind++;
			}
			else {
				m_MapTodMuStar[i] = m_MapVDW[VDWParam];
				m_dMuStarIndices[m_MapVDW[VDWParam]].push_back(i);
			}
		}

		printf("Optimization: %d --> %d\n", NN, m_MapFromdMuStar.size());
	}

	clock_t tbeg = clock();

	{

		vector<double> muStarInit = m_MuStar;

		for (int i = 0; i<muStarInit.size(); ++i) {
			if (m_TPS->Particles()[i].Statistics() == -1 && muStarInit[i] > m_TPS->Particles()[i].Mass())
				muStarInit[i] = 0.98 * m_TPS->Particles()[i].Mass();
		}


		m_MuStar = SearchSolutionsSingleBroydenOptimized(muStarInit);
	}
	

	printf("Solution time = %lf ms\n", (clock() - tbeg) / (double)(CLOCKS_PER_SEC) * 1.e3);

	tbeg = clock();

	for(int i=0;i<m_TPS->Particles().size();++i) 
		m_DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_MuStar[i], 0.);

	MatrixXd densMatrix(NN, NN);
	VectorXd solVector(NN), xVector(NN);

	for(int i=0;i<NN;++i)
		for(int j=0;j<NN;++j) {
			densMatrix(i,j) = m_Virial[j][i] * m_DensitiesId[i];
			if (i==j) densMatrix(i,j) += 1.;
		}

	PartialPivLU<MatrixXd> decomp(densMatrix);

	for(int i=0;i<NN;++i) xVector[i] = m_DensitiesId[i];
	solVector = decomp.solve(xVector);
	for(int i=0;i<NN;++i) m_densities[i] = solVector[i];

	// TODO: Scalar density properly
	for(int i=0;i<NN;++i) xVector[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ScalarDensity, m_UseWidth, m_MuStar[i], 0.);
	solVector = decomp.solve(xVector);
	m_scaldens.resize(m_densities.size());
	for(int i=0;i<NN;++i) m_scaldens[i] = solVector[i];


	tbeg = clock();

	CalculateFeeddown();

	m_Calculated = true;
}

void ThermalModelVDWFull::CalculateDensitiesNew() {
	m_FluctuationsCalculated = false;

	std::map< std::vector<double>, int> m_MapVDW;

	int NN = m_densities.size();

	{
		m_MapTodMuStar.resize(NN);
		m_MapFromdMuStar.clear();
		m_MapVDW.clear();
		m_dMuStarIndices.clear();

		int tind = 0;
		for (int i = 0; i < NN; ++i) {
			std::vector<double> VDWParam(0);
			for (int j = 0; j < NN; ++j) {
				VDWParam.push_back(m_Virial[i][j]);
			}
			for (int j = 0; j < NN; ++j) {
				VDWParam.push_back(m_Attr[i][j] + m_Attr[j][i]);
			}

			if (m_MapVDW.count(VDWParam) == 0) {
				m_MapVDW[VDWParam] = tind;
				m_MapTodMuStar[i] = tind;
				m_MapFromdMuStar.push_back(i);
				m_dMuStarIndices.push_back(std::vector<int>(1, i));
				tind++;
			}
			else {
				m_MapTodMuStar[i] = m_MapVDW[VDWParam];
				m_dMuStarIndices[m_MapVDW[VDWParam]].push_back(i);
			}
		}
	}

	clock_t tbeg = clock();

	SearchSolutions();

	tbeg = clock();

	for (int i = 0; i<m_TPS->Particles().size(); ++i)
		m_DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_MuStar[i], 0.);

	
	int NNdmu = m_MapFromdMuStar.size();

	MatrixXd densMatrix(NNdmu, NNdmu);
	VectorXd solVector(NNdmu), xVector(NNdmu);

	for (int i = 0; i < NNdmu; ++i) {
		for (int j = 0; j < NNdmu; ++j) {
			densMatrix(i, j) = 0.;
			if (i == j)
				densMatrix(i, j) += 1.;

			for (int m = 0; m < m_dMuStarIndices[i].size(); ++m) {
				densMatrix(i, j) += m_Virial[m_MapFromdMuStar[j]][m_dMuStarIndices[i][m]] * m_DensitiesId[m_dMuStarIndices[i][m]];
			}
		}
	}

	PartialPivLU<MatrixXd> decomp(densMatrix);

	for (int kp = 0; kp < NNdmu; ++kp) {
		xVector[kp] = 0.;
		for (int m = 0; m < m_dMuStarIndices[kp].size(); ++m) {
			xVector[kp] += m_DensitiesId[m_dMuStarIndices[kp][m]];
		}
	}

	solVector = decomp.solve(xVector);

	vector<double> ntildek(NNdmu, 0.);
	for (int i = 0; i < NNdmu; ++i)
		ntildek[i] = solVector[i];

	//vector<double> np(m_densities.size(), 0.);
	for (int i = 0; i < NN; ++i) {
		m_densities[i] = 0.;
		for (int k = 0; k < NNdmu; ++k) {
			m_densities[i] += m_Virial[m_MapFromdMuStar[k]][i] * solVector[k];
		}
		m_densities[i] = (1. - m_densities[i]) * m_DensitiesId[i];
	}

	// TODO: Scalar density properly
	m_scaldens = m_densities;

	tbeg = clock();

	CalculateFeeddown();

	m_Calculated = true;
}

vector<double> ThermalModelVDWFull::CalculateChargeFluctuations(const std::vector<double> &chgs, int order) {
	vector<double> ret(order + 1, 0.);
	
	// chi1
	for(int i=0;i<m_densities.size();++i)
		ret[0] += chgs[i] * m_densities[i];

	ret[0] /= pow(m_Parameters.T * xMath::GeVtoifm(), 3);

	if (order<2) return ret;
	// Preparing matrix for system of linear equations
	int NN = m_densities.size();
	MatrixXd densMatrix(2*NN, 2*NN);
	VectorXd solVector(2*NN), xVector(2*NN);

	vector<double> chi2id(m_densities.size());
	for(int i=0;i<NN;++i) 
		chi2id[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_MuStar[i], 0.);

	for(int i=0;i<NN;++i)
		for(int j=0;j<NN;++j) {
			densMatrix(i,j) = m_Virial[j][i] * m_DensitiesId[i];
			if (i==j) densMatrix(i,j) += 1.;
		}

	for(int i=0;i<NN;++i)
		for(int j=0;j<NN;++j) 
			densMatrix(i,NN+j) = 0.;

	for(int i=0;i<NN;++i) {
		densMatrix(i,NN+i) = 0.;
		for(int k=0;k<NN;++k) {
			densMatrix(i,NN+i) += m_Virial[k][i] * m_densities[k];
		}
		densMatrix(i,NN+i) = (densMatrix(i,NN+i) - 1.) * chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T;
	}

	for(int i=0;i<NN;++i)
		for(int j=0;j<NN;++j) {
			densMatrix(NN+i,j) = -(m_Attr[i][j] + m_Attr[j][i]);
		}
	
	for(int i=0;i<NN;++i)
		for(int j=0;j<NN;++j) {
			densMatrix(NN+i,NN+j) = m_Virial[i][j] * m_DensitiesId[j];
			if (i==j) densMatrix(NN+i,NN+j) += 1.;
		}


	PartialPivLU<MatrixXd> decomp(densMatrix);

	// chi2
	vector<double> dni(NN, 0.), dmus(NN, 0.);

	for(int i=0;i<NN;++i) {
		xVector[i]    = 0.;
		xVector[NN+i] = chgs[i];
	}

	solVector = decomp.solve(xVector);

	for(int i=0;i<NN;++i) {
		dni[i]  = solVector[i];
		dmus[i] = solVector[NN+i];
	}

	for(int i=0;i<NN;++i)
		ret[1] += chgs[i] * dni[i];

	ret[1] /= pow(m_Parameters.T, 2) * pow(xMath::GeVtoifm(), 3);

	if (order<3) return ret;
	// chi3
	vector<double> d2ni(NN, 0.), d2mus(NN, 0.);

	vector<double> chi3id(m_densities.size());
	for(int i=0;i<NN;++i) 
		chi3id[i] = m_TPS->Particles()[i].chi(3, m_Parameters, m_UseWidth, m_MuStar[i], 0.);

	for(int i=0;i<NN;++i) {
		xVector[i]    = 0.;

		double tmp = 0.;
		for(int j=0;j<NN;++j) tmp += m_Virial[j][i] * dni[j];
		tmp = -2. * tmp * chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[i];
		xVector[i] += tmp;

		tmp = 0.;
		for(int j=0;j<NN;++j) tmp += m_Virial[j][i] * m_densities[j];
		tmp = -(tmp - 1.) * chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * dmus[i] * dmus[i];
		xVector[i] += tmp;
	}
	for(int i=0;i<NN;++i) {
		xVector[NN+i]    = 0.;

		double tmp = 0.;
		for(int j=0;j<NN;++j) tmp += -m_Virial[i][j] * dmus[j] * chi2id[j] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[j];
		
		xVector[NN+i] = tmp;
	}

	solVector = decomp.solve(xVector);

	for(int i=0;i<NN;++i) {
		d2ni[i]  = solVector[i];
		d2mus[i] = solVector[NN+i];
	}

	for(int i=0;i<NN;++i)
		ret[2] += chgs[i] * d2ni[i];

	ret[2] /= m_Parameters.T * pow(xMath::GeVtoifm(), 3);


	if (order<4) return ret;

	// chi4
	vector<double> d3ni(NN, 0.), d3mus(NN, 0.);

	vector<double> chi4id(m_densities.size());
	for(int i=0;i<NN;++i) 
		chi4id[i] = m_TPS->Particles()[i].chi(4, m_Parameters, m_UseWidth, m_MuStar[i], 0.);

	vector<double> dnis(NN, 0.);
	for(int i=0;i<NN;++i) {
		dnis[i] = chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[i];
	}

	vector<double> d2nis(NN, 0.);
	for(int i=0;i<NN;++i) {
		d2nis[i] = chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * dmus[i] * dmus[i] + 
			chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * d2mus[i];
	}

	for(int i=0;i<NN;++i) {
		xVector[i]    = 0.;

		double tmp = 0.;
		for(int j=0;j<NN;++j) tmp += m_Virial[j][i] * dni[j];
		tmp = -3. * tmp * d2nis[i];
		xVector[i] += tmp;

		tmp = 0.;
		for(int j=0;j<NN;++j) tmp += m_Virial[j][i] * d2ni[j];
		tmp = -3. * tmp * dnis[i];
		xVector[i] += tmp;

		double tmps = 0.;
		for(int j=0;j<NN;++j) tmps += m_Virial[j][i] * m_densities[j];

		tmp = -(tmps - 1.) * chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * d2mus[i] * 3. * dmus[i];
		xVector[i] += tmp;

		tmp = -(tmps - 1.) * chi4id[i] * pow(xMath::GeVtoifm(), 3) * dmus[i] * dmus[i] * dmus[i];
		xVector[i] += tmp;
	}
	for(int i=0;i<NN;++i) {
		xVector[NN+i]    = 0.;

		double tmp = 0.;
		for(int j=0;j<NN;++j) tmp += -2. * m_Virial[i][j] * d2mus[j] * dnis[j];
		xVector[NN+i] += tmp;

		tmp = 0.;
		for(int j=0;j<NN;++j) tmp += -m_Virial[i][j] * dmus[j] * d2nis[j];
		xVector[NN+i] += tmp;
	}

	solVector = decomp.solve(xVector);

	for(int i=0;i<NN;++i) {
		d3ni[i]  = solVector[i];
		d3mus[i] = solVector[NN+i];
	}

	for(int i=0;i<NN;++i)
		ret[3] += chgs[i] * d3ni[i];

	ret[3] /= pow(xMath::GeVtoifm(), 3);

	return ret;
}

// TODO include correlations
std::vector< std::vector<double> >  ThermalModelVDWFull::CalculateFluctuations(int order) {
	if (order<1) return m_chi;
	
	vector<double> chgs(m_densities.size());
	vector<double> chis;

	// Baryon charge
	for(int i=0;i<chgs.size();++i)
		chgs[i] = m_TPS->Particles()[i].BaryonCharge();
	chis = CalculateChargeFluctuations(chgs, order);
	for(int i=0;i<order;++i) m_chi[i][0] = chis[i];

	// Electric charge
	for(int i=0;i<chgs.size();++i)
		chgs[i] = m_TPS->Particles()[i].ElectricCharge();
	chis = CalculateChargeFluctuations(chgs, order);
	for(int i=0;i<order;++i) m_chi[i][1] = chis[i];

	// Strangeness charge
	for(int i=0;i<chgs.size();++i)
		chgs[i] = m_TPS->Particles()[i].Strangeness();
	chis = CalculateChargeFluctuations(chgs, order);
	for(int i=0;i<order;++i) m_chi[i][2] = chis[i];

	// Arbitrary charge
	for(int i=0;i<chgs.size();++i)
		chgs[i] = m_TPS->Particles()[i].ArbitraryCharge();
	chis = CalculateChargeFluctuations(chgs, order);
	for(int i=0;i<order;++i) m_chiarb[i] = chis[i];

	return m_chi;
}

void ThermalModelVDWFull::CalculateTwoParticleCorrelations()
{
	int NN = m_densities.size();

	m_PrimCorrel.resize(NN);
	for (int i = 0; i < NN; ++i)
		m_PrimCorrel[i].resize(NN);
	m_TotalCorrel = m_PrimCorrel;

	MatrixXd densMatrix(2 * NN, 2 * NN);
	VectorXd solVector(2 * NN), xVector(2 * NN);

	vector<double> chi2id(m_densities.size());
	for (int i = 0; i<NN; ++i)
		chi2id[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_MuStar[i], 0.);

	for (int i = 0; i<NN; ++i)
		for (int j = 0; j<NN; ++j) {
			densMatrix(i, j) = m_Virial[j][i] * m_DensitiesId[i];
			if (i == j) densMatrix(i, j) += 1.;
		}

	for (int i = 0; i<NN; ++i)
		for (int j = 0; j<NN; ++j)
			densMatrix(i, NN + j) = 0.;

	for (int i = 0; i<NN; ++i) {
		densMatrix(i, NN + i) = 0.;
		for (int k = 0; k<NN; ++k) {
			densMatrix(i, NN + i) += m_Virial[k][i] * m_densities[k];
		}
		densMatrix(i, NN + i) = (densMatrix(i, NN + i) - 1.) * chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T;
	}

	for (int i = 0; i<NN; ++i)
		for (int j = 0; j<NN; ++j) {
			densMatrix(NN + i, j) = -(m_Attr[i][j] + m_Attr[j][i]);
		}

	for (int i = 0; i<NN; ++i)
		for (int j = 0; j<NN; ++j) {
			densMatrix(NN + i, NN + j) = m_Virial[i][j] * m_DensitiesId[j];
			if (i == j) densMatrix(NN + i, NN + j) += 1.;
		}

	PartialPivLU<MatrixXd> decomp(densMatrix);

	for (int k = 0; k < NN; ++k) {
		vector<double> dni(NN, 0.), dmus(NN, 0.);

		for (int i = 0; i < NN; ++i) {
			xVector[i] = 0.;
			xVector[NN + i] = static_cast<int>(i == k);
		}

		solVector = decomp.solve(xVector);

		for (int i = 0; i < NN; ++i) {
			dni[i]  = solVector[i];
			dmus[i] = solVector[NN + i];
		}

		for (int j = 0; j < NN; ++j) {
			m_PrimCorrel[j][k] = dni[j];
		}
	}

	for (int i = 0; i < NN; ++i) {
		m_wprim[i] = m_PrimCorrel[i][i];
		if (m_densities[i] > 0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
		else m_wprim[i] = 1.;
	}

	CalculateSusceptibilityMatrix();
	CalculateTwoParticleFluctuationsDecays();
	CalculateProxySusceptibilityMatrix();
}

void ThermalModelVDWFull::CalculateFluctuations()
{
	CalculateTwoParticleCorrelations();

	for (int i = 0; i < m_wprim.size(); ++i) {
		m_skewprim[i] = 1.;
		m_kurtprim[i] = 1.;
		m_skewtot[i] = 1.;
		m_kurttot[i] = 1.;
	}

	m_FluctuationsCalculated = true;
}


double ThermalModelVDWFull::CalculateEnergyDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for(int i=0;i<m_densities.size();++i) 
		if (m_densities[i]>0.) 
			ret += m_densities[i] / m_DensitiesId[i] * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_MuStar[i], 0.);
	for(int i=0;i<m_densities.size();++i)
		for(int j=0;j<m_densities.size();++j)
			ret += -m_Attr[i][j] * m_densities[i] * m_densities[j];

	if (m_TemperatureDependentAB) {
		for (int i = 0; i < m_densities.size(); ++i) {
			double tPid = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i], 0.);
			for (int j = 0; j < m_densities.size(); ++j) {
				ret += -tPid * m_densities[j] * m_Parameters.T * m_VirialdT[j][i];
				ret += m_Parameters.T * m_AttrdT[i][j] * m_densities[i] * m_densities[j];
			}
		}
	}

	return ret;
}

double ThermalModelVDWFull::CalculateEntropyDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for(int i=0;i<m_densities.size();++i) 
		if (m_densities[i]>0.) 
			ret += m_densities[i] / m_DensitiesId[i] * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_MuStar[i], 0.);
	
	if (m_TemperatureDependentAB) {
		for (int i = 0; i < m_densities.size(); ++i) {
			double tPid = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i], 0.);
			for (int j = 0; j < m_densities.size(); ++j) {
				ret += -tPid * m_densities[j] * m_VirialdT[j][i];
				ret += m_AttrdT[i][j] * m_densities[i] * m_densities[j];
			}
		}
	}
	
	return ret;
}

// Dummy
double ThermalModelVDWFull::CalculateBaryonMatterEntropyDensity() {
    double ret = 0.;
    return ret;
}
double ThermalModelVDWFull::CalculateMesonMatterEntropyDensity() {
    double ret = 0.;
    return ret;
}

double ThermalModelVDWFull::CalculatePressure() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for(int i=0;i<m_densities.size();++i) ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i], 0.);
	for(int i=0;i<m_densities.size();++i)
		for(int j=0;j<m_densities.size();++j)
			ret += -m_Attr[i][j] * m_densities[i] * m_densities[j];
	return ret;
}


double ThermalModelVDWFull::ParticleScalarDensity(int part) {
	if (!m_Calculated) CalculateDensities();

	return m_scaldens[part];
}

double ThermalModelVDWFull::MuShift(int id)
{
	if (id >= 0. && id < m_Virial.size())
		return m_MuStar[id] - m_Chem[id];
	else
		return 0.0;
}

double ThermalModelVDWFull::VirialCoefficient(int i, int j) const
{
	if (i<0 || i >= m_Virial.size() || j<0 || j>=m_Virial.size())
		return 0.;
	return m_Virial[i][j];
}

double ThermalModelVDWFull::AttractionCoefficient(int i, int j) const
{
	if (i<0 || i >= m_Attr.size() || j<0 || j>=m_Attr.size())
		return 0.;
	return m_Attr[i][j];
}

double ThermalModelVDWFull::VirialCoefficientdT(int i, int j) const
{
	if (i<0 || i >= m_VirialdT.size() || j<0 || j>=m_VirialdT.size())
		return 0.;
	return m_VirialdT[i][j];
}

double ThermalModelVDWFull::AttractionCoefficientdT(int i, int j) const
{
	if (i<0 || i >= m_AttrdT.size() || j<0 || j>=m_AttrdT.size())
		return 0.;
	return m_AttrdT[i][j];
}