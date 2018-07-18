#include "HRGBase/ThermalModelBase.h"

#include <cstdio>
#include <algorithm>

#include <Eigen/Dense>

#include "HRGBase/ThermalParticleSystem.h"

using namespace Eigen;

using namespace std;


namespace ThermalModelBaseNamespace {
	ThermalModelBase *gThM;

	void broyden2(vector<double> &xin, vector<double>(*func)(const vector<double>&, ThermalModelBase*), ThermalModelBase *ThM) {
		const double TOLF = 1.0e-8, EPS = 1.0e-8;
		const int MAXITS = 200;
		bool fl = 0;
		if (abs(xin[0]) < 1e-8 && abs(xin[1]) < 1e-8) fl = 1;
		vector<double> h = xin, xh1 = xin, xh2 = xin;
		for (int i = 0; i < xin.size(); ++i) {
			h[i] = EPS*abs(h[i]);
			if (h[i] == 0.0) h[i] = EPS;
			h[i] = max(EPS, h[i]);
			if (i == 0) xh1[i] = xin[i] + h[i];
			else xh2[i] = xin[i] + h[i];
		}
		vector<double> r1 = func(xin, ThM), r21 = func(xh1, ThM), r22 = func(xh2, ThM);
		if (abs(r1[0]) < TOLF && abs(r1[1]) < TOLF) return;
		double J[2][2];
		J[0][0] = (r21[0] - r1[0]) / h[0];
		J[0][1] = (r22[0] - r1[0]) / h[1];
		J[1][0] = (r21[1] - r1[1]) / h[0];
		J[1][1] = (r22[1] - r1[1]) / h[1];
		double Jinv[2][2];
		double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
		if (det == 0.0) {
			printf("**WARNING** singular Jacobian in broyden2");
			return;
		}
		Jinv[0][0] = J[1][1] / det;
		Jinv[0][1] = -J[0][1] / det;
		Jinv[1][0] = -J[1][0] / det;
		Jinv[1][1] = J[0][0] / det;
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				if (J[i][j] > 1.e8) return;
		vector<double> xold = xin;
		vector<double> rold = r1;
		double Jinvold[2][2];
		vector<double> rprevten = r1;
		for (int i1 = 0; i1 < 2; ++i1)
			for (int i2 = 0; i2 < 2; ++i2) Jinvold[i1][i2] = Jinv[i1][i2];
		for (int iter = 1; iter <= MAXITS; ++iter) {
			xin[0] = xold[0] - Jinv[0][0] * rold[0] - Jinv[0][1] * rold[1];
			xin[1] = xold[1] - Jinv[1][0] * rold[0] - Jinv[1][1] * rold[1];
			r1 = func(xin, ThM);
			if (abs(r1[0]) < TOLF && abs(r1[1]) < TOLF) {
				return;
			}
			double JF[2];
			double DF[2];
			double dx[2];
			DF[0] = (r1[0] - rold[0]);
			DF[1] = (r1[1] - rold[1]);
			JF[0] = Jinv[0][0] * DF[0] + Jinv[0][1] * DF[1];
			JF[1] = Jinv[1][0] * DF[0] + Jinv[1][1] * DF[1];
			dx[0] = xin[0] - xold[0];
			dx[1] = xin[1] - xold[1];
			double znam = dx[0] * JF[0] + dx[1] * JF[1];
			if (znam == 0.0) {
				printf("**WARNING** singular Jacobian in broyden2");
				return;
			}
			double xJ[2];
			JF[0] = dx[0] - JF[0];
			JF[1] = dx[1] - JF[1];
			xJ[0] = dx[0] * Jinv[0][0] + dx[1] * Jinv[1][0];
			xJ[1] = dx[0] * Jinv[0][1] + dx[1] * Jinv[1][1];
			Jinv[0][0] = Jinv[0][0] + JF[0] * xJ[0] / znam;
			Jinv[0][1] = Jinv[0][1] + JF[0] * xJ[1] / znam;
			Jinv[1][0] = Jinv[1][0] + JF[1] * xJ[0] / znam;
			Jinv[1][1] = Jinv[1][1] + JF[1] * xJ[1] / znam;
			xold = xin;
			rold = r1;
		}
		printf("**WARNING** reached maximum number of interations in broyden");
	}

	vector<double> function2(const vector<double> &xin, ThermalModelBase *ThM) {
		vector<double> ret(2);
		ThM->SetElectricChemicalPotential(xin[0]);
		ThM->SetStrangenessChemicalPotential(xin[1]);
		ThM->FillChemicalPotentials();
		ThM->CalculateDensities();
		double fBd = ThM->CalculateBaryonDensity();
		double fCd = ThM->CalculateChargeDensity();
		double fSd = ThM->CalculateStrangenessDensity();
		double fASd = ThM->CalculateAbsoluteStrangenessDensity();
		ret[0] = (fCd / fBd - ThM->QoverB()) / ThM->QoverB();
		ret[1] = fSd / fASd;
		return ret;
	}

	VectorXd broydenEigenFunc(const VectorXd &xin, ThermalModelBase *ThM) {
		int i1 = 0;
		if (ThM->ConstrainMuQ()) { ThM->SetElectricChemicalPotential(xin[i1]); i1++; }
		if (ThM->ConstrainMuS()) { ThM->SetStrangenessChemicalPotential(xin[i1]); i1++; }
		if (ThM->ConstrainMuC()) { ThM->SetCharmChemicalPotential(xin[i1]); i1++; }
		ThM->FillChemicalPotentials();
		ThM->CalculateDensities();

		double fBd		= ThM->CalculateBaryonDensity();
		double fQd		= ThM->CalculateChargeDensity();
		double fSd		= ThM->CalculateStrangenessDensity();
		double fASd = ThM->CalculateAbsoluteStrangenessDensity();
		double fCd		= ThM->CalculateCharmDensity();
		double fACd = ThM->CalculateAbsoluteCharmDensity();

		int NNN = 0;
		if (ThM->ConstrainMuQ()) NNN++;
		if (ThM->ConstrainMuS()) NNN++;
		if (ThM->ConstrainMuC()) NNN++;
		VectorXd ret(NNN);

		i1 = 0;
		// Charge derivatives
		if (ThM->ConstrainMuQ()) {
			ret[i1] = (fQd / fBd - ThM->QoverB()) / ThM->QoverB();

			i1++;
		}


		// Strangeness derivatives
		if (ThM->ConstrainMuS()) {
			ret[i1] = fSd / fASd;

			i1++;
		}


		// Charm derivatives
		if (ThM->ConstrainMuC()) {
			ret[i1] = fCd / fACd;

			i1++;
		}

		return ret;
	}

	MatrixXd broydenEigenJacobian(const VectorXd &xin, ThermalModelBase *ThM) {
		int i1 = 0;
		if (ThM->ConstrainMuQ()) { ThM->SetElectricChemicalPotential(xin[i1]); i1++; }
		if (ThM->ConstrainMuS()) { ThM->SetStrangenessChemicalPotential(xin[i1]); i1++; }
		if (ThM->ConstrainMuC()) { ThM->SetCharmChemicalPotential(xin[i1]); i1++; }
		ThM->FillChemicalPotentials();
		ThM->CalculateDensities();

		double fBd = ThM->CalculateBaryonDensity();
		double fQd = ThM->CalculateChargeDensity();
		double fSd = ThM->CalculateStrangenessDensity();
		double fASd = ThM->CalculateAbsoluteStrangenessDensity();
		double fCd = ThM->CalculateCharmDensity();
		double fACd = ThM->CalculateAbsoluteCharmDensity();

		vector<double> m_wprim;
		m_wprim.resize(ThM->Densities().size());
		for (int i = 0; i < m_wprim.size(); ++i) 
			m_wprim[i] = ThM->CalculateParticleScaledVariance(i);

		int NNN = 0;
		if (ThM->ConstrainMuQ()) NNN++;
		if (ThM->ConstrainMuS()) NNN++;
		if (ThM->ConstrainMuC()) NNN++;
		MatrixXd ret(NNN, NNN);

		i1 = 0;
		// Charge derivatives
		if (ThM->ConstrainMuQ()) {
			int i2 = 0;

			double d1 = 0., d2 = 0.;

			if (ThM->ConstrainMuQ()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).ElectricCharge() * ThM->TPS()->Particle(i).ElectricCharge() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).BaryonCharge() * ThM->TPS()->Particle(i).ElectricCharge() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2) / ThM->QoverB();

				i2++;
			}


			if (ThM->ConstrainMuS()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).ElectricCharge() * ThM->TPS()->Particle(i).Strangeness() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).BaryonCharge() * ThM->TPS()->Particle(i).Strangeness() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2) / ThM->QoverB();

				i2++;
			}


			if (ThM->ConstrainMuC()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).ElectricCharge() * ThM->TPS()->Particle(i).Charm() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).BaryonCharge() * ThM->TPS()->Particle(i).Charm() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = (d1 / fBd - fQd / fBd / fBd * d2) / ThM->QoverB();

				i2++;
			}

			i1++;
		}


		// Strangeness derivatives
		if (ThM->ConstrainMuS()) {
			int i2 = 0;

			double d1 = 0., d2 = 0.;

			if (ThM->ConstrainMuQ()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).Strangeness()    * ThM->TPS()->Particle(i).ElectricCharge() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).AbsoluteStrangeness() * ThM->TPS()->Particle(i).ElectricCharge() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = d1 / fASd - fSd / fASd / fASd * d2;

				i2++;
			}


			if (ThM->ConstrainMuS()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).Strangeness()    * ThM->TPS()->Particle(i).Strangeness() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).AbsoluteStrangeness() * ThM->TPS()->Particle(i).Strangeness() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = d1 / fASd - fSd / fASd / fASd * d2;

				i2++;
			}


			if (ThM->ConstrainMuC()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).Strangeness()    * ThM->TPS()->Particle(i).Charm() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).AbsoluteStrangeness() * ThM->TPS()->Particle(i).Charm() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = d1 / fASd - fSd / fASd / fASd * d2;

				i2++;
			}

			i1++;
		}


		// Charm derivatives
		if (ThM->ConstrainMuC()) {
			int i2 = 0;

			double d1 = 0., d2 = 0.;

			if (ThM->ConstrainMuQ()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).Charm() * ThM->TPS()->Particle(i).ElectricCharge() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).AbsoluteCharm()  * ThM->TPS()->Particle(i).ElectricCharge() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = d1 / fACd - fCd / fACd / fACd * d2;

				i2++;
			}


			if (ThM->ConstrainMuS()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).Charm() * ThM->TPS()->Particle(i).Strangeness() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).AbsoluteCharm()  * ThM->TPS()->Particle(i).Strangeness() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = d1 / fACd - fCd / fACd / fACd * d2;

				i2++;
			}


			if (ThM->ConstrainMuC()) {
				d1 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d1 += ThM->TPS()->Particle(i).Charm() * ThM->TPS()->Particle(i).Charm() * ThM->Densities()[i] * m_wprim[i];
				d1 /= ThM->Parameters().T;

				d2 = 0.;
				for (int i = 0; i < m_wprim.size(); ++i)
					d2 += ThM->TPS()->Particle(i).AbsoluteCharm()  * ThM->TPS()->Particle(i).Charm() * ThM->Densities()[i] * m_wprim[i];
				d2 /= ThM->Parameters().T;

				ret(i1, i2) = d1 / fACd - fCd / fACd / fACd * d2;

				i2++;
			}

			i1++;
		}

		return ret;
	}

	void broydenEigen(vector<double> &xin, ThermalModelBase *ThM) {
		const double TOLF = 1.0e-8, EPS = 1.0e-8;
		const int MAXITS = 200;
		bool fl = 0;

		int NNN = 0;
		if (ThM->ConstrainMuQ()) NNN++;
		if (ThM->ConstrainMuS()) NNN++;
		if (ThM->ConstrainMuC()) NNN++;
		if (NNN == 0) {
			ThM->FillChemicalPotentials();
			return;
		}
		VectorXd xold(NNN), xnew(NNN), xdelta(NNN);
		VectorXd fold(NNN), fnew(NNN), fdelta(NNN);
		MatrixXd Jac(NNN, NNN), Jinv(NNN, NNN);

		NNN = 0;
		if (ThM->ConstrainMuQ()) { xold[NNN] = xin[0]; NNN++; }
		if (ThM->ConstrainMuS()) { xold[NNN] = xin[1]; NNN++; }
		if (ThM->ConstrainMuC()) { xold[NNN] = xin[2]; NNN++; }



		Jac = broydenEigenJacobian(xold, ThM);


		bool constrmuQ = ThM->ConstrainMuQ();
		bool constrmuS = ThM->ConstrainMuS();
		bool constrmuC = ThM->ConstrainMuC();
		bool repeat = false;
		NNN = 0;
		if (ThM->ConstrainMuQ()) {
			for (int j = 0; j < Jac.rows(); ++j)
				if (Jac(NNN, j) > 1.e8) { repeat = true; ThM->ConstrainMuQ(false); }
			double nQ = ThM->CalculateChargeDensity();
			double nB = ThM->CalculateBaryonDensity();
			if (abs(nQ) < 1.e-25 || abs(nB) < 1.e-25) { repeat = true; ThM->ConstrainMuQ(false); }
			NNN++;
		}
		if (ThM->ConstrainMuS()) {
			for (int j = 0; j < Jac.rows(); ++j)
				if (Jac(NNN, j) > 1.e8) { repeat = true; ThM->ConstrainMuS(false); }

			double nS = ThM->CalculateAbsoluteStrangenessDensity();
			if (abs(nS) < 1.e-25) { repeat = true; ThM->ConstrainMuS(false); }
			NNN++;
		}
		if (ThM->ConstrainMuC()) {
			for (int j = 0; j < Jac.rows(); ++j)
				if (Jac(NNN, j) > 1.e8) { repeat = true; ThM->ConstrainMuC(false); }
			double nC = ThM->CalculateAbsoluteCharmDensity();
			if (abs(nC) < 1.e-25) { repeat = true; ThM->ConstrainMuC(false); }
			NNN++;
		}
		if (repeat) {
			broydenEigen(xin, ThM);
			ThM->ConstrainMuQ(constrmuQ);
			ThM->ConstrainMuS(constrmuS);
			ThM->ConstrainMuC(constrmuC);
			return;
		}

		fold = broydenEigenFunc(xold, ThM);

		xnew = xold;
		fnew = fold;

		if (Jac.determinant() == 0.0)
		{
			printf("**WARNING** Singular Jacobian in BroydenEigen\n");
			return;
			//throw("singular Jacobian in broydenEigen");
		}
		Jinv = Jac.inverse();

		for (int iter = 1; iter < MAXITS; ++iter) {
			xnew = xold - Jinv * fold;
			fnew = broydenEigenFunc(xnew, ThM);
			xdelta = xnew - xold;
			fdelta = fnew - fold;
			double norm = 0.;
			for (int i = 0; i < NNN; ++i)
				for (int j = 0; j < NNN; ++j)
					norm += xdelta[i] * Jinv(i, j) * fdelta[j];
			//TVectorD p1(NNN);
			VectorXd p1(NNN);
			p1 = (xdelta - Jinv * fdelta);
			for (int i = 0; i < NNN; ++i) p1[i] *= 1. / norm;
			Jinv = Jinv + (p1 * xdelta.transpose()) * Jinv;
			/*Jac = broydenEigenJacobian(xnew, ThM);
			if (Jac.determinant()==0.0) throw("singular Jacobian in broydenEigen");
			Jinv = Jac.inverse();*/
			double maxdiff = 0.;
			for (int i = 0; i < fnew.size(); ++i) {
				maxdiff = max(maxdiff, fabs(fnew[i]));
			}
			if (maxdiff < TOLF) {
				return;
			}
			xold = xnew;
			fold = fnew;
		}
		printf("**WARNING** Reached maximum number of iterations in BroydenEigen\n");
		//throw("exceed MAXITS in broyden");
	}


	VectorXd broydenEigenFunc2(const VectorXd &xin, const vector<int> &vConstr, const vector<int> &vType, const vector<double> &vTotals, ThermalModelBase *ThM) {
		int NNN = 0;
		for (int i = 0; i < 4; ++i) NNN += vConstr[i];

		int i1 = 0;

		for (int i = 0; i < 4; ++i) {
			if (vConstr[i]) {
				if (i == 0) ThM->SetBaryonChemicalPotential(xin[i1]);
				if (i == 1) ThM->SetElectricChemicalPotential(xin[i1]);
				if (i == 2) ThM->SetStrangenessChemicalPotential(xin[i1]);
				if (i == 3) ThM->SetCharmChemicalPotential(xin[i1]);
				i1++;
			}
		}

		ThM->FillChemicalPotentials();
		ThM->CalculateDensities();

		vector<double> dens(4, 0.), absdens(4, 0.);
		if (vConstr[0]) {
			dens[0] = ThM->CalculateBaryonDensity();
			absdens[0] = ThM->CalculateAbsoluteBaryonDensity();
		}
		if (vConstr[1]) {
			dens[1] = ThM->CalculateChargeDensity();
			absdens[1] = ThM->CalculateAbsoluteChargeDensity();
		}
		if (vConstr[2]) {
			dens[2] = ThM->CalculateStrangenessDensity();
			absdens[2] = ThM->CalculateAbsoluteStrangenessDensity();
		}
		if (vConstr[3]) {
			dens[3] = ThM->CalculateCharmDensity();
			absdens[3] = ThM->CalculateAbsoluteCharmDensity();
		}

		VectorXd ret(NNN);

		i1 = 0;

		for (int i = 0; i < 4; ++i) {
			if (vConstr[i]) {
				if (vType[i] == 0)
					//ret[i1] = (dens[i1] * ThM->Parameters().V - vTotals[i1]) / vTotals[i1];
					ret[i1] = (dens[i] * ThM->Parameters().V - vTotals[i]) / vTotals[i];
				else
					//ret[i1] = dens[i1] / absdens[i1];
					ret[i1] = dens[i] / absdens[i];
				i1++;
			}
		}

		return ret;
	}

	MatrixXd broydenEigenJacobian2(const VectorXd &xin, const vector<int> &vConstr, const vector<int> &vType, const vector<double> &vTotals, ThermalModelBase *ThM) {
		int NNN = 0;
		for (int i = 0; i < 4; ++i) NNN += vConstr[i];

		int i1 = 0;

		for (int i = 0; i < 4; ++i) {
			if (vConstr[i]) {
				if (i == 0) ThM->SetBaryonChemicalPotential(xin[i1]);
				if (i == 1) ThM->SetElectricChemicalPotential(xin[i1]);
				if (i == 2) ThM->SetStrangenessChemicalPotential(xin[i1]);
				if (i == 3) ThM->SetCharmChemicalPotential(xin[i1]);
				i1++;
			}
		}

		vector<double> tfug(4, 0.);
		tfug[0] = exp(ThM->Parameters().muB / ThM->Parameters().T);
		tfug[1] = exp(ThM->Parameters().muQ / ThM->Parameters().T);
		tfug[2] = exp(ThM->Parameters().muS / ThM->Parameters().T);
		tfug[3] = exp(ThM->Parameters().muC / ThM->Parameters().T);

		ThM->FillChemicalPotentials();
		ThM->CalculateDensities();

		vector<double> dens(4, 0.), absdens(4, 0.);
		if (vConstr[0]) {
			dens[0] = ThM->CalculateBaryonDensity();
			absdens[0] = ThM->CalculateAbsoluteBaryonDensity();
		}
		if (vConstr[1]) {
			dens[1] = ThM->CalculateChargeDensity();
			absdens[1] = ThM->CalculateAbsoluteChargeDensity();
		}
		if (vConstr[2]) {
			dens[2] = ThM->CalculateStrangenessDensity();
			absdens[2] = ThM->CalculateAbsoluteStrangenessDensity();
		}
		if (vConstr[3]) {
			dens[3] = ThM->CalculateCharmDensity();
			absdens[3] = ThM->CalculateAbsoluteCharmDensity();
		}

		vector<double> m_wprim;
		m_wprim.resize(ThM->Densities().size());
		for (int i = 0; i < m_wprim.size(); ++i) m_wprim[i] = ThM->CalculateParticleScaledVariance(i);

		vector< vector<double> > deriv(4, vector<double>(4)), derivabs(4, vector<double>(4));
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j) {
				deriv[i][j] = 0.;
				for (int part = 0; part < m_wprim.size(); ++part)
					deriv[i][j] += ThM->TPS()->Particles()[part].GetCharge(i) * ThM->TPS()->Particles()[part].GetCharge(j) * ThM->Densities()[part] * m_wprim[part];
				deriv[i][j] /= ThM->Parameters().T;

				derivabs[i][j] = 0.;
				for (int part = 0; part < m_wprim.size(); ++part)
					derivabs[i][j] += ThM->TPS()->Particles()[part].GetAbsCharge(i) * ThM->TPS()->Particles()[part].GetCharge(j) * ThM->Densities()[part] * m_wprim[part];
				derivabs[i][j] /= ThM->Parameters().T;
			}


		MatrixXd ret(NNN, NNN);

		i1 = 0;

		for (int i = 0; i < 4; ++i) {
			if (vConstr[i]) {
				int i2 = 0;
				for (int j = 0; j < 4; ++j)
					if (vConstr[j]) {
						ret(i1, i2) = 0.;
						if (vType[i] == 0)
							//ret(i1, i2) = deriv[i1][i2] * ThM->Parameters().V / vTotals[i1];
							ret(i1, i2) = deriv[i][j] * ThM->Parameters().V / vTotals[i];
						else
							//ret(i1, i2) = deriv[i1][i2] / absdens[i1] - dens[i1] / absdens[i1] / absdens[i1] * derivabs[i1][i2];
							ret(i1, i2) = deriv[i][j] / absdens[i] - dens[i] / absdens[i] / absdens[i] * derivabs[i][j];
						i2++;
					}
				i1++;
			}
		}

		return ret;
	}

	void broydenEigen2(vector<double> &xin, const vector<int> &vConstr, const vector<int> &vType, const vector<double> &vTotals, ThermalModelBase *ThM) {
		const double TOLF = 1.0e-8, EPS = 1.0e-8;
		const int MAXITS = 400;
		bool fl = 0;


		int NNN = 0;
		for (int i = 0; i < 4; ++i) NNN += vConstr[i];
		if (NNN == 0) return;
		VectorXd xold(NNN), xnew(NNN), xdelta(NNN);
		VectorXd fold(NNN), fnew(NNN), fdelta(NNN);
		MatrixXd Jac(NNN, NNN), Jinv(NNN, NNN);

		int i1 = 0;

		for (int i = 0; i < 4; ++i) {
			if (vConstr[i]) {
				xold[i1] = xin[i];
				i1++;
			}
		}

		Jac = broydenEigenJacobian2(xold, vConstr, vType, vTotals, ThM);

		fold = broydenEigenFunc2(xold, vConstr, vType, vTotals, ThM);

		xnew = xold;
		fnew = fold;

		if (Jac.determinant() == 0.0) throw("singular Jacobian in broydenEigen");
		Jinv = Jac.inverse();

		for (int iter = 1; iter < MAXITS; ++iter) {
			xnew = xold - Jinv * fold;
			fnew = broydenEigenFunc2(xnew, vConstr, vType, vTotals, ThM);
			xdelta = xnew - xold;
			fdelta = fnew - fold;
			double norm = 0.;
			for (int i = 0; i < NNN; ++i)
				for (int j = 0; j < NNN; ++j)
					norm += xdelta[i] * Jinv(i, j) * fdelta[j];
			VectorXd p1(NNN);
			p1 = (xdelta - Jinv * fdelta);
			for (int i = 0; i < NNN; ++i) p1[i] *= 1. / norm;
			Jinv = Jinv + (p1 * xdelta.transpose()) * Jinv;
			//Jac  = broydenEigenJacobian2(xnew, vConstr, vType, vTotals, ThM);
			//if (Jac.determinant()==0.0) throw("singular Jacobian in broydenEigen");
			//Jinv = Jac.inverse();
			double maxdiff = 0.;
			for (int i = 0; i < fnew.size(); ++i) {
				maxdiff = max(maxdiff, fabs(fnew[i]));
			}
			if (maxdiff < TOLF) {
				return;
			}
			xold = xnew;
			fold = fnew;
		}
		throw("exceed MAXITS in broyden");
	}

}

using namespace ThermalModelBaseNamespace;


ThermalModelBase::ThermalModelBase(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
	m_TPS(TPS_), 
	m_Parameters(params),
	m_UseWidth(false),
	m_Calculated(false),
	m_FluctuationsCalculated(false),
	m_GCECalculated(false),
	m_NormBratio(false),
	m_QuantumStats(true),
	m_MaxDiff(0.),
	m_useOpenMP(0)
{
	m_QBgoal = 0.4;
	m_Chem.resize(m_TPS->Particles().size());
	m_Volume = params.V;
	m_densities.resize(m_TPS->Particles().size());
	m_densitiestotal.resize(m_TPS->Particles().size());
	m_densitiestotalweak.resize(m_TPS->Particles().size());
	m_wprim.resize(m_TPS->Particles().size());
	m_wtot.resize(m_TPS->Particles().size());
	m_skewprim.resize(m_TPS->Particles().size());
	m_skewtot.resize(m_TPS->Particles().size());
	m_kurtprim.resize(m_TPS->Particles().size());
	m_kurttot.resize(m_TPS->Particles().size());

	m_ConstrainMuB = m_ConstrainMuC = m_ConstrainMuQ = m_ConstrainMuS = true;

	m_Susc.resize(4);
	for (int i = 0; i < 4; ++i) m_Susc[i].resize(4);

	m_NormBratio = false;
	
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		const ThermalParticle &tpart = m_TPS->Particles()[i];
		for (int j = 0; j < tpart.Decays().size(); ++j) {
			if (tpart.Decays()[j].mBratio != tpart.DecaysOriginal()[j].mBratio)
				m_NormBratio = true;
		}
	}

	m_Ensemble = GCE;
	m_InteractionModel = Ideal;

	m_ValidityLog = "";
}


void ThermalModelBase::SetUseWidth(bool useWidth)
{
	if (!useWidth && m_TPS->ResonanceWidthIntegrationType() == ThermalParticle::eBW) {
		m_TPS->SetResonanceWidthIntegrationType(ThermalParticle::TwoGamma);
		m_TPS->ProcessDecays();
	}
	m_UseWidth = useWidth;
}

void ThermalModelBase::SetUseWidth(ThermalParticle::ResonanceWidthIntegration type)
{
	m_UseWidth = (type != ThermalParticle::ZeroWidth);
	m_TPS->SetResonanceWidthIntegrationType(type);
}

void ThermalModelBase::SetNormBratio(bool normBratio) {
	if (normBratio != m_NormBratio) {
		m_NormBratio = normBratio;
		if (m_NormBratio) {
			m_TPS->NormalizeBranchingRatios();
		}
		else {
			m_TPS->RestoreBranchingRatios();
		}
	}
}


void ThermalModelBase::ResetChemicalPotentials() {
	m_Parameters.muS = m_Parameters.muB / 5.;
	m_Parameters.muQ = -m_Parameters.muB / 50.;
	m_Parameters.muC = 0.;
}


void ThermalModelBase::SetParameters(const ThermalModelParameters& params) {
	m_Parameters = params;
	m_Volume = m_Parameters.V;
	m_Calculated = false;
}

void ThermalModelBase::SetTemperature(double T)
{
	m_Parameters.T = T;
	m_Calculated = false;
}

void ThermalModelBase::SetBaryonChemicalPotential(double muB)
{
	m_Parameters.muB = muB;
	FillChemicalPotentials();
	m_Calculated = false;
}

void ThermalModelBase::SetElectricChemicalPotential(double muQ)
{
	m_Parameters.muQ = muQ;
	FillChemicalPotentials();
	m_Calculated = false;
}

void ThermalModelBase::SetStrangenessChemicalPotential(double muS)
{
	m_Parameters.muS = muS;
	FillChemicalPotentials();
	m_Calculated = false;
}

void ThermalModelBase::SetCharmChemicalPotential(double muC)
{
	m_Parameters.muC = muC;
	FillChemicalPotentials();
	m_Calculated = false;
}

void ThermalModelBase::SetGammaS(double gammaS)
{
	m_Parameters.gammaS = gammaS;
	m_Calculated = false;
}

void ThermalModelBase::SetGammaC(double gammaC)
{
	m_Parameters.gammaC = gammaC;
	m_Calculated = false;
}

void ThermalModelBase::SetBaryonCharge(int B)
{
	m_Parameters.B = B;
	m_Calculated = false;
}

void ThermalModelBase::SetElectricCharge(int Q)
{
	m_Parameters.Q = Q;
	m_Calculated = false;
}

void ThermalModelBase::SetStrangeness(int S)
{
	m_Parameters.S = S;
	m_Calculated = false;
}

void ThermalModelBase::SetCharm(int C)
{
	m_Parameters.C = C;
	m_Calculated = false;
}

void ThermalModelBase::SetGammaq(double gammaq)
{
	m_Parameters.gammaq = gammaq;
	m_Calculated = false;
}


void ThermalModelBase::ChangeTPS(ThermalParticleSystem *TPS_) {
	m_TPS = TPS_;
	m_Chem.resize(m_TPS->Particles().size());
	m_densities.resize(m_TPS->Particles().size());
	m_densitiestotal.resize(m_TPS->Particles().size());
	m_densitiestotalweak.resize(m_TPS->Particles().size());
	m_wprim.resize(m_TPS->Particles().size());
	m_wtot.resize(m_TPS->Particles().size());
	m_skewprim.resize(m_TPS->Particles().size());
	m_skewtot.resize(m_TPS->Particles().size());
	m_kurtprim.resize(m_TPS->Particles().size());
	m_kurttot.resize(m_TPS->Particles().size());
	m_Calculated = false;
}

void ThermalModelBase::SetStatistics(bool stats) {
	m_QuantumStats = stats;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		m_TPS->Particle(i).UseStatistics(stats);
}

void ThermalModelBase::SetResonanceWidthIntegrationType(ThermalParticle::ResonanceWidthIntegration type)
{
	if (!m_UseWidth) {
		printf("**WARNING** ThermalModelBase::SetResonanceWidthIntegrationType: Using resonance widths is switched off!\n");
		m_TPS->SetResonanceWidthIntegrationType(ThermalParticle::TwoGamma);
	}
	else
		m_TPS->SetResonanceWidthIntegrationType(type);
}

void ThermalModelBase::FillChemicalPotentials() {
	m_Chem.resize(m_TPS->Particles().size());
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		m_Chem[i] = m_TPS->Particles()[i].BaryonCharge() * m_Parameters.muB + m_TPS->Particles()[i].Strangeness() * m_Parameters.muS + m_TPS->Particles()[i].ElectricCharge() * m_Parameters.muQ + m_TPS->Particles()[i].Charm() * m_Parameters.muC;
}

void ThermalModelBase::SetChemicalPotentials(const std::vector<double>& chem)
{
	if (chem.size() != m_TPS->Particles().size()) {
		printf("**WARNING** %s::SetChemicalPotentials(const std::vector<double> & chem): size of chem does not match number of hadrons in the list", m_TAG.c_str());
		return;
	}
	m_Chem = chem;
}


void ThermalModelBase::CalculateFeeddown() {
	if (m_UseWidth && m_TPS->ResonanceWidthIntegrationType() == ThermalParticle::eBW) {
		for (int i = 0; i < m_TPS->Particles().size(); ++i) {
			m_TPS->Particle(i).CalculateThermalBranchingRatios(m_Parameters, m_UseWidth, m_Chem[i], MuShift(i));
		}
		m_TPS->ProcessDecays();
	}

	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		m_densitiestotal[i] = m_densities[i];
		for (int j = 0; j < m_TPS->Particles()[i].DecayContributions().size(); ++j)
			if (i != m_TPS->Particles()[i].DecayContributions()[j].second) m_densitiestotal[i] += m_TPS->Particles()[i].DecayContributions()[j].first * m_densities[m_TPS->Particles()[i].DecayContributions()[j].second];
	}

	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		m_densitiestotalweak[i] = m_densities[i];
		for (int j = 0; j < m_TPS->Particles()[i].WeakDecayContributions().size(); ++j)
			if (i != m_TPS->Particles()[i].WeakDecayContributions()[j].second) m_densitiestotalweak[i] += m_TPS->Particles()[i].WeakDecayContributions()[j].first * m_densities[m_TPS->Particles()[i].WeakDecayContributions()[j].second];
	}

}


void ThermalModelBase::FixParameters() {
	if (fabs(m_Parameters.muB) < 1e-6) {
		if (m_ConstrainMuS)
			m_Parameters.muS = 0.;
		if (m_ConstrainMuQ)
			m_Parameters.muQ = 0.;
		if (m_ConstrainMuC)
			m_Parameters.muC = 0.;
		//m_Parameters.muS = m_Parameters.muQ = m_Parameters.muC = 0.;
		FillChemicalPotentials();
		CalculateDensities();
		return;
	}
	double suppr = 10;
	if (m_Parameters.muB > 0.150) suppr = 8.;
	if (m_Parameters.muB > 0.300) suppr = 7.;
	if (m_Parameters.muB > 0.450) suppr = 6.;
	if (m_Parameters.muB > 0.600) suppr = 6.;
	if (m_Parameters.muB > 0.750) suppr = 5.;
	if (m_Parameters.muB > 0.900) suppr = 4.;
	if (m_Parameters.muB > 1.000) suppr = 3.;
	if (m_ConstrainMuS)
		m_Parameters.muS = m_Parameters.muB / suppr;
	if (m_ConstrainMuQ)
		m_Parameters.muQ = -m_Parameters.muB / suppr / 10.;
	//m_Parameters.muC = 0.;
	if (m_ConstrainMuC)
		m_Parameters.muC = -m_Parameters.muS;

	FixParametersNoReset();
}

void ThermalModelBase::FixParametersNoReset() {
	//printf("FixQB: %lf\n", m_Parameters.muB);
	if (fabs(m_Parameters.muB) < 1e-6) {
		m_Parameters.muS = m_Parameters.muQ = m_Parameters.muC = 0.;
		FillChemicalPotentials();
		CalculateDensities();
		return;
	}

	m_ConstrainMuQ &= (m_TPS->hasCharged() && m_TPS->hasBaryons());
	m_ConstrainMuS &= m_TPS->hasStrange();
	m_ConstrainMuC &= m_TPS->hasCharmed();

	//printf("%d %d %d\n", (int)m_ConstrainMuQ, (int)m_ConstrainMuS, (int)m_ConstrainMuC);

	vector<double> x22(3);
	x22[0] = m_Parameters.muQ;
	x22[1] = m_Parameters.muS;
	x22[2] = m_Parameters.muC;
	vector<double> x2(3), xinit(3);
	xinit[0] = x2[0] = m_Parameters.muQ;
	xinit[1] = x2[1] = m_Parameters.muS;
	xinit[2] = x2[2] = m_Parameters.muC;
	int iter = 0, iterMAX = 2;
	//m_ConstrainMuS = 0;
	while (iter < iterMAX) {
		try {
			if (0 && m_ConstrainMuQ && m_ConstrainMuS && !m_ConstrainMuC)
				broyden2(x22, ThermalModelBaseNamespace::function2, this);
			else
				broydenEigen(x22, this);
		}
		catch (...) {
		}
		break;
		iter++;
	}
}

void ThermalModelBase::SolveChemicalPotentials(double totB, double totQ, double totS, double totC,
	double muBinit, double muQinit, double muSinit, double muCinit,
	bool ConstrMuB, bool ConstrMuQ, bool ConstrMuS, bool ConstrMuC) {
	m_Parameters.muB = muBinit;
	m_Parameters.muS = muSinit;
	m_Parameters.muQ = muQinit;
	m_Parameters.muC = muCinit;
	if (totB == 0.0 && totQ == 0.0 && totS == 0.0 && totC == 0.0) {
		m_Parameters.muB = 0.;
		m_Parameters.muS = 0.;
		m_Parameters.muQ = 0.;
		m_Parameters.muC = 0.;
		FillChemicalPotentials();
		CalculateDensities();
		return;
	}
	vector<int> vConstr(4, 1);
	vector<int> vType(4, 0);

	vConstr[0] = m_TPS->hasBaryons() && ConstrMuB;
	vConstr[1] = m_TPS->hasCharged() && ConstrMuQ;
	vConstr[2] = m_TPS->hasStrange() && ConstrMuS;
	vConstr[3] = m_TPS->hasCharmed() && ConstrMuC;

	vType[0] = (int)(totB == 0.0);
	vType[1] = (int)(totQ == 0.0);
	vType[2] = (int)(totS == 0.0);
	vType[3] = (int)(totC == 0.0);

	vector<double> vTotals(4);
	vTotals[0] = totB;
	vTotals[1] = totQ;
	vTotals[2] = totS;
	vTotals[3] = totC;

	vector<double> xin(4, 0.);
	xin[0] = muBinit;
	xin[1] = muQinit;
	xin[2] = muSinit;
	xin[3] = muCinit;

	try {
		broydenEigen2(xin, vConstr, vType, vTotals, this);
	}
	catch (...) {
	}
}

void ThermalModelBase::ValidateCalculation()
{
	m_ValidityLog = "";

	char cc[1000];

	m_LastCalculationSuccessFlag = true;
	for (int i = 0; i < m_densities.size(); ++i) {
		if (m_densities[i] != m_densities[i]) {
			m_LastCalculationSuccessFlag = false;
			
			sprintf(cc, "**WARNING** Density for particle %d (%s) is NaN!\n\n", m_TPS->Particle(i).PdgId(), m_TPS->Particle(i).Name().c_str());
			printf("%s", cc);

			m_ValidityLog.append(cc);
		}
		//m_LastCalculationSuccessFlag &= (m_densities[i] == m_densities[i]);
	}
}

void ThermalModelBase::FixParameters(double QB) {
	m_QBgoal = QB;
	FixParameters();
}

double ThermalModelBase::GetParticlePrimordialDensity(unsigned int part) {
	if (!m_Calculated) CalculateDensities();
	if (part >= m_densities.size()) return 0.;
	return m_densities[part];
}

double ThermalModelBase::GetParticleTotalDensity(unsigned int part) {
	if (!m_Calculated) CalculateDensities();
	if (part >= m_densitiestotal.size()) return 0.;
	return m_densitiestotal[part];
}

double ThermalModelBase::CalculateHadronDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_densities[i];

	return ret;
}

double ThermalModelBase::CalculateBaryonDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_TPS->Particles()[i].BaryonCharge() * m_densities[i];

	return ret;
}

double ThermalModelBase::CalculateChargeDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_TPS->Particles()[i].ElectricCharge() * m_densities[i];

	return ret;
}

double ThermalModelBase::CalculateStrangenessDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_TPS->Particles()[i].Strangeness() * m_densities[i];

	return ret;
}

double ThermalModelBase::CalculateCharmDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_TPS->Particles()[i].Charm() * m_densities[i];
	return ret;
}

double ThermalModelBase::CalculateAbsoluteBaryonDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += fabs((double)m_TPS->Particles()[i].BaryonCharge()) * m_densities[i];
	return ret;
}

double ThermalModelBase::CalculateAbsoluteChargeDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += fabs((double)m_TPS->Particles()[i].ElectricCharge()) * m_densities[i];
	return ret;
}

double ThermalModelBase::CalculateAbsoluteStrangenessDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_TPS->Particles()[i].AbsoluteStrangeness() * m_densities[i];
	return ret;
}

double ThermalModelBase::CalculateAbsoluteCharmDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_TPS->Particles()[i].AbsoluteCharm() * m_densities[i];
	return ret;
}

double ThermalModelBase::CalculateArbitraryChargeDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;
	for (int i = 0; i < m_TPS->Particles().size(); ++i)
		ret += m_TPS->Particles()[i].ArbitraryCharge() * m_densities[i];
	return ret;
}


double ThermalModelBase::GetDensity(int PDGID, int feeddown) {
	std::vector<double> *dens;
	if (feeddown == 0) dens = &m_densities;
	else if (feeddown == 1) dens = &m_densitiestotal;
	else dens = &m_densitiestotalweak;

	if (m_TPS->PdgToId(PDGID) != -1)
		return dens->operator[](m_TPS->PdgToId(PDGID));

	// 1 - Npart
	if (PDGID == 1) return CalculateBaryonDensity();

	// 33340 - \Omega + \Omegabar
	if (PDGID == 33340 && m_TPS->PdgToId(3334) != -1 && m_TPS->PdgToId(-3334) != -1)
		return dens->operator[](m_TPS->PdgToId(3334)) + dens->operator[](m_TPS->PdgToId(-3334));

	// 22120 - nucleons
	if (PDGID == 22120 && m_TPS->PdgToId(2212) != -1 && m_TPS->PdgToId(2112)  != -1)
		return  dens->operator[](m_TPS->PdgToId(2212)) + dens->operator[](m_TPS->PdgToId(2112));

	printf("**WARNING** %s: Density with PDG ID %d not found!\n", m_TAG.c_str(), PDGID);

	return 0.;
}


std::vector<double> ThermalModelBase::GetIdealGasDensities() const {
	std::vector<double> ret = m_densities;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		ret[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i], 0.);
	}
	return ret;
}

double ThermalModelBase::ChargedMultiplicity(int type)
{
	if (!m_Calculated) CalculateDensities();
	double ret = 0.0;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		int tQ = m_TPS->Particles()[i].ElectricCharge();
		bool fl = false;
		if (type == 0  && tQ != 0)
			fl = true;
		if (type == 1  && tQ > 0)
			fl = true;
		if (type == -1 && tQ < 0)
			fl = true;
		if (fl)
			ret += m_densities[i];
	}
	return ret * Volume();
}

double ThermalModelBase::ChargedScaledVariance(int type)
{
	if (!m_FluctuationsCalculated) {
		printf("**WARNING** %s: ChargedScaledVariance(int): Fluctuations were not calculated", m_TAG.c_str());
		return 1.;
	}
	double ret = 0.0;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		int tQ = m_TPS->Particles()[i].ElectricCharge();
		bool fl = false;
		if (type == 0 && tQ != 0)
			fl = true;
		if (type == 1 && tQ > 0)
			fl = true;
		if (type == -1 && tQ < 0)
			fl = true;
		if (fl) {
			for (int j = 0; j < m_TPS->Particles().size(); ++j) {
				int tQ2 = m_TPS->Particles()[j].ElectricCharge();
				bool fl2 = false;
				if (type == 0 && tQ2 != 0)
					fl2 = true;
				if (type == 1 && tQ2 > 0)
					fl2 = true;
				if (type == -1 && tQ2 < 0)
					fl2 = true;

				if (fl2) {
					ret += m_PrimCorrel[i][j];
				}
			}
		}
	}
	return ret * m_Parameters.T * Volume() / ChargedMultiplicity(type);
}

double ThermalModelBase::ChargedMultiplicityFinal(int type)
{
	if (!m_Calculated) CalculateDensities();

	int op = type;
	if (type == -1)
		op = 2;

	double ret = 0.0;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		ret += m_densities[i] * m_TPS->Particles()[i].Nch()[op];
	}
	return ret * Volume();
}

double ThermalModelBase::ChargedScaledVarianceFinal(int type)
{
	if (!m_FluctuationsCalculated) {
		printf("**WARNING** %s: ChargedScaledVarianceFinal(int): Fluctuations were not calculated", m_TAG.c_str());
		return 1.;
	}
	int op = type;
	if (type == -1)
		op = 2;
	double ret = 0.0;
	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		ret += m_densities[i] * Volume() * m_TPS->Particles()[i].DeltaNch()[op];
		for (int j = 0; j < m_TPS->Particles().size(); ++j) {
			ret += m_PrimCorrel[i][j] * m_Parameters.T * Volume() * m_TPS->Particles()[i].Nch()[op] * m_TPS->Particles()[j].Nch()[op];
		}
	}
	return ret / ChargedMultiplicityFinal(type);
}

void ThermalModelBase::CalculateTwoParticleCorrelations() {
	printf("**WARNING** %s: Calculation of two-particle correlations and fluctuations is not implemented", m_TAG.c_str());
}

void ThermalModelBase::CalculateTwoParticleFluctuationsDecays()
{
	int NN = m_densities.size();

	// Only diagonal for now
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
}

void ThermalModelBase::CalculateSusceptibilityMatrix()
{
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

void ThermalModelBase::CalculateFluctuations() {
	printf("**WARNING** %s: Calculation of fluctuations is not implemented", m_TAG.c_str());
}
