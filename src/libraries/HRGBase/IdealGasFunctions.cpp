#include "HRGBase/IdealGasFunctions.h"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cfloat>
#include <vector>

#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"

using namespace std;

namespace IdealGasFunctions {

	double BoltzmannDensity(double T, double mu, double m, double deg) {
		return deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::BesselKexp(2, m / T) * exp((mu - m) / T) * xMath::GeVtoifm3();
	}

	double BoltzmannPressure(double T, double mu, double m, double deg) {
		return T * BoltzmannDensity(T, mu, m, deg);
	}

	double BoltzmannEnergyDensity(double T, double mu, double m, double deg) {
		return (3 * T + m * xMath::BesselK1exp(m / T) / xMath::BesselKexp(2, m / T)) * BoltzmannDensity(T, mu, m, deg);
	}

	double BoltzmannEntropyDensity(double T, double mu, double m, double deg) {
		return (BoltzmannPressure(T, mu, m, deg) + BoltzmannEnergyDensity(T, mu, m, deg) - mu * BoltzmannDensity(T, mu, m, deg)) / T;
	}

	double BoltzmannScalarDensity(double T, double mu, double m, double deg) {
		return deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::BesselKexp(1, m / T) * exp((mu - m) / T) * xMath::GeVtoifm3();
	}

	double BoltzmannTdndmu(int N, double T, double mu, double m, double deg)
	{
		return BoltzmannDensity(T, mu, m, deg);
	}

	double BoltzmannChiN(int N, double T, double mu, double m, double deg)
	{
		return BoltzmannTdndmu(N - 1, T, mu, m, deg) / pow(T, 3) / xMath::GeVtoifm3();
	}

	double QuantumClusterExpansionDensity(int statistics, double T, double mu, double m, double deg, int order)
	{
		double sign = 1.;
		bool signchange = true;
		if (statistics == 1) //Fermi
			signchange = true;
		else if (statistics == -1) //Bose
			signchange = false;
		else
			return BoltzmannDensity(T, mu, m, deg);

		double tfug  = exp((mu - m) / T);
		double cfug  = tfug;
		double moverT = m / T;
		double ret = 0.;
		for (int i = 1; i <= order; ++i) {
			ret += sign * xMath::BesselKexp(2, i*moverT) * cfug / static_cast<double>(i);
			cfug *= tfug;
			if (signchange) sign = -sign;
		}
		ret *= deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
		return ret;
	}

	double QuantumClusterExpansionPressure(int statistics, double T, double mu, double m, double deg, int order)
	{
		double sign = 1.;
		bool signchange = true;
		if (statistics == 1) //Fermi
			signchange = true;
		else if (statistics == -1) //Bose
			signchange = false;
		else
			return BoltzmannPressure(T, mu, m, deg);

		double tfug = exp((mu - m) / T);
		double cfug = tfug;
		double moverT = m / T;
		double ret = 0.;
		for (int i = 1; i <= order; ++i) {
			ret  += sign * xMath::BesselKexp(2, i*moverT) * cfug / static_cast<double>(i) / static_cast<double>(i);
			cfug *= tfug;
			if (signchange) sign = -sign;
		}
		ret *= deg * m * m * T * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
		return ret;
	}

	double QuantumClusterExpansionEnergyDensity(int statistics, double T, double mu, double m, double deg, int order)
	{
		double sign = 1.;
		bool signchange = true;
		if (statistics == 1) //Fermi
			signchange = true;
		else if (statistics == -1) //Bose
			signchange = false;
		else
			return BoltzmannEnergyDensity(T, mu, m, deg);

		double tfug = exp((mu - m) / T);
		double cfug = tfug;
		double moverT = m / T;
		double ret = 0.;
		for (int i = 1; i <= order; ++i) {
			ret += sign * (xMath::BesselKexp(1, i*moverT) + 3. * xMath::BesselKexp(2, i*moverT) / moverT / static_cast<double>(i)) * cfug / static_cast<double>(i);
			cfug *= tfug;
			if (signchange) sign = -sign;
		}
		ret *= deg * m * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
		return ret;
	}

	double QuantumClusterExpansionEntropyDensity(int statistics, double T, double mu, double m, double deg, int order)
	{
		return (QuantumClusterExpansionPressure(statistics, T, mu, m, deg, order) + QuantumClusterExpansionEnergyDensity(statistics, T, mu, m, deg, order) - mu * QuantumClusterExpansionDensity(statistics, T, mu, m, deg, order)) / T;
	}

	double QuantumClusterExpansionScalarDensity(int statistics, double T, double mu, double m, double deg, int order)
	{
		double sign = 1.;
		bool signchange = true;
		if (statistics == 1) //Fermi
			signchange = true;
		else if (statistics == -1) //Bose
			signchange = false;
		else
			return BoltzmannScalarDensity(T, mu, m, deg);

		double tfug = exp((mu - m) / T);
		double cfug = tfug;
		double moverT = m / T;
		double ret = 0.;
		for (int i = 1; i <= order; ++i) {
			ret += sign * xMath::BesselKexp(1, i*moverT) * cfug / static_cast<double>(i);
			cfug *= tfug;
			if (signchange) sign = -sign;
		}
		ret *= deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
		return ret;
	}

	double QuantumClusterExpansionTdndmu(int N, int statistics, double T, double mu, double m, double deg, int order)
	{
		double sign = 1.;
		bool signchange = true;
		if (statistics == 1) //Fermi
			signchange = true;
		else if (statistics == -1) //Bose
			signchange = false;
		else
			return BoltzmannTdndmu(N, T, mu, m, deg);

		double tfug = exp((mu - m) / T);
		double cfug = tfug;
		double moverT = m / T;
		double ret = 0.;
		for (int i = 1; i <= order; ++i) {
			ret += sign * xMath::BesselKexp(2, i*moverT) * cfug * pow(static_cast<double>(i), N - 1);
			cfug *= tfug;
			if (signchange) sign = -sign;
		}
		ret *= deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
		return ret;
	}

	double QuantumClusterExpansionChiN(int N, int statistics, double T, double mu, double m, double deg, int order)
	{
		return QuantumClusterExpansionTdndmu(N - 1, statistics, T, mu, m, deg, order) / pow(T, 3) / xMath::GeVtoifm3();
	}


	// Gauss-Legendre 32-point quadrature for [0,1] interval
	const double *legx32 = NumericalIntegration::coefficients_xleg32_zeroone;
	const double *legw32 = NumericalIntegration::coefficients_wleg32_zeroone;
	// Gauss-Laguerre 32-point quadrature for [0,\infty] interval
	const double *lagx32 = NumericalIntegration::coefficients_xlag32;
	const double *lagw32 = NumericalIntegration::coefficients_wlag32;

	double QuantumNumericalIntegrationDensity(int statistics, double T, double mu, double m, double deg)
	{
		if (statistics == 0)           return BoltzmannDensity(T, mu, m, deg);
		if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuDensity(T, mu, m, deg);
		if (statistics == -1 && mu > m) {
			printf("**WARNING** QuantumNumericalIntegrationDensity: Bose-Einstein condensation\n");
			return 0.;
		}

		double ret = 0.;
		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = lagx32[i];
			ret += lagw32[i] * T * tx * T * tx * T / (exp(sqrt(tx*tx + moverT*moverT) - muoverT) + statistics);
		}

		ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret;
	}

	double QuantumNumericalIntegrationPressure(int statistics, double T, double mu, double m, double deg)
	{
		if (statistics == 0)           return BoltzmannPressure(T, mu, m, deg);
		if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuPressure(T, mu, m, deg);
		if (statistics == -1 && mu > m) {
			printf("**WARNING** QuantumNumericalIntegrationPressure: Bose-Einstein condensation\n");
			return 0.;
		}

		double ret = 0.;
		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = lagx32[i];
			double x2 = tx * T * tx * T;
			double E = sqrt(tx*tx + moverT*moverT);
			ret += lagw32[i] * T * x2 * x2 / E / T / (exp(E - muoverT) + statistics);
		}

		ret *= deg / 6. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret;
	}

	double QuantumNumericalIntegrationEnergyDensity(int statistics, double T, double mu, double m, double deg)
	{
		if (statistics == 0)           return BoltzmannEnergyDensity(T, mu, m, deg);
		if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuEnergyDensity(T, mu, m, deg);
		if (statistics == -1 && mu > m) {
			printf("**WARNING** QuantumNumericalIntegrationEnergyDensity: Bose-Einstein condensation\n");
			return 0.;
		}

		double ret = 0.;
		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = lagx32[i];
			ret += lagw32[i] * T * tx * T * tx * T * sqrt(tx*tx + moverT*moverT) * T / (exp(sqrt(tx*tx + moverT*moverT) - muoverT) + statistics);
		}

		ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret;
	}

	double QuantumNumericalIntegrationEntropyDensity(int statistics, double T, double mu, double m, double deg)
	{
		return (QuantumNumericalIntegrationPressure(statistics, T, mu, m, deg) + QuantumNumericalIntegrationEnergyDensity(statistics, T, mu, m, deg) - mu * QuantumNumericalIntegrationDensity(statistics, T, mu, m, deg)) / T;
	}

	double QuantumNumericalIntegrationScalarDensity(int statistics, double T, double mu, double m, double deg)
	{
		if (statistics == 0)           return BoltzmannScalarDensity(T, mu, m, deg);
		if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuScalarDensity(T, mu, m, deg);
		if (statistics == -1 && mu > m) {
			printf("**WARNING** QuantumNumericalIntegrationScalarDensity: Bose-Einstein condensation\n");
			return 0.;
		}

		double ret = 0.;
		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = lagx32[i];
			ret += lagw32[i] * T * tx * T * tx * T * moverT / sqrt(tx*tx + moverT*moverT) / (exp(sqrt(tx*tx + moverT*moverT) - muoverT) + statistics);
		}

		ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret;
	}

	double QuantumNumericalIntegrationT1dn1dmu1(int statistics, double T, double mu, double m, double deg)
	{
		if (statistics == 0)           return BoltzmannTdndmu(1, T, mu, m, deg);
		if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuT1dn1dmu1(T, mu, m, deg);
		if (statistics == -1 && mu > m) {
			printf("**WARNING** QuantumNumericalIntegrationT1dn1dmu1: Bose-Einstein condensation\n");
			return 0.;
		}

		double ret = 0.;
		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = lagx32[i];
			double Eexp = exp(sqrt(tx*tx + moverT*moverT) - muoverT);
			ret += lagw32[i] * T * tx * T * tx * T * Eexp / (Eexp + statistics) / (Eexp + statistics);
		}

		ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret;
	}

	double QuantumNumericalIntegrationT2dn2dmu2(int statistics, double T, double mu, double m, double deg)
	{
		if (statistics == 0)           return BoltzmannTdndmu(2, T, mu, m, deg);
		if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuT2dn2dmu2(T, mu, m, deg);
		if (statistics == -1 && mu > m) {
			printf("**WARNING** QuantumNumericalIntegrationT2dn2dmu2: Bose-Einstein condensation\n");
			return 0.;
		}

		double ret = 0.;
		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = lagx32[i];
			double Eexp = exp(sqrt(tx*tx + moverT*moverT) - muoverT);
			ret += lagw32[i] * T * tx * T * tx * T * (Eexp*Eexp - statistics * Eexp) / (Eexp + statistics) / (Eexp + statistics) / (Eexp + statistics);
		}

		ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret;
	}

	double QuantumNumericalIntegrationT3dn3dmu3(int statistics, double T, double mu, double m, double deg)
	{
		if (statistics == 0)           return BoltzmannTdndmu(3, T, mu, m, deg);
		if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuT3dn3dmu3(T, mu, m, deg);
		if (statistics == -1 && mu > m) {
			printf("**WARNING** QuantumNumericalIntegrationT3dn3dmu3: Bose-Einstein condensation\n");
			return 0.;
		}

		double ret = 0.;
		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = lagx32[i];
			double Eexp = exp(sqrt(tx*tx + moverT*moverT) - muoverT);
			ret += lagw32[i] * T * tx * T * tx * T * (Eexp*Eexp*Eexp - 4.*statistics*Eexp*Eexp + statistics*statistics*Eexp) / (Eexp + statistics) / (Eexp + statistics) / (Eexp + statistics) / (Eexp + statistics);
		}

		ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret;
	}

	double QuantumNumericalIntegrationTdndmu(int N, int statistics, double T, double mu, double m, double deg)
	{
		if (N < 0 || N>3) {
			printf("**WARNING** QuantumNumericalIntegrationTdndmu: N must be between 0 and 3!\n");
			exit(1);
		}
		if (N == 0)
			return QuantumNumericalIntegrationDensity(statistics, T, mu, m, deg);

		if (N == 1)
			return QuantumNumericalIntegrationT1dn1dmu1(statistics, T, mu, m, deg);

		if (N == 2)
			return QuantumNumericalIntegrationT2dn2dmu2(statistics, T, mu, m, deg);

		return QuantumNumericalIntegrationT3dn3dmu3(statistics, T, mu, m, deg);
	}

	double QuantumNumericalIntegrationChiN(int N, int statistics, double T, double mu, double m, double deg)
	{
		return QuantumNumericalIntegrationTdndmu(N - 1, statistics, T, mu, m, deg) / pow(T, 3) / xMath::GeVtoifm3();
	}

	double psi(double x)
	{
		double x2      = x*x;
		double tmpsqrt = sqrt(1. + x2);
		return (1. + x2 / 2.) * tmpsqrt - x2*x2 / 2. * log((1. + tmpsqrt) / x);
	}

	double psi2(double x)
	{
		double x2      = x*x;
		double tmpsqrt = sqrt(1. + x2);
		return 2. * tmpsqrt + 2. * x2 * log((1. + tmpsqrt) / x);
	}

	double FermiNumericalIntegrationLargeMuDensity(double T, double mu, double m, double deg)
	{
		if (mu <= m)
			return QuantumNumericalIntegrationDensity(1, T, mu, m, deg);

		double pf = sqrt(mu*mu - m*m);
		double ret1 = 0.;
		for (int i = 0; i < 32; i++) {
			ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf / (exp(-(sqrt(legx32[i] * legx32[i] * pf * pf + m*m) - mu) / T) + 1.);
		}

		double moverT  = m  / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = pf / T + lagx32[i];
			ret1 += lagw32[i] * T * tx * T * tx * T / (exp(sqrt(tx*tx + moverT*moverT) - muoverT) + 1.);
			//ret += lagw32[i] * lagx32[i] * lagx32[i] / (exp(tMu*(sqrt(lagx32[i]*lagx32[i]+tMass*tMass) - 20.)) + fStatistics);
			//std::cout << lagx32[i] * T * tMu << " " << sqrt(T*T*tMu*tMu*20.*20.-mass*mass*tMu*tMu*T*T*20.*20.) << "\n";
		}

		ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		double ret2 = 0.;
		ret2 += deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() * pf * pf * pf / 3.;

		return ret1 + ret2;
	}

	double FermiNumericalIntegrationLargeMuPressure(double T, double mu, double m, double deg)
	{
		if (mu <= m)
			return QuantumNumericalIntegrationPressure(1, T, mu, m, deg);

		double pf = sqrt(mu*mu - m*m);
		double ret1 = 0.;
		for (int i = 0; i < 32; i++) {
			double x2 = legx32[i] * pf * legx32[i] * pf;
			double E = sqrt(legx32[i] * legx32[i] * pf*pf + m*m);
			ret1 += -legw32[i] * pf * x2 * x2 / E / (exp(-(E - mu) / T) + 1.);
		}

		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = pf / T + lagx32[i];
			double x2 = tx * T * tx * T;
			double E = sqrt(tx*tx + moverT*moverT);
			ret1 += lagw32[i] * T * x2 * x2 / E / T / (exp(E - muoverT) + 1.);
		}

		ret1 *= deg / 6. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		double ret2 = 0.;
		ret2 += mu * pf * pf * pf;
		ret2 += -3. / 4. * pf * pf * pf * pf * psi(m / pf);
		ret2 *= deg / 6. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret1 + ret2;
	}

	double FermiNumericalIntegrationLargeMuEnergyDensity(double T, double mu, double m, double deg)
	{
		if (mu <= m)
			return QuantumNumericalIntegrationEnergyDensity(1, T, mu, m, deg);

		double pf = sqrt(mu*mu - m*m);
		double ret1 = 0.;
		for (int i = 0; i < 32; i++) {
			ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * sqrt(legx32[i] * legx32[i] * pf*pf + m*m) / (exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m*m) - mu) / T) + 1.);
		}

		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = pf / T + lagx32[i];
			ret1 += lagw32[i] * T * tx * T * tx * T * sqrt(tx*tx + moverT*moverT) * T / (exp(sqrt(tx*tx + moverT*moverT) - muoverT) + 1.);
		}

		ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		double ret2 = 0.;
		ret2 += deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() * pf * pf * pf * pf / 4. * psi(m / pf);

		return ret1 + ret2;
	}

	double FermiNumericalIntegrationLargeMuEntropyDensity(double T, double mu, double m, double deg)
	{
		return (FermiNumericalIntegrationLargeMuPressure(T, mu, m, deg) + FermiNumericalIntegrationLargeMuEnergyDensity(T, mu, m, deg) - mu * FermiNumericalIntegrationLargeMuDensity(T, mu, m, deg)) / T;
	}

	double FermiNumericalIntegrationLargeMuScalarDensity(double T, double mu, double m, double deg)
	{
		if (mu <= m)
			return QuantumNumericalIntegrationScalarDensity(1, T, mu, m, deg);

		double pf = sqrt(mu*mu - m*m);
		double ret1 = 0.;
		for (int i = 0; i < 32; i++) {
			ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * m / sqrt(legx32[i] * legx32[i] * pf*pf + m*m) / (exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m*m) - mu) / T) + 1.);
		}

		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = pf / T + lagx32[i];
			ret1 += lagw32[i] * T * tx * T * tx * T * moverT / sqrt(tx*tx + moverT*moverT) / (exp(sqrt(tx*tx + moverT*moverT) - muoverT) + 1.);
		}

		ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		double ret2 = 0.;
		ret2 += deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() * m * pf * (mu - pf / 4. * psi2(m / pf));

		return ret1 + ret2;
	}

	double FermiNumericalIntegrationLargeMuT1dn1dmu1(double T, double mu, double m, double deg)
	{
		if (mu <= m)
			return QuantumNumericalIntegrationT1dn1dmu1(1, T, mu, m, deg);

		double pf = sqrt(mu*mu - m*m);
		double ret1 = 0.;
		for (int i = 0; i < 32; i++) {
			double Eexp = exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m*m) - mu) / T);
			ret1 += legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * Eexp / (Eexp + 1.) / (Eexp + 1.);
		}

		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = pf / T + lagx32[i];
			double Eexp = exp(sqrt(tx*tx + moverT*moverT) - muoverT);
			ret1 += lagw32[i] * T * tx * T * tx * T * Eexp / (Eexp + 1.) / (Eexp + 1.);
		}

		ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret1;
	}

	double FermiNumericalIntegrationLargeMuT2dn2dmu2(double T, double mu, double m, double deg)
	{
		if (mu <= m)
			return QuantumNumericalIntegrationT2dn2dmu2(1, T, mu, m, deg);

		double pf = sqrt(mu*mu - m*m);
		double ret1 = 0.;
		for (int i = 0; i < 32; i++) {
			double Eexp = exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m*m) - mu) / T);
			ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * (Eexp*Eexp - Eexp) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.);
		}

		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = pf / T + lagx32[i];
			double Eexp = exp(sqrt(tx*tx + moverT*moverT) - muoverT);
			ret1 += lagw32[i] * T * tx * T * tx * T * (Eexp*Eexp - Eexp) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.);
		}

		ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret1;
	}

	double FermiNumericalIntegrationLargeMuT3dn3dmu3(double T, double mu, double m, double deg)
	{
		if (mu <= m)
			return QuantumNumericalIntegrationT3dn3dmu3(1, T, mu, m, deg);

		double pf = sqrt(mu*mu - m*m);
		double ret1 = 0.;
		for (int i = 0; i < 32; i++) {
			double Eexp = exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m*m) - mu) / T);
			ret1 += legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * (Eexp*Eexp*Eexp - 4.*Eexp*Eexp + Eexp) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.);
		}

		double moverT = m / T;
		double muoverT = mu / T;
		for (int i = 0; i < 32; i++) {
			double tx = pf / T + lagx32[i];
			double Eexp = exp(sqrt(tx*tx + moverT*moverT) - muoverT);
			ret1 += lagw32[i] * T * tx * T * tx * T * (Eexp*Eexp*Eexp - 4.*Eexp*Eexp + Eexp) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.);
		}

		ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

		return ret1;
	}

	double FermiNumericalIntegrationLargeMuTdndmu(int N, double T, double mu, double m, double deg)
	{
		if (N < 0 || N>3) {
			printf("**ERROR** FermiNumericalIntegrationLargeMuTdndmu: N < 0 or N > 3\n");
			exit(1);
		}
		if (N == 0)
			return FermiNumericalIntegrationLargeMuDensity(T, mu, m, deg);

		if (N == 1)
			return FermiNumericalIntegrationLargeMuT1dn1dmu1(T, mu, m, deg);

		if (N == 2)
			return FermiNumericalIntegrationLargeMuT2dn2dmu2(T, mu, m, deg);
		
		return FermiNumericalIntegrationLargeMuT3dn3dmu3(T, mu, m, deg);
	}

	double FermiNumericalIntegrationLargeMuChiN(int N, double T, double mu, double m, double deg)
	{
		return FermiNumericalIntegrationLargeMuTdndmu(N - 1, T, mu, m, deg) / pow(T, 3) / xMath::GeVtoifm3();
	}

	double IdealGasQuantity(Quantity quantity, QStatsCalculationType calctype, int statistics, double T, double mu, double m, double deg, int order)
	{
		if (statistics == 0) {
			if (quantity == ParticleDensity)
				return BoltzmannDensity(T, mu, m, deg);
			if (quantity == Pressure)
				return BoltzmannPressure(T, mu, m, deg);
			if (quantity == EnergyDensity)
				return BoltzmannEnergyDensity(T, mu, m, deg);
			if (quantity == EntropyDensity)
				return BoltzmannEntropyDensity(T, mu, m, deg);
			if (quantity == ScalarDensity)
				return BoltzmannScalarDensity(T, mu, m, deg);
			if (quantity == chi2)
				return BoltzmannChiN(2, T, mu, m, deg);
			if (quantity == chi3)
				return BoltzmannChiN(3, T, mu, m, deg);
			if (quantity == chi4)
				return BoltzmannChiN(4, T, mu, m, deg);
		}
		else {
			if (calctype == ClusterExpansion) {
				if (quantity == ParticleDensity)
					return QuantumClusterExpansionDensity(statistics, T, mu, m, deg, order);
				if (quantity == Pressure)
					return QuantumClusterExpansionPressure(statistics, T, mu, m, deg, order);
				if (quantity == EnergyDensity)
					return QuantumClusterExpansionEnergyDensity(statistics, T, mu, m, deg, order);
				if (quantity == EntropyDensity)
					return QuantumClusterExpansionEntropyDensity(statistics, T, mu, m, deg, order);
				if (quantity == ScalarDensity)
					return QuantumClusterExpansionScalarDensity(statistics, T, mu, m, deg, order);
				if (quantity == chi2)
					return QuantumClusterExpansionChiN(2, statistics, T, mu, m, deg, order);
				if (quantity == chi3)
					return QuantumClusterExpansionChiN(3, statistics, T, mu, m, deg, order);
				if (quantity == chi4)
					return QuantumClusterExpansionChiN(4, statistics, T, mu, m, deg, order);
			}
			else {
				if (quantity == ParticleDensity)
					return QuantumNumericalIntegrationDensity(statistics, T, mu, m, deg);
				if (quantity == Pressure)
					return QuantumNumericalIntegrationPressure(statistics, T, mu, m, deg);
				if (quantity == EnergyDensity)
					return QuantumNumericalIntegrationEnergyDensity(statistics, T, mu, m, deg);
				if (quantity == EntropyDensity)
					return QuantumNumericalIntegrationEntropyDensity(statistics, T, mu, m, deg);
				if (quantity == ScalarDensity)
					return QuantumNumericalIntegrationScalarDensity(statistics, T, mu, m, deg);
				if (quantity == chi2)
					return QuantumNumericalIntegrationChiN(2, statistics, T, mu, m, deg);
				if (quantity == chi3)
					return QuantumNumericalIntegrationChiN(3, statistics, T, mu, m, deg);
				if (quantity == chi4)
					return QuantumNumericalIntegrationChiN(4, statistics, T, mu, m, deg);
			}
		}
		printf("**WARNING** IdealGasFunctions::IdealGasQuantity: Unknown quantity\n");
		return 0.;
	}




}
