/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef IDEALGASFUNCTIONS_H
#define IDEALGASFUNCTIONS_H

namespace thermalfist {

  namespace IdealGasFunctions {
    enum Quantity { ParticleDensity, EnergyDensity, EntropyDensity, Pressure, chi2, chi3, chi4, ScalarDensity };
    enum QStatsCalculationType { ClusterExpansion, Quadratures };

    double BoltzmannDensity(double T, double mu, double m, double deg);
    double BoltzmannPressure(double T, double mu, double m, double deg);
    double BoltzmannEnergyDensity(double T, double mu, double m, double deg);
    double BoltzmannEntropyDensity(double T, double mu, double m, double deg);
    double BoltzmannScalarDensity(double T, double mu, double m, double deg);  // TODO: Check for correctness
    double BoltzmannTdndmu(int N, double T, double mu, double m, double deg);
    double BoltzmannChiN(int N, double T, double mu, double m, double deg);

    double QuantumClusterExpansionDensity(int statistics, double T, double mu, double m, double deg, int order = 1);
    double QuantumClusterExpansionPressure(int statistics, double T, double mu, double m, double deg, int order = 1);
    double QuantumClusterExpansionEnergyDensity(int statistics, double T, double mu, double m, double deg, int order = 1);  // TODO: Check for correctness
    double QuantumClusterExpansionEntropyDensity(int statistics, double T, double mu, double m, double deg, int order = 1);
    double QuantumClusterExpansionScalarDensity(int statistics, double T, double mu, double m, double deg, int order = 1);   // TODO: Check for correctness
    double QuantumClusterExpansionTdndmu(int N, int statistics, double T, double mu, double m, double deg, int order = 1);
    double QuantumClusterExpansionChiN(int N, int statistics, double T, double mu, double m, double deg, int order = 1);

    double QuantumNumericalIntegrationDensity(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationPressure(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationEnergyDensity(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationEntropyDensity(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationScalarDensity(int statistics, double T, double mu, double m, double deg);  // TODO: Check for correctness
    double QuantumNumericalIntegrationT1dn1dmu1(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationT2dn2dmu2(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationT3dn3dmu3(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationTdndmu(int N, int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationChiN(int N, int statistics, double T, double mu, double m, double deg);

    double psi(double x);
    double psi2(double x);
    double FermiNumericalIntegrationLargeMuDensity(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuPressure(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuEnergyDensity(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuEntropyDensity(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuScalarDensity(double T, double mu, double m, double deg);  // TODO: Check for correctness
    double FermiNumericalIntegrationLargeMuT1dn1dmu1(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuT2dn2dmu2(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuT3dn3dmu3(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuTdndmu(int N, double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuChiN(int N, double T, double mu, double m, double deg);

    double IdealGasQuantity(Quantity quantity, QStatsCalculationType calctype, int statistics, double T, double mu, double m, double deg, int order = 1);
  }

} // namespace thermalfist

#endif // IDEALGASFUNCTIONS_H
