/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <limits.h>
#include "HRGBase/xMath.h"
#include "HRGBase/IdealGasFunctions.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

	TEST(BoltzmannTest, Density) {
		// For zero degeneracy one should get zero density
		EXPECT_EQ(IdealGasFunctions::BoltzmannDensity(0.100, 0.100, 0.938, 0), 0.);
		// Cross-checking accuracy of the Boltzmann distribution functions, using Mathematica as a reference
		// Relative error of at least 10^-5
		double accuracy = 1.e-5;

		// Nucleons
		double MathematicaRef = 0.0000563767;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.160, 0.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.68958e-4;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.160, 0.250, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.0000118172;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.160, -0.250, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.0156317;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.160, 0.900, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.03326e-7;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.160, -0.900, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 15.128;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.160, 2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.10096e-10;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.160, -2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.02516e-11;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.050, 0.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 4.76692e6;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.050, 2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 8.60358e-29;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.050, -2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 4.31427e-45;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.010, 0.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 3.11748e42;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.010, 2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 5.9705e-132;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.010, -2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 3.17596e-109;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.001, 0.700, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.294942440341e-22;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.001, 0.900, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 4.457756771677e108;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.001, 1.200, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);

		// pions
		MathematicaRef = 3.2752e-24;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.001, 0.100, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.36665e63;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.001, 0.300, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 1.71014e150;
		EXPECT_LT(abs(IdealGasFunctions::BoltzmannDensity(0.001, 0.500, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
	}


	TEST(QSClusterExpansionTest, Density) {
		// For zero degeneracy one should get zero density
		EXPECT_EQ(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.100, 0.100, 0.938, 0, 5), 0.);
		// Testing Cluster Expansion Quantum Statistics calculations
		// Requiring accuracy of at least 1% compared to the Mathematica reference
		double accuracy = 1.e-2;

		// Nucleons
		double MathematicaRef = 4.3952e-196;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.001, 0.500, 0.938, 4, 5) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.29494e-22;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.001, 0.900, 0.938, 4, 5) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.45219e-9;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.001, 0.930, 0.938, 4, 5) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.45219e-9;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.001, 0.930, 0.938, 4, 5) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.0000922762;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.010, 0.930, 0.938, 4, 5) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.00163886;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.050, 0.920, 0.938, 4, 5) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.0185169;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.200, 0.880, 0.938, 4, 5) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.116894;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(1, 0.500, 0.800, 0.938, 4, 5) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);

		// pions
		MathematicaRef = 1.2184e-67;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(-1, 0.001, 0.000, 0.138, 1, 10) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.0000749895;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(-1, 0.100, 0.000, 0.138, 1, 10) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.120262;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(-1, 1.000, 0.000, 0.138, 1, 10) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 121.767;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(-1, 10.000, 0.000, 0.138, 1, 10) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		
		MathematicaRef = 3.50045e-11;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(-1, 0.001, 0.130, 0.138, 1, 10) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.0000812976;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(-1, 0.050, 0.130, 0.138, 1, 10) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.140011;
		EXPECT_LT(abs(IdealGasFunctions::QuantumClusterExpansionDensity(-1, 1.000, 0.110, 0.138, 1, 10) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
	}


	TEST(QSQuadratures, Density) {
		// For zero degeneracy one should get zero density
		EXPECT_EQ(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.100, 0.100, 0.938, 0), 0.);
		// Testing Quantum Statistics with Numerical Integraion (Quadratures) calculations
		// Requiring accuracy of at least 1% compared to the Mathematica reference
		double accuracy = 1.e-2;

		// Nucleons
		double MathematicaRef = 4.3952e-196;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.001, 0.500, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.29494e-22;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.001, 0.900, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.45219e-9;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.001, 0.930, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 2.45219e-9;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.001, 0.930, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.00281445;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.001, 1.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		
		MathematicaRef = 0.00292343;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.010, 1.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.372418;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.010, 2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);

		MathematicaRef = 0.0111384;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.100, 1.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.385726;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.100, 2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);

		MathematicaRef = 0.771986;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 1.000, 1.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 1.74909;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 1.000, 2.000, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);


		MathematicaRef = 0.0000922762;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.010, 0.930, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.00163886;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.050, 0.920, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.0185169;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.200, 0.880, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.116894;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(1, 0.500, 0.800, 0.938, 4) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);

		// pions
		MathematicaRef = 1.2184e-67;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(-1, 0.001, 0.000, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.0000749895;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(-1, 0.100, 0.000, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.120262;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(-1, 1.000, 0.000, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 121.767;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(-1, 10.000, 0.000, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);

		MathematicaRef = 4.46552e-8;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(-1, 0.001, 0.137, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.000113812;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(-1, 0.050, 0.137, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
		MathematicaRef = 0.146728;
		EXPECT_LT(abs(IdealGasFunctions::QuantumNumericalIntegrationDensity(-1, 1.000, 0.137, 0.138, 1) / xMath::GeVtoifm3() - MathematicaRef) / MathematicaRef, accuracy);
	}

}