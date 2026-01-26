/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <string>
#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelParameters.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGVDW/ThermalModelVDW.h"
#include "HRGRealGas/ThermalModelRealGas.h"
#include "HRGRealGas/ExcludedVolumeModelsMulti.h"
#include "HRGRealGas/MeanFieldModelsMulti.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  // Helper function to check if a value is NaN or Inf
  bool isValidNumber(double x) {
    return !std::isnan(x) && !std::isinf(x);
  }

  // Helper function to set up QvdW-HRG model with standard baryon-baryon interactions
  void SetupQvdWHRG(ThermalModelVDW& model, double a = 0.329, double b = 3.42) {
    // Set QvdW interactions for baryon-baryon and antibaryon-antibaryon pairs
    // as in https://arxiv.org/abs/1609.03975
    SetVDWHRGInteractionParameters(&model, a, b);
  }

  //============================================================================
  // QvdW-HRG Stress Tests
  //============================================================================

  class QvdWHRGStressTest : public ::testing::Test {
  protected:
    void SetUp() override {
      // Use PDG2020 list
      std::string inputDir = std::string(ThermalFIST_INPUT_FOLDER);
      TPS = new ThermalParticleSystem(inputDir + "/list/PDG2020/list.dat", inputDir + "/list/PDG2020/decays.dat");
    }

    void TearDown() override {
      delete TPS;
    }

    ThermalParticleSystem* TPS;
  };

  // Test case that previously caused NaN: T=400 MeV, muB=350 MeV
  TEST_F(QvdWHRGStressTest, HighTemperatureHighMuB) {
    ThermalModelParameters params;
    params.T = 0.400;    // 400 MeV
    params.muB = 0.350;  // 350 MeV
    params.muS = 0.0;
    params.muQ = 0.0;
    params.muC = 0.0;
    params.V = 1.0;

    ThermalModelVDW model(TPS, params);
    model.SetStatistics(false);  // Maxwell-Boltzmann
    SetupQvdWHRG(model);

    model.CalculatePrimordialDensities();

    // Check thermodynamic quantities
    double pressure = model.CalculatePressure();
    double energyDensity = model.CalculateEnergyDensity();
    double entropyDensity = model.CalculateEntropyDensity();

    EXPECT_TRUE(isValidNumber(pressure)) << "Pressure is NaN/Inf";
    EXPECT_TRUE(isValidNumber(energyDensity)) << "Energy density is NaN/Inf";
    EXPECT_TRUE(isValidNumber(entropyDensity)) << "Entropy density is NaN/Inf";

    // Check physical constraints
    EXPECT_GT(pressure, 0.0) << "Pressure should be positive";
    EXPECT_GT(energyDensity, 0.0) << "Energy density should be positive";
    EXPECT_GT(entropyDensity, 0.0) << "Entropy density should be positive";

    // Verify solution succeeded
    EXPECT_TRUE(model.IsLastSolutionOK()) << "Broyden solver did not converge";
  }

  // Scan over a range of T and muB values that may be challenging
  TEST_F(QvdWHRGStressTest, TemperatureMuBScan) {
    std::vector<double> temperatures = {0.100, 0.150, 0.200, 0.300, 0.400, 0.500};
    std::vector<double> muBValues = {0.0, 0.200, 0.350, 0.500, 0.700};

    for (double T : temperatures) {
      for (double muB : muBValues) {
        ThermalModelParameters params;
        params.T = T;
        params.muB = muB;
        params.muS = 0.0;
        params.muQ = 0.0;
        params.muC = 0.0;
        params.V = 1.0;

        ThermalModelVDW model(TPS, params);
        model.SetStatistics(false);  // Maxwell-Boltzmann
        SetupQvdWHRG(model);

        model.CalculatePrimordialDensities();

        double pressure = model.CalculatePressure();
        double energyDensity = model.CalculateEnergyDensity();
        double entropyDensity = model.CalculateEntropyDensity();

        EXPECT_TRUE(isValidNumber(pressure))
          << "Pressure is NaN/Inf at T=" << T << " GeV, muB=" << muB << " GeV";
        EXPECT_TRUE(isValidNumber(energyDensity))
          << "Energy density is NaN/Inf at T=" << T << " GeV, muB=" << muB << " GeV";
        EXPECT_TRUE(isValidNumber(entropyDensity))
          << "Entropy density is NaN/Inf at T=" << T << " GeV, muB=" << muB << " GeV";

        EXPECT_TRUE(model.IsLastSolutionOK())
          << "Broyden solver failed at T=" << T << " GeV, muB=" << muB << " GeV";
      }
    }
  }

  // Test with quantum statistics at various T and muB
  TEST_F(QvdWHRGStressTest, QuantumStatisticsScan) {
    std::vector<double> temperatures = {0.100, 0.140, 0.160, 0.200};
    std::vector<double> muBValues = {0.0, 0.100, 0.200, 0.300};

    for (double T : temperatures) {
      for (double muB : muBValues) {
        ThermalModelParameters params;
        params.T = T;
        params.muB = muB;
        params.muS = 0.0;
        params.muQ = 0.0;
        params.muC = 0.0;
        params.V = 1.0;

        ThermalModelVDW model(TPS, params);
        model.SetStatistics(true);  // Quantum statistics
        SetupQvdWHRG(model);

        model.CalculatePrimordialDensities();

        double pressure = model.CalculatePressure();
        double energyDensity = model.CalculateEnergyDensity();

        EXPECT_TRUE(isValidNumber(pressure))
          << "Pressure is NaN/Inf with quantum stats at T=" << T << " GeV, muB=" << muB << " GeV";
        EXPECT_TRUE(isValidNumber(energyDensity))
          << "Energy density is NaN/Inf with quantum stats at T=" << T << " GeV, muB=" << muB << " GeV";
        EXPECT_TRUE(model.IsLastSolutionOK())
          << "Broyden solver failed with quantum stats at T=" << T << " GeV, muB=" << muB << " GeV";
      }
    }
  }

  // Test at very high temperature
  TEST_F(QvdWHRGStressTest, VeryHighTemperature) {
    ThermalModelParameters params;
    params.T = 0.800;  // 800 MeV - very high temperature
    params.muB = 0.0;
    params.muS = 0.0;
    params.muQ = 0.0;
    params.muC = 0.0;
    params.V = 1.0;

    ThermalModelVDW model(TPS, params);
    model.SetStatistics(false);
    SetupQvdWHRG(model);

    model.CalculatePrimordialDensities();

    double pressure = model.CalculatePressure();
    EXPECT_TRUE(isValidNumber(pressure)) << "Pressure is NaN/Inf at very high T";
    EXPECT_TRUE(model.IsLastSolutionOK()) << "Broyden solver failed at very high T";
  }

  // Test near nuclear matter density
  TEST_F(QvdWHRGStressTest, NuclearMatterRegion) {
    ThermalModelParameters params;
    params.T = 0.010;  // 10 MeV - low temperature
    params.muB = 0.920;  // Close to nucleon mass
    params.muS = 0.0;
    params.muQ = 0.0;
    params.muC = 0.0;
    params.V = 1.0;

    ThermalModelVDW model(TPS, params);
    model.SetStatistics(false);
    SetupQvdWHRG(model);

    model.CalculatePrimordialDensities();

    double pressure = model.CalculatePressure();
    EXPECT_TRUE(isValidNumber(pressure)) << "Pressure is NaN/Inf in nuclear matter region";
    // Note: solution may not always converge in this challenging region,
    // but we at least check for NaN
  }

  // Test with different VDW parameters
  TEST_F(QvdWHRGStressTest, DifferentVDWParameters) {
    std::vector<std::pair<double, double>> abPairs = {
      {0.0, 0.0},      // No interactions
      {0.329, 3.42},   // Standard parameters
      {0.5, 5.0},      // Stronger interactions
      {1.0, 8.0},      // Very strong interactions
    };

    for (const auto& ab : abPairs) {
      double a = ab.first;
      double b = ab.second;

      ThermalModelParameters params;
      params.T = 0.160;
      params.muB = 0.200;
      params.muS = 0.0;
      params.muQ = 0.0;
      params.muC = 0.0;
      params.V = 1.0;

      ThermalModelVDW model(TPS, params);
      model.SetStatistics(false);
      SetupQvdWHRG(model, a, b);

      model.CalculatePrimordialDensities();

      double pressure = model.CalculatePressure();
      EXPECT_TRUE(isValidNumber(pressure))
        << "Pressure is NaN/Inf with a=" << a << ", b=" << b;
    }
  }

  //============================================================================
  // Real Gas HRG Stress Tests
  //============================================================================

  class RealGasHRGStressTest : public ::testing::Test {
  protected:
    void SetUp() override {
      // Use PDG2020 list
      std::string inputDir = std::string(ThermalFIST_INPUT_FOLDER);
      TPS = new ThermalParticleSystem(inputDir + "/list/PDG2020/list.dat", inputDir + "/list/PDG2020/decays.dat");
    }

    void TearDown() override {
      delete TPS;
    }

    // Helper to create diagonal excluded volume parameters for baryons
    std::vector<double> CreateBaryonEVParams(double radius) {
      std::vector<double> b(TPS->ComponentsNumber(), 0.);
      double v0 = (4./3.) * 3.14159265 * radius * radius * radius;  // eigenvolume
      for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
        if (TPS->Particle(i).BaryonCharge() != 0)
          b[i] = v0;
      }
      return b;
    }

    ThermalParticleSystem* TPS;
  };

  // Test RealGas model at high T and muB with diagonal EV
  TEST_F(RealGasHRGStressTest, HighTemperatureHighMuB) {
    ThermalModelParameters params;
    params.T = 0.400;    // 400 MeV
    params.muB = 0.350;  // 350 MeV
    params.muS = 0.0;
    params.muQ = 0.0;
    params.muC = 0.0;
    params.V = 1.0;

    ThermalModelRealGas model(TPS, params);
    model.SetStatistics(false);  // Maxwell-Boltzmann

    // Set up diagonal excluded volume model for baryons
    std::vector<double> b = CreateBaryonEVParams(0.3);  // 0.3 fm radius
    ExcludedVolumeModelDiagonalVDW* evmod = new ExcludedVolumeModelDiagonalVDW(b);
    model.SetExcludedVolumeModel(evmod);

    model.CalculatePrimordialDensities();

    double pressure = model.CalculatePressure();
    double energyDensity = model.CalculateEnergyDensity();
    double entropyDensity = model.CalculateEntropyDensity();

    EXPECT_TRUE(isValidNumber(pressure)) << "Pressure is NaN/Inf";
    EXPECT_TRUE(isValidNumber(energyDensity)) << "Energy density is NaN/Inf";
    EXPECT_TRUE(isValidNumber(entropyDensity)) << "Entropy density is NaN/Inf";

    EXPECT_GT(pressure, 0.0) << "Pressure should be positive";
    EXPECT_TRUE(model.IsLastSolutionOK()) << "Broyden solver did not converge";
    // Note: model takes ownership of evmod and deletes it in destructor
  }

  // Scan over T and muB for RealGas model
  TEST_F(RealGasHRGStressTest, TemperatureMuBScan) {
    std::vector<double> temperatures = {0.100, 0.150, 0.200, 0.300, 0.400};
    std::vector<double> muBValues = {0.0, 0.200, 0.350, 0.500};

    std::vector<double> b = CreateBaryonEVParams(0.3);

    for (double T : temperatures) {
      for (double muB : muBValues) {
        ThermalModelParameters params;
        params.T = T;
        params.muB = muB;
        params.muS = 0.0;
        params.muQ = 0.0;
        params.muC = 0.0;
        params.V = 1.0;

        ThermalModelRealGas model(TPS, params);
        model.SetStatistics(false);

        ExcludedVolumeModelDiagonalVDW* evmod = new ExcludedVolumeModelDiagonalVDW(b);
        model.SetExcludedVolumeModel(evmod);

        model.CalculatePrimordialDensities();

        double pressure = model.CalculatePressure();
        double energyDensity = model.CalculateEnergyDensity();

        EXPECT_TRUE(isValidNumber(pressure))
          << "RealGas: Pressure is NaN/Inf at T=" << T << " GeV, muB=" << muB << " GeV";
        EXPECT_TRUE(isValidNumber(energyDensity))
          << "RealGas: Energy density is NaN/Inf at T=" << T << " GeV, muB=" << muB << " GeV";
        EXPECT_TRUE(model.IsLastSolutionOK())
          << "RealGas: Broyden solver failed at T=" << T << " GeV, muB=" << muB << " GeV";
        // Note: model takes ownership of evmod and deletes it in destructor
      }
    }
  }

  // Test RealGas with quantum statistics
  TEST_F(RealGasHRGStressTest, QuantumStatisticsScan) {
    std::vector<double> temperatures = {0.100, 0.140, 0.160, 0.200};
    std::vector<double> muBValues = {0.0, 0.100, 0.200};

    std::vector<double> b = CreateBaryonEVParams(0.3);

    for (double T : temperatures) {
      for (double muB : muBValues) {
        ThermalModelParameters params;
        params.T = T;
        params.muB = muB;
        params.muS = 0.0;
        params.muQ = 0.0;
        params.muC = 0.0;
        params.V = 1.0;

        ThermalModelRealGas model(TPS, params);
        model.SetStatistics(true);  // Quantum statistics

        ExcludedVolumeModelDiagonalVDW* evmod = new ExcludedVolumeModelDiagonalVDW(b);
        model.SetExcludedVolumeModel(evmod);

        model.CalculatePrimordialDensities();

        double pressure = model.CalculatePressure();
        double energyDensity = model.CalculateEnergyDensity();

        EXPECT_TRUE(isValidNumber(pressure))
          << "RealGas QS: Pressure is NaN/Inf at T=" << T << " GeV, muB=" << muB << " GeV";
        EXPECT_TRUE(isValidNumber(energyDensity))
          << "RealGas QS: Energy density is NaN/Inf at T=" << T << " GeV, muB=" << muB << " GeV";
        EXPECT_TRUE(model.IsLastSolutionOK())
          << "RealGas QS: Broyden solver failed at T=" << T << " GeV, muB=" << muB << " GeV";
        // Note: model takes ownership of evmod and deletes it in destructor
      }
    }
  }

  // Test RealGas with crossterms excluded volume
  TEST_F(RealGasHRGStressTest, CrosstermsEV) {
    ThermalModelParameters params;
    params.T = 0.160;
    params.muB = 0.200;
    params.muS = 0.0;
    params.muQ = 0.0;
    params.muC = 0.0;
    params.V = 1.0;

    ThermalModelRealGas model(TPS, params);
    model.SetStatistics(false);

    // Create crossterms excluded volume matrix for baryons
    int N = TPS->ComponentsNumber();
    double radius = 0.3;
    double v0 = (4./3.) * 3.14159265 * radius * radius * radius;
    std::vector<std::vector<double>> bij(N, std::vector<double>(N, 0.));
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        if (TPS->Particle(i).BaryonCharge() != 0 && TPS->Particle(j).BaryonCharge() != 0)
          bij[i][j] = v0;
      }
    }

    ExcludedVolumeModelCrosstermsVDW* evmod = new ExcludedVolumeModelCrosstermsVDW(bij);
    model.SetExcludedVolumeModel(evmod);

    model.CalculatePrimordialDensities();

    double pressure = model.CalculatePressure();
    double energyDensity = model.CalculateEnergyDensity();

    EXPECT_TRUE(isValidNumber(pressure)) << "Crossterms EV: Pressure is NaN/Inf";
    EXPECT_TRUE(isValidNumber(energyDensity)) << "Crossterms EV: Energy density is NaN/Inf";
    EXPECT_TRUE(model.IsLastSolutionOK()) << "Crossterms EV: Broyden solver failed";
    // Note: model takes ownership of evmod and deletes it in destructor
  }

  // Test RealGas with mean field
  TEST_F(RealGasHRGStressTest, MeanFieldModel) {
    ThermalModelParameters params;
    params.T = 0.160;
    params.muB = 0.200;
    params.muS = 0.0;
    params.muQ = 0.0;
    params.muC = 0.0;
    params.V = 1.0;

    ThermalModelRealGas model(TPS, params);
    model.SetStatistics(false);

    // Set up diagonal excluded volume
    std::vector<double> b = CreateBaryonEVParams(0.3);
    ExcludedVolumeModelDiagonalVDW* evmod = new ExcludedVolumeModelDiagonalVDW(b);
    model.SetExcludedVolumeModel(evmod);

    // Set up mean field model (attractive interactions for baryons)
    int N = TPS->ComponentsNumber();
    std::vector<std::vector<double>> aij(N, std::vector<double>(N, 0.));
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        if (TPS->Particle(i).BaryonCharge() != 0 && TPS->Particle(j).BaryonCharge() != 0) {
          if (TPS->Particle(i).BaryonCharge() * TPS->Particle(j).BaryonCharge() > 0)
            aij[i][j] = 0.329;  // Same sign baryons attract
        }
      }
    }
    MeanFieldModelMultiVDW* mfmod = new MeanFieldModelMultiVDW(aij);
    model.SetMeanFieldModel(mfmod);

    model.CalculatePrimordialDensities();

    double pressure = model.CalculatePressure();
    double energyDensity = model.CalculateEnergyDensity();

    EXPECT_TRUE(isValidNumber(pressure)) << "Mean field: Pressure is NaN/Inf";
    EXPECT_TRUE(isValidNumber(energyDensity)) << "Mean field: Energy density is NaN/Inf";
    EXPECT_TRUE(model.IsLastSolutionOK()) << "Mean field: Broyden solver failed";
    // Note: model takes ownership of evmod and mfmod and deletes them in destructor
  }

  // Test RealGas at very high temperature
  TEST_F(RealGasHRGStressTest, VeryHighTemperature) {
    ThermalModelParameters params;
    params.T = 0.800;  // 800 MeV
    params.muB = 0.0;
    params.muS = 0.0;
    params.muQ = 0.0;
    params.muC = 0.0;
    params.V = 1.0;

    ThermalModelRealGas model(TPS, params);
    model.SetStatistics(false);

    std::vector<double> b = CreateBaryonEVParams(0.3);
    ExcludedVolumeModelDiagonalVDW* evmod = new ExcludedVolumeModelDiagonalVDW(b);
    model.SetExcludedVolumeModel(evmod);

    model.CalculatePrimordialDensities();

    double pressure = model.CalculatePressure();
    EXPECT_TRUE(isValidNumber(pressure)) << "RealGas: Pressure is NaN/Inf at very high T";
    EXPECT_TRUE(model.IsLastSolutionOK()) << "RealGas: Broyden solver failed at very high T";
    // Note: model takes ownership of evmod and deletes it in destructor
  }

}
