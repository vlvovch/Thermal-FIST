/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

/**
 * Unit tests for Effective Mass Model (EMM) temperature derivatives
 * in ThermalModelRealGas with differentiated EV/MF parameters.
 *
 * Tests:
 *  - Thermodynamic identity: e + P = T*s + sum(mu_i * n_i)
 *  - Temperature derivatives: ds/dT, de/dT vs numerical finite differences
 *  - Chemical potential derivatives: nB = dP/dmuB, nQ = dP/dmuQ
 *  - Both normal and BEC phases
 */

#include <cmath>
#include <vector>
#include <string>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelParameters.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGRealGas/ThermalModelRealGas.h"
#include "HRGRealGas/ExcludedVolumeModelsMulti.h"
#include "HRGRealGas/MeanFieldModelsMulti.h"
#include "HRGEV/ExcludedVolumeHelper.h"
#include "CosmicEos/EffectiveMassModel.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  // Accuracy threshold for all comparisons
  const double kAccuracy = 1.e-6;

  // EV/MF parameters
  const double kRMeson  = 0.15;  // fm
  const double kRBaryon = 0.4;   // fm
  const double kAMeson  = 0.05;  // GeV fm^3
  const double kABaryon = 0.3;   // GeV fm^3

  // Helper: set up a RealGas model with given EV/MF and optionally EMM on pions
  void SetupModel(ThermalModelRealGas& model, ThermalParticleSystem& tps,
                  const std::vector<std::vector<double>>& bij,
                  const std::vector<std::vector<double>>& aij,
                  bool useEMM,
                  double T, double muB, double muQ, double muS) {
    model.SetExcludedVolumeModel(new ExcludedVolumeModelCrosstermsVDW(bij));
    model.SetMeanFieldModel(new MeanFieldModelMultiVDW(aij));

    if (useEMM) {
      int pdgs[] = {211, 111, -211};
      for (int pdg : pdgs) {
        tps.ParticleByPDG(pdg).SetGeneralizedDensity(
          new EffectiveMassModel(
            tps.ParticleByPDG(pdg),
            new EMMFieldPressureChPT(tps.ParticleByPDG(pdg).Mass(), 0.133)));
      }
    }

    ThermalModelParameters params;
    params.T = T;
    params.muB = muB;
    params.muQ = muQ;
    params.muS = muS;
    params.muC = 0.;
    model.SetParameters(params);
    model.SetStatistics(true);
    model.ConstrainMuB(false);
    model.ConstrainMuQ(false);
    model.ConstrainMuS(false);
    model.ConstrainMuC(false);
  }

  // Helper: compute relative difference (or absolute if both values are tiny)
  double RelDiff(double val, double ref) {
    double scale = std::max(std::abs(val), std::abs(ref));
    if (scale > 1.e-10)
      return std::abs(val - ref) / scale;
    return std::abs(val - ref);  // both essentially zero
  }

  //============================================================================
  // Test fixture: builds the bij and aij matrices once
  //============================================================================
  class EMMDerivativesTest : public ::testing::Test {
  protected:
    void SetUp() override {
      std::string listpath = std::string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list.dat";
      TPS = new ThermalParticleSystem(listpath);
      int N = TPS->Particles().size();

      bij.resize(N, std::vector<double>(N, 0.));
      aij.resize(N, std::vector<double>(N, 0.));

      for (int i = 0; i < N; ++i) {
        int Bi = TPS->Particles()[i].BaryonCharge();
        double ri = (Bi == 0) ? kRMeson : kRBaryon;
        double ai_coeff = (Bi == 0) ? kAMeson : kABaryon;
        for (int j = 0; j < N; ++j) {
          int Bj = TPS->Particles()[j].BaryonCharge();
          double rj = (Bj == 0) ? kRMeson : kRBaryon;
          double aj_coeff = (Bj == 0) ? kAMeson : kABaryon;
          bij[i][j] = CuteHRGHelper::brr(ri, rj);
          aij[i][j] = std::sqrt(ai_coeff * aj_coeff);
        }
      }
    }

    void TearDown() override {
      delete TPS;
    }

    // Run all checks at given (T, muB, muQ, muS) with or without EMM
    void RunChecks(double T, double muB, double muQ, double muS, bool useEMM) {
      ThermalParticleSystem TPScopy(*TPS);
      ThermalModelRealGas model(&TPScopy);
      SetupModel(model, TPScopy, bij, aij, useEMM, T, muB, muQ, muS);
      model.CalculateDensities();

      // --- Thermodynamic identity: e + P = T*s + sum(mu_i * n_i) ---
      double P = model.CalculatePressure();
      double e = model.CalculateEnergyDensity();
      double s = model.CalculateEntropyDensity();

      double sum_mu_n = 0.;
      for (size_t i = 0; i < model.TPS()->Particles().size(); ++i) {
        double mu_i = model.ChemicalPotential(i);
        double n_i  = model.Densities()[i];
        sum_mu_n += mu_i * n_i;
      }

      double lhs = e + P;
      double rhs = T * s + sum_mu_n;
      EXPECT_LT(RelDiff(rhs, lhs), kAccuracy)
        << "Thermodynamic identity e+P = Ts+sum(mu*n) failed"
        << " at T=" << T << ", muB=" << muB << ", muQ=" << muQ
        << ", EMM=" << useEMM;

      // --- Temperature derivatives via CalculateTemperatureDerivatives ---
      model.CalculateTemperatureDerivatives();
      double dsdT_model = model.CalculateEntropyDensityDerivativeT();
      double dedT_model = model.CalculateEnergyDensityDerivativeT();

      // --- Numerical finite differences for T-derivatives ---
      double dT = 1.e-4 * T;

      ThermalParticleSystem TPSp(*TPS);
      ThermalModelRealGas model_p(&TPSp);
      SetupModel(model_p, TPSp, bij, aij, useEMM, T + 0.5 * dT, muB, muQ, muS);
      model_p.CalculateDensities();
      double s_p = model_p.CalculateEntropyDensity();
      double e_p = model_p.CalculateEnergyDensity();

      ThermalParticleSystem TPSm(*TPS);
      ThermalModelRealGas model_m(&TPSm);
      SetupModel(model_m, TPSm, bij, aij, useEMM, T - 0.5 * dT, muB, muQ, muS);
      model_m.CalculateDensities();
      double s_m = model_m.CalculateEntropyDensity();
      double e_m = model_m.CalculateEnergyDensity();

      double dsdT_fd = (s_p - s_m) / dT;
      double dedT_fd = (e_p - e_m) / dT;

      EXPECT_LT(RelDiff(dsdT_model, dsdT_fd), kAccuracy)
        << "ds/dT mismatch at T=" << T << ", muB=" << muB << ", muQ=" << muQ
        << ", EMM=" << useEMM
        << " (model=" << dsdT_model << ", FD=" << dsdT_fd << ")";

      EXPECT_LT(RelDiff(dedT_model, dedT_fd), kAccuracy)
        << "de/dT mismatch at T=" << T << ", muB=" << muB << ", muQ=" << muQ
        << ", EMM=" << useEMM
        << " (model=" << dedT_model << ", FD=" << dedT_fd << ")";

      // --- dP/dmuB = nB ---
      double nB_model = model.CalculateBaryonDensity();
      double dmu = 1.e-4 * T;

      {
        ThermalParticleSystem TPS_Bp(*TPS), TPS_Bm(*TPS);
        ThermalModelRealGas m_Bp(&TPS_Bp), m_Bm(&TPS_Bm);
        SetupModel(m_Bp, TPS_Bp, bij, aij, useEMM, T, muB + 0.5 * dmu, muQ, muS);
        SetupModel(m_Bm, TPS_Bm, bij, aij, useEMM, T, muB - 0.5 * dmu, muQ, muS);
        m_Bp.CalculateDensities();
        m_Bm.CalculateDensities();
        double dPdmuB_fd = (m_Bp.CalculatePressure() - m_Bm.CalculatePressure()) / dmu;

        EXPECT_LT(RelDiff(nB_model, dPdmuB_fd), kAccuracy)
          << "nB != dP/dmuB at T=" << T << ", muB=" << muB << ", muQ=" << muQ
          << ", EMM=" << useEMM
          << " (nB=" << nB_model << ", dP/dmuB=" << dPdmuB_fd << ")";
      }

      // --- dP/dmuQ = nQ ---
      double nQ_model = model.CalculateChargeDensity();

      {
        ThermalParticleSystem TPS_Qp(*TPS), TPS_Qm(*TPS);
        ThermalModelRealGas m_Qp(&TPS_Qp), m_Qm(&TPS_Qm);
        SetupModel(m_Qp, TPS_Qp, bij, aij, useEMM, T, muB, muQ + 0.5 * dmu, muS);
        SetupModel(m_Qm, TPS_Qm, bij, aij, useEMM, T, muB, muQ - 0.5 * dmu, muS);
        m_Qp.CalculateDensities();
        m_Qm.CalculateDensities();
        double dPdmuQ_fd = (m_Qp.CalculatePressure() - m_Qm.CalculatePressure()) / dmu;

        EXPECT_LT(RelDiff(nQ_model, dPdmuQ_fd), kAccuracy)
          << "nQ != dP/dmuQ at T=" << T << ", muB=" << muB << ", muQ=" << muQ
          << ", EMM=" << useEMM
          << " (nQ=" << nQ_model << ", dP/dmuQ=" << dPdmuQ_fd << ")";
      }
    }

    ThermalParticleSystem* TPS;
    std::vector<std::vector<double>> bij;
    std::vector<std::vector<double>> aij;
  };

  //============================================================================
  // Normal phase tests (no EMM)
  //============================================================================

  TEST_F(EMMDerivativesTest, NormalPhase_NoEMM_T150_muB0) {
    RunChecks(0.150, 0.0, 0.0, 0.0, false);
  }

  TEST_F(EMMDerivativesTest, NormalPhase_NoEMM_T150_muB300) {
    RunChecks(0.150, 0.300, 0.0, 0.0, false);
  }

  TEST_F(EMMDerivativesTest, NormalPhase_NoEMM_T120_muB0) {
    RunChecks(0.120, 0.0, 0.0, 0.0, false);
  }

  //============================================================================
  // Normal phase tests (with EMM on pions)
  //============================================================================

  TEST_F(EMMDerivativesTest, NormalPhase_EMM_T150_muB0) {
    RunChecks(0.150, 0.0, 0.0, 0.0, true);
  }

  TEST_F(EMMDerivativesTest, NormalPhase_EMM_T150_muB300) {
    RunChecks(0.150, 0.300, 0.0, 0.0, true);
  }

  TEST_F(EMMDerivativesTest, NormalPhase_EMM_T120_muB0) {
    RunChecks(0.120, 0.0, 0.0, 0.0, true);
  }

  TEST_F(EMMDerivativesTest, NormalPhase_EMM_T120_muB200) {
    RunChecks(0.120, 0.200, 0.0, 0.0, true);
  }

  //============================================================================
  // BEC phase tests (EMM + large muQ pushes pi+ into BEC)
  //============================================================================

  TEST_F(EMMDerivativesTest, BECPhase_EMM_T100_muQ200) {
    RunChecks(0.100, 0.0, 0.200, 0.0, true);
  }

  TEST_F(EMMDerivativesTest, BECPhase_EMM_T80_muQ200) {
    RunChecks(0.080, 0.0, 0.200, 0.0, true);
  }

  TEST_F(EMMDerivativesTest, NearBEC_EMM_T120_muQ150) {
    RunChecks(0.120, 0.0, 0.150, 0.0, true);
  }

  TEST_F(EMMDerivativesTest, BECPhase_EMM_T100_muB300_muQ200) {
    RunChecks(0.100, 0.300, 0.200, 0.0, true);
  }

}
