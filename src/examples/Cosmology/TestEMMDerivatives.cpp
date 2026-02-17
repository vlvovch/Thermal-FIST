/*
 * Thermal-FIST package
 *
 * Test: EMM temperature derivatives with differentiated EV/MF parameters
 *       r_meson = 0.15 fm, r_baryon = 0.4 fm, small mean-field attraction
 *       Tests thermodynamic identity (e+P = Ts + sum mu_i n_i)
 *       and temperature derivatives (ds/dT, de/dT) vs numerical FD
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "HRGBase.h"
#include "HRGRealGas.h"
#include "HRGEV.h"
#include "ThermalFISTConfig.h"
#include "HRGEV/ExcludedVolumeHelper.h"
#include "CosmicEos/EffectiveMassModel.h"

using namespace std;
using namespace thermalfist;

// Helper: set up a RealGas model with given EV/MF and optionally EMM on pions
void setupModel(ThermalModelRealGas& model, ThermalParticleSystem& tps,
                const vector<vector<double>>& bij, const vector<vector<double>>& aij,
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

int main() {
  // Load particle list
  string listpath = string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list.dat";
  ThermalParticleSystem TPS(listpath);
  int N = TPS.Particles().size();

  // EV radii
  double r_meson  = 0.15; // fm
  double r_baryon = 0.4;  // fm

  // Mean-field attraction parameter (GeV fm^3)
  double a_meson  = 0.05;
  double a_baryon = 0.3;

  // Build bij and aij matrices
  vector<vector<double>> bij(N, vector<double>(N, 0.));
  vector<vector<double>> aij(N, vector<double>(N, 0.));

  for (int i = 0; i < N; ++i) {
    int Bi = TPS.Particles()[i].BaryonCharge();
    double ri = (Bi == 0) ? r_meson : r_baryon;
    double ai_coeff = (Bi == 0) ? a_meson : a_baryon;
    for (int j = 0; j < N; ++j) {
      int Bj = TPS.Particles()[j].BaryonCharge();
      double rj = (Bj == 0) ? r_meson : r_baryon;
      double aj_coeff = (Bj == 0) ? a_meson : a_baryon;
      bij[i][j] = CuteHRGHelper::brr(ri, rj);
      aij[i][j] = sqrt(ai_coeff * aj_coeff);  // geometric mean
    }
  }

  // Helper lambda: run test at given T, muB, muQ, muS
  auto runTest = [&](double T, double muB, double muQ, double muS, bool useEMM, const string& label) {
    cout << "=== " << label << " ===" << endl;
    cout << "T = " << T << " GeV, muB = " << muB << " GeV, muQ = " << muQ << " GeV, muS = " << muS << " GeV" << endl;
    cout << "EMM on pions: " << (useEMM ? "YES" : "NO") << endl;

    ThermalParticleSystem TPScopy(TPS);
    ThermalModelRealGas model(&TPScopy);
    setupModel(model, TPScopy, bij, aij, useEMM, T, muB, muQ, muS);
    model.CalculateDensities();

    // Check BEC status of pions
    if (useEMM) {
      int pdgs[] = {211, 111, -211};
      const char* names[] = {"pi+", "pi0", "pi-"};
      for (int ip = 0; ip < 3; ++ip) {
        int idx = TPScopy.PdgToId(pdgs[ip]);
        auto* gd = TPScopy.Particles()[idx].GetGeneralizedDensity();
        if (gd) {
          cout << "  " << names[ip] << ": mu=" << model.ChemicalPotential(idx)
               << " m*=" << gd->EffectiveMass() << " BEC=" << gd->IsBECPhase() << endl;
        }
      }
    }

    // --- Thermodynamic identity: e + P = T*s + sum(mu_i * n_i) ---
    double P = model.CalculatePressure();
    double e = model.CalculateEnergyDensity();
    double s = model.CalculateEntropyDensity();

    double sum_mu_n = 0.;
    for (int i = 0; i < model.TPS()->Particles().size(); ++i) {
      double mu_i = model.ChemicalPotential(i);
      double n_i  = model.Densities()[i]; // primordial
      sum_mu_n += mu_i * n_i;
    }

    double lhs = e + P;
    double rhs = T * s + sum_mu_n;
    double rel_identity = (lhs != 0.) ? std::abs(lhs - rhs) / std::abs(lhs) : std::abs(lhs - rhs);
    cout << fixed << setprecision(8);
    cout << "  e + P           = " << lhs << endl;
    cout << "  T*s + sum(mu*n) = " << rhs << endl;
    cout << "  Rel. diff       = " << scientific << rel_identity << endl;

    // --- Temperature derivatives via CalculateTemperatureDerivatives ---
    model.CalculateTemperatureDerivatives();
    double dsdT_model = model.CalculateEntropyDensityDerivativeT();
    double dedT_model = model.CalculateEnergyDensityDerivativeT();

    // --- Numerical finite differences ---
    double dT = 1.e-4 * T;

    // s(T+dT/2)
    ThermalParticleSystem TPSp(TPS);
    ThermalModelRealGas model_p(&TPSp);
    setupModel(model_p, TPSp, bij, aij, useEMM, T + 0.5 * dT, muB, muQ, muS);
    model_p.CalculateDensities();
    double s_p = model_p.CalculateEntropyDensity();
    double e_p = model_p.CalculateEnergyDensity();

    // s(T-dT/2)
    ThermalParticleSystem TPSm(TPS);
    ThermalModelRealGas model_m(&TPSm);
    setupModel(model_m, TPSm, bij, aij, useEMM, T - 0.5 * dT, muB, muQ, muS);
    model_m.CalculateDensities();
    double s_m = model_m.CalculateEntropyDensity();
    double e_m = model_m.CalculateEnergyDensity();

    double dsdT_fd = (s_p - s_m) / dT;
    double dedT_fd = (e_p - e_m) / dT;

    double rel_dsdT = (dsdT_fd != 0.) ? std::abs(dsdT_model - dsdT_fd) / std::abs(dsdT_fd) : std::abs(dsdT_model - dsdT_fd);
    double rel_dedT = (dedT_fd != 0.) ? std::abs(dedT_model - dedT_fd) / std::abs(dedT_fd) : std::abs(dedT_model - dedT_fd);

    cout << fixed << setprecision(8);
    cout << "  ds/dT (model) = " << dsdT_model << "  (FD) = " << dsdT_fd
         << "  rel.diff = " << scientific << rel_dsdT << endl;
    cout << "  de/dT (model) = " << fixed << dedT_model << "  (FD) = " << dedT_fd
         << "  rel.diff = " << scientific << rel_dedT << endl;

    // --- Check dP/dmuB = nB and dP/dmuQ = nQ ---
    double nB_model = model.CalculateBaryonDensity();
    double nQ_model = model.CalculateChargeDensity();

    double dmu = 1.e-4 * T;  // small shift in chemical potential

    // dP/dmuB via FD
    {
      ThermalParticleSystem TPS_Bp(TPS), TPS_Bm(TPS);
      ThermalModelRealGas m_Bp(&TPS_Bp), m_Bm(&TPS_Bm);
      setupModel(m_Bp, TPS_Bp, bij, aij, useEMM, T, muB + 0.5 * dmu, muQ, muS);
      setupModel(m_Bm, TPS_Bm, bij, aij, useEMM, T, muB - 0.5 * dmu, muQ, muS);
      m_Bp.CalculateDensities();
      m_Bm.CalculateDensities();
      double dPdmuB_fd = (m_Bp.CalculatePressure() - m_Bm.CalculatePressure()) / dmu;
      double rel_nB = (std::abs(nB_model) > 1.e-20) ?
        std::abs(nB_model - dPdmuB_fd) / std::abs(nB_model) : std::abs(nB_model - dPdmuB_fd);
      cout << fixed << setprecision(8);
      cout << "  nB (model)    = " << nB_model << "  dP/dmuB (FD) = " << dPdmuB_fd
           << "  rel.diff = " << scientific << rel_nB << endl;
    }

    // dP/dmuQ via FD
    {
      ThermalParticleSystem TPS_Qp(TPS), TPS_Qm(TPS);
      ThermalModelRealGas m_Qp(&TPS_Qp), m_Qm(&TPS_Qm);
      setupModel(m_Qp, TPS_Qp, bij, aij, useEMM, T, muB, muQ + 0.5 * dmu, muS);
      setupModel(m_Qm, TPS_Qm, bij, aij, useEMM, T, muB, muQ - 0.5 * dmu, muS);
      m_Qp.CalculateDensities();
      m_Qm.CalculateDensities();
      double dPdmuQ_fd = (m_Qp.CalculatePressure() - m_Qm.CalculatePressure()) / dmu;
      double rel_nQ = (std::abs(nQ_model) > 1.e-20) ?
        std::abs(nQ_model - dPdmuQ_fd) / std::abs(nQ_model) : std::abs(nQ_model - dPdmuQ_fd);
      cout << fixed << setprecision(8);
      cout << "  nQ (model)    = " << nQ_model << "  dP/dmuQ (FD) = " << dPdmuQ_fd
           << "  rel.diff = " << scientific << rel_nQ << endl;
    }

    cout << endl;
  };

  cout << "====================================================" << endl;
  cout << "  EMM T-derivative test with differentiated EV/MF" << endl;
  cout << "  r_meson = " << r_meson << " fm, r_baryon = " << r_baryon << " fm" << endl;
  cout << "  a_meson = " << a_meson << ", a_baryon = " << a_baryon << " GeV fm^3" << endl;
  cout << "====================================================" << endl << endl;

  // Normal phase tests
  runTest(0.150, 0.0,   0.0,  0.0, false, "T=150 MeV, muB=0, no EMM");
  runTest(0.150, 0.0,   0.0,  0.0, true,  "T=150 MeV, muB=0, with EMM");
  runTest(0.150, 0.300, 0.0,  0.0, false, "T=150 MeV, muB=300 MeV, no EMM");
  runTest(0.150, 0.300, 0.0,  0.0, true,  "T=150 MeV, muB=300 MeV, with EMM");
  runTest(0.120, 0.0,   0.0,  0.0, false, "T=120 MeV, muB=0, no EMM");
  runTest(0.120, 0.0,   0.0,  0.0, true,  "T=120 MeV, muB=0, with EMM");
  runTest(0.120, 0.200, 0.0,  0.0, true,  "T=120 MeV, muB=200 MeV, with EMM");

  // BEC phase tests: large muQ pushes pi+ into BEC (mu_pi+ = muQ > m_pi)
  runTest(0.100, 0.0,   0.200, 0.0, true,  "T=100 MeV, muQ=200 MeV, with EMM (BEC)");
  runTest(0.080, 0.0,   0.200, 0.0, true,  "T=80 MeV, muQ=200 MeV, with EMM (BEC)");
  runTest(0.120, 0.0,   0.150, 0.0, true,  "T=120 MeV, muQ=150 MeV, with EMM (near BEC)");
  runTest(0.100, 0.300, 0.200, 0.0, true,  "T=100 MeV, muB=300, muQ=200 MeV, with EMM (BEC)");

  return 0;
}

/**
 * \example TestEMMDerivatives.cpp
 *
 * Tests the thermodynamic consistency of the RealGas HRG model
 * with effective mass model (EMM) for pions and differentiated
 * excluded-volume/mean-field parameters for mesons and baryons.
 *
 * The program verifies:
 *   1. The thermodynamic identity: e + P = T*s + sum(mu_i * n_i).
 *   2. Analytical temperature derivatives (ds/dT, de/dT) against numerical finite differences.
 *   3. The relations dP/dmuB = nB and dP/dmuQ = nQ.
 *
 * Tests are run at various temperatures and chemical potentials,
 * including configurations in the Bose-Einstein condensation (BEC) phase.
 */
