/*
 * Thermal-FIST package
 *
 * Copyright (c) 2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "HRGBase.h"
#include "HRGEV.h"
#include "HRGRealGas.h"
#include "HRGVDW.h"
#include "HRGEV/ExcludedVolumeHelper.h"
#include "ThermalFISTConfig.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

// Computes zero-temperature EoS in the QvdW setup for either:
//   - ThermalModelVDW (useRG=0)
//   - ThermalModelRealGas with CS EV + VDW mean field (useRG=1)
//
// Usage:
//   ZeroTemperatureQvdW <useRG> <a> <b> <muBmin> <muBmax> <dmuB> <outputfile>
//
// Defaults:
//   useRG  = 0      (0: ThermalModelVDW, 1: ThermalModelRealGas)
//   a      = 0.329  GeV*fm^3
//   b      = 3.42   fm^3
//   muBmin = 0.90   GeV
//   muBmax = 1.70   GeV
//   dmuB   = 0.01   GeV
int main(int argc, char *argv[])
{
  bool useRG = false;
  if (argc > 1) useRG = (atoi(argv[1]) != 0);

  double a = 0.329;
  if (argc > 2) a = atof(argv[2]);

  double b = 3.42;
  if (argc > 3) b = atof(argv[3]);

  double muBmin = 0.90;
  if (argc > 4) muBmin = atof(argv[4]);

  double muBmax = 1.70;
  if (argc > 5) muBmax = atof(argv[5]);

  double dmuB = 0.01;
  if (argc > 6) dmuB = atof(argv[6]);

  string outputfile = string("ZeroT-EoS-") + (useRG ? "RG-HRG" : "QvdW-HRG") + ".dat";
  if (argc > 7) outputfile = argv[7];

  if (dmuB <= 0.0) {
    cerr << "Error: dmuB must be positive." << endl;
    return 1;
  }
  if (muBmax < muBmin) {
    cerr << "Error: muBmax must be >= muBmin." << endl;
    return 1;
  }

  ThermalParticleSystem tps(string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2025/list.dat");

  unique_ptr<ThermalModelBase> model;
  if (!useRG) {
    auto *vdw = new ThermalModelVDW(&tps);
    SetVDWHRGInteractionParameters(vdw, a, b);
    model.reset(vdw);
    cout << "Model: ThermalModelVDW (QvdW-HRG)" << endl;
  } else {
    auto *rg = new ThermalModelRealGas(&tps);

    if (b != 0.0) {
      ExcludedVolumeModelBase *evmod = new ExcludedVolumeModelCS();
      auto vdWb = GetBaryonBaryonInteractionMatrix(rg->TPS(), b);
      rg->SetExcludedVolumeModel(new ExcludedVolumeModelCrosstermsGeneralized(evmod, vdWb));
    }
    if (a != 0.0) {
      auto vdWa = GetBaryonBaryonInteractionMatrix(rg->TPS(), a);
      rg->SetMeanFieldModel(new MeanFieldModelMultiVDW(vdWa));
    }

    model.reset(rg);
    cout << "Model: ThermalModelRealGas (CS + VDW mean field)" << endl;
  }

  model->SetUseWidth(false);
  model->SetStatistics(true);
  model->SetCalculationType(IdealGasFunctions::Quadratures);

  // Zero-temperature setup with fixed charge chemical potentials.
  model->SetTemperature(0.0);
  model->SetElectricChemicalPotential(0.0);
  model->SetStrangenessChemicalPotential(0.0);
  model->SetCharmChemicalPotential(0.0);

  ofstream fout(outputfile.c_str());
  if (!fout.is_open()) {
    cerr << "Error: cannot open output file " << outputfile << endl;
    return 1;
  }

  const int tab = 20;
  fout << setw(tab) << "muB[GeV]" << " "
       << setw(tab) << "P[GeV/fm3]" << " "
       << setw(tab) << "e[GeV/fm3]" << " "
       << setw(tab) << "s[fm-3]" << " "
       << setw(tab) << "rhoB[fm-3]" << " "
       << setw(tab) << "rhoQ[fm-3]" << " "
       << setw(tab) << "rhoS[fm-3]" << " "
       << setw(tab) << "(1/3-p/e)" << " "
       << setw(tab) << "cT2(B)" << " "
       << setw(tab) << "cs2(B)" << " "
       << setw(tab) << "cT2(BQ)" << " "
       << setw(tab) << "cs2(BQ)" << " "
       << setw(tab) << "cT2(BS)" << " "
       << setw(tab) << "cs2(BS)" << " "
       << setw(tab) << "cT2(BQS)" << " "
       << setw(tab) << "cs2(BQS)" << " "
       << setw(tab) << "cV(B)[fm-3]" << " "
       << setw(tab) << "cV(BQ)[fm-3]" << " "
       << setw(tab) << "cV(BS)[fm-3]" << " "
       << setw(tab) << "cV(BQS)[fm-3]" << " "
       << setw(tab) << "cMu[fm-3]" << " "
       << setw(tab) << "(ds/dT)_mu_num@1e-4" << " "
       << setw(tab) << "(ds/dT)_mu_num@1e-3" << " "
       << setw(tab) << "(ds/dT)_mu_num@1e-2" << " "
       << setw(tab) << "(ds/dT)_mu_an@1e-4" << " "
       << setw(tab) << "(ds/dT)_mu_an@1e-3" << " "
       << setw(tab) << "(ds/dT)_mu_an@1e-2" << " "
       << setw(tab) << "(ds/dT)_mu_an" << "\n";

  const double wt1 = get_wall_time();
  int iters = 0;

  for (double muB = muBmin; muB <= muBmax + 0.1 * dmuB; muB += dmuB) {
    model->SetBaryonChemicalPotential(muB);
    model->CalculatePrimordialDensities();
    model->CalculateFluctuations();
    model->CalculateTemperatureDerivatives();

    const double p = model->Pressure();
    const double e = model->EnergyDensity();
    const double s = model->EntropyDensity();
    const double rhoB = model->BaryonDensity();
    const double rhoQ = model->ElectricChargeDensity();
    const double rhoS = model->StrangenessDensity();
    const double trace = (e != 0.0) ? (1.0 / 3.0 - p / e) : 0.0;
    const double eps = 1.e-12;
    const bool hasB = (std::abs(rhoB) > eps);
    const bool hasQ = (std::abs(rhoQ) > eps);
    const bool hasS = (std::abs(rhoS) > eps);

    const auto safe_cT2 = [&](bool b, bool q, bool ss, double fallback = 0.0) {
      const double v = model->cT2(b, q, ss, false);
      return std::isfinite(v) ? v : fallback;
    };
    const auto safe_cs2 = [&](bool b, bool q, bool ss, double fallback = 0.0) {
      const double v = model->cs2(b, q, ss, false);
      return std::isfinite(v) ? v : fallback;
    };

    double cT2_B = 0.0, cs2_B = 0.0;
    double cT2_BQ = 0.0, cs2_BQ = 0.0;
    double cT2_BS = 0.0, cs2_BS = 0.0;
    double cT2_BQS = 0.0, cs2_BQS = 0.0;

    if (hasB) {
      cT2_B = safe_cT2(true, false, false, 0.0);
      cs2_B = safe_cs2(true, false, false, 0.0);

      // At T=0, if extra conserved charges vanish, constraints reduce to B-only.
      cT2_BQ = hasQ ? safe_cT2(true, true, false, cT2_B) : cT2_B;
      cs2_BQ = hasQ ? safe_cs2(true, true, false, cs2_B) : cs2_B;
      cT2_BS = hasS ? safe_cT2(true, false, true, cT2_B) : cT2_B;
      cs2_BS = hasS ? safe_cs2(true, false, true, cs2_B) : cs2_B;
      cT2_BQS = (hasQ && hasS) ? safe_cT2(true, true, true, cT2_B) : (hasS ? cT2_BS : cT2_BQ);
      cs2_BQS = (hasQ && hasS) ? safe_cs2(true, true, true, cs2_B) : (hasS ? cs2_BS : cs2_BQ);
    }

    // At exactly T=0 these heat capacities vanish by definition.
    const double cV_B = 0.0;
    const double cV_BQ = 0.0;
    const double cV_BS = 0.0;
    const double cV_BQS = 0.0;

    const double cMu = model->HeatCapacityMu();
    double dsdT_mu_an = model->CalculateEntropyDensityDerivativeT();

    // Numerical (ds/dT)_mu cross-check from one-sided derivatives at several dT.
    // For a degenerate Fermi system, corrections are ~T^2 (Sommerfeld), so fit vs dT^2.
    const std::vector<double> dTscan = {1.e-4, 1.e-3, 1.e-2};
    std::vector<double> dsdT_scan, dsdT_scan_an;
    dsdT_scan.reserve(dTscan.size());
    for (double dTnum : dTscan) {
      model->SetTemperature(dTnum);
      model->CalculatePrimordialDensities();
      const double sT = model->EntropyDensity();
      dsdT_scan.push_back((sT - s) / dTnum);

      double cVmu = model->HeatCapacityMu();
      dsdT_scan_an.push_back(cVmu / dTnum);
    }

    const double dsdT_mu_num_1e4 = dsdT_scan[0];
    const double dsdT_mu_num_1e3 = dsdT_scan[1];
    const double dsdT_mu_num_1e2 = dsdT_scan[2];

    const double dsdT_mu_an_1e4 = dsdT_scan_an[0];
    const double dsdT_mu_an_1e3 = dsdT_scan_an[1];
    const double dsdT_mu_an_1e2 = dsdT_scan_an[2];
    


    // Restore the baseline zero-temperature state.
    model->SetTemperature(0.0);
    model->CalculatePrimordialDensities();
    model->CalculateFluctuations();
    model->CalculateTemperatureDerivatives();

    fout << setw(tab) << muB << " "
         << setw(tab) << p << " "
         << setw(tab) << e << " "
         << setw(tab) << s << " "
         << setw(tab) << rhoB << " "
         << setw(tab) << rhoQ << " "
         << setw(tab) << rhoS << " "
         << setw(tab) << trace << " "
         << setw(tab) << cT2_B << " "
         << setw(tab) << cs2_B << " "
         << setw(tab) << cT2_BQ << " "
         << setw(tab) << cs2_BQ << " "
         << setw(tab) << cT2_BS << " "
         << setw(tab) << cs2_BS << " "
         << setw(tab) << cT2_BQS << " "
         << setw(tab) << cs2_BQS << " "
         << setw(tab) << cV_B << " "
         << setw(tab) << cV_BQ << " "
         << setw(tab) << cV_BS << " "
         << setw(tab) << cV_BQS << " "
         << setw(tab) << cMu << " "
         << setw(tab) << dsdT_mu_num_1e4 << " "
         << setw(tab) << dsdT_mu_num_1e3 << " "
         << setw(tab) << dsdT_mu_num_1e2 << " "
         << setw(tab) << dsdT_mu_an_1e4 << " "
         << setw(tab) << dsdT_mu_an_1e3 << " "
         << setw(tab) << dsdT_mu_an_1e2 << " "
         << setw(tab) << dsdT_mu_an << "\n";

    ++iters;
  }

  fout.close();

  const double wt2 = get_wall_time();
  cout << "Output file: " << outputfile << endl;
  cout << setw(30) << "Running time: " << (wt2 - wt1) << " s" << endl;
  if (iters > 0)
    cout << setw(30) << "Time per point: " << (wt2 - wt1) / iters << " s" << endl;

  return 0;
}

/**
 * \example ZeroTemperatureQvdW.cpp
 *
 * Zero-temperature equation-of-state scan for QvdW-HRG.
 * Demonstrates both ThermalModelVDW and ThermalModelRealGas setups.
 */
