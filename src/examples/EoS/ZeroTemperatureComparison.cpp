/*
 * Thermal-FIST package
 *
 * Copyright (c) 2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

/**
 * \example ZeroTemperatureComparison.cpp
 *
 * Prints all key thermodynamic quantities side-by-side at
 *   T = 0,  1e-4,  1e-3,  1e-2  GeV
 * for a fixed muB, for each of four models:
 *   Ideal HRG, EV-Diagonal HRG, QvdW HRG, RealGas (CS+VDW) HRG.
 *
 * Purpose: quickly see which quantities are computed reliably at
 * very low T and which suffer from numerical issues.
 */
#include <cmath>
#include <cstdlib>
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

// -------------------------------------------------------------------
//  Collect all quantities at one (T, muB) point
// -------------------------------------------------------------------
struct ThermoRow {
  double T;
  double P, e, s, rhoB, rhoQ;
  double cs2_B, cT2_B;
  double cs2_BQ, cT2_BQ;
  double cMu;        // C_mu  = T ds/dT|_mu  (heat capacity at fixed mu)
  double dsdT_mu;     // ds/dT|_mu  (= C_mu/T for T>0, Sommerfeld for T=0)
  double chi2BB;      // baryon number susceptibility d^2P/dmuB^2 [GeV^-2]
  double chi3BB;      // baryon number susceptibility d^3P/dmuB^3 [GeV^-3]
  double chi4BB;      // baryon number susceptibility d^4P/dmuB^4 [GeV^-4]
};

static ThermoRow Evaluate(ThermalModelBase* model) {
  ThermoRow r;
  r.T = model->Parameters().T;

  model->CalculatePrimordialDensities();
  model->CalculateFluctuations();
  model->CalculateTemperatureDerivatives();

  r.P    = model->Pressure();
  r.e    = model->EnergyDensity();
  r.s    = model->EntropyDensity();
  r.rhoB = model->BaryonDensity();
  r.rhoQ = model->ElectricChargeDensity();

  auto safe = [](double v) { return std::isfinite(v) ? v : std::numeric_limits<double>::quiet_NaN(); };
  r.cs2_B  = safe(model->cs2(true, false, false, false));
  r.cT2_B  = safe(model->cT2(true, false, false, false));
  r.cs2_BQ = safe(model->cs2(true, true, false, false));
  r.cT2_BQ = safe(model->cT2(true, true, false, false));

  r.cMu = model->HeatCapacityMu();

  if (r.T > 0.)
    r.dsdT_mu = r.cMu / r.T;
  else {
    try   { r.dsdT_mu = model->CalculateEntropyDensityDerivativeTZeroTemperature(); }
    catch (...) { r.dsdT_mu = std::numeric_limits<double>::quiet_NaN(); }
  }

  r.chi2BB = model->SusceptibilityDimensionfull(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge)
             * xMath::GeVtoifm3();

  // chi3BB via CalculateChargeFluctuations with dimensionfull=true
  // Returns d^3P/dmuB^3 in [GeV^{-2}] units (natural units)
  {
    vector<double> chargesB(model->TPS()->ComponentsNumber());
    for (int i = 0; i < model->TPS()->ComponentsNumber(); ++i)
      chargesB[i] = model->TPS()->Particles()[i].BaryonCharge();

    auto chis = model->CalculateChargeFluctuations(chargesB, 4, true);
    r.chi3BB = chis[2] * xMath::GeVtoifm3();
    r.chi4BB = chis[3] * xMath::GeVtoifm3();
  }

  return r;
}


// -------------------------------------------------------------------
//  Print table for one model
// -------------------------------------------------------------------
static void PrintTable(const string& modelName,
                       ThermalModelBase* model,
                       double muB,
                       const vector<double>& Tvals)
{
  model->SetBaryonChemicalPotential(muB);
  model->SetElectricChemicalPotential(0.);
  model->SetStrangenessChemicalPotential(0.);
  model->SetCharmChemicalPotential(0.);

  // Collect rows
  vector<ThermoRow> rows;
  for (double T : Tvals) {
    model->SetTemperature(T);
    rows.push_back(Evaluate(model));
  }

  // Print
  const int w = 16;
  const int wlab = 22;
  cout << "\n";
  cout << "================================================================\n";
  cout << "  " << modelName << "    muB = " << muB << " GeV\n";
  cout << "================================================================\n";

  // Header
  cout << setw(wlab) << left << "Quantity";
  for (auto& r : rows)
    cout << setw(w) << right << (r.T == 0. ? "T=0" :
                                  r.T == 1.e-6 ? "T=0.001MeV" :
                                  r.T == 1.e-5 ? "T=0.01MeV" :
                                  r.T == 1.e-4 ? "T=0.1MeV" :
                                  r.T == 1.e-3 ? "T=1MeV" :
                                  r.T == 1.e-2 ? "T=10MeV" : "T=?");
  cout << "\n";
  cout << string(wlab + w * rows.size(), '-') << "\n";

  auto line = [&](const char* label, auto getter) {
    cout << setw(wlab) << left << label;
    for (auto& r : rows)
      cout << setw(w) << right << scientific << setprecision(6) << getter(r);
    cout << "\n";
  };

  line("P [GeV/fm3]",      [](const ThermoRow& r){ return r.P; });
  line("e [GeV/fm3]",      [](const ThermoRow& r){ return r.e; });
  line("s [fm-3]",         [](const ThermoRow& r){ return r.s; });
  line("rhoB [fm-3]",      [](const ThermoRow& r){ return r.rhoB; });
  line("rhoQ [fm-3]",      [](const ThermoRow& r){ return r.rhoQ; });
  line("chi2_BB [fm-3]",    [](const ThermoRow& r){ return r.chi2BB; });
  line("chi3_BB [fm-3GeV-1]",[](const ThermoRow& r){ return r.chi3BB; });
  line("chi4_BB [fm-3GeV-2]",[](const ThermoRow& r){ return r.chi4BB; });
  line("cs2(B)",           [](const ThermoRow& r){ return r.cs2_B; });
  line("cT2(B)",           [](const ThermoRow& r){ return r.cT2_B; });
  line("cs2(BQ)",          [](const ThermoRow& r){ return r.cs2_BQ; });
  line("cT2(BQ)",          [](const ThermoRow& r){ return r.cT2_BQ; });
  line("C_mu [fm-3]",      [](const ThermoRow& r){ return r.cMu; });
  line("ds/dT|mu [fm-3/GeV]", [](const ThermoRow& r){ return r.dsdT_mu; });

  // Relative deviations from T=0
  if (rows.size() > 1) {
    const auto& ref = rows[0];
    cout << "\n  Relative deviations from T=0:\n";
    cout << string(wlab + w * rows.size(), '-') << "\n";

    auto relline = [&](const char* label, auto getter) {
      double refv = getter(ref);
      cout << setw(wlab) << left << label;
      for (size_t i = 0; i < rows.size(); ++i) {
        if (i == 0) { cout << setw(w) << right << "---"; continue; }
        double v = getter(rows[i]);
        if (std::isfinite(v) && std::isfinite(refv) && std::abs(refv) > 1.e-30)
          cout << setw(w) << right << scientific << setprecision(3) << (v - refv) / refv;
        else
          cout << setw(w) << right << "n/a";
      }
      cout << "\n";
    };

    relline("dP/P",         [](const ThermoRow& r){ return r.P; });
    relline("de/e",         [](const ThermoRow& r){ return r.e; });
    relline("drhoB/rhoB",   [](const ThermoRow& r){ return r.rhoB; });
    relline("dchi3/chi3(B)", [](const ThermoRow& r){ return r.chi3BB; });
    relline("dchi4/chi4(B)", [](const ThermoRow& r){ return r.chi4BB; });
    relline("dcs2/cs2(B)",  [](const ThermoRow& r){ return r.cs2_B; });
    relline("dcT2/cT2(B)",  [](const ThermoRow& r){ return r.cT2_B; });
    relline("d(dsdT)/dsdT", [](const ThermoRow& r){ return r.dsdT_mu; });
  }

  cout << "\n";
}

// -------------------------------------------------------------------
//  main
// -------------------------------------------------------------------
int main(int argc, char *argv[])
{
  double muB = 1.2;
  if (argc > 1) muB = atof(argv[1]);

  cout << "Zero-temperature comparison at muB = " << muB << " GeV\n";
  cout << "Temperatures: 0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2 GeV\n";

  string inputDir = string(ThermalFIST_INPUT_FOLDER);
  ThermalParticleSystem tps(inputDir + "/list/PDG2025/list.dat");

  const vector<double> Tvals = {0., 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2};

  // --- 1) Ideal HRG ---
  {
    ThermalModelIdeal model(&tps);
    model.SetStatistics(true);
    model.SetUseWidth(false);
    model.SetCalculationType(IdealGasFunctions::Quadratures);
    PrintTable("Ideal HRG", &model, muB, Tvals);
  }

  // --- 2) EV-Diagonal HRG ---
  {
    ThermalModelEVDiagonal model(&tps);
    model.SetStatistics(true);
    model.SetUseWidth(false);
    model.SetCalculationType(IdealGasFunctions::Quadratures);
    double b = 1.0;
    for (int i = 0; i < tps.ComponentsNumber(); ++i)
      for (int j = 0; j < tps.ComponentsNumber(); ++j)
        if (tps.Particle(i).BaryonCharge() != 0 && tps.Particle(j).BaryonCharge() != 0)
          model.SetVirial(i, j, b);
    PrintTable("EV-Diagonal HRG (b=1 fm3)", &model, muB, Tvals);
  }

  // --- 3) QvdW HRG ---
  {
    ThermalModelVDW model(&tps);
    model.SetStatistics(true);
    model.SetUseWidth(false);
    model.SetCalculationType(IdealGasFunctions::Quadratures);
    SetVDWHRGInteractionParameters(&model, 0.329, 3.42);
    PrintTable("QvdW HRG (a=0.329, b=3.42)", &model, muB, Tvals);
  }

  // --- 4) RealGas HRG (CS + VDW MF) ---
  {
    ThermalModelRealGas model(&tps);
    model.SetStatistics(true);
    model.SetUseWidth(false);
    model.SetCalculationType(IdealGasFunctions::Quadratures);

    double a = 0.329, b = 3.42;
    {
      ExcludedVolumeModelBase *evmod = new ExcludedVolumeModelCS();
      auto bij = GetBaryonBaryonInteractionMatrix(model.TPS(), b);
      model.SetExcludedVolumeModel(new ExcludedVolumeModelCrosstermsGeneralized(evmod, bij));
    }
    {
      auto aij = GetBaryonBaryonInteractionMatrix(model.TPS(), a);
      model.SetMeanFieldModel(new MeanFieldModelMultiVDW(aij));
    }
    PrintTable("RealGas HRG (CS+VDW, a=0.329, b=3.42)", &model, muB, Tvals);
  }

  return 0;
}
