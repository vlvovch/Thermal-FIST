/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdio>
#include <vector>
#include <string>

#include "HRGBase.h"
#include "HRGEV.h"
#include "HRGVDW.h"
#include "HRGFit.h"

#include "ThermalFISTConfig.h"

#include "CosmicEos/CosmicEoS.h"


using namespace std;
using namespace thermalfist;


// Compute chemical potentials as a function of temperature in the Early Universe
// Reference: V. Vovchenko et al., Phys. Rev. Lett. 126 (2021) 012701, https://inspirehep.net/literature/1815329
int main(int argc, char *argv[])
{
  // Baryon asymmetry
  double b = 8.6e-11;
  // Electric charge asymmetry
  double q = 0.0;

  // Electron flavor asymmetry (read from command line)
  double le = 0.0;
  if (argc > 1)
    le = atof(argv[1]);

  // Muon flavor asymmetry (read from command line)
  double lmu = 0.0;
  if (argc > 2)
    lmu = atof(argv[2]);

  // Tau flavor asymmetry (read from command line)
  double ltau = 0.0;
  if (argc > 3)
    ltau = atof(argv[3]);

  // If very small, set to zero
  if (fabs(le) < 1.e-13)
    le = 0.0;
  if (fabs(lmu) < 1.e-13)
    lmu = 0.0;
  if (fabs(ltau) < 1.e-13)
    ltau = 0.0;


  // Particle list for the HRG model
  ThermalParticleSystem parts(string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list.dat");

  // HRG model
  ThermalModelIdeal modelHRG(&parts);
  
  // For interactin HRG use e.g.
  // ThermalModelVDW model(&parts);
  // double a = 0.329, b = 3.42;
  // SetVDWHRGInteractionParameters(&model, a, b);

  // Resonance widths
  bool useWidth  = false;
  modelHRG.SetUseWidth(useWidth);
  if (useWidth)
    modelHRG.SetUseWidth(ThermalParticle::eBWconstBR);

  // Quantum statistics
  bool useQStats = true;
  modelHRG.SetStatistics(useQStats);

  // Effective mass model for pions (needed for BEC)
  bool interactingpions = true;

  // Cosmic EoS object
  CosmicEoS cosmos(&modelHRG, interactingpions);




  cout << "le + lmu = " << setw(15) << le + lmu << endl;
  cout << "le - lmu = " << setw(15) << le - lmu << endl;
  cout << "le   = " << setw(15) << le << endl;
  cout << "lmu  = " << setw(15) << lmu << endl;
  cout << "ltau = " << setw(15) << ltau << endl;

  cosmos.SetAsymmetries(vector<double>({ b, 0., le, lmu, ltau }));

  // Range of temperatures
  double Tmin = 0.010; // GeV
  double Tmax = 0.180; // GeV
  double dT = 0.001;   // GeV
  vector<double> Temps;
  for (double tT = Tmin; tT <= Tmax + 0.1 * dT; tT += dT) {
    Temps.push_back(tT);
  }

  // Output file (write le + lmu since it's the main ingredient for pion condersation)
  string filename = "CosmicTrajectory";
  {
    char cc[500];
    sprintf(cc, "le+lmu.%lf", le + lmu);
    filename += "." + string(cc);
  }
  filename += ".dat";
  
  ofstream fout(filename);

  // We will output the trajectories and the various properties of the EoS
  fout << setw(15) << "T[MeV]" << " ";
  fout << setw(15) << "muB[MeV]" << " ";
  fout << setw(15) << "muQ[MeV]" << " ";
  fout << setw(15) << "mue[MeV]" << " ";
  fout << setw(15) << "mum[MeV]" << " ";
  fout << setw(15) << "mut[MeV]" << " ";
  fout << setw(15) << "pion_bec" << " ";
  fout << setw(15) << "nI/T3" << " ";
  fout << setw(15) << "pT4" << " ";
  fout << setw(15) << "eT4" << " ";
  fout << setw(15) << "IT4" << " ";
  fout << setw(15) << "pT4_QCD" << " ";
  fout << setw(15) << "eT4_QCD" << " ";
  fout << setw(15) << "IT4_QCD" << " ";
  fout << setw(15) << "sT3" << " ";
  fout << setw(15) << "sT3_QCD" << " ";
  fout << setw(15) << "IT4_e" << " ";
  fout << setw(15) << "IT4_mu" << " ";
  fout << setw(15) << "IT4_tau" << " ";
  fout << setw(15) << "nQ_QCD/T3" << " ";
  fout << setw(15) << "npi_QCD/T3" << " ";
  fout << setw(15) << "ne/T3" << " ";
  fout << setw(15) << "nmu/T3" << " ";
  fout << setw(15) << "ntau/T3" << " ";
  fout << endl;

  // On-screen
  cout << setw(15) << "T[MeV]" << " ";
  cout << setw(15) << "muB[MeV]" << " ";
  cout << setw(15) << "muQ[MeV]" << " ";
  cout << setw(15) << "mue[MeV]" << " ";
  cout << setw(15) << "mum[MeV]" << " ";
  cout << setw(15) << "mut[MeV]" << " ";
  cout << setw(15) << "pion_bec" << endl;

  // Timer
  double wt1 = get_wall_time();
  int iters = 0; // number of iterations


  // Initial guess for the chemical potentials (muB, muQ, mue, mumu, mutau)
  vector<double> prev = vector<double>({ 0.700, -1.e-7, -1.e-7, -1.e-7, -1.e-7 });

  // Loop over the temperatures
  for (auto&& T : Temps) {
    // Solve the equations for the chemical potentials
    vector<double> chems = cosmos.SolveChemicalPotentials(T, prev);
    iters++;
    // Will be the initial guess for the next temperature
    prev = chems;

    // Check for pion condensation
    if (!interactingpions && abs(chems[1]) > 0.139)
      break;

    cout << setw(15) << T * 1.e3 << " ";
    cout << setw(15) << chems[0] * 1.e3 << " ";
    cout << setw(15) << chems[1] * 1.e3 << " ";
    cout << setw(15) << chems[2] * 1.e3 << " ";
    cout << setw(15) << chems[3] * 1.e3 << " ";
    cout << setw(15) << chems[4] * 1.e3 << " ";
    cout << setw(15) << cosmos.InPionCondensedPhase() << endl;


    // Write to file
    fout << setw(15) << T * 1.e3 << " ";
    fout << setw(15) << chems[0] * 1.e3 << " ";
    fout << setw(15) << chems[1] * 1.e3 << " ";
    fout << setw(15) << chems[2] * 1.e3 << " ";
    fout << setw(15) << chems[3] * 1.e3 << " ";
    fout << setw(15) << chems[4] * 1.e3 << " ";
    fout << setw(15) << cosmos.InPionCondensedPhase() << " ";

    fout << setw(15) << cosmos.IsospinChargeDensity() / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";

    fout << setw(15) << cosmos.Pressure() / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << cosmos.EnergyDensity() / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << (cosmos.EnergyDensity() - 3. * cosmos.Pressure()) / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << cosmos.PressureHRG() / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << cosmos.EnergyDensityHRG() / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << (cosmos.EnergyDensityHRG() - 3. * cosmos.PressureHRG()) / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << cosmos.EntropyDensity() / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << cosmos.EntropyDensityHRG() / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";

    //fout << setw(15) << (cosmos.EnergyDensity() - 3. * cosmos.Pressure()) / pow(T, 4) / pow(xMath::GeVtoifm(), 3);
    fout << setw(15) << (cosmos.EnergyDensityChargedLepton(0) - 3. * cosmos.PressureChargedLepton(0)) / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << (cosmos.EnergyDensityChargedLepton(1) - 3. * cosmos.PressureChargedLepton(1)) / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << (cosmos.EnergyDensityChargedLepton(2) - 3. * cosmos.PressureChargedLepton(2)) / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";

    fout << setw(15) << cosmos.ElectricChargeDensityHRG() / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << (cosmos.ElectricChargeDensityHRG() - cosmos.HRGModel()->ElectricChargeDensity()) / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << cosmos.NetDensityChargedLepton(0) / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << cosmos.NetDensityChargedLepton(1) / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << setw(15) << cosmos.NetDensityChargedLepton(2) / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";
    fout << endl;
  }

  // Optional cross-check (run backwards in temperature and make sure the trajectory is the same)
  if (0) {
    for(int iT = Temps.size() - 1; 0 && iT >= 0; iT--) {
      double T = Temps[iT];
      vector<double> chems = cosmos.SolveChemicalPotentials(T, prev);
      iters++;
      prev = chems;

      if (!interactingpions && abs(chems[1]) > 0.139)
        break;

      cout << setw(15) << T * 1.e3 << " ";
      cout << setw(15) << chems[0] * 1.e3 << " ";
      cout << setw(15) << chems[1] * 1.e3 << " ";
      cout << setw(15) << chems[2] * 1.e3 << " ";
      cout << setw(15) << chems[3] * 1.e3 << " ";
      cout << setw(15) << chems[4] * 1.e3 << " ";
      cout << setw(15) << cosmos.InPionCondensedPhase() << endl;

      fout << setw(15) << T * 1.e3 << " ";
      fout << setw(15) << chems[0] * 1.e3 << " ";
      fout << setw(15) << chems[1] * 1.e3 << " ";
      fout << setw(15) << chems[2] * 1.e3 << " ";
      fout << setw(15) << chems[3] * 1.e3 << " ";
      fout << setw(15) << chems[4] * 1.e3 << " ";
      fout << setw(15) << cosmos.InPionCondensedPhase() << " ";
      fout << setw(15) << cosmos.Pressure() / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
      fout << setw(15) << cosmos.EnergyDensity() / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
      fout << setw(15) << (cosmos.EnergyDensity() - 3. * cosmos.Pressure()) / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
      fout << setw(15) << cosmos.PressureHRG() / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
      fout << setw(15) << cosmos.EnergyDensityHRG() / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
      fout << setw(15) << (cosmos.EnergyDensityHRG() - 3. * cosmos.PressureHRG()) / pow(T, 4) / pow(xMath::GeVtoifm(), 3) << " ";
      fout << setw(15) << cosmos.EntropyDensity() / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";
      fout << setw(15) << cosmos.EntropyDensityHRG() / pow(T, 3) / pow(xMath::GeVtoifm(), 3) << " ";
      fout << endl;
    }
  }

  double wt2 = get_wall_time();

  printf("%30s %lf s\n", "Running time:", (wt2 - wt1));
  printf("%30s %lf s\n", "Time per single calculation:", (wt2 - wt1) / iters);

  return 0;
}

/**
 * \example CosmicTrajectory.cpp
 * 
 * An example of calculating the chemical potentials as a function of temperature
 * in the Early Universe.
 * 
 * The program does the following:
 * 
 *   1. Reads the electron, muon, and tau flavor asymmetries from the command line.
 *   2. Initializes the HRG model and the Cosmic EoS object.
 *   3. Sets the asymmetries and the range of temperatures.
 *   4. Solves the equations for the chemical potentials at each temperature.
 *   5. Outputs the results to a file and the console.
 * 
 * Usage:
 * ~~~.bash
 * CosmicTrajectory <le> <lmu> <ltau>
 * ~~~
 * 
 * Where <le>, <lmu>, and <ltau> are the electron, muon, and tau flavor asymmetries, respectively.
 * 
 * Reference: V. Vovchenko et al., Phys. Rev. Lett. 126 (2021) 012701, https://inspirehep.net/literature/1815329
 */
