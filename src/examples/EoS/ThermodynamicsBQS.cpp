/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cstdio>
#include <cassert>

#include "HRGBase.h"
#include "HRGEV.h"
#include "HRGFit.h"
#include "HRGVDW.h"
#include "HRGRealGas.h"

#include "HRGEV/ExcludedVolumeHelper.h"

#include "ThermalFISTConfig.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

// The range in T, muB, muQ, muS to scan
map<string, vector<double>> param_range = {
  { "T", { 0.100, 0.180, 0.001 } },  // GeV
  { "muB", { 0.000, 0.000, 0.01 } }, // GeV
  { "muQ", { 0.000, 0.000, 0.01 } }, // GeV
  { "muS", { 0.000, 0.000, 0.01 } }  // GeV
};

// Read the parameter range from file
void ReadParameterRangeFromFile(const std::string &filename)
{
  if (filename == "")
    return;
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    std::cerr << "Error: cannot open file " << filename << " using default parameter range" << std::endl;
    return;
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') { continue; } // skip empty lines and comments
    std::istringstream iss(line);
    std::string param;
    if (!(iss >> param)) continue;
    if (param[0] == '#') { continue; } // skip comments
    double min, max, step;
    if (!(iss >> min >> max >> step)) { break; } // error
    param_range[param] = { min, max, step };
  }

  fin.close();
}

// Calculates the HRG model equation of state properties
// for a given range of T, muB, muQ, muS
// Usage: ThermodynamicsBQS <param_range_file> <outputfile> <a> <b>
// * <a> -- parameter a for the QvdW model (GeV fm^3)
// * <b> -- parameter b for the QvdW model (fm^3)
// * useCS -- Use real gas model and Carnahan-Starling EV
// * <param_range_file> -- file with the parameter range
// * <outputfile> -- file to write the results to
int main(int argc, char *argv[])
{
  // van der Waals attraction
  double a = 0.;
  if (argc > 1)
    a = atof(argv[1]);
  
  // van der Waals repulsion
  double b = 0.;
  if (argc > 2)
    b = atof(argv[2]);

  // Use real gas
  bool useRG = false;
  if (argc > 3)
    useRG = atoi(argv[3]);
  
  // Parameter range from file
  string param_range_file = "";
  if (argc > 4)
    param_range_file = argv[4];
  ReadParameterRangeFromFile(param_range_file);

  // Model type
  // 0 - Ideal HRG, 1 - EV HRG (no attraction), 2 - QvdW HRG (from 1609.03975), 3 - Real gas HRG
  int ModelType = 0;
  string ModelPrefix = "Id-HRG";
  if (a == 0. && b == 0.) {
    ModelType = 0;
    ModelPrefix = "Id-HRG";
  }
  else if (a == 0.) {
    ModelType = 1;
    ModelPrefix = "EV-HRG";
  }
  else {
    ModelType = 2;
    ModelPrefix = "QvdW-HRG";
  }

  if (useRG) {
    ModelType = 3;
    ModelPrefix = "RG-HRG";
  }
  
  std::string outputfile = "Thermodynamics-" + ModelPrefix + "-output.dat";
  if (argc > 5)
    outputfile = argv[5];
  
  
  // Create the hadron list instance and read the list from file
  ThermalParticleSystem TPS(string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list.dat");  // <-- Default list, no light nuclei
  //ThermalParticleSystem TPS(string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list-withnuclei.dat");  // <-- Default list, with light nuclei

  // Create the ThermalModel instance
  // Choose the class which fits the required variant of HRG model
  ThermalModelBase *model;
  assert(ModelType >= 0 && ModelType <= 3);
  if (ModelType == 0) {
    model = new ThermalModelIdeal(&TPS);
    cout << "Performing calculation in the Id-HRG model..." << endl;
  }
  else if (ModelType == 1) {
    model = new ThermalModelEVCrossterms(&TPS);
    SetEVHRGInteractionParameters(model, b);
    cout << "Performing calculation in the EV-HRG model..." << endl;
  }
  else if (ModelType == 2) {
    model = new ThermalModelVDW(&TPS);
    SetVDWHRGInteractionParameters(model, a, b);
    cout << "Performing calculation in the QvdW-HRG model..." << endl;
  }
  else if (ModelType == 3) {
    model = new ThermalModelRealGas(&TPS);

    // Excluded volume model
    if (b != 0.) {
      // ExcludedVolumeModelBase *evmod = new ExcludedVolumeModelVDW();
      ExcludedVolumeModelBase *evmod = new ExcludedVolumeModelCS();
      auto vdWb = GetBaryonBaryonInteractionMatrix(model->TPS(), b);
      static_cast<ThermalModelRealGas *>(model)->SetExcludedVolumeModel(
              new ExcludedVolumeModelCrosstermsGeneralized(evmod, vdWb));
    }

    // Mean field model
    if (a != 0.) {
      auto vdWa = GetBaryonBaryonInteractionMatrix(model->TPS(), a);
      static_cast<ThermalModelRealGas *>(model)->SetMeanFieldModel(new MeanFieldModelMultiVDW(vdWa));
    }
    cout << "Performing calculation in the RG-HRG model..." << endl;
  }

  // Use (or not) finite resonance width
  bool useWidth = false;
  model->SetUseWidth(useWidth);
  if (useWidth)
    model->SetUseWidth(ThermalParticle::eBWconstBR);

  // Include (or not) quantum statistics
  bool useQStats = true;
  model->SetStatistics(useQStats);
  if (useQStats)
    model->SetCalculationType(IdealGasFunctions::Quadratures);

  // Prepare the output
  ofstream fout(outputfile);
  const int tabsize = 20;
  fout << std::setw(tabsize) << "T[GeV]" << " ";
  fout << std::setw(tabsize) << "muB[GeV]" << " ";
  fout << std::setw(tabsize) << "muQ[GeV]" << " ";
  fout << std::setw(tabsize) << "muS[GeV]" << " ";
  fout << std::setw(tabsize) << "P[GeV/fm3]" << " ";
  fout << std::setw(tabsize) << "e[GeV/fm3]" << " ";
  fout << std::setw(tabsize) << "s[fm-3]" << " ";
  fout << std::setw(tabsize) << "(1/3-p/e)" << " ";    // Another defition of the trace anomaly
  fout << std::setw(tabsize) << "rhoB[fm-3]" << " ";
  fout << std::setw(tabsize) << "rhoQ[fm-3]" << " ";
  fout << std::setw(tabsize) << "rhoS[fm-3]" << " ";
  // Temperature normalized quantities (a la lattice QCD)
  fout << std::setw(tabsize) << "P/T^4" << " ";
  fout << std::setw(tabsize) << "e/T^4" << " ";
  fout << std::setw(tabsize) << "s/T^3" << " ";
  fout << std::setw(tabsize) << "(e-3p)/T^3" << " ";
  fout << std::setw(tabsize) << std::endl;

  // Timer
  double wt1 = get_wall_time();
  int iters = 0; // number of iterations (phase diagram points)

  // Loop over the thermodynamic variables
  double Tmin = param_range["T"][0];
  double Tmax = param_range["T"][1];
  double dT   = param_range["T"][2];
  for(double T = Tmin; T <= Tmax + 0.1 * dT; T += dT) {
    model->SetTemperature(T);

    // Output the outer loop
    cout << "Calculating the temperature " << T << " GeV" << endl;

    double muBmin = param_range["muB"][0];
    double muBmax = param_range["muB"][1];
    double dmuB   = param_range["muB"][2];
    for(double muB = muBmin; muB <= muBmax + 0.1 * dmuB; muB += dmuB) {
      model->SetBaryonChemicalPotential(muB);

      double muQmin = param_range["muQ"][0];
      double muQmax = param_range["muQ"][1];
      double dmuQ   = param_range["muQ"][2];
      for(double muQ = muQmin; muQ <= muQmax + 0.1 * dmuQ; muQ += dmuQ) {
        model->SetElectricChemicalPotential(muQ);

        double muSmin = param_range["muS"][0];
        double muSmax = param_range["muS"][1];
        double dmuS   = param_range["muS"][2];
        for(double muS = muSmin; muS <= muSmax + 0.1 * dmuS; muS += dmuS) {
          model->SetStrangenessChemicalPotential(muS);

          // Perform the calculation
          model->CalculatePrimordialDensities();
          iters++;
          
          // Collect the output
          double p = model->Pressure();       // Pressure in GeV/fm3
          double e = model->EnergyDensity();  // Energy density in GeV/fm3
          double s = model->EntropyDensity(); // Entropy density in fm-3
          double rhoB = model->BaryonDensity();         // Baryon density in fm-3
          double rhoQ = model->ElectricChargeDensity(); // Electric charge density in fm-3
          double rhoS = model->StrangenessDensity();    // Strangeness density in fm-3
          double trace_anomaly = (1./3. - p/e);

          double pT4 = p / pow(T, 4) / xMath::GeVtoifm3();
          double eT4 = e / pow(T, 4) / xMath::GeVtoifm3();
          double sT3 = s / pow(T, 3) / xMath::GeVtoifm3();
          double IT4 = (e - 3 * p) / pow(T, 4) / xMath::GeVtoifm3();

          // Print to file
          fout << setw(tabsize) << T << " ";
          fout << setw(tabsize) << muB << " ";
          fout << setw(tabsize) << muQ << " ";
          fout << setw(tabsize) << muS << " ";
          fout << setw(tabsize) << p << " ";
          fout << setw(tabsize) << e << " ";
          fout << setw(tabsize) << s << " ";
          fout << setw(tabsize) << trace_anomaly << " ";
          fout << setw(tabsize) << rhoB << " ";
          fout << setw(tabsize) << rhoQ << " ";
          fout << setw(tabsize) << rhoS << " ";
          fout << setw(tabsize) << pT4 << " ";
          fout << setw(tabsize) << eT4 << " ";
          fout << setw(tabsize) << sT3 << " ";
          fout << setw(tabsize) << IT4 << " ";
          fout << setw(tabsize) << endl;
        }
      }
    }
  }
  
  fout.close();

  double wt2 = get_wall_time();
  cout << setw(30) <<  "Running time: " << (wt2 - wt1) << " s" << endl;
  cout << setw(30) <<  "Time per single calculation: " << (wt2 - wt1) / iters << " s" << endl;

  // Cleanup
  delete model;

  return 0;
}

/**
 * \example ThermodynamicsBQS.cpp
 * 
 * An example of calculating various thermodynamic quantities over a range
 * of T, \f$\mu_B\f$, \f$\mu_Q\f$, and \f$\mu_S\f$.
 * 
 * For each specified value of T, \f$\mu_B\f$, \f$\mu_Q\f$, and \f$\mu_S\f$,
 * the following quantities are calculated:
 *   - Pressure
 *   - Energy density
 *   - Entropy density
 *   - Trace anomaly
 *   - Baryon density
 *   - Electric charge density
 *   - Strangeness density
 *   - Temperature normalized quantities (a la lattice QCD)
 * 
 * Usage:
 * ~~~.bash
 * ThermodynamicsBQS <a> <b> <param_range_file> <outputfile>
 * ~~~
 * 
 * Where:
 * - <a> is the parameter a for the QvdW model (GeV fm^3)
 * - <b> is the parameter b for the QvdW model (fm^3)
 * - <param_range_file> is the file with the parameter range
 * - <outputfile> is the file to write the results to
 */