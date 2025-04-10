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
  { "T",   { 0.100, 0.180, 0.001 } },  // GeV
  { "muB", { 0.000, 0.000, 0.01  } }, // GeV
  { "muQ", { 0.000, 0.000, 0.01  } }, // GeV
  { "muS", { 0.000, 0.000, 0.01  } }  // GeV
};

// Read the parameter range from file
void ReadParameterRangeFromFile(const std::string &filename)
{
  if (filename == "")
    return;
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    std::cerr << "Error: cannot open file " << filename << "using default parameter range" << std::endl;
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

// Calculates the HRG model 2nd order susceptibilities and their temperature derivatives
// for a given range of T, muB, muQ, muS
// Usage: SusceptibilitiesBQS <param_range_file> <outputfile> <a> <b> <useRG>
// * <a> -- parameter a for the QvdW model (GeV fm^3)
// * <b> -- parameter b for the QvdW model (fm^3)
// * <useRG> -- Use real gas model and Carnahan-Starling EV
// * <param_range_file> -- file with the parameter range
// * <outputfile> -- file to write the results to
int main(int argc, char *argv[])
{
  // van der Waals attraction
  double a = 0.329;
  if (argc > 1)
    a = atof(argv[1]);
  
  // van der Waals repulsion
  double b = 3.42;
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

  // Compute isothermal speed of sound
  bool compute_cT2 = false;
  if (argc > 5)
    compute_cT2 = atoi(argv[5]);

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
  
  std::string outputfile = "Susceptibilities-" + ModelPrefix + "-output.dat";
  if (argc > 6)
    outputfile = argv[6];
  
  
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
  // Speed of sound squared
  fout << std::setw(tabsize) << "cs2" << " ";
  if (compute_cT2)
    fout << std::setw(tabsize) << "cT2" << " ";
  // Heat capacity
  fout << std::setw(tabsize) << "cV/T^3" << " ";
  // Susceptibilities
  fout << std::setw(tabsize) << "chi2B" << " ";
  fout << std::setw(tabsize) << "chi2Q" << " ";
  fout << std::setw(tabsize) << "chi2S" << " ";
  fout << std::setw(tabsize) << "chi11BQ" << " ";    
  fout << std::setw(tabsize) << "chi11BS" << " ";
  fout << std::setw(tabsize) << "chi11QS" << " ";
  // Temperature derivatives
  fout << std::setw(tabsize) << "T*chi2B'" << " ";
  fout << std::setw(tabsize) << "T*chi2Q'" << " ";
  fout << std::setw(tabsize) << "T*chi2S'" << " ";
  fout << std::setw(tabsize) << "T*chi11BQ'" << " ";
  fout << std::setw(tabsize) << "T*chi11BS'" << " ";
  fout << std::setw(tabsize) << "T*chi11QS'" << " ";
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
          model->CalculateFluctuations();
          iters++;
          
          
          // Susceptibilities
          double chi2B = model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge);
          double chi2Q = model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
          double chi2S = model->Susc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge);
          double chi11BQ = model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::ElectricCharge);
          double chi11BS = model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge);
          double chi11QS = model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge);

          model->CalculateTemperatureDerivatives();
          // Temperature derivatives
          double Tchi2Bpr = T * model->dSuscdT(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge);
          double Tchi2Qpr = T * model->dSuscdT(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
          double Tchi2Spr = T * model->dSuscdT(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge);
          double Tchi11BQpr = T * model->dSuscdT(ConservedCharge::BaryonCharge, ConservedCharge::ElectricCharge);
          double Tchi11BSpr = T * model->dSuscdT(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge);
          double Tchi11QSpr = T * model->dSuscdT(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge);

          // Speed of sound
          double cs2 = model->cs2();

          // Heat capacity
          double cV  = model->HeatCapacity();
          double cVT3 = cV / pow(T, 3.) / xMath::GeVtoifm3();

          // Print to file
          fout << setw(tabsize) << T << " ";
          fout << setw(tabsize) << muB << " ";
          fout << setw(tabsize) << muQ << " ";
          fout << setw(tabsize) << muS << " ";
          fout << setw(tabsize) << cs2 << " ";
          if (compute_cT2) {
            double cT2 = model->cT2();
            fout << setw(tabsize) << cT2 << " ";
          }
          fout << setw(tabsize) << cVT3 << " ";
          fout << setw(tabsize) << chi2B << " ";
          fout << setw(tabsize) << chi2Q << " ";
          fout << setw(tabsize) << chi2S << " ";
          fout << setw(tabsize) << chi11BQ << " ";
          fout << setw(tabsize) << chi11BS << " ";
          fout << setw(tabsize) << chi11QS << " ";
          fout << setw(tabsize) << Tchi2Bpr << " ";
          fout << setw(tabsize) << Tchi2Qpr << " ";
          fout << setw(tabsize) << Tchi2Spr << " ";
          fout << setw(tabsize) << Tchi11BQpr << " ";
          fout << setw(tabsize) << Tchi11BSpr << " ";
          fout << setw(tabsize) << Tchi11QSpr << " ";
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
 * \example SusceptibilitiesBQS.cpp
 * 
 * An example of calculating 2nd order susceptibilities and their temperature derivatives
 * for a given range of T, \f$\mu_B\f$, \f$\mu_Q\f$, and \f$\mu_S\f$.
 * 
 * For each specified value of T, \f$\mu_B\f$, \f$\mu_Q\f$, and \f$\mu_S\f$ does the following:
 * 
 *   1. Sets the temperature and chemical potentials.
 *   2. Calculates the following quantities:
 *      - Susceptibilities \f$\chi_2^B\f$, \f$\chi_2^Q\f$, \f$\chi_2^S\f$, \f$\chi_{11}^{BQ}\f$, \f$\chi_{11}^{BS}\f$, \f$\chi_{11}^{QS}\f$
 *      - Temperature derivatives of susceptibilities \f$T \chi_2^B'\f$, \f$T \chi_2^Q'\f$, \f$T \chi_2^S'\f$, \f$T \chi_{11}^{BQ}'\f$, \f$T \chi_{11}^{BS}'\f$, \f$T \chi_{11}^{QS}'\f$
 *      - Adiabatic \f$c_s^2\f$ and isothermal \f$c_T^2\f$ squared sound velocity
 *      - Heat capacity \f$c_V/T^3\f$
 * 
 * Usage:
 * ~~~.bash
 * SusceptibilitiesBQS <a> <b> <useRG> <param_range_file> <cT2flag> <outputfile>
 * ~~~
 * 
 * Where:
 * - <a> is the parameter a for the QvdW model (GeV fm^3)
 * - <b> is the parameter b for the QvdW model (fm^3)
 * - <useRG> is a flag to use the real gas model and Carnahan-Starling EV
 * - <param_range_file> is the file with the parameter range
 * - <cT2flag> is a flag to compute the isothermal speed of sound squared
 * - <outputfile> is the file to write the results to
 */
