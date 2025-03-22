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

// The range in muB to scan
map<string, double> params = {
	{ "muBmin", 0.940 }, // GeV
  { "muBmax", 1.500 }, // GeV
  { "muBstep", 0.010 }, // GeV
  { "include_leptons", 1}, // Include leptons in the calculation
  {"a", 0.160}, // van der Waals attraction   (default value from https://arxiv.org/pdf/2109.06799)
  {"b", 2.236} // Carnahan-Starling repulsion (default value from https://arxiv.org/pdf/2109.06799)
};

// Read the parameters from file
void ReadParametersFromFile(const std::string &filename)
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
		double val;
		if (!(iss >> val)) { break; } // error
    params[param] = val;
	}

	fin.close();
}

// Calculates the neutron star matter EoS in beta-equilibrium
// for a given range of muB
// Usage: NeutronStars-CSHRG <param_file> <outputfile>
// * <param_file> -- file with the parameters
// * <outputfile> -- file to write the results to
int main(int argc, char *argv[])
{
  // Parameter range from file
  string params_file = "";
  if (argc > 1)
    params_file = argv[1];
  ReadParametersFromFile(params_file);

  bool include_leptons = (params["include_leptons"] != 0.);

	std::string outputfile = "NSMatter-CSHRG";
  if (include_leptons)
    outputfile += "-leptons";
  outputfile += "-output.dat";
	if (argc > 2)
		outputfile = argv[2];


  vector<string> lists = {string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list.dat"};

  // If including leptons, add the list of charged leptons
  if (include_leptons)
    lists.push_back(string(ThermalFIST_INPUT_FOLDER) + "/list/electroweak/list-charged-leptons.dat");
	
	// Create the hadron list instance and read the list from file
	ThermalParticleSystem TPS(lists);

	// Create the ThermalModel instance, we will use real gas
	// Choose the class which fits the required variant of HRG model
	ThermalModelRealGas *model;
	model = new ThermalModelRealGas(&TPS);
  {
    ExcludedVolumeModelBase *evmod = new ExcludedVolumeModelCS();
    double b = params["b"];
    if (b != 0.) {
      auto vdWb = GetBaryonBaryonInteractionMatrix(model->TPS(), b);
      model->SetExcludedVolumeModel(new ExcludedVolumeModelCrosstermsGeneralized(evmod, vdWb));
    }

    double a = params["a"];
    if (a != 0.) {
      auto vdWa = GetBaryonBaryonInteractionMatrix(model->TPS(), a);
      static_cast<ThermalModelRealGas *>(model)->SetMeanFieldModel(new MeanFieldModelMultiVDW(vdWa));
    }
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
  fout << std::setw(tabsize) << "nB[n0]" << " ";
  fout << std::setw(tabsize) << "nQ[n0]" << " ";     // Electri charge density (should be zero)
  fout << std::setw(tabsize) << "P[GeV/fm3]" << " "; // Pressure
  fout << std::setw(tabsize) << "e[GeV/fm3]" << " "; // Energy density
  fout << std::setw(tabsize) << "s[GeV/fm3]" << " "; // Entropy density (should be zero)
  fout << std::setw(tabsize) << "(1/3-p/e)" << " "; // trace anomaly
  // fout << std::setw(tabsize) << "vs2" << " "; // sound velocity squared
  fout << std::setw(tabsize) << "vs2" << " "; // adiabatic sound velocity squared
  fout << std::setw(tabsize) << "vT2" << " "; // isothermal sound velocity squared
  fout << std::setw(tabsize) << "Yp" << " ";  // proton fraction
  fout << std::setw(tabsize) << "Yn" << " ";  // neutron fraction
  fout << std::setw(tabsize) << "YSig-" << " ";  // proton fraction
  fout << std::setw(tabsize) << "YLambda" << " ";  // proton fraction
  fout << std::setw(tabsize) << "Ye" << " ";  // electron-to-baryon
  fout << std::setw(tabsize) << "Ymu" << " ";  // muon-to-baryon
	fout << std::setw(tabsize) << std::endl;

  const double n0 = 0.16; // fm^-3

	// Timer
  double wt1 = get_wall_time();
  int iters = 0; // number of iterations

  // Beta-equilibrium conditions
  model->SetTemperature(0.);
  model->SetQoverB(0.);
  model->ConstrainMuS(false);
  model->SetStrangenessChemicalPotential(0.);

  // "Good" initial guess for muQ for muB = 940 MeV
  model->SetElectricChemicalPotential(-0.018);

  // For vs2
  double pprev = 0., eprev = 0.;

	// Loop over the chemical potential
  double muBmin = params["muBmin"];
  double muBmax = params["muBmax"];
  double dmuB   = params["muBstep"];
  for(double muB = muBmin; muB <= muBmax + 0.1 * dmuB; muB += dmuB) {
    model->SetBaryonChemicalPotential(muB);

    // Solve for beta-equilibrium (use previous muQ as initial guess)
    model->ConstrainChemicalPotentials(false);

    // Sometimes Broyden does not converge from first try (try to repeat and improve)
    {
      double rhoQ = model->ElectricChargeDensity() / n0;
      int max_repeats = 5;
      for(int i = 0; i < max_repeats && abs(rhoQ) > 1.e-10; i++) {
        model->ConstrainChemicalPotentials(false);
        rhoQ = model->ElectricChargeDensity() / n0;
      }
    }

    // Compute the densities
    model->CalculatePrimordialDensities();
    iters++;

    // Collect the output

    double T = model->Parameters().T; // Temperature in GeV
    double muQ = model->Parameters().muQ; // Electric charge chemical potential in GeV
    double muS = model->Parameters().muS; // Strangeness chemical potential in GeV

    double p = model->Pressure();       // Pressure in GeV/fm3
    double e = model->EnergyDensity();  // Energy density in GeV/fm3
    double s = model->EntropyDensity(); // Entropy density in fm-3
    double rhoB = model->BaryonDensity();         // Baryon density in fm-3
    double rhoQ = model->ElectricChargeDensity(); // Electric charge density in fm-3
    double rhoS = model->StrangenessDensity();    // Strangeness density in fm-3
    double trace_anomaly = (1./3. - p/e);         // Trace anomaly

    double Yp = model->GetDensity(2212, Feeddown::Primordial) / rhoB; // proton fraction
    double Yn = model->GetDensity(2112, Feeddown::Primordial) / rhoB; // neutron fraction
    double YSig = model->GetDensity(3112, Feeddown::Primordial) / rhoB; // Sigma- fraction
    double YLam = model->GetDensity(3122, Feeddown::Primordial) / rhoB; // Lambda fraction
    double Ye = 0.;
    double Ymu = 0.;
    if (include_leptons) {
      Ye  = model->GetDensity(11, Feeddown::Primordial) / rhoB; // electron-to-baryon
      Ymu = model->GetDensity(13, Feeddown::Primordial) / rhoB; // muon-to-baryon
    }


    // vs2 = rhoB / muB * dmub / drhob
    // For strangeness-neutral: vs2 = rhoB / muB / (chi2B + chi11BQ * dmuQ/dmuB)
    model->CalculateFluctuations();
    double chi2B = model->SusceptibilityDimensionfull(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge);
    double chi11BQ = model->SusceptibilityDimensionfull(ConservedCharge::BaryonCharge, ConservedCharge::ElectricCharge);
    double chi2Q = model->SusceptibilityDimensionfull(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    // double vs2 = (rhoB / xMath::GeVtoifm3()) / muB / (chi2B - chi11BQ * chi11BQ / chi2Q); // Cross-check
    // if (chi11BQ == 0.)
    //   vs2 = (rhoB / xMath::GeVtoifm3()) / muB / (chi2B);
    // // Finite difference (backward) for vs2
    // double vs2FD = (p - pprev) / (e - eprev);

    // Speed of sound
    double vs2fct = model->cs2(true, true, false); // vs2 through the standard function
    double vT2fct = model->cT2(true, true, false); // vT2 through the standard function

    // Print to file
    fout << setw(tabsize) << T << " ";
    fout << setw(tabsize) << muB << " ";
    fout << setw(tabsize) << muQ << " ";
    fout << setw(tabsize) << muS << " ";
    fout << setw(tabsize) << rhoB / n0 << " ";
    fout << setw(tabsize) << rhoQ / n0 << " ";
    fout << setw(tabsize) << p << " ";
    fout << setw(tabsize) << e << " ";
    fout << setw(tabsize) << s << " ";
    fout << setw(tabsize) << trace_anomaly << " ";
    // fout << setw(tabsize) << vs2 << " "; // vs2 direct
    fout << setw(tabsize) << vs2fct << " "; // vs2 standard
    fout << setw(tabsize) << vT2fct << " ";
    fout << setw(tabsize) << Yp << " ";
    fout << setw(tabsize) << Yn << " ";
    fout << setw(tabsize) << YSig << " ";
    fout << setw(tabsize) << YLam << " ";
    fout << setw(tabsize) << Ye << " ";
    fout << setw(tabsize) << Ymu << " ";
    fout << setw(tabsize) << endl;

    pprev = p;
    eprev = e;
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
 * \example NeutronStars-CSHRG.cpp
 * 
 * An example of calculating the neutron star matter equation of state (EoS)
 * in beta-equilibrium over a range of baryon chemical potentials (\f$\mu_B\f$).
 * 
 * The EoS is calculated using an interacting HRG model with mean field attraction
 * and Carnahan-Starling excluded volume repulsion. 
 * Model parameters are based on Fujimoto et al., Phys. Lett. B 835 (2022) 137524  (https://arxiv.org/pdf/2109.06799)
 * 
 * For each specified value of \f$\mu_B\f$, the following quantities are calculated:
 *   - Temperature
 *   - Baryon chemical potential
 *   - Electric charge chemical potential
 *   - Strangeness chemical potential
 *   - Baryon density
 *   - Electric charge density
 *   - Pressure
 *   - Energy density
 *   - Entropy density
 *   - Trace anomaly
 *   - Sound velocity squared (adiabatic)
 *   - Sound velocity squared (isothermal)
 *   - Proton fraction
 *   - Neutron fraction
 *   - Sigma- fraction
 *   - Lambda fraction
 *   - Electron-to-baryon ratio
 *   - Muon-to-baryon ratio
 * 
 * Usage:
 * ~~~.bash
 * NeutronStars-CSHRG <param_file> <outputfile>
 * ~~~
 * 
 * Where:
 * - <param_file> is the file with the input parameters
 * - <outputfile> is the file to write the results to
 */