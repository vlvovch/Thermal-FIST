#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdio>

#include "HRGBase.h"
#include "HRGEV.h"
#include "HRGFit.h"

#include "ThermalFISTConfig.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

// Fits within Bag Model, as in 1512.08046
// First the global fit
// Then, the chi2 profile
// Usage: BagModelFit <rp> <Tmin> <Tmax> <dT>
int main(int argc, char *argv[])
{
	// Proton radius
	double radProton = 0.50;
	if (argc>1) 
		radProton = atof(argv[1]);

	// Minimum temperature for chi2 profile calculation [GeV]
	double Tmin = 0.100;
	if (argc>2)
		Tmin = atof(argv[2]);

	// Maximum temperature for chi2 profile calculation [GeV]
	double Tmax = 0.301;
	if (argc>3)
		Tmin = atof(argv[3]);

	// Step in temperature for chi2 profile calculation [GeV]
	double dT = 0.001;
	if (argc>4)
		dT = atof(argv[4]);

	// Create the hadron list instance and read the list from file
	ThermalParticleSystem TPS(string(INPUT_FOLDER) + "/list/thermus23/list.dat");
	//ThermalParticleSystem TPS(string(INPUT_FOLDER) + "/list/PDG2014/list.dat");

	// Create the EV-HRG ThermalModel instance
	ThermalModelEVDiagonal model(&TPS);

	// Use finite resonance width
	model.SetUseWidth(true);

	// Include quantum statistics
	model.SetStatistics(true);

	// All chemical potentials are zero
	model.SetBaryonChemicalPotential(0.);
	model.SetStrangenessChemicalPotential(0.);
	model.SetElectricChemicalPotential(0.);
	model.SetCharmChemicalPotential(0.);

	// Loop over all particles and set their EV parameters in accordance with the bag model
	for (int i = 0; i < model.TPS()->Particles().size(); ++i) {
		const ThermalParticle &part = model.TPS()->Particles()[i];
		double mass = part.Mass(); // Mass of particle i
		double rad  = radProton * pow(mass / xMath::mnucleon(), 1. / 3.); // Radius parameter of particle i
		model.SetRadius(i, rad); // Sets the radius
	}

	// Load the experimental data
	vector<FittedQuantity> quantities = ThermalModelFit::loadExpDataFromFile(string(INPUT_FOLDER) + "/data/ALICE-PbPb2.76TeV-0-5-1512.08046.dat");
		

	// Global fit
	ThermalModelFit pfit(&model);

	// By default T, muB, and R parameters are fitted, while others (gammaS etc.) are fixed
	// Initial parameters values are taken from those currently set in a ThermalModel object
	// Here we do not fit muB, which is set to zero
	pfit.SetParameterFitFlag("muB", false); // Do not fit muB
	
	// R is fitted by default
	// We can still specify the initial value, the initial delta used by minuit,
	// and the lower and upper limits using the SetParameter function
	double Rinit  = 10.0;
	double Rdelta = 1.0;
	double Rmin = 0.0;
	double Rmax = 30.0;
	pfit.SetParameter("R", Rinit, Rdelta, Rmin, Rmax);

	pfit.SetQuantities(quantities); // Provide the data to be fitted

	pfit.PerformFit();  // Performs the fit
	pfit.PrintFitLog(); // Print the fit result


	// Now perform the chi-square profile scan
	printf("%15s%15s%15s\n", "T [GeV]", "R [fm]", "chi2");
	pfit.SetParameterFitFlag("T", false); // We are not fitting T anymore

	// Iterate over all the T values
	for (double T = Tmin; T <= Tmax; T += dT) {
		pfit.SetParameterValue("T", T); // Set the temperature
		ThermalModelFitParameters result = pfit.PerformFit(false);  // We still have to fit the radius, the argument suppresses the output during minimization  
		printf("%15lf%15lf%15lf\n", T, result.R.value, result.chi2);
	}

	return 0;
}
