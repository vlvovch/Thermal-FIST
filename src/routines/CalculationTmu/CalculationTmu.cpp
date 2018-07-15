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


// Calculates equation of state, particle yields and fluctuations at a given set of T-mu values
// Usage: CalculationTmu
int main(int argc, char *argv[])
{
	// Fill the T-mu values where calculations should be performed
	vector<double> Tvalues, muvalues;

	// Here done by hand
	// Alternatively one can read those from external file, or populate in a loop, etc.
	// Note that all energy units are in GeV!
	// 1
	Tvalues.push_back(0.100); muvalues.push_back(0.600);
	// 2
	Tvalues.push_back(0.130); muvalues.push_back(0.500);
	// 3
	Tvalues.push_back(0.160); muvalues.push_back(0.000);


	// Create the hadron list instance and read the list from file
	//ThermalParticleSystem TPS(string(INPUT_FOLDER) + "/list/thermus23mod/list.dat");
	ThermalParticleSystem TPS(string(INPUT_FOLDER) + "/list/PDG2014/list.dat");


	// Create the ThermalModel instance
	// Choose the class which fits the required variant of HRG model
	ThermalModelIdeal model(&TPS);
	//ThermalModelCanonical model(&TPS);

	// Use (or not) finite resonance width
	model.SetUseWidth(true);

	// Include (or not) quantum statistics
	model.SetStatistics(true);

	// Output, here on screen, to write into file use, e.g., fprintf
	printf("%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", 
		"T[GeV]", "muB[GeV]", 
		"P[GeV/fm3]", "e[GeV/fm3]", "s[fm-3]", 
		"<K+>", "<pi+>", "<K+>/<pi+>", 
		"w[K+]", "w[pi+]",
		"<N->", "w[N-]");

	// Iterate over all the T-muB pair values
	for (int i = 0; i < Tvalues.size(); ++i) {
		double T   = Tvalues[i];
		double muB = muvalues[i];


		// Set temperature and baryon chemical potential
		model.SetTemperature(T);
		model.SetBaryonChemicalPotential(muB);

		// Constrain muB from strangeness neutrality condition
		model.ConstrainMuS(true);
		// Alternatively set the muS value manually
		// model.ConstrainMuS(false);
		// model.SetStrangenessChemicalPotential(0.);

		// Constrain muq from Q/B = 0.4 condition
		model.ConstrainMuQ(true);
		model.SetQoverB(0.4);
		// Alternatively set the muQ value manually
		//model.ConstrainMuQ(false);
		//model.SetElectricChemicalPotential(0.);

		// Chemical non-equilbrium parameters
		model.SetGammaq(1.);
		model.SetGammaS(1.);
		
		// Set volume
		model.SetVolumeRadius(3.); // System radius R in fm, volume is V = (4/3) * \pi * R^3
		//model.SetVolume(5000.); //<-- Alternative, volume V in fm^3


		// Determine muS and/or muQ from constraints, if there are any
		model.FixParameters();

		// Calculate all hadron densities, both primordial and final
		model.CalculateDensities();

		// Calculate fluctuations
		model.CalculateFluctuations();

		// Equation of state parameters
		double p = model.CalculatePressure();       // Pressure in GeV/fm3
		double e = model.CalculateEnergyDensity();  // Energy density in GeV/fm3
		double s = model.CalculateEntropyDensity(); // Entropy density in fm-3

		// Calculate final yields, 
		// Usage: model.GetDensity(pdgid, feeddown)
		// pdgid -- PDG code for the desired particle species
		// feeddown: 0 - primordial, 1 - final, 2 - final with additional feeddown from weak decays
		// yield = density * volume
		double yieldKplus  = model.GetDensity(321, 1) * model.Volume();
		double yieldpiplus = model.GetDensity(211, 1) * model.Volume();

		// Scaled variance of final state particle number fluctuations
		double wKplus  = model.ScaledVarianceTotal( model.TPS()->PdgToId(321) );
		double wpiplus = model.ScaledVarianceTotal( model.TPS()->PdgToId(211) );

		// Charged particle mean multplicities, after decays
		// Argument: 0 - all charged, 1 - positively charged, 2 - negatively charged
		double Nch    = model.ChargedMultiplicityFinal(0);
		double Nplus  = model.ChargedMultiplicityFinal(1);
		double Nminus = model.ChargedMultiplicityFinal(2);

		// Scaled variance for charged particle multplicity distribution, after decays 
		double wNch    = model.ChargedScaledVarianceFinal(0);
		double wNplus  = model.ChargedScaledVarianceFinal(1);
		double wNminus = model.ChargedScaledVarianceFinal(2);

		printf("%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf\n",
			T,
			muB,
			p,
			e,
			s,
			yieldKplus,
			yieldpiplus,
			yieldKplus/yieldpiplus,
			wKplus,
			wpiplus,
			Nminus,
			wNminus);
	}

	return 0;
}
