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


// Calculates particle yields at a given T-mu values
// Usage: CalculationTmu
int main(int argc, char *argv[])
{
	// Fill the T-mu values where calculations should be performed
	vector<double> Tvalues, muvalues;

	// Here done by hand
	// Alternatively one can read those from external file
	// Note that all energy units are in GeV!
	// 1
	Tvalues.push_back(0.124425); muvalues.push_back(0.535912);
	// 2
	Tvalues.push_back(0.13344); muvalues.push_back(0.468765);
	// 3
	Tvalues.push_back(0.139899); muvalues.push_back(0.40872);
	// 4
	Tvalues.push_back(0.143476); muvalues.push_back(0.368669);
	// 5
	Tvalues.push_back(0.149036); muvalues.push_back(0.289944);
	// 6
	Tvalues.push_back(0.152653); muvalues.push_back(0.218109);
	// 7
	Tvalues.push_back(0.152849); muvalues.push_back(0.213352);


	// Create the hadron list instance and read the list from file
	//ThermalParticleSystem TPS(string(INPUT_FOLDER) + "/list/thermus23mod/list.dat");
	ThermalParticleSystem TPS(string(INPUT_FOLDER) + "/list/PDG2014/list.dat");


	// Create the ThermalModelIdeal instance (ideal HRG model)
	//ThermalModelIdeal model(&TPS);
	ThermalModelCanonical model(&TPS);

	// Use (or not) finite resonance width
	model.SetUseWidth(true);

	// Include (or not) quantum statistics
	model.SetStatistics(true);
	for (int i = 0; i < model.TPS()->Particles().size(); ++i) {
		ThermalParticle &part = model.TPS()->Particle(i);
		if (part.PdgId() != 211 && part.PdgId() != 111 && part.PdgId() != -211) {
			part.UseStatistics(false);
		}
	}
	model.CalculateQuantumNumbersRange();

	// Output, to write into file use, e.g., fprintf
	printf("%15s%15s%15s%15s%15s%15s%15s\n", 
		"T[GeV]", "muB[GeV]", "<K+>", "<pi+>", "<K+>/<pi+>",
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
		// Alternatively set the muS value
		//model.SetStrangenessChemicalPotential(0.);

		// Constrain muq from Q/B = 0.4 condition
		model.ConstrainMuQ(true);
		model.SetQoverB(0.4);
		// Alternatively set the muS value
		//model.SetElectricChemicalPotential(0.);

		// Chemical non-equilbrium parameters
		model.SetGammaq(1.);
		model.SetGammaS(1.);
		
		// Set volume
		model.SetVolumeRadius(3.); // in fm
		//model.SetVolume(5000.); //<-- Alternative, in fm^3


		// Determine muS and muQ
		model.FixParameters();

		// Calculate all hadron densities
		model.CalculateDensities();

		// Fluctuations
		model.CalculateFluctuations();

		// Calculate final yields, 2nd argument equal to 1 in GetDensity means that we take final yields, with feeddown
		double yieldKplus  = model.GetDensity(321, 1) * model.Volume();
		double yieldpiplus = model.GetDensity(211, 1) * model.Volume();

		double Nch    = model.ChargedMultiplicityFinal(0);
		double Nplus  = model.ChargedMultiplicityFinal(1);
		double Nminus = model.ChargedMultiplicityFinal(2);

		double wNch    = model.ChargedScaledVarianceFinal(0);
		double wNplus  = model.ChargedScaledVarianceFinal(1);
		double wNminus = model.ChargedScaledVarianceFinal(2);

		printf("%15lf%15lf%15lf%15lf%15lf%15lf%15lf\n",
			T,
			muB,
			yieldKplus,
			yieldpiplus,
			yieldKplus/yieldpiplus,
			Nminus,
			wNminus);
	}

	return 0;
}
