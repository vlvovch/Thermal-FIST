/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdio>

#include "HRGBase.h"
#include "HRGEV.h"
#include "HRGFit.h"
#include "HRGVDW.h"

#include "ThermalFISTConfig.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

// Calculates equation of state, particle yields and fluctuations at a given set of T-mu values
// Usage: CalculationTmu
int main(int argc, char *argv[])
{
	int ModelType = 1; // 0 - Ideal HRG, 1 - QvdW HRG (from 1609.03975)
	if (argc > 1)
		ModelType = atoi(argv[1]);

	std::string prefix = "QvdW-HRG";
	if (ModelType != 1)
		prefix = "IdealHRG";
	
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


	//ThermalParticleSystem TPS(string(ThermalFIST_INPUT_FOLDER) + "/list/thermus23mod/list.dat"); // <-- modified THERMUS-2.3 list
	ThermalParticleSystem TPS(string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2014/list.dat");  // <-- Default list, no light nuclei
	//ThermalParticleSystem TPS(string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2014/list-withnuclei.dat");  // <-- Default list, with light nuclei

	// Create the ThermalModel instance
	// Choose the class which fits the required variant of HRG model
	//ThermalModelIdeal model(&TPS);
	//ThermalModelCanonical model(&TPS);
	ThermalModelBase *model;
	if (ModelType == 1) // QvdW-HRG
	{
		model = new ThermalModelVDWFull(&TPS);

		// Set the QvdW interaction parameters
		// As in 1609.03975, here we consider (anti)baryon-(anti)baryon interactions only
		double a = 0.329;
		double b = 3.42;

		// Loop over all hadron-hadron pairs to set a and b for each of these pairs
		for (int i1 = 0; i1 < model->TPS()->Particles().size(); ++i1) {
			for (int i2 = 0; i2 < model->TPS()->Particles().size(); ++i2) {
				const ThermalParticle &part1 = model->TPS()->Particles()[i1];
				const ThermalParticle &part2 = model->TPS()->Particles()[i2];

				int B1 = part1.BaryonCharge();
				int B2 = part2.BaryonCharge();

				// Or use pdgid's to identify the two particles
				//int pdgid1 = part1.PdgId();
				//int pdgid2 = part2.PdgId();

				// Meson-meson, meson-baryon, baryon-antibaryon non-interacting
				if (!(B1 * B2 > 0)) {

					model->SetVirial(i1, i2, 0.); // No repulsion

					model->SetAttraction(i1, i2, 0.); // No attraction
					continue;
				}
				else {
					// BB excluded volume
					model->SetVirial(i1, i2, b);

					// BB attraction
					model->SetAttraction(i1, i2, a);
				}
			}
		}
	}
	else {
		model = new ThermalModelIdeal(&TPS);
	}

	// Use (or not) finite resonance width
	model->SetUseWidth(true);

	// Include (or not) quantum statistics
	model->SetStatistics(true);

	// Output, here on screen, to write into file use, e.g., fprintf
	printf("%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", 
		"T[GeV]", "muB[GeV]", 
		"P[GeV/fm3]", "e[GeV/fm3]", "s[fm-3]", 
		"<K+>", "<pi+>", "<K+>/<pi+>", 
		"w[K+]", "w[pi+]",
		"<N->", "w[N-]",
		"chi3B/chi2B", "chi4B/chi2B",
		"chi3Q/chi2Q", "chi4Q/chi2Q",
		"chi3S/chi2S", "chi4S/chi2S");

	// The same output to file
	std::string filename = prefix + ".CalculationTmu.dat";
	FILE *f = fopen(filename.c_str(), "w");
	fprintf(f, "%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n",
		"T[GeV]", "muB[GeV]",
		"P[GeV/fm3]", "e[GeV/fm3]", "s[fm-3]",
		"<K+>", "<pi+>", "<K+>/<pi+>",
		"w[K+]", "w[pi+]",
		"<N->", "w[N-]",
		"chi3B/chi2B", "chi4B/chi2B",
		"chi3Q/chi2Q", "chi4Q/chi2Q",
		"chi3S/chi2S", "chi4S/chi2S");

	// Iterate over all the T-muB pair values
	for (int i = 0; i < Tvalues.size(); ++i) {
		double T   = Tvalues[i];
		double muB = muvalues[i];


		// Set temperature and baryon chemical potential
		model->SetTemperature(T);
		model->SetBaryonChemicalPotential(muB);

		// Constrain muB from strangeness neutrality condition
		model->ConstrainMuS(true);
		// Alternatively set the muS value manually
		// model->ConstrainMuS(false);
		// model->SetStrangenessChemicalPotential(0.);

		// Constrain muq from Q/B = 0.4 condition
		model->ConstrainMuQ(true);
		model->SetQoverB(0.4);
		// Alternatively set the muQ value manually
		//model->ConstrainMuQ(false);
		//model->SetElectricChemicalPotential(0.);

		// Chemical non-equilbrium parameters
		model->SetGammaq(1.);
		model->SetGammaS(1.);
		
		// Set volume
		model->SetVolumeRadius(3.); // System radius R in fm, volume is V = (4/3) * \pi * R^3
		//model->SetVolume(5000.); //<-- Alternative, volume V in fm^3


		// Determine muS and/or muQ from constraints, if there are any
		model->FixParameters();

		// Calculate all hadron densities, both primordial and final
		model->CalculateDensities();

		// Calculate fluctuations
		model->CalculateFluctuations();

		// Equation of state parameters
		double p = model->CalculatePressure();       // Pressure in GeV/fm3
		double e = model->CalculateEnergyDensity();  // Energy density in GeV/fm3
		double s = model->CalculateEntropyDensity(); // Entropy density in fm-3

		// Calculate final yields, 
		// Usage: model->GetDensity(pdgid, feeddown)
		// pdgid -- PDG code for the desired particle species
		// feeddown: 0 - primordial, 1 - final, 2 - final with additional feeddown from weak decays
		// yield = density * volume
		double yieldKplus  = model->GetDensity(321, Feeddown::StabilityFlag) * model->Volume();
		double yieldpiplus = model->GetDensity(211, Feeddown::StabilityFlag) * model->Volume();

		// Scaled variance of final state particle number fluctuations
		double wKplus  = model->ScaledVarianceTotal( model->TPS()->PdgToId(321) );
		double wpiplus = model->ScaledVarianceTotal( model->TPS()->PdgToId(211) );

		// Charged particle mean multplicities, after decays
		// Argument: 0 - all charged, 1 - positively charged, 2 - negatively charged
		double Nch    = model->ChargedMultiplicityFinal(0);
		double Nplus  = model->ChargedMultiplicityFinal(1);
		double Nminus = model->ChargedMultiplicityFinal(2);

		// Scaled variance for charged particle multplicity distribution, after decays 
		double wNch    = model->ChargedScaledVarianceFinal(0);
		double wNplus  = model->ChargedScaledVarianceFinal(1);
		double wNminus = model->ChargedScaledVarianceFinal(2);

		// Higher-order fluctuations of conserved charges B, Q, S

		// Array of charges of all particles in the list
		// E.g., if the ith particle has baryon charge 1, then chargesB[i] = 1 etc. 
		vector<double> chargesB(model->Densities().size()), chargesQ(model->Densities().size()), chargesS(model->Densities().size());

		// Array with the values of the calculated susceptibilities
		vector<double> chchis;

		// Baryon number
		for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
			chargesB[i] = model->TPS()->Particles()[i].BaryonCharge();
		}

		// Calculation of susceptibilities chi1-chi4
		chchis = model->CalculateChargeFluctuations(chargesB, 4);
		double chi2B = chchis[1];
		double chi3B = chchis[2];
		double chi4B = chchis[3];

		// Electric charge, same procedure
		for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
			chargesQ[i] = model->TPS()->Particles()[i].ElectricCharge();
		}

		chchis = model->CalculateChargeFluctuations(chargesQ, 4);
		double chi2Q = chchis[1];
		double chi3Q = chchis[2];
		double chi4Q = chchis[3];

		// Strangeness
		for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
			chargesS[i] = model->TPS()->Particles()[i].Strangeness();
		}

		chchis = model->CalculateChargeFluctuations(chargesS, 4);
		double chi2S = chchis[1];
		double chi3S = chchis[2];
		double chi4S = chchis[3];


		printf("%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf\n",
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
			wNminus,
			chi3B / chi2B,
			chi4B / chi2B,
			chi3Q / chi2Q,
			chi4Q / chi2Q,
			chi3S / chi2S,
			chi4S / chi2S);

		fprintf(f, "%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf\n",
			T,
			muB,
			p,
			e,
			s,
			yieldKplus,
			yieldpiplus,
			yieldKplus / yieldpiplus,
			wKplus,
			wpiplus,
			Nminus,
			wNminus,
			chi3B / chi2B,
			chi4B / chi2B,
			chi3Q / chi2Q,
			chi4Q / chi2Q,
			chi3S / chi2S,
			chi4S / chi2S);
	}

	fclose(f);

	delete model;

	return 0;
}

/**
 * \example CalculationTmu.cpp
 * 
 * An example of calculating various thermodynamic quantities at fixed
 * T and \f$\mu_B\f$.
 * 
 * For each specified value of a \f$T-\mu_B\f$ pair does the following.
 * 
 *   1. Constrains the chemical potentials \f$\mu_Q\f$ and \f$\mu_S\f$ from
 *      the conditions Q/B = 0.4 and S = 0
 *   2. Calculates the following quantities:
 *      - Pressure
 *      - Energy density
 *      - Entropy density
 *      - K+ yield (with feeddown)
 *      - pi+ yield (with feeddown)
 *      - K+/pi+ ratio (with feeddown)
 *      - Scaled variances of K+ and pi+
 *      - Multiplicity and scaled variance of negatively charged particles
 *      - Skewness and kurtosis of net baryon, charge, and strangeness fluctuations
 * 
 * Usage:
 * ~~~.bash
 * CalculationTmu <ModelType>
 * ~~~
 * 
 * Where <ModelType> is 0 for Ideal HRG, and 1 for QvdW-HRG (from [arXiv:1609.03975](https://arxiv.org/abs/1609.03975))
 */