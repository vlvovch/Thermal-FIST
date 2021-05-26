/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2020 Volodymyr Vovchenko
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
#include "HRGPCE.h"

#include "ThermalFISTConfig.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

// This is an example of doing PCE-HRG model calculations at the LHC energies using Thermal-FIST
// Usage: PCE-Saha-LHC
int main(int argc, char *argv[])
{
	// The default particle list. As of version 1.3 this is PDG2020 list including light nuclei
	ThermalParticleSystem parts(ThermalFIST_DEFAULT_LIST_FILE);
	
	// To include excited nuclei use the following line instead
	//ThermalParticleSystem parts(string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list-withexcitednuclei.dat");

	// To reproduce arXiv:1903.10024 use the PDG2014 list
	//ThermalParticleSystem TPS(string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2014/list-withnuclei.dat"); 

	// Use ideal HRG model
	ThermalModelIdeal model(&parts);
	
	// PCE-HRG model
	ThermalModelPCE modelpce(&model);
	modelpce.UseSahaForNuclei(true); // Light nuclei evaluated using the Saha equation (arXiv:1903.10024)
	modelpce.FreezeLonglivedResonances(false); // All strongly decaying resonance are in partial equilibrium

	// Chemical freeze-out conditions: 2.76 TeV 0-10% Pb-Pb collisions
	ThermalModelParameters params_chemical_freezeout;
	params_chemical_freezeout.T = 0.155; // Temperature in GeV
	params_chemical_freezeout.muB = 0.;
	params_chemical_freezeout.V = 4700.; // Volume in fm^3

	model.SetParameters(params_chemical_freezeout);

	// For finite baryon density: constrain muQ and muS
	model.ConstrainChemicalPotentials();
	params_chemical_freezeout = model.Parameters();

	model.FillChemicalPotentials(); // Fills chemical potentials for all species at Tch

	// Set the chemical freeze-out as an "initial" condition for PCE
	modelpce.SetChemicalFreezeout(params_chemical_freezeout);

	// The list of chemical potentials for output, coded by the pdg code
	vector<long long> pdgcodes_stable;
	pdgcodes_stable.push_back(211);  // pions (pi+)
	pdgcodes_stable.push_back(321);  // kaons (K+)
	pdgcodes_stable.push_back(2212); // protons (p+)
	pdgcodes_stable.push_back(3122); // Lambdas
	pdgcodes_stable.push_back(3222); // Sigma+
	pdgcodes_stable.push_back(3312); // Xi-
	pdgcodes_stable.push_back(3334); // Omega

	// The list of yield ratios to output
	vector< pair<long long, long long> > ratios;
	// First the nuclei
	ratios.push_back(make_pair(1000010020, 2212)); // d/p
	ratios.push_back(make_pair(1000020030, 2212)); // He3/p
	ratios.push_back(make_pair(1000010030, 2212)); // H3/p
	ratios.push_back(make_pair(1000020040, 2212)); // He4/p
	ratios.push_back(make_pair(1010010030, 2212)); // Hypertriton/p
	ratios.push_back(make_pair(1010010040, 2212)); // HyperHydrogen4/p
	// Now the resonances
	ratios.push_back(make_pair(313, -321));    // K^*0 / K^-
	ratios.push_back(make_pair(113, 211));     // rho^0/ pi^+
	ratios.push_back(make_pair(3124, 3122));   // \Lambda(1520)/\Lambda
	ratios.push_back(make_pair(9010221, 211)); // f0(980) / pi^+
	ratios.push_back(make_pair(2224, 2212));   // \Delta(1232)++/p

	// Preparing the output files
	// The file to output the parameters (volume, entropy, chemical potentials)
	FILE* fout_params = fopen("PCE.LHC.Parameters.dat", "w");
	fprintf(fout_params, "%15s %15s %15s ", "T[MeV]", "V/Vch", "S/Sch");
	for (int i = 0; i < pdgcodes_stable.size(); ++i) {
		fprintf(fout_params, "%15s ", ("mu_" + string(parts.ParticleByPDG(pdgcodes_stable[i]).Name())).c_str());
	}
	fprintf(fout_params, "\n");

	// The file to output the yield ratios
	FILE* fout_ratios = fopen("PCE.LHC.Ratios.dat", "w");
	fprintf(fout_ratios, "%15s ", "T[MeV]");
	for (int i = 0; i < ratios.size(); ++i) {
		fprintf(fout_ratios, "%15s ", (parts.ParticleByPDG(ratios[i].first).Name() + "/" + parts.ParticleByPDG(ratios[i].second).Name()).c_str());
	}
	fprintf(fout_ratios, "\n");

	// The temperature scan
	double T0 = params_chemical_freezeout.T;
	double dT = 0.001;   // steps of 1 MeV
	double Tmin = 0.070; // Down to 70 MeV

	// Store the value of the total entropy at the chemical freeze-out
	double entropy_chemical_freezeout = modelpce.ThermalModel()->EntropyDensity() * params_chemical_freezeout.V;

	// Loop over temperatures
	for (double T = T0; T >= Tmin - 1.e-9; T -= dT) {
		printf("T = %lf MeV\n", T * 1.e3);

		// Compute the PCE chemical potentials at a given temperature
		modelpce.CalculatePCE(T);

		// Output the parameters at the current temperature
		fprintf(fout_params, "%15lf %15lf %15lf ",
			T * 1.e3,
			modelpce.ThermalModel()->Volume() / params_chemical_freezeout.V,
			modelpce.ThermalModel()->EntropyDensity() * modelpce.ThermalModel()->Volume() / entropy_chemical_freezeout
			);

		for (int i = 0; i < pdgcodes_stable.size(); ++i) {
			fprintf(fout_params, "%15lf ",
				modelpce.ChemicalPotentials()[ parts.PdgToId(pdgcodes_stable[i]) ]
			);
		}
		fprintf(fout_params, "\n");

		// Output the yield ratios at the current temperature
		fprintf(fout_ratios, "%15lf ", T * 1.e3);
		for (int i = 0; i < ratios.size(); ++i) {
			fprintf(fout_ratios, "%15E ",
				modelpce.ThermalModel()->GetYield(ratios[i].first, Feeddown::Electromagnetic) /
				modelpce.ThermalModel()->GetYield(ratios[i].second, Feeddown::Electromagnetic));
		}
		fprintf(fout_ratios, "\n");
		
	}

	fclose(fout_params);
	fclose(fout_ratios);

	return 0;
}

/**
 * \example PCE-Saha-LHC.cpp
 *
 * An example of doing partial chemical equilibrium HRG model calculations at the LHC energies using Thermal-FIST
 *
 * Calculates the evolution of the non-equilibrium chemical potentials (fugacities) and various particle ratios
 * in the hadronic phase of 0-10% central 2.76 TeV Pb-Pb collisions at the LHC.
 *
 * The calculations closely correspond to the results published in [arXiv:1903.10024](https://arxiv.org/abs/1903.10024)
 *
 * Calculations start at T<sub>ch</sub> = 155 MeV and go down to a specified temperature (by default down to 70 MeV in steps of 1 MeV).
 * The values of the chemical potentials, as well as of the system volume relative to the volume at the freeze-out, are 
 * written to a file `PCE.LHC.Parameters.dat'
 *
 * The particle yield ratios at each temperature are written to a file `PCE.LHC.Ratios.dat'.
 * 
 * The abundances of light nuclei are calculated using the Saha equation.
 *
 * The source code can be modified to obtain other particle yields, to change the particle list or
 * or the HRG model type (e.g. an excluded volume HRG instead of an ideal HRG), or to explore other collision energies.
 *
 * Usage:
 * ~~~.bash
 * ./PCE-Saha-LHC
 * ~~~
 *
 */