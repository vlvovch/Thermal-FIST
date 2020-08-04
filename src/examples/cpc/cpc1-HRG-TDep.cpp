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


// Temperature dependence of HRG thermodynamics at \mu = 0
// Three variants of the HRG model: 
// 1. Ideal HRG: <config> = 0
// 2. EV-HRG with constant radius parameter r = 0.3 fm for all hadrons (as in 1412.5478): <config> = 1
// 3. QvdW-HRG with a and b for baryons only, fixed to nuclear ground state (as in 1609.03975): <config> = 2
// Usage: cpc1HRGTDep <config>
int main(int argc, char *argv[])
{
  // Particle list file
  // Here we will use the list from THERMUS-2.3, for comparing the results with THERMUS-2.3
  string listname = string(ThermalFIST_INPUT_FOLDER) + "/list/thermus23/list.dat";

  // Alternative: use the default PDG2014 list
  //string listname = string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2014/list.dat";
  
  // Create the hadron list instance and read the list from file
  ThermalParticleSystem TPS(listname);

  // Which variant of the HRG model to use
  int config = 0;

  // Read config from command line
  if (argc > 1)
    config = atoi(argv[1]);

  
  string modeltype; // For output
  

  // Pointer to the thermal model instance used in calculations
  ThermalModelBase *model;

  

  if (config == 0) // Ideal HRG
  {
    model = new ThermalModelIdeal(&TPS);
    
    printf("#Calculating thermodynamics at \\mu = 0 in Id-HRG model\n");

    modeltype = "Id-HRG";
  }
  else if (config == 1) // EV-HRG, r = 0.3 fm, to reproduce 1412.5478
  {
    model = new ThermalModelEVDiagonal(&TPS);
    // Set r = 0.3 fm for each hadron in the list
    double rad = 0.3;
    for (int i = 0; i < model->TPS()->ComponentsNumber(); ++i)
      model->SetRadius(i, rad);

    printf("#Calculating thermodynamics at \\mu = 0 in EV-HRG model with r = %lf fm\n", rad);

    modeltype = "EV-HRG";
  }
  else if (config == 2) // QvdW-HRG, to reproduce 1609.03975
  {
    model = new ThermalModelVDWFull(&TPS);

    // vdW parameters, for baryon-baryon, antibaryon-antibaryon ONLY, otherwise zero
    double a = 0.329; // In GeV*fm3
    double b = 3.42;  // In fm3

    // Iterate over all pairs of hadron species and set a and b
    for (int i = 0; i < model->TPS()->ComponentsNumber(); ++i) {
      for (int j = 0; j < model->TPS()->ComponentsNumber(); ++j) {
        int B1 = model->TPS()->Particle(i).BaryonCharge(); // Baryon number of 1st species
        int B2 = model->TPS()->Particle(j).BaryonCharge(); // Baryon number of 2nd species

        if ((B1 > 0 && B2 > 0) || (B1 < 0 && B2 < 0))  // baryon-baryon or antibaryon-antibaryon
        {
          model->SetAttraction(i, j, a); // Set QvdW repulsion
          model->SetVirial(i, j, b);     // Set QvdW attraction
        }
        else // no vdW interactions for meson-meson, meson-baryon or baryon-anitbaryon pairs
        {
          model->SetAttraction(i, j, 0.);
          model->SetVirial(i, j, 0.);
        }
      }
    }

    printf("#Calculating thermodynamics at \\mu = 0 in QvdW-HRG model\n");

    modeltype = "QvdW-HRG";
  }
  else // Ideal HRG by default
  {
    model = new ThermalModelIdeal(&TPS);

    modeltype = "Id-HRG";
  }

  // Use quantum statistics
  model->SetStatistics(true);
  //model->SetStatistics(false);

  // Use mass integration over Breit-Wigner shapes in +-2Gamma interval, as in THERMUS
  model->SetUseWidth(ThermalParticle::BWTwoGamma);
  //model->SetUseWidth(ThermalParticle::ZeroWidth);

  // Set chemical potentials to zero
  model->SetBaryonChemicalPotential(0.0);
  model->SetElectricChemicalPotential(0.0);
  model->SetStrangenessChemicalPotential(0.0);
  model->SetCharmChemicalPotential(0.0);
  model->FillChemicalPotentials();

  
  // Prepare output
  char tmpc[1000];
  sprintf(tmpc, "cpc1.%s.TDep.out", modeltype.c_str());
  FILE *fout = fopen(tmpc, "w");

  // Output to screen
  printf("%15s%15s%15s%15s%15s%15s%15s\n",
    "T[MeV]", // Temperature
    "p/T^4",  // Scaled pressure
    "e/T^4",  // Scaled energy density
    "s/T^3",   // Scaled entropy density
    "chi2B",   // Baryon number susceptibility
    "chi4B",   // 4th order baryon cumulant
    "chi2B-chi4B"   // Difference of 2nd and 4th order baryon susceptibilities
  );

  // Output to file
  fprintf(fout, "%15s%15s%15s%15s%15s%15s%15s\n",
    "T[MeV]", // Temperature
    "p/T^4",  // Scaled pressure
    "e/T^4",  // Scaled energy density
    "s/T^3",   // Scaled entropy density
    "chi2B",   // Baryon number susceptibility
    "chi4B",   // 4th order baryon cumulant
    "chi2B-chi4B"   // Difference of 2nd and 4th order baryon susceptibilities
  );


  double wt1 = get_wall_time(); // Timing

  int iters = 0; // Number of data points

  // Temperature interval, in GeV
  double Tmin = 0.020;
  double Tmax = 0.2001;
  double dT   = 0.001;

  for (double T = Tmin; T <= Tmax; T += dT) {
    model->SetTemperature(T);

    // Calculates densities, solves all necessary transcendental equations, if necessary
    model->CalculateDensities();

    // Output temperature, scale to get the unit of MeV
    printf("%15lf", T * 1000.);
    fprintf(fout, "%15lf", T * 1000.);

    // Pressure, in [GeV/fm3]
    double p = model->CalculatePressure();
    // Scaled pressure, pow(xMath::GeVtoifm(), 3) needed to make it dimensionless
    double pT4 = p / pow(T, 4) / pow(xMath::GeVtoifm(), 3);
    printf("%15lf", pT4);
    fprintf(fout, "%15lf", pT4);

    // Energy density
    double e = model->CalculateEnergyDensity();
    double eT4 = e / pow(T, 4) / pow(xMath::GeVtoifm(), 3);
    printf("%15lf", eT4);
    fprintf(fout, "%15lf", eT4);

    // Entropy density
    double s = model->CalculateEntropyDensity();
    double sT3 = s / pow(T, 3) / pow(xMath::GeVtoifm(), 3);
    printf("%15lf", sT3);
    fprintf(fout, "%15lf", sT3);


    // Baryon number fluctuations

    // Vector containing baryon charge of each species
    vector<double> baryon_charges;
    for (int i = 0; i < model->TPS()->ComponentsNumber(); ++i) {
      baryon_charges.push_back(static_cast<double>(model->TPS()->Particle(i).BaryonCharge()));
    }

    // Calculate baryon number cumulants up to 4th order
    vector<double> chiB = model->CalculateChargeFluctuations(baryon_charges, 4);

    if (chiB.size() >= 4) {
      // chi2B
      printf("%15E", chiB[1]);
      fprintf(fout, "%15E", chiB[1]);

      // chi4B
      printf("%15E", chiB[3]);
      fprintf(fout, "%15E", chiB[3]);

      // chi2B - chi4B
      printf("%15E", chiB[1] - chiB[3]);
      fprintf(fout, "%15E", chiB[1] - chiB[3]);
    }

    printf("\n");
    fprintf(fout, "\n");

    iters++;

  }

  fclose(fout);

  delete model;
  

  double wt2 = get_wall_time();

  printf("%30s %lf s\n", "Running time:", (wt2 - wt1));
  printf("%30s %lf s\n", "Time per single calculation:", (wt2 - wt1) / iters);

  return 0;
}

/**
 * \example cpc1-HRG-TDep.cpp
 * 
 * Calculates the temperature dependence of HRG thermodynamics
 * at zero chemical potentials.
 * 
 * The calculated quantities include scaled pressure, scaled energy density,
 * scaled entropy density, the 2nd and 4th order baryon number susceptibilities.
 * 
 * Calculations can be done within three different models:
 *   - <config> = 0: Ideal HRG model
 *   - <config> = 1: EV-HRG model with a constant radius parameter r = 0.3 fm for all hadrons (as in 1412.5478) 
 *   - <config> = 2: QvdW-HRG with QvdW interactions between baryons only,
 *     with parameters being fixed to the nuclear ground state (as in 1609.03975) 
 * 
 * Usage:
 * ~~~.bash
 * cpc1HRGTDep <config>
 * ~~~
 * 
 */