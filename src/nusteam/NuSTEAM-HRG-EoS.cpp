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

// Function which sets the a and b parameters of the van der Waals HRG model 
// for baryon-baryon and antibaryon-antibaryon pairs
void SetVDWInteractionParameters(ThermalModelBase *model, double a, double b)
{
  // Iterate over all hadron pairs
  for (int i1 = 0; i1 < model->TPS()->Particles().size(); ++i1) {
    for (int i2 = 0; i2 < model->TPS()->Particles().size(); ++i2) {
      // Baryon charge of first species
      int B1 = model->TPS()->Particles()[i1].BaryonCharge();
      // Baryon charge of second species
      int B2 = model->TPS()->Particles()[i2].BaryonCharge();

      // TODO part:
      // Comment the lines below and write code to
      // set the vdW parameters for baryon-baryon and antibaryon-antibaryon pairs
      // based on the values (signs) of B1 and B2
      // Use model->SetAttraction(i1, i2, a) for attraction and
      // model->SetRepulsion(i1, i2, b) for repulsion
      cout << "Implement setting the vdW parameters for baryon-baryon and antibaryon-antibaryon interaction" << endl;
      exit(1);
    }
  }
}


// Temperature dependence of HRG thermodynamics and baryon susceptibilities at \mu = 0
// Calculations within van der Waals HRG model with parameters a and b taken as command line input
// For a = 0 and b = 0, the model reduces to the ideal HRG
// Usage: NuSTEAM-HRG-EoS <a> <b>
// a is in GeV * fm^3 and b is in fm^3
int main(int argc, char *argv[])
{
  // Particle list file
  // Here we will use the default PDG2020 list without light nuclei
  string listname = string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list.dat";
  
  // Create the hadron list instance and read the list from file
  ThermalParticleSystem TPS(listname);

  // van der Waals interaction parameters (zero by default)
  double a = 0; // GeV * fm^3
  // Read the value from command line
  if (argc > 1)
    a = atof(argv[1]);

  double b = 0; // fm^3
  // Read the value from command line
  if (argc > 2)
    b = atof(argv[2]);

  // Pointer to the thermal model instance used in calculations
  ThermalModelBase *model;
  

  if (a == 0. && b == 0.) // Ideal HRG
  {
    model = new ThermalModelIdeal(&TPS);
    
    cout << "#Calculating thermodynamics at \\mu = 0 in Id-HRG model" << endl;
  }
  else // vdW-HRG, to reproduce 1609.03975
  {
    model = new ThermalModelVDW(&TPS);

    SetVDWInteractionParameters(model, a, b);

    cout << "#Calculating thermodynamics at \\mu = 0 in vdW-HRG model" << endl;
  }

  // Use quantum statistics
  model->SetStatistics(true);
  //model->SetStatistics(false);

  // Use zero resonance widths
  model->SetUseWidth(ThermalParticle::ZeroWidth);
  // To use finite resonance widths (energy-dependent Breit-Wigner) uncomment below
  //model->SetUseWidth(ThermalParticle::eBW);

  // Set chemical potentials to zero
  model->SetBaryonChemicalPotential(0.0);
  model->SetElectricChemicalPotential(0.0);
  model->SetStrangenessChemicalPotential(0.0);
  model->SetCharmChemicalPotential(0.0);
  model->FillChemicalPotentials();

  // Prepare output
  const int tabsize = 15;
  cout << setw(tabsize) << "T[MeV]"  // Temperature
       << setw(tabsize) << "p/T^4"   // Scaled pressure
       << setw(tabsize) << "e/T^4"   // Scaled energy density
       << setw(tabsize) << "s/T^3"   // Scaled entropy density
       << setw(tabsize) << "chi2B"   // Baryon number susceptibility
       << setw(tabsize) << "chi4B"   // 4th order baryon cumulant
       << setw(tabsize) << "chi4B/chi2B" // Baryon number kurtosis
       << endl;

  double wt1 = get_wall_time(); // Timing

  int iters = 0; // Number of calculation points (to be incremented)

  // Temperature interval, in GeV
  double Tmin = 0.100;
  double Tmax = 0.2001;
  double dT   = 0.001;

  // Loop over temperature values
  for (double T = Tmin; T <= Tmax; T += dT) {
    // Set the desired tempetarure
    model->SetTemperature(T);

    // Calculates hadron densities, solving all necessary non-linear equations as necessary
    model->CalculateDensities();

    // Output temperature, scale to get the unit of MeV
    cout << setw(tabsize) << T * 1000.;

    // Pressure, in [GeV/fm3]
    double p = model->CalculatePressure();
    // Scaled pressure, pow(xMath::GeVtoifm(), 3) needed to make it dimensionless
    double pT4 = p / pow(T, 4) / pow(xMath::GeVtoifm(), 3);
    cout << setw(tabsize) << pT4;

    // Energy density
    double e = model->CalculateEnergyDensity();
    double eT4 = e / pow(T, 4) / pow(xMath::GeVtoifm(), 3);
    cout << setw(tabsize) << eT4;

    // Entropy density
    double s = model->CalculateEntropyDensity();
    double sT3 = s / pow(T, 3) / pow(xMath::GeVtoifm(), 3);
    cout << setw(tabsize) << sT3;


    // Baryon number fluctuations

    // Vector containing baryon charge of each species
    vector<double> baryon_charges = model->TPS()->GetConservedChargesVector(ConservedCharge::BaryonCharge);
    // Calculate baryon number cumulants up to 4th order
    vector<double> chiBs = model->CalculateChargeFluctuations(baryon_charges, 4);

    double chi2B = chiBs[1];
    double chi4B = chiBs[3];

    // Output baryon susceptibilities
    cout << setw(tabsize) << chi2B;
    cout << setw(tabsize) << chi4B;
    cout << setw(tabsize) << chi4B / chi2B;

    cout << endl;

    iters++;

  }

  // Clean-up
  delete model;
  

  // Time metrics
  double wt2 = get_wall_time();

  printf("%30s %lf s\n", "Running time:", (wt2 - wt1));
  printf("%30s %lf s\n", "Time per single calculation:", (wt2 - wt1) / iters);

  return 0;
}

/**
 * \example NuSTEAM-HRG-EoS.cpp
 * 
 * Calculates the temperature dependence of HRG thermodynamics
 * at zero chemical potentials.
 * 
 * The calculated quantities include scaled pressure, scaled energy density,
 * scaled entropy density, the 2nd and 4th order baryon number susceptibilities.
 * 
 * Calculations are done within vdW-HRG model (https://arxiv.org/abs/1609.03975) with interaction parameters from command-line
 *   <a> - vdW attraction parameter in GeV * fm^3 for baryons (GeV * fm^3)
 *   <b> - vdW repulsion parameter in GeV * fm^3 for baryons (fm^3)
 * 
 * Usage:
 * ~~~.bash
 * NuSTEAM-HRG-EoS <a> <b>
 * ~~~
 * 
 */