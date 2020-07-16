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


// Temperature dependence of the fit to ALICE 2.76 TeV data, 0-5% centrality, as in 1512.08046
// Four variants of HRG model: 
// 1. Ideal HRG: <config> = 0
// 2. Diagonal EV-HRG with bag model parametrization r = r_p * (m/m_p)^1/3, where r_p = 0.5 is proton radius parameter (as in 1512.08046): <config> = 1
// 3. Diagonal EV-HRG with constant radius parameter r = 0.3 fm for all baryons and r = 0 for all mesons (as in 1201.0693): <config> = 2
// 4. QvdW-HRG with a and b for baryons only, fixed to nuclear ground state (as in 1609.03975): <config> = 3
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

    printf("#Fitting 2.76 TeV ALICE data at \\mu = 0 in Id-HRG model\n");

    modeltype = "Id-HRG";
  }
  else if (config == 1) // EV-HRG, r = 0.3 fm for baryons, r = 0 for mesons, no light nuclei, as in 1512.08046
  {
    model = new ThermalModelEVDiagonal(&TPS);
    double rad = 0.3;
    for (int i = 0; i < model->TPS()->ComponentsNumber(); ++i) {
      if (model->TPS()->Particle(i).BaryonCharge() == 0) model->SetRadius(i, 0.);
      else model->SetRadius(i, rad);
    }

    printf("#Fitting 2.76 TeV ALICE data at \\mu = 0 in EV-HRG model with r = %lf fm for baryons, and r = 0 for mesons\n", rad);

    modeltype = "EV-HRG-TwoComponent";
  }
  else if (config == 2) // EV-HRG, Bag Model with r_p = 0.5 fm, as in 1512.08046
  {
    model = new ThermalModelEVDiagonal(&TPS);
    double radProton = 0.5;
    for (int i = 0; i < model->TPS()->ComponentsNumber(); ++i) {
      model->SetRadius(i, radProton * pow(model->TPS()->Particle(i).Mass() / 0.938, 1/3.));
    }

    printf("#Fitting 2.76 TeV ALICE data at \\mu = 0 in Bag Model EV-HRG model with proton r = %lf fm\n", radProton);

    modeltype = "EV-HRG-BagModel";
  }
  else if (config == 3) // QvdW-HRG, to reproduce 1609.03975
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

    printf("#Fitting 2.76 TeV ALICE data at \\mu = 0 in QvdW-HRG model\n");

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

  // Prepare fitter
  ThermalModelFit fitter(model);

  // By default T, muB, and R parameters are fitted, while others (gammaS etc.) are fixed
  // Initial parameters values are taken from those currently set in the ThermalModel object
  // Here we do not fit muB, which is set to zero above
  fitter.SetParameterFitFlag("muB", false); // Do not fit muB

  

  // R is fitted by default
  // We can still specify the initial value, the initial delta used by minuit,
  // and the lower and upper limits using the SetParameter function
  double Rinit  = 10.0;
  double Rdelta = 1.0;
  double Rmin   = 0.0;
  double Rmax   = 30.0;
  fitter.SetParameter("R", Rinit, Rdelta, Rmin, Rmax);

  // Load the experimental data
  vector<FittedQuantity> quantities = ThermalModelFit::loadExpDataFromFile(string(ThermalFIST_INPUT_FOLDER) + "/data/ALICE-PbPb2.76TeV-0-5-1512.08046.dat");
  fitter.SetQuantities(quantities);


  printf("%15s%15s%15s%15s\n",
    "T[MeV]",   // Temperature in MeV
    "R[fm]",    // System radius in fm
    "chi2",     // chi_2
    "chi2_dof"  // Reduced chi2
  );

  
  // Prepare for output to file
  char tmpc[1000];
  sprintf(tmpc, "cpc2.%s.ALICE2_76.chi2.TDep.out", modeltype.c_str());
  FILE *fout = fopen(tmpc, "w");

  fprintf(fout, "%15s%15s%15s%15s\n",
    "T[MeV]",   // Temperature in MeV
    "R[fm]",    // System radius in fm
    "chi2",     // chi_2
    "chi2_dof"  // Reduced chi2
  );

  double wt1 = get_wall_time(); // Timing

  int iters = 0; // Number of data points

  // Temperature interval, in GeV
  double Tmin = 0.100;
  double Tmax = 0.2501;
  double dT   = 0.002;

  if (config == 0)
    dT = 0.001;

  if (config == 2) {
    dT = 0.005;
    Tmax = 0.4001;
  }

  for (double T = Tmin; T <= Tmax; T += dT) {
    // We also do not fit T, but fix it at each iteration to a given value
    fitter.SetParameterFitFlag("T", false);
    fitter.SetParameterValue("T", T); // Set the temperature

    ThermalModelFitParameters result = fitter.PerformFit(false);  // We still have to fit the radius, the argument suppresses the output during minimization  

    double Rfit = result.R.value;
    double chi2 = result.chi2;
    double chi2dof = result.chi2ndf;

    printf("%15lf%15lf%15lf%15lf\n", T * 1000., Rfit, chi2, chi2 / (result.ndf - 1.));

    fprintf(fout, "%15lf%15lf%15lf%15lf\n", T * 1000., Rfit, chi2, chi2 / (result.ndf - 1.));

    iters++;

    fflush(stdout);

  }

  fclose(fout);

  delete model;


  double wt2 = get_wall_time();

  printf("%30s %lf s\n", "Running time:", (wt2 - wt1));
  printf("%30s %lf s\n", "Time per single calculation:", (wt2 - wt1) / iters);

  return 0;
}


/**
 * \example cpc2-chi2-vs-T.cpp
 * 
 * Calculates the temperature profile of \f$ \chi^2 \f$ of a fit to
 * the ALICE 2.76 TeV data, 0-5% centrality, as in 1512.08046
 * 
 * Calculations can be done within four variants of the HRG model:
 *   - <config> = 0: Ideal HRG model
 *   - <config> = 1: Diagonal EV-HRG model with the bag model radii parametrization r = r_p * (m/m_p)^1/3, 
 *     where r_p = 0.5 is proton radius parameter (as in 1512.08046) 
 *   - <config> = 2: Diagonal EV-HRG model with a constant radius parameter r = 0.3 fm for all baryons 
 *     and r = 0 for all mesons (as in 1201.0693)
 *   - <config> = 3: QvdW-HRG with QvdW interactions between baryons only,
 *     with parameters being fixed to the nuclear ground state (as in 1609.03975) 
 * 
 * Usage:
 * ~~~.bash
 * cpc2chi2 <config>
 * ~~~
 * 
 */