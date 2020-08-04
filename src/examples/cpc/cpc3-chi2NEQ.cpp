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


// Thermal fits to NA49 and LHC data within the Ideal HRG model in
// 1. Chemical equilibrium scenario (<config> = 0);
// 2. Chemical non-equilibrium scenario (<config> = 1).
// 
// Usage: cpc3chi2NEQ <config>
int main(int argc, char *argv[])
{
  // Particle list file
  // Here we will use the list from THERMUS-2.3, for comparing the results with THERMUS-2.3
  //string listname = string(ThermalFIST_INPUT_FOLDER) + "/list/thermus23/list.dat";

  // Alternative: use the default PDG2014 list
  string listname = string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2014/list.dat";
  //string listname = string(ThermalFIST_INPUT_FOLDER) + "/../../input/list/PDG2014update/list.dat";

  // Create the hadron list instance and read the list from file
  ThermalParticleSystem TPS(listname);

  // Which variant of the HRG model to use
  int config = 0;

  // Read config from command line
  if (argc > 1)
    config = atoi(argv[1]);


  string fittype; // For output
  if (config == 1)
    fittype = "NEQ";
  else
    fittype = "EQ";

                    // Pointer to the thermal model instance used in calculations
  ThermalModelBase *model;

  model = new ThermalModelIdeal(&TPS);


  // Use quantum statistics
  model->SetStatistics(true);
  model->SetCalculationType(IdealGasFunctions::Quadratures); // Use quadratures for better accuracy, needed due to Bose-Einstein condensation effects for pions
  //model->SetStatistics(false);

  // Use mass integration over Breit-Wigner shapes in +-2Gamma interval, as in THERMUS
  model->SetUseWidth(ThermalParticle::BWTwoGamma);
  //model->SetUseWidth(ThermalParticle::ZeroWidth);

  // Prepare for output
  char tmpc[1000];
  sprintf(tmpc, "cpc3.%s.chi2.out", fittype.c_str());
  FILE *fout = fopen(tmpc, "w");


  // Prepare fitter
  ThermalModelFit fitter(model);


  


  printf("%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n",
    "Dataset",  // What data are fitted
    "T[MeV]",   // Temperature in MeV
    "T_err[MeV]",   // Temperature error in MeV
    "muB[MeV]", // Baryon chemical potential in MeV
    "muB_err[MeV]", // Baryon chemical potential error in MeV
    "R[fm]",    // System radius in fm
    "R_err[fm]",    // System radius error in fm
    "gammaq",    // gamma_q
    "gammaq_err",    // gamma_q
    "gammaS",    // gamma_S
    "gammaS_err",    // gamma_S
    "chi2",     // chi_2
    "chi2/dof"  // Reduced chi2
  );

  fprintf(fout, "%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n",
    "Dataset",  // What data are fitted
    "T[MeV]",   // Temperature in MeV
    "muB[MeV]", // Baryon chemical potential in MeV
    "R[fm]",    // System radius in fm
    "gammaq",    // gamma_q
    "gammaS",    // gamma_S
    "chi2",     // chi_2
    "chi2/dof",  // Reduced chi2
    "Q/B",      // Electric-to-baryon charge ratio, must be 0.4
    "S/|S|"     // Strangeness-to-absolutestrangeness ratio, must be very close to 0
  );


  // Experimental data sets
  vector<string> names;
  vector<string> filenames;

  names.push_back("NA49-30GeV-4pi");
  filenames.push_back(string(ThermalFIST_INPUT_FOLDER) + "/data/NA49/NA49-PbPb30AGeV-4pi.dat");

  names.push_back("NA49-40GeV-4pi");
  filenames.push_back(string(ThermalFIST_INPUT_FOLDER) + "/data/NA49/NA49-PbPb40AGeV-4pi.dat");

  names.push_back("NA49-80GeV-4pi");
  filenames.push_back(string(ThermalFIST_INPUT_FOLDER) + "/data/NA49/NA49-PbPb80AGeV-4pi.dat");

  names.push_back("NA49-158GeV-4pi");
  filenames.push_back(string(ThermalFIST_INPUT_FOLDER) + "/data/NA49/NA49-PbPb158AGeV-4pi.dat");

  names.push_back("ALICE-2_76-0-5");
  filenames.push_back(string(ThermalFIST_INPUT_FOLDER) + "/data/ALICE-PbPb2.76TeV-0-5-1512.08046.dat");





  double wt1 = get_wall_time(); // Timing

  int iters = 0; // Number of data sets

  for(size_t ind = 0; ind < names.size(); ++ind)
  {
    // Load the data to be fitted
    vector<FittedQuantity> quantities = ThermalModelFit::loadExpDataFromFile(filenames[ind]);
    fitter.SetQuantities(quantities);

    // Temperature is fitted by default
    // We can specify initial guess, as well as the bounds
    double Tinit  = 0.140;
    double Tdelta = 0.030;
    double Tmin   = 0.050;
    double Tmax   = 0.170; // If Tmax is larger, than fit converges to a local minimum with higher chi2 for NA49-158 GeV
    fitter.SetParameter("T", Tinit, Tdelta, Tmin, Tmax);

    // muB is fitted by default
    // We can specify initial guess, as well as the bounds
    double muBinit  = 0.250;
    double muBdelta = 0.150;
    double muBmin   = -0.050;
    double muBmax   = 1.000;
    fitter.SetParameter("muB", muBinit, muBdelta, muBmin, muBmax);


    // R is fitted by default
    // We can specify initial guess, as well as the bounds
    double Rinit = 10.0;
    double Rdelta = 1.0;
    double Rmin = 0.0;
    double Rmax = 30.0;
    fitter.SetParameter("R", Rinit, Rdelta, Rmin, Rmax);


    // gammaq is NOT fitted by default
    // We can specify initial guess, as well as the bounds
    double gqinit  = 1.0;
    double gqdelta = 0.6;
    double gqmin   = 0.001;
    double gqmax   = 1.800;
    if (config != 0) {
      fitter.SetParameterFitFlag("gammaq", true);
      fitter.SetParameter("gammaq", gqinit, gqdelta, gqmin, gqmax);
    }


    // gammaS is NOT fitted by default, the default value is 1
    // We can specify initial guess, as well as the bounds
    if (config != 0) {
    double gSinit  = 1.0;
    double gSdelta = 0.6;
    double gSmin   = 0.001;
    double gSmax   = 3.000;
      fitter.SetParameterFitFlag("gammaS", true);
      fitter.SetParameter("gammaS", gSinit, gSdelta, gSmin, gSmax);
    }


    ThermalModelFitParameters result = fitter.PerformFit(true);  // The argument is set to show additional (verbose) output during minimization  

    double Tfit = result.T.value;
    double Terr = result.T.error;
    double muBfit = result.muB.value;
    double muBerr = result.muB.error;
    double Rfit = result.R.value;
    double Rerr = result.R.error;
    double gqfit = result.gammaq.value;
    double gqerr = result.gammaq.error;
    double gSfit = result.gammaS.value;
    double gSerr = result.gammaS.error;
    double chi2 = result.chi2;
    double chi2dof = result.chi2ndf;

    printf("%15s%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf", 
      names[ind].c_str(),
      Tfit * 1000., 
      Terr * 1000.,
      muBfit * 1000.,
      muBerr * 1000.,
      Rfit, 
      Rerr,
      gqfit,
      gqerr,
      gSfit,
      gSerr,
      chi2, 
      chi2dof);

    printf("\n");


    fprintf(fout, "%15s%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15E%15E\n",
      names[ind].c_str(),
      Tfit * 1000.,
      muBfit * 1000.,
      Rfit,
      gqfit,
      gSfit,
      chi2,
      chi2dof,
      model->CalculateChargeDensity() / model->CalculateBaryonDensity(),
      model->CalculateStrangenessDensity() / model->CalculateAbsoluteStrangenessDensity());

    fflush(fout);

    iters++;

  }

  delete model;


  double wt2 = get_wall_time();

  printf("%30s %lf s\n", "Running time:", (wt2 - wt1));
  printf("%30s %lf s\n", "Time per single calculation:", (wt2 - wt1) / iters);

  return 0;
}

/**
 * \example cpc3-chi2NEQ.cpp
 * 
 * Thermal fits to NA49 and LHC data within the Ideal HRG model in  
 *   1. Chemical equilibrium scenario (<config> = 0)
 *   2. Chemical non-equilibrium scenario (<config> = 1)
 * 
 * Usage:
 * ~~~.bash
 * cpc3chi2NEQ <config>
 * ~~~
 * 
 */