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
#include "HRGEventGenerator.h"

#include "ThermalFISTConfig.h"

using namespace std;

#ifdef ThermalFIST_USENAMESPACE
using namespace thermalfist;
#endif

double muBss(double ss) {
  return 1.308 / (1. + 0.273 * ss);
}


double Tss(double ss) {
  double tmpmuB = muBss(ss);
  return 0.166 - 0.139 * tmpmuB * tmpmuB - 0.053 * tmpmuB * tmpmuB * tmpmuB * tmpmuB;
}


// Calculates the (mixed) susceptibility using the average decay procedure of 1402.1238
// by introducing auxiliary chemical potential for all decaying resonances
// ch1,2 -- vector of charges that all final state hadrons contribute to the observable 1,2
double CalculateAveragedDecaysChi2(ThermalModelBase *model, const std::vector<double> & ch1, const std::vector<double> & ch2)
{
  double ret = 0.;
  for (int r = 0; r < model->TPS()->ComponentsNumber(); ++r) {
    //const ThermalParticle &part_r = model->TPS()->Particle(r);

    const ThermalParticleSystem::ResonanceFinalStatesDistribution& decayDistributions = model->TPS()->ResonanceFinalStatesDistributions()[r];

    double mn1 = model->ScaledVariancePrimordial(r) * model->Densities()[r] / pow(model->Parameters().T, 3) / pow(xMath::GeVtoifm(), 3);

    double mn2 = 0., mn3 = 0.;
    for (int j = 0; j < model->TPS()->ComponentsNumber(); ++j) {
      //const ThermalParticle &part_j = model->TPS()->Particle(j);

      double njavr = 0.;
      for (size_t i = 0; i < decayDistributions.size(); ++i) {
        njavr += decayDistributions[i].first * decayDistributions[i].second[j];
      }

      mn2 += ch1[j] * njavr;
      mn3 += ch2[j] * njavr;
    }

    ret += mn1 * mn2 * mn3;
  }

  return ret;
}

// Collision energy dependence of 2nd order susceptibilities of
// (proxy) conserved charges, computed within the Ideal HRG model 
// along the phenomenological chemical freeze-out curve
// Comparison of analytic and Monte Carlo calculations
// Usage: cpc4mcHRG <withMonteCarlo> <nevents>
// where <withMonteCarlo> flag determines whether Monte Carlo calculations
// are performed and <nevents> is the number of Monte Carlo events per
// single collision energy
int main(int argc, char *argv[])
{
  int withMonteCarlo = 1;
  if (argc > 1)
    withMonteCarlo = atoi(argv[1]);

  int nevents = 100000;
  if (argc >2)
    nevents = atoi(argv[2]);
  
  string listname = string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2014/list.dat";
  ThermalParticleSystem parts(listname);

  // Disable K0 and K0bar decays to avoid strangeness non-conservation
  parts.ParticleByPDG(311).SetStable(true);
  parts.ParticleByPDG(-311).SetStable(true);

  // Optionally switch off some decays
  //for (size_t i = 0; i < parts.Particles().size(); ++i) {
  //  ThermalParticle& part = parts.Particle(i);
  //  
  //  // Switch off |S| = 1 hyperon decays
  //  if (abs(part.BaryonCharge()) == 1 && abs(part.Strangeness()) == 1)
  //    part.SetStable(true);
  //  
  //  // Switch off K*0 like decays
  //  if (abs(part.BaryonCharge()) == 0 && abs(part.Strangeness()) == 1 && part.ElectricCharge() == 0)
  //    part.SetStable(true);

  //  // Switch off all decays
  //  part.SetStable(true);
  //}

  ThermalModelParameters params;
  params.V   = 1000.;
  
  ThermalModelBase *model;

  model = new ThermalModelIdeal(&parts);

  model->SetParameters(params);
  model->SetStatistics(false);

  model->SetUseWidth(ThermalParticle::BWTwoGamma);
  //model->SetUseWidth(ThermalParticle::ZeroWidth);

  // Normalize branching ratios to avoid charge non-conservation in decays
  model->SetNormBratio(true);
  

  // For calculating net-kaon scaled variance using "average decays"
  vector<double> coefsk(parts.Particles().size(), 0.);
  coefsk[parts.PdgToId(321)]  =  1.;
  coefsk[parts.PdgToId(-321)] = -1.;


  double smin = 3.;
  double smax = 3000.;
  int iterss = 100;

  double wt1 = get_wall_time();
  int iters = 0;

  double musp = 0., muqp = 0., mucp = 0.;

  // First analytic calculations
  printf("Analytic\n");

  FILE *f = fopen("cpc4.analyt.dat", "w");

  fprintf(f, "%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n",
    "sqrts[GeV]",
    "T[MeV]",
    "muB[MeV]",
    "muQ[MeV]",
    "muS[MeV]",
    "c2k/c1k",
    "C_B,S",
    "C_p,k_fd",
    "C_Q,S",
    "C_Q,k_fd",
    "C_B,Q",
    "C_p,Q_fd",
    "c2k_av/c1k"
  );

  vector<double> ens, Ts, muBs, muQs, muSs;

  double logSmin = log(smin);
  double logSmax = log(smax);
  double dlogS = (logSmax - logSmin) / iterss;
  for (int is = 0; is <= iterss; ++is) {
    double ss = exp(logSmin + is * dlogS);
    double T = Tss(ss);
    double muB = muBss(ss);

    model->SetTemperature(T);
    model->ConstrainMuS(true);
    model->ConstrainMuQ(true);
    model->ConstrainMuC(true);
    model->SetBaryonChemicalPotential(muB);
    model->SetStrangenessChemicalPotential(musp);
    model->SetElectricChemicalPotential(muqp);
    model->SetCharmChemicalPotential(mucp);
    model->ConstrainChemicalPotentials();

    musp = model->Parameters().muS;
    muqp = model->Parameters().muQ;
    mucp = model->Parameters().muC;

    ens.push_back(ss);
    Ts.push_back(T);
    muBs.push_back(muB);
    muQs.push_back(muqp);
    muSs.push_back(musp);

    model->CalculateFluctuations();


    double chi1B = model->CalculateBaryonDensity() / pow(T, 3) / pow(xMath::GeVtoifm(), 3);
    double chi1Q = model->CalculateChargeDensity() / pow(T, 3) / pow(xMath::GeVtoifm(), 3);
    double chi1S = model->CalculateStrangenessDensity() / pow(T, 3) / pow(xMath::GeVtoifm(), 3);

    double chi1p = (model->GetDensity(2212, Feeddown::StabilityFlag) - model->GetDensity(-2212, Feeddown::StabilityFlag)) / pow(T, 3) / pow(xMath::GeVtoifm(), 3);
    double chi1k = (model->GetDensity(321, Feeddown::StabilityFlag) - model->GetDensity(-321, Feeddown::StabilityFlag)) / pow(T, 3) / pow(xMath::GeVtoifm(), 3);

    double chi2B   = model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge);
    double chi2Q   = model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    double chi2S   = model->Susc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge);
    double chi11BS = model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge);
    double chi11QS = model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge);
    double chi11BQ = model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::BaryonCharge);

    double chi2p   = model->ProxySusc(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge);
    double chi2k   = model->ProxySusc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge);
    double chi11pk = model->ProxySusc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge);
    double chi11Qk = model->ProxySusc(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge);
    double chi11pQ = model->ProxySusc(ConservedCharge::ElectricCharge, ConservedCharge::BaryonCharge);

    double chi2kav = CalculateAveragedDecaysChi2(model, coefsk, coefsk);

    printf("%15lf%15lf%15lf%15lf%15lf%15lf%15lf\n", 
      ss, 
      chi11BS / chi2S, 
      chi11QS / chi2S, 
      chi11BQ / chi2B, 
      chi11pk / chi2k, 
      chi11Qk / chi2k, 
      chi11pQ / chi2p);

    fprintf(f, "%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf\n",
      ss,
      T * 1.e3,
      muB * 1.e3,
      muqp * 1.e3,
      musp * 1.e3,
      chi2k / chi1k,
      chi11BS / chi2S,
      chi11pk / chi2k,
      chi11QS / chi2S,
      chi11Qk / chi2k,
      chi11BQ / chi2B,
      chi11pQ / chi2p,
      chi2kav / chi1k);

    fflush(f);

    iters++;
  }

  fclose(f);

  double wt2 = get_wall_time();

  printf("%30s %lf s\n", "Running time:", (wt2 - wt1));
  printf("%30s %lf s\n", "Time per single calculation:", (wt2 - wt1) / iters);

  printf("\n\n");

  if (!withMonteCarlo)
    return 0;

  printf("Monte Carlo\n");

  f = fopen("cpc4.montecarlo.dat", "w");

  fprintf(f, "%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n",
    "sqrts[GeV]",
    "T[MeV]",
    "muB[MeV]",
    "muQ[MeV]",
    "muS[MeV]",
    "c2k/c1k",
    "C_B,S",
    "C_p,k_fd",
    "C_Q,S",
    "C_Q,k_fd",
    "C_B,Q",
    "C_p,Q_fd"
  );

  // Now Monte Carlo
  iters = 0;
  wt1 = get_wall_time();
  RandomGenerators::SetSeed(1);

  for (size_t ind = 0; ind < ens.size(); ind++) {
    model->SetTemperature(Ts[ind]);
    model->SetBaryonChemicalPotential(muBs[ind]);
    model->SetElectricChemicalPotential(muQs[ind]);
    model->SetStrangenessChemicalPotential(muSs[ind]);

    double mc_B = 0., mc_Q = 0., mc_S = 0., mc_p = 0., mc_k = 0.;
    double mc_B2 = 0., mc_Q2 = 0., mc_S2 = 0., mc_BQ = 0., mc_QS = 0., mc_BS = 0.;
    double mc_p2 = 0., mc_k2 = 0., mc_pk = 0., mc_pQ = 0., mc_Qk = 0.;

    // Setup the configuration for event generator
    EventGeneratorConfiguration config;
    config.fEnsemble = EventGeneratorConfiguration::GCE;
    config.fModelType = EventGeneratorConfiguration::PointParticle;
    config.CFOParameters = model->Parameters();

    SphericalBlastWaveEventGenerator generator(model->TPS(), config, 0.100, 0.5);
    double wsum = 0.;
    for (int i = 0; i < nevents; ++i) {
      SimpleEvent ev = generator.GetEvent();
      int mcev_B = 0, mcev_Q = 0, mcev_S = 0, mcev_p = 0, mcev_k = 0;

      for (size_t part = 0; part < ev.Particles.size(); ++part) {
        long long pdgid = ev.Particles[part].PDGID;
        const ThermalParticle &thermpart = parts.ParticleByPDG(pdgid);
        mcev_B += thermpart.BaryonCharge();
        mcev_Q += thermpart.ElectricCharge();
        mcev_S += thermpart.Strangeness();

        if (pdgid == 2212)
          mcev_p += 1;
        else if (pdgid == -2212)
          mcev_p += (-1);

        if (pdgid == 321)
          mcev_k += 1;
        else if (pdgid == -321)
          mcev_k += (-1);
      }

      mc_B += ev.weight * mcev_B;
      mc_Q += ev.weight * mcev_Q;
      mc_S += ev.weight * mcev_S;
      mc_p += ev.weight * mcev_p;
      mc_k += ev.weight * mcev_k;

      mc_B2 += ev.weight * mcev_B * mcev_B;
      mc_Q2 += ev.weight * mcev_Q * mcev_Q;
      mc_S2 += ev.weight * mcev_S * mcev_S;
      mc_p2 += ev.weight * mcev_p * mcev_p;
      mc_k2 += ev.weight * mcev_k * mcev_k;

      mc_BQ += ev.weight * mcev_B * mcev_Q;
      mc_QS += ev.weight * mcev_Q * mcev_S;
      mc_BS += ev.weight * mcev_B * mcev_S;
      mc_pQ += ev.weight * mcev_p * mcev_Q;
      mc_pk += ev.weight * mcev_p * mcev_k;
      mc_Qk += ev.weight * mcev_Q * mcev_k;

      wsum += ev.weight;
    }

    double chi1B = mc_B / wsum;
    double chi1Q = mc_Q / wsum;
    double chi1S = mc_S / wsum;

    double chi1p = mc_p / wsum;
    double chi1k = mc_k / wsum;

    double chi2B = mc_B2 / wsum - mc_B / wsum * mc_B / wsum;
    double chi2Q = mc_Q2 / wsum - mc_Q / wsum * mc_Q / wsum;
    double chi2S = mc_S2 / wsum - mc_S / wsum * mc_S / wsum;
    double chi11BS = mc_BS / wsum - mc_B / wsum * mc_S / wsum;
    double chi11QS = mc_QS / wsum - mc_Q / wsum * mc_S / wsum;
    double chi11BQ = mc_BQ / wsum - mc_B / wsum * mc_Q / wsum;

    double chi2p = mc_p2 / wsum - mc_p / wsum * mc_p / wsum;
    double chi2k = mc_k2 / wsum - mc_k / wsum * mc_k / wsum;
    double chi11pk = mc_pk / wsum - mc_p / wsum * mc_k / wsum;
    double chi11Qk = mc_Qk / wsum - mc_Q / wsum * mc_k / wsum;
    double chi11pQ = mc_pQ / wsum - mc_p / wsum * mc_Q / wsum;

    printf("%15lf%15lf%15lf%15lf%15lf%15lf%15lf\n", ens[ind], chi11BS / chi2S, chi11QS / chi2S, chi11BQ / chi2B, chi11pk / chi2k, chi11Qk / chi2k, chi11pQ / chi2p);
  

    fprintf(f, "%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf\n",
      ens[ind],
      Ts[ind] * 1.e3,
      muBs[ind] * 1.e3,
      muQs[ind] * 1.e3,
      muSs[ind] * 1.e3,
      chi2k / chi1k,
      chi11BS / chi2S,
      chi11pk / chi2k,
      chi11QS / chi2S,
      chi11Qk / chi2k,
      chi11BQ / chi2B,
      chi11pQ / chi2p);

    fflush(f);
    iters++;
  }

  fclose(f);

  wt2 = get_wall_time();

  printf("%30s %lf s\n", "Running time:", (wt2 - wt1));
  printf("%30s %lf s\n", "Time per single calculation:", (wt2 - wt1) / iters);

  delete model;

  return 0;
}

/**
 * \example cpc4-mcHRG.cpp
 * 
 * Calculates the collision energy dependence of 2nd order susceptibilities of
 * (proxy) conserved charges, computed within the Ideal HRG model 
 * along the phenomenological chemical freeze-out curve.
 * 
 * Calculations are done in two steps:
 *   1. Analytically
 *   2. With Monte Carlo (if <withMonteCarlo> = 1)
 * 
 * 
 * 
 * Usage:
 * ~~~.bash
 * cpc4mcHRG <withMonteCarlo> <nevents>
 * ~~~
 * 
 * where <withMonteCarlo> flag determines whether Monte Carlo calculations
 * are performed and <nevents> is the number of Monte Carlo events per
 * single collision energy
 * 
 */