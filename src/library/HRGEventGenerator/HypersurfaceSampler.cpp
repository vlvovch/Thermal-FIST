/*
 * Thermal-FIST package
 *
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

#include <iostream>
#include <stdexcept>

#include "HRGEventGenerator/SimpleParticle.h"
#include "HRGEventGenerator/ParticleDecaysMC.h"
#include "HRGEventGenerator/HypersurfaceSampler.h"
#include "HRGBase/NumericalIntegration.h"
#include "HRGEV/ExcludedVolumeHelper.h"

using namespace std;

namespace thermalfist {


  RandomGenerators::VolumeElementSampler::VolumeElementSampler(const ParticlizationHypersurface* Hypersurface)
  {
    FillProbabilities(Hypersurface);
  }


 RandomGenerators::VolumeElementSampler::VolumeElementSampler(const std::vector<double>& Weights)
  {
    FillProbabilities(Weights);
  }

  void RandomGenerators::VolumeElementSampler::FillProbabilities(const ParticlizationHypersurface* Hypersurface)
  {
    std::vector<double> Weights;
    for (auto& element : (*Hypersurface)) {
      double dVeff = 0.;
      for (int mu = 0; mu < 4; ++mu)
        dVeff += element.dsigma[mu] * element.u[mu];
      Weights.push_back(dVeff);
    }
    FillProbabilities(Weights);
  }

  void RandomGenerators::VolumeElementSampler::FillProbabilities(const std::vector<double>& Weights)
  {
    m_CumulativeProbabilities = std::vector<double>(Weights.size(), 0.);
    double totalWeight = 0.;
    for (int i = 0; i < Weights.size(); ++i) {
      totalWeight += Weights[i];
      m_CumulativeProbabilities[i] += Weights[i];
      if (i > 0)
        m_CumulativeProbabilities[i] += m_CumulativeProbabilities[i - 1];
    }
    for (double& prob : m_CumulativeProbabilities) {
      prob /= totalWeight;
    }
  }

  int RandomGenerators::VolumeElementSampler::SampleVolumeElement(MTRand& rangen) const
  {
    double prob = rangen.rand();
    const auto& it = std::lower_bound(m_CumulativeProbabilities.begin(), m_CumulativeProbabilities.end(), prob);
    int tind = std::distance(m_CumulativeProbabilities.begin(), it);
    if (tind < 0) tind = 0;
    if (tind >= static_cast<int>(m_CumulativeProbabilities.size())) tind = m_CumulativeProbabilities.size() - 1;
    return tind;
  }


  RandomGenerators::HypersurfaceMomentumGenerator::HypersurfaceMomentumGenerator
  (const ParticlizationHypersurface* hypersurface,
    const ThermalParticle* particle,
    const VolumeElementSampler* positionsampler,
    const HypersurfaceMomentumGeneratorConfiguration& config) :
    m_ParticlizationHypersurface(hypersurface),
    m_Particle(particle),
    m_VolumeElementSampler(positionsampler),
    m_EtaSmear(config.etaSmear),
    m_ShearCorrection(config.shearCorrection),
    m_BulkCorrection(config.bulkCorrection),
    m_SpeedOfSoundSquared(config.speedOfSoundSquared)
  {

  }

  std::vector<double> RandomGenerators::HypersurfaceMomentumGenerator::GetMomentum(double mass) const
  {
    if (m_VolumeElementSampler == NULL || m_ParticlizationHypersurface == NULL) {
      throw std::runtime_error("RandomGenerators::HypersurfaceMomentumGenerator::GetMomentum(double mass): Hypersurface not initialized!");
    }

    if (mass < 0.)
      mass = Mass();

    int VolumeElementIndex = m_VolumeElementSampler->SampleVolumeElement();

    const ParticlizationHypersurfaceElement& elem = (*m_ParticlizationHypersurface)[VolumeElementIndex];

    return SamplePhaseSpaceCoordinateFromElement(&elem, m_Particle, mass, EtaSmear(), ShearCorrection(), BulkCorrection(), SpeedOfSoundSquared());
  }

  HypersurfaceEventGenerator::HypersurfaceEventGenerator(ThermalParticleSystem* TPS, const EventGeneratorConfiguration& config, const ParticlizationHypersurface* hypersurface, double etasmear, bool shear_correction, bool bulk_correction, double speed_of_sound_squared) :
    EventGeneratorBase()
  {
    SetConfiguration(TPS, config);
    SetHypersurface(hypersurface);
    SetEtaSmear(etasmear);
    SetRescaleTmu();
    SetShearCorrection(shear_correction);
    SetBulkCorrection(bulk_correction);
    SetSpeedOfSoundSquared(speed_of_sound_squared);
    //SetParameters(hypersurface, m_THM, etasmear);
  }


  HypersurfaceEventGenerator::HypersurfaceEventGenerator(ThermalParticleSystem *TPS,
                                                         const EventGeneratorConfiguration &config,
                                                         const ParticlizationHypersurface *hypersurface,
                                                         const RandomGenerators::HypersurfaceMomentumGenerator::HypersurfaceMomentumGeneratorConfiguration &configMomentumGenerator
                                                         ) : EventGeneratorBase(), m_MomentumGeneratorConfig(configMomentumGenerator)
  {
    SetConfiguration(TPS, config);
    SetHypersurface(hypersurface);
    SetMomentumGeneratorConfig(configMomentumGenerator);
    SetRescaleTmu();
  }

  std::vector<double> HypersurfaceEventGenerator::GCEMeanYields() const
  {
    return m_FullSpaceYields;
  }

  std::vector<double>& HypersurfaceEventGenerator::GCEMeanYields()
  {
    return m_FullSpaceYields;
  }

  //void HypersurfaceEventGenerator::SetParameters(const ParticlizationHypersurface* hypersurface, ThermalModelBase* model, double etasmear)
  void HypersurfaceEventGenerator::SetParameters()
  {
    if (m_RescaleTmu) {
      // Rescale T, mu's, and P for all cells
      RescaleHypersurfaceParametersEdens(const_cast<ParticlizationHypersurface *>(m_ParticlizationHypersurface),
                                         m_THM,
                                         m_edens);
    }
    ProcessVolumeElements();
    SetMomentumGenerators();
    m_ParametersSet = true;
  }

  void HypersurfaceEventGenerator::ProcessVolumeElements()
  {
    // Densities for the volume element sampling
    vector<vector<double>> allweights(m_THM->TPS()->ComponentsNumber(), vector<double>(m_ParticlizationHypersurface->size(), 0.));

    m_FullSpaceYields = vector<double>(m_THM->TPS()->ComponentsNumber(), 0.);
    m_Tav = 0.;
    m_Musav = vector<double>(m_THM->TPS()->ComponentsNumber(), 0.);
    double Veff = 0.;

    // For the standard sampling
    vector<double> FullDensities(m_THM->TPS()->ComponentsNumber(), 0.), FullDensitiesIdeal(m_THM->TPS()->ComponentsNumber(), 0.);

    cout << "Processing " << m_ParticlizationHypersurface->size() << " volume elements" << endl;


    // If rescaling T & mu, keep track how much of the energy-weighted volume correspond to given energy density
    // This is to catch the case when the energy density is not uniform or mismatched in the input parameters
    double rescaleTmu_EMatch = 0., rescaleTmu_Etot = 0.;

    // Process all the hypersurface elements
    for (size_t ielem = 0; ielem < m_ParticlizationHypersurface->size(); ++ielem) {
      if (ielem % 10000 == 0) {
        cout << ielem << " ";
        cout.flush();
      }
      const auto& elem = m_ParticlizationHypersurface->operator[](ielem);

      m_THM->SetTemperature(elem.T);
      m_THM->SetBaryonChemicalPotential(elem.muB);
      m_THM->SetElectricChemicalPotential(elem.muQ);
      m_THM->SetStrangenessChemicalPotential(elem.muS);

      if (m_Config.CFOParameters.gammaq != 1.0)
        m_THM->SetGammaq(m_Config.CFOParameters.gammaq);

      if (m_Config.CFOParameters.gammaS != 1.0)
        m_THM->SetGammaS(m_Config.CFOParameters.gammaS);

      if (m_Config.CFOParameters.gammaC != 1.0)
        m_THM->SetGammaC(m_Config.CFOParameters.gammaC);

      double dVeff = 0.;
      for (int mu = 0; mu < 4; ++mu)
        dVeff += elem.dsigma[mu] * elem.u[mu];

      if (dVeff <= 0.) {
        continue;
      }

      if (abs(elem.edens - m_edens) <= 1.e-3)
        rescaleTmu_EMatch += dVeff * elem.edens;
//      else
//        cout << "Energy density mismatch: " << elem.edens << " vs " << m_edens << endl;
      rescaleTmu_Etot += dVeff * elem.edens;

      Veff += dVeff;
      m_Tav += elem.T * dVeff;

      m_THM->CalculatePrimordialDensities();

      std::vector<double>* densitiesIdeal = &m_THM->Densities();
      std::vector<double> tdens;
      if (m_THM->TAG() != "ThermalModelIdeal") {
        tdens = m_THM->GetIdealGasDensities();
        densitiesIdeal = &tdens;
      }

      for (size_t ipart = 0; ipart < m_THM->TPS()->ComponentsNumber(); ++ipart) {
        allweights[ipart][ielem] = m_THM->Densities()[ipart] * dVeff;

        m_Musav[ipart] += m_THM->ChemicalPotential(ipart) * dVeff;
        FullDensities[ipart] += m_THM->Densities()[ipart] * dVeff;
        FullDensitiesIdeal[ipart] += densitiesIdeal->operator[](ipart) * dVeff;

        //Npart += m_THM->TPS()->Particle(ipart).BaryonCharge() * m_THM->Densities()[ipart] * dVeff;
      }
    }

    // Free memory just in case
    std::vector<RandomGenerators::VolumeElementSampler>().swap(m_VolumeElementSamplers);
    m_VolumeElementSamplers.clear();
    for (size_t ipart = 0; ipart < m_THM->TPS()->ComponentsNumber(); ++ipart) {
      m_VolumeElementSamplers.push_back(RandomGenerators::VolumeElementSampler(allweights[ipart]));
      // Free memory
      vector<double>().swap(allweights[ipart]);
    }

    m_FullSpaceYields = FullDensities;

    m_Tav /= Veff;
    for (size_t ipart = 0; ipart < m_THM->TPS()->ComponentsNumber(); ++ipart) {
      m_Musav[ipart] /= Veff;
      FullDensities[ipart] /= Veff;
      FullDensitiesIdeal[ipart] /= Veff;
    }

    cout << endl;
    cout << "V     = " << Veff << endl;
    cout << "<T>   = " << m_Tav << endl;
    cout << "<muB> = " << m_Musav[m_THM->TPS()->PdgToId(2112)] << endl;

    // Check for energy density mismatch for T-Mu rescaling
    if (m_RescaleTmu){
      if (rescaleTmu_EMatch/rescaleTmu_Etot < 0.90) {
        printf("**WARNING** Energy density mismatch for T-Mu rescaling: %lf\n", rescaleTmu_EMatch/rescaleTmu_Etot);
      }
      cout << "Edens match fraction = " << rescaleTmu_EMatch/rescaleTmu_Etot << endl;
    }

    // The canonical ensemble
    // Make the canonical ensemble integers
    // Rescale all the means to get integer baryon number
    // Then round the other conserved charges to the nearest integer
    if (m_Config.fEnsemble == EventGeneratorConfiguration::CE) {
      double totB = 0., totQ = 0., totS = 0., totC = 0.;
      for (size_t ipart = 0; ipart < m_THM->TPS()->ComponentsNumber(); ++ipart) {
        totB += m_THM->TPS()->Particle(ipart).BaryonCharge() * m_FullSpaceYields[ipart];
        totQ += m_THM->TPS()->Particle(ipart).ElectricCharge() * m_FullSpaceYields[ipart];
        totS += m_THM->TPS()->Particle(ipart).Strangeness() * m_FullSpaceYields[ipart];
        totC += m_THM->TPS()->Particle(ipart).Charm() * m_FullSpaceYields[ipart];
      }

      cout << endl;
      cout << "Integrated values of conserved charges:" << endl;
      cout << "B = " << totB << endl;
      cout << "Q = " << totQ << endl;
      cout << "S = " << totS << endl;
      cout << "C = " << totC << endl;

      double Vcorr = 1.;
      if (m_Config.CanonicalB) {

        if (abs(totB) > 1.e-6) {
          Vcorr = round(totB) / totB;
          Veff *= Vcorr;
          cout << "Volume rescaling factor Vcorr = " << Vcorr << endl;
          cout << "B: " << totB << " -> " << round(totB) << endl;
        }

        m_Config.B = round(totB);
      }

      if (m_Config.CanonicalQ) {
        cout << "Q: " << totQ * Vcorr << " -> " << round(totQ * Vcorr) << endl;
        m_Config.Q = round(totQ * Vcorr);
      }

      if (m_Config.CanonicalS && m_Config.S != 0) {
        cout << "S: " << totS * Vcorr << " -> " << round(totS * Vcorr) << endl;
        m_Config.S = round(totS * Vcorr);
      }

      if (m_Config.CanonicalC && m_Config.C != 0) {
        cout << "C: " << totC * Vcorr << " -> " << round(totC * Vcorr) << endl;
        m_Config.C = round(totC * Vcorr);
      }
    }

    cout.flush();


    m_THM->SetVolume(Veff);
    m_THM->SetCanonicalVolume(Veff);
    m_THM->Densities() = FullDensities;
    m_DensitiesIdeal = FullDensitiesIdeal;
    if (m_Config.fEnsemble != EventGeneratorConfiguration::GCE)
      PrepareMultinomials();
  }

  std::vector<std::vector<double>> HypersurfaceEventGenerator::CalculateTMuMap(ThermalModelBase* model, double edens, double rhomin, double rhomax, double drho)
  {
    cout << "Remapping T and mu along e = " << edens << " surface..." << endl;
    vector<double> rhos, Ts, muBs, muSs, muQs, Ps;
    vector<double> Tmusini = { 0.150, 0.100, 0.033, -0.05 };
    for (double rho = rhomin; rho < rhomax + 1.e-9; rho += drho) {
      rhos.push_back(rho);

      if (rho == rhomin || rho - drho == rhomin) {
        model->SetTemperature(Tmusini[0]);
        model->SetBaryonChemicalPotential(Tmusini[1]);
        model->SetStrangenessChemicalPotential(Tmusini[2]);
        model->SetElectricChemicalPotential(Tmusini[3]);
      }

      auto Tmus = MatchEnergyBaryonDensities(model, edens, rho);
      Ts.push_back(Tmus[0]);
      muBs.push_back(Tmus[1]);
      muSs.push_back(Tmus[2]);
      muQs.push_back(Tmus[3]);
      Ps.push_back(Tmus[4]);
    }

    return { rhos, Ts, muBs, muSs, muQs, Ps };
  }


  void HypersurfaceEventGenerator::SetRescaleTmu(bool rescale, double edens)
  {
    m_RescaleTmu = rescale;
    m_edens = edens;
  }

  SimpleEvent HypersurfaceEventGenerator::GetEvent(bool PerformDecays) const
  {
    const_cast<HypersurfaceEventGenerator*>(this)->CheckSetParameters();
    return EventGeneratorBase::GetEvent(PerformDecays);
  }

  void HypersurfaceEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);

    if (m_THM != NULL) {
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        const ThermalParticle& part = m_THM->TPS()->Particles()[i];
        m_MomentumGens.push_back(new RandomGenerators::HypersurfaceMomentumGenerator(
          m_ParticlizationHypersurface,
          &m_THM->TPS()->Particle(i),
          &m_VolumeElementSamplers[i],
          m_MomentumGeneratorConfig
        ));

        // Should not be used
        double T = m_Tav;
        double Mu = m_Musav[i];
        if (m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW || m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBWconstBR)
          m_BWGens.push_back(new RandomGenerators::ThermalEnergyBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        else
          m_BWGens.push_back(new RandomGenerators::ThermalBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
      }
    }
  }

  void HypersurfaceEventGenerator::RescaleHypersurfaceParametersEdens(ParticlizationHypersurface *hypersurface,
                                                                      ThermalModelBase *model, double edens,
                                                                      double rhomin, double rhomax, double drho,
                                                                      double rhocrit) {
    auto TMuMap = CalculateTMuMap(model, edens, rhomin, rhomax, drho);
    vector<SplineFunction> SplinesTMu(5);
    for(int i = 0; i < 5; ++i)
      SplinesTMu[i].fill(TMuMap[0], TMuMap[1+i]);

    for(auto& elem : *hypersurface) {
      if (abs(elem.edens - edens) > 1.e-3 || elem.rhoB < 0.0 || elem.rhoB > 0.25) {
        // Do nothing
      }
      else {
        elem.T = SplinesTMu[0].f(elem.rhoB);
        elem.muB = SplinesTMu[1].f(elem.rhoB);
        elem.muS = SplinesTMu[2].f(elem.rhoB);
        elem.muQ = SplinesTMu[3].f(elem.rhoB);
        elem.P = SplinesTMu[4].f(elem.rhoB);
      }
    }
  }


  RandomGenerators::BoostInvariantHypersurfaceMomentumGenerator::BoostInvariantHypersurfaceMomentumGenerator(
    const ParticlizationHypersurface* hypersurface,
    VolumeElementSampler* positionsampler, double Tkin, double etamax, double mass, int statistics, double mu) :
    m_ParticlizationHypersurface(hypersurface),
    m_VolumeElementSampler(positionsampler),
    m_Tkin(Tkin), m_EtaMax(etamax), m_Mass(mass),
    m_Generator(mass, statistics, Tkin, mu)
  {

  }

  std::vector<double> RandomGenerators::BoostInvariantHypersurfaceMomentumGenerator::GetMomentum(double mass) const
  {
    if (m_VolumeElementSampler == NULL || m_ParticlizationHypersurface == NULL) {
      throw std::runtime_error("RandomGenerators::BoostInvariantHypersurfaceMomentumGenerator::GetMomentum(double mass): Hypersurface not initialized!");
    }

    if (mass < 0.)
      mass = Mass();

    std::vector<double> ret(3, 0.);
    int VolumeElementIndex = m_VolumeElementSampler->SampleVolumeElement();

    const ParticlizationHypersurfaceElement& elem = (*m_ParticlizationHypersurface)[VolumeElementIndex];

    double etaF = 0.5 * log((elem.u[0] + elem.u[3]) / (elem.u[0] - elem.u[3]));


    double vx = elem.u[1] / elem.u[0];// / cosheta;
    double vy = elem.u[2] / elem.u[0];// / cosheta;
    double vz = tanh(etaF);// tanh(eta);

    SimpleParticle part(0., 0., 0., mass, 0);

    // dsigma^\mu in the local rest frame
    std::vector<double> dsigma_loc =
      LorentzBoost(
        { elem.dsigma[0], -elem.dsigma[1], -elem.dsigma[2], -elem.dsigma[3] },
        vx, vy, vz);

    double maxWeight = 1. + std::abs(dsigma_loc[1] / dsigma_loc[0]) + std::abs(dsigma_loc[2] / dsigma_loc[0]) + std::abs(dsigma_loc[3] / dsigma_loc[0]);

    while (1) {
      double tp = m_Generator.GetP(mass);
      double tphi = 2. * xMath::Pi() * RandomGenerators::randgenMT.rand();
      double cthe = 2. * RandomGenerators::randgenMT.rand() - 1.;
      double sthe = sqrt(1. - cthe * cthe);
      part.px = tp * cos(tphi) * sthe;
      part.py = tp * sin(tphi) * sthe;
      part.pz = tp * cthe;
      part.p0 = sqrt(mass * mass + tp * tp);

      double p0LRF = part.p0;

      double dsigmamu_pmu_loc = dsigma_loc[0] * part.p0
        - dsigma_loc[1] * part.px - dsigma_loc[2] * part.py - dsigma_loc[3] * part.pz;


      double dsigmamu_umu_loc = dsigma_loc[0];

      double dumu_pmu_loc = p0LRF;

      double Weight = dsigmamu_pmu_loc / dsigmamu_umu_loc / dumu_pmu_loc / maxWeight;

      if (Weight > 1.) {
        printf("**WARNING** BoostInvariantHypersurfaceMomentumGenerator::GetMomentum: Weight exceeds unity by %E\n",
          Weight - 1.);
      }

      if (RandomGenerators::randgenMT.rand() < Weight)
        break;
    }

    if (vx != 0.0 || vy != 0.0 || vz != 0.0)
      part = ParticleDecaysMC::LorentzBoostMomentumOnly(part, -vx, -vy, -vz);

    double eta = elem.eta;
    double cosheta = cosh(eta);
    double sinheta = sinh(eta);
    // Smearing in eta
    if (EtaMax() > 0.0) {
      double deta = -EtaMax() + 2. * EtaMax() * RandomGenerators::randgenMT.rand();

      double tpz = part.GetMt() * std::sinh(part.GetY() + deta);
      double tp0 = sqrt(part.m * part.m + part.px * part.px + part.py * part.py + tpz * tpz);

      part.pz = tpz;
      part.p0 = tp0;

      eta += deta;
      cosheta = cosh(eta);
      sinheta = sinh(eta);

      //vz = tanh(eta);
      //if (vz != 0.0)
      //  part = ParticleDecaysMC::LorentzBoost(part, 0., 0., -vz);
    }


    ret[0] = part.px;
    ret[1] = part.py;
    ret[2] = part.pz;

    // Space-time coordinates
    double tau = elem.tau;
    double r0 = tau * cosheta;
    double rz = tau * sinheta;

    double rx = elem.x;
    double ry = elem.y;



    ret.push_back(r0);
    ret.push_back(rx);
    ret.push_back(ry);
    ret.push_back(rz);

    return ret;
  }

  void fillBoostMatrix(double vx, double vy, double vz, double boostMatrix[4][4])
  // Lorentz boost matrix
  // here in boostMatrix [0]=t, [1]=x, [2]=y, [3]=z
  {
    const double vv [3] = {vx, vy, vz} ;
    const double v2 = vx*vx+vy*vy+vz*vz ;
    const double gamma = 1.0/sqrt(1.0-v2) ;
    if(std::isinf(gamma)||std::isnan(gamma)){ throw std::runtime_error("boost vector invalid"); }
    boostMatrix[0][0] = gamma ;
    boostMatrix[0][1] = boostMatrix[1][0] = vx*gamma ;
    boostMatrix[0][2] = boostMatrix[2][0] = vy*gamma ;
    boostMatrix[0][3] = boostMatrix[3][0] = vz*gamma ;
    if(v2>0.0){
    for(int i=1; i<4; i++)
    for(int j=1; j<4; j++)
    boostMatrix[i][j] = (gamma-1.0)*vv[i-1]*vv[j-1]/v2 ;
    }else{
    for(int i=1; i<4; i++)
    for(int j=1; j<4; j++)
    boostMatrix[i][j] = 0.0 ;
    }
    for(int i=1; i<4; i++) boostMatrix[i][i] += 1.0 ;
  }

  int index44(const int &i, const int &j){
    // index44: returns an index of pi^{mu nu} mu,nu component in a plain 1D array
    if(i>3 || j>3 || i<0 || j<0) {std::cout<<"index44: i j " <<i<<" "<<j<<endl ; exit(1) ; }
    if(j<i) return (i*(i+1))/2 + j ;
    else return (j*(j+1))/2 + i ;
  }

  std::vector<double> RandomGenerators::HypersurfaceMomentumGenerator::SamplePhaseSpaceCoordinateFromElement(const ParticlizationHypersurfaceElement* elem, const ThermalParticle* particle, const double& mass, const double& etasmear, const bool shear_correction, const bool bulk_correction, const double speed_of_sound_squared) 
  {
    if (particle == NULL) {
      throw std::runtime_error("HypersurfaceMomentumGenerator::SamplePhaseSpaceCoordinateFromElement(): Unknown particle species!");
    }

    double etaF = 0.5 * log((elem->u[0] + elem->u[3]) / (elem->u[0] - elem->u[3]));

    double vx = elem->u[1] / elem->u[0];
    double vy = elem->u[2] / elem->u[0];
    double vz = tanh(etaF);

    SimpleParticle part(0., 0., 0., mass, 0);

    // dsigma^\mu in the local rest frame
    std::vector<double> dsigma_loc =
      LorentzBoost(
        { elem->dsigma[0], -elem->dsigma[1], -elem->dsigma[2], -elem->dsigma[3] },
        vx, vy, vz);

    // Maximum weight for the rejection sampling of the momentum
    double maxWeight = 1. + std::abs(dsigma_loc[1] / dsigma_loc[0]) + std::abs(dsigma_loc[2] / dsigma_loc[0]) + std::abs(dsigma_loc[3] / dsigma_loc[0]);

    // Regulating linear shear corrections
    double Wvisc_min = 0.1, Wvisc_max = 10.0;  // Lower bound consistent with smash-hadron-sampler, upper bound large enough so that it is not reached in practice
    // Wvisc_min = 0.; Wvisc_max = 2.0; // |delta f| < feq, following https://arxiv.org/pdf/1912.08271.pdf

    double T = elem->T;
    double mu = particle->BaryonCharge() * elem->muB + particle->ElectricCharge() * elem->muQ + particle->Strangeness() * elem->muS;
    ThermalMomentumGenerator Generator(mass, particle->Statistics(), T, mu);

    // Shear correction based on smash-hadron-sampler (ideal gas type), see https://github.com/smash-transport/smash-hadron-sampler/blob/main/src/gen.cpp#L297
    const double gmumu[4] = {1., -1., -1., -1.};
    const double C_Feq = pow(0.5*xMath::GeVtoifm() / xMath::Pi(), 3);
    int stat = -particle->Statistics(); // Opposite convention to smash-hadron-sampler for quantum statistics
    // Ideal gas distribution function f0 = feq, excluded volume not included
    // Note that feq plays no role for Maxwell-Boltzmann statistics (stat = 0)
    const double feq = C_Feq / (exp((part.p0-mu)/T) - stat); // TODO: The possible excluded volume factor not included, effect likely small and disappears in the absence of quantum statistics
    double pi_lrf[10];
    double boostMatrix[4][4];
    if (shear_correction){
      maxWeight *= Wvisc_max;
      // boost pi^{mu nu} into the local rest frame
      fillBoostMatrix(-vx, -vy, -vz, boostMatrix);
      for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
          pi_lrf[index44(i,j)] = 0.0;
          for (int k=0; k<4; k++){
            for (int l=0; l<4; l++){
              pi_lrf[index44(i,j)] += elem->pi[index44(k,l)] * boostMatrix[i][k] * boostMatrix[j][l];
            }
          }
        }
      }
    }
    
    while (true) {
      double tp = Generator.GetP(mass);
      double tphi = 2. * xMath::Pi() * RandomGenerators::randgenMT.rand();
      double cthe = 2. * RandomGenerators::randgenMT.rand() - 1.;
      double sthe = sqrt(1. - cthe * cthe);
      part.px = tp * cos(tphi) * sthe;
      part.py = tp * sin(tphi) * sthe;
      part.pz = tp * cthe;
      part.p0 = sqrt(mass * mass + tp * tp);

      double p0LRF = part.p0;

      double dsigmamu_pmu_loc = dsigma_loc[0] * part.p0
        - dsigma_loc[1] * part.px - dsigma_loc[2] * part.py - dsigma_loc[3] * part.pz;

      double dsigmamu_umu_loc = dsigma_loc[0];

      double dumu_pmu_loc = p0LRF;

      double Weight = dsigmamu_pmu_loc / dsigmamu_umu_loc / dumu_pmu_loc / maxWeight;

      double Weight_visc = 1.0;
      if (shear_correction){
        double mom[4] = {part.p0, part.px, part.py, part.pz};
        double pipp = 0.0;
        for (int i=0; i<4; i++){
          for (int j=0; j<4; j++){
            pipp += mom[i] * mom[j] * gmumu[i] * gmumu[j] * pi_lrf[index44(i,j)];
          }
        }
        // this is in principle the ansatz which is also used in https://github.com/smash-transport/smash-hadron-sampler
        // from this paper Phys.Rev.C 73 (2006) 064903
        Weight_visc += ((1.0 + stat * feq) * pipp / (2.0 * T * T * (elem->edens + elem->P)));
        
      }
      if (bulk_correction){
         // this is in principle the ansatz which is also used in https://github.com/smash-transport/smash-hadron-sampler
         // Eq. (7) of https://arxiv.org/pdf/1509.06738 plus Eq. (4) from https://arxiv.org/pdf/1403.0962
         // see also https://github.com/smash-transport/smash-hadron-sampler/files/14011233/bulk_corrections_note.pdf
         double mom[4] = {part.p0, part.px, part.py, part.pz};
         Weight_visc -= (1.0+stat*feq)*elem->Pi
                 *(mass*mass/(3*mom[0])-mom[0]*(1.0/3.0-speed_of_sound_squared))
                 /(15*(1.0/3.0-speed_of_sound_squared)*(1.0/3.0-speed_of_sound_squared)*T*(elem->edens + elem->P))  ;
      }
      if (bulk_correction || shear_correction){
        if (Weight_visc<Wvisc_min) Weight_visc = Wvisc_min;
        if (Weight_visc>Wvisc_max) Weight_visc = Wvisc_max;
        // update weight with viscosity factor
        Weight *= Weight_visc;
      }
      
      

      if (Weight > 1.) {
        printf("**WARNING** HypersurfaceSampler::GetMomentum: Weight exceeds unity by %E\n",
          Weight - 1.);
      }
      if (RandomGenerators::randgenMT.rand() < Weight)
        break;
    }

    if (vx != 0.0 || vy != 0.0 || vz != 0.0)
      part = ParticleDecaysMC::LorentzBoostMomentumOnly(part, -vx, -vy, -vz);

    double eta = elem->eta;
    double cosheta = cosh(eta);
    double sinheta = sinh(eta);

    // Smearing in eta
    if (etasmear > 0.0) {
      double deta = -0.5 * etasmear + 1. * etasmear * RandomGenerators::randgenMT.rand();

      double tpz = part.GetMt() * std::sinh(part.GetY() + deta);
      double tp0 = sqrt(part.m * part.m + part.px * part.px + part.py * part.py + tpz * tpz);

      part.pz = tpz;
      part.p0 = tp0;

      eta += deta;
      cosheta = cosh(eta);
      sinheta = sinh(eta);
    }

    // Space-time coordinates
    double tau = elem->tau;
    double r0 = tau * cosheta;
    double rz = tau * sinheta;

    double rx = elem->x;
    double ry = elem->y;

    return { part.px, part.py, part.pz, r0, rx, ry, rz };
  }

  BoostInvariantHypersurfaceEventGenerator::BoostInvariantHypersurfaceEventGenerator(ThermalParticleSystem* TPS, const EventGeneratorConfiguration& config, double etamax, const ParticlizationHypersurface* hypersurface)
    : m_VolumeElementSampler(NULL)
  {
    SetConfiguration(TPS, config);

    SetParameters(etamax, hypersurface);
  }

  BoostInvariantHypersurfaceEventGenerator::~BoostInvariantHypersurfaceEventGenerator()
  {
    if (m_VolumeElementSampler != NULL)
      delete m_VolumeElementSampler;
  }

  void BoostInvariantHypersurfaceEventGenerator::SetParameters(double etamax, const ParticlizationHypersurface* hypersurface)
  {
    m_EtaMax = etamax;
    m_ParticlizationHypersurface = hypersurface;
    SetMomentumGenerators();
  }

  void BoostInvariantHypersurfaceEventGenerator::SetMomentumGenerators()
  {
    ClearMomentumGenerators();
    m_BWGens.resize(0);

    if (m_VolumeElementSampler != NULL)
      delete m_VolumeElementSampler;

    m_VolumeElementSampler =
      new RandomGenerators::VolumeElementSampler(m_ParticlizationHypersurface);


    if (m_THM != NULL) {
      for (size_t i = 0; i < m_THM->TPS()->Particles().size(); ++i) {
        const ThermalParticle& part = m_THM->TPS()->Particles()[i];
        m_MomentumGens.push_back(new RandomGenerators::BoostInvariantHypersurfaceMomentumGenerator(
          m_ParticlizationHypersurface,
          m_VolumeElementSampler,
          m_Config.CFOParameters.T,
          GetEtaMax(),
          part.Mass(),
          part.Statistics(),
          m_THM->FullIdealChemicalPotential(i)
        ));

        double T = m_THM->Parameters().T;
        double Mu = m_THM->FullIdealChemicalPotential(i);
        if (m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW || m_THM->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBWconstBR)
          m_BWGens.push_back(new RandomGenerators::ThermalEnergyBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
        else
          m_BWGens.push_back(new RandomGenerators::ThermalBreitWignerGenerator(&m_THM->TPS()->Particle(i), T, Mu));
      }
    }
  }

  SimpleEvent HypersurfaceEventGeneratorEVHRG::GetEvent(bool DoDecays) const
  {
    const_cast<HypersurfaceEventGeneratorEVHRG*>(this)->CheckSetParameters();

    std::pair< std::vector<int>, double > yieldsW = SampleYields();
    std::pair<int, int> NBBbar = ComputeNBNBbar(yieldsW.first);
    double wgtB    = EVHRGWeight( NBBbar.first,  m_MeanB, m_VEff, m_b);
    double wgtBbar = EVHRGWeight(NBBbar.second, m_MeanAB, m_VEff, m_b);
    while (RandomGenerators::randgenMT.rand() > wgtB || RandomGenerators::randgenMT.rand() > wgtBbar) {
      yieldsW = SampleYields();
      NBBbar = ComputeNBNBbar(yieldsW.first);
      wgtB = EVHRGWeight(NBBbar.first, m_MeanB, m_VEff, m_b);
      wgtBbar = EVHRGWeight(NBBbar.second, m_MeanAB, m_VEff, m_b);
    }

    std::vector<int>& yields = yieldsW.first;

    SimpleEvent ret = SampleParticles(yields);

    if (DoDecays)
      return PerformDecays(ret, m_THM->TPS());
    else
      return ret;
  }

  double HypersurfaceEventGeneratorEVHRG::EVHRGWeight(int sampledN, double meanN, double V, double b)
  {
    double nev = meanN / V;
    return pow((1. - b * sampledN / V) / (1. - b * nev), sampledN)
      * exp(b * nev / (1. - b * nev) * (sampledN - meanN));
  }

  void HypersurfaceEventGeneratorEVHRG::SetParameters()
  {
    HypersurfaceEventGenerator::SetParameters();
    m_MeanB = 0.0; m_MeanAB = 0.0;
    m_VEff = m_THM->Volume();

    for (int ipart = 0; ipart < m_THM->TPS()->ComponentsNumber(); ++ipart) {
      const ThermalParticle& part = m_THM->TPS()->Particle(ipart);
      if (part.BaryonCharge() == 1)
        m_MeanB += FullSpaceYields()[ipart];
      if (part.BaryonCharge() == -1)
        m_MeanAB += FullSpaceYields()[ipart];
    }

    cout << "nB  = " << m_MeanB  / m_VEff << endl;
    cout << "naB = " << m_MeanAB / m_VEff << endl;
  }

  SimpleEvent HypersurfaceEventGeneratorEVHRG::SampleParticles(const std::vector<int>& yields) const
  {
    SimpleEvent ret;

    // Hard-core radius for (anti)baryons in the sampling procedure
    double radB = m_rad;
    if (radB < 0.0)
      radB = CuteHRGHelper::rv(m_b);

    // Reshuffle the order of particles to be sampled
    std::vector<int> idsM, idsB, idsaB;
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
      if (m_THM->TPS()->Particles()[i].BaryonCharge() != 1 && m_THM->TPS()->Particles()[i].BaryonCharge() != -1)
        for (int part = 0; part < yields[i]; ++part)
          idsM.push_back(i);
    //std::random_shuffle(idsM.begin(), idsM.end()); // Removed in C++17
    std::shuffle(idsM.begin(), idsM.end(), RandomGenerators::rng_std);
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
      if (m_THM->TPS()->Particles()[i].BaryonCharge() == 1)
        for (int part = 0; part < yields[i]; ++part)
          idsB.push_back(i);
    //std::random_shuffle(idsB.begin(), idsB.end()); // Removed in C++17
    std::shuffle(idsB.begin(), idsB.end(), RandomGenerators::rng_std);
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
      if (m_THM->TPS()->Particles()[i].BaryonCharge() == -1)
        for (int part = 0; part < yields[i]; ++part)
          idsaB.push_back(i);
    //std::random_shuffle(idsaB.begin(), idsaB.end()); // Removed in C++17
    std::shuffle(idsaB.begin(), idsaB.end(), RandomGenerators::rng_std);

    std::vector<int> ids;
    ids.insert(ids.end(), idsM.begin(), idsM.end());
    ids.insert(ids.end(), idsB.begin(), idsB.end());
    ids.insert(ids.end(), idsaB.begin(), idsaB.end());

    int idBstart = idsM.size();
    int idaBstart = idsM.size() + idsB.size();

    ret.Particles.resize(ids.size());
    for (int ip = 0; ip < ids.size(); ++ip) {
      int pid = ids[ip];
      const ThermalParticle& species = m_THM->TPS()->Particles()[pid];
      SimpleParticle part = SampleParticle(pid);

      if (radB > 0.0) {
        if (species.BaryonCharge() == 1 || species.BaryonCharge() == -1) {
          bool flOverlap = false;
          int tip = ip - 1;
          while (tip >= 0 && m_THM->TPS()->Particles()[ids[tip]].BaryonCharge() == species.BaryonCharge()) {
            double dist2 = ParticleDecaysMC::ParticleDistanceSquared(ret.Particles[tip], part);
            flOverlap |= (dist2 <= 4. * radB * radB);
            tip--;
          }

          if (flOverlap) {
            if (EVUseSPR()) {
              ip--;
            }
            else {
              if (species.BaryonCharge() == 1)
                ip = idBstart - 1;
              else
                ip = idaBstart - 1;
            }
            continue;
          }
        }
      }

      ret.Particles[ip] = part;
    }

    ret.AllParticles = ret.Particles;

    ret.DecayMap.resize(ret.Particles.size());
    fill(ret.DecayMap.begin(), ret.DecayMap.end(), -1);

    ret.DecayMapFinal.resize(ret.Particles.size());
    for (int i = 0; i < ret.DecayMapFinal.size(); ++i)
      ret.DecayMapFinal[i] = i;

    return ret;
  }

  std::pair<int, int> HypersurfaceEventGeneratorEVHRG::ComputeNBNBbar(const std::vector<int>& yields) const
  {
    int NB = 0, NBbar = 0;
    for (int ipart = 0; ipart < m_THM->TPS()->ComponentsNumber(); ++ipart) {
      const ThermalParticle& part = m_THM->TPS()->Particle(ipart);
      if (part.BaryonCharge() == 1)
        NB += yields[ipart];
      if (part.BaryonCharge() == -1)
        NBbar += yields[ipart];
    }
    return std::make_pair(NB, NBbar);
  }

} // namespace thermalfist
