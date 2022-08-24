/*
 * Thermal-FIST package
 *
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

#include <iostream>

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
    double etasmear) :
    m_ParticlizationHypersurface(hypersurface),
    m_Particle(particle),
    m_VolumeElementSampler(positionsampler),
    m_EtaSmear(etasmear)
  {

  }

  std::vector<double> RandomGenerators::HypersurfaceMomentumGenerator::GetMomentum(double mass) const
  {
    if (m_VolumeElementSampler == NULL || m_ParticlizationHypersurface == NULL) {
      printf("**ERROR** in RandomGenerators::HypersurfaceMomentumGenerator::GetMomentum(double mass): Hypersurface not initialized!\n");
      return { 0., 0., 0., 0., 0., 0., 0. };
    }

    if (mass < 0.)
      mass = Mass();

    int VolumeElementIndex = m_VolumeElementSampler->SampleVolumeElement();

    const ParticlizationHypersurfaceElement& elem = (*m_ParticlizationHypersurface)[VolumeElementIndex];

    return SamplePhaseSpaceCoordinateFromElement(&elem, m_Particle, mass, EtaSmear());
  }

  HypersurfaceEventGenerator::HypersurfaceEventGenerator(ThermalParticleSystem* TPS, const EventGeneratorConfiguration& config, const ParticlizationHypersurface* hypersurface, double etasmear) :
    EventGeneratorBase()
  {
    SetConfiguration(TPS, config);
    SetHypersurface(hypersurface);
    SetEtaSmear(etasmear);
    SetRescaleTmu();
    //SetParameters(hypersurface, m_THM, etasmear);
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
      std::vector<std::vector<double>> TMuMap = CalculateTMuMap(m_THM, m_edens);
      m_SplinesTMu.resize(4);
      for(int i = 0; i < 4; ++i)
        m_SplinesTMu[i].fill(TMuMap[0], TMuMap[1+i]);
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

    //for (const auto& elem : *m_ParticlizationHypersurface) {
    for (size_t ielem = 0; ielem < m_ParticlizationHypersurface->size(); ++ielem) {
      if (ielem % 10000 == 0) {
        cout << ielem << " ";
        cout.flush();
      }
      const auto& elem = m_ParticlizationHypersurface->operator[](ielem);

      if (!m_RescaleTmu || abs(elem.edens - m_edens) > 1.e-3 || elem.rhoB < 0.0 || elem.rhoB > 0.25) {
        m_THM->SetTemperature(elem.T);
        m_THM->SetBaryonChemicalPotential(elem.muB);
        m_THM->SetElectricChemicalPotential(elem.muQ);
        m_THM->SetStrangenessChemicalPotential(elem.muS);
      }
      else {
        m_THM->SetTemperature(m_SplinesTMu[0].f(elem.rhoB));
        m_THM->SetBaryonChemicalPotential(m_SplinesTMu[1].f(elem.rhoB));
        m_THM->SetElectricChemicalPotential(m_SplinesTMu[3].f(elem.rhoB));
        m_THM->SetStrangenessChemicalPotential(m_SplinesTMu[2].f(elem.rhoB));
      }

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
    vector<double> rhos, Ts, muBs, muSs, muQs;
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
      //cout << rho << " : " << Tmus[0] << " " << Tmus[1] << " " << Tmus[2] << " " << Tmus[3] << " ";
    }

    return { rhos, Ts, muBs, muSs, muQs };
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
          GetEtaSmear()
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
      printf("**ERROR** in RandomGenerators::BoostInvariantHypersurfaceMomentumGenerator::GetMomentum(double mass): Hypersurface not initialized!\n");
      return { 0., 0., 0., 0., 0., 0., 0. };
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

  std::vector<double> RandomGenerators::HypersurfaceMomentumGenerator::SamplePhaseSpaceCoordinateFromElement(const ParticlizationHypersurfaceElement* elem, const ThermalParticle* particle, const double& mass, const double& etasmear)
  {
    if (particle == NULL) {
      printf("**ERROR** in HypersurfaceMomentumGenerator::SamplePhaseSpaceCoordinateFromElement(): Unknown particle species!\n");
      return { 0., 0., 0., 0., 0., 0., 0. };
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

    double T = elem->T;
    double mu = particle->BaryonCharge() * elem->muB + particle->ElectricCharge() * elem->muQ + particle->Strangeness() * elem->muS;
    ThermalMomentumGenerator Generator(mass, particle->Statistics(), T, mu);

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

      if (Weight > 1.) {
        printf("**WARNING** BoostInvariantHypersurfaceMomentumGenerator::GetMomentum: Weight exceeds unity by %E\n",
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
    std::random_shuffle(idsM.begin(), idsM.end());
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
      if (m_THM->TPS()->Particles()[i].BaryonCharge() == 1)
        for (int part = 0; part < yields[i]; ++part)
          idsB.push_back(i);
    std::random_shuffle(idsB.begin(), idsB.end());
    for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
      if (m_THM->TPS()->Particles()[i].BaryonCharge() == -1)
        for (int part = 0; part < yields[i]; ++part)
          idsaB.push_back(i);
    std::random_shuffle(idsaB.begin(), idsaB.end());

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


  
  //SimpleEvent HypersurfaceEventGeneratorEVHRG::SampleParticles(const std::vector<int>& yields) const
  //{
  //  SimpleEvent ret;

  //  // Hard-core radius for (anti)baryons in the sampling procedure
  //  double radB = m_rad;
  //  if (radB < 0.0)
  //      radB = CuteHRGHelper::rv(m_b);

  //  // Reshuffle the order of particles to be sampled
  //  std::vector<int> ids;
  //  for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
  //    for (int part = 0; part < yields[i]; ++part)
  //      ids.push_back(i);
  //  std::random_shuffle(ids.begin(), ids.end());

  //  // Sample the particles
  //  bool flOverlap = true;

  //  while (flOverlap) {
  //    int sampled = 0;
  //    std::vector<int> idBaryons, idAntiBaryons;
  //    ret.Particles.clear();
  //    while (sampled < ids.size()) {
  //      flOverlap = false;
  //      int i = ids[sampled];
  //      const ThermalParticle& species = m_THM->TPS()->Particles()[i];
  //      SimpleParticle part = SampleParticle(i);

  //      // Reject the (anti)baryon if it overlaps with another (anti)baryon
  //      if (radB > 0.0) {
  //        //bool flOverlap = false;

  //        if (species.BaryonCharge() == 1) {
  //          for (int ip = 0; ip < idBaryons.size(); ++ip) {
  //            double dist2 = ParticleDecaysMC::ParticleDistanceSquared(ret.Particles[idBaryons[ip]], part);
  //            flOverlap |= (dist2 <= 4. * radB * radB);
  //          }
  //        }

  //        if (species.BaryonCharge() == -1) {
  //          for (int ip = 0; ip < idAntiBaryons.size(); ++ip) {
  //            double dist2 = ParticleDecaysMC::ParticleDistanceSquared(ret.Particles[idAntiBaryons[ip]], part);
  //            flOverlap |= (dist2 <= 4. * radB * radB);
  //          }
  //        }
  //        if (flOverlap) {
  //          if (!EVFastMode())
  //            break;
  //          else
  //            continue;
  //          //printf("Reject baryon  ");
  //          //continue;
  //        }
  //      }

  //      ret.Particles.push_back(part);
  //      if (species.BaryonCharge() == 1)
  //        idBaryons.push_back(sampled);
  //      if (species.BaryonCharge() == -1)
  //        idAntiBaryons.push_back(sampled);
  //      sampled++;
  //    }
  //  }

  //  ret.AllParticles = ret.Particles;

  //  ret.DecayMap.resize(ret.Particles.size());
  //  fill(ret.DecayMap.begin(), ret.DecayMap.end(), -1);

  //  ret.DecayMapFinal.resize(ret.Particles.size());
  //  for (int i = 0; i < ret.DecayMapFinal.size(); ++i)
  //    ret.DecayMapFinal[i] = i;

  //  return ret;
  //}
  //

  //SimpleEvent HypersurfaceEventGeneratorEVHRG::SampleParticles(const std::vector<int>& yields) const
  //{
  //  SimpleEvent ret;

  //  // Hard-core radius for (anti)baryons in the sampling procedure
  //  double radB = m_rad;
  //  if (radB < 0.0)
  //    radB = CuteHRGHelper::rv(m_b);

  //  // Reshuffle the order of particles to be sampled
  //  std::vector<int> ids;
  //  for (int i = 0; i < m_THM->TPS()->Particles().size(); ++i)
  //    for (int part = 0; part < yields[i]; ++part)
  //      ids.push_back(i);
  //  std::random_shuffle(ids.begin(), ids.end());

  //  // Sample the particles

  //  // First the mesons (and light nuclei if applicable)
  //  int sampled = 0;
  //  ret.Particles.resize(ids.size());
  //  for (int ip = 0; ip < ids.size(); ++ip) {
  //    int pid = ids[ip];
  //    const ThermalParticle& species = m_THM->TPS()->Particles()[pid];
  //    if (species.BaryonCharge() != 1 && species.BaryonCharge() != -1) {
  //      ret.Particles[sampled] = SampleParticle(pid);
  //      sampled++;
  //    }
  //  }

  //  // Then, the (anti)baryons
  //  int sampledSoFar = sampled;

  //  //std::vector<int> idBaryons, idAntiBaryons;
  //  //for (int ip = 0; ip < ids.size(); ++ip) {
  //  //  int pid = ids[ip];
  //  //  const ThermalParticle& species = m_THM->TPS()->Particles()[pid];
  //  //  if (species.BaryonCharge() == 1)
  //  //    idBaryons.push_back(sampled);
  //  //  if (species.BaryonCharge() == -1)
  //  //    idAntiBaryons.push_back(sampled);

  //  //  sampled++;
  //  //}

  //  for (int idB = 0; idB < 2; idB++) {

  //    int Bcharge = 1;
  //    if (idB == 1)
  //      Bcharge = -1;

  //    sampled = sampledSoFar;

  //    std::vector<int> idBaryons;
  //    for (int ip = 0; ip < ids.size(); ++ip) {
  //      int pid = ids[ip];
  //      const ThermalParticle& species = m_THM->TPS()->Particles()[pid];
  //      if (species.BaryonCharge() == Bcharge)
  //        idBaryons.push_back(sampled);

  //      sampled++;
  //    }

  //    bool flOverlap = true;

  //    while (flOverlap) {
  //      sampled = sampledSoFar;
  //      for (int ip = 0; ip < ids.size(); ++ip) {
  //        flOverlap = false;
  //        int pid = ids[ip];
  //        const ThermalParticle& species = m_THM->TPS()->Particles()[pid];
  //        if (species.BaryonCharge() != Bcharge)
  //          continue;

  //        SimpleParticle part = SampleParticle(pid);

  //        // Reject the (anti)baryon if it overlaps with another (anti)baryon
  //        if (radB > 0.0) {

  //          for (int iip = 0; iip < idBaryons.size(); ++iip) {
  //            if (idBaryons[iip] >= sampled)
  //              break;
  //            double dist2 = ParticleDecaysMC::ParticleDistanceSquared(ret.Particles[idBaryons[iip]], part);
  //            flOverlap |= (dist2 <= 4. * radB * radB);
  //          }

  //          if (flOverlap) {
  //            //printf("Reject baryon  ");
  //            if (!EVFastMode()) {
  //              break;
  //            }
  //            else {
  //              ip--;
  //              continue;
  //            }
  //          }

  //        }

  //        ret.Particles[sampled] = part;
  //        sampled++;
  //      }
  //    }

  //    sampledSoFar = sampled;
  //  }

  //  ret.AllParticles = ret.Particles;

  //  ret.DecayMap.resize(ret.Particles.size());
  //  fill(ret.DecayMap.begin(), ret.DecayMap.end(), -1);

  //  ret.DecayMapFinal.resize(ret.Particles.size());
  //  for (int i = 0; i < ret.DecayMapFinal.size(); ++i)
  //    ret.DecayMapFinal[i] = i;

  //  return ret;
  //}

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



