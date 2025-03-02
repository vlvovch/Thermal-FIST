/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2024 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ExtraParticles.h"

#include <algorithm>
#include <cstdlib>

#include "HRGBase/Utility.h"

using namespace std;

namespace thermalfist {

  namespace ExtraParticles {
    static std::vector<ThermalParticle> Particles;
    static std::map<long long, int> PdgIdMap;
    static bool isInitialized = []() {
      return Init();
    }();

    const ThermalParticle& Particle(int id) {
      if (id < 0 || id >= Particles.size()) {
        throw std::out_of_range("ExtraParticles::Particle: id is out of bounds!");
      }
      return Particles[id];
    }

    const ThermalParticle& ParticleByPdg(long long pdgid)
    {
      int tid = PdgToId(pdgid);
      if (tid == -1) {
        throw std::invalid_argument("ExtraParticles::ParticleByPdg: pdgid " + std::to_string(pdgid) + " is unknown");
      }
      return Particle(tid);
    }
    int PdgToId(long long pdgid)
    {
      return (PdgIdMap.count(pdgid) > 0) ? PdgIdMap[pdgid] : -1;
    }

    // Initializes the particles and their PDG ID mappings.
    // This function is used for static initialization.
    bool Init()
    {
      Particles.clear();
      PdgIdMap.clear();
      
      int tsz = 0;
      // photons
      Particles.push_back(ThermalParticle(true, "gamma", 22, 2., -1, 0.));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;

      // electrons
      Particles.push_back(ThermalParticle(true, "e-", 11, 2., 1, 5.109989461E-04, 0, 0, -1));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;
      Particles.push_back(ThermalParticle(true, "e+", -11, 2., 1, 5.109989461E-04, 0, 0, 1));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;

      // muons
      Particles.push_back(ThermalParticle(true, "mu-", 13, 2., 1, 1.056583745E-01, 0, 0, -1));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;
      Particles.push_back(ThermalParticle(true, "mu+", -13, 2., 1, 1.056583745E-01, 0, 0, 1));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;

      // tauons
      Particles.push_back(ThermalParticle(true, "tau-", 15, 2., 1, 1.77686E+00, 0, 0, -1));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;
      Particles.push_back(ThermalParticle(true, "tau+", -15, 2., 1, 1.77686E+00, 0, 0, 1));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;

      // nu(e)
      Particles.push_back(ThermalParticle(true, "nu(e)", 12, 1., 1, 0.));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;
      Particles.push_back(ThermalParticle(true, "anti-nu(e)", -12, 1., 1, 0.));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;

      // nu(mu)
      Particles.push_back(ThermalParticle(true, "nu(mu)", 14, 1., 1, 0.));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;
      Particles.push_back(ThermalParticle(true, "anti-nu(mu)", -14, 1., 1, 0.));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;

      // nu(tau)
      Particles.push_back(ThermalParticle(true, "nu(tau)", 16, 1., 1, 0.));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;
      Particles.push_back(ThermalParticle(true, "anti-nu(tau)", -16, 1., 1, 0.));
      PdgIdMap[Particles[tsz].PdgId()] = tsz;
      tsz++;

      return true;
    }

    std::string NameByPdg(long long pdg)
    {
      int tid = PdgToId(pdg);
      if (tid != -1)
        return Particle(tid).Name();
      return std::string("???");
    }
  } // namespace ExtraParticles

} // namespace thermalfist