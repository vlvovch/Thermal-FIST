/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2024 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EXTRAPARTICLES_H
#define EXTRAPARTICLES_H

#include <map>
#include <vector>
#include <set>
#include <fstream>

#include "HRGBase/ThermalParticle.h"

namespace thermalfist {

  /// Contains properties of non-QCD particles such as photons and leptons
  namespace ExtraParticles {
    const ThermalParticle& Particle(int id);
    const ThermalParticle& ParticleByPdg(long long pdg);
    int PdgToId(long long pdg);
    bool Init();
    std::string NameByPdg(long long pdg);
  }

  namespace DecayLifetimes {
    // In units of ctau (fm)
    double GetLifetime(long long pdg);
  }

} // namespace thermalfist

#endif
