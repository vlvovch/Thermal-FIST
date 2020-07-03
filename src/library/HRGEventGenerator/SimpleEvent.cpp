/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/SimpleEvent.h"

#include <iomanip>

namespace thermalfist {

  void SimpleEvent::writeToFile(std::ofstream & fout, int eventnumber)
  {
    fout << "Event " << eventnumber << std::endl;
    fout << "Weight: " << weight << std::endl;

    fout << std::setw(20) << "pdgid"
      << std::setw(20) << "px[GeV]"
      << std::setw(20) << "py[GeV]"
      << std::setw(20) << "pz[GeV]"
      << std::setw(20) << "motherpdgid"
      << std::endl;

    for (size_t i = 0; i < Particles.size(); ++i) {
      fout << std::setw(20) << Particles[i].PDGID
        << std::setw(20) << Particles[i].px
        << std::setw(20) << Particles[i].py
        << std::setw(20) << Particles[i].pz
        << std::setw(20) << Particles[i].MotherPDGID
        << std::endl;
    }
    fout << std::endl;
  }

  void SimpleEvent::RapidityBoost(double dY)
  {
    for (size_t i = 0; i < Particles.size(); ++i)
      Particles[i].RapidityBoost(dY);
  }

  SimpleEvent SimpleEvent::MergeEvents(const SimpleEvent & evt1, const SimpleEvent & evt2)
  {
    SimpleEvent ret;
    ret.Particles.reserve(evt1.Particles.size() + evt2.Particles.size());
    ret.Particles.insert(ret.Particles.end(), evt1.Particles.begin(), evt1.Particles.end());
    ret.Particles.insert(ret.Particles.end(), evt2.Particles.begin(), evt2.Particles.end());

    ret.AllParticles.reserve(evt1.AllParticles.size() + evt2.AllParticles.size());
    ret.AllParticles.insert(ret.AllParticles.end(), evt1.AllParticles.begin(), evt1.AllParticles.end());
    ret.AllParticles.insert(ret.AllParticles.end(), evt2.AllParticles.begin(), evt2.AllParticles.end());

    ret.DecayMap.reserve(evt1.DecayMap.size() + evt2.DecayMap.size());
    ret.DecayMap.insert(ret.DecayMap.end(), evt1.DecayMap.begin(), evt1.DecayMap.end());
    ret.DecayMap.insert(ret.DecayMap.end(), evt2.DecayMap.begin(), evt2.DecayMap.end());

    ret.DecayMapFinal.reserve(evt1.DecayMapFinal.size() + evt2.DecayMapFinal.size());
    ret.DecayMapFinal.insert(ret.DecayMapFinal.end(), evt1.DecayMapFinal.begin(), evt1.DecayMapFinal.end());
    ret.DecayMapFinal.insert(ret.DecayMapFinal.end(), evt2.DecayMapFinal.begin(), evt2.DecayMapFinal.end());
    
    // TODO: check if proper to combine weights like that
    ret.weight = evt1.weight * evt2.weight;
    ret.logweight = evt1.logweight + evt2.logweight;

    return ret;
  }

} // namespace thermalfist
