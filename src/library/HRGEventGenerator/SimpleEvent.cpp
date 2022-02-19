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

  void SimpleEvent::writeToFile(std::ofstream& fout, const EventOutputConfig& config, int eventnumber) const
  {
    fout << "Event " << eventnumber << std::endl;

    if (config.printWeight)
      fout << "Weight: " << weight << std::endl;

    fout << std::setw(20) << "pdgid";

    if (config.printCoordinates)
      fout << std::setw(20) << "r0[fm/c]"
      << std::setw(20) << "rx[fm]"
      << std::setw(20) << "ry[fm]"
      << std::setw(20) << "rz[fm]";

    if (config.printEnergy)
      fout << std::setw(20) << "p0[GeV/c2]";

    fout << std::setw(20) << "px[GeV/c]"
      << std::setw(20) << "py[GeV/c]"
      << std::setw(20) << "pz[GeV/c]";

    if (config.printMotherPdg)
      fout << std::setw(20) << "mother_pdgid";

    if (config.printDecayEpoch)
      fout << std::setw(20) << "decay_epoch";

    fout << std::endl;

    fout.precision(10);
    fout << std::scientific;

    for (size_t i = 0; i < Particles.size(); ++i) {
      fout << std::setw(20) << Particles[i].PDGID;

      if (config.printCoordinates)
        fout << std::setw(20) << Particles[i].r0
        << std::setw(20) << Particles[i].rx
        << std::setw(20) << Particles[i].ry
        << std::setw(20) << Particles[i].rz;


      if (config.printEnergy)
        fout << std::setw(20) << Particles[i].p0;

      fout << std::setw(20) << Particles[i].px
        << std::setw(20) << Particles[i].py
        << std::setw(20) << Particles[i].pz;

      if (config.printMotherPdg)
        fout << std::setw(20) << Particles[i].MotherPDGID;

      if (config.printDecayEpoch)
        fout << std::setw(20) << Particles[i].epoch;

      fout << std::endl;
    }


    if (config.printPhotonsLeptons) {
      for (size_t i = 0; i < PhotonsLeptons.size(); ++i) {
        fout << std::setw(20) << PhotonsLeptons[i].PDGID;

        if (config.printCoordinates)
          fout << std::setw(20) << PhotonsLeptons[i].r0
          << std::setw(20) << PhotonsLeptons[i].rx
          << std::setw(20) << PhotonsLeptons[i].ry
          << std::setw(20) << PhotonsLeptons[i].rz;

        if (config.printEnergy)
          fout << std::setw(20) << PhotonsLeptons[i].p0;

        fout << std::setw(20) << PhotonsLeptons[i].px
          << std::setw(20) << PhotonsLeptons[i].py
          << std::setw(20) << PhotonsLeptons[i].pz;


        if (config.printMotherPdg)
          fout << std::setw(20) << PhotonsLeptons[i].MotherPDGID;

        if (config.printDecayEpoch)
          fout << std::setw(20) << PhotonsLeptons[i].epoch;

        fout << std::endl;
      }
    }
    fout << std::endl;
    fout << std::fixed;
  }

  void SimpleEvent::writeToFileForUrqmd(std::ofstream& fout) const
  {
    fout << "# " << Particles.size() << std::endl;

    fout.precision(16);
    fout << std::scientific;

    const int tabsize = 23;

    for (size_t i = 0; i < Particles.size(); ++i) {
      fout << std::setw(12) << Particles[i].PDGID << " ";

      fout << std::setw(tabsize) << Particles[i].r0 << " "
        << std::setw(tabsize) << Particles[i].rx << " "
        << std::setw(tabsize) << Particles[i].ry << " "
        << std::setw(tabsize) << Particles[i].rz << " ";


      fout << std::setw(tabsize) << Particles[i].p0 << " ";
      fout << std::setw(tabsize) << Particles[i].px << " "
        << std::setw(tabsize) << Particles[i].py << " "
        << std::setw(tabsize) << Particles[i].pz << " ";

      fout << std::endl;
    }

    fout << std::fixed;
  }

  void SimpleEvent::RapidityBoost(double dY)
  {
    for (size_t i = 0; i < Particles.size(); ++i)
      Particles[i].RapidityBoost(dY);
    for (size_t i = 0; i < AllParticles.size(); ++i)
      AllParticles[i].RapidityBoost(dY);
  }

  SimpleEvent SimpleEvent::MergeEvents(const SimpleEvent& evt1, const SimpleEvent& evt2)
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
    int offset = evt1.DecayMap.size();
    for (int i = evt1.DecayMap.size(); i < ret.DecayMap.size(); i++)
      if (ret.DecayMap[i] != -1)
        ret.DecayMap[i] += offset;

    ret.DecayMapFinal.reserve(evt1.DecayMapFinal.size() + evt2.DecayMapFinal.size());
    ret.DecayMapFinal.insert(ret.DecayMapFinal.end(), evt1.DecayMapFinal.begin(), evt1.DecayMapFinal.end());
    ret.DecayMapFinal.insert(ret.DecayMapFinal.end(), evt2.DecayMapFinal.begin(), evt2.DecayMapFinal.end());
    for (int i = evt1.DecayMapFinal.size(); i < ret.DecayMapFinal.size(); i++)
      ret.DecayMapFinal[i] += offset;

    // TODO: check if proper to combine weights like that
    ret.weight = evt1.weight * evt2.weight;
    ret.logweight = evt1.logweight + evt2.logweight;

    return ret;
  }

} // namespace thermalfist
