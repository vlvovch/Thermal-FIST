/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
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

    for (int i = 0; i < Particles.size(); ++i) {
      fout << std::setw(20) << Particles[i].PDGID
        << std::setw(20) << Particles[i].px
        << std::setw(20) << Particles[i].py
        << std::setw(20) << Particles[i].pz
        << std::setw(20) << Particles[i].MotherPDGID
        << std::endl;
    }
    fout << std::endl;
  }

} // namespace thermalfist
