/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef SIMPLEEVENT_H
#define SIMPLEEVENT_H

#include <cmath>
#include <vector>
#include <fstream>

#include "HRGEventGenerator/SimpleParticle.h"

namespace thermalfist {

  struct SimpleEvent {
    double weight, logweight;
    std::vector<SimpleParticle> Particles;
    SimpleEvent() { Particles.resize(0); weight = 1.; }
    void writeToFile(std::ofstream & fout, int eventnumber = 1);
  };

} // namespace thermalfist

#endif
