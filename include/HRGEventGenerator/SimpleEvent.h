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
  /// Structure holding information about a single event in the event generator.
  struct SimpleEvent {
    /// Event weight factor
    double weight;

    /// Log of the event weight factor
    double logweight;

    /// Vector of all final particles in the event
    std::vector<SimpleParticle> Particles;

     /// Vector of all particles which ever appeared in the event (including those that decay)
    std::vector<SimpleParticle> AllParticles;

    /// Default constructor, empty event
    SimpleEvent() { Particles.resize(0); AllParticles.resize(0); weight = 1.; }

    /// Writes the event to an output file stream
    void writeToFile(std::ofstream & fout, int eventnumber = 1);

    /// Rapidity boost for all particles
    void RapidityBoost(double dY);

    static SimpleEvent MergeEvents(const SimpleEvent &evt1, const SimpleEvent &evt2);
  };

} // namespace thermalfist

#endif
