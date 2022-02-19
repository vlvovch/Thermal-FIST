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

    /// Vector of all particles which ever appeared in the event (including those that decay and photons/leptons)
    std::vector<SimpleParticle> AllParticles;

    /// Vector of all decay photons/leptons
    std::vector<SimpleParticle> PhotonsLeptons;

    /// Vector for each AllParticles element pointing to the index of the mother resonance.
    /// If the element corresponds to a primordial particle, the id is -1.
    std::vector<int> DecayMap;

    /// Vector for each Particles element pointing to the index of the primordial resonance from which this particle originated.
    std::vector<int> DecayMapFinal;

    /// Default constructor, empty event
    SimpleEvent() { Particles.resize(0); AllParticles.resize(0); PhotonsLeptons.resize(0); weight = 1.; logweight = 0.; }

    /// Rapidity boost by Y -> Y + dY for all particles
    void RapidityBoost(double dY);

    /// Merge particles from two events (e.g. two patches, two canonical volumes, etc.)
    static SimpleEvent MergeEvents(const SimpleEvent &evt1, const SimpleEvent &evt2);

    /// Configuration for the event output
    struct EventOutputConfig {

      /// Output the particle's energy in addition to its 3-momentum
      bool printEnergy;

      /// Output the pdg code of the mother particle
      bool printMotherPdg;

      /// Output photons and leptons, if any
      bool printPhotonsLeptons;

      /// Print the number of successive decays before the particle was produced
      bool printDecayEpoch;

      /// Print the space-time coordinates of the particles
      bool printCoordinates;

      /// Print the event weight for importance sampling
      bool printWeight;

      EventOutputConfig() :
        printEnergy(true), printMotherPdg(false), printPhotonsLeptons(false), printDecayEpoch(false), printCoordinates(false), printWeight(true) { }
    };

    /// Writes the event to an output file stream
    void writeToFile(std::ofstream& fout, const EventOutputConfig& config = EventOutputConfig(), int eventnumber = 1) const;

    /// Writes the event to an output file stream
    void writeToFile(std::ofstream& fout, int eventnumber = 1) const { writeToFile(fout, EventOutputConfig(), eventnumber); }

    /// Writes the event in a format suitable for UrQMD afterburner, as described here https://github.com/jbernhard/urqmd-afterburner
    void writeToFileForUrqmd(std::ofstream& fout) const;
  };

} // namespace thermalfist

#endif
