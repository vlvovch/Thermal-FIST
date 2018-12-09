#ifndef PARTICLEDECAYS_H
#define PARTICLEDECAYS_H

#define DEBUGDECAYS

#include <cmath>
#include <vector>
#include <map>

#include "HRGEventGenerator/SimpleParticle.h"

namespace thermalfist {

  namespace ParticleDecays {
    double ThreeBodym12F2(double m12, double M, double m1, double m2, double m3);
    double TernaryThreeBodym12Maximum(double M, double m1_, double m2_, double m3_);
    extern int threebodysucc, threebodytot;
    double GetRandomThreeBodym12(double M, double m1_, double m2_, double m3_, double fm12max);

    SimpleParticle LorentzBoost(const SimpleParticle &part, double vx, double vy, double vz);
    std::vector<SimpleParticle> TwoBodyDecay(const SimpleParticle & Mother, double m1, int pdg1, double m2, int pdg2);
    std::vector<SimpleParticle> ManyBodyDecay(const SimpleParticle & Mother, std::vector<double> masses, std::vector<int> pdgs); // TODO: proper implementation for 4+ - body decays
  }

} // namespace thermalfist

#endif
