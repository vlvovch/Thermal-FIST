/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef PARTICLEDECAYS_H
#define PARTICLEDECAYS_H

#include <cmath>
#include <vector>
#include <map>

#include "HRGEventGenerator/SimpleParticle.h"

namespace thermalfist {

  /// Contains functions for Monte Carlo generation of decays
  namespace ParticleDecaysMC {
    /**
     * \brief Square of the \f$m_{12}\f$ probability density function 
     *        used in a three-body decay (unnormalized)
     * 
     * \param m12 Invariant mass of the leading two daughter particles in a three-body decay
     * \param M   Mass of the decaying particle
     * \param m1  Mass of the first daughter particle
     * \param m2  Mass of the first daughter particle
     * \param m3  Mass of the first daughter particle
     * \return    Square of the probability density for \f$m_{12}\f$
     */
    double ThreeBodym12F2(double m12, double M, double m1, double m2, double m3);

    /**
     * \brief Determines the maximum of the \f$m_{12}\f$ probability density function 
     *        used in a three-body decay.
     * 
     * Uses ternary search.
     * 
     * \param M   Mass of the decaying particle
     * \param m1  Mass of the first daughter particle
     * \param m2  Mass of the first daughter particle
     * \param m3  Mass of the first daughter particle
     * \return    Maximum of the \f$m_{12}\f$ probability density function
     */
    double TernaryThreeBodym12Maximum(double M, double m1, double m2, double m3);

    /// Used for debugging the succes rate in the rejection sampling of \f$m_{12}\f$.
    extern int threebodysucc, threebodytot;

    /**
     * \brief Sample the invariant mass \f$m_{12}\f$ of the leading two
     *        daughter particles in a three-body decay.
     * 
     * \param M   Mass of the decaying particle
     * \param m1  Mass of the first daughter particle
     * \param m2  Mass of the first daughter particle
     * \param m3  Mass of the first daughter particle
     * \param fm12max The maximum of the \f$m_{12}\f$ probability density (precomputed with TernaryThreeBodym12Maximum())
     * \return    The sampled value of \f$m_{12}\f$
     */
    double GetRandomThreeBodym12(double M, double m1, double m2, double m3, double fm12max);

    /**
     * \brief Lorentz boost of the 4-momentum and 4-coordinate of a particle
     * 
     * \param[in] part Input particle with a 4-momentum
     * \param[in] vx   Lorentz boost velocity x component
     * \param[in] vy   Lorentz boost velocity y component
     * \param[in] vz   Lorentz boost velocity z component
     * \return SimpleParticle Particle with Lorentz-boosted 4-momentum and 4-coordinate
     */
    SimpleParticle LorentzBoost(const SimpleParticle& part, double vx, double vy, double vz);

    /**
     * \brief Lorentz boost of the 4-momentum of a particle
     *
     *        Does *not* boost the coordinates!
     *
     * \param[in] part Input particle with a 4-momentum
     * \param[in] vx   Lorentz boost velocity x component
     * \param[in] vy   Lorentz boost velocity y component
     * \param[in] vz   Lorentz boost velocity z component
     * \return SimpleParticle Particle with Lorentz-boosted 4-momentum
     */
    SimpleParticle LorentzBoostMomentumOnly(const SimpleParticle& part, double vx, double vy, double vz);

    /**
     * \brief Lorentz boost of the 4-coordinate and 4-momentum of a particle
     *
     * \param[in] part Input particle with a 4-momentum
     * \param[in] vx   Lorentz boost velocity x component
     * \param[in] vy   Lorentz boost velocity y component
     * \param[in] vz   Lorentz boost velocity z component
     * \return SimpleParticle Particle with Lorentz-boosted 4-momentum
     */
    SimpleParticle LorentzBoostMomentaAndCoordinates(const SimpleParticle& part, double vx, double vy, double vz);

    /**
     * \brief Samples the decay products of a two-body decay.
     * 
     * \param Mother The decaying particle
     * \param m1     Mass of the first daughter (in GeV)
     * \param pdg1   Pdg code of the first daughter
     * \param m2     Mass of the second daughter (in GeV)
     * \param pdg2   Pdg code of the second daughter
     * \return std::vector<SimpleParticle> Two-component vector of decay products
     */
    std::vector<SimpleParticle> TwoBodyDecay(const SimpleParticle & Mother, double m1, long long pdg1, double m2, long long pdg2);
    
    /**
     * \brief Samples the decay products of a many-body decay.
     * 
     * ---
     * **NOTE**
     * 
     * While the decay kinematics of a three-body decay here are exact,
     * decay kinematics for (4+)-body decays are approximate and should
     * be taken into consideration.
     * 
     * ---
     * 
     * \param Mother The decaying particle
     * \param masses Masses of the decay products (in GeV)
     * \param pdgs   Pdg codes of the decay products
     * \return std::vector<SimpleParticle> 
     */
    std::vector<SimpleParticle> ManyBodyDecay(const SimpleParticle & Mother, std::vector<double> masses, std::vector<long long> pdgs); // TODO: proper implementation for 4+ - body decays


    /**
     * \brief Shuffles the decay products.
     *        
     * Currently used for four+ body decays to avoid having asymmetry in 
     * momentum distributions of decay products due to non-exact current implementation of these decays.
     *
     * \param masses Masses of the decay products (in GeV)
     * \param pdgs   Pdg codes of the decay products
     */
    void ShuffleDecayProducts(std::vector<double> &masses, std::vector<long long> &pdgs); // TODO: proper implementation for 4+ - body decays
  
    
    /**
     * \brief Return the square of the distance between particles at equal time.
     *
     * Makes Lorentz boost into the pair rest frame, then propagates the earlier particle
     * to other particle's time, and computed the square of the distance between particles at that time.
     *
     * \param part1  First particle
     * \param part2  Second particle
     */
    double ParticleDistanceSquared(const SimpleParticle& part1, const SimpleParticle& part2);
  }

} // namespace thermalfist

#endif
