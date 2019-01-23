/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef PARTICLEDECAY_H
#define PARTICLEDECAY_H

/**
 * \file ParticleDecay.h
 * \brief Contains structures related to particle decays.
 * 
 */

#include <string>
#include <vector>
#include <cmath>

#include "HRGBase/ThermalModelParameters.h"
#include "HRGBase/IdealGasFunctions.h"
#include "HRGBase/xMath.h"

namespace thermalfist {
  /**
   * \brief An auxiliary struct containing the list of feeddown flags.
   * 
   */
  struct Feeddown {
    /**
     * List of feeddown flags.
     * 
     */
    enum Type { 
      Primordial = 0,      ///< No feeddown.
      StabilityFlag = 1,   ///< Feeddown from all particles marked as unstable.
      Weak = 2,            ///< Feeddown from strong, electromagnetic, and weak decays.
      Electromagnetic = 3, ///< Feeddown from strong and electromagnetic decays.
      Strong = 4           ///< Feeddown from strong decays.
      };
    static const int NumberOfTypes = 5;
  };

  /**
  * \brief An auxiliary struct containing the list of decay types.
  *
  */
  struct ParticleDecayType
  {
    /**
     * \brief Type of particle's decay.
     * 
     * Stable -- does not decay.
     * Default -- not used.
     * Weak -- weakly decaying.
     * Electromagnetic -- electromagnetically decaying.
     * Strong -- strongly decaying.
     */
    enum DecayType { 
      Stable = 0,           ///< Does not decay
      Default = 1,          ///< Not used
      Weak = 2,             ///< Weakly decaying
      Electromagnetic = 3,  ///< Electromagnetically decaying
      Strong = 4            ///< Strongly decaying
      };
    static const int NumberOfDecayTypes = 5;
  };

  /**
   * \brief Structure containing information about a single decay channel of a particle.
   */
  struct ParticleDecayChannel {
    double mBratio;                      /**< Branching ratio */
    std::vector<long long> mDaughters;   /**< PDGID numbers of daughter particles */
    double mM0;                          /**< Sum of masses of decay products */

    double mPole;
    double mL;                      /**< Orbital angular momentum for decay, used in the eBW scheme */
    std::vector<double> mBratioVsM; /**< Mass-dependent branching ratios */
    double mBratioAverage;          /**< Average branching ratios after integrating with the Boltzmann factor*/

    /// Name of the decay channel. Not used.
    std::string mChannelName;

    /**
     * \brief Construct a new ParticleDecay object.
     * 
     * \param bratio Branching ratio of the decay. Between 0 and 1.
     * \param daughters A vector of PDG ID numbers of all daughter products.
     */
    ParticleDecayChannel(double bratio = 0., const std::vector<long long> &daughters = std::vector<long long>(0)) :
      mBratio(bratio), mDaughters(daughters), mM0(0.), mPole(0.), mL(0.),
      mBratioVsM(std::vector<double>(0)), mBratioAverage(bratio), mChannelName("decay1") {
    }

    /**
     * \brief Energy depedent modified branching ratio.
     * 
     * Energy dependent modification of the branching ratio.
     * Used in the eBW scheme.
     * 
     * \param m Energy [GeV].
     * \return Modification of the branching ratio.
     */
    double ModifiedWidth(double m) const;

    bool operator==(const ParticleDecayChannel &rhs) const; // TODO: improve
    bool operator!=(const ParticleDecayChannel &rhs) const { return !(*this == rhs); }
  };  

} // namespace thermalfist

#endif