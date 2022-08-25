/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef CRACOWFREEZEOUTEVENTGENERATOR_H
#define CRACOWFREEZEOUTEVENTGENERATOR_H

#include <cstdlib>

#include "HRGEventGenerator/EventGeneratorBase.h"
#include "HRGEventGenerator/RandomGenerators.h"

namespace thermalfist {

  /// \brief Class implementing the Thermal Event Generator for
  ///        the Cracow freeze-out model scenario
  class CracowFreezeoutEventGenerator : public EventGeneratorBase
  {
  public:

    CracowFreezeoutEventGenerator();

    /**
     * \brief Construct a new CylindricalBlastWaveEventGenerator object
     *
     * \param TPS        A pointer to the particle list
     * \param config     Event generator configuration
     * \param Tkin       The kinetic freeze-out temperature (in GeV)
     * \param RoverTauH  The r_max/tau_H ratio parameter
     * \param etamax     The longitudinal space-time rapidity cut-off
     */
    CracowFreezeoutEventGenerator(ThermalParticleSystem* TPS,
      const EventGeneratorConfiguration& config,
      double T = 0.120,
      double RoverTauH = 1.5,
      double etamax = 0.5);

    ~CracowFreezeoutEventGenerator() { }

    /// Sets the momentum distribution parameters
    void SetParameters(double T, double RoverTauH, double etamax);


    /// Sets up the random generators of particle momenta
    /// and resonances masses
    virtual void SetMomentumGenerators();

    double GetTkin() const { return m_T; }
    double GetRoverTauH() const { return m_RoverTauH; }
    double GetEtaMax() const { return m_EtaMax; }
  private:
    double m_T, m_RoverTauH, m_EtaMax;
  };

} // namespace thermalfist

#endif
