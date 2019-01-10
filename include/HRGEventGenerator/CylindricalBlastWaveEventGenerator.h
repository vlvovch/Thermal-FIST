/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef CYLINDRICALBLASTWAVEEVENTGENERATOR_H
#define CYLINDRICALBLASTWAVEEVENTGENERATOR_H

#include <cstdlib>

#include "HRGEventGenerator/EventGeneratorBase.h"
#include "HRGEventGenerator/RandomGenerators.h"

namespace thermalfist {

  /// \brief Class implementing the Thermal Event Generator for
  ///        the longitudinally symmetric blast-wave scenario
  class CylindricalBlastWaveEventGenerator : public EventGeneratorBase
  {
  public:
    
    CylindricalBlastWaveEventGenerator();

    /**
     * \brief Construct a new CylindricalBlastWaveEventGenerator object
     * 
     * \param TPS    A pointer to the particle list
     * \param config Event generator configuration
     * \param Tkin   The kinetic freeze-out temperature (in GeV) 
     * \param beta   The transverse flow velocity
     * \param etamax The longitudinal space-time rapidity cut-off
     * \param npow   The power in the transverse flow profile function
     */
    CylindricalBlastWaveEventGenerator(ThermalParticleSystem *TPS, 
                                       const EventGeneratorConfiguration& config, 
                                       double T = 0.120, 
                                       double beta = 0.5, 
                                       double etamax = 0.5, 
                                       double npow = 1.);
    
    /// \deprecated
    /// \brief Old constructor. Included for backward compatibility.
    CylindricalBlastWaveEventGenerator(ThermalModelBase *THM, double T = 0.120, double beta = 0.5, double etamax = 0.5, double npow = 1., bool onlyStable = false, EventGeneratorConfiguration::ModelType EV = EventGeneratorConfiguration::PointParticle, ThermalModelBase *THMEVVDW = NULL);
    
    ~CylindricalBlastWaveEventGenerator() { }

    /// Sets the momentum distribution parameters
    void SetParameters(double T, double beta, double etamax, double npow = 1.);

    /// Sets up the random generators of particle momenta
    /// and resonances masses
    void SetMomentumGenerators();
  private:
    double m_T, m_Beta, m_EtaMax, m_n;
  };

  /// For backward compatibility
  typedef CylindricalBlastWaveEventGenerator SSHEventGenerator;
} // namespace thermalfist

#endif
