/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef SPHERICALBLASTWAVEEVENTGENERATOR_H
#define SPHERICALBLASTWAVEEVENTGENERATOR_H

#include <cstdlib>

#include "HRGEventGenerator/EventGeneratorBase.h"
#include "HRGEventGenerator/RandomGenerators.h"

namespace thermalfist {


  /// \brief Class implementing the Thermal Event Generator for
  ///        the isotropic blast-wave scenario
  ///
  /// Momentum distribution is based on the Siemens-Rasmussen formula
  ///        
  class SphericalBlastWaveEventGenerator : public EventGeneratorBase
  {
  public:

    SphericalBlastWaveEventGenerator();

    /**
     * \brief Construct a new SphericalBlastWaveEventGenerator object
     * 
     * \param TPS    A pointer to the particle list
     * \param config Event generator configuration
     * \param Tkin   The kinetic freeze-out temperature (in GeV)      
     * \param beta   Radial flow velocity
     */
    SphericalBlastWaveEventGenerator(ThermalParticleSystem *TPS, 
                                     const EventGeneratorConfiguration& config,
                                     double Tkin = 0.120, 
                                     double beta = 0.5);

    /// \deprecated
    /// \brief Old constructor. Included for backward compatibility.
    SphericalBlastWaveEventGenerator(ThermalModelBase *THM, double T = 0.120, double beta = 0.5, bool onlyStable = false, EventGeneratorConfiguration::ModelType EV = EventGeneratorConfiguration::PointParticle, ThermalModelBase *THMEVVDW = NULL);
    
    ~SphericalBlastWaveEventGenerator() { }

    /// Sets the momentum distribution parameters
    void SetBlastWaveParameters(double Tkin, double beta);

    /// Sets up the random generators of particle momenta
    /// and resonances masses
    void SetMomentumGenerators();

    double GetTkin() { return m_T; }
    double GetBeta() { return m_Beta; }

    //virtual void SetParameters();

  private:
    double m_T;
    double m_Beta;
  };

  /// For backward compatibility
  typedef SphericalBlastWaveEventGenerator SiemensRasmussenEventGenerator;
} // namespace thermalfist

#endif
