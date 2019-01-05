/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef SPHERICALBLASTWAVEEVENTGENERATOR_H
#define SPHERICALBLASTWAVEEVENTGENERATOR_H

#include <cstdlib>

#include "HRGEventGenerator/EventGeneratorBase.h"
#include "HRGEventGenerator/RandomGenerators.h"

namespace thermalfist {


  // Class for generating events from Thermal Model with Siemens-Rasmussen momentum distribution
  class SphericalBlastWaveEventGenerator : public EventGeneratorBase
  {
  public:
    SphericalBlastWaveEventGenerator();
    SphericalBlastWaveEventGenerator(ThermalModelBase *THM, double T = 0.120, double beta = 0.5, bool onlyStable = false, EventGeneratorConfiguration::ModelType EV = EventGeneratorConfiguration::PointParticle, ThermalModelBase *THMEVVDW = NULL);
    ~SphericalBlastWaveEventGenerator() { }

    void SetParameters(double T, double beta);

    void SetThermalModel(ThermalModelBase *THM, bool regen = true);

    void SetMode(bool OnlyStable) { m_OnlyStable = OnlyStable; }

    void SetMomentumGenerators();

  private:
    double m_T;
    double m_Beta;
  };

  /// For backward compatibility
  typedef SphericalBlastWaveEventGenerator SiemensRasmussenEventGenerator;
} // namespace thermalfist

#endif
