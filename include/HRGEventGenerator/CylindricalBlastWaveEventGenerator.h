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

  // Class for generating events from Thermal Model with Schnedermann-Sollfrank-Heinz (Blast-Wave) prescription for momentum distribution
  class CylindricalBlastWaveEventGenerator : public EventGeneratorBase
  {
  public:
    CylindricalBlastWaveEventGenerator();
    CylindricalBlastWaveEventGenerator(ThermalModelBase *THM, double T = 0.120, double beta = 0.5, double etamax = 0.5, double npow = 1., bool onlyStable = false, EventGeneratorConfiguration::ModelType EV = EventGeneratorConfiguration::PointParticle, ThermalModelBase *THMEVVDW = NULL);
    ~CylindricalBlastWaveEventGenerator() { }

    void SetParameters(double T, double beta, double etamax);

    void SetThermalModel(ThermalModelBase *THM, bool regen = true);

    void SetMode(bool OnlyStable) { m_OnlyStable = OnlyStable; }

    void SetMomentumGenerators();
  private:
    double m_T, m_Beta, m_EtaMax, m_n;
  };

  /// For backward compatibility
  typedef CylindricalBlastWaveEventGenerator SSHEventGenerator;
} // namespace thermalfist

#endif
