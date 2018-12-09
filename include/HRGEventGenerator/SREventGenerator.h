#ifndef SREVENTGENERATOR_H
#define SREVENTGENERATOR_H

#include <cstdlib>

#include "HRGEventGenerator/EventGeneratorBase.h"
#include "HRGEventGenerator/RandomGenerators.h"

namespace thermalfist {

  //class ExcludedVolumeModel;


  // Class for generating events from Thermal Model with Siemens-Rasmussen momentum distribution
  class SiemensRasmussenEventGenerator : public EventGeneratorBase
  {
  public:
    SiemensRasmussenEventGenerator();
    SiemensRasmussenEventGenerator(ThermalModelBase *THM, double T = 0.120, double beta = 0.5, bool onlyStable = false, int EV = 0, ThermalModelBase *THMEVVDW = NULL);
    ~SiemensRasmussenEventGenerator() { }

    void SetParameters(double T, double beta);

    void SetThermalModel(ThermalModelBase *THM, bool regen = true);

    void SetMode(bool OnlyStable) { m_OnlyStable = OnlyStable; }

    void SetMomentumGenerators();

  private:
    double m_T;
    double m_Beta;
  };

} // namespace thermalfist

#endif
