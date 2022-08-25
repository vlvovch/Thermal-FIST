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
    
    //CylindricalBlastWaveEventGenerator();

    /**
     * \brief Construct a new CylindricalBlastWaveEventGenerator object
     * 
     * \param TPS    A pointer to the particle list
     * \param config Event generator configuration
     * \param Tkin   The kinetic freeze-out temperature (in GeV) 
     * \param betas  The transverse flow velocity at the surface
     * \param etamax The longitudinal space-time rapidity cut-off
     * \param npow   The power in the transverse flow profile function
     */
    CylindricalBlastWaveEventGenerator(ThermalParticleSystem *TPS = NULL, 
                                       const EventGeneratorConfiguration& config = EventGeneratorConfiguration(),
                                       double T = 0.120, 
                                       double betas = 0.5, 
                                       double etamax = 0.5, 
                                       double npow = 1.,
                                       double Rperp = 6.5);
    
    /// \deprecated
    /// \brief Old constructor. Included for backward compatibility.
    CylindricalBlastWaveEventGenerator(ThermalModelBase *THM, 
      double T = 0.120, 
      double betas = 0.5, 
      double etamax = 0.5, 
      double npow = 1., 
      bool onlyStable = false, 
      EventGeneratorConfiguration::ModelType EV = EventGeneratorConfiguration::PointParticle, 
      ThermalModelBase *THMEVVDW = NULL);
    
    ~CylindricalBlastWaveEventGenerator() { }

    /// Sets the momentum distribution parameters
    void SetParameters(double T, double betas, double etamax, double npow = 1.);

    /**
     * \brief Set the mean transverse flow velocity.
     *
     * Surface velocity is calculated as \\beta_s = \\langle \\beta_T \\rangle (2+n)/2 
     * 
     * \param betaT  The mean transverse flow velocity
     */
    void SetMeanBetaT(double betaT);

    /// Sets up the random generators of particle momenta
    /// and resonances masses
    virtual void SetMomentumGenerators();

    void RecalculateTotalConservedNumbers();

    //virtual void SetParameters();

    double GetTkin() const { return m_T; }
    double GetBetaSurface() const { return m_BetaS; }
    double GetNPow() const { return m_n; }
    double GetEtaMax() const { return m_EtaMax; }
    double GetRperp() const { return m_Rperp; }
  private:
    double GetVeffIntegral() const;

    double m_T, m_BetaS, m_EtaMax, m_n, m_Rperp;
  };

  /// For backward compatibility
  typedef CylindricalBlastWaveEventGenerator SSHEventGenerator;

} // namespace thermalfist

#endif
