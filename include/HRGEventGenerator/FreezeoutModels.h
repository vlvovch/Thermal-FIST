/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef FREEZEOUTMODELS_H
#define FREEZEOUTMODELS_H

//#include "HRGEventGenerator/RandomGenerators.h"
#include <cmath>

namespace thermalfist {

  /**
   * \brief Base class implementing a longitudinally boost-invariant
   *        azimuthally symmetric freeze-out parametrization.
   *
   *
   * The hypersurface is parametrized using variables:
   *  - Space-time rapidity \eta
   *  - Azimuthal angle \phi
   *  - Variable \zeta defining which parametrizes the hypersurface in \tau-R_T coordinates
   *
   * This base class implements a static fireball (no transverse flow) of 1 fm transverse radius at \tau = 1 fm/c proper time
   */
  class BoostInvariantFreezeoutParametrization {
  public:

    BoostInvariantFreezeoutParametrization() : m_ProbabilityMaximumComputed(false), m_ProbabilityMaximum(1.) {}
    virtual ~BoostInvariantFreezeoutParametrization() {}

    /**
     * \brief Transverse radius vs \zeta
     */
    virtual double Rfunc(double zeta) const { return zeta; }

    /**
     * \brief dR/d\zeta
     */
    virtual double dRdZeta(double zeta) const { return (Rfunc(zeta + dzeta) - Rfunc(zeta)) / dzeta; }

    /**
     * \brief Proper time \tau vs \zeta
     */
    virtual double taufunc(double zeta) const { return 1.; }

    /**
     * \brief d\tau/d\zeta
     */
    virtual double dtaudZeta(double zeta) const { return (taufunc(zeta + dzeta) - taufunc(zeta)) / dzeta; }

    /**
     * \brief Transverse flow rapidity as a function of \zeta
     */
    virtual double etaperp(double zeta) const { return 0.; }

    virtual double sinhetaperp(double zeta) const { return sinh(etaperp(zeta)); }
    virtual double coshetaperp(double zeta) const { return cosh(etaperp(zeta)); }
    virtual double tanhetaperp(double zeta) const { return sinhetaperp(zeta) / coshetaperp(zeta); }

    /**
     * \brief Proportional to probability of having given \zeta value
     *
     * Given by d\Sigma_\mu u^\mu.
     * Used by Monte Carlo sampler.
     */
    virtual double ZetaProbability(double zeta) const;// { return zeta; }

    virtual double ProbabilityMaximum();

    /**
     * \brief Samples zeta for use in Monte Carlo event generator.
     *
     * Uses rejection sampling.
     *
     */
    //virtual double GetRandomZeta(MTRand& rangen = RandomGenerators::randgenMT);

    /**
    *  \brief Whether explicit inverse function of \zeta variable distribution is available.
    */
    virtual bool InverseZetaDistributionIsExplicit() const { return false; }

    /**
    *  \brief Inverse function of \zeta variable distribution used in random number generation.
    */
    virtual double InverseZetaDistribution(double xi) const { return 0.; }


  protected:
    /**
    * \brief Computes and sets the maximum of the \zeta probability density.
    */
    virtual double ComputeProbabilitydMaximum();


  private:
    /**
    * \brief Searches for the maximum of ZetaProbability(zeta) in \zeta range [zetaMin,zetaMax].
    *
    * Uses ternary search.
    */
    virtual double TernarySearchForIntegrandMaximum(double zetaMin = 0., double zetaMax = 1.) const;

    static const double dzeta;

    bool m_ProbabilityMaximumComputed;
    double m_ProbabilityMaximum;
  };

  /**
   * \brief Implements the cylindrically symmetric blast-wave model parametrization.
   *
   *
   * Reference: E. Schnedermann, J. Sollfrank, U. Heinz, Phys. Rev. C 48, 2462 (1993)
   */
  class CylindricalBlastWaveParametrization : public BoostInvariantFreezeoutParametrization {
  public:

    /**
     * \param betaSurface Transverse flow velocity at the surface
     * \param nPower      The power in the transverse flow profile function
     * \param tau         Freeze-out proper time (fm/c)
     * \param Rmax        Transverse radius (fm)
     *
     * The values of parameters tau and Rmax do not influence the shape
     * of the momenum distribution. Their values are irrelevant for the sampling of the momenta
     */
    CylindricalBlastWaveParametrization(double betaSurface = 0.5, double nPower = 1., double tau = 10., double Rmax = 6.);

    virtual ~CylindricalBlastWaveParametrization() {}

    // Override functions begin

    virtual double Rfunc(double zeta) const { return zeta * m_R; }

    virtual double dRdZeta(double zeta) const { return m_R; }

    virtual double taufunc(double zeta) const { return m_tau; }

    virtual double dtaudZeta(double zeta) const { return 0.; }

    virtual double etaperp(double zeta) const { return atanh(m_BetaS * pow(zeta, m_n)); }


    virtual double sinhetaperp(double zeta) const { return m_BetaS * pow(zeta, m_n) / sqrt(1. - m_BetaS * m_BetaS * pow(zeta, 2. * m_n)); }
    virtual double coshetaperp(double zeta) const { return 1. / sqrt(1. - m_BetaS * m_BetaS * pow(zeta, 2. * m_n)); }
    virtual double tanhetaperp(double zeta) const { return m_BetaS * pow(zeta, m_n); }

    virtual double ZetaProbability(double zeta) const;

  protected:
    virtual double ComputeProbabilitydMaximum() { return ZetaProbability(1.); }

    // Override functions end

  private:
    double m_BetaS, m_n;
    double m_tau, m_R;
  };

  /**
   * \brief Implements the Cracow (Hubble-like) freeze-out model parametrization.
   *
   *
   * Reference: W. Broniowski, W. Florkowski, Phys. Rev. Lett. 87, 272302 (2001)
   */
  class CracowFreezeoutParametrization : public BoostInvariantFreezeoutParametrization {
  public:

    /**
     * \param RoverTauH   Ratio of maximum transverse radius over the freeze-out (Hubble) proper time
     * \param tauH        The freeze-out (Hubble) proper time
     *
     * The value of parameter tauH does not influence the shape
     * of the momenum distribution. Its value is irrelevant for the sampling of the momenta
     */
    CracowFreezeoutParametrization(double RoverTauH = 1., double tauH = 10.);

    virtual ~CracowFreezeoutParametrization() {}

    // Override functions begin

    virtual double Rfunc(double zeta) const { return zeta * Rmax(); }

    virtual double dRdZeta(double zeta) const { return Rmax(); }

    virtual double taufunc(double zeta) const { return m_tauH * coshetaperp(zeta); }

    virtual double dtaudZeta(double zeta) const { return Rmax() * sinhetaperp(zeta) / coshetaperp(zeta); }

    virtual double etaperp(double zeta) const { return asinh(zeta * m_RoverTauH); }


    virtual double sinhetaperp(double zeta) const { return zeta * m_RoverTauH; }
    virtual double coshetaperp(double zeta) const { return sqrt(1. + zeta * zeta * m_RoverTauH * m_RoverTauH); }

    virtual double ZetaProbability(double zeta) const;

    //virtual double GetRandomZeta(MTRand& rangen = RandomGenerators::randgenMT);

  protected:
    virtual double ComputeProbabilitydMaximum() { return ZetaProbability(1.); }

    virtual bool InverseZetaDistributionIsExplicit() const { return true; }

    virtual double InverseZetaDistribution(double xi) const { return sqrt(xi); }

    // Override functions end

  private:
    double Rmax() const { return  m_RoverTauH * m_tauH; }

    double m_RoverTauH, m_tauH;
  };

} // namespace thermalfist

#endif
