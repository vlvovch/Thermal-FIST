/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef MOMENTUMDISTRIBUTION_H
#define MOMENTUMDISTRIBUTION_H

#include <cmath>
#include <vector>

#include "HRGBase/SplineFunction.h"
#include "HRGEventGenerator/Acceptance.h"
#include "HRGEventGenerator/FreezeoutModels.h"

namespace thermalfist {

  /// \brief Class implementing the primordial 3-momentum distribution function
  ///        of certain particle species.
  ///
  /// Assumes that the distribution is azimuthally symmetric
  class MomentumDistributionBase {
  public:
    /**
     * \param pdgid PDG code of particle
     * \param mass  Mass of particle (in GeV)
     */
    MomentumDistributionBase(int pdgid = 0, double mass = 0.) : m_PDGID(pdgid), m_Mass(mass), m_Normalized(false) { }
    
    /// Destructor
    virtual ~MomentumDistributionBase() { }

    /// Normalizes the momentum distribution to unity
    virtual void Normalize() = 0;

    /// Distribution density over the absolute value of the 3-momentum
    virtual double dndp(double p) const = 0;

    /// Distribution density over the longitudinal rapidity
    virtual double dndy(double y) const = 0;

    /// Transverse mass distribution
    virtual double dnmtdmt(double mt) const = 0;

    /// 2D distribution density in rapidity and transverse momentum
    virtual double d2ndptdy(double pt, double y) const = 0;

    /// Whether the distribution has been normalized to unity
    bool isNormalized() const { return m_Normalized; }

    /// Acceptance function in y-pT. If set, the y-pT distribution is
    /// multiplied by the acceptance
    void SetAcceptance(Acceptance::AcceptanceFunction *acc_, double ycm_ = 0.) {
      m_acc = acc_;
      m_useacc = true;
      m_ycm = ycm_;
    }

  protected:
    int m_PDGID;       ///< PDG code of a particle
    double m_Mass;     ///< Mass of a particle
    bool m_Normalized; ///< \copydoc isNormalized()
    Acceptance::AcceptanceFunction *m_acc; ///< Pointer to acceptance function
    double m_ycm;      ///< Center-of-mass rapidity for the acceptance function
    bool m_useacc;     ///< Whether the acceptance functions are used
  };


  /**
   * \brief Class implementing the momentum distribution 
   *        in the spherically symmetric Blast-Wave model of Siemens and Rasmussen
   * 
   * Reference: P. Siemens, J. Rasmussen, Phys. Rev. Lett. 42, 880 (1979)
   */
  class SiemensRasmussenDistribution : public MomentumDistributionBase {
  public:
    /**
     * \copydoc MomentumDistributionBase()
     * \param T The kinetic temperature (in GeV)
     * \param beta The radial flow velocity
     */
    SiemensRasmussenDistribution(int pdgid = 0, double mass = 0., double T = 0.100, double beta = 0.5) :
      MomentumDistributionBase(pdgid, mass),
      m_T(T), m_Beta(beta)
    {
      m_Gamma = 1. / sqrt(1. - m_Beta * m_Beta);
      Normalize();
      m_useacc = false;
    }

    virtual ~SiemensRasmussenDistribution() { }

    /**
     * \brief Set the parameters of the Siemens-Rasmussen distribution
     * 
     * \param T     Kinetic temperature (in GeV)
     * \param beta  Radial flow velocity
     * \param mass  Particle mass (in GeV)
     * \param pdgid Particle PDG code
     */
    void SetParameters(double T, double beta, double mass, int pdgid = 0) {
      m_T = T;
      m_Beta = beta;
      m_Mass = mass;
      if (pdgid != 0) m_PDGID = pdgid;
      m_Gamma = 1. / sqrt(1. - m_Beta * m_Beta);
      Normalize();
    }
 
    // Override functions begin

    void Normalize();

    virtual double dndp(double p) const;

    virtual double dndy(double y) const;

    virtual double dnmtdmt(double mt) const;

    virtual double d2ndptdy(double pt, double y) const;

    // Override functions end

  private:
    double w(double p) const {
      return sqrt(p*p + m_Mass * m_Mass);
    }

    double alpha(double p) const {
      return m_Gamma * m_Beta * p / m_T;
    }

    double PAv() const;

    double m_T;
    double m_Beta;
    double m_Gamma;
    double m_Norm;
    std::vector<double> m_xlag, m_wlag;
  };


  /**
   * \brief Class implementing the momentum distribution
   *        of boost-invariant, azimuthally symmetric freeze-out models
   *        using Maxwell-Boltzmann statistics.
   *
   */
  class BoostInvariantMomentumDistribution : public MomentumDistributionBase {
  public:
    /**
     * \copydoc MomentumDistributionBase()
     * \param FreezeoutModel Pointer to a BoostInvariantFreezeoutParametrization object. Will be deleted on destruction!
     * \param T      The kinetic temperature (in GeV)
     * \param etamax The longitudinal space-time rapidity cut-off
     * \param npow   The power in the transverse flow profile function
     * \param norm   Whether the momentum distribution should be normalized to unity
     */
    BoostInvariantMomentumDistribution(BoostInvariantFreezeoutParametrization* freezeoutModel = NULL, int pdgid = 0, double mass = 0., double T = 0.100, double etamax = 0.5, bool norm = false) :
      MomentumDistributionBase(pdgid, mass),
      m_FreezeoutModel(freezeoutModel),
      m_T(T), m_EtaMax(etamax)
    {
      if (m_FreezeoutModel == NULL) {
        m_FreezeoutModel = new BoostInvariantFreezeoutParametrization();
      }
      m_NormY = m_NormPt = m_Norm = 1.;
      if (norm) Normalize();
      else Initialize();
      m_useacc = false;
    }

    virtual ~BoostInvariantMomentumDistribution();

    void Normalize();

    ///// Rapidity distribution at fixed pT
    //virtual double dndy(double y, double pt) const;

    ///// Rapidity distribution of a single fireball at fixed pT
    //virtual double dndysingle(double y, double pt) const;

    ///// The pT distribution function
    //virtual double dndpt(double pt) const;

    // Override functions begin

    virtual double dndp(double /*p*/) const { return 0.; }

    virtual double dndy(double y) const;

    virtual double dnmtdmt(double mt) const;

    virtual double d2ndptdy(double pt, double y) const;

    // Override functions end

  protected:
    virtual double ZetaIntegrandpTYSingleFireball(double zeta, double pt, double y) const;
    virtual double ZetaIntegrandpT(double zeta, double pt) const;

  private:
    void Initialize();

    virtual double dndysingle(double y) const;

    virtual double dndptsingle(double pt, double y) const;

    /// The pT distribution function
    virtual double dndpt(double pt) const { return pt * dnmtdmt(sqrt(pt * pt + m_Mass * m_Mass)); }

    virtual double dndpt(double pt, double y) const;

    BoostInvariantFreezeoutParametrization* m_FreezeoutModel;

    double m_T;
    double m_EtaMax;
    double m_NormY, m_NormPt, m_Norm;
    std::vector<double> m_xlag, m_wlag;
    std::vector<double> m_xlegT, m_wlegT;
    std::vector<double> m_xlegY, m_wlegY;
    std::vector<double> m_xlegeta, m_wlegeta;

    SplineFunction m_dndy, m_dndyint;
  };

} // namespace thermalfist


/** \defgroup deprecatedEventGenerators Deprecated objects from event generator module.
 *  \deprecated The source code below is no longer used.
 *  @{
 */

namespace thermalfist {

  /**
   * \brief Class implementing the momentum distribution
   *        in the longitudinally symmetric Blast-Wave model
   *
   * \deprecared Superseded by BoostInvariantMomentumDistribution class
   *
   * Reference: E. Schnedermann, J. Sollfrank, U. Heinz, Phys. Rev. C 48, 2462 (1993)
   */
  class SSHDistribution : public MomentumDistributionBase {
  public:
    /**
     * \copydoc MomentumDistributionBase()
     * \param T      The kinetic temperature (in GeV)
     * \param beta   The transverse flow velocity at the surface
     * \param etamax The longitudinal space-time rapidity cut-off
     * \param npow   The power in the transverse flow profile function
     * \param norm   Whether the momentum distribution should be normalized to unity
     */
    SSHDistribution(int pdgid = 0, double mass = 0., double T = 0.100, double betas = 0.5, double etamax = 0.5, double npow = 1., bool norm = false) :
      MomentumDistributionBase(pdgid, mass),
      m_T(T), m_BetaS(betas), m_EtaMax(etamax), m_n(npow)
    {
      m_NormY = m_NormPt = m_Norm = 1.;
      if (norm) Normalize();
      else Initialize();
      m_useacc = false;
    }

    virtual ~SSHDistribution() { }

    /**
     * \brief Set the parameters of the longitudinal blast-wave distribution
     *
     * \param T      The kinetic temperature (in GeV)
     * \param betas  The transverse flow velocity at the surface
     * \param etamax The longitudinal space-time rapidity cut-off
     * \param npow   The power in the transverse flow profile function
     * \param mass   Particle mass (in GeV)
     * \param pdgid  Particle PDG code
     * \param norm   Whether the momentum distribution should be normalized to unity
     */
    void SetParameters(double T, double betas, double etamax, double npow, double mass, int pdgid = 0, bool norm = true) {
      m_T = T;
      m_BetaS = betas;
      m_EtaMax = etamax;
      m_n = npow;
      m_Mass = mass;
      if (pdgid != 0) m_PDGID = pdgid;
      m_NormY = m_NormPt = m_Norm = 1.;
      m_Normalized = false;
      if (norm) Normalize();
      else Initialize();
    }

    /**
     * \brief Set the mean transverse flow velocity.
     *
     * Surface velocity is calculated as \\beta_s = \\langle \\beta_T \\rangle (2+n)/2
     *
     * \param betaT  The mean transverse flow velocity
     */
    void SetMeanBetaT(double betaT) {
      m_BetaS = (2. + m_n) / 2. * betaT;
    }

    void Normalize();

    /// Rapidity distribution at fixed pT
    virtual double dndy(double y, double pt) const;

    /// Rapidity distribution of a single fireball at fixed pT
    virtual double dndysingle(double y, double pt) const;

    /// The pT distribution function
    virtual double dndpt(double pt) const;

    // Override functions begin

    virtual double dndp(double /*p*/) const { return 0.; }

    virtual double dndy(double y) const;

    virtual double dnmtdmt(double mt) const;

    virtual double d2ndptdy(double pt, double y) const;

    // Override functions end

  private:
    void Initialize();

    virtual double dndysingle(double y) const;

    virtual double dndptsingle(double pt, double y) const;


    virtual double dndpt(double pt, double y) const;

    double w(double p) const {
      return sqrt(p * p + m_Mass * m_Mass);
    }

    double asinh(double x) const {
      return log(x + sqrt(1. + x * x));
    }

    double atanh(double x) const {
      return 0.5 * log((1. + x) / (1. - x));
    }

    double betar(double r) const {
      if (m_n == 1.)
        return m_BetaS * r;
      else if (m_n == 2.)
        return m_BetaS * r * r;
      else
        return m_BetaS * pow(r, m_n);
    }

    double rho(double r) const { return atanh(betar(r)); }

    double MtAv() const;

    double y2Av() const;

    double m_T;
    double m_BetaS, m_EtaMax;
    double m_NormY, m_NormPt, m_Norm;
    double m_n;
    std::vector<double> m_xlag, m_wlag;
    std::vector<double> m_xlegT, m_wlegT;
    std::vector<double> m_xlegY, m_wlegY;
    std::vector<double> m_xlegeta, m_wlegeta;

    SplineFunction m_dndy, m_dndyint;
  };

} // namespace thermalfist

 /** @}*/

#endif
