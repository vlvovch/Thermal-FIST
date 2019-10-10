/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/FreezeoutModels.h"

#include <fstream>
#include <iostream>

namespace thermalfist {

  const double BoostInvariantFreezeoutParametrization::dzeta = 0.001;

  double BoostInvariantFreezeoutParametrization::ZetaProbability(double zeta) const
  {
    return Rfunc(zeta) * taufunc(zeta) * (coshetaperp(zeta) * dRdZeta(zeta) - sinhetaperp(zeta) * dtaudZeta(zeta));
  }

  double BoostInvariantFreezeoutParametrization::ProbabilityMaximum()
  {
    if (!m_ProbabilityMaximumComputed) {
      m_ProbabilityMaximum = ComputeProbabilitydMaximum();
      m_ProbabilityMaximumComputed = true;
    }
    return m_ProbabilityMaximum;
  }

  double BoostInvariantFreezeoutParametrization::ComputeProbabilitydMaximum()
  {
    // Global ternary search
    m_ProbabilityMaximum = TernarySearchForIntegrandMaximum(0., 1.);

    // Look for a possibility that the ternary search has produced a local minimum instead of the global one
    double tmax = 0., tzetamax = 0.;
    double dzeta = 0.01;
    for (double tzeta = 0.; tzeta <= 1. + 1.e9; tzeta += dzeta) {
      double tprob = ZetaProbability(tzeta);
      if (tprob > tmax) {
        tmax = tprob;
        tzetamax = tzeta;
      }
    }

    if (tmax > m_ProbabilityMaximum) {
      m_ProbabilityMaximum = TernarySearchForIntegrandMaximum(tzetamax - dzeta, tzetamax + dzeta);
    }

    return m_ProbabilityMaximum;
  }

  double BoostInvariantFreezeoutParametrization::TernarySearchForIntegrandMaximum(double zetaMin, double zetaMax) const
  {
    double eps = 1e-8;
    double l = zetaMin, r = zetaMax;

    if (l < 0.)
      l = 0.;
    if (r > 1.)
      r = 1.;

    double m1 = l + (r - l) / 3.;
    double m2 = r - (r - l) / 3.;
    int MAXITERS = 200;
    int iter = 0;
    while (fabs(m2 - m1) > eps && iter < MAXITERS) {
      if (ZetaProbability(m1) < ZetaProbability(m2)) {
        l = m1;
      }
      else {
        r = m2;
      }
      m1 = l + (r - l) / 3.;
      m2 = r - (r - l) / 3.;
      iter++;
    }
    return ZetaProbability((m1 + m2) / 2.);
  }

  CylindricalBlastWaveParametrization::CylindricalBlastWaveParametrization(double betaSurface, double nPower, double tau, double Rmax) :
    BoostInvariantFreezeoutParametrization(),
    m_BetaS(betaSurface),
    m_n(nPower),
    m_tau(tau),
    m_R(Rmax)
  {
    if (tau <= 0. || Rmax <= 0. || m_n < 0. || (m_BetaS < 0. || m_BetaS > 1.)) {
      std::cerr << "**ERROR** CylindricalBlastWaveParametrization::CylindricalBlastWaveParametrization: invalid parameter values!" << std::endl;
      exit(1);
    }
  }

  double CylindricalBlastWaveParametrization::ZetaProbability(double zeta) const
  {
    return m_R * zeta * m_tau * coshetaperp(zeta) * m_R;
  }

  CracowFreezeoutParametrization::CracowFreezeoutParametrization(double RoverTauH, double tauH) :
    BoostInvariantFreezeoutParametrization(),
    m_RoverTauH(RoverTauH),
    m_tauH(tauH)
  {
    if (tauH <= 0. || RoverTauH <= 0.) {
      std::cerr << "**ERROR** CracowFreezeoutParametrization::CracowFreezeoutParametrization: invalid parameter values!" << std::endl;
      exit(1);
    }
  }

  double CracowFreezeoutParametrization::ZetaProbability(double zeta) const
  {
    return Rmax() * zeta * m_tauH * Rmax();
  }

} // namespace thermalfist