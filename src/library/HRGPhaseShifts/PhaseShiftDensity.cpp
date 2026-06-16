/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGPhaseShifts/PhaseShiftDensity.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "HRGBase/NumericalIntegration.h"
#include "HRGBase/xMath.h"

namespace thermalfist {

  PhaseShiftDensity::PhaseShiftDensity(const std::vector<PhaseShiftPartialWave>& waves,
                                       double m1, double m2, double Mmax,
                                       int statistics, int quadratureNodes)
    : m_waves(waves), m_m1(m1), m_m2(m2), m_Mmax(Mmax), m_stat(statistics), m_nodes(quadratureNodes), m_enabled(true)
  {
    if (m1 <= 0. || m2 <= 0.)
      throw std::invalid_argument("PhaseShiftDensity: constituent masses must be positive");
    if (Mmax <= m1 + m2)
      throw std::invalid_argument("PhaseShiftDensity: Mmax must exceed the threshold m1 + m2");
    if (quadratureNodes <= 0)
      throw std::invalid_argument("PhaseShiftDensity: quadratureNodes must be positive");
    m_qmax = qfromM(m_Mmax);
  }

  double PhaseShiftDensity::Mfromq(double q) const {
    return std::sqrt(q * q + m_m1 * m_m1) + std::sqrt(q * q + m_m2 * m_m2);
  }

  double PhaseShiftDensity::qfromM(double M) const {
    // Two-body CM momentum: q = sqrt([M^2-(m1+m2)^2][M^2-(m1-m2)^2]) / (2M).
    double s = M * M;
    double sp = (m_m1 + m_m2) * (m_m1 + m_m2);
    double sm = (m_m1 - m_m2) * (m_m1 - m_m2);
    double arg = (s - sp) * (s - sm);
    return (arg > 0.) ? std::sqrt(arg) / (2. * M) : 0.;
  }

  double PhaseShiftDensity::SpectralWeight(double q) const {
    // (1/pi) sum_l (2J_l+1) ddelta_l/dq.
    // If a wave provides an analytic d(delta)/dM, use it via ddelta/dq =
    // (ddelta/dM)(dM/dq). Otherwise fall back to a central finite difference in q
    // (delta(M(q)) is smooth for q > 0; the q=0 kink is never sampled by the
    // interior Gauss-Legendre nodes). The FD step is clamped to stay at q > 0.
    const double M = Mfromq(q);
    const double dMdq = q / std::sqrt(q * q + m_m1 * m_m1)
                      + q / std::sqrt(q * q + m_m2 * m_m2);

    bool needFD = false;
    for (size_t i = 0; i < m_waves.size(); ++i)
      if (!m_waves[i].ddelta_dM) { needFD = true; break; }
    const double h = std::min(1.0e-5, 0.25 * q);
    const double Mp = needFD ? Mfromq(q + h) : 0.;
    const double Mm = needFD ? Mfromq(q - h) : 0.;

    double w = 0.;
    for (size_t i = 0; i < m_waves.size(); ++i) {
      double ddelta_dq;
      if (m_waves[i].ddelta_dM)
        ddelta_dq = m_waves[i].ddelta_dM(M) * dMdq;                          // analytic
      else
        ddelta_dq = (m_waves[i].delta(Mp) - m_waves[i].delta(Mm)) / (2. * h); // finite diff
      w += m_waves[i].twoJplus1 * ddelta_dq;
    }
    return w / xMath::Pi();
  }

  double PhaseShiftDensity::Integrate(IdealGasFunctions::Quantity quantity, int statistics, double T, double mu) const {
    if (!m_enabled) return 0.;   // disabled channel contributes nothing
    std::vector<double> xleg, wleg;
    NumericalIntegration::GetCoefsIntegrateLegendre(m_nodes, 0., m_qmax, &xleg, &wleg);

    double ret = 0.;
    for (int i = 0; i < m_nodes; ++i) {
      double q = xleg[i];
      double M = Mfromq(q);
      double Xideal = IdealGasFunctions::IdealGasQuantity(
        quantity, IdealGasFunctions::Quadratures, statistics, T, mu, M, /*deg*/ 1., /*order*/ 1);
      ret += wleg[i] * SpectralWeight(q) * Xideal;
    }
    return ret;
  }

  double PhaseShiftDensity::Quantity(IdealGasFunctions::Quantity quantity, double T, double mu) {
    return Integrate(quantity, m_stat, T, mu);
  }

  double PhaseShiftDensity::QuantityCluster(int n, IdealGasFunctions::Quantity quantity, double T, double mu) {
    if (n < 1) return 0.;
    // n-th term of the Boltzmann cluster (fugacity) expansion: the spectral
    // integral with Boltzmann statistics at temperature T/n. The mass-M ideal
    // term n_id(T/n, mu) carries fugacity exp(n mu/T), so the sum over n
    // reconstructs the quantum (m_stat) Quantity() for the fugacity-linear
    // quantities the canonical ensemble requests (ParticleDensity, Pressure,
    // EnergyDensity). For other quantities (entropy, susceptibilities, T
    // derivatives) the T/n term is NOT the correct cluster term, so we keep only
    // the leading n==1 Boltzmann term - matching the base-class semantics - rather
    // than add wrong higher-order terms. The cluster sign (fermionic clusters) is
    // applied by ThermalParticle::DensityCluster.
    if (n > 1 && quantity != IdealGasFunctions::ParticleDensity
              && quantity != IdealGasFunctions::Pressure
              && quantity != IdealGasFunctions::EnergyDensity)
      return 0.;
    return Integrate(quantity, /*Boltzmann*/ 0, T / static_cast<double>(n), mu);
  }

} // namespace thermalfist
