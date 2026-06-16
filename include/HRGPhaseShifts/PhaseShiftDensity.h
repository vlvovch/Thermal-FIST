/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef PHASESHIFTDENSITY_H
#define PHASESHIFTDENSITY_H

/**
 * \file PhaseShiftDensity.h
 *
 * \brief S-matrix / Beth-Uhlenbeck effective degree of freedom: a two-body
 *        scattering channel whose thermal contribution is weighted by the
 *        energy derivative of the phase shift. The contribution can be NEGATIVE
 *        for repulsive interactions (ddelta/dM < 0).
 */

#include <vector>
#include <functional>

#include "HRGBase/IdealGasFunctions.h"

namespace thermalfist {

  /**
   * \brief One partial wave of a scattering channel: its degeneracy (2J+1) and
   *        the elastic phase shift delta(M) [radians] as a function of the
   *        two-body invariant mass M [GeV].
   */
  struct PhaseShiftPartialWave {
    int twoJplus1;                       ///< (2J+1) partial-wave/spin degeneracy
    std::function<double(double)> delta; ///< delta(M) in radians, M in GeV
    /// Optional analytic d(delta)/dM [rad/GeV]. If empty, the engine falls back
    /// to a central finite difference of delta() in q (the default).
    std::function<double(double)> ddelta_dM;
    PhaseShiftPartialWave() : twoJplus1(1) {}
    PhaseShiftPartialWave(int g, std::function<double(double)> d,
                          std::function<double(double)> dd = std::function<double(double)>())
      : twoJplus1(g), delta(std::move(d)), ddelta_dM(std::move(dd)) {}
  };

  /**
   * \brief GeneralizedDensity implementing the S-matrix / Beth-Uhlenbeck
   *        contribution of a two-body scattering channel.
   *
   * The effective spectral weight is
   *   \f$ \rho(M)\,dM = \frac{1}{\pi}\sum_l (2J_l+1)\,\frac{d\delta_l}{dM}\,dM \f$,
   * and any thermodynamic quantity of the channel is
   *   \f$ X = \int dM\,\rho(M)\, X_\mathrm{ideal}(T,\mu,M) \f$.
   *
   * The integral is performed in the centre-of-mass momentum \f$q\f$ (where the
   * integrand is regular at threshold, the \f$dM/dq\f$ Jacobian cancelling the
   * \f$1/q\f$ singularity of \f$d\delta/dM\f$). The result is NEGATIVE for
   * repulsive channels. \f$d\delta/dq\f$ is evaluated by central finite
   * differences in \f$q\f$ (no external autodiff dependency).
   *
   * Attach to a ThermalParticle via SetGeneralizedDensity(); the particle's
   * quantum numbers (Q, B, S, ...) then route the contribution into all
   * densities and susceptibilities automatically.
   */
  class PhaseShiftDensity : public GeneralizedDensity {
  public:
    /**
     * \param waves       Partial waves (degeneracy + delta(M)) of the channel.
     * \param m1,m2       Masses [GeV] of the two scattering constituents.
     * \param Mmax        Upper invariant-mass limit [GeV] of the parametrization
     *                    (the integration cutoff). Must exceed m1+m2.
     * \param statistics  Quantum statistics of the effective cluster
     *                    (0 Boltzmann, +1 Fermi, -1 Bose; default Bose).
     * \param quadratureNodes  Number of Gauss-Legendre nodes for the q-integral.
     */
    PhaseShiftDensity(const std::vector<PhaseShiftPartialWave>& waves,
                      double m1, double m2, double Mmax,
                      int statistics = -1, int quadratureNodes = 64);

    double Quantity(IdealGasFunctions::Quantity quantity, double T, double mu) override;

    /// n-th Boltzmann cluster term of the Beth-Uhlenbeck integral (for the
    /// canonical ensemble): the spectral integral with Boltzmann statistics at
    /// temperature T/n. Summing over n reconstructs the quantum Quantity() for the
    /// fugacity-linear quantities the canonical ensemble sums (ParticleDensity,
    /// Pressure, EnergyDensity); for other quantities (entropy, susceptibilities)
    /// only the n==1 Boltzmann term is kept (see GeneralizedDensity::QuantityCluster).
    /// The cluster sign is applied by the caller.
    double QuantityCluster(int n, IdealGasFunctions::Quantity quantity, double T, double mu) override;

    /// Two-body invariant mass M [GeV] from the CM momentum q [GeV].
    double Mfromq(double q) const;
    /// CM momentum q [GeV] from the invariant mass M [GeV].
    double qfromM(double M) const;

    double Threshold() const { return m_m1 + m_m2; } ///< M threshold = m1 + m2.
    double Mmax() const { return m_Mmax; }
    double qMax() const { return m_qmax; }

    /// Effective spectral weight in q-space: (1/pi) sum_l (2J_l+1) ddelta_l/dq.
    double SpectralWeight(double q) const;

    /// Enable/disable the contribution. When disabled, Quantity() returns 0 for
    /// every quantity, so the channel contributes nothing to densities,
    /// susceptibilities or feeddown without removing it from the particle list.
    void SetEnabled(bool enabled) { m_enabled = enabled; }
    bool IsEnabled() const { return m_enabled; }

    /// Phase-shift channels are managed independently of the thermal-model
    /// configuration, so ClearDensityModels() (which resets the EMM pion models)
    /// must not delete them.
    bool IsClearable() const override { return false; }

  private:
    /// Spectral integral of the ideal-gas quantity with the given quantum
    /// statistics at temperature T (shared by Quantity and QuantityCluster).
    double Integrate(IdealGasFunctions::Quantity quantity, int statistics, double T, double mu) const;

    std::vector<PhaseShiftPartialWave> m_waves;
    double m_m1, m_m2;
    double m_Mmax, m_qmax;
    int    m_stat;
    int    m_nodes;
    bool   m_enabled;
  };

} // namespace thermalfist

#endif
