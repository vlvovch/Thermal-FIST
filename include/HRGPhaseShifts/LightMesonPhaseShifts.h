/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef LIGHTMESONPHASESHIFTS_H
#define LIGHTMESONPHASESHIFTS_H

/**
 * \file LightMesonPhaseShifts.h
 *
 * \brief Catalog of empirical elastic phase-shift parametrizations for light
 *        meson-meson scattering, for use with PhaseShiftDensity (S-matrix HRG).
 *
 * Currently implemented:
 *   - pi-pi, I=0 (attractive, the sigma/f0(500)): S-wave (J=0), conformal CFD
 *     parametrization of Garcia-Martin et al., Phys.Rev. D83 (2011) 074004
 *     [arXiv:1102.2183], Appendix A.1 (Eqs. A1-A3, Table V). The phase passes
 *     through 90 deg (sigma), so it is branch-tracked; elastic below the K-Kbar
 *     threshold.
 *   - pi-pi, I=2 (repulsive): S-wave (J=0) and D-wave (J=2),
 *     conformal CFD parametrization of Garcia-Martin et al.,
 *     Phys.Rev. D83 (2011) 074004 [arXiv:1102.2183], Tables VI, IX.
 *   - pi-K, I=3/2 (repulsive) and I=1/2 (attractive, the kappa/K0*(700)):
 *     S-wave (J=0), conformal CFD parametrization of Pelaez, Rodas,
 *     Phys.Rev. D93 (2016) 074025 [arXiv:1602.08404], Eqs. (11), (13),
 *     Tables I, II. The I=1/2 wave is elastic only below the K-eta threshold.
 */

#include <vector>

#include "HRGPhaseShifts/PhaseShiftDensity.h"

namespace thermalfist {

  namespace PhaseShifts {

    /// Charged-pion mass [GeV] used by the pi-pi parametrizations.
    constexpr double PionMass() { return 0.13957; }
    /// Upper validity limit sqrt(s)_max [GeV] of the pi-pi I=2 parametrization.
    constexpr double PiPiI2_Mmax() { return 1.420; }
    /// pi-pi I=0 S-wave is treated as elastic up to the K-Kbar threshold (2 M_K),
    /// which is the physical limit of the phase-shift-only Beth-Uhlenbeck term and
    /// covers the sigma/f0(500). Raise to 1.420 for the full Garcia-Martin range
    /// (then the f0(980) is included and must not be double-counted in the list).
    constexpr double PiPiI0_Mmax() { return 2. * 0.496; }   // 2 M_K

    /// pi-pi I=0 S-wave phase shift delta_0^{(0)}(M) [radians]. Attractive (> 0),
    /// branch-tracked through 90 deg (sigma); full Garcia-Martin A1-A3 form.
    double PiPi_delta_I0_S(double M);

    /// The part of delta_0^{(0)} ABOVE the K-Kbar threshold (the f0(980) region),
    /// offset to 0 at 2 M_K so it carries only the high-mass spectral weight (the
    /// sigma part below 2 M_K is the separate PiPi_delta_I0_S wave). [radians]
    double PiPi_delta_I0_f0980_S(double M);

    /// S-wave (2J+1=1) of the pi-pi I=0 channel (the sigma). d(delta)/dM uses the
    /// engine's finite-difference fallback.
    std::vector<PhaseShiftPartialWave> PiPi_I0_Waves();

    /// S-wave (2J+1=1) of the pi-pi I=0 f0(980) region (above the K-Kbar threshold).
    std::vector<PhaseShiftPartialWave> PiPi_I0_f0980_Waves();

    /// pi-pi I=2 S-wave phase shift delta_0^{(2)}(M) [radians]. Repulsive (< 0).
    double PiPi_delta_I2_S(double M);
    /// pi-pi I=2 D-wave phase shift delta_2^{(2)}(M) [radians]. Repulsive (< 0).
    double PiPi_delta_I2_D(double M);

    /// Analytic d(delta)/dM [rad/GeV] for the pi-pi I=2 S- and D-waves.
    double PiPi_ddelta_I2_S_dM(double M);
    double PiPi_ddelta_I2_D_dM(double M);

    /// Partial waves of the pi-pi I=2 channel: S (2J+1=1) and D (2J+1=5),
    /// each carrying its exact analytic d(delta)/dM.
    std::vector<PhaseShiftPartialWave> PiPi_I2_Waves();

    /// Isospin-averaged kaon mass [GeV] used by the pi-K parametrizations.
    constexpr double KaonMass() { return 0.496; }
    /// Eta mass [GeV]; sets the K-eta inelastic threshold for the pi-K I=1/2 wave.
    constexpr double EtaMass()  { return 0.5478; }
    /// Upper validity limit sqrt(s)_max [GeV] of the pi-K I=3/2 S-wave (elastic).
    constexpr double PiK_I32_Mmax() { return 1.740; }
    /// pi-K I=1/2 S-wave is elastic only below the K-eta threshold m_K + m_eta.
    constexpr double PiK_I12_Mmax() { return KaonMass() + EtaMass(); }

    /// pi-K I=3/2 S-wave phase shift delta_0^{3/2}(M) [radians]. Repulsive (< 0).
    double PiK_delta_I32_S(double M);
    /// pi-K I=1/2 S-wave phase shift delta_0^{1/2}(M) [radians]. Attractive (> 0),
    /// the kappa/K0*(700); returns 0 at/above the K-eta threshold.
    double PiK_delta_I12_S(double M);

    /// S-wave (2J+1=1) of the pi-K I=3/2 (repulsive) and I=1/2 (attractive)
    /// channels. d(delta)/dM uses the engine's finite-difference fallback.
    std::vector<PhaseShiftPartialWave> PiK_I32_Waves();
    std::vector<PhaseShiftPartialWave> PiK_I12_Waves();

  } // namespace PhaseShifts

} // namespace thermalfist

#endif
