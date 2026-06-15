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
 *   - pi-pi, I=2 (repulsive): S-wave (J=0) and D-wave (J=2),
 *     conformal CFD parametrization of Garcia-Martin et al.,
 *     Phys.Rev. D83 (2011) 074004 [arXiv:1102.2183], Tables VI, IX.
 */

#include <vector>

#include "HRGPhaseShifts/PhaseShiftDensity.h"

namespace thermalfist {

  namespace PhaseShifts {

    /// Charged-pion mass [GeV] used by the pi-pi parametrizations.
    constexpr double PionMass() { return 0.13957; }
    /// Upper validity limit sqrt(s)_max [GeV] of the pi-pi I=2 parametrization.
    constexpr double PiPiI2_Mmax() { return 1.420; }

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

  } // namespace PhaseShifts

} // namespace thermalfist

#endif
