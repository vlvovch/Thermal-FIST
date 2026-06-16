/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef MESONBARYONPHASESHIFTS_H
#define MESONBARYONPHASESHIFTS_H

/**
 * \file MesonBaryonPhaseShifts.h
 *
 * \brief Catalog of empirical elastic phase-shift parametrizations for
 *        meson-baryon scattering (the baryon-sector companion of
 *        LightMesonPhaseShifts.h), for use with PhaseShiftDensity (S-matrix HRG).
 *
 * Currently implemented:
 *   pi-N Roy-Steiner waves of Hoferichter, Ruiz de Elvira, Kubis, Meissner,
 *   Phys.Rept. 625 (2016) 1 [arXiv:1510.06039], Appendix D (Tables D.1-D.2):
 *   - P33 = the Delta(1232): conformal parametrization (Eq. D.3), resonant
 *     (branch-tracked through 90 deg), elastic across the resonance up to the
 *     matching point sqrt(s) = 1.38 GeV.
 *   - S31, S11: S-wave Schenk parametrization (Eq. D.1). S31 repulsive, S11
 *     attractive.
 *   - P31, P11, P13: P-wave Schenk parametrization (Eq. D.2), non-resonant.
 *   S31 and P13 are elastic up to the matching point 1.38 GeV; S11, P31 and P11
 *   are elastic only below the pi-pi-N threshold m_N + 2 m_pi ~ 1.22 GeV (above
 *   it the inelasticity is non-negligible), so their integration is cut there.
 *
 * See input/list/phaseshifts/piN_RoySteiner_reference.md for the full coefficient
 * tables and a numerical validation of all waves.
 */

#include <vector>

#include "HRGPhaseShifts/PhaseShiftDensity.h"
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"   // PionMass()

namespace thermalfist {

  namespace PhaseShifts {

    /// Nucleon mass [GeV] used by the pi-N parametrizations (reproduces the
    /// published Roy-Steiner P33 curve, e.g. delta = 90 deg at the Delta).
    constexpr double NucleonMass() { return 0.938272; }

    /// Roy-Steiner matching point sqrt(s)_max [GeV] (the validity limit). The
    /// S31, P13 and P33(Delta) waves are elastic up to here.
    constexpr double PiN_Mmax_elastic() { return 1.380; }
    /// pi-pi-N inelastic threshold m_N + 2 m_pi [GeV] (~1.217); the S11, P31 and
    /// P11 waves are elastic only below it, so their integration is cut here.
    constexpr double PiN_Mmax_inelastic() { return NucleonMass() + 2. * PionMass(); }
    /// Upper validity limit for the Delta(1232) P33 wave (the matching point).
    constexpr double PiN_Delta_Mmax() { return PiN_Mmax_elastic(); }

    /// pi-N I=3/2 P-wave phase shift delta_1+^{3/2}(M) [radians] - the Delta(1232).
    /// Resonant (branch-tracked through 90 deg at the Delta), conformal Roy-Steiner
    /// parametrization (Eq. D.3 of arXiv:1510.06039).
    double PiN_delta_P33(double M);

    /// pi-N S-wave phase shifts [radians], Schenk parametrization (Eq. D.1).
    double PiN_delta_S31(double M);   ///< I=3/2 (repulsive, < 0)
    double PiN_delta_S11(double M);   ///< I=1/2 (attractive, > 0)
    /// pi-N P-wave phase shifts [radians], Schenk parametrization (Eq. D.2);
    /// non-resonant (the Delta is the separate conformal P33 wave).
    double PiN_delta_P31(double M);   ///< I=3/2, J=1/2
    double PiN_delta_P11(double M);   ///< I=1/2, J=1/2 (the Roper N(1440) tail)
    double PiN_delta_P13(double M);   ///< I=1/2, J=3/2

    /// P-wave (2J+1=4, J=3/2) of the pi-N I=3/2 channel - the Delta(1232).
    std::vector<PhaseShiftPartialWave> PiN_Delta_Waves();
    /// Single-wave sets for the non-resonant pi-N background waves. S/P J=1/2
    /// waves have 2J+1=2; the P13 J=3/2 wave has 2J+1=4.
    std::vector<PhaseShiftPartialWave> PiN_S31_Waves();
    std::vector<PhaseShiftPartialWave> PiN_S11_Waves();
    std::vector<PhaseShiftPartialWave> PiN_P31_Waves();
    std::vector<PhaseShiftPartialWave> PiN_P11_Waves();
    std::vector<PhaseShiftPartialWave> PiN_P13_Waves();

  } // namespace PhaseShifts

} // namespace thermalfist

#endif
