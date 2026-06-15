/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef PHASESHIFTMODEL_H
#define PHASESHIFTMODEL_H

/**
 * \file PhaseShiftModel.h
 *
 * \brief Selectable phase-shift models, at the granularity of a single partial
 *        wave. The dynamics of one wave is a PhaseShiftPartialWave (2J+1 +
 *        delta(M) [+ analytic ddelta/dM]); a channel's dynamics is just a set of
 *        such per-wave models, so each wave can use a different parametrization
 *        (analytic, or a tabulated phase shift) and coupled-channel models can
 *        later be attached at the same per-wave granularity.
 */

#include <string>

#include "HRGPhaseShifts/PhaseShiftDensity.h"

namespace thermalfist {

  namespace PhaseShifts {

    /**
     * \brief One partial wave of a named analytic parametrization.
     *
     * The model registry is per-wave, and the model NAME is wave-specific: each
     * (channel, param) resolves to exactly one wave's delta(M) (the returned wave
     * carries its own 2J+1). There is no name that stands for "all waves" - two
     * waves of one paper get two names (they have different phase shifts).
     *
     * Known (channel, param):
     *   - ("pipi_I2", "GarciaMartin2011_S") : Garcia-Martin et al. (2011) I=2
     *     S-wave (2J+1=1).
     *   - ("pipi_I2", "GarciaMartin2011_D") : Garcia-Martin et al. (2011) I=2
     *     D-wave (2J+1=5).
     *
     * \param channel    Channel name, e.g. "pipi_I2".
     * \param param      Wave-specific parametrization name, e.g. "GarciaMartin2011_S".
     * \throws std::invalid_argument if the (channel, param) is unknown.
     */
    PhaseShiftPartialWave AnalyticWave(const std::string& channel, const std::string& param);

    /**
     * \brief One tabulated partial wave: delta(M) [radians] read from a
     *        whitespace-separated two-column text file "M[GeV]  delta[rad]"
     *        (lines starting with '#' ignored), natural-cubic-spline interpolated;
     *        ddelta/dM is the spline's analytic derivative.
     *
     * \param twoJplus1  Spin degeneracy 2J+1 of the wave.
     * \param file       Path to the (M, delta) table.
     * \throws std::runtime_error if the file cannot be read or has < 2 points.
     */
    PhaseShiftPartialWave TabulatedWave(int twoJplus1, const std::string& file);

  } // namespace PhaseShifts

} // namespace thermalfist

#endif
