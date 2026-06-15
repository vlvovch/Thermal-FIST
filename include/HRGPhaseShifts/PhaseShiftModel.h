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
#include <vector>
#include <map>
#include <utility>

#include "HRGPhaseShifts/PhaseShiftDensity.h"

namespace thermalfist {

  namespace PhaseShifts {

    /**
     * \brief One partial wave of a named analytic parametrization.
     *
     * Known names:
     *   - "pipi_I2:GarciaMartin2011" : Garcia-Martin et al. (2011), I=2 (S, D).
     *
     * \param name       Registry key (channel:parametrization).
     * \param twoJplus1  Spin degeneracy 2J+1 of the wanted wave (e.g. 1=S, 5=D).
     * \param params     Optional named parameter overrides (ignored otherwise).
     * \throws std::invalid_argument if the name is unknown or has no such wave.
     */
    PhaseShiftPartialWave AnalyticWave(const std::string& name, int twoJplus1,
                                       const std::map<std::string, double>& params = std::map<std::string, double>());

    /// All partial waves of a named analytic parametrization (whole-channel
    /// convenience, e.g. the full S+D set of "pipi_I2:GarciaMartin2011").
    std::vector<PhaseShiftPartialWave> AnalyticWaves(const std::string& name,
                                       const std::map<std::string, double>& params = std::map<std::string, double>());

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

    /// Several tabulated waves at once (one per (2J+1, file) entry).
    std::vector<PhaseShiftPartialWave> TabulatedWaves(
                                       const std::vector<std::pair<int, std::string> >& waveFiles);

  } // namespace PhaseShifts

} // namespace thermalfist

#endif
