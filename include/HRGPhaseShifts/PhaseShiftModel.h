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
 * \brief A selectable phase-shift model: the partial waves (2J+1 + delta(M)
 *        [+ analytic ddelta/dM]) of one scattering channel. The model is
 *        decoupled from the channel's physical structure so the same channel
 *        can be evaluated with different dynamics (a named analytic
 *        parametrization, or a tabulated phase shift read from file).
 */

#include <string>
#include <vector>
#include <map>
#include <utility>

#include "HRGPhaseShifts/PhaseShiftDensity.h"

namespace thermalfist {

  namespace PhaseShifts {

    /**
     * \brief A phase-shift model for one channel: its partial waves, each
     *        carrying 2J+1, delta(M) and (when available) an analytic ddelta/dM.
     */
    struct PhaseShiftModel {
      std::string name;                          ///< human-readable model tag
      std::vector<PhaseShiftPartialWave> waves;  ///< partial waves (delta providers)
    };

    /**
     * \brief Build a named analytic phase-shift model from the built-in registry.
     *
     * Known names:
     *   - "pipi_I2:GarciaMartin2011" : Garcia-Martin et al. (2011) I=2 S+D waves.
     *
     * \param name    Registry key (channel:parametrization).
     * \param params  Optional named parameter overrides for parametrizations
     *                that expose tunable parameters (ignored otherwise).
     * \throws std::invalid_argument if the name is not in the registry.
     */
    PhaseShiftModel AnalyticModel(const std::string& name,
                                  const std::map<std::string, double>& params = std::map<std::string, double>());

    /**
     * \brief Build a tabulated phase-shift model.
     *
     * Each entry (twoJplus1, file) supplies one partial wave whose delta(M)
     * [radians] is read from a whitespace-separated two-column text file
     * "M[GeV]  delta[rad]" (lines starting with '#' ignored) and interpolated
     * with a natural cubic spline; ddelta/dM is the spline's analytic derivative.
     *
     * \param name       Human-readable model tag.
     * \param waveFiles  List of (2J+1, file path) pairs, one per partial wave.
     * \throws std::runtime_error if a file cannot be read or has < 2 points.
     */
    PhaseShiftModel TabulatedModel(const std::string& name,
                                   const std::vector<std::pair<int, std::string> >& waveFiles);

  } // namespace PhaseShifts

} // namespace thermalfist

#endif
