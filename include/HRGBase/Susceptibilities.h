/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2024 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef SUSCEPTIBILITIES_H
#define SUSCEPTIBILITIES_H

/**
 * \file Susceptibilities.h
 * 
 * \brief Contains some helper functions for calculating (mixed) susceptibilities.
 * 
 */

#include <string>
#include <vector>
#include <map>
#include "ThermalModelBase.h"

namespace thermalfist {

  /// Calculates mixed susceptibilities of Baryon (B), Charge (Q), and Strangeness (S) charges for the given model up to the given order (fourth-order max)
  /// The susceptibilities of third and fourth order are calculated using the central finite difference method with step size dmuTnum
  std::map<std::vector<int>, double> ComputeBQSSusceptibilities(ThermalModelBase *model, int order = 4, double dmuTnum = 0.001);

  /// Calculates temperature derivatives of the mixed susceptibilities of Baryon (B), Charge (Q), and Strangeness (S) charges for the given model up to the given order (fourth-order max)
  /// The susceptibilities of third and fourth order are calculated using the central finite difference method with step size dmuTnum
  std::map<std::vector<int>, double> ComputeBQSSusceptibilitiesDerivativeT(ThermalModelBase *model, int order = 4, double dmuTnum = 0.001);

} // namespace thermalfist

#endif