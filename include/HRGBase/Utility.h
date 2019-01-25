/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef UTILITY_H
#define UTILITY_H

/**
 * \file Utility.h
 * 
 * \brief Contains some helper functions.
 * 
 */

#include <string>

namespace thermalfist {

  class Disclaimer {

  public:
    static bool PrintDisclaimer();

    static bool DisclaimerPrinted;
  
  };


  // For C99 compatibility
  long long stringToLongLong(const std::string &str);

  double get_wall_time();

  double get_cpu_time();

} // namespace thermalfist

#endif