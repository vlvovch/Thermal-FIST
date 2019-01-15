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

namespace thermalfist {

  class Disclaimer {

  public:
    static bool PrintDisclaimer();

    static bool DisclaimerPrinted;
  
  };

} // namespace thermalfist

#endif