/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019-2025 Volodymyr Vovchenko
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
#include <map>

namespace thermalfist {

  class Disclaimer {

  public:
    static bool PrintDisclaimer();

    static bool DisclaimerPrinted;

  };

  // Read parameters from a file
  /**
   *  Records the parameters from a file into a map.
   * 
   *  The file should contain the parameters in the format:
   *  
   *  # comment
   *  var1 val1  # comment
   *  var2 val2
   *  ...
   * 
   *  where var1, var2, ... are the parameter names and val1, val2, ... are the parameter values.
   * 
   *  The function returns the map with the parameters.
   *  It overrides the parameters in @params with the parameters from the file.
   * 
   *  @param filename The name of the file with the parameters.
   *  @param params The map with the existing/default parameters.
   * 
   */
  std::map<std::string, std::string> ReadParametersFromFile(const std::string& filename, const std::map<std::string, std::string>& params = {});

  /**
   *  Records the parameters from command line arguments into a map.
   * 
   *  --var1 val1 --var2 val2 ...
   * 
   *  where var1, var2, ... are the parameter names and val1, val2, ... are the parameter values.
   * 
   *  The function writes the parameters to the map.
   * 
   *  @param argc The number of command line arguments.
   *  @param argv The array of command line arguments.
   *  @param params The map with the existing/default parameters.
   * 
   */
  void ParametersFromArgs(int argc, char* argv[], std::map<std::string, std::string>& params);

  // For C99 compatibility
  long long stringToLongLong(const std::string &str);

  double get_wall_time();

  double get_cpu_time();

} // namespace thermalfist

#endif