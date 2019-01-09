/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <vector>
#include <string>

#include "HRGBase/BilinearSplineFunction.h"

namespace thermalfist {

  /// \brief Routines used for modeling the acceptance effects.
  namespace Acceptance {

    /**
     *  \brief Structure which contains the binomial probabilities
     *         for particle with given y and pt to be accepted.
     * 
     *  Assumes that all bins have the same width in rapidity and pT.
     *  Acceptance function is \f$p(y,p_T)\f$ where p is the acceptance probability.
     */
    struct AcceptanceFunction {
      bool init;

      double dy;  ///< Rapidity width of a bin
      double dpt; ///< pT width of a bin

      std::vector<double> ys;       ///< Vector of bin rapidities. One element per bin.
      std::vector<double> pts;      ///< Vector of bin pT values. One element per bin.
      std::vector<double> probs;    ///< Vector of acceptance probabilities for each bin.
      BilinearSplineFunction sfunc; ///< 2D spline interpolation of the acceptance function

      AcceptanceFunction() : dy(), dpt(), ys(), pts(), probs(), sfunc() { init = false; }

      void setSpline() { sfunc.setData(ys, pts, probs); init = true; }

      /// Binomial acceptance for the given values of y and pt
      double getAcceptance(const double & y, const double & pt) const;
    };



    /**
     *  Read the acceptance function from file.
     */
    int ReadAcceptanceFunction(AcceptanceFunction & func, std::string filename);
  }

} // namespace thermalfist

#endif
