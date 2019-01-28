/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELFITQUANTITIES_H
#define THERMALMODELFITQUANTITIES_H

#include <string>

#include "HRGBase/ParticleDecay.h"

/**
 * \file  ThermalModelFitQuantities.h
 * \brief Contains structures describing the yields and ratios used in thermal fits
 * 
 */

namespace thermalfist {
  
  /**
   * \brief Structure containing the experimental yield
   *        (multiplicity) to be fitted.
   * 
   */
  struct ExperimentMultiplicity {
    /// PDG code of the particle yield
    long long fPDGID;

    /// Experimental value
    double fValue;
    
    /// Experimental error
    double fError;

    /// The feeddown contributions to be included
    Feeddown::Type fFeedDown; 

    /**
     * \brief Construct a new ExperimentMultiplicity object.
     * 
     * \param PDGID \copydoc fPDGID
     * \param value \copydoc fValue
     * \param error \copydoc fError
     * \param fd    \copydoc fFeedDown
     */
    ExperimentMultiplicity(long long PDGID = -211, double value = 300., double error = 20., Feeddown::Type fd = Feeddown::StabilityFlag) :
      fPDGID(PDGID), fValue(value), fError(error), fFeedDown(fd) { }

    /// Adds a relative systematic error as a fraction of the total yield
    void addSystematicError(double percent) {
      fError = sqrt(fError*fError + percent * percent * fValue * fValue);
    }
  };

  /**
   * \brief Structure containing the experimental ratio
   *        of yields to be fitted.
   * 
   */
  struct ExperimentRatio {
    /// PDG code of the particle yield in the numerator
    long long fPDGID1;

    /// PDG code of the particle yield in the denominator
    long long fPDGID2;

    /// Experimental value of the yield ratio
    double fValue;

    /// Experimental error of the yield ratio
    double fError;

    /// The feeddown contributions to be included for the yield in the numerator
    Feeddown::Type fFeedDown1;

    /// The feeddown contributions to be included for the yield in the denominator
    Feeddown::Type fFeedDown2; 

    /**
     * \brief Construct a new ExperimentRatio object.
     * 
     * \param PDGID1 \copydoc fPDGID1
     * \param PDGID2 \copydoc fPDGID2
     * \param value  \copydoc fValue
     * \param error  \copydoc fError
     * \param fd1    \copydoc fFeedDown1
     * \param fd2    \copydoc fFeedDown2
     */
    ExperimentRatio(long long PDGID1 = 211, long long PDGID2 = -211, double value = 1., double error = 0.1, Feeddown::Type fd1 = Feeddown::StabilityFlag, Feeddown::Type fd2 = Feeddown::StabilityFlag) :
      fPDGID1(PDGID1), fPDGID2(PDGID2), fValue(value), fError(error), fFeedDown1(fd1), fFeedDown2(fd2) { }
    
    /**
     * \brief  Construct a new ExperimentRatio object two individual yields and their errors.
     * 
     * The error of the ratio is computed through standard
     * propagation of uncertainties in the individual yields.
     * 
     * \param PDGID1 \copydoc fPDGID1
     * \param PDGID2 \copydoc fPDGID2
     * \param value1 Experimental value of the yield in the numerator
     * \param error1 Experimental error of the yield in the numerator
     * \param value2 Experimental value of the yield in the denominator
     * \param error2 Experimental error of the yield in the denominator
     * \param fd1    \copydoc fFeedDown1
     * \param fd2    \copydoc fFeedDown2
     */
    ExperimentRatio(long long PDGID1, long long PDGID2, double value1, double error1, double value2, double error2, Feeddown::Type fd1 = Feeddown::StabilityFlag, Feeddown::Type fd2 = Feeddown::StabilityFlag) :
      fPDGID1(PDGID1), fPDGID2(PDGID2), fFeedDown1(fd1), fFeedDown2(fd2) {
      fValue = value1 / value2;
      fError = sqrt(error1*error1 / value2 / value2 + value1 * value1 / value2 / value2 / value2 / value2 * error2 * error2);
    }
  };

  /**
   * \brief Structure describing the measurement to be fitted
   *        or compared to model
   * 
   */
  struct FittedQuantity {
    /// Yield (multiplicity) or ratio
    enum FittedQuantityType { 
        Multiplicity = 0, 
        Ratio = 1 
    };

    /// Whether it is a yield (multiplicity) or a ratio
    FittedQuantityType type;

    /// Whether this quantity contributes to the \f$ \chi^2 \f$ of a fit
    bool toFit;

    /// The yield data. Used if type is FittedQuantityType::Multiplicity
    ExperimentMultiplicity mult;

    /// The ratio data. Used if type is FittedQuantityType::Ratio
    ExperimentRatio ratio;

    /// Default constructor
    FittedQuantity() {
      toFit = true;
      type = FittedQuantity::Multiplicity;
      mult = ExperimentMultiplicity(-211, 10., 1.);
    }

    /// Constructs a yield measurement
    /// \param[in] op Yield measurement data
    FittedQuantity(const ExperimentMultiplicity & op) {
      toFit = true;
      type = FittedQuantity::Multiplicity;
      mult = op;
    }

    /// Constructs a yield ratio measurement
    /// \param[in] op Yield ratio measurement data
    FittedQuantity(const ExperimentRatio & op) {
      toFit = true;
      type = FittedQuantity::Ratio;
      ratio = op;
    }

    /// Value of the measurement
    double Value() const {
      if (type == Multiplicity)
        return mult.fValue;
      else
        return ratio.fValue;
    }

    /// Error of the measurement
    double ValueError() const {
      if (type == Multiplicity)
        return mult.fError;
      else
        return ratio.fError;
    }
  };

} // namespace thermalfist

#endif