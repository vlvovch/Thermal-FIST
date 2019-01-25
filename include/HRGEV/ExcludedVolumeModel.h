/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2016-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EXCLUDEDVOLUMEMODEL_H
#define EXCLUDEDVOLUMEMODEL_H

#include <cmath>

namespace thermalfist {

  /** \file ExcludedVolumeModel.h
      \brief Header with helper excluded-volume class

  *   A header with classes that contain implementation of auxiliary excluded-volume functions
  *   needed for multi-component mean-field approach.
  *   Contains ``Diagonal'' van der Waals model (default, base class ExcludedVolumeModel) and
  *   ``Diagonal'' Carnahan-Starling model (class ExcludedVolumeModelCS, overrides base class).
  *   There are not yet fully used in the package.
  */

  /**
  *   \brief Base class implementing auxiliary excluded-volume functions
  *          needed for multi-component mean-field approach.
  *   Contains van der Waals functions.
  */
  class ExcludedVolumeModel
  {
  public:
    ExcludedVolumeModel()
    {
    }
    virtual ~ExcludedVolumeModel() { }
    virtual double xMax() const { return 1.; }
    virtual double f(double x) const { return 1. / (1. - x); }
    virtual double g1(double x) const { return 1. / (1. - x); }
    virtual double g2(double x) const { return -log(1. - x); }
    virtual double g3(double /*x*/) const { return 1.; }
    virtual double Dg3(double /*x*/) const { return 0.; }
  };

  /**
  *   \brief Derived class implementing auxiliary excluded-volume functions
  *          for multi-component mean-field approach from the Carnahan-Starling model.
  */
  class ExcludedVolumeModelCS : public ExcludedVolumeModel
  {
  public:
    ExcludedVolumeModelCS() : ExcludedVolumeModel()
    {
    }
    virtual ~ExcludedVolumeModelCS() { }
    virtual double xMax() const { return 4.; }
    virtual double f(double x) const {
      double tx = x / 4.;
      return (1 + tx + tx * tx - tx * tx*tx) / (1. - tx) / (1. - tx) / (1. - tx);
    }
    virtual double g1(double x) const {
      double tx = x / 4.;
      return (1 - tx / 2.) / (1. - tx) / (1. - tx) / (1. - tx);
    }
    virtual double g2(double x) const {
      double tx = x / 4.;
      return (1 - 3.*x / 16.) * x / (1. - tx) / (1. - tx);
    }
    virtual double g3(double x) const {
      double tx = x / 4.;
      return  (1 - tx / 2.) / (1 + tx + tx * tx - tx * tx*tx);
    }
    virtual double Dg3(double x) const {
      double tx = x / 4.;
      double zn = (1 + tx + tx * tx - tx * tx*tx);
      return -1. / 8. / zn - (1. - tx / 2.) / zn / zn * (1. / 4. + 2. * tx - 3.*tx*tx);
    }
  };

} // namespace thermalfist

#endif

