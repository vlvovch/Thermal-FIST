/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2016-2022 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EXCLUDEDVOLUMEMODELS_H
#define EXCLUDEDVOLUMEMODELS_H

#include <cmath>

namespace thermalfist {

  /** \file ExcludedVolumeModel.h
      \brief Header with helper excluded-volume class

  *   A header with classes that contain implementation of auxiliary excluded-volume functions
  *   needed for real gas interacting HRG implementation.
  */

  /**
  *   \brief Base class implementing the ideal gas.
  */
  class ExcludedVolumeModelBase
  {
  public:
    ExcludedVolumeModelBase()
    {
    }
    virtual ~ExcludedVolumeModelBase() { }
    virtual double EtaMax() const { return -1.; }
    virtual double f(double eta) const { return 1.; }
    virtual double df(int n, double eta) const { return (n == 0) ? 1. : 0.; }
    virtual double etasol(double etatil) const { return etatil; }
    virtual double etasolBinarySearch(double etatil) const;
    //virtual double etasolBroyden(double etatil) const { return etatil; } // TODO: later
  };

  /**
  *   \brief Class implementing auxiliary functions for the van der Waals excluded volume model.
  */
  class ExcludedVolumeModelVDW : public ExcludedVolumeModelBase
  {
  public:
    ExcludedVolumeModelVDW()
    {
    }
    virtual ~ExcludedVolumeModelVDW() { }
    virtual double EtaMax() const { return 0.25; }
    virtual double f(double eta) const { return 1. - 4 * eta; }
    virtual double df(int n, double eta) const { 
      if (n == 0)
        return f(eta);
      if (n == 1)
        return -4.;
      return 0.; 
    }
    virtual double etasol(double etatil) const { return etatil / (1. + 4. * etatil); }
  };

  /**
  *   \brief Class implementing auxiliary functions for the Carnahan-Starling excluded volume model.
  */
  class ExcludedVolumeModelCS : public ExcludedVolumeModelBase
  {
  public:
    ExcludedVolumeModelCS() : ExcludedVolumeModelBase()
    {
    }
    virtual ~ExcludedVolumeModelCS() { }
    virtual double EtaMax() const { return 1.; }
    virtual double f(double eta) const {
      return exp(-(4. - 3. * eta) * eta / (1. - eta) / (1. - eta));
    }
    double d1f(double eta) const {
      return -f(eta) * 2. * (2. - eta) / pow(1. - eta, 3);
    }
    double d2f(double eta) const {
      return f(eta) * 2. * (3. + eta * (4. + eta * (-7. + 2. * eta))) / pow(1. - eta, 6);
    }
    double d3f(double eta) const {
      return f(eta) * 4. * (5. + eta * (-24. + eta * (12. + eta * (17. + 3. * eta * (-5. + eta))))) / pow(1. - eta, 9);
    }
    double d4f(double eta) const {
      double tret = 2. * eta;
      tret += -13.; tret *= eta;
      tret += 22.; tret *= 6. * eta;
      tret += 85.; tret *= eta;
      tret += -440.; tret *= eta;
      tret += 396.; tret *= eta;
      tret += -104.; tret *= eta;
      tret += 1.; tret *= 4.;
      return f(eta) * tret / pow((1. - eta), 12);
    }
    virtual double df(int n, double eta) const;

    virtual double etasol(double etatil) const { return etasolBinarySearch(etatil); }

  };


  /**
  *   \brief Class implementing auxiliary functions for the VDW excluded volume model truncated at n^2.
  */
  class ExcludedVolumeModelVirial : public ExcludedVolumeModelBase
  {
  public:
    ExcludedVolumeModelVirial() : ExcludedVolumeModelBase()
    {
    }
    virtual ~ExcludedVolumeModelVirial() { }
    virtual double EtaMax() const { return 1.e9; }
    virtual double f(double eta) const {
      return exp(-4. * eta);
    }virtual double df(int n, double eta) const {
      if (n == 0)
        return f(eta);

      double mult = 1.;
      if (n & 1)
        mult = -1.;
      return mult * pow(4., n) * f(eta);
    }

    virtual double etasol(double etatil) const { return etasolBinarySearch(etatil); }
  };

  /**
  *   \brief Class implementing auxiliary functions for the VDW excluded volume model truncated at n^3. 
  *   This corresponds to the trivial model (TVM) from https://arxiv.org/abs/1909.02276
  */
  class ExcludedVolumeModelTVM : public ExcludedVolumeModelBase
  {
  public:
    ExcludedVolumeModelTVM() : ExcludedVolumeModelBase()
    {
    }
    virtual ~ExcludedVolumeModelTVM() { }
    virtual double EtaMax() const { return 1.e9; }
    virtual double f(double eta) const {
      return exp(-4. * eta + 8. * eta * eta);
    }
    double d1f(double eta) const {
      return -f(eta) * 4. * (1. + 4. * eta);
    }
    double d2f(double eta) const {
      return f(eta) * 128. * eta * (1. + 2. * eta);
    }
    double d3f(double eta) const {
      return -f(eta) * 128. * eta * (-1. + 8. * eta * eta * (3. + 4. * eta));
    }
    double d4f(double eta) const {
      return f(eta) * 512. * eta * (-1. + 16. * eta * (-1. + 8. * eta * eta * (1. + eta)));
    }
    virtual double df(int n, double eta) const;

    virtual double etasol(double etatil) const { return etasolBinarySearch(etatil); }
  };

} // namespace thermalfist

#endif

