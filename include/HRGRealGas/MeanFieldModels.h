/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2022 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef MEANFIELDMODELS_H
#define MEANFIELDMODELS_H

#include <cmath>
#include <vector>

namespace thermalfist {

  /** \file MeanFieldModel.h
      \brief Header with helper mean-field class

  *   A header with classes that contain implementation of auxiliary functions
  *   describing density-dependent mean fields.
  */

  /**
  *   \brief Base class implementing the ideal gas.
  */
  class MeanFieldModelBase
  {
  public:
    MeanFieldModelBase()
    {
    }
    virtual ~MeanFieldModelBase() { }

    // Mean field in units of GeV/fm^3
    virtual double v(double n) const { return 0.; }

    // Density derivatives of the mean field in units of GeV * fm^{3 * (order-1)}
    virtual double dv(int order, double n) const { return 0.; }

    // Temperature derivative of the mean field in units of fm^{-3}
    virtual double dvdT(double n) const { return 0.; }
  };

  /**
  *   \brief Class implementing auxiliary functions for the linear (van der Waals) mean-field model.
  */
  class MeanFieldModelVDW : public MeanFieldModelBase
  {
  public:
    MeanFieldModelVDW(double a, double dadT = 0.) : m_a(a), m_dadT(dadT)
    {
    }
    virtual ~MeanFieldModelVDW() { }
    virtual double v(double n) const { return -m_a * n * n; }
    virtual double dv(int order, double n) const {
      if (order == 0)
        return v(n);
      if (order == 1)
        return -2. * m_a * n;
      if (order == 2)
        return -2. * m_a;
      return 0.;
    }
    virtual double dvdT(double n) const { return -m_dadT * n * n; }
  private:
    double m_a, m_dadT;
  };

  /**
  *   \brief Class implementing auxiliary functions for the Clausius mean-field model (https://arxiv.org/pdf/1701.06524.pdf).
  */
  class MeanFieldModelClausius : public MeanFieldModelBase
  {
  public:
    MeanFieldModelClausius(double a, double c, double dadT = 0., double dcdT = 0.) : 
      m_a(a), m_c(c), m_dadT(dadT), m_dcdT(dcdT)
    {
    }
    virtual ~MeanFieldModelClausius() { }
    virtual double v(double n) const { return -m_a * n * n / (1. + m_c * n); }
    virtual double dv(int order, double n) const {
      if (order == 0)
        return v(n);
      if (order == 1)
        return - m_a * n * (2. + m_c * n) / pow(1. + m_c * n, 2);

      double mult = 1.;
      if (order & 1)
        mult = -1.;

      //double ret = m_c / (1. + m_c * n);
      double ret = 1. / (1. + m_c * n);
      double tret = ret;
      double k = 1.;
      for (int i = 1; i <= order; ++i) {
        k *= i;
        ret *= tret * k;
        if (i == 2)
          tret *= m_c;
      }

      ret *= -mult * m_a;

      return ret;
    }
    virtual double dvdT(double n) const { 
      return -m_dadT * n * n / (1. + m_c * n) 
        + m_a * n * n / (1. + m_c * n) / (1. + m_c * n) * m_dcdT * n;
    }
  private:
    double m_a, m_c;
    double m_dadT, m_dcdT;
  };

  /**
  *   \brief Class implementing auxiliary functions for the Clausius mean-field model (https://arxiv.org/pdf/1701.06524.pdf).
  */
  class MeanFieldModelSkyrme : public MeanFieldModelBase
  {
  public:
    MeanFieldModelSkyrme(
      double alpha, 
      double beta, 
      double gam = 2.0,
      double n0 = 0.16,
      double dalphadT = 0.,
      double dbetadT = 0.
    ) : m_alpha(alpha), m_beta(beta), m_gam(gam), m_n0(n0), m_dalphadT(dalphadT), m_dbetadT(dbetadT)
    {
    }

    virtual ~MeanFieldModelSkyrme() { }
    virtual double v(double n) const { 
      return m_alpha * n * n / m_n0 
        + m_beta * n * pow(n/m_n0, m_gam); 
    }
    virtual double dv(int order, double n) const {
      if (order == 0)
        return v(n);
      if (order == 1)
        return 2. * m_alpha * n / m_n0 + (1. + m_gam) * m_beta * pow(n/m_n0, m_gam);
      if (order == 2)
        return 2. * m_alpha / m_n0 + (1. + m_gam) * m_gam * m_beta * pow(n, m_gam - 1.) / pow(m_n0, m_gam);

      if (order >= 3) {
        double gam_mult = 1.0;
        for(int i = 0; i < order; ++i) {
          gam_mult *= (1. + m_gam - i);
        }
        return gam_mult * m_beta * pow(n, m_gam + 1 - order) / pow(m_n0, m_gam);
      }

      // if (order == 3)
      //   return 6. * m_beta / m_n0 / m_n0;
      return 0.;
    }
    virtual double dvdT(double n) const {
      return m_dalphadT * n * n / m_n0 
        + m_dbetadT * n * pow(n/m_n0, m_gam); 
      // return m_dalphadT * n * n / m_n0
      //   + m_dbetadT * n * n * n / m_n0 / m_n0;
    }
  private:
    double m_alpha, m_beta, m_gam, m_n0;
    double m_dalphadT, m_dbetadT;
  };

  /**
  *   \brief Class implementing auxiliary functions for the vector density functional (VDF) model from https://arxiv.org/pdf/2011.06635.
  */
  class MeanFieldModelVDF : public MeanFieldModelBase
  {
  public:
    MeanFieldModelVDF(
      int N,
      std::vector<double> Ck,
      std::vector<double> bk,
      double n0 = 0.,
      std::vector<double> dCkdT = std::vector<double>(),
      std::vector<double> dbkdT = std::vector<double>()
    ) : m_N(N), m_Ck(Ck), m_bk(bk), m_dCkdT(dCkdT), m_dbkdT(dbkdT), m_n0(n0)
    {
      if (dCkdT.size() == 0)
        dCkdT = std::vector<double>(Ck.size(), 0.);
      if (dbkdT.size() == 0)
        dbkdT = std::vector<double>(bk.size(), 0.);
    }

    virtual ~MeanFieldModelVDF() { }
    virtual double v(double n) const;
    virtual double dv(int order, double n) const;
    virtual double dvdT(double n) const;
  private:
    int m_N;
    std::vector<double> m_Ck, m_bk;
    std::vector<double> m_dCkdT, m_dbkdT;
    double m_n0;
  };

} // namespace thermalfist

#endif

