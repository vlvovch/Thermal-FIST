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

  /** \file MeanFieldModels.h
      \brief Header with helper mean-field classes

  *   A header with classes that contain implementation of auxiliary functions
  *   describing density-dependent mean fields.
  */

  /**
  *   \brief Base class implementing the ideal gas (no mean field).
  *
  *   This class serves as the base for all mean field models and implements
  *   the trivial case of an ideal gas with no mean field.
  */
  class MeanFieldModelBase
  {
  public:
    /**
     * \brief Constructor for the MeanFieldModelBase class.
     */
    MeanFieldModelBase()
    {
    }
    
    /**
     * \brief Destructor for the MeanFieldModelBase class.
     */
    virtual ~MeanFieldModelBase() { }

    /**
     * \brief Calculates the mean field value at a given density.
     * 
     * \param n Particle number density in fm^-3.
     * \return Mean field value in units of GeV/fm^3.
     */
    virtual double v(double n) const { return 0.; }

    /**
     * \brief Calculates the density derivatives of the mean field.
     * 
     * \param order Order of the derivative.
     * \param n Particle number density in fm^-3.
     * \return Derivative of the mean field in units of GeV * fm^{3 * (order-1)}.
     */
    virtual double dv(int order, double n) const { return 0.; }

    /**
     * \brief Calculates the temperature derivative of the mean field.
     * 
     * \param n Particle number density in fm^-3.
     * \return Temperature derivative of the mean field in units of fm^{-3}.
     */
    virtual double dvdT(double n) const { return 0.; }
  };

  /**
  *   \brief Class implementing auxiliary functions for the linear (van der Waals) mean-field model.
  *
  *   This class implements the linear mean-field model, also known as the van der Waals model,
  *   where the mean field is proportional to the particle number density.
  */
  class MeanFieldModelVDW : public MeanFieldModelBase
  {
  public:
    /**
     * \brief Constructor for the MeanFieldModelVDW class.
     * 
     * \param a Attraction parameter in units of GeV*fm^3.
     * \param dadT Temperature derivative of the attraction parameter (optional).
     */
    MeanFieldModelVDW(double a, double dadT = 0.) : m_a(a), m_dadT(dadT)
    {
    }
    
    /**
     * \brief Destructor for the MeanFieldModelVDW class.
     */
    virtual ~MeanFieldModelVDW() { }

    /**
     * \brief Calculates the mean field value at a given density.
     * 
     * \param n Particle number density in fm^-3.
     * \return Mean field value in units of GeV/fm^3.
     */
    virtual double v(double n) const { return -m_a * n * n; }

    /**
     * \brief Calculates the density derivatives of the mean field.
     * 
     * \param order Order of the derivative.
     * \param n Particle number density in fm^-3.
     * \return Derivative of the mean field in units of GeV * fm^{3 * (order-1)}.
     */
    virtual double dv(int order, double n) const {
      if (order == 0)
        return v(n);
      if (order == 1)
        return -2. * m_a * n;
      if (order == 2)
        return -2. * m_a;
      return 0.;
    }

    /**
     * \brief Calculates the temperature derivative of the mean field.
     * 
     * \param n Particle number density in fm^-3.
     * \return Temperature derivative of the mean field in units of fm^{-3}.
     */
    virtual double dvdT(double n) const { return -m_dadT * n * n; }
  private:
    double m_a;
    double m_dadT;
  };

  /**
  *   \brief Class implementing auxiliary functions for the Clausius mean-field model.
  *
  *   This class implements the Clausius mean-field model as described in
  *   https://arxiv.org/pdf/1701.06524.pdf, which includes a denominator term
  *   to account for saturation effects at high densities.
  */
  class MeanFieldModelClausius : public MeanFieldModelBase
  {
  public:
    /**
     * \brief Constructor for the MeanFieldModelClausius class.
     * 
     * \param a Attraction parameter in units of GeV*fm^3.
     * \param c Parameter controlling the saturation of the mean field at high densities.
     * \param dadT Temperature derivative of the attraction parameter (optional).
     * \param dcdT Temperature derivative of the c parameter (optional).
     */
    MeanFieldModelClausius(double a, double c, double dadT = 0., double dcdT = 0.) : 
      m_a(a), m_c(c), m_dadT(dadT), m_dcdT(dcdT)
    {
    }
    
    /**
     * \brief Destructor for the MeanFieldModelClausius class.
     */
    virtual ~MeanFieldModelClausius() { }

    /**
     * \brief Calculates the mean field value at a given density.
     * 
     * \param n Particle number density in fm^-3.
     * \return Mean field value in units of GeV/fm^3.
     */
    virtual double v(double n) const { return -m_a * n * n / (1. + m_c * n); }

    /**
     * \brief Calculates the density derivatives of the mean field.
     * 
     * \param order Order of the derivative.
     * \param n Particle number density in fm^-3.
     * \return Derivative of the mean field in units of GeV * fm^{3 * (order-1)}.
     */
    virtual double dv(int order, double n) const {
      if (order == 0)
        return v(n);
      if (order == 1)
        return - m_a * n * (2. + m_c * n) / pow(1. + m_c * n, 2);

      double ret = -2. * m_a / pow(1. + m_c * n, 3);

      if (order == 2)
        return ret;

      for(int i = 3; i <= order; ++i) {
        ret *= -i * m_c;
        ret *= 1. / (1. + m_c * n);
      }

      return ret;
    }

    /**
     * \brief Calculates the temperature derivative of the mean field.
     * 
     * \param n Particle number density in fm^-3.
     * \return Temperature derivative of the mean field in units of fm^{-3}.
     */
    virtual double dvdT(double n) const { 
      return -m_dadT * n * n / (1. + m_c * n) 
        + m_a * n * n / (1. + m_c * n) / (1. + m_c * n) * m_dcdT * n;
    }
  private:
    double m_a, m_c;
    double m_dadT, m_dcdT;
  };

  /**
  *   \brief Class implementing auxiliary functions for the Skyrme mean-field model.
  *
  *   This class implements the Skyrme mean-field model, which is commonly used
  *   in nuclear physics to describe the interaction between nucleons.
  */
  class MeanFieldModelSkyrme : public MeanFieldModelBase
  {
  public:
    /**
     * \brief Constructor for the MeanFieldModelSkyrme class.
     * 
     * \param alpha Parameter for the quadratic term in units of GeV*fm^3.
     * \param beta Parameter for the higher-order term in units of GeV*fm^3.
     * \param gam Exponent for the higher-order term, typically 2.0.
     * \param n0 Saturation density in fm^-3, typically 0.16.
     * \param dalphadT Temperature derivative of the alpha parameter (optional).
     * \param dbetadT Temperature derivative of the beta parameter (optional).
     */
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
    
    /**
     * \brief Destructor for the MeanFieldModelSkyrme class.
     */
    virtual ~MeanFieldModelSkyrme() { }

    /**
     * \brief Calculates the mean field value at a given density.
     * 
     * \param n Particle number density in fm^-3.
     * \return Mean field value in units of GeV/fm^3.
     */
    virtual double v(double n) const { 
      return m_alpha * n * n / m_n0 
        + m_beta * n * pow(n/m_n0, m_gam); 
    }

    /**
     * \brief Calculates the density derivatives of the mean field.
     * 
     * \param order Order of the derivative.
     * \param n Particle number density in fm^-3.
     * \return Derivative of the mean field in units of GeV * fm^{3 * (order-1)}.
     */
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

    /**
     * \brief Calculates the temperature derivative of the mean field.
     * 
     * \param n Particle number density in fm^-3.
     * \return Temperature derivative of the mean field in units of fm^{-3}.
     */
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
  *   \brief Class implementing auxiliary functions for the vector density functional (VDF) model.
  *
  *   This class implements the vector density functional (VDF) model from
  *   https://arxiv.org/pdf/2011.06635, which provides a flexible framework
  *   for describing mean-field interactions.
  */
  class MeanFieldModelVDF : public MeanFieldModelBase
  {
  public:
    /**
     * \brief Constructor for the MeanFieldModelVDF class.
     * 
     * \param N Number of terms in the VDF model.
     * \param Ck Vector of coefficients for each term.
     * \param bk Vector of density exponents for each term.
     * \param n0 Reference density (optional).
     * \param dCkdT Vector of temperature derivatives of the Ck coefficients (optional).
     * \param dbkdT Vector of temperature derivatives of the bk exponents (optional).
     */
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
    
    /**
     * \brief Destructor for the MeanFieldModelVDF class.
     */
    virtual ~MeanFieldModelVDF() { }

    /**
     * \brief Calculates the mean field value at a given density.
     * 
     * \param n Particle number density in fm^-3.
     * \return Mean field value in units of GeV/fm^3.
     */
    virtual double v(double n) const;

    /**
     * \brief Calculates the density derivatives of the mean field.
     * 
     * \param order Order of the derivative.
     * \param n Particle number density in fm^-3.
     * \return Derivative of the mean field in units of GeV * fm^{3 * (order-1)}.
     */
    virtual double dv(int order, double n) const;

    /**
     * \brief Calculates the temperature derivative of the mean field.
     * 
     * \param n Particle number density in fm^-3.
     * \return Temperature derivative of the mean field in units of fm^{-3}.
     */
    virtual double dvdT(double n) const;
  private:
    int m_N;
    std::vector<double> m_Ck, m_bk;
    std::vector<double> m_dCkdT, m_dbkdT;
    double m_n0;
  };

} // namespace thermalfist

#endif
