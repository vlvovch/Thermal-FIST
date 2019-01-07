/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef BILINEARSPLINEFUNCTION_H
#define BILINEARSPLINEFUNCTION_H

#include "HRGBase/SplineFunction.h"

namespace thermalfist {

  /// A class implementing a bilinear spline.
  /**
  * Implementation of bilinear spline function f(x,y) of two arguments.
  * Uses a 1D spline function tp model f(y;x) at each discrete y
  * Requires that (x,y) values form a grid.
  */
  class BilinearSplineFunction
  {
  public:
    /**
    * Default constructror. Empty function.
    */
    BilinearSplineFunction(void) : m_xs(), m_xspls() {
      m_xs.resize(0);
      m_xspls.resize(0);
    }

    /**
    * Constructor which sets the data from the provided vectors.
    * \param x A vector of x values.
    * \param y A vector of y values.
    * \param vals A vector of f(x,y) values.
    * Data must be sorted in non-descending order of x values, and then in
    * non-descending order of y values for equal x values.
    * Provided (x,y) set must form a grid.
    */
    BilinearSplineFunction(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &vals)
      : m_xs(), m_xspls() {
      setData(x, y, vals);
    }

    /**
    * Method which sets the data from the provided vectors.
    * \param x A vector of x values.
    * \param y A vector of y values.
    * \param vals A vector of f(x,y) values.
    * Data must be sorted in non-descending order of x values, and then in
    * non-descending order of y values for equal x values.
    * Provided (x,y) set must form a grid.
    */
    void setData(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &vals) {
      if (x.size() > 0) {
        m_xs.resize(0);
        m_xspls.resize(0);
        double cx = -1e50;
        for (unsigned int i = 0; i < x.size(); ++i) {
          if (fabs(x[i] - cx) > 1e-6) {
            m_xspls.push_back(SplineFunction());
            m_xs.push_back(x[i]);
            m_xspls[m_xspls.size() - 1].add_val(y[i], vals[i]);
            cx = x[i];
          }
          else {
            m_xspls[m_xspls.size() - 1].add_val(y[i], vals[i]);
          }
        }
      }
    }

    /// Evaluates interpolated f(x,y)
    double Eval(double x, double y) const {
      if (m_xs.size() < 2) return -1.;
      unsigned int indx = 0;
      std::vector< double >::const_iterator it = std::lower_bound(m_xs.begin(), m_xs.end(), x);
      indx = std::distance(m_xs.begin(), it);
      int ind1 = 0, ind2 = 0;
      if (indx == 0) {
        ind1 = 0;
        ind2 = 1;
      }
      else if (indx == m_xs.size()) {
        ind1 = indx - 2;
        ind2 = indx - 1;
      }
      else {
        ind1 = indx - 1;
        ind2 = indx;
      }
      double f1v = m_xspls[ind1].f(y);
      double f2v = m_xspls[ind2].f(y);
      return f1v + (x - m_xs[ind1]) * (f2v - f1v) / (m_xs[ind2] - m_xs[ind1]);
    }

    /// Destructor.
    ~BilinearSplineFunction(void) { }

  private:
    std::vector<double> m_xs;
    std::vector<SplineFunction> m_xspls;
  };

} // namespace thermalfist

#endif
