/*
 * * GNU General Public License (GPLv3 or later)
 */
#ifndef SPLINEFUNCTION_H
#define SPLINEFUNCTION_H

#include <vector>
#include <algorithm>
#include <cmath>

namespace thermalfist {

  /// Class implementing a simple linear spline.
  /**
   * This class implements a simple linear spline function f(x) of one argument.
   * It is used to model a function f(x) which is defined by a set of (x,y) pairs.
   * The function is evaluated using linear interpolation between the provided
   * (x,y) pairs.
   */
  class SplineFunction
  {
  public:
    /**
     * Default constructor. Empty function.
     */
    SplineFunction() : m_vals(0) {
    }

    /**
     * Constructor which sets the data from the provided vectors.
     * \param x A vector of x values.
     * \param y A vector of y values.
     * Data must be sorted in non-descending order of x values.
     */
    SplineFunction(std::vector<double> x, std::vector<double> y) : m_vals(0) {
      for (unsigned int i = 0; i < x.size(); ++i)
      {
        m_vals.push_back(std::make_pair(x[i], y[i]));
      }
      sort(m_vals.begin(), m_vals.end());
    }

    /// Adds a new pair of x,y values.
    /**
     * \param x A new x value.
     * \param val A new y value.
     * Adds a new (x,y) pair to the function. The function is sorted after
     * adding the new value.
     */
    void add_val(double x, double val)
    {
      m_vals.push_back(std::make_pair(x, val));
      sort(m_vals.begin(), m_vals.end());
    }

    /// Evaluates interpolated function at x = arg.
    /**
     * \param arg The x value at which to evaluate the function.
     * \return The value of the function at x = arg.
     * Evaluates the function at x = arg using linear interpolation between
     * the provided (x,y) pairs.
     */
    double f(double arg) const
    {
      std::pair<double, double> op = std::make_pair(arg, 0.);
      std::vector< std::pair<double, double> >::const_iterator it = std::lower_bound(m_vals.begin(), m_vals.end(), op);
      unsigned int ind = distance(m_vals.begin(), it);

      if (ind == 0) return m_vals[0].second +
        (arg - m_vals[0].first) *
        (m_vals[1].second - m_vals[0].second) / (m_vals[1].first - m_vals[0].first);

      if (ind == m_vals.size()) return m_vals[ind - 2].second +
        (arg - m_vals[ind - 2].first) *
        (m_vals[ind - 1].second - m_vals[ind - 2].second) / (m_vals[ind - 1].first - m_vals[ind - 2].first);

      return m_vals[ind - 1].second +
        (arg - m_vals[ind - 1].first) *
        (m_vals[ind].second - m_vals[ind - 1].second) / (m_vals[ind].first - m_vals[ind - 1].first);
    }

    /// Evaluates slope (derivative) at x = arg.
    /**
     * \param arg The x value at which to evaluate the derivative.
     * \return The value of the derivative at x = arg.
     * Evaluates the derivative of the function at x = arg using linear
     * interpolation between the provided (x,y) pairs.
     */
    double df(double arg) const {
      unsigned int ind = 0;
      std::pair<double, double> op = std::make_pair(arg, 0.);
      std::vector< std::pair<double, double> >::const_iterator it = std::lower_bound(m_vals.begin(), m_vals.end(), op);
      ind = std::distance(m_vals.begin(), it);
      if (ind == 0)
        return (m_vals[1].second - m_vals[0].second) / (m_vals[1].first - m_vals[0].first);
      if (ind == m_vals.size())
        return (m_vals[ind - 1].second - m_vals[ind - 2].second) / (m_vals[ind - 1].first - m_vals[ind - 2].first);
      return (m_vals[ind].second - m_vals[ind - 1].second) / (m_vals[ind].first - m_vals[ind - 1].first);
    }

    /// Evaluates f(arg)^2.
    /**
     * \param arg The x value at which to evaluate the function.
     * \return The value of the function at x = arg squared.
     * Evaluates the function at x = arg using linear interpolation between
     * the provided (x,y) pairs and returns the square of the result.
     */
    double fsquare(double arg) {
      double ret = f(arg);
      return ret * ret;
    }

    /**
     * \return The value of the function at x = arg.
     * Clears all data and fills the function with a constant function
     * f(x) = 0.
     */
    void clear() {
      m_vals.resize(2);
      m_vals[0].first = 0.;
      m_vals[0].second = 0.;
      m_vals[1].first = 1.;
      m_vals[1].second = 0.;
    }

    /// Just clear all data.
    void clearall() {
      m_vals.resize(0);
    }

    /// Fill (x,y) pairs from provided vectors.
    /**
     * \param x A vector of x values.
     * \param y A vector of y values.
     * Data must be sorted in non-descending order of x values.
     */
    void fill(std::vector<double> x, std::vector<double> y) {
      m_vals.resize(0);
      for (unsigned int i = 0; i < x.size(); ++i)
      {
        m_vals.push_back(std::make_pair(x[i], y[i]));
      }
      sort(m_vals.begin(), m_vals.end());
    }

    /// Models constant f(x) == val function.
    /**
     * \param val The value of the constant function.
     * Clears all data and fills the function with a constant function
     * f(x) = val.
     */
    void setConstant(double val) {
      m_vals.resize(0);
      m_vals.push_back(std::make_pair(0., val));
      m_vals.push_back(std::make_pair(1., val));
    }

    // TODO: Read (x,y) pairs from file.
    //void loadFromFile(const char *file);

  private:
    std::vector< std::pair<double, double> > m_vals;    /**< A set of x and y values of a discrete function */
  };

} // namespace thermalfist

#endif // SPLINEFUNCTION_H
