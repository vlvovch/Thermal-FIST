/*
/* * GNU General Public License (GPLv3 or later)
 */
#ifndef SPLINEFUNCTION_H
#define SPLINEFUNCTION_H

#include <vector>
#include <algorithm>
#include <cmath>

namespace thermalfist {

  /// Class implementing a simple linear spline.
  class SplineFunction
  {
  public:
    SplineFunction() : m_vals(0) {
      m_vals.resize(0);
    }

    SplineFunction(std::vector<double> x, std::vector<double> y) : m_vals(0) {
      for (unsigned int i = 0; i < x.size(); ++i)
      {
        m_vals.push_back(std::make_pair(x[i], y[i]));
      }
      sort(m_vals.begin(), m_vals.end());
    }

    /// Adds a new pair of x,y values.
    void add_val(double x, double val)
    {
      m_vals.push_back(std::make_pair(x, val));
      sort(m_vals.begin(), m_vals.end());
    }

    /// Evaluates interpolated function at x = arg.
    double f(double arg) const
    {
      unsigned int ind = 0;
      std::pair<double, double> op = std::make_pair(arg, 0.);
      std::vector< std::pair<double, double> >::const_iterator it = std::lower_bound(m_vals.begin(), m_vals.end(), op);
      ind = distance(m_vals.begin(), it);

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
    double fsquare(double arg) {
      double ret = f(arg);
      return ret * ret;
    }

    /// Clear all data and refill with zero function.
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
    void fill(std::vector<double> x, std::vector<double> y) {
      m_vals.resize(0);
      for (unsigned int i = 0; i < x.size(); ++i)
      {
        m_vals.push_back(std::make_pair(x[i], y[i]));
      }
      sort(m_vals.begin(), m_vals.end());
    }

    /// Models constnat f(x) == val function.
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
