/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGPhaseShifts/PhaseShiftModel.h"
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <cmath>

namespace thermalfist {

  namespace PhaseShifts {

    namespace {

      // Natural cubic spline of (x,y) data with analytic first derivative.
      // Outside the tabulated range the spline is extrapolated linearly using
      // the end-segment slope, which keeps delta(M) and ddelta/dM finite.
      class CubicSpline {
      public:
        CubicSpline(std::vector<double> x, std::vector<double> y)
          : m_x(std::move(x)), m_y(std::move(y))
        {
          const size_t n = m_x.size();
          if (n < 2)
            throw std::runtime_error("CubicSpline: need at least 2 points");
          // ensure ascending x
          for (size_t i = 1; i < n; ++i)
            if (m_x[i] <= m_x[i - 1])
              throw std::runtime_error("CubicSpline: x values must be strictly increasing");

          m_y2.assign(n, 0.);
          std::vector<double> u(n, 0.);
          // natural boundary conditions: y2[0] = y2[n-1] = 0
          for (size_t i = 1; i + 1 < n; ++i) {
            double sig = (m_x[i] - m_x[i - 1]) / (m_x[i + 1] - m_x[i - 1]);
            double p = sig * m_y2[i - 1] + 2.;
            m_y2[i] = (sig - 1.) / p;
            double d = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i])
                     - (m_y[i] - m_y[i - 1]) / (m_x[i] - m_x[i - 1]);
            u[i] = (6. * d / (m_x[i + 1] - m_x[i - 1]) - sig * u[i - 1]) / p;
          }
          for (size_t k = n - 1; k-- > 0; )
            m_y2[k] = m_y2[k] * m_y2[k + 1] + u[k];
        }

        double f(double x) const {
          const size_t n = m_x.size();
          if (x <= m_x[0])
            return m_y[0] + slope(0) * (x - m_x[0]);            // linear extrapolation
          if (x >= m_x[n - 1])
            return m_y[n - 1] + slope(n - 2) * (x - m_x[n - 1]);
          size_t lo = interval(x);
          double h = m_x[lo + 1] - m_x[lo];
          double a = (m_x[lo + 1] - x) / h, b = (x - m_x[lo]) / h;
          return a * m_y[lo] + b * m_y[lo + 1]
               + ((a * a * a - a) * m_y2[lo] + (b * b * b - b) * m_y2[lo + 1]) * (h * h) / 6.;
        }

        double df(double x) const {
          const size_t n = m_x.size();
          if (x <= m_x[0])    return slope(0);
          if (x >= m_x[n - 1]) return slope(n - 2);
          size_t lo = interval(x);
          double h = m_x[lo + 1] - m_x[lo];
          double a = (m_x[lo + 1] - x) / h, b = (x - m_x[lo]) / h;
          return (m_y[lo + 1] - m_y[lo]) / h
               - (3. * a * a - 1.) / 6. * h * m_y2[lo]
               + (3. * b * b - 1.) / 6. * h * m_y2[lo + 1];
        }

      private:
        // slope of the cubic at the left node of segment i (= secant for endpoints
        // since y2 vanishes there under natural BC, giving clean linear extrapolation)
        double slope(size_t i) const {
          double h = m_x[i + 1] - m_x[i];
          return (m_y[i + 1] - m_y[i]) / h
               - h / 6. * (2. * m_y2[i] + m_y2[i + 1]);
        }
        size_t interval(double x) const {
          // largest lo with m_x[lo] <= x, x strictly inside the range here
          size_t lo = std::upper_bound(m_x.begin(), m_x.end(), x) - m_x.begin();
          return lo - 1;
        }
        std::vector<double> m_x, m_y, m_y2;
      };

      // Read a two-column "M delta" text file (skips blank and '#' lines).
      void readTable(const std::string& path, std::vector<double>& M, std::vector<double>& d) {
        std::ifstream fin(path.c_str());
        if (!fin.is_open())
          throw std::runtime_error("TabulatedModel: cannot open file " + path);
        std::string line;
        while (std::getline(fin, line)) {
          size_t p = line.find_first_not_of(" \t\r\n");
          if (p == std::string::npos || line[p] == '#') continue;
          std::istringstream iss(line);
          double mm, dd;
          if (iss >> mm >> dd) { M.push_back(mm); d.push_back(dd); }
        }
        if (M.size() < 2)
          throw std::runtime_error("TabulatedModel: file " + path + " has fewer than 2 points");
      }

    } // anonymous namespace

    PhaseShiftPartialWave AnalyticWave(const std::string& channel, const std::string& param) {
      // Per-wave registry: the model NAME is wave-specific, so (channel, param)
      // maps to exactly one wave. The wave's 2J+1 comes from the catalog wave set.
      std::vector<PhaseShiftPartialWave> waves;
      int twoJplus1 = -1;
      if (channel == "pipi_I2") {
        waves = PiPi_I2_Waves();
        if      (param == "GarciaMartin2011_S") twoJplus1 = 1;
        else if (param == "GarciaMartin2011_D") twoJplus1 = 5;
      }
      else if (channel == "pipi_I0") {
        waves = PiPi_I0_Waves();
        if (param == "GarciaMartin2011_S") twoJplus1 = 1;
      }
      else if (channel == "pipi_I0_f0980") {
        waves = PiPi_I0_f0980_Waves();
        if (param == "GarciaMartin2011_S") twoJplus1 = 1;
      }
      else if (channel == "piK_I32") {
        waves = PiK_I32_Waves();
        if (param == "PelaezRodas2016_S") twoJplus1 = 1;
      }
      else if (channel == "piK_I12") {
        waves = PiK_I12_Waves();
        if (param == "PelaezRodas2016_S") twoJplus1 = 1;
      }
      if (twoJplus1 < 0)
        throw std::invalid_argument("AnalyticWave: unknown analytic model '"
                                    + channel + ":" + param + "'");
      for (size_t i = 0; i < waves.size(); ++i)
        if (waves[i].twoJplus1 == twoJplus1) return waves[i];
      throw std::invalid_argument("AnalyticWave: model '" + channel + ":" + param
                                  + "' wave (2J+1=" + std::to_string(twoJplus1) + ") unavailable");
    }

    PhaseShiftPartialWave TabulatedWave(int twoJplus1, const std::string& file) {
      std::vector<double> M, d;
      readTable(file, M, d);
      std::shared_ptr<CubicSpline> sp = std::make_shared<CubicSpline>(M, d);
      return PhaseShiftPartialWave(
        twoJplus1,
        [sp](double m) { return sp->f(m); },     // delta(M)
        [sp](double m) { return sp->df(m); });   // analytic ddelta/dM
    }

  } // namespace PhaseShifts

} // namespace thermalfist
