/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/NumericalIntegration.h"

#include <map>
#include <mutex>
#include <cmath>
#include <utility>
#include <stdexcept>

namespace thermalfist {

  namespace NumericalIntegration {

    namespace {
      // Shared, mutex-guarded caches of base quadrature rules keyed by node count.
      // One lock per rule-setup call (negligible next to the n function evaluations
      // that follow); std::map nodes are stable, so returned references stay valid.
      std::mutex g_cacheMutex;

      typedef std::pair<std::vector<double>, std::vector<double>> Rule;

      const Rule& CachedLegendre(int n) {
        static std::map<int, Rule> cache;
        std::lock_guard<std::mutex> lock(g_cacheMutex);
        auto it = cache.find(n);
        if (it == cache.end()) {
          LegendrePolynomial poly;
          it = cache.emplace(n, GaussianQuadrature(n, &poly)).first;
        }
        return it->second;
      }

      const Rule& CachedLegendreZeroOne(int n) {
        static std::map<int, Rule> cache;
        std::lock_guard<std::mutex> lock(g_cacheMutex);
        auto it = cache.find(n);
        if (it == cache.end()) {
          // Compute the [-1,1] base rule directly here (not via CachedLegendre)
          // to avoid re-locking the non-recursive mutex.
          LegendrePolynomial poly;
          Rule base = GaussianQuadrature(n, &poly);
          Rule r;
          r.first.resize(n);
          r.second.resize(n);
          for (int k = 0; k < n; ++k) {
            r.first[k]  = 0.5 * (base.first[k] + 1.0);
            r.second[k] = 0.5 * base.second[k];
          }
          it = cache.emplace(n, std::move(r)).first;
        }
        return it->second;
      }

      const Rule& CachedLaguerre(int n) {
        static std::map<int, Rule> cache;
        std::lock_guard<std::mutex> lock(g_cacheMutex);
        auto it = cache.find(n);
        if (it == cache.end()) {
          LaguerrePolynomial poly;
          Rule r = GaussianQuadrature(n, &poly);   // standard weights (sum = 1)
          for (int k = 0; k < n; ++k)
            r.second[k] *= std::exp(r.first[k]);    // modified weights w_k * e^{x_k}
          it = cache.emplace(n, std::move(r)).first;
        }
        return it->second;
      }
    } // anonymous namespace

    // ---- Cached base rules -------------------------------------------------

    const std::vector<double>& GetGaussLegendreNodes(int n)          { return CachedLegendre(n).first; }
    const std::vector<double>& GetGaussLegendreWeights(int n)        { return CachedLegendre(n).second; }
    const std::vector<double>& GetGaussLegendreNodesZeroOne(int n)   { return CachedLegendreZeroOne(n).first; }
    const std::vector<double>& GetGaussLegendreWeightsZeroOne(int n) { return CachedLegendreZeroOne(n).second; }
    const std::vector<double>& GetGaussLaguerreNodes(int n)          { return CachedLaguerre(n).first; }
    const std::vector<double>& GetGaussLaguerreWeights(int n)        { return CachedLaguerre(n).second; }

    // ---- Parameterized 1D coefficient helpers ------------------------------

    void GetCoefsIntegrateLegendre(int n, double a, double b, std::vector<double> *xp, std::vector<double> *wp) {
      if (n <= 0)
        throw std::invalid_argument("GetCoefsIntegrateLegendre: number of nodes n must be positive");
      std::vector<double> &x = *xp;
      std::vector<double> &w = *wp;
      x.resize(n);
      w.resize(n);

      const std::vector<double>& xlego = GetGaussLegendreNodes(n);
      const std::vector<double>& wlego = GetGaussLegendreWeights(n);

      for (int i = 0; i < n; i++) {
        w[i] = (b - a) / 2. * wlego[i];
        x[i] = (b - a) / 2. * xlego[i] + (b + a) / 2.;
      }
    }

    void GetCoefsIntegrateLaguerre(int n, std::vector<double> *xp, std::vector<double> *wp) {
      if (n <= 0)
        throw std::invalid_argument("GetCoefsIntegrateLaguerre: number of nodes n must be positive");
      std::vector<double> &x = *xp;
      std::vector<double> &w = *wp;
      x.resize(n);
      w.resize(n);

      const std::vector<double>& xlago = GetGaussLaguerreNodes(n);
      const std::vector<double>& wlago = GetGaussLaguerreWeights(n);

      for (int i = 0; i < n; i++) {
        x[i] = xlago[i];
        w[i] = wlago[i];   // modified weights (e^{-x} folded in)
      }
    }

    // ---- 2D helpers --------------------------------------------------------

    double Integrate2DLaguerre32Legendre32(double(*func)(double, double), double ay, double by)
    {
      const std::vector<double>& xleg = GetGaussLegendreNodes(32);
      const std::vector<double>& wleg = GetGaussLegendreWeights(32);
      const std::vector<double>& xlag = GetGaussLaguerreNodes(32);
      const std::vector<double>& wlag = GetGaussLaguerreWeights(32);

      double sum = 0.;
      for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) {
          double x = (by - ay) / 2.*xleg[j] + (by + ay) / 2.;
          double w = (by - ay) / 2.*wleg[j];
          sum += wlag[i] * w * func(xlag[i], x);
        }
      }
      return sum;
    }

    void GetCoefs2DLaguerre32Legendre32(double ay, double by,
      std::vector<double> *xlagp, std::vector<double> *wlagp,
      std::vector<double> *xlegp, std::vector<double> *wlegp) {
      GetCoefsIntegrateLegendre(32, ay, by, xlegp, wlegp);
      GetCoefsIntegrateLaguerre(32, xlagp, wlagp);
    }

    void GetCoefs2DLegendre32Legendre32(double ay, double by, double a2y, double b2y,
      std::vector<double> *xlegp1, std::vector<double> *wlegp1,
      std::vector<double> *xlegp2, std::vector<double> *wlegp2) {
      GetCoefsIntegrateLegendre(32, ay, by, xlegp1, wlegp1);
      GetCoefsIntegrateLegendre(32, a2y, b2y, xlegp2, wlegp2);
    }

    // ---- Fixed-order wrappers (preserved API) ------------------------------

    void GetCoefsIntegrateLegendre32(double a, double b, std::vector<double> *xp, std::vector<double> *wp) {
      GetCoefsIntegrateLegendre(32, a, b, xp, wp);
    }

    void GetCoefsIntegrateLegendre10(double a, double b, std::vector<double> *xp, std::vector<double> *wp) {
      GetCoefsIntegrateLegendre(10, a, b, xp, wp);
    }

    void GetCoefsIntegrateLegendre5(double a, double b, std::vector<double> *xp, std::vector<double> *wp) {
      GetCoefsIntegrateLegendre(5, a, b, xp, wp);
    }

    void GetCoefsIntegrateLegendre40(double a, double b, std::vector<double> *xp, std::vector<double> *wp) {
      GetCoefsIntegrateLegendre(40, a, b, xp, wp);
    }

    void GetCoefsIntegrateLaguerre32(std::vector<double> *xp, std::vector<double> *wp) {
      GetCoefsIntegrateLaguerre(32, xp, wp);
    }

  } // namespace NumericalIntegration

} // namespace thermalfist
