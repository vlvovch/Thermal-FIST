/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

/**
 * \file NumericalIntegration.h
 *
 * Gauss-Legendre and Gauss-Laguerre quadratures used in numerical integrations.
 *
 * Nodes and weights are computed on demand (and cached) by the in-house
 * Golub-Welsch solver in NumericalIntegration/GaussianQuadrature.h, rather than
 * stored as hard-coded tables. The interface below therefore supports an
 * arbitrary number of nodes: the historical fixed-order helpers
 * (GetCoefsIntegrateLegendre5/10/32/40, GetCoefsIntegrateLaguerre32) are thin
 * wrappers over the parameterized GetCoefsIntegrateLegendre(n, ...) /
 * GetCoefsIntegrateLaguerre(n, ...).
 */

#include <vector>
#include <stdexcept>

#include "HRGBase/NumericalIntegration/OrthogonalPolynomials.h"
#include "HRGBase/NumericalIntegration/GaussianQuadrature.h"

namespace thermalfist {

  /// \brief Contains various Gauss-Legendre and Gauss-Laguerre quadratures
  ///        used in numerical integrations.
  namespace NumericalIntegration {

    // NOTE (API change): the hard-coded coefficient tables that used to live in
    // this header --
    //     coefficients_xleg5/10/32/40, coefficients_wleg5/10/32/40,
    //     coefficients_xleg32_zeroone, coefficients_wleg32_zeroone,
    //     coefficients_xlag32, coefficients_wlag32
    // -- have been removed. Nodes and weights are now computed and cached on
    // demand. Code that referenced those arrays should switch to the cached
    // getters below (e.g. coefficients_xleg32  -> GetGaussLegendreNodes(32),
    // coefficients_wlag32 -> GetGaussLaguerreWeights(32), the _zeroone arrays ->
    // GetGaussLegendre{Nodes,Weights}ZeroOne(32)) or to GetCoefsIntegrate*.

    /** \name Cached base quadrature rules on the reference intervals.
     *  Each rule is computed once per node count and cached (thread-safe).
     *  The returned references are stable for the lifetime of the program.
     *  @{ */

    /// Gauss-Legendre nodes on [-1, 1] (sorted ascending).
    const std::vector<double>& GetGaussLegendreNodes(int n);
    /// Gauss-Legendre weights on [-1, 1].
    const std::vector<double>& GetGaussLegendreWeights(int n);

    /// Gauss-Legendre nodes mapped to [0, 1].
    const std::vector<double>& GetGaussLegendreNodesZeroOne(int n);
    /// Gauss-Legendre weights mapped to [0, 1].
    const std::vector<double>& GetGaussLegendreWeightsZeroOne(int n);

    /// Gauss-Laguerre nodes on [0, inf) (weight e^{-x}).
    const std::vector<double>& GetGaussLaguerreNodes(int n);
    /// Gauss-Laguerre *modified* weights w_k * e^{x_k}, i.e. such that
    /// \f$ \int_0^\infty f(x)\,dx \approx \sum_k w_k f(x_k) \f$.
    const std::vector<double>& GetGaussLaguerreWeights(int n);

    /** @} */

    /**
     * Populates nodes and weights for integrating f(x) on [a, b] using an
     * n-point Gauss-Legendre quadrature.
     */
    void GetCoefsIntegrateLegendre(int n, double a, double b, std::vector<double>* x, std::vector<double>* w);

    /**
     * Populates nodes and weights for integrating f(x) on [0, inf) using an
     * n-point Gauss-Laguerre quadrature. The weights are the modified weights
     * (see GetGaussLaguerreWeights), i.e. the e^{-x} factor is folded in.
     */
    void GetCoefsIntegrateLaguerre(int n, std::vector<double>* x, std::vector<double>* w);

    /**
     * Integrates a function f(x,y)
     * in range 0 < x < \infty, ay <= y <= by
     * using the combined 32-point
     * Gauss-Laguerre and the 32-point Gauss-Legendre quadrature.
     *
     * \param func A pointer to a function to be integrated.
     * \param ay Left limit of integration for variable y.
     * \param by Right limit of integration for variable y.
     * \return Result of the integration.
     */
    double Integrate2DLaguerre32Legendre32(double(*func)(double, double), double ay, double by);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x,y)
     * in the range 0 < x < \infty, ay <= y <= by
     * using the combined 32-point
     * Gauss-Laguerre and the 32-point Gauss-Legendre quadrature.
     *
     * \param [in] ay Left limit of integration for variable y.
     * \param [in] by Right limit of integration for variable y.
     * \param [out] xlag Gauss-Laguerre nodes for the variable x.
     * \param [out] wlag Gauss-Laguerre weights for the variable x.
     * \param [out] xleg Gauss-Legendre nodes for the variable y.
     * \param [out] wleg Gauss-Legendre weights for the variable y.
     */
    void GetCoefs2DLaguerre32Legendre32(double ay, double by, std::vector<double> *xlag, std::vector<double> *wlag, std::vector<double> *xleg, std::vector<double> *wleg);

    /**
     * Populates the nodes and weights for integrating
     * a function f(x,y)
     * in the range ax < x < bx, ay <= y <= by
     * using two 32-point Gauss-Legendre quadratures.
     *
     * \param [in] ax Left limit of integration for variable x.
     * \param [in] bx Right limit of integration for variable x.
     * \param [in] ay Left limit of integration for variable y.
     * \param [in] by Right limit of integration for variable y.
     * \param [out] xlag Gauss-Legendre nodes for the variable x.
     * \param [out] wlag Gauss-Legendre weights for the variable x.
     * \param [out] xleg Gauss-Legendre nodes for the variable y.
     * \param [out] wleg Gauss-Legendre weights for the variable y.
     */
    void GetCoefs2DLegendre32Legendre32(double ax, double bx, double ay, double by, std::vector<double> *xleg1, std::vector<double> *wleg1, std::vector<double> *xleg2, std::vector<double> *wleg2);

    /**
     * Populates the nodes and weights for integrating f(x) on [a, b] using the
     * 32-point Gauss-Legendre quadrature.
     */
    void GetCoefsIntegrateLegendre32(double a, double b, std::vector<double> *x, std::vector<double> *w);

    /**
     * Populates the nodes and weights for integrating f(x) on [a, b] using the
     * 10-point Gauss-Legendre quadrature.
     */
    void GetCoefsIntegrateLegendre10(double a, double b, std::vector<double> *x, std::vector<double> *w);

    /**
     * Populates the nodes and weights for integrating f(x) on [a, b] using the
     * 5-point Gauss-Legendre quadrature.
     */
    void GetCoefsIntegrateLegendre5(double a, double b, std::vector<double> *x, std::vector<double> *w);

    /**
     * Populates the nodes and weights for integrating f(x) on [a, b] using the
     * 40-point Gauss-Legendre quadrature.
     */
    void GetCoefsIntegrateLegendre40(double a, double b, std::vector<double> *x, std::vector<double> *w);

    /**
     * Populates the nodes and weights for integrating f(x) on [0, inf) using the
     * 32-point Gauss-Laguerre quadrature.
     */
    void GetCoefsIntegrateLaguerre32(std::vector<double> *x, std::vector<double> *w);


    /**
     * Template version.
     * Populates the nodes and weights for integrating f(x) on [a, b] using an
     * n-point Gauss-Legendre quadrature, for an arbitrary scalar type T
     * (e.g. an automatic-differentiation type).
     */
    template <typename T>
    void GetCoefsIntegrateLegendreTemplate(int n, T a, T b, std::vector<T> *xp, std::vector<T> *wp) {
      if (n <= 0)
        throw std::invalid_argument("GetCoefsIntegrateLegendreTemplate: number of nodes n must be positive");
      std::vector<T> &x = *xp;
      std::vector<T> &w = *wp;
      x.resize(n);
      w.resize(n);

      const std::vector<double>& xlego = GetGaussLegendreNodes(n);
      const std::vector<double>& wlego = GetGaussLegendreWeights(n);

      for (int i = 0; i < n; i++) {
        w[i] = (b - a) / 2. * wlego[i];
        x[i] = (b - a) / 2. * xlego[i] + (b + a) / 2.;
      }
    }

    /**
     * Template version.
     * Populates the nodes and weights for integrating f(x) on [a, b] using the
     * 32-point Gauss-Legendre quadrature.
     */
    template <typename T>
    void GetCoefsIntegrateLegendre32Template(T a, T b, std::vector<T> *xp, std::vector<T> *wp) {
      GetCoefsIntegrateLegendreTemplate<T>(32, a, b, xp, wp);
    }
  }

} // namespace thermalfist

#endif
