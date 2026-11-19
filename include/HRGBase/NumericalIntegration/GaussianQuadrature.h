/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 *
 * Vendored from GaussQuad (https://github.com/vlvovch/gaussquad-cpp), v1.1.0,
 * commit a1dc387. Local changes: namespace gaussquad -> thermalfist::NumericalIntegration,
 * include paths and include guards, and InitialGuessLaguerre/InitialGuessHermite
 * use root-by-root Numerical Recipes guesses (the upstream asymptotic seeds
 * converge to duplicate roots). Keep in sync with upstream.
 */
/**
 * @file GaussianQuadrature.h
 * @brief Header file containing functions to compute nodes and weights for Gaussian quadrature.
 *
 * This file implements various methods for computing Gaussian quadrature rules, which are
 * numerical integration schemes that integrate polynomials exactly up to a certain degree.
 * The quadrature rules are based on orthogonal polynomials and their zeros.
 *
 * The orthogonal polynomials are defined by a recurrence relation:
 * p_{n+1}(x) = (k_{n+1} x - a_n)p_n(x) - b_n p_{n-1}(x)
 *
 * Two main algorithms are implemented:
 * 1. Newton-Raphson method with initial guesses
 * 2. Golub-Welsch algorithm via an in-house tridiagonal QL eigensolver
 *
 * Supported polynomial families:
 * - Legendre polynomials (interval [-1,1], weight w(x)=1)
 * - Laguerre polynomials (interval [0,∞), weight w(x)=exp(-x))
 * - Generalized Laguerre polynomials (interval [0,∞), weight w(x)=x^α*exp(-x))
 * - Hermite polynomials (interval (-∞,∞), weight w(x)=exp(-x²))
 * - Jacobi polynomials (interval [-1,1], weight w(x)=(1-x)^α*(1+x)^β)
 *
 * This is a header-only, dependency-free library: it requires only the C++
 * standard library (no Eigen, no third-party packages). Just #include it.
 */

#ifndef NUMERICALINTEGRATION_GAUSSIANQUADRATURE_H
#define NUMERICALINTEGRATION_GAUSSIANQUADRATURE_H
#include <vector>
#include <utility>
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <string>
#include "HRGBase/NumericalIntegration/OrthogonalPolynomials.h"

namespace thermalfist {
namespace NumericalIntegration {

  /**
   * @brief Compute nodes and weights for Gaussian quadrature using Newton-Raphson method.
   *
   * This function refines the initial guess for the nodes using the Newton-Raphson method.
   * The nodes are the zeros of the orthogonal polynomial of degree n, and the weights are
   * computed using the derivative of the polynomial at the nodes.
   *
   * @param n Number of quadrature nodes (polynomial degree)
   * @param x0 Initial guess for the nodes (must have size n)
   * @param poly Pointer to the orthogonal polynomial object
   * @param TOL Tolerance for the Newton-Raphson convergence (default: 1e-15)
   * @param MAX_ITER Maximum number of iterations for Newton-Raphson (default: 100)
   *
   * @return Pair of vectors: first contains nodes, second contains weights
   * @throws std::invalid_argument if x0.size() != n
   * @throws std::runtime_error if Newton-Raphson fails to converge for any node
   */
  std::pair<std::vector<double>,std::vector<double>> GaussianQuadraturesNewton(int n,
                                                                               const std::vector<double> &x0,
                                                                               const OrthogonalPolynomial* poly,
                                                                               const double TOL = 1e-15,
                                                                               const int MAX_ITER = 100);

  /**
   * @brief Generate initial guess for Legendre polynomial nodes.
   *
   * Uses the asymptotic formula for Legendre polynomial zeros as initial estimates.
   * The approximation is: x_k ≈ cos(π(4k+3)/(4n+2)) for k = 0, 1, ..., n-1
   *
   * @param n Number of nodes
   * @return Vector of initial node estimates
   */
  std::vector<double> InitialGuessLegendre(int n);

  /**
   * @brief Generate initial guess for Laguerre polynomial nodes.
   *
   * Uses asymptotic approximations for Laguerre polynomial zeros.
   * The approximation involves trigonometric functions and corrections.
   *
   * @param n Number of nodes
   * @return Vector of initial node estimates
   */
  std::vector<double> InitialGuessLaguerre(int n);

  /**
   * @brief Generate initial guess for Hermite polynomial nodes.
   *
   * Uses asymptotic approximations for Hermite polynomial zeros.
   * The approximation is based on the eigenvalues of the Hermite operator.
   *
   * @param n Number of nodes
   * @return Vector of initial node estimates
   */
  std::vector<double> InitialGuessHermite(int n);

  /**
   * @brief Compute nodes using the Golub-Welsch algorithm.
   *
   * Convenience wrapper that returns just the nodes from GolubWelschNodesWeights
   * (the eigenvalues of the symmetric tridiagonal Jacobi matrix).
   *
   * @param n Number of nodes
   * @param poly Pointer to the orthogonal polynomial object
   * @return Vector of quadrature nodes (sorted)
   */
  std::vector<double> NodesGolubWelsch(int n, const OrthogonalPolynomial* poly);

  /**
   * @brief Compute nodes and weights using the full Golub-Welsch algorithm.
   *
   * Builds the diagonal/subdiagonal bands of the symmetric tridiagonal Jacobi
   * matrix from the recurrence coefficients and solves it with an in-house
   * implicit-shift QL iteration that tracks only the first row of the
   * eigenvector matrix (all the Golub-Welsch weight formula needs). The
   * eigenvalues are the quadrature nodes; the weights are w_k = mu0 * v_k[0]^2,
   * where mu0 = Norm(0) is the integral of the weight function.
   *
   * Working directly on the bands makes this O(n^2) time and O(n) memory with no
   * external dependency. Unlike the Newton-Raphson weight formula it never
   * evaluates the orthogonal polynomial, so it does not overflow for high-order
   * Laguerre/Hermite quadrature, and the sum rule holds to machine precision by
   * the orthonormality of the eigenvectors.
   *
   * If poly->HasEvenWeight() is true (Legendre, Hermite, Jacobi with
   * alpha == beta), the rule is made exactly symmetric about the origin:
   * x[k] == -x[n-1-k] and w[k] == w[n-1-k], and the middle node is exactly 0
   * for odd n.
   *
   * @param n Number of nodes
   * @param poly Pointer to the orthogonal polynomial object
   * @return Pair of vectors: first contains nodes (sorted), second the weights
   */
  std::pair<std::vector<double>,std::vector<double>> GolubWelschNodesWeights(int n, const OrthogonalPolynomial* poly);

  /**
   * @brief Compute Gauss-Legendre quadrature nodes and weights.
   *
   * Convenience function specifically for Legendre polynomials.
   * Uses Newton-Raphson method with asymptotic initial guesses.
   *
   * @param n Number of quadrature nodes
   * @return Pair of vectors: first contains nodes, second contains weights
   */
  std::pair<std::vector<double>,std::vector<double>> GaussianQuadraturesLegendre(int n);

  /**
   * @brief Compute Gaussian quadrature nodes and weights for any orthogonal polynomial.
   *
   * This is the main entry point. It uses the full Golub-Welsch algorithm
   * (eigen-decomposition of the Jacobi matrix, see GolubWelschNodesWeights):
   * nodes are the eigenvalues, weights come from the eigenvectors. This
   * supersedes the previous "Golub-Welsch seed + Newton-Raphson refinement"
   * approach, which moved the nodes by only ~machine epsilon yet overflowed in
   * its weight computation for high-order Laguerre/Hermite. The eigenvector
   * weights are both more accurate and robust to large n.
   *
   * @param n Number of quadrature nodes
   * @param poly Pointer to the orthogonal polynomial object
   * @return Pair of vectors: first contains nodes, second contains weights
   */
  std::pair<std::vector<double>,std::vector<double>> GaussianQuadrature(int n, const OrthogonalPolynomial* poly);


  // ===========================================================================
  //  Inline implementations (header-only)
  // ===========================================================================

  namespace detail {

    /// sqrt(a^2 + b^2) computed so as to avoid spurious overflow/underflow.
    inline double pythag(double a, double b) {
      double absa = std::fabs(a), absb = std::fabs(b);
      if (absa > absb) { double t = absb / absa; return absa * std::sqrt(1.0 + t * t); }
      if (absb == 0.0) return 0.0;
      double t = absa / absb; return absb * std::sqrt(1.0 + t * t);
    }

    /// |a| with the sign of b.
    inline double sign_of(double a, double b) { return b >= 0.0 ? std::fabs(a) : -std::fabs(a); }

    // Make a rule for an even weight function exactly symmetric about the
    // origin. The nodes must be sorted in ascending order. The eigensolver
    // returns such a rule symmetric only to round-off (~1e-15, with a middle
    // node of ~1e-16 instead of 0 for odd n). Replacing each mirror pair by its
    // average gives x_k == -x_{n-1-k} and w_k == w_{n-1-k} bitwise, and since
    // it is an average it does not increase the error.
    inline void Symmetrize(std::vector<double>& x, std::vector<double>& w) {
      const int n = static_cast<int>(x.size());
      for (int lo = 0, hi = n - 1; lo < hi; ++lo, --hi) {
        const double xm = 0.5 * (x[hi] - x[lo]);
        x[lo] = -xm;
        x[hi] =  xm;
        const double wm = 0.5 * (w[lo] + w[hi]);
        w[lo] = wm;
        w[hi] = wm;
      }
      if (n % 2 == 1) x[n / 2] = 0.0;
    }

    // Implicit-shift QL eigensolver for a symmetric tridiagonal matrix, tracking
    // only the FIRST ROW of the eigenvector matrix -- all the Golub-Welsch weight
    // formula needs (w_k = mu0 * v_k[0]^2). This is the classic Golub-Welsch
    // device (cf. Golub & Welsch 1969; Numerical Recipes "tqli"; Gautschi OPQ).
    // Working directly on the (d, e) bands makes it O(n^2) time and O(n) memory.
    //
    //  d[0..n-1] : diagonal      (in) -> eigenvalues, unsorted (out)
    //  e[0..n-2] : subdiagonal, e[i] couples d[i] and d[i+1]; e[n-1] = 0 (in) -> destroyed
    //  z[0..n-1] : first eigenvector-row accumulator, init (1,0,...,0) (in)
    //              -> z[k] = first component of the eigenvector for d[k] (out)
    inline void TridiagonalQL(std::vector<double>& d, std::vector<double>& e, std::vector<double>& z) {
      const int n = static_cast<int>(d.size());
      const int MAX_ITER = 50;
      const double eps = std::numeric_limits<double>::epsilon();

      for (int l = 0; l < n; ++l) {
        int iter = 0;
        int m;
        do {
          // Look for a small subdiagonal element to split the matrix.
          for (m = l; m < n - 1; ++m) {
            double dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
            if (std::fabs(e[m]) <= eps * dd) break;
          }
          if (m != l) {
            if (iter++ == MAX_ITER)
              throw std::runtime_error("TridiagonalQL: too many iterations for eigenvalue " + std::to_string(l));
            // Form the Wilkinson shift.
            double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            double r = pythag(g, 1.0);
            g = d[m] - d[l] + e[l] / (g + sign_of(r, g));
            double s = 1.0, c = 1.0, p = 0.0;
            int i;
            for (i = m - 1; i >= l; --i) {
              double f = s * e[i];
              double b = c * e[i];
              e[i + 1] = (r = pythag(f, g));
              if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
              s = f / r;
              c = g / r;
              g = d[i + 1] - p;
              r = (d[i] - g) * s + 2.0 * c * b;
              d[i + 1] = g + (p = s * r);
              g = c * r - b;
              // Apply the Givens rotation to the tracked first eigenvector row.
              f = z[i + 1];
              z[i + 1] = s * z[i] + c * f;
              z[i]     = c * z[i] - s * f;
            }
            if (r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
          }
        } while (m != l);
      }
    }

  } // namespace detail

  inline std::pair<std::vector<double>, std::vector<double>>
  GaussianQuadraturesNewton(int n, const std::vector<double> &x0, const OrthogonalPolynomial* poly,
                            const double TOL, const int MAX_ITER) {
    if (n <= 0) {
      throw std::invalid_argument("GaussianQuadraturesNewton: number of nodes n must be positive");
    }
    if (static_cast<int>(x0.size()) != n) {
      throw std::invalid_argument("Initial guess size must be equal to the number of nodes.");
    }

    std::vector<double> x(n);
    std::vector<double> w(n);

    x = x0;

    std::vector<std::vector<double>> allcoeffs(n);
    for (int i = 0; i < n; ++i) {
      allcoeffs[i] = poly->Recurrence(i);
    }

    // Newton-Raphson method to find the roots
    for (int k = 0; k < n; ++k) {
      double p0 = 0, p1 = 1;
      double dp0 = 0, dp1 = 0;
      double dx = 1.0;
      for (int iter = 0; iter < MAX_ITER; ++iter) {
        p0 = 0;
        p1 = 1;
        dp0 = 0;
        dp1 = 0;
        for(int j = 0; j < n; ++j) {
          double np = (allcoeffs[j][0] * x[k] - allcoeffs[j][1]) * p1 - allcoeffs[j][2] * p0;
          double ndp = allcoeffs[j][0] * p1
                       + (allcoeffs[j][0] * x[k] - allcoeffs[j][1]) * dp1 - allcoeffs[j][2] * dp0;
          p0 = p1;
          p1 = np;
          dp0 = dp1;
          dp1 = ndp;
        }
        if (dp1 == 0.0 || !std::isfinite(dp1)) {
          throw std::runtime_error("Newton-Raphson: zero or non-finite derivative for root " + std::to_string(k));
        }
        dx = -p1/dp1;
        x[k] += dx;
        if (std::abs(dx) < TOL) {
          break;
        }
      }
      if (std::abs(dx) > 1.e2 * TOL) {
        throw std::runtime_error("Newton-Raphson method did not converge for root " + std::to_string(k));
      }

      w[k] = allcoeffs[n-1][0] * poly->Norm(n-1) / (p0 * dp1);
    }

    return std::make_pair(x, w);
  }

  inline std::vector<double> InitialGuessLegendre(int n) {
    if (n <= 0)
      throw std::invalid_argument("InitialGuessLegendre: number of nodes n must be positive");
    std::vector<double> x0(n);
    for (int i = 0; i < n; ++i) {
      x0[n - i - 1] = std::cos(detail::pi * (4 * i + 3) / (4 * n + 2));
    }
    return x0;
  }

  namespace detail {
    // p_n(z) and p_n'(z) of the family, up to a common positive factor, from
    // the three-term recurrence. The rescaling keeps high-order Hermite and
    // Laguerre values from overflowing; it changes neither signs nor p/p'.
    inline void EvalPolyAndDerivative(int n, double z, const OrthogonalPolynomial* poly,
                                      double& p, double& dp) {
      double p0 = 0., p1 = 1., dp0 = 0., dp1 = 0.;
      for (int j = 0; j < n; ++j) {
        const std::vector<double> c = poly->Recurrence(j);
        const double np  = (c[0] * z - c[1]) * p1 - c[2] * p0;
        const double ndp = c[0] * p1 + (c[0] * z - c[1]) * dp1 - c[2] * dp0;
        p0 = p1; p1 = np; dp0 = dp1; dp1 = ndp;
        if (std::max(std::abs(p1), std::abs(dp1)) > 1e150) {
          p0 *= 1e-150; p1 *= 1e-150; dp0 *= 1e-150; dp1 *= 1e-150;
        }
      }
      p = p1; dp = dp1;
    }

    // Newton's method safeguarded by bisection on a bracket [a, b] with
    // p_n(a) and p_n(b) of opposite signs, so it cannot leave the bracket and
    // always returns the root inside it.
    inline double SafeNewtonRoot(int n, double a, double b, double z,
                                 const OrthogonalPolynomial* poly, const int MAX_ITER = 200) {
      double pa, dpa, pz, dpz;
      if (a > b) std::swap(a, b);
      EvalPolyAndDerivative(n, a, poly, pa, dpa);
      if (!(z > a && z < b)) z = 0.5 * (a + b);
      for (int iter = 0; iter < MAX_ITER; ++iter) {
        EvalPolyAndDerivative(n, z, poly, pz, dpz);
        if (pz == 0.0) return z;
        if ((pz < 0.) == (pa < 0.)) { a = z; pa = pz; }
        else b = z;
        double znew = (dpz != 0.0) ? z - pz / dpz : 0.5 * (a + b);
        if (!(znew > a && znew < b)) znew = 0.5 * (a + b);
        if (std::abs(znew - z) <= 1e-15 * std::max(1., std::abs(znew)) || b - a <= 1e-15 * std::max(1., std::abs(z)))
          return znew;
        z = znew;
      }
      return z;
    }

    // The root next to a known root xprev in the direction of the estimate
    // zguess. Steps away from xprev by a quarter of the estimated spacing until
    // p_n changes sign, which brackets exactly the neighbouring root as long as
    // the step is below the root spacing, then refines it inside the bracket.
    inline double NextRoot(int n, double xprev, double zguess, const OrthogonalPolynomial* poly) {
      const double d = zguess - xprev;
      double a = xprev + 1e-6 * d, pa, dpa, pb, dpb;
      EvalPolyAndDerivative(n, a, poly, pa, dpa);
      double b = a;
      for (int k = 0; k < 10000; ++k) {
        const double bn = b + 0.25 * d;
        EvalPolyAndDerivative(n, bn, poly, pb, dpb);
        if ((pb < 0.) != (pa < 0.) || pb == 0.)
          return SafeNewtonRoot(n, b, bn, zguess, poly);
        b = bn;
      }
      throw std::runtime_error("NextRoot: could not bracket the next root");
    }

    // Plain Newton from a guess, for the first (outermost) root.
    inline double NewtonRoot(int n, double z, const OrthogonalPolynomial* poly, const int MAX_ITER = 200) {
      for (int iter = 0; iter < MAX_ITER; ++iter) {
        double p, dp;
        EvalPolyAndDerivative(n, z, poly, p, dp);
        if (dp == 0.0 || !std::isfinite(dp) || !std::isfinite(p))
          throw std::runtime_error("NewtonRoot: zero or non-finite derivative");
        const double dx = -p / dp;
        z += dx;
        if (std::abs(dx) <= 1e-15 * std::max(1., std::abs(z))) return z;
      }
      throw std::runtime_error("NewtonRoot: Newton-Raphson did not converge");
    }
  } // namespace detail

  // Seeds for GaussianQuadraturesNewton. The asymptotic formulas used upstream
  // (GaussQuad v1.1.0) do not isolate the roots: from them Newton converges to
  // duplicate roots (Hermite already at n = 4, Laguerre at n = 32), giving an
  // invalid rule. Here each root is guessed from the previously found ones with
  // the empirical formulas of Numerical Recipes (gaulag, gauher; after Stroud and
  // Secrest). Each later root is bracketed next to the previous one by a sign
  // change of p_n and refined by Newton safeguarded with bisection, so no root
  // can be found twice. No Golub-Welsch step is involved.
  inline std::vector<double> InitialGuessLaguerre(int n) {
    if (n <= 0)
      throw std::invalid_argument("InitialGuessLaguerre: number of nodes n must be positive");
    LaguerrePolynomial poly;
    std::vector<double> x(n);
    double z = 0.;
    for (int i = 0; i < n; ++i) {
      if (i == 0)
        z = 3. / (1. + 2.4 * n);
      else if (i == 1)
        z += 15. / (1. + 2.5 * n);
      else {
        const double ai = i - 1.;
        z += (1. + 2.55 * ai) / (1.9 * ai) * (z - x[i - 2]);
      }
      if (i > 0 && !(z > x[i - 1]))  // the estimate must lie above the last root
        z = x[i - 1] + (i > 1 ? x[i - 1] - x[i - 2] : x[i - 1]);
      z = (i == 0) ? detail::NewtonRoot(n, z, &poly) : detail::NextRoot(n, x[i - 1], z, &poly);
      x[i] = z;
    }
    return x;
  }

  inline std::vector<double> InitialGuessHermite(int n) {
    if (n <= 0)
      throw std::invalid_argument("InitialGuessHermite: number of nodes n must be positive");
    HermitePolynomial poly;
    std::vector<double> x(n);
    const int m = (n + 1) / 2;
    double z = 0.;
    // Positive roots from the largest down; the rest follow by symmetry.
    for (int i = 0; i < m; ++i) {
      if (i == 0)
        z = std::sqrt(2. * n + 1.) - 1.85575 * std::pow(2. * n + 1., -0.16667);
      else if (i == 1)
        z -= 1.14 * std::pow(static_cast<double>(n), 0.426) / z;
      else if (i == 2)
        z = 1.86 * z - 0.86 * x[n - 1];
      else if (i == 3)
        z = 1.91 * z - 0.91 * x[n - 2];
      else
        z = 2. * z - x[n - i + 1];
      if (n % 2 == 1 && i == m - 1) break;  // the middle root is 0
      if (i > 0 && !(z < x[n - i]))  // the estimate must lie below the last root
        z = x[n - i] - (i > 1 ? x[n - i + 1] - x[n - i] : 0.5 * x[n - i]);
      z = (i == 0) ? detail::NewtonRoot(n, z, &poly) : detail::NextRoot(n, x[n - i], z, &poly);
      x[n - 1 - i] = z;
      x[i] = -z;
    }
    if (n % 2 == 1) x[n / 2] = 0.;
    return x;
  }

  inline std::pair<std::vector<double>, std::vector<double>> GaussianQuadraturesLegendre(int n) {
    auto x0 = InitialGuessLegendre(n);
    LegendrePolynomial poly;
    return GaussianQuadraturesNewton(n, x0, &poly);
  }

  inline std::pair<std::vector<double>, std::vector<double>>
  GolubWelschNodesWeights(int n, const OrthogonalPolynomial *poly) {
    if (n <= 0)
      throw std::invalid_argument("GolubWelschNodesWeights: number of nodes n must be positive");
    // Build the (d, e) bands of the symmetric tridiagonal Jacobi matrix directly:
    //   diagonal      d[i]   = a_i / k_i
    //   subdiagonal   e[i-1] = sqrt(b_i / (k_{i-1} k_i))   (couples d[i-1], d[i])
    std::vector<double> d(n), e(n, 0.0), z(n, 0.0);
    std::vector<double> rec = poly->Recurrence(0);
    d[0] = rec[1] / rec[0];
    z[0] = 1.0; // first row of the identity eigenvector accumulator
    for (int i = 1; i < n; ++i) {
      std::vector<double> recPrev = rec;        // Recurrence(i-1)
      rec = poly->Recurrence(i);                // Recurrence(i)
      d[i] = rec[1] / rec[0];
      e[i - 1] = std::sqrt(rec[2] / recPrev[0] / rec[0]);
    }

    const double mu0 = poly->Norm(0);
    detail::TridiagonalQL(d, e, z);

    // Pair (node, weight) and sort by node (ascending).
    std::vector<std::pair<double, double>> nw(n);
    for (int k = 0; k < n; ++k) nw[k] = std::make_pair(d[k], mu0 * z[k] * z[k]);
    std::sort(nw.begin(), nw.end(),
              [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                return a.first < b.first;
              });

    std::vector<double> x(n), w(n);
    for (int k = 0; k < n; ++k) { x[k] = nw[k].first; w[k] = nw[k].second; }

    // Even weight -> symmetric rule; impose it exactly (see detail::Symmetrize).
    if (poly->HasEvenWeight()) detail::Symmetrize(x, w);

    return std::make_pair(x, w);
  }

  inline std::vector<double> NodesGolubWelsch(int n, const OrthogonalPolynomial *poly) {
    // Nodes are the eigenvalues; reuse the full Golub-Welsch routine and drop
    // the weights. (The extra first-eigenvector-row tracking is only O(n).)
    return GolubWelschNodesWeights(n, poly).first;
  }

  inline std::pair<std::vector<double>, std::vector<double>> GaussianQuadrature(int n, const OrthogonalPolynomial *poly) {
    // Full Golub-Welsch: nodes from eigenvalues, weights from eigenvectors.
    // No Newton-Raphson refinement -- it moved the nodes by only ~machine
    // epsilon (the tridiagonal eigenvalues are already accurate) while its
    // polynomial-evaluation weight formula overflowed for high-order
    // Laguerre/Hermite. See GolubWelschNodesWeights.
    return GolubWelschNodesWeights(n, poly);
  }
} // namespace NumericalIntegration
} // namespace thermalfist

#endif