/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef NUMERICALINTEGRATION_GAUSSIANQUADRATURE_H
#define NUMERICALINTEGRATION_GAUSSIANQUADRATURE_H

/**
 * \file GaussianQuadrature.h
 *
 * Computation of Gaussian quadrature nodes and weights for the orthogonal
 * polynomial families in OrthogonalPolynomials.h.
 *
 * The default method is the full Golub-Welsch algorithm implemented with an
 * in-house implicit-shift tridiagonal QL eigensolver that tracks only the first
 * row of the eigenvector matrix (all the weight formula needs). It is O(n^2)
 * time / O(n) memory, dependency-free, and accurate to machine precision for
 * high orders. A Newton-Raphson path with asymptotic seeds is also provided.
 *
 * This is a header-only component (C++ standard library only).
 */

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

      // Implicit-shift QL eigensolver for a symmetric tridiagonal matrix, tracking
      // only the FIRST ROW of the eigenvector matrix -- all the Golub-Welsch weight
      // formula needs (w_k = mu0 * v_k[0]^2). Classic Golub-Welsch device
      // (cf. Golub & Welsch 1969; Numerical Recipes "tqli"; Gautschi OPQ).
      // Working directly on the (d, e) bands makes it O(n^2) time / O(n) memory.
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

    /**
     * \brief Compute nodes and weights via the full Golub-Welsch algorithm.
     *
     * Builds the diagonal/subdiagonal bands of the symmetric tridiagonal Jacobi
     * matrix from the recurrence coefficients and solves it with the in-house QL
     * iteration. Nodes are the eigenvalues; weights are w_k = mu0 * v_k[0]^2 with
     * mu0 = Norm(0) the integral of the weight function. O(n^2) time / O(n)
     * memory; no polynomial evaluation, so no overflow at high order, and the sum
     * rule holds to machine precision by orthonormality of the eigenvectors.
     *
     * \param n    Number of nodes
     * \param poly Orthogonal polynomial family
     * \return Pair (nodes [sorted ascending], weights)
     */
    inline std::pair<std::vector<double>, std::vector<double>>
    GolubWelschNodesWeights(int n, const OrthogonalPolynomial* poly) {
      if (n <= 0)
        throw std::invalid_argument("GolubWelschNodesWeights: number of nodes n must be positive");
      // Bands of the symmetric tridiagonal Jacobi matrix:
      //   diagonal      d[i]   = a_i / k_i
      //   subdiagonal   e[i-1] = sqrt(b_i / (k_{i-1} k_i))   (couples d[i-1], d[i])
      std::vector<double> d(n), e(n, 0.0), z(n, 0.0);
      std::vector<double> rec = poly->Recurrence(0);
      d[0] = rec[1] / rec[0];
      z[0] = 1.0;
      for (int i = 1; i < n; ++i) {
        std::vector<double> recPrev = rec;
        rec = poly->Recurrence(i);
        d[i] = rec[1] / rec[0];
        e[i - 1] = std::sqrt(rec[2] / recPrev[0] / rec[0]);
      }

      const double mu0 = poly->Norm(0);
      detail::TridiagonalQL(d, e, z);

      std::vector<std::pair<double, double>> nw(n);
      for (int k = 0; k < n; ++k) nw[k] = std::make_pair(d[k], mu0 * z[k] * z[k]);
      std::sort(nw.begin(), nw.end(),
                [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                  return a.first < b.first;
                });

      std::vector<double> x(n), w(n);
      for (int k = 0; k < n; ++k) { x[k] = nw[k].first; w[k] = nw[k].second; }
      return std::make_pair(x, w);
    }

    /// \brief Quadrature nodes only (eigenvalues of the Jacobi matrix), sorted.
    inline std::vector<double> NodesGolubWelsch(int n, const OrthogonalPolynomial* poly) {
      return GolubWelschNodesWeights(n, poly).first;
    }

    /**
     * \brief Main entry point: nodes and weights for any orthogonal family.
     *
     * Uses the full Golub-Welsch algorithm (see GolubWelschNodesWeights).
     */
    inline std::pair<std::vector<double>, std::vector<double>>
    GaussianQuadrature(int n, const OrthogonalPolynomial* poly) {
      return GolubWelschNodesWeights(n, poly);
    }

    /**
     * \brief Refine quadrature nodes with the Newton-Raphson method.
     *
     * Refines initial guesses x0 to the zeros of the degree-n orthogonal
     * polynomial and computes the weights from the polynomial derivative.
     * Note: the polynomial evaluation can overflow for high-order Laguerre/Hermite;
     * GolubWelschNodesWeights is preferred for those. Kept for the Legendre path.
     *
     * \throws std::invalid_argument if x0.size() != n
     * \throws std::runtime_error    if Newton-Raphson fails to converge
     */
    inline std::pair<std::vector<double>, std::vector<double>>
    GaussianQuadraturesNewton(int n, const std::vector<double>& x0, const OrthogonalPolynomial* poly,
                              const double TOL = 1e-15, const int MAX_ITER = 100) {
      if (n <= 0)
        throw std::invalid_argument("GaussianQuadraturesNewton: number of nodes n must be positive");
      if (static_cast<int>(x0.size()) != n)
        throw std::invalid_argument("Initial guess size must be equal to the number of nodes.");

      std::vector<double> x(x0), w(n);

      std::vector<std::vector<double>> allcoeffs(n);
      for (int i = 0; i < n; ++i) allcoeffs[i] = poly->Recurrence(i);

      for (int k = 0; k < n; ++k) {
        double p0 = 0, p1 = 1, dp0 = 0, dp1 = 0, dx = 1.0;
        for (int iter = 0; iter < MAX_ITER; ++iter) {
          p0 = 0; p1 = 1; dp0 = 0; dp1 = 0;
          for (int j = 0; j < n; ++j) {
            double np = (allcoeffs[j][0] * x[k] - allcoeffs[j][1]) * p1 - allcoeffs[j][2] * p0;
            double ndp = allcoeffs[j][0] * p1
                         + (allcoeffs[j][0] * x[k] - allcoeffs[j][1]) * dp1 - allcoeffs[j][2] * dp0;
            p0 = p1; p1 = np; dp0 = dp1; dp1 = ndp;
          }
          dx = -p1 / dp1;
          x[k] += dx;
          if (std::abs(dx) < TOL) break;
        }
        if (std::abs(dx) > 1.e2 * TOL)
          throw std::runtime_error("Newton-Raphson method did not converge for root " + std::to_string(k));
        w[k] = allcoeffs[n - 1][0] * poly->Norm(n - 1) / (p0 * dp1);
      }
      return std::make_pair(x, w);
    }

    /// \brief Asymptotic initial guesses for Legendre zeros.
    inline std::vector<double> InitialGuessLegendre(int n) {
      if (n <= 0)
        throw std::invalid_argument("InitialGuessLegendre: number of nodes n must be positive");
      std::vector<double> x0(n);
      for (int i = 0; i < n; ++i)
        x0[n - i - 1] = std::cos(detail::pi * (4 * i + 3) / (4 * n + 2));
      return x0;
    }

    /// \brief Gauss-Legendre nodes and weights via Newton refinement of asymptotic seeds.
    inline std::pair<std::vector<double>, std::vector<double>> GaussianQuadraturesLegendre(int n) {
      LegendrePolynomial poly;
      return GaussianQuadraturesNewton(n, InitialGuessLegendre(n), &poly);
    }

  } // namespace NumericalIntegration

} // namespace thermalfist

#endif
