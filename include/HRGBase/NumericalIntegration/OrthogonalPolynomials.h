/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef NUMERICALINTEGRATION_ORTHOGONALPOLYNOMIALS_H
#define NUMERICALINTEGRATION_ORTHOGONALPOLYNOMIALS_H

/**
 * \file OrthogonalPolynomials.h
 *
 * Orthogonal polynomial families used to build Gaussian quadrature rules.
 *
 * Each family is characterized by the three-term recurrence
 *   p_{n+1}(x) = (k_n x - a_n) p_n(x) - b_n p_{n-1}(x),
 * a weight function w(x), and a normalization constant h_n = \int w p_n^2 dx.
 * This is a header-only, dependency-free component (C++ standard library only).
 */

#include <vector>
#include <cmath>
#include <stdexcept>

namespace thermalfist {

  namespace NumericalIntegration {

    namespace detail {
      // Local pi constant. M_PI is not standard C++ and is unavailable on some
      // toolchains (e.g. MSVC without _USE_MATH_DEFINES), so it is not used here.
      constexpr double pi = 3.141592653589793238462643383279502884;
    }

    /// \brief Abstract base class for an orthogonal polynomial family.
    class OrthogonalPolynomial {
    public:
      virtual ~OrthogonalPolynomial() = default;

      /// Recurrence coefficients [k_n, a_n, b_n] for
      /// p_{n+1} = (k_n x - a_n) p_n - b_n p_{n-1}.
      virtual std::vector<double> Recurrence(int n) const = 0;

      /// Squared norm h_n = \int w(x) p_n(x)^2 dx.
      virtual double Norm(int n) const = 0;

      /// Weight function w(x).
      virtual double Weight(double x) const = 0;
    };

    /// \brief Legendre polynomials: interval [-1,1], weight w(x) = 1.
    class LegendrePolynomial : public OrthogonalPolynomial {
    public:
      std::vector<double> Recurrence(int n) const override {
        return { (2.*n + 1.) / (n + 1.), 0., n / (n + 1.) };
      }
      double Norm(int n) const override { return 2. / (2.*n + 1.); }
      double Weight(double /*x*/) const override { return 1.; }
    };

    /// \brief Laguerre polynomials: interval [0,inf), weight w(x) = exp(-x).
    class LaguerrePolynomial : public OrthogonalPolynomial {
    public:
      std::vector<double> Recurrence(int n) const override {
        return { -1. / (n + 1.), -(2.*n + 1.) / (n + 1.), n / (n + 1.) };
      }
      double Norm(int /*n*/) const override { return 1.; }
      double Weight(double x) const override { return exp(-x); }
    };

    /// \brief Generalized Laguerre polynomials: weight w(x) = x^alpha exp(-x), alpha > -1.
    class GeneralizedLaguerrePolynomial : public LaguerrePolynomial {
      double m_alpha;
    public:
      GeneralizedLaguerrePolynomial(double alpha) : m_alpha(alpha) {
        if (alpha <= -1.)
          throw std::invalid_argument("GeneralizedLaguerrePolynomial: alpha must be > -1");
      }
      std::vector<double> Recurrence(int n) const override {
        return { -1. / (n + 1.), -(2.*n + 1. + m_alpha) / (n + 1.), (n + m_alpha) / (n + 1.) };
      }
      double Norm(int n) const override {
        return tgamma(n + m_alpha + 1.) / tgamma(n + 1.);
      }
      double Weight(double x) const override { return exp(-x) * pow(x, m_alpha); }
    };

    /// \brief Hermite polynomials: interval (-inf,inf), weight w(x) = exp(-x^2).
    class HermitePolynomial : public OrthogonalPolynomial {
    public:
      std::vector<double> Recurrence(int n) const override {
        return { 2., 0., 2.*n };
      }
      double Norm(int n) const override {
        double ret = sqrt(detail::pi);
        for (int i = 1; i <= n; i++) ret *= 2. * i;
        return ret;
      }
      double Weight(double x) const override { return exp(-x * x); }
    };

    /// \brief Jacobi polynomials: interval [-1,1], weight w(x) = (1-x)^alpha (1+x)^beta.
    class JacobiPolynomial : public OrthogonalPolynomial {
      double m_alpha, m_beta;
    public:
      JacobiPolynomial(double alpha, double beta) : m_alpha(alpha), m_beta(beta) {
        if (alpha <= -1. || beta <= -1.)
          throw std::invalid_argument("JacobiPolynomial: alpha and beta must be > -1");
      }
      std::vector<double> Recurrence(int n) const override {
        // The general formula has a removable 0/0 at n == 0 when alpha+beta == 0
        // (e.g. Legendre (0,0)) or alpha+beta == -1 (e.g. Chebyshev (-1/2,-1/2)).
        // Use the closed form P_1(x) = (a+b+2)/2 x + (a-b)/2 directly.
        if (n == 0) {
          return { 0.5 * (m_alpha + m_beta + 2.),   // k_0
                   -0.5 * (m_alpha - m_beta),         // a_0
                   0. };                              // b_0
        }
        double a = n + 1. + m_alpha;
        double b = n + 1. + m_beta;
        double c = a + b;
        double den = 2. * (n + 1.) * (c - n - 1.) * (c - 2);
        double an = (c - 1.) * c * (c - 2.) / den;
        double aln = -(a - b) * (c - 2.*n - 2.) / den;
        double betn = 2. * (a - 1.) * (b - 1.) * c / den;
        return { an, aln, betn };
      }
      double Norm(int n) const override {
        // h_0 = 2^(a+b+1) Gamma(a+1) Gamma(b+1) / Gamma(a+b+2); the closed form
        // avoids the (2n+a+b+1)/Gamma(n+a+b+1) singularity at n == 0, a+b == -1.
        if (n == 0) {
          return pow(2., m_alpha + m_beta + 1.)
            * tgamma(m_alpha + 1.) * tgamma(m_beta + 1.)
            / tgamma(m_alpha + m_beta + 2.);
        }
        return pow(2., m_alpha + m_beta + 1.) / (2.*n + m_alpha + m_beta + 1.)
          * tgamma(n + m_alpha + 1.) / tgamma(n + 1.)
          * tgamma(n + m_beta + 1.) / tgamma(n + m_alpha + m_beta + 1.);
      }
      double Weight(double x) const override {
        return pow(1. - x, m_alpha) * pow(1. + x, m_beta);
      }
    };

  } // namespace NumericalIntegration

} // namespace thermalfist

#endif
