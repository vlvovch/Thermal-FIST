/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 *
 * Vendored from GaussQuad (https://github.com/vlvovch/gaussquad-cpp), v1.1.0,
 * commit a1dc387. Local changes: namespace gaussquad -> thermalfist::NumericalIntegration,
 * include paths and include guards. Keep in sync with upstream.
 */
/**
 * @file OrthogonalPolynomials.h
 * @brief Header file defining orthogonal polynomial families for Gaussian quadrature.
 *
 * This file contains the abstract base class OrthogonalPolynomial and concrete
 * implementations for various families of orthogonal polynomials commonly used
 * in numerical integration via Gaussian quadrature.
 *
 * Each orthogonal polynomial family is characterized by:
 * - A three-term recurrence relation
 * - A weight function w(x)
 * - A normalization constant h_n
 * - An interval of orthogonality
 *
 * The recurrence relation has the form:
 * p_{n+1}(x) = (k_n x - a_n)p_n(x) - b_n p_{n-1}(x)
 * where k_n, a_n, b_n are the recurrence coefficients.
 */

#ifndef NUMERICALINTEGRATION_ORTHOGONALPOLYNOMIALS_H
#define NUMERICALINTEGRATION_ORTHOGONALPOLYNOMIALS_H
#include <vector>
#include <cmath>
#include <stdexcept>

namespace thermalfist {
namespace NumericalIntegration {

  namespace detail {
    // Local pi constant. M_PI is not standard C++ and is unavailable on some
    // toolchains (e.g. MSVC without _USE_MATH_DEFINES), so it is not used here.
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double ln2 = 0.693147180559945309417232121458176568;

    // Remainder of Stirling's series,
    //   lgamma(x) = (x - 1/2) log(x) - x + log(2 pi)/2 + StirlingRemainder(x),
    // summed through the x^-13 term. For x >= 10 the truncation error is
    // below 3e-17.
    inline double StirlingRemainder(double x) {
      const double r = 1. / x, r2 = r * r;
      return r * (1. / 12. + r2 * (-1. / 360. + r2 * (1. / 1260. + r2 * (-1. / 1680.
           + r2 * (1. / 1188. + r2 * (-691. / 360360. + r2 * (1. / 156.)))))));
    }

    // mu0 = 2^(a+b+1) Gamma(a+1) Gamma(b+1) / Gamma(a+b+2), the integral of
    // the Jacobi weight (1-x)^a (1+x)^b over [-1, 1], for a, b > -1.
    //
    // With p = min(a,b)+1 and q = max(a,b)+1:
    //  - p, q < 10: the direct product of tgamma values. The arguments are
    //    small, where tgamma is accurate to a few ulps and cannot overflow.
    //  - otherwise, Stirling's series with the large terms cancelled
    //    analytically. A plain sum of lgamma values would lose accuracy for
    //    large exponents: the terms grow like a log(a) and cancel almost
    //    completely (for a = b, mu0 ~ sqrt(pi/a)), leaving their rounding
    //    errors in the result (~1e-12 relative at a = b = 1000, ~1% at
    //    a = b = 1e12). The direct tgamma product, besides overflowing, is
    //    also less accurate here, since tgamma itself errs by up to a few
    //    hundred ulps at large non-integer arguments (glibc 2.39).
    //
    //    First, while p < 10, p is raised by one using
    //      mu0(p, q) = mu0(p+1, q) * (p+q) / (2p),
    //    which follows from B(p, q) = B(p+1, q) (p+q)/p; each factor exceeds 1
    //    because q > p, and is exact to half an ulp. Then, with s = p + q,
    //    d = (q - p)/s and R = StirlingRemainder,
    //      log mu0 = log(2 pi / s)/2 + (s-1)/2 log(1 - d^2) + s d atanh(d)
    //                + R(p) + R(q) - R(s),
    //    the (s-1) log 2 and (x - 1/2) log x terms of Stirling's series having
    //    cancelled exactly. This is accurate to a few ulps of log mu0, which is
    //    the best that exp() of a computed logarithm can give.
    inline double JacobiWeightIntegral(double a, double b) {
      double p = (a < b ? a : b) + 1.;
      const double q = (a < b ? b : a) + 1.;
      const double big = 10.;
      if (q < big)
        return std::pow(2., a + b + 1.) * std::tgamma(a + 1.) * std::tgamma(b + 1.)
             / std::tgamma(a + b + 2.);
      double scale = 1.;
      int k = 0;
      for (; p < big; p += 1., ++k) scale *= (p + q) / (2. * p);
      // h = s/2, formed from halves so that it stays finite when a + b
      // overflows (mu0 ~ sqrt(pi/a) is still representable there), and used
      // in place of s throughout: (s-1)/2 = h - 1/2 and log(2 pi / s) =
      // log(pi / h).
      const double h = 0.5 * a + 0.5 * b + 1. + 0.5 * k;
      // q - p, from |a - b| rather than from p and q, so that it stays exact
      // when the exponents are huge and nearly equal.
      const double diff = std::fabs(a - b) - k;
      const double d = 0.5 * diff / h;
      // For small d, log1p and atanh keep the O(s d^2) terms accurate. For
      // d >= 1/2, using 1 - d^2 = (p/h)(q/h) and atanh(d) = log(q/p)/2, the
      // body regroups exactly as (p - 1/2) log(p/h) + (q - 1/2) log(q/h),
      // which avoids computing 1 - d with cancellation. In that form each
      // term stays finite for any p, q up to the double range (|x log(x/h)|
      // <= h/e for x <= h, and q < 2h), so no intermediate overflows.
      const double body = d < 0.5
          ? (h - 0.5) * std::log1p(-d * d) + diff * std::atanh(d)
          : (p - 0.5) * std::log(p / h) + (q - 0.5) * std::log(q / h);
      const double logMu0 = 0.5 * std::log(pi / h) + body
                          + StirlingRemainder(p) + StirlingRemainder(q) - StirlingRemainder(2. * h);
      return std::exp(logMu0) * scale;
    }
  }

  /**
   * @brief Abstract base class for orthogonal polynomials.
   *
   * This class defines the interface for orthogonal polynomial families.
   * Each derived class must implement the recurrence relation coefficients,
   * normalization constant, and weight function.
   *
   * The orthogonal polynomials satisfy the orthogonality condition:
   * ∫[a,b] w(x) p_m(x) p_n(x) dx = h_n δ_{mn}
   *
   * where [a,b] is the interval of orthogonality, w(x) is the weight function,
   * h_n is the normalization constant, and δ_{mn} is the Kronecker delta.
   */

  class OrthogonalPolynomial {
  public:
    /**
     * @brief Virtual destructor so derived instances can be safely deleted
     *        through a base-class pointer.
     */
    virtual ~OrthogonalPolynomial() = default;

    /**
     * @brief Get the recurrence relation coefficients for polynomial of order n.
     *
     * Returns the coefficients [k_n, a_n, b_n] for the three-term recurrence relation:
     * p_{n+1}(x) = (k_n x - a_n)p_n(x) - b_n p_{n-1}(x)
     *
     * @param n Order of the polynomial (n ≥ 0)
     * @return Vector containing [k_n, a_n, b_n] coefficients
     */
    virtual std::vector<double> Recurrence(int n) const = 0;

    /**
     * @brief Get the normalization constant for polynomial of order n.
     *
     * Returns h_n = ∫[a,b] w(x) p_n(x)² dx, the squared norm of the polynomial.
     * This is used in computing the quadrature weights.
     *
     * @param n Order of the polynomial (n ≥ 0)
     * @return Normalization constant h_n
     */
    virtual double Norm(int n) const = 0;

    /**
     * @brief Evaluate the weight function at point x.
     *
     * Returns w(x), the weight function associated with this polynomial family.
     * The weight function defines the measure for the inner product and integration.
     *
     * @param x Point at which to evaluate the weight function
     * @return Value of the weight function w(x)
     */
    virtual double Weight(double x) const = 0;

    /**
     * @brief Whether the weight function is even on an interval symmetric about zero.
     *
     * When this returns true, i.e. w(-x) = w(x) on (-a, a), the Gauss rule is
     * symmetric about the origin, and the Golub-Welsch solver imposes that
     * symmetry exactly: x_k == -x_{n-1-k} and w_k == w_{n-1-k} bitwise, with
     * the middle node exactly 0 for odd n. Otherwise the eigensolver would
     * leave the symmetry broken at the round-off level.
     *
     * The default is false, so families that do not override this method are
     * left untouched.
     *
     * @return true if the weight function is even about x = 0
     */
    virtual bool HasEvenWeight() const { return false; }
  };

  /**
   * @brief Implementation of Legendre polynomials.
   *
   * Legendre polynomials are orthogonal on the interval [-1, 1] with weight function w(x) = 1.
   * They are solutions to Legendre's differential equation and are widely used in
   * physics and engineering applications.
   *
   * Properties:
   * - Interval: [-1, 1]
   * - Weight function: w(x) = 1
   * - Recurrence: (n+1) P_{n+1}(x) = ((2n+1)x - 0)P_n(x) - n P_{n-1}(x)
   * - Normalization: h_n = 2/(2n+1)
   */
  class LegendrePolynomial : public OrthogonalPolynomial {
  public:
    std::vector<double> Recurrence(int n) const override {
      return { (2.*n + 1.)/(n+1.),
                0.,
                n/(n+1.) };
    }

    double Norm(int n) const override {
      return 2. / (2.*n + 1.);
    }

    double Weight(double /*x*/) const override {
      return 1.;
    }

    bool HasEvenWeight() const override { return true; }
  };

  /**
   * @brief Implementation of Laguerre polynomials.
   *
   * Laguerre polynomials are orthogonal on the interval [0, ∞) with weight function w(x) = exp(-x).
   * They are solutions to Laguerre's differential equation and are commonly used in
   * quantum mechanics and probability theory.
   *
   * Properties:
   * - Interval: [0, ∞)
   * - Weight function: w(x) = exp(-x)
   * - Recurrence: (n+1) L_{n+1}(x) = (-x + (2n+1))L_n(x) - n L_{n-1}(x)
   * - Normalization: h_n = 1
   */
  class LaguerrePolynomial : public OrthogonalPolynomial {
  public:
    virtual std::vector<double> Recurrence(int n) const override {
      return { -1. / (n + 1.),
                -(2. * n + 1) / (n + 1.),
                n / (n + 1.) };
    }

    virtual double Norm(int /*n*/) const override {
      return 1.;
    }

    virtual double Weight(double x) const override {
      return std::exp(-x);
    }
  };

  /**
   * @brief Implementation of generalized Laguerre polynomials.
   *
   * Generalized Laguerre polynomials L_n^(α)(x) are orthogonal on the interval [0, ∞)
   * with weight function w(x) = x^α exp(-x), where α > -1.
   * They generalize the standard Laguerre polynomials (α = 0).
   *
   * Properties:
   * - Interval: [0, ∞)
   * - Weight function: w(x) = x^α exp(-x)
   * - Parameter: α > -1
   * - Recurrence: (n+1) L_{n+1}^(α)(x) = (-x + (2n+1+α))L_n^(α)(x) - (n+α) L_{n-1}^(α)(x)
   * - Normalization: h_n = Γ(n+α+1)/Γ(n+1)
   */
  class GeneralizedLaguerrePolynomial : public LaguerrePolynomial {
    double m_alpha;
  public:
    explicit GeneralizedLaguerrePolynomial(double alpha) : m_alpha(alpha) {
      if (alpha <= -1.)
        throw std::invalid_argument("GeneralizedLaguerrePolynomial: alpha must be > -1");
    }
    std::vector<double> Recurrence(int n) const override {
      return { -1. / (n + 1.),
               -(2. * n + 1 + m_alpha) / (n + 1.),
               (n + m_alpha) / (n + 1.) };
    }

    double Norm(int n) const override {
      // h_n = Gamma(n+alpha+1)/Gamma(n+1). Both gammas overflow to +inf once
      // n + alpha exceeds ~170, long before the ratio itself does (inf/inf =
      // NaN). The direct ratio is kept wherever it is a normal number;
      // elsewhere the log-gamma difference is used, which overflows only if
      // h_n itself is not representable. The direct ratio is the more
      // accurate of the two where both apply: the difference of two large
      // lgammas loses relative accuracy in proportion to their size (measured
      // with glibc 2.39 for alpha < 171: median 0.8-10 ulp for the ratio
      // against 2-139 ulp for the lgamma difference, worst cases 416 against
      // 1214 ulp).
      const double h = std::tgamma(n + m_alpha + 1.) / std::tgamma(n + 1.);
      if (std::isnormal(h)) return h;
      return std::exp(std::lgamma(n + m_alpha + 1.) - std::lgamma(n + 1.));
    }

    double Weight(double x) const override {
      return std::exp(-x) * std::pow(x, m_alpha);
    }
  };

  /**
   * @brief Implementation of Hermite polynomials.
   *
   * Hermite polynomials are orthogonal on the interval (-∞, ∞) with weight function w(x) = exp(-x²).
   * They are solutions to Hermite's differential equation and are fundamental in quantum mechanics,
   * particularly in the quantum harmonic oscillator.
   *
   * Properties:
   * - Interval: (-∞, ∞)
   * - Weight function: w(x) = exp(-x²)
   * - Recurrence: H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x)
   * - Normalization: h_n = √π * 2^n * n!
   */
  class HermitePolynomial : public OrthogonalPolynomial {
  public:
    std::vector<double> Recurrence(int n) const override {
      return { 2.,
                0.,
                2. * n };
    }

    double Norm(int n) const override {
      // h_n = sqrt(pi) * 2^n * n!. This grows past the double range around
      // n = 150, so the returned value is genuinely +inf for large n (only the
      // Newton weight path reads Norm(n>0); Golub-Welsch uses Norm(0) only).
      double ret = std::sqrt(detail::pi);
      for (int i = 1; i <= n; i++) {
        ret *= 2. * i;
      }
      return ret;
    }

    double Weight(double x) const override {
      return std::exp(-x * x);
    }

    bool HasEvenWeight() const override { return true; }
  };

  /**
   * @brief Implementation of Jacobi polynomials.
   *
   * Jacobi polynomials P_n^(α,β)(x) are orthogonal on the interval [-1, 1] with weight function
   * w(x) = (1-x)^α (1+x)^β, where α, β > -1. They are the most general class of classical
   * orthogonal polynomials, with Legendre polynomials as a special case (α = β = 0).
   *
   * Properties:
   * - Interval: [-1, 1]
   * - Weight function: w(x) = (1-x)^α (1+x)^β
   * - Parameters: α, β > -1
   * - Recurrence: Complex three-term recurrence relation
   * - Normalization: h_n = 2^(α+β+1) * Γ(n+α+1) * Γ(n+β+1) / ((2n+α+β+1) * Γ(n+1) * Γ(n+α+β+1))
   */
  class JacobiPolynomial : public OrthogonalPolynomial {
    double m_alpha, m_beta;
  public:
    JacobiPolynomial(double alpha, double beta) : m_alpha(alpha), m_beta(beta) {
      if (alpha <= -1. || beta <= -1.)
        throw std::invalid_argument("JacobiPolynomial: alpha and beta must be > -1");
    }
    std::vector<double> Recurrence(int n) const override {
      // The general formula below has a removable 0/0 at n == 0 whenever
      // alpha + beta == 0 (e.g. Legendre (0,0)) or alpha + beta == -1 (e.g.
      // Chebyshev (-1/2,-1/2)). Use the closed form for P_1 directly:
      //   P_1(x) = (alpha+beta+2)/2 * x + (alpha-beta)/2
      // i.e. p_1 = (k_0 x - a_0) p_0 with k_0, a_0 below and b_0 = 0.
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
      double aln = -(c - 1.) * (a - b) * (c - 2. * n - 2.) / den;
      double betn = 2. * (a - 1.) * (b - 1.) * c / den;

      return { an,
               aln,
               betn };
    }

    double Norm(int n) const override {
      // h_0 = integral of the weight = 2^(a+b+1) * B(a+1,b+1)
      //     = 2^(a+b+1) * Gamma(a+1) Gamma(b+1) / Gamma(a+b+2).
      // The general-n formula below has a (2n+a+b+1) denominator and a
      // Gamma(n+a+b+1) factor that are singular at n == 0 when a+b == -1;
      // the closed form for h_0 avoids both.
      //
      // h_0 scales every Golub-Welsch weight, and is computed by
      // detail::JacobiWeightIntegral: directly for small exponents, and
      // otherwise in a form that is accurate to a few ulps of log h_0 for any
      // exponent size.
      //
      // For n > 0, which only the Newton-Raphson weight formula reads, the
      // products of tgamma values overflow once an argument exceeds ~170,
      // although h_n itself may be of moderate size. The direct form is kept
      // where it is a normal number; otherwise the same expression is
      // evaluated through lgamma. That fallback is simpler, but loses accuracy
      // as the gamma arguments grow, roughly 1e-16 * N log(N) relative with
      // N = n + alpha + beta (~1e-12 at N ~ 2000). All gamma arguments are
      // positive for alpha, beta > -1, so lgamma is the logarithm of a
      // positive Gamma.
      const double a = m_alpha, b = m_beta;
      if (n == 0) return detail::JacobiWeightIntegral(a, b);
      const double h = std::pow(2., a + b + 1.) / (2. * n + a + b + 1.)
        * std::tgamma(n + a + 1.) / std::tgamma(n + 1.)
        * std::tgamma(n + b + 1.) / std::tgamma(n + a + b + 1.);
      if (std::isnormal(h)) return h;
      return std::exp((a + b + 1.) * detail::ln2 - std::log(2. * n + a + b + 1.)
                      + std::lgamma(n + a + 1.) - std::lgamma(n + 1.)
                      + std::lgamma(n + b + 1.) - std::lgamma(n + a + b + 1.));
    }

    double Weight(double x) const override {
      return std::pow(1. - x, m_alpha) * std::pow(1. + x, m_beta);
    }

    bool HasEvenWeight() const override { return m_alpha == m_beta; }
  };
} // namespace NumericalIntegration
} // namespace thermalfist

#endif