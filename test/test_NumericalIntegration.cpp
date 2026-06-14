/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <vector>
#include <cmath>
#include <utility>
#include "HRGBase/NumericalIntegration.h"
#include "gtest/gtest.h"

using namespace thermalfist;
namespace NI = thermalfist::NumericalIntegration;

namespace {

  // The quadrature nodes/weights used to be hard-coded tables (12 figures).
  // They are now computed by the in-house Golub-Welsch solver; these reference
  // arrays pin the computed values to the historical constants so the switch is
  // provably behavior-preserving. Tolerance 1e-10 absorbs the 12-digit rounding
  // of the references while still catching any real algorithmic regression.
  const double kTol = 1e-10;

  // ---- Reference constants (the original hard-coded tables) ----------------
  const double ref_xleg5[5] = { -0.906179845939, -0.538469310106, 0., 0.538469310106, 0.906179845939 };
  const double ref_wleg5[5] = { 0.236926885056, 0.478628670499, 0.568888888889, 0.478628670499, 0.236926885056 };

  const double ref_xleg10[10] = { -0.973906528517, -0.865063366689, -0.679409568299, -0.433395394129,
    -0.148874338982, 0.148874338982, 0.433395394129, 0.679409568299, 0.865063366689, 0.973906528517 };
  const double ref_wleg10[10] = { 0.0666713443087, 0.149451349151, 0.219086362516, 0.269266719310,
    0.295524224715, 0.295524224715, 0.269266719310, 0.219086362516, 0.149451349151, 0.0666713443087 };

  const double ref_xleg32[32] = { -0.997263861849, -0.985611511545, -0.964762255588, -0.934906075938,
    -0.896321155766, -0.849367613733, -0.794483795968, -0.732182118740, -0.663044266930, -0.587715757241,
    -0.506899908932, -0.421351276131, -0.331868602282, -0.239287362252, -0.144471961583, -0.048307665688,
    0.048307665688, 0.144471961583, 0.239287362252, 0.331868602282, 0.421351276131, 0.506899908932,
    0.587715757241, 0.663044266930, 0.732182118740, 0.794483795968, 0.849367613733, 0.896321155766,
    0.934906075938, 0.964762255588, 0.985611511545, 0.997263861849 };
  const double ref_wleg32[32] = { 0.00701861000947, 0.0162743947309, 0.0253920653093, 0.0342738629130,
    0.0428358980222, 0.0509980592624, 0.0586840934785, 0.0658222227764, 0.0723457941088, 0.0781938957871,
    0.0833119242269, 0.0876520930044, 0.0911738786958, 0.0938443990808, 0.0956387200793, 0.0965400885147,
    0.0965400885147, 0.0956387200793, 0.0938443990808, 0.0911738786958, 0.0876520930044, 0.0833119242269,
    0.0781938957871, 0.0723457941088, 0.0658222227764, 0.0586840934785, 0.0509980592624, 0.0428358980222,
    0.0342738629130, 0.0253920653093, 0.0162743947309, 0.00701861000947 };

  const double ref_xlag32[32] = { 0.0444893658333, 0.234526109520, 0.576884629302, 1.07244875382,
    1.72240877644, 2.52833670643, 3.49221327302, 4.61645676975, 5.90395850417, 7.35812673319,
    8.98294092421, 10.7830186325, 12.7636979867, 14.9311397555, 17.2924543367, 19.8558609403,
    22.6308890132, 25.6286360225, 28.8621018163, 32.3466291540, 36.1004948058, 40.1457197715,
    44.5092079958, 49.2243949873, 54.3337213334, 59.8925091621, 65.9753772879, 72.6876280907,
    80.1874469779, 88.7353404179, 98.8295428683, 111.751398098 };
  const double ref_wlag32[32] = { 0.114187105768, 0.266065216898, 0.418793137325, 0.572532846500,
    0.727648788381, 0.884536719340, 1.04361887589, 1.20534927415, 1.37022133852, 1.53877725647,
    1.71161935269, 1.88942406345, 2.07295934025, 2.26310663400, 2.46088907249, 2.66750812640,
    2.88439209292, 3.11326132704, 3.35621769260, 3.61586985648, 3.89551304495, 4.19939410471,
    4.53311497853, 4.90427028761, 5.32350097202, 5.80633321423, 6.37661467416, 7.07352658071,
    7.96769350930, 9.20504033128, 11.1630130908, 15.3901804153 };

  void expectMatch(const std::vector<double>& got, const double* ref, int n, const char* what) {
    ASSERT_EQ(static_cast<int>(got.size()), n) << what;
    for (int i = 0; i < n; ++i)
      EXPECT_NEAR(got[i], ref[i], kTol) << what << " index " << i;
  }

  // ---- Regression: computed base rules match the historical tables ---------

  TEST(NumericalIntegrationReference, Legendre5)  { expectMatch(NI::GetGaussLegendreNodes(5),   ref_xleg5, 5, "xleg5");
                                                    expectMatch(NI::GetGaussLegendreWeights(5),  ref_wleg5, 5, "wleg5"); }
  TEST(NumericalIntegrationReference, Legendre10) { expectMatch(NI::GetGaussLegendreNodes(10),  ref_xleg10, 10, "xleg10");
                                                    expectMatch(NI::GetGaussLegendreWeights(10), ref_wleg10, 10, "wleg10"); }
  TEST(NumericalIntegrationReference, Legendre32) { expectMatch(NI::GetGaussLegendreNodes(32),  ref_xleg32, 32, "xleg32");
                                                    expectMatch(NI::GetGaussLegendreWeights(32), ref_wleg32, 32, "wleg32"); }
  TEST(NumericalIntegrationReference, Laguerre32) { expectMatch(NI::GetGaussLaguerreNodes(32),  ref_xlag32, 32, "xlag32");
                                                    expectMatch(NI::GetGaussLaguerreWeights(32), ref_wlag32, 32, "wlag32 (modified)"); }

  // ---- GetCoefs* interface ------------------------------------------------

  TEST(NumericalIntegrationInterface, Legendre32MatchesReferenceOnUnitInterval) {
    std::vector<double> x, w;
    NI::GetCoefsIntegrateLegendre32(-1., 1., &x, &w);
    expectMatch(x, ref_xleg32, 32, "GetCoefsIntegrateLegendre32 nodes");
    expectMatch(w, ref_wleg32, 32, "GetCoefsIntegrateLegendre32 weights");
  }

  TEST(NumericalIntegrationInterface, Laguerre32MatchesReference) {
    std::vector<double> x, w;
    NI::GetCoefsIntegrateLaguerre32(&x, &w);
    expectMatch(x, ref_xlag32, 32, "GetCoefsIntegrateLaguerre32 nodes");
    expectMatch(w, ref_wlag32, 32, "GetCoefsIntegrateLaguerre32 weights");
  }

  TEST(NumericalIntegrationInterface, LegendreSumRules) {
    std::vector<double> x, w;
    for (int n : {5, 10, 32, 40, 64}) {
      NI::GetCoefsIntegrateLegendre(n, 0., 3., &x, &w);
      double s = 0; for (double v : w) s += v;
      EXPECT_NEAR(s, 3.0, 1e-12) << "n=" << n;        // length of [0,3]
    }
  }

  TEST(NumericalIntegrationInterface, LaguerreModifiedWeightSumRule) {
    std::vector<double> x, w;
    NI::GetCoefsIntegrateLaguerre(32, &x, &w);
    double s = 0; for (size_t i = 0; i < w.size(); ++i) s += w[i] * std::exp(-x[i]);
    EXPECT_NEAR(s, 1.0, 1e-12);                        // sum w_k e^{-x_k} == 1
  }

  // ---- Accuracy / arbitrary node count ------------------------------------

  TEST(NumericalIntegrationAccuracy, PolynomialExactness) {
    // n-point Gauss-Legendre is exact for polynomials up to degree 2n-1.
    std::vector<double> x, w;
    NI::GetCoefsIntegrateLegendre(7, -1., 1., &x, &w);
    double s = 0; for (size_t i = 0; i < x.size(); ++i) s += w[i] * std::pow(x[i], 6);
    EXPECT_NEAR(s, 2.0 / 7.0, 1e-13);                 // \int_{-1}^1 x^6 dx
  }

  TEST(NumericalIntegrationAccuracy, LaguerreCosineIntegral) {
    std::vector<double> x, w;
    NI::GetCoefsIntegrateLaguerre(32, &x, &w);
    double s = 0; for (size_t i = 0; i < x.size(); ++i) s += w[i] * std::exp(-x[i]) * std::cos(x[i]);
    EXPECT_NEAR(s, 0.5, 1e-12);                        // \int_0^inf e^-x cos x dx
  }

  TEST(NumericalIntegrationAccuracy, ArbitraryNodeCountRunge) {
    std::vector<double> x, w;
    NI::GetCoefsIntegrateLegendre(128, -1., 1., &x, &w);
    double s = 0; for (size_t i = 0; i < x.size(); ++i) s += w[i] / (1. + 25. * x[i] * x[i]);
    EXPECT_NEAR(s, 0.4 * std::atan(5.0), 1e-12);      // hard Runge integral, n=128
  }

  // ---- Input validation: non-positive node counts must throw ---------------

  TEST(NumericalIntegrationValidation, NonPositiveNodeCountThrows) {
    std::vector<double> x, w;
    NI::LegendrePolynomial leg;
    for (int n : {0, -1, -100}) {
      EXPECT_THROW(NI::GetCoefsIntegrateLegendre(n, 0., 1., &x, &w), std::invalid_argument) << "n=" << n;
      EXPECT_THROW(NI::GetCoefsIntegrateLaguerre(n, &x, &w),         std::invalid_argument) << "n=" << n;
      EXPECT_THROW(NI::GetGaussLegendreNodes(n),                     std::invalid_argument) << "n=" << n;
      EXPECT_THROW(NI::GetGaussLaguerreWeights(n),                   std::invalid_argument) << "n=" << n;
      EXPECT_THROW(NI::GaussianQuadrature(n, &leg),                  std::invalid_argument) << "n=" << n;
      EXPECT_THROW(NI::GaussianQuadraturesLegendre(n),               std::invalid_argument) << "n=" << n;
      // Template entry point must validate before resize() mangles a negative n.
      EXPECT_THROW(NI::GetCoefsIntegrateLegendreTemplate<double>(n, 0., 1., &x, &w),
                   std::invalid_argument) << "n=" << n;
    }
  }

  // ---- Jacobi edge cases & parameter-domain validation --------------------

  TEST(NumericalIntegrationJacobi, ReducesToLegendreAndHandlesSingularAlphaBeta) {
    // Jacobi(0,0) reproduces Legendre (alpha+beta==0 removable singularity).
    NI::LegendrePolynomial leg;
    NI::JacobiPolynomial jac00(0.0, 0.0);
    auto L = NI::GaussianQuadrature(32, &leg);
    auto J = NI::GaussianQuadrature(32, &jac00);
    for (int k = 0; k < 32; ++k) {
      EXPECT_NEAR(J.first[k],  L.first[k],  1e-12) << "k=" << k;
      EXPECT_NEAR(J.second[k], L.second[k], 1e-12) << "k=" << k;
    }
    // alpha+beta == 0 (0.5,-0.5) and alpha+beta == -1 (Chebyshev) both -> mu0 = pi.
    for (auto ab : { std::make_pair(0.5, -0.5), std::make_pair(-0.5, -0.5) }) {
      NI::JacobiPolynomial jac(ab.first, ab.second);
      auto q = NI::GaussianQuadrature(64, &jac);
      double s = 0; for (double w : q.second) { EXPECT_TRUE(std::isfinite(w)); s += w; }
      EXPECT_NEAR(s, std::acos(-1.0), 1e-12) << "alpha=" << ab.first << " beta=" << ab.second;
    }
  }

  TEST(NumericalIntegrationValidation, InvalidPolynomialParametersThrow) {
    EXPECT_THROW(NI::GeneralizedLaguerrePolynomial(-1.0), std::invalid_argument);
    EXPECT_NO_THROW(NI::GeneralizedLaguerrePolynomial(-0.5));
    EXPECT_THROW(NI::JacobiPolynomial(-1.0, 0.0), std::invalid_argument);
    EXPECT_THROW(NI::JacobiPolynomial(0.0, -1.0), std::invalid_argument);
    EXPECT_NO_THROW(NI::JacobiPolynomial(-0.5, -0.5));
  }

}
