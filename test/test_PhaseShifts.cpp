/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <vector>
#include <cmath>
#include "HRGPhaseShifts/PhaseShiftDensity.h"
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"
#include "HRGBase/IdealGasFunctions.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  // Build the pi-pi I=2 channel S-matrix d.o.f. (repulsive, S+D waves).
  // Uses the catalog waves, which carry the exact analytic d(delta)/dM.
  PhaseShiftDensity makePiPiI2(int nodes = 64) {
    return PhaseShiftDensity(PhaseShifts::PiPi_I2_Waves(),
                             PhaseShifts::PionMass(), PhaseShifts::PionMass(),
                             PhaseShifts::PiPiI2_Mmax(), /*Bose*/ -1, nodes);
  }

  // Same channel but WITHOUT analytic derivatives -> engine uses the
  // central-finite-difference fallback for d(delta)/dq.
  PhaseShiftDensity makePiPiI2FiniteDiff(int nodes = 64) {
    std::vector<PhaseShiftPartialWave> waves;
    waves.push_back(PhaseShiftPartialWave(1, &PhaseShifts::PiPi_delta_I2_S)); // no ddelta_dM
    waves.push_back(PhaseShiftPartialWave(5, &PhaseShifts::PiPi_delta_I2_D)); // no ddelta_dM
    return PhaseShiftDensity(waves,
                             PhaseShifts::PionMass(), PhaseShifts::PionMass(),
                             PhaseShifts::PiPiI2_Mmax(), /*Bose*/ -1, nodes);
  }

  // chi_2^Q contribution of the pi-pi I=2 multiplet at mu_Q = 0:
  // Iz = +-2 (Q=+-2), +-1 (Q=+-1), 0 (Q=0) -> sum Q^2 = 2*(4+1) = 10 times chi2(channel).
  double chi2QPiPiI2(PhaseShiftDensity& psd, double T) {
    return 10.0 * psd.Quantity(IdealGasFunctions::chi2, T, 0.0);
  }

  TEST(PhaseShifts, PiPiI2IsRepulsive) {
    PhaseShiftDensity psd = makePiPiI2();
    // The I=2 phase shift is negative and decreasing -> spectral weight < 0,
    // and it must be finite everywhere on (0, qmax) including near threshold.
    for (double q : {0.02, 0.05, 0.1, 0.2, 0.3}) {
      double w = psd.SpectralWeight(q);
      EXPECT_TRUE(std::isfinite(w)) << "q=" << q;
      EXPECT_LT(w, 0.0) << "q=" << q;   // repulsive
    }
  }

  TEST(PhaseShifts, PhaseShiftSignAndThreshold) {
    // delta_I2 < 0 above threshold; exactly 0 at/below threshold.
    double Mthr = 2.0 * PhaseShifts::PionMass();
    EXPECT_DOUBLE_EQ(PhaseShifts::PiPi_delta_I2_S(Mthr), 0.0);
    EXPECT_LT(PhaseShifts::PiPi_delta_I2_S(0.5), 0.0);
    EXPECT_LT(PhaseShifts::PiPi_delta_I2_S(0.8), 0.0);
    EXPECT_LT(PhaseShifts::PiPi_delta_I2_D(0.8), 0.0);
  }

  TEST(PhaseShifts, Chi2QContributionIsNegative) {
    // A repulsive channel suppresses electric-charge fluctuations.
    PhaseShiftDensity psd = makePiPiI2();
    double c2 = chi2QPiPiI2(psd, 0.150);
    EXPECT_TRUE(std::isfinite(c2));
    EXPECT_LT(c2, 0.0);
  }

  TEST(PhaseShifts, AnalyticMatchesFiniteDiff) {
    // The analytic d(delta)/dM overload and the finite-difference fallback must
    // agree to ~FD accuracy across the whole q range (away from the q=0 kink).
    PhaseShiftDensity psA = makePiPiI2();           // analytic
    PhaseShiftDensity psF = makePiPiI2FiniteDiff(); // finite diff
    for (double q : {0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7}) {
      double wA = psA.SpectralWeight(q);
      double wF = psF.SpectralWeight(q);
      EXPECT_NEAR(wA, wF, 1e-6 * std::fabs(wA) + 1e-9) << "q=" << q;
    }
    // and the integrated thermodynamics must match to quadrature precision.
    double cA = chi2QPiPiI2(psA, 0.150);
    double cF = chi2QPiPiI2(psF, 0.150);
    EXPECT_NEAR(cA, cF, 1e-6 * std::fabs(cA));
  }

  TEST(PhaseShifts, SpectralWeightReferenceValues) {
    // Analytic spectral weight (1/pi) sum_l (2J_l+1) ddelta_l/dq vs an
    // independent mpmath high-precision reference (S+D, charged-pion mass).
    PhaseShiftDensity psd = makePiPiI2();
    EXPECT_NEAR(psd.SpectralWeight(0.05), -1.430522040942e-01, 1e-9);
    EXPECT_NEAR(psd.SpectralWeight(0.20), -3.366924650159e-01, 1e-9);
    EXPECT_NEAR(psd.SpectralWeight(0.50), -4.885358082049e-01, 1e-9);
  }

  TEST(PhaseShifts, QuadratureConvergence) {
    // The q-integral must be converged: 32 -> 64 -> 128 nodes should agree well.
    double T = 0.150;
    PhaseShiftDensity p32 = makePiPiI2(32);
    PhaseShiftDensity p64 = makePiPiI2(64);
    PhaseShiftDensity p128 = makePiPiI2(128);
    double c32  = chi2QPiPiI2(p32,  T);
    double c64  = chi2QPiPiI2(p64,  T);
    double c128 = chi2QPiPiI2(p128, T);
    ASSERT_LT(std::fabs(c128), 1.0);          // sanity: O(<1) contribution
    EXPECT_NEAR(c64,  c128, 1e-4 * std::fabs(c128));
    EXPECT_NEAR(c32,  c128, 1e-2 * std::fabs(c128));
  }

}
