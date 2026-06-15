/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"

#include <cmath>

namespace thermalfist {

  namespace PhaseShifts {

    namespace {
      const double mpi = PionMass();

      // ---- pi-pi I=2 S-wave (J=0), Garcia-Martin et al. (2011), Eqs. A5/A6, Table VI ----
      const double S_B0  = -79.4;
      const double S_B1  = -63.0;
      const double S_z2  = 0.1435;        // GeV (Adler zero at s = 2 z2^2)
      const double S_sl  = 1.05 * 1.05;   // GeV^2 (low-energy conformal point)
      const double S_sh  = 1.42 * 1.42;   // GeV^2 (intermediate conformal point)
      const double S_sM  = 0.85 * 0.85;   // GeV^2 (matching point A5 <-> A6)
      const double S_Bh2 = 32.;

      inline double conformal_w(double s, double s0) {
        if (s >= s0) return 1.;
        double sqrts = std::sqrt(s), sqrtd = std::sqrt(s0 - s);
        return (sqrts - sqrtd) / (sqrts + sqrtd);
      }

      // d w/ds of the conformal variable above (0 for s >= s0).
      inline double conformal_dw(double s, double s0) {
        if (s >= s0) return 0.;
        double a = std::sqrt(s), b = std::sqrt(s0 - s);
        return s0 / (a * b * (a + b) * (a + b));
      }

      // Precomputed matching coefficients for Eq. A6.
      const double S_wl_sM = conformal_w(S_sM, S_sl);
      const double S_wh_sM = conformal_w(S_sM, S_sh);
      const double S_Bh0   = S_B0 + S_B1 * S_wl_sM;
      const double S_Bh1   = [] {
        double sqrtsM = std::sqrt(S_sM);
        double sqrt_sh_sM = std::sqrt(S_sh - S_sM), sqrt_sl_sM = std::sqrt(S_sl - S_sM);
        return S_B1 * (S_sl * sqrt_sh_sM) / (S_sh * sqrt_sl_sM)
               * std::pow((sqrtsM + sqrt_sh_sM) / (sqrtsM + sqrt_sl_sM), 2);
      }();

      // ---- pi-pi I=2 D-wave (J=2), Eq. A12, Table IX ----
      const double D_B0  = 4.1e3;
      const double D_B1  = 8.6e3;
      const double D_B2  = 25.5e3;
      const double D_Delta = 0.233;          // GeV
      const double D_s0  = 1.45 * 1.45;      // GeV^2
      const double D_sAdler = 4. * (mpi * mpi + D_Delta * D_Delta);
    } // anonymous namespace

    double PiPi_delta_I2_S(double M) {
      double s = M * M;
      double arg = s - 4. * mpi * mpi;
      if (arg <= 0.) return 0.;
      double k = 0.5 * std::sqrt(arg);
      if (k < 1.e-12) return 0.;
      double prefactor = (M / (2. * k)) * (mpi * mpi / (s - 2. * S_z2 * S_z2));
      double cotd;
      if (s <= S_sM) {
        cotd = prefactor * (S_B0 + S_B1 * conformal_w(s, S_sl));
      } else {
        double dw = conformal_w(s, S_sh) - S_wh_sM;
        cotd = prefactor * (S_Bh0 + S_Bh1 * dw + S_Bh2 * dw * dw);
      }
      // |delta| << pi/2 throughout the I=2 S-wave, so atan(1/cotd) is single-branch.
      return std::atan(1. / cotd);
    }

    double PiPi_delta_I2_D(double M) {
      double s = M * M;
      double arg = s - 4. * mpi * mpi;
      if (arg <= 0.) return 0.;
      double k = 0.5 * std::sqrt(arg);
      if (k < 1.e-12) return 0.;
      double k2 = k * k, k5 = k2 * k2 * k;
      double prefactor = (M / (2. * k5)) * (mpi * mpi * mpi * mpi * s) / (D_sAdler - s);
      double wval = conformal_w(s, D_s0);
      double cotd = prefactor * (D_B0 + D_B1 * wval + D_B2 * wval * wval);
      return std::atan(1. / cotd);
    }

    // Analytic d(delta)/dM via the logarithmic derivative of cot(delta):
    //   delta = atan(1/u),  ddelta/dM = -2M (du/ds) / (1 + u^2),
    //   du/ds = u * d ln u / ds   (u = product of factors -> sum of log-derivatives).
    // Verified against an mpmath high-precision reference to ~1e-16 over the whole
    // validity range, including the S-wave A5/A6 match and the D-wave Adler point.
    double PiPi_ddelta_I2_S_dM(double M) {
      double s = M * M;
      double arg = s - 4. * mpi * mpi;
      if (arg <= 0.) return 0.;
      double k = 0.5 * std::sqrt(arg);
      if (k < 1.e-12) return 0.;
      double F1 = M / (2. * k);
      double F2 = mpi * mpi / (s - 2. * S_z2 * S_z2);
      double P, dPds;
      if (s <= S_sM) {
        double w = conformal_w(s, S_sl), dw = conformal_dw(s, S_sl);
        P = S_B0 + S_B1 * w;  dPds = S_B1 * dw;
      } else {
        double dwv = conformal_w(s, S_sh) - S_wh_sM, dwh = conformal_dw(s, S_sh);
        P = S_Bh0 + S_Bh1 * dwv + S_Bh2 * dwv * dwv;  dPds = (S_Bh1 + 2. * S_Bh2 * dwv) * dwh;
      }
      double u = F1 * F2 * P;
      double duds = u * (0.5 * (1. / s - 1. / arg) - 1. / (s - 2. * S_z2 * S_z2) + dPds / P);
      return -2. * M * duds / (1. + u * u);
    }

    double PiPi_ddelta_I2_D_dM(double M) {
      double s = M * M;
      double arg = s - 4. * mpi * mpi;
      if (arg <= 0.) return 0.;
      double k = 0.5 * std::sqrt(arg);
      if (k < 1.e-12) return 0.;
      double k2 = k * k, k5 = k2 * k2 * k;
      double G1 = M / (2. * k5);
      double G2 = mpi * mpi * mpi * mpi * s / (D_sAdler - s);
      double w = conformal_w(s, D_s0), dw = conformal_dw(s, D_s0);
      double P = D_B0 + D_B1 * w + D_B2 * w * w;
      double dPds = (D_B1 + 2. * D_B2 * w) * dw;
      double u = G1 * G2 * P;
      double duds = u * ((0.5 / s - 2.5 / arg) + (1. / s + 1. / (D_sAdler - s)) + dPds / P);
      return -2. * M * duds / (1. + u * u);
    }

    std::vector<PhaseShiftPartialWave> PiPi_I2_Waves() {
      std::vector<PhaseShiftPartialWave> waves;
      waves.push_back(PhaseShiftPartialWave(1, &PiPi_delta_I2_S, &PiPi_ddelta_I2_S_dM)); // J=0 -> 2J+1 = 1
      waves.push_back(PhaseShiftPartialWave(5, &PiPi_delta_I2_D, &PiPi_ddelta_I2_D_dM)); // J=2 -> 2J+1 = 5
      return waves;
    }

  } // namespace PhaseShifts

} // namespace thermalfist
