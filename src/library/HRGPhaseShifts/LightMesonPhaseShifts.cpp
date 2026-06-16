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

    // ---- pi-K S-wave phase shifts -------------------------------------------
    // Pelaez, Rodas, Phys.Rev. D93 (2016) 074025 [arXiv:1602.08404], Eqs.(11),(13),
    // Tables I (I=3/2) and II (I=1/2), CFD. Conformal variable
    //   omega(y) = (sqrt y - alpha sqrt(y0 - y)) / (sqrt y + alpha sqrt(y0 - y)),
    //   y(s) = ((s - Delta_Kpi)/(s + Delta_Kpi))^2,  Delta_Kpi = m_K^2 - m_pi^2.
    //   cot delta = -+ (sqrt s / (2q)) * 1/(s - sAdler) * {B0 + B1 omega + ...}
    // (minus = repulsive I=3/2; plus = attractive I=1/2).
    namespace {
      const double mK_piK    = KaonMass();    // 0.496 GeV (isospin-averaged)
      const double meta_piK  = EtaMass();     // 0.5478 GeV
      const double SigmaKpi  = mK_piK * mK_piK + mpi * mpi;   // m_K^2 + m_pi^2
      const double DeltaKpi  = mK_piK * mK_piK - mpi * mpi;   // m_K^2 - m_pi^2
      const double sthr_piK  = (mK_piK + mpi) * (mK_piK + mpi);
      const double sKeta_piK = (mK_piK + meta_piK) * (mK_piK + meta_piK);

      // S^{3/2} parameters (Table I, CFD)
      const double B0_K32 = 2.27, B1_K32 = 3.94, B2_K32 = 3.36;
      const double sAdler_K32 = SigmaKpi;            // ChPT LO
      const double alpha_K32 = 1.4, s0_K32 = 1.84 * 1.84;
      // S^{1/2} elastic parameters (Table II, CFD)
      const double B0_K12 = 0.411, B1_K12 = 0.162;
      const double sAdler_K12 = 0.236;               // GeV^2, ChPT LO (Eq.14)
      const double alpha_K12 = 1.15, s0_K12 = 1.1 * 1.1;
      // P^{3/2} (Eq.26, Table V) and D^{3/2} (Eq.35, Table VIII), CFD. Both use
      // alpha=1.45, s0=(1.84 GeV)^2 (same s0 as S^{3/2}, hence y0_K32).
      const double alpha_PD32 = 1.45;
      const double B0_P32 = -15.6, B1_P32 = -2.2;
      const double B0_D32 = -1.67, B1_D32 = -7.0, B2_D32 = -38.;
      // I=1/2 P-wave = K*(892) (Eq.28, Table VI elastic), CFD. alpha=1.15,
      // s0=(1.1 GeV)^2 (same as S^{1/2}, hence y0_K12). The (m_r^2 - s) factor
      // makes the phase cross 90 deg at the K*(892).
      const double B0_K892 = 0.97, B1_K892 = 0.55, B2_K892 = 0.75;
      const double mr_K892 = 0.8957;                 // GeV (Table VI)

      inline double yofs_piK(double s) {
        double r = (s - DeltaKpi) / (s + DeltaKpi);
        return r * r;
      }
      const double y0_K32 = yofs_piK(s0_K32);
      const double y0_K12 = yofs_piK(s0_K12);

      inline double qCM_piK(double M) {
        double s = M * M;
        double a1 = s - sthr_piK;
        double a2 = s - (mK_piK - mpi) * (mK_piK - mpi);
        return (a1 > 0.) ? std::sqrt(a1 * a2) / (2. * M) : 0.;
      }
      inline double omega_piK(double y, double alpha, double y0) {
        if (y >= y0) return 1.;
        double sy = std::sqrt(y), sd = std::sqrt(y0 - y);
        return (sy - alpha * sd) / (sy + alpha * sd);
      }
    } // anonymous namespace

    double PiK_delta_I32_S(double M) {
      double s = M * M;
      if (s <= sthr_piK) return 0.;
      double q = qCM_piK(M);
      if (q < 1.e-12) return 0.;
      double w = omega_piK(yofs_piK(s), alpha_K32, y0_K32);
      // minus sign: repulsive channel (delta < 0)
      double cotd = -(M / (2. * q * (s - sAdler_K32))) * (B0_K32 + B1_K32 * w + B2_K32 * w * w);
      return std::atan(1. / cotd);
    }

    double PiK_delta_I12_S(double M) {
      double s = M * M;
      if (s <= sthr_piK || s >= sKeta_piK) return 0.;   // elastic region only
      double q = qCM_piK(M);
      if (q < 1.e-12) return 0.;
      double w = omega_piK(yofs_piK(s), alpha_K12, y0_K12);
      // plus sign: attractive channel (delta > 0, kappa/K0*(700))
      double cotd = +(M / (2. * q * (s - sAdler_K12))) * (B0_K12 + B1_K12 * w);
      return std::atan(1. / cotd);
    }

    double PiK_delta_I32_P(double M) {       // Eq.26, non-resonant (small, repulsive)
      double s = M * M;
      if (s <= sthr_piK) return 0.;
      double q = qCM_piK(M);
      if (q < 1.e-12) return 0.;
      double w = omega_piK(yofs_piK(s), alpha_PD32, y0_K32);
      double cotd = std::sqrt(s) / (2. * q * q * q) * (B0_P32 + B1_P32 * w);
      return std::atan(1. / cotd);
    }

    double PiK_delta_I32_D(double M) {       // Eq.35, non-resonant (small)
      double s = M * M;
      if (s <= sthr_piK) return 0.;
      double q = qCM_piK(M);
      if (q < 1.e-12) return 0.;
      double q5 = q * q * q * q * q;
      double w = omega_piK(yofs_piK(s), alpha_PD32, y0_K32);
      double cotd = std::sqrt(s) / (2. * q5) * (B0_D32 + B1_D32 * w + B2_D32 * w * w);
      return std::atan(1. / cotd);
    }

    double PiK_delta_I12_P(double M) {       // K*(892), Eq.28, resonant
      double s = M * M;
      if (s <= sthr_piK || s >= sKeta_piK) return 0.;   // elastic region only
      double q = qCM_piK(M);
      if (q < 1.e-12) return 0.;
      double w = omega_piK(yofs_piK(s), alpha_K12, y0_K12);
      double cotd = std::sqrt(s) / (2. * q * q * q) * (mr_K892 * mr_K892 - s)
                  * (B0_K892 + B1_K892 * w + B2_K892 * w * w);
      return std::atan2(1.0, cotd);          // branch-tracked through 90 deg at K*(892)
    }

    std::vector<PhaseShiftPartialWave> PiK_I32_Waves() {
      std::vector<PhaseShiftPartialWave> w;
      w.push_back(PhaseShiftPartialWave(1, &PiK_delta_I32_S));   // S (2J+1=1)
      w.push_back(PhaseShiftPartialWave(3, &PiK_delta_I32_P));   // P (2J+1=3)
      w.push_back(PhaseShiftPartialWave(5, &PiK_delta_I32_D));   // D (2J+1=5)
      return w;
    }
    std::vector<PhaseShiftPartialWave> PiK_I12_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(1, &PiK_delta_I12_S));
    }
    std::vector<PhaseShiftPartialWave> PiK_K892_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(3, &PiK_delta_I12_P));  // P-wave
    }

    // ---- pi-pi I=0 S-wave (the sigma/f0(500)) -------------------------------
    // Garcia-Martin et al., PRD83 (2011) 074004 [1102.2183], Appendix A.1
    // (Eqs. A1-A3, Table V CFD). Low energy (s <= sM): conformal expansion of
    // cot delta. Intermediate (sM < s < (1.42 GeV)^2): matched polynomial in the
    // kaon/eta momenta (in degrees). The phase rises through 90 deg (sigma), so
    // the low-energy branch uses atan2(1, cot delta) for a continuous delta(s).
    namespace {
      const double mK_I0   = 0.496;       // M_K  [GeV]
      const double meta_I0 = 0.54751;     // M_eta [GeV]
      const double deg2rad = 3.14159265358979323846 / 180.;

      // low-energy conformal parameters (Table V, CFD)
      const double S0_B0 = 7.14, S0_B1 = -25.3, S0_B2 = -33.2, S0_B3 = -26.2;
      const double S0_z0 = mpi;                       // z0 = M_pi
      const double S0_s0conf = 4. * mK_I0 * mK_I0;    // conformal s0 = 4 M_K^2
      const double S0_sM = 0.85 * 0.85;               // matching point [GeV^2]
      const double S0_4MK2   = 4. * mK_I0 * mK_I0;
      const double S0_4Meta2 = 4. * meta_I0 * meta_I0;
      // intermediate parameters (Table V, CFD), converted from degrees
      const double S0_d0 = 226.5 * deg2rad;
      const double S0_c  = -81.0 * deg2rad;
      const double S0_B  =  93.3 * deg2rad;
      const double S0_C  =  48.7 * deg2rad;
      const double S0_D  = -88.3 * deg2rad;

      inline double w_I0(double s) {                  // conformal variable, Eq.A2
        double a = std::sqrt(s), b = std::sqrt(S0_s0conf - s);
        return (a - b) / (a + b);
      }
      inline double cotd_I0_lowE(double s) {          // Eq.A1
        double k = std::sqrt(s / 4. - mpi * mpi);     // CM momentum
        double w = w_I0(s);
        double poly = S0_z0 * S0_z0 / (mpi * std::sqrt(s))
                    + S0_B0 + S0_B1 * w + S0_B2 * w * w + S0_B3 * w * w * w;
        return std::sqrt(s) * mpi * mpi / (2. * k * (s - 0.5 * S0_z0 * S0_z0)) * poly;
      }
      inline double delta_I0_lowE(double s) {         // radians, branch-tracked in (0, pi)
        return std::atan2(1.0, cotd_I0_lowE(s));
      }
      // matching-point values (computed once from the low-energy parametrization)
      const double S0_kM = std::sqrt(mK_I0 * mK_I0 - S0_sM / 4.);   // |k2(sM)|
      const double S0_deltaM = delta_I0_lowE(S0_sM);                // radians
      const double S0_deltaMp = (delta_I0_lowE(S0_sM + 1e-5)
                               - delta_I0_lowE(S0_sM - 1e-5)) / (2e-5);  // d delta/ds
    } // anonymous namespace

    double PiPi_delta_I0_S(double M) {
      double s = M * M;
      if (s <= 4. * mpi * mpi) return 0.;
      if (s <= S0_sM) return delta_I0_lowE(s);        // low energy (Eq.A1), branch-tracked
      if (s < S0_4MK2) {                              // Eq.A3, sM < s < 4 M_K^2
        double k2 = std::sqrt(mK_I0 * mK_I0 - s / 4.);   // |k2| (real below K-Kbar)
        double r  = k2 / S0_kM;
        return S0_d0 * (1. - r) * (1. - r)
             + S0_deltaM * r * (2. - r)
             + k2 * (S0_kM - k2) * (8. * S0_deltaMp + S0_c * (S0_kM - k2) / (mK_I0 * mK_I0 * mK_I0));
      }
      // Eq.A3, 4 M_K^2 < s < (1.42 GeV)^2 (above the elastic limit; only reached
      // if the channel Mmax is raised past 2 M_K)
      double u = (s / 4. - mK_I0 * mK_I0) / (mK_I0 * mK_I0);   // k2^2 / M_K^2
      double d = S0_d0 + S0_B * u + S0_C * u * u;
      if (s > S0_4Meta2)
        d += S0_D * (s / 4. - meta_I0 * meta_I0) / (meta_I0 * meta_I0);  // eta-momentum term
      return d;
    }

    double PiPi_delta_I0_f0980_S(double M) {
      double s = M * M;
      if (s <= S0_4MK2) return 0.;          // below K-Kbar: handled by the sigma wave
      // offset so the wave starts at 0 at 2 M_K (S0_d0 = delta_0^0 there); only
      // d(delta)/dM matters for the spectral weight, the offset is cosmetic.
      return PiPi_delta_I0_S(M) - S0_d0;
    }

    std::vector<PhaseShiftPartialWave> PiPi_I0_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(1, &PiPi_delta_I0_S));
    }

    std::vector<PhaseShiftPartialWave> PiPi_I0_f0980_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(1, &PiPi_delta_I0_f0980_S));
    }

    // ---- pi-pi I=1 P-wave (the rho(770)) -----------------------------------
    // Garcia-Martin et al., PRD83 (2011) 074004 [1102.2183], Eq.A7, Table VII CFD:
    //   cot delta_1^1(s) = (sqrt s / 2k^3)(M_rho^2 - s)
    //                      [ 2 M_pi^3 / (M_rho^2 sqrt s) + B0 + B1 w(s) ],
    //   w(s) = (sqrt s - sqrt(s0 - s))/(sqrt s + sqrt(s0 - s)), s0 = (1.05 GeV)^2.
    // (M_rho^2 - s) makes cot delta cross 0 at the rho (delta = 90 deg), so the
    // phase is branch-tracked with atan2. Elastic up to the K-Kbar threshold.
    namespace {
      const double Mrho_P = 0.7736;        // GeV (Garcia-Martin value, 773.6 MeV)
      const double P_B0 = 1.043, P_B1 = 0.19;   // Table VII, CFD
      const double P_s0 = 1.05 * 1.05;     // GeV^2 (conformal pivot)
    }

    double PiPi_delta_I1_P(double M) {
      double s = M * M;
      if (s <= 4. * mpi * mpi) return 0.;
      double k = std::sqrt(s / 4. - mpi * mpi);   // CM momentum
      if (k < 1.e-12) return 0.;
      double w = conformal_w(s, P_s0);
      double cotd = std::sqrt(s) / (2. * k * k * k) * (Mrho_P * Mrho_P - s)
                  * (2. * mpi * mpi * mpi / (Mrho_P * Mrho_P * std::sqrt(s)) + P_B0 + P_B1 * w);
      return std::atan2(1.0, cotd);               // branch-tracked through the rho (90 deg)
    }

    std::vector<PhaseShiftPartialWave> PiPi_I1_Waves() {
      // P-wave: J=1 -> 2J+1 = 3.
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(3, &PiPi_delta_I1_P));
    }

    // ---- pi-pi I=1 F-wave ---------------------------------------------------
    // Garcia-Martin et al., PRD83 (2011) 074004 [1102.2183], Eq.A14, Table X CFD:
    //   cot delta_3^1(s) = (sqrt s / 2k^7) M_pi^6 [ 2 lambda M_pi / sqrt s + B0 + B1 w(s) ],
    //   w(s) = (sqrt s - sqrt(s0 - s))/(sqrt s + sqrt(s0 - s)), s0 = (1.45 GeV)^2.
    // Non-resonant below 1.42 GeV (rho3(1690) is higher) and small (k^7), so the
    // single-branch atan is fine; the paper neglects the inelasticity to 1.42 GeV.
    namespace {
      const double F_B0 = 1.09e5, F_B1 = 1.41e5, F_lambda = 0.051e5;  // Table X, CFD
      const double F_s0 = 1.45 * 1.45;     // GeV^2 (conformal pivot)
    }

    double PiPi_delta_I1_F(double M) {
      double s = M * M;
      if (s <= 4. * mpi * mpi) return 0.;
      double k = std::sqrt(s / 4. - mpi * mpi);   // CM momentum
      if (k < 1.e-12) return 0.;
      double k7 = k * k * k * k * k * k * k;
      double mpi6 = mpi * mpi * mpi * mpi * mpi * mpi;
      double w = conformal_w(s, F_s0);
      double cotd = std::sqrt(s) / (2. * k7) * mpi6
                  * (2. * F_lambda * mpi / std::sqrt(s) + F_B0 + F_B1 * w);
      return std::atan(1. / cotd);                // small, attractive (non-resonant)
    }

    std::vector<PhaseShiftPartialWave> PiPi_I1_F_Waves() {
      // F-wave: J=3 -> 2J+1 = 7.
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(7, &PiPi_delta_I1_F));
    }

  } // namespace PhaseShifts

} // namespace thermalfist
