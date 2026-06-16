/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGPhaseShifts/MesonBaryonPhaseShifts.h"
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"   // PionMass()

#include <cmath>

namespace thermalfist {

  namespace PhaseShifts {

    namespace {
      const double mpi = PionMass();
      const double mN  = NucleonMass();
      const double sthr_piN = (mN + mpi) * (mN + mpi);          // pi-N threshold s_+
      const double sdif_piN = (mN - mpi) * (mN - mpi);

      // Delta(1232) = P33 conformal parameters (Roy-Steiner, arXiv:1510.06039,
      // Table D.2, I_s = 3/2, 1+ wave). All in appropriate powers of GeV.
      const double Atil_P33 = 7.781e1;     // 1/Atil is the CONSTANT term of the bracket
      const double Btil_P33 = -3.986e-2;
      const double Ctil_P33 = -0.3098;
      const double s1_P33   = 0.4509;      // position parameter [GeV^2]
      const double sbar_P33 = 1.540 * 1.540; // conformal branch point sbar = (1.540 GeV)^2

      inline double qCM_piN(double M) {
        double s = M * M;
        double a1 = s - sthr_piN;
        double a2 = s - sdif_piN;
        return (a1 > 0.) ? std::sqrt(a1 * a2) / (2. * M) : 0.;
      }
      // Conformal variable w(s) = (sqrt(s) - sqrt(sbar - s)) / (sqrt(s) + sqrt(sbar - s)).
      inline double wofs_piN(double s) {
        double sd = sbar_P33 - s;
        if (sd <= 0.) return 1.;                 // s >= sbar (outside the elastic range)
        double rs = std::sqrt(s), rd = std::sqrt(sd);
        return (rs - rd) / (rs + rd);
      }

      // Schenk S-wave (Eq. D.1): tan d = |q| (A + B q^2 + C q^4 + D q^6 + E q^8)
      //                                      (s+ - s0)/(s - s0).   Non-resonant.
      inline double schenkS(double M, double A, double B, double C, double D, double E, double s0) {
        double s = M * M;
        if (s <= sthr_piN) return 0.;
        double q = qCM_piN(M);
        if (q < 1.e-12) return 0.;
        double q2 = q * q;
        double poly = A + q2 * (B + q2 * (C + q2 * (D + q2 * E)));
        return std::atan(q * poly * (sthr_piN - s0) / (s - s0));
      }
      // Schenk P-wave (Eq. D.2): tan d = |q|^3 (A + B q^2 + C q^4 + D q^6)
      //                                        (s+ - s1)/(s - s1).   Non-resonant.
      inline double schenkP(double M, double A, double B, double C, double D, double s1) {
        double s = M * M;
        if (s <= sthr_piN) return 0.;
        double q = qCM_piN(M);
        if (q < 1.e-12) return 0.;
        double q2 = q * q;
        double poly = A + q2 * (B + q2 * (C + q2 * D));
        return std::atan(q * q * q * poly * (sthr_piN - s1) / (s - s1));
      }
    } // anonymous namespace

    double PiN_delta_P33(double M) {            // Eq. D.3, resonant (the Delta)
      double s = M * M;
      if (s <= sthr_piN) return 0.;
      double q = qCM_piN(M);
      if (q < 1.e-12) return 0.;
      double x = wofs_piN(s) - wofs_piN(sthr_piN);
      // cot d = (1/q^3) (s - s1)/(s+ - s1) [ 1/Atil + Btil x + Ctil x^2 ]
      double cotd = (1. / (q * q * q)) * (s - s1_P33) / (sthr_piN - s1_P33)
                  * (1. / Atil_P33 + Btil_P33 * x + Ctil_P33 * x * x);
      return std::atan2(1.0, cotd);             // branch-tracked through 90 deg at the Delta
    }

    // S-waves (Eq. D.1). Coefficients from Table D.1 (in powers of GeV); the
    // leading A is the scattering length. S31 repulsive, S11 attractive.
    double PiN_delta_S31(double M) {
      return schenkS(M, -0.6183, -1.831e1, 3.090e2, -2.846e3, 9.529e3, -1.809e3);
    }
    double PiN_delta_S11(double M) {
      return schenkS(M, 1.217, -1.879e1, 1.958e2, -1.235e3, 3.350e3, 2.494);
    }
    // P-waves (Eq. D.2). Coefficients from Table D.2 (in powers of GeV).
    double PiN_delta_P31(double M) {
      return schenkP(M, -1.477e1, 1.467e2, -1.633e3, 6.508e3, 0.4081);
    }
    double PiN_delta_P11(double M) {
      return schenkP(M, -2.569e1, 8.062e2, -4.214e3, 3.986e4, 0.9340);
    }
    double PiN_delta_P13(double M) {
      return schenkP(M, -1.085e1, -1.145e1, 3.651e2, -1.052e3, 0.9639);
    }

    std::vector<PhaseShiftPartialWave> PiN_Delta_Waves() {
      // J = 3/2 -> 2J+1 = 4 (the spectral-weight degeneracy of the Delta wave).
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(4, &PiN_delta_P33));
    }
    // S/P J=1/2 waves -> 2J+1 = 2; the P13 J=3/2 wave -> 2J+1 = 4.
    std::vector<PhaseShiftPartialWave> PiN_S31_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(2, &PiN_delta_S31));
    }
    std::vector<PhaseShiftPartialWave> PiN_S11_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(2, &PiN_delta_S11));
    }
    std::vector<PhaseShiftPartialWave> PiN_P31_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(2, &PiN_delta_P31));
    }
    std::vector<PhaseShiftPartialWave> PiN_P11_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(2, &PiN_delta_P11));
    }
    std::vector<PhaseShiftPartialWave> PiN_P13_Waves() {
      return std::vector<PhaseShiftPartialWave>(1, PhaseShiftPartialWave(4, &PiN_delta_P13));
    }

  } // namespace PhaseShifts

} // namespace thermalfist
