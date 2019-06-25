/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/xMath.h"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cfloat>

using namespace std;

namespace thermalfist {

  // Implementation of the special functions below adapted from
  // CERN-ROOT package: https://root.cern.ch/

  //______________________________________________________________________________
  double xMath::BesselI0(double x)
  {
    // Compute the modified Bessel function I_0(x) for any real x.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht

    // Parameters of the polynomial approximation
    const double p1 = 1.0, p2 = 3.5156229, p3 = 3.0899424,
      p4 = 1.2067492, p5 = 0.2659732, p6 = 3.60768e-2, p7 = 4.5813e-3;

    const double q1 = 0.39894228, q2 = 1.328592e-2, q3 = 2.25319e-3,
      q4 = -1.57565e-3, q5 = 9.16281e-3, q6 = -2.057706e-2,
      q7 = 2.635537e-2, q8 = -1.647633e-2, q9 = 3.92377e-3;

    const double k1 = 3.75;
    double ax = fabs(x);

    double y = 0, result = 0;

    if (ax < k1) {
      double xx = x / k1;
      y = xx * xx;
      result = p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7)))));
    }
    else {
      y = k1 / ax;
      result = (exp(ax) / sqrt(ax))*(q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))))))));
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselK0(double x)
  {
    // Compute the modified Bessel function K_0(x) for positive real x.
    //
    //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
    //     Applied Mathematics Series vol. 55 (1964), Washington.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht

    // Parameters of the polynomial approximation
    const double p1 = -0.57721566, p2 = 0.42278420, p3 = 0.23069756,
      p4 = 3.488590e-2, p5 = 2.62698e-3, p6 = 1.0750e-4, p7 = 7.4e-6;

    const double q1 = 1.25331414, q2 = -7.832358e-2, q3 = 2.189568e-2,
      q4 = -1.062446e-2, q5 = 5.87872e-3, q6 = -2.51540e-3, q7 = 5.3208e-4;

    if (x <= 0) {
      //Error("xMath::BesselK0", "*K0* Invalid argument x = %g\n",x);
      return 0;
    }

    double y = 0, result = 0;

    if (x <= 2) {
      y = x * x / 4;
      result = (-log(x / 2.)*xMath::BesselI0(x)) + (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7))))));
    }
    else {
      y = 2 / x;
      result = (exp(-x) / sqrt(x))*(q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * q7))))));
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselI1(double x)
  {
    // Compute the modified Bessel function I_1(x) for any real x.
    //
    //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
    //     Applied Mathematics Series vol. 55 (1964), Washington.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht

    // Parameters of the polynomial approximation
    const double p1 = 0.5, p2 = 0.87890594, p3 = 0.51498869,
      p4 = 0.15084934, p5 = 2.658733e-2, p6 = 3.01532e-3, p7 = 3.2411e-4;

    const double q1 = 0.39894228, q2 = -3.988024e-2, q3 = -3.62018e-3,
      q4 = 1.63801e-3, q5 = -1.031555e-2, q6 = 2.282967e-2,
      q7 = -2.895312e-2, q8 = 1.787654e-2, q9 = -4.20059e-3;

    const double k1 = 3.75;
    double ax = fabs(x);

    double y = 0, result = 0;

    if (ax < k1) {
      double xx = x / k1;
      y = xx * xx;
      result = x * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7))))));
    }
    else {
      y = k1 / ax;
      result = (exp(ax) / sqrt(ax))*(q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))))))));
      if (x < 0) result = -result;
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselK1(double x)
  {
    // Compute the modified Bessel function K_1(x) for positive real x.
    //
    //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
    //     Applied Mathematics Series vol. 55 (1964), Washington.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht

    // Parameters of the polynomial approximation
    const double p1 = 1., p2 = 0.15443144, p3 = -0.67278579,
      p4 = -0.18156897, p5 = -1.919402e-2, p6 = -1.10404e-3, p7 = -4.686e-5;

    const double q1 = 1.25331414, q2 = 0.23498619, q3 = -3.655620e-2,
      q4 = 1.504268e-2, q5 = -7.80353e-3, q6 = 3.25614e-3, q7 = -6.8245e-4;

    if (x <= 0) {
      //Error("xMath::BesselK1", "*K1* Invalid argument x = %g\n",x);
      return 0;
    }

    double y = 0, result = 0;

    if (x <= 2) {
      y = x * x / 4;
      result = (log(x / 2.)*xMath::BesselI1(x)) + (1. / x)*(p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7))))));
    }
    else {
      y = 2 / x;
      result = (exp(-x) / sqrt(x))*(q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * q7))))));
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselK(int n, double x)
  {
    if (n < 0) n = -n;
    // Compute the Integer Order Modified Bessel function K_n(x)
     // for n=0,1,2,... and positive real x.
     //
     //--- NvE 12-mar-2000 UU-SAP Utrecht

    if (x <= 0 || n < 0) {
      //Error("xMath::BesselK", "*K* Invalid argument(s) (n,x) = (%d, %g)\n",n,x);
      return 0;
    }

    if (n == 0) return xMath::BesselK0(x);
    if (n == 1) return xMath::BesselK1(x);

    // Perform upward recurrence for all x
    double tox = 2 / x;
    double bkm = xMath::BesselK0(x);
    double bk = xMath::BesselK1(x);
    double bkp = 0;
    for (int j = 1; j < n; j++) {
      bkp = bkm + double(j)*tox*bk;
      bkm = bk;
      bk = bkp;
    }
    return bk;
  }

  //______________________________________________________________________________
  double xMath::BesselI(int n, double x)
  {
    // Compute the Integer Order Modified Bessel function I_n(x)
    // for n=0,1,2,... and any real x.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht
    if (n < 0) n = -n;

    int iacc = 40; // Increase to enhance accuracy
    const double kBigPositive = 1.e10;
    const double kBigNegative = 1.e-10;

    if (n < 0) {
      //Error("xMath::BesselI", "*I* Invalid argument (n,x) = (%d, %g)\n",n,x);
      return 0;
    }

    if (n == 0) return xMath::BesselI0(x);
    if (n == 1) return xMath::BesselI1(x);

    if (x == 0) return 0;
    if (fabs(x) > kBigPositive) return 0;

    double tox = 2 / fabs(x);
    double bip = 0, bim = 0;
    double bi = 1;
    double result = 0;
    int m = 2 * ((n + int(sqrt(float(iacc*n)))));
    for (int j = m; j >= 1; j--) {
      bim = bip + double(j)*tox*bi;
      bip = bi;
      bi = bim;
      // Renormalise to prevent overflows
      if (fabs(bi) > kBigPositive) {
        result *= kBigNegative;
        bi *= kBigNegative;
        bip *= kBigNegative;
      }
      if (j == n) result = bip;
    }

    result *= xMath::BesselI0(x) / bi; // Normalise with BesselI0(x)
    if ((x < 0) && (n % 2 == 1)) result = -result;

    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselJ0(double x)
  {
    // Returns the Bessel function J0(x) for any real x.

    double ax, z;
    double xx, y, result, result1, result2;
    const double p1 = 57568490574.0, p2 = -13362590354.0, p3 = 651619640.7;
    const double p4 = -11214424.18, p5 = 77392.33017, p6 = -184.9052456;
    const double p7 = 57568490411.0, p8 = 1029532985.0, p9 = 9494680.718;
    const double p10 = 59272.64853, p11 = 267.8532712;

    const double q1 = 0.785398164;
    const double q2 = -0.1098628627e-2, q3 = 0.2734510407e-4;
    const double q4 = -0.2073370639e-5, q5 = 0.2093887211e-6;
    const double q6 = -0.1562499995e-1, q7 = 0.1430488765e-3;
    const double q8 = -0.6911147651e-5, q9 = 0.7621095161e-6;
    const double q10 = 0.934935152e-7, q11 = 0.636619772;

    if ((ax = fabs(x)) < 8) {
      y = x * x;
      result1 = p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * p6))));
      result2 = p7 + y * (p8 + y * (p9 + y * (p10 + y * (p11 + y))));
      result = result1 / result2;
    }
    else {
      z = 8 / ax;
      y = z * z;
      xx = ax - q1;
      result1 = 1 + y * (q2 + y * (q3 + y * (q4 + y * q5)));
      result2 = q6 + y * (q7 + y * (q8 + y * (q9 - y * q10)));
      result = sqrt(q11 / ax)*(cos(xx)*result1 - z * sin(xx)*result2);
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselJ1(double x)
  {
    // Returns the Bessel function J1(x) for any real x.

    double ax, z;
    double xx, y, result, result1, result2;
    const double p1 = 72362614232.0, p2 = -7895059235.0, p3 = 242396853.1;
    const double p4 = -2972611.439, p5 = 15704.48260, p6 = -30.16036606;
    const double p7 = 144725228442.0, p8 = 2300535178.0, p9 = 18583304.74;
    const double p10 = 99447.43394, p11 = 376.9991397;

    const double q1 = 2.356194491;
    const double q2 = 0.183105e-2, q3 = -0.3516396496e-4;
    const double q4 = 0.2457520174e-5, q5 = -0.240337019e-6;
    const double q6 = 0.04687499995, q7 = -0.2002690873e-3;
    const double q8 = 0.8449199096e-5, q9 = -0.88228987e-6;
    const double q10 = 0.105787412e-6, q11 = 0.636619772;

    if ((ax = fabs(x)) < 8) {
      y = x * x;
      result1 = x * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * p6)))));
      result2 = p7 + y * (p8 + y * (p9 + y * (p10 + y * (p11 + y))));
      result = result1 / result2;
    }
    else {
      z = 8 / ax;
      y = z * z;
      xx = ax - q1;
      result1 = 1 + y * (q2 + y * (q3 + y * (q4 + y * q5)));
      result2 = q6 + y * (q7 + y * (q8 + y * (q9 + y * q10)));
      result = sqrt(q11 / ax)*(cos(xx)*result1 - z * sin(xx)*result2);
      if (x < 0) result = -result;
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselY0(double x)
  {
    // Returns the Bessel function Y0(x) for positive x.

    double z, xx, y, result, result1, result2;
    const double p1 = -2957821389., p2 = 7062834065.0, p3 = -512359803.6;
    const double p4 = 10879881.29, p5 = -86327.92757, p6 = 228.4622733;
    const double p7 = 40076544269., p8 = 745249964.8, p9 = 7189466.438;
    const double p10 = 47447.26470, p11 = 226.1030244, p12 = 0.636619772;

    const double q1 = 0.785398164;
    const double q2 = -0.1098628627e-2, q3 = 0.2734510407e-4;
    const double q4 = -0.2073370639e-5, q5 = 0.2093887211e-6;
    const double q6 = -0.1562499995e-1, q7 = 0.1430488765e-3;
    const double q8 = -0.6911147651e-5, q9 = 0.7621095161e-6;
    const double q10 = -0.934945152e-7, q11 = 0.636619772;

    if (x < 8) {
      y = x * x;
      result1 = p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * p6))));
      result2 = p7 + y * (p8 + y * (p9 + y * (p10 + y * (p11 + y))));
      result = (result1 / result2) + p12 * xMath::BesselJ0(x)*log(x);
    }
    else {
      z = 8 / x;
      y = z * z;
      xx = x - q1;
      result1 = 1 + y * (q2 + y * (q3 + y * (q4 + y * q5)));
      result2 = q6 + y * (q7 + y * (q8 + y * (q9 + y * q10)));
      result = sqrt(q11 / x)*(sin(xx)*result1 + z * cos(xx)*result2);
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselY1(double x)
  {
    // Returns the Bessel function Y1(x) for positive x.

    double z, xx, y, result, result1, result2;
    const double p1 = -0.4900604943e13, p2 = 0.1275274390e13;
    const double p3 = -0.5153438139e11, p4 = 0.7349264551e9;
    const double p5 = -0.4237922726e7, p6 = 0.8511937935e4;
    const double p7 = 0.2499580570e14, p8 = 0.4244419664e12;
    const double p9 = 0.3733650367e10, p10 = 0.2245904002e8;
    const double p11 = 0.1020426050e6, p12 = 0.3549632885e3;
    const double p13 = 0.636619772;
    const double q1 = 2.356194491;
    const double q2 = 0.183105e-2, q3 = -0.3516396496e-4;
    const double q4 = 0.2457520174e-5, q5 = -0.240337019e-6;
    const double q6 = 0.04687499995, q7 = -0.2002690873e-3;
    const double q8 = 0.8449199096e-5, q9 = -0.88228987e-6;
    const double q10 = 0.105787412e-6, q11 = 0.636619772;

    if (x < 8) {
      y = x * x;
      result1 = x * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * p6)))));
      result2 = p7 + y * (p8 + y * (p9 + y * (p10 + y * (p11 + y * (p12 + y)))));
      result = (result1 / result2) + p13 * (xMath::BesselJ1(x)*log(x) - 1 / x);
    }
    else {
      z = 8 / x;
      y = z * z;
      xx = x - q1;
      result1 = 1 + y * (q2 + y * (q3 + y * (q4 + y * q5)));
      result2 = q6 + y * (q7 + y * (q8 + y * (q9 + y * q10)));
      result = sqrt(q11 / x)*(sin(xx)*result1 + z * cos(xx)*result2);
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::StruveH0(double x)
  {
    // Struve Functions of Order 0
    //
    // Converted from CERNLIB M342 by Rene Brun.

    const int n1 = 15;
    const int n2 = 25;
    const double c1[16] = { 1.00215845609911981, -1.63969292681309147,
                              1.50236939618292819, -.72485115302121872,
                               .18955327371093136, -.03067052022988,
                               .00337561447375194, -2.6965014312602e-4,
                              1.637461692612e-5,   -7.8244408508e-7,
                              3.021593188e-8,      -9.6326645e-10,
                              2.579337e-11,        -5.8854e-13,
                              1.158e-14,           -2e-16 };
    const double c2[26] = { .99283727576423943, -.00696891281138625,
                              1.8205103787037e-4,  -1.063258252844e-5,
                              9.8198294287e-7,     -1.2250645445e-7,
                              1.894083312e-8,      -3.44358226e-9,
                              7.1119102e-10,       -1.6288744e-10,
                              4.065681e-11,        -1.091505e-11,
                              3.12005e-12,         -9.4202e-13,
                              2.9848e-13,          -9.872e-14,
                              3.394e-14,           -1.208e-14,
                              4.44e-15,            -1.68e-15,
                              6.5e-16,             -2.6e-16,
                              1.1e-16,             -4e-17,
                              2e-17,               -1e-17 };

    const double c0 = 2 / xMath::Pi();

    int i;
    double alfa, h, r, y, b0, b1, b2;
    double v = fabs(x);

    v = fabs(x);
    if (v < 8) {
      y = v / 8;
      h = 2 * y*y - 1;
      alfa = h + h;
      b0 = 0;
      b1 = 0;
      b2 = 0;
      for (i = n1; i >= 0; --i) {
        b0 = c1[i] + alfa * b1 - b2;
        b2 = b1;
        b1 = b0;
      }
      h = y * (b0 - h * b2);
    }
    else {
      r = 1 / v;
      h = 128 * r*r - 1;
      alfa = h + h;
      b0 = 0;
      b1 = 0;
      b2 = 0;
      for (i = n2; i >= 0; --i) {
        b0 = c2[i] + alfa * b1 - b2;
        b2 = b1;
        b1 = b0;
      }
      h = xMath::BesselY0(v) + r * c0*(b0 - h * b2);
    }
    if (x < 0)  h = -h;
    return h;
  }

  //______________________________________________________________________________
  double xMath::StruveH1(double x)
  {
    // Struve Functions of Order 1
    //
    // Converted from CERNLIB M342 by Rene Brun.

    const int n3 = 16;
    const int n4 = 22;
    const double c3[17] = { .5578891446481605,   -.11188325726569816,
                             -.16337958125200939,   .32256932072405902,
                             -.14581632367244242,   .03292677399374035,
                             -.00460372142093573,  4.434706163314e-4,
                             -3.142099529341e-5,   1.7123719938e-6,
                             -7.416987005e-8,      2.61837671e-9,
                             -7.685839e-11,        1.9067e-12,
                             -4.052e-14,           7.5e-16,
                             -1e-17 };
    const double c4[23] = { 1.00757647293865641,  .00750316051248257,
                             -7.043933264519e-5,   2.66205393382e-6,
                             -1.8841157753e-7,     1.949014958e-8,
                             -2.6126199e-9,        4.236269e-10,
                             -7.955156e-11,        1.679973e-11,
                             -3.9072e-12,          9.8543e-13,
                             -2.6636e-13,          7.645e-14,
                             -2.313e-14,           7.33e-15,
                             -2.42e-15,            8.3e-16,
                             -3e-16,               1.1e-16,
                             -4e-17,               2e-17,-1e-17 };

    const double c0 = 2 / xMath::Pi();
    const double cc = 2 / (3 * xMath::Pi());

    int i, i1;
    double alfa, h, r, y, b0, b1, b2;
    double v = fabs(x);

    if (v == 0) {
      h = 0;
    }
    else if (v <= 0.3) {
      y = v * v;
      r = 1;
      h = 1;
      i1 = (int)(-8. / log10(v));
      for (i = 1; i <= i1; ++i) {
        h = -h * y / ((2 * i + 1)*(2 * i + 3));
        r += h;
      }
      h = cc * y*r;
    }
    else if (v < 8) {
      h = v * v / 32 - 1;
      alfa = h + h;
      b0 = 0;
      b1 = 0;
      b2 = 0;
      for (i = n3; i >= 0; --i) {
        b0 = c3[i] + alfa * b1 - b2;
        b2 = b1;
        b1 = b0;
      }
      h = b0 - h * b2;
    }
    else {
      h = 128 / (v*v) - 1;
      alfa = h + h;
      b0 = 0;
      b1 = 0;
      b2 = 0;
      for (i = n4; i >= 0; --i) {
        b0 = c4[i] + alfa * b1 - b2;
        b2 = b1;
        b1 = b0;
      }
      h = xMath::BesselY1(v) + c0 * (b0 - h * b2);
    }
    return h;
  }


  //______________________________________________________________________________
  double xMath::StruveL0(double x)
  {
    // Modified Struve Function of Order 0.
    // By Kirill Filimonov.

    const double pi = xMath::Pi();

    double s = 1.0;
    double r = 1.0;

    double a0, sl0, a1, bi0;

    int km, i;

    if (x <= 20.) {
      a0 = 2.0*x / pi;
      for (i = 1; i <= 60; i++) {
        r *= (x / (2 * i + 1))*(x / (2 * i + 1));
        s += r;
        if (fabs(r / s) < 1.e-12) break;
      }
      sl0 = a0 * s;
    }
    else {
      km = int(5 * (x + 1.0));
      if (x >= 50.0)km = 25;
      for (i = 1; i <= km; i++) {
        r *= (2 * i - 1)*(2 * i - 1) / x / x;
        s += r;
        if (fabs(r / s) < 1.0e-12) break;
      }
      a1 = exp(x) / sqrt(2 * pi*x);
      r = 1.0;
      bi0 = 1.0;
      for (i = 1; i <= 16; i++) {
        r = 0.125*r*(2.0*i - 1.0)*(2.0*i - 1.0) / (i*x);
        bi0 += r;
        if (fabs(r / bi0) < 1.0e-12) break;
      }

      bi0 = a1 * bi0;
      sl0 = -2.0 / (pi*x)*s + bi0;
    }
    return sl0;
  }

  //______________________________________________________________________________
  double xMath::StruveL1(double x)
  {
    // Modified Struve Function of Order 1.
    // By Kirill Filimonov.

    const double pi = xMath::Pi();
    double a1, sl1, bi1, s;
    double r = 1.0;
    int km, i;

    if (x <= 20.) {
      s = 0.0;
      for (i = 1; i <= 60; i++) {
        r *= x * x / (4.0*i*i - 1.0);
        s += r;
        if (fabs(r) < fabs(s)*1.e-12) break;
      }
      sl1 = 2.0 / pi * s;
    }
    else {
      s = 1.0;
      km = int(0.5*x);
      if (x > 50.0)km = 25;
      for (i = 1; i <= km; i++) {
        r *= (2 * i + 3)*(2 * i + 1) / x / x;
        s += r;
        if (fabs(r / s) < 1.0e-12) break;
      }
      sl1 = 2.0 / pi * (-1.0 + 1.0 / (x*x) + 3.0*s / (x*x*x*x));
      a1 = exp(x) / sqrt(2 * pi*x);
      r = 1.0;
      bi1 = 1.0;
      for (i = 1; i <= 16; i++) {
        r = -0.125*r*(4.0 - (2.0*i - 1.0)*(2.0*i - 1.0)) / (i*x);
        bi1 += r;
        if (fabs(r / bi1) < 1.0e-12) break;
      }
      sl1 += a1 * bi1;
    }
    return sl1;
  }



  //______________________________________________________________________________
  double xMath::BesselK0exp(double x)
  {
    // Compute the modified Bessel function K_0(x) for positive real x.
    //
    //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
    //     Applied Mathematics Series vol. 55 (1964), Washington.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht

    // Parameters of the polynomial approximation
    const double p1 = -0.57721566, p2 = 0.42278420, p3 = 0.23069756,
      p4 = 3.488590e-2, p5 = 2.62698e-3, p6 = 1.0750e-4, p7 = 7.4e-6;

    const double q1 = 1.25331414, q2 = -7.832358e-2, q3 = 2.189568e-2,
      q4 = -1.062446e-2, q5 = 5.87872e-3, q6 = -2.51540e-3, q7 = 5.3208e-4;

    if (x <= 0) {
      //Error("xMath::BesselK0", "*K0* Invalid argument x = %g\n",x);
      return 0;
    }

    double y = 0, result = 0;

    if (x <= 2) {
      y = x * x / 4;
      result = exp(x)*((-log(x / 2.)*xMath::BesselI0(x)) + (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7)))))));
    }
    else {
      y = 2 / x;
      result = (1. / sqrt(x))*(q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * q7))))));
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselK1exp(double x)
  {
    // Compute the modified Bessel function K_1(x) for positive real x.
    //
    //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
    //     Applied Mathematics Series vol. 55 (1964), Washington.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht

    // Parameters of the polynomial approximation
    const double p1 = 1., p2 = 0.15443144, p3 = -0.67278579,
      p4 = -0.18156897, p5 = -1.919402e-2, p6 = -1.10404e-3, p7 = -4.686e-5;

    const double q1 = 1.25331414, q2 = 0.23498619, q3 = -3.655620e-2,
      q4 = 1.504268e-2, q5 = -7.80353e-3, q6 = 3.25614e-3, q7 = -6.8245e-4;

    if (x <= 0) {
      //Error("xMath::BesselK1", "*K1* Invalid argument x = %g\n",x);
      return 0;
    }

    double y = 0, result = 0;

    if (x <= 2) {
      y = x * x / 4;
      result = exp(x) * ((log(x / 2.)*xMath::BesselI1(x)) + (1. / x)*(p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7)))))));
    }
    else {
      y = 2 / x;
      result = (1. / sqrt(x))*(q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * q7))))));
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselKexp(int n, double x)
  {
    if (n < 0) n = -n;
    // Compute the Integer Order Modified Bessel function K_n(x)
     // for n=0,1,2,... and positive real x.
     //
     //--- NvE 12-mar-2000 UU-SAP Utrecht

    if (x <= 0 || n < 0) {
      //Error("xMath::BesselK", "*K* Invalid argument(s) (n,x) = (%d, %g)\n",n,x);
      return 0;
    }

    if (n == 0) return xMath::BesselK0exp(x);
    if (n == 1) return xMath::BesselK1exp(x);

    // Perform upward recurrence for all x
    double tox = 2 / x;
    double bkm = xMath::BesselK0exp(x);
    double bk = xMath::BesselK1exp(x);
    double bkp = 0;
    for (int j = 1; j < n; j++) {
      bkp = bkm + double(j)*tox*bk;
      bkm = bk;
      bk = bkp;
    }
    return bk;
  }

  //______________________________________________________________________________
  double xMath::BesselI0exp(double x)
  {
    // Compute the modified Bessel function I_0(x) for any real x.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht

    // Parameters of the polynomial approximation
    const double p1 = 1.0, p2 = 3.5156229, p3 = 3.0899424,
      p4 = 1.2067492, p5 = 0.2659732, p6 = 3.60768e-2, p7 = 4.5813e-3;

    const double q1 = 0.39894228, q2 = 1.328592e-2, q3 = 2.25319e-3,
      q4 = -1.57565e-3, q5 = 9.16281e-3, q6 = -2.057706e-2,
      q7 = 2.635537e-2, q8 = -1.647633e-2, q9 = 3.92377e-3;

    const double k1 = 3.75;
    double ax = fabs(x);

    double y = 0, result = 0;

    if (ax < k1) {
      double xx = x / k1;
      y = xx * xx;
      result = exp(-ax) * ( p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7))))) );
    }
    else {
      y = k1 / ax;
      result = (1. / sqrt(ax))*(q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))))))));
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselI1exp(double x)
  {
    // Compute the modified Bessel function I_1(x) for any real x.
    //
    //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
    //     Applied Mathematics Series vol. 55 (1964), Washington.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht

    // Parameters of the polynomial approximation
    const double p1 = 0.5, p2 = 0.87890594, p3 = 0.51498869,
      p4 = 0.15084934, p5 = 2.658733e-2, p6 = 3.01532e-3, p7 = 3.2411e-4;

    const double q1 = 0.39894228, q2 = -3.988024e-2, q3 = -3.62018e-3,
      q4 = 1.63801e-3, q5 = -1.031555e-2, q6 = 2.282967e-2,
      q7 = -2.895312e-2, q8 = 1.787654e-2, q9 = -4.20059e-3;

    const double k1 = 3.75;
    double ax = fabs(x);

    double y = 0, result = 0;

    if (ax < k1) {
      double xx = x / k1;
      y = xx * xx;
      result = exp(-ax) * (x * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7)))))));
    }
    else {
      y = k1 / ax;
      result = (1. / sqrt(ax))*(q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))))))));
      if (x < 0) result = -result;
    }
    return result;
  }

  //______________________________________________________________________________
  double xMath::BesselIexp(int n, double x)
  {
    // Compute the Integer Order Modified Bessel function I_n(x)
    // for n=0,1,2,... and any real x.
    //
    //--- NvE 12-mar-2000 UU-SAP Utrecht
    if (n < 0) n = -n;

    int iacc = 40; // Increase to enhance accuracy
    const double kBigPositive = 1.e10;
    const double kBigNegative = 1.e-10;

    if (n < 0) {
      //Error("xMath::BesselI", "*I* Invalid argument (n,x) = (%d, %g)\n",n,x);
      return 0;
    }

    if (n == 0) return xMath::BesselI0exp(x);
    if (n == 1) return xMath::BesselI1exp(x);

    if (x == 0) return 0;
    if (fabs(x) > kBigPositive) return 0;

    double tox = 2 / fabs(x);
    double bip = 0, bim = 0;
    double bi = 1;
    double result = 0;
    int m = 2 * ((n + int(sqrt(float(iacc*n)))));
    for (int j = m; j >= 1; j--) {
      bim = bip + double(j)*tox*bi;
      bip = bi;
      bi = bim;
      // Renormalise to prevent overflows
      if (fabs(bi) > kBigPositive) {
        result *= kBigNegative;
        bi *= kBigNegative;
        bip *= kBigNegative;
      }
      if (j == n) result = bip;
    }

    result *= xMath::BesselI0exp(x) / bi; // Normalise with BesselI0exp(x)
    if ((x < 0) && (n % 2 == 1)) result = -result;

    return result;
  }


  double xMath::Gamma
  (
    double x    // We require x > 0
  )
  {
    if (x <= 0.0)
    {
      std::stringstream os;
      os << "Invalid input argument " << x << ". Argument must be positive.";
      throw std::invalid_argument(os.str());
    }

    // Split the function domain into three intervals:
    // (0, 0.001), [0.001, 12), and (12, infinity)

    ///////////////////////////////////////////////////////////////////////////
    // First interval: (0, 0.001)
  //
  // For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
  // So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
  // The relative error over this interval is less than 6e-7.

    const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

    if (x < 0.001)
      return 1.0 / (x*(1.0 + gamma * x));

    ///////////////////////////////////////////////////////////////////////////
    // Second interval: [0.001, 12)

    if (x < 12.0)
    {
      // The algorithm directly approximates gamma over (1,2) and uses
      // reduction identities to reduce other arguments to this interval.

      double y = x;
      int n = 0;
      bool arg_was_less_than_one = (y < 1.0);

      // Add or subtract integers as necessary to bring y into (1,2)
      // Will correct for this below
      if (arg_was_less_than_one)
      {
        y += 1.0;
      }
      else
      {
        n = static_cast<int> (floor(y)) - 1;  // will use n later
        y -= n;
      }

      // numerator coefficients for approximation over the interval (1,2)
      static const double p[] =
      {
          -1.71618513886549492533811E+0,
           2.47656508055759199108314E+1,
          -3.79804256470945635097577E+2,
           6.29331155312818442661052E+2,
           8.66966202790413211295064E+2,
          -3.14512729688483675254357E+4,
          -3.61444134186911729807069E+4,
           6.64561438202405440627855E+4
      };

      // denominator coefficients for approximation over the interval (1,2)
      static const double q[] =
      {
          -3.08402300119738975254353E+1,
           3.15350626979604161529144E+2,
          -1.01515636749021914166146E+3,
          -3.10777167157231109440444E+3,
           2.25381184209801510330112E+4,
           4.75584627752788110767815E+3,
          -1.34659959864969306392456E+5,
          -1.15132259675553483497211E+5
      };

      double num = 0.0;
      double den = 1.0;
      int i;

      double z = y - 1;
      for (i = 0; i < 8; i++)
      {
        num = (num + p[i])*z;
        den = den * z + q[i];
      }
      double result = num / den + 1.0;

      // Apply correction if argument was not initially in (1,2)
      if (arg_was_less_than_one)
      {
        // Use identity gamma(z) = gamma(z+1)/z
        // The variable "result" now holds gamma of the original y + 1
        // Thus we use y-1 to get back the orginal y.
        result /= (y - 1.0);
      }
      else
      {
        // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
        for (i = 0; i < n; i++)
          result *= y++;
      }

      return result;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Third interval: [12, infinity)

    if (x > 171.624)
    {
      // Correct answer too large to display. Force +infinity.
      double temp = DBL_MAX;
      return temp * 2.0;
    }

    return exp(LogGamma(x));
  }

  double xMath::LogGamma
  (
    double x    // x must be positive
  )
  {
    if (x <= 0.0)
    {
      std::stringstream os;
      os << "Invalid input argument " << x << ". Argument must be positive.";
      throw std::invalid_argument(os.str());
    }

    if (x < 12.0)
    {
      return log(fabs(Gamma(x)));
    }

    // Abramowitz and Stegun 6.1.41
      // Asymptotic series should be good to at least 11 or 12 figures
      // For error analysis, see Whittiker and Watson
      // A Course in Modern Analysis (1927), page 252

    static const double c[8] =
    {
     1.0 / 12.0,
    -1.0 / 360.0,
    1.0 / 1260.0,
    -1.0 / 1680.0,
    1.0 / 1188.0,
    -691.0 / 360360.0,
    1.0 / 156.0,
    -3617.0 / 122400.0
    };
    double z = 1.0 / (x*x);
    double sum = c[7];
    for (int i = 6; i >= 0; i--)
    {
      sum *= z;
      sum += c[i];
    }
    double series = sum / x;

    static const double halfLogTwoPi = 0.91893853320467274178032973640562;
    double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;
    return logGamma;
  }

} // namespace thermalfist