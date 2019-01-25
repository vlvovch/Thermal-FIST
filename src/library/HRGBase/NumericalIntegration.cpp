/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/NumericalIntegration.h"

namespace thermalfist {

  namespace NumericalIntegration {

    double Integrate2DLaguerre32Legendre32(double(*func)(double, double), double ay, double by)
    {
      // Integrate 2D function from 0 to infinity using Gauss-Laguerre integration
      // with 32 points and from ay to by using Gauss-Legendre integration with 32 points
      //

      const double *xleg = coefficients_xleg32;
      const double *wleg = coefficients_wleg32;
      double x[32];
      double w[32];


      const double *xlag = coefficients_xlag32;
      const double *wlag = coefficients_wlag32;

      double sum = 0.;

      for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) {
          x[j] = (by - ay) / 2.*xleg[j] + (by + ay) / 2.;
          w[j] = (by - ay) / 2.*wleg[j];

          sum += wlag[i] * w[j] * func(xlag[i], x[j]);
        }
      }
      return sum;
    }

    void GetCoefs2DLaguerre32Legendre32(double ay, double by,
      std::vector<double> *xlagp, std::vector<double> *wlagp,
      std::vector<double> *xlegp, std::vector<double> *wlegp) {
      std::vector<double> &xlag = *xlagp;
      std::vector<double> &wlag = *wlagp;
      std::vector<double> &xleg = *xlegp;
      std::vector<double> &wleg = *wlegp;

      xlag.resize(32);
      wlag.resize(32);
      xleg.resize(32);
      wleg.resize(32);

      const double *xlego = coefficients_xleg32;
      const double *wlego = coefficients_wleg32;
      const double *xlago = coefficients_xlag32;
      const double *wlago = coefficients_wlag32;

      for (int j = 0; j < 32; j++) {
        xleg[j] = (by - ay) / 2.*xlego[j] + (by + ay) / 2.;
        wleg[j] = (by - ay) / 2.*wlego[j];
        xlag[j] = xlago[j];
        wlag[j] = wlago[j];
      }
    }

    void GetCoefs2DLegendre32Legendre32(double ay, double by, double a2y, double b2y,
      std::vector<double> *xlegp1, std::vector<double> *wlegp1,
      std::vector<double> *xlegp2, std::vector<double> *wlegp2) {
      std::vector<double> &xleg1 = *xlegp1;
      std::vector<double> &wleg1 = *wlegp1;
      std::vector<double> &xleg2 = *xlegp2;
      std::vector<double> &wleg2 = *wlegp2;

      xleg1.resize(32);
      wleg1.resize(32);
      xleg2.resize(32);
      wleg2.resize(32);

      const double *xlego = coefficients_xleg32;
      const double *wlego = coefficients_wleg32;

      for (int j = 0; j < 32; j++) {
        xleg1[j] = (by - ay) / 2.*xlego[j] + (by + ay) / 2.;
        wleg1[j] = (by - ay) / 2.*wlego[j];
        xleg2[j] = (b2y - a2y) / 2.*xlego[j] + (b2y + a2y) / 2.;
        wleg2[j] = (b2y - a2y) / 2.*wlego[j];
      }
    }

    void GetCoefsIntegrateLegendre32(double a, double b, std::vector<double> *xp, std::vector<double> *wp)
    {
      // Integrate function from a to b using Legendre-Gaussian integration
      // with 32 points.
      //
      std::vector<double> &x = *xp;
      std::vector<double> &w = *wp;

      x.resize(32);
      w.resize(32);

      const double *xlego = coefficients_xleg32;
      const double *wlego = coefficients_wleg32;


      for (int i = 0; i < 32; i++) {
        w[i] = (b - a) / 2.*wlego[i];
        x[i] = (b - a) / 2.*xlego[i] + (b + a) / 2.;
      }
    }

    void GetCoefsIntegrateLegendre10(double a, double b, std::vector<double> *xp, std::vector<double> *wp)
    {
      // Integrate function from a to b using Legendre-Gaussian integration
      // with 32 points.
      //
      std::vector<double> &x = *xp;
      std::vector<double> &w = *wp;

      x.resize(10);
      w.resize(10);

      const double *xlego = coefficients_xleg10;
      const double *wlego = coefficients_wleg10;

      for (int i = 0; i < 10; i++) {
        w[i] = (b - a) / 2.*wlego[i];
        x[i] = (b - a) / 2.*xlego[i] + (b + a) / 2.;
      }
    }

    void GetCoefsIntegrateLegendre5(double a, double b, std::vector<double> *xp, std::vector<double> *wp)
    {
      // Integrate function from a to b using Legendre-Gaussian integration
      // with 32 points.
      //
      std::vector<double> &x = *xp;
      std::vector<double> &w = *wp;

      x.resize(5);
      w.resize(5);

      const double *xlego = coefficients_xleg5;
      const double *wlego = coefficients_wleg5;

      for (int i = 0; i < 5; i++) {
        w[i] = (b - a) / 2.*wlego[i];
        x[i] = (b - a) / 2.*xlego[i] + (b + a) / 2.;
      }
    }

    void GetCoefsIntegrateLegendre40(double a, double b, std::vector<double> *xp, std::vector<double> *wp)
    {
      // Integrate function from a to b using Legendre-Gaussian integration
      // with 40 points.
      //
      std::vector<double> &x = *xp;
      std::vector<double> &w = *wp;

      x.resize(40);
      w.resize(40);


      const double *xlego = coefficients_xleg40;
      const double *wlego = coefficients_wleg40;

      for (int i = 0; i < 40; i++) {
        w[i] = (b - a) / 2.*wlego[i];
        x[i] = (b - a) / 2.*xlego[i] + (b + a) / 2.;
      }
    }

    void GetCoefsIntegrateLaguerre32(std::vector<double> *xp, std::vector<double> *wp)
    {
      // Integrate function from 0 to infinity using Gauss-Laguerre integration
      // with 32 points
      //
      std::vector<double> &x = *xp;
      std::vector<double> &w = *wp;

      x.resize(32);
      w.resize(32);

      const double *xlago = coefficients_xlag32;
      const double *wlago = coefficients_wlag32;

      for (int i = 0; i < 32; i++) {
        w[i] = wlago[i];
        x[i] = xlago[i];
      }
    }

  } // namespace NumericalIntegration

} // namespace thermalfist