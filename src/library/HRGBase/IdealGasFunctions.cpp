/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2017-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/IdealGasFunctions.h"

#include <stdio.h>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cfloat>
#include <vector>

#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"

using namespace std;

namespace thermalfist {

  namespace IdealGasFunctions {

    /// Calculate all the spin values -S, S+1,...,S-1,S
    std::vector<double> GetSpins(double degSpin) {
      // Check that spin degeneracy is integer
      assert(abs(degSpin - round(degSpin)) < 1.e-6);
      int spinDegeneracy = round(degSpin);
      std::vector<double> ret;
      for(int isz = 0; isz < spinDegeneracy; ++isz) {
        ret.push_back(-(spinDegeneracy - 1.)/2. + isz);
      }
      return ret;
    }

    bool calculationHadBECIssue = false;

    double BoltzmannDensity(double T, double mu, double m, double deg,
                            const IdealGasFunctionsExtraConfig& extraConfig) {
      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);
        // Use Eq. (7) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));
            ret += e0 * xMath::BesselKexp(1, e0 / T) * exp((mu - e0) / T);
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      if (m == 0.)
        return deg * T * T * T / 2. / xMath::Pi() / xMath::Pi() * 2. * exp(mu/ T) * xMath::GeVtoifm3();
      return deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::BesselKexp(2, m / T) * exp((mu - m) / T) * xMath::GeVtoifm3();
    }

    double BoltzmannPressure(double T, double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig) {
      return T * BoltzmannDensity(T, mu, m, deg, extraConfig);
    }

    double BoltzmannEnergyDensity(double T, double mu, double m, double deg,
                                  const IdealGasFunctionsExtraConfig& extraConfig) {
      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        // TODO: Check that it's done correctly
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);
        // Use e = -p + T (dp/dT) + mu (dp/mu)
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));
            ret += (T + e0 * xMath::BesselKexp(0, e0 / T) / xMath::BesselKexp(1, e0 / T)) * e0 * xMath::BesselKexp(1, e0 / T) * exp((mu - e0) / T);
          }
        }
        ret *= Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
        return ret;
      }

      // No magnetic field
      if (m == 0.)
        return 3 * T * BoltzmannDensity(T, mu, m, deg);
      return (3 * T + m * xMath::BesselK1exp(m / T) / xMath::BesselKexp(2, m / T)) * BoltzmannDensity(T, mu, m, deg);
    }

    double BoltzmannEntropyDensity(double T, double mu, double m, double deg,
                                   const IdealGasFunctionsExtraConfig& extraConfig) {
      return (BoltzmannPressure(T, mu, m, deg, extraConfig) +
      BoltzmannEnergyDensity(T, mu, m, deg, extraConfig) -
      mu * BoltzmannDensity(T, mu, m, deg, extraConfig)) / T;
    }

    double BoltzmannScalarDensity(double T, double mu, double m, double deg,
                                  const IdealGasFunctionsExtraConfig& extraConfig) {
      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        // TODO: Check that it's done correctly
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);
        // Use nsig = -dp/dm
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));
            ret += m * xMath::BesselKexp(0, e0 / T) * exp((mu - e0) / T);
          }
        }
        ret *= Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
        return ret;
      }

      // No magnetic field
      if (m == 0.)
        return 0.;
      return deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::BesselKexp(1, m / T) * exp((mu - m) / T) * xMath::GeVtoifm3();
    }

    double BoltzmannTdndmu(int /*N*/, double T, double mu, double m, double deg,
                           const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return BoltzmannDensity(T, mu, m, deg, extraConfig);
    }

    double BoltzmannChiN(int N, double T, double mu, double m, double deg,
                         const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return BoltzmannTdndmu(N - 1, T, mu, m, deg, extraConfig) / pow(T, 3) / xMath::GeVtoifm3();
    }

    double BoltzmannChiNDimensionfull(int N, double T, double mu, double m, double deg,
                                      const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return BoltzmannTdndmu(N - 1, T, mu, m, deg, extraConfig) / pow(T, N - 1) / xMath::GeVtoifm3();
    }

    double QuantumClusterExpansionDensity(int statistics, double T, double mu, double m, double deg, int order,
                                          const IdealGasFunctionsExtraConfig& extraConfig)
    {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return BoltzmannDensity(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);


        // Use Eq. (7) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));

            // Sum over clusters
            double tfug = exp((mu - e0) / T);
            double EoverT = e0 / T;
            double cfug = tfug;
            double sign = 1.;
            for (int i = 1; i <= order; ++i) {
              ret += e0 * sign * xMath::BesselKexp(1, i*EoverT) * cfug;
              cfug *= tfug;
              if (signchange) sign = -sign;
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double moverT = m / T;
      double cfug = tfug;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        ret += sign * xMath::BesselKexp(2, i*moverT) * cfug / static_cast<double>(i);
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      ret *= deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      return ret;
    }

    double QuantumClusterExpansionPressure(int statistics, double T, double mu, double m, double deg, int order,
                                           const IdealGasFunctionsExtraConfig& extraConfig)
    {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return BoltzmannPressure(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);


        // Use Eq. (7) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));

            // Sum over clusters
            double tfug = exp((mu - e0) / T);
            double EoverT = e0 / T;
            double cfug = tfug;
            double sign = 1.;
            for (int i = 1; i <= order; ++i) {
              ret += e0 * sign * xMath::BesselKexp(1, i*EoverT) * cfug * T / static_cast<double>(i);
              cfug *= tfug;
              if (signchange) sign = -sign;
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double cfug = tfug;
      double moverT = m / T;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        ret += sign * xMath::BesselKexp(2, i*moverT) * cfug / static_cast<double>(i) / static_cast<double>(i);
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      ret *= deg * m * m * T * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      return ret;
    }

    double QuantumClusterExpansionEnergyDensity(int statistics, double T, double mu, double m, double deg, int order,
                                                const IdealGasFunctionsExtraConfig& extraConfig)
    {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return BoltzmannEnergyDensity(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);


        // Use Eq. (7) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));

            // Sum over clusters
            double tfug = exp((mu - e0) / T);
            double EoverT = e0 / T;
            double cfug = tfug;
            double sign = 1.;
            for (int i = 1; i <= order; ++i) {
              ret += sign * e0 * T / static_cast<double>(i)  * xMath::BesselKexp(1, i*EoverT) * cfug;
              ret += sign * e0 * e0 * xMath::BesselKexp(0, i*EoverT) * cfug;
              cfug *= tfug;
              if (signchange) sign = -sign;
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double cfug = tfug;
      double moverT = m / T;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        ret += sign * (xMath::BesselKexp(1, i*moverT) + 3. * xMath::BesselKexp(2, i*moverT) / moverT / static_cast<double>(i)) * cfug / static_cast<double>(i);
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      ret *= deg * m * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      return ret;
    }

    double QuantumClusterExpansionEntropyDensity(int statistics, double T, double mu, double m, double deg, int order,
                                                 const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return (QuantumClusterExpansionPressure(statistics, T, mu, m, deg, order, extraConfig) +
      QuantumClusterExpansionEnergyDensity(statistics, T, mu, m, deg, order, extraConfig) -
      mu * QuantumClusterExpansionDensity(statistics, T, mu, m, deg, order, extraConfig)) / T;
    }

    double QuantumClusterExpansionScalarDensity(int statistics, double T, double mu, double m, double deg, int order,
                                                const IdealGasFunctionsExtraConfig& extraConfig)
    {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return BoltzmannScalarDensity(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);


        // Use Eq. (7) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));

            // Sum over clusters
            double tfug = exp((mu - e0) / T);
            double EoverT = e0 / T;
            double cfug = tfug;
            double sign = 1.;
            for (int i = 1; i <= order; ++i) {
              ret += m * sign * xMath::BesselKexp(0, i*EoverT) * cfug;
              cfug *= tfug;
              if (signchange) sign = -sign;
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double cfug = tfug;
      double moverT = m / T;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        ret += sign * xMath::BesselKexp(1, i*moverT) * cfug / static_cast<double>(i);
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      ret *= deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      return ret;
    }

    double QuantumClusterExpansionTdndmu(int N, int statistics, double T, double mu, double m, double deg, int order,
                                         const IdealGasFunctionsExtraConfig& extraConfig)
    {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return BoltzmannTdndmu(N, T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);


        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));

            // Sum over clusters
            double tfug = exp((mu - e0) / T);
            double EoverT = e0 / T;
            double cfug = tfug;
            double sign = 1.;
            for (int i = 1; i <= order; ++i) {
              ret += e0 * sign * xMath::BesselKexp(1, i*EoverT) * cfug * pow(static_cast<double>(i), N);
              cfug *= tfug;
              if (signchange) sign = -sign;
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double cfug = tfug;
      double moverT = m / T;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        ret += sign * xMath::BesselKexp(2, i*moverT) * cfug * pow(static_cast<double>(i), N - 1);
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      ret *= deg * m * m * T / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      return ret;
    }

    double QuantumClusterExpansionChiN(int N, int statistics, double T, double mu, double m, double deg, int order,
                                       const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return QuantumClusterExpansionTdndmu(N - 1, statistics, T, mu, m, deg, order, extraConfig) /
        pow(T, 3) / xMath::GeVtoifm3();
    }


    double QuantumClusterExpansionChiNDimensionfull(int N, int statistics, double T, double mu, double m, double deg, int order,
                                                    const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return QuantumClusterExpansionTdndmu(N - 1, statistics, T, mu, m, deg, order, extraConfig) /
      pow(T, N - 1) / xMath::GeVtoifm3();
    }

    // Gauss-Legendre 32-point quadrature for [0,1] interval
    const double* legx32 = NumericalIntegration::coefficients_xleg32_zeroone;
    const double* legw32 = NumericalIntegration::coefficients_wleg32_zeroone;
    // Gauss-Laguerre 32-point quadrature for [0,\infty] interval
    const double* lagx32 = NumericalIntegration::coefficients_xlag32;
    const double* lagw32 = NumericalIntegration::coefficients_wlag32;

    double QuantumNumericalIntegrationDensity(int statistics, double T, double mu, double m, double deg,
                                              const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0)           return BoltzmannDensity(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return FermiZeroTDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuDensity(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        printf("**WARNING** QuantumNumericalIntegrationDensity: Bose-Einstein condensation, mass = %lf, mu = %lf\n", m, mu);
        calculationHadBECIssue = true;
        return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);

        // Use Eq. (5) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {

            // Numerical integration over pz
            double moverT = m / T;
            double muoverT = mu / T;
            double BoverT2 = extraConfig.MagneticField.B / T / T;
            for (int i = 0; i < 32; i++) {
              double tx = lagx32[i];
              double EoverT = sqrt(tx*tx + moverT*moverT + Qmod*BoverT2*(2.*l + 1. - 2.*sz));
              ret += lagw32[i] * T / (exp(EoverT - muoverT) + statistics);
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        ret += lagw32[i] * T * tx * T * tx * T / (exp(sqrt(tx*tx + moverT * moverT) - muoverT) + statistics);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationPressure(int statistics, double T, double mu, double m, double deg,
                                               const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0)           return BoltzmannPressure(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return FermiZeroTPressure(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuPressure(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        printf("**WARNING** QuantumNumericalIntegrationPressure: Bose-Einstein condensation\n");
        calculationHadBECIssue = true;
        return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);

        // Use Eq. (5) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {

            // Numerical integration over pz
            double moverT = m / T;
            double muoverT = mu / T;
            double BoverT2 = extraConfig.MagneticField.B / T / T;
            for (int i = 0; i < 32; i++) {
              double tx = lagx32[i];
              double x2 = tx * T * tx * T;
              double EoverT = sqrt(tx*tx + moverT*moverT + Qmod*BoverT2*(2.*l + 1. - 2.*sz));
              ret += lagw32[i] * T * x2 / EoverT / T / (exp(EoverT - muoverT) + statistics);
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double x2 = tx * T * tx * T;
        double E = sqrt(tx*tx + moverT * moverT);
        ret += lagw32[i] * T * x2 * x2 / E / T / (exp(E - muoverT) + statistics);
      }

      ret *= deg / 6. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationEnergyDensity(int statistics, double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0)           return BoltzmannEnergyDensity(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return FermiZeroTEnergyDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuEnergyDensity(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        printf("**WARNING** QuantumNumericalIntegrationEnergyDensity: Bose-Einstein condensation\n");
        calculationHadBECIssue = true;
        return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);

        // Use Eq. (5) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {

            // Numerical integration over pz
            double moverT = m / T;
            double muoverT = mu / T;
            double BoverT2 = extraConfig.MagneticField.B / T / T;
            for (int i = 0; i < 32; i++) {
              double tx = lagx32[i];
              double EoverT = sqrt(tx*tx + moverT*moverT + Qmod*BoverT2*(2.*l + 1. - 2.*sz));
              ret += lagw32[i] * T * EoverT * T / (exp(EoverT - muoverT) + statistics);
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        ret += lagw32[i] * T * tx * T * tx * T * sqrt(tx*tx + moverT * moverT) * T / (exp(sqrt(tx*tx + moverT * moverT) - muoverT) + statistics);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationEntropyDensity(int statistics, double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (T == 0.)
        return 0.;
      return (QuantumNumericalIntegrationPressure(statistics, T, mu, m, deg, extraConfig)
      + QuantumNumericalIntegrationEnergyDensity(statistics, T, mu, m, deg, extraConfig)
      - mu * QuantumNumericalIntegrationDensity(statistics, T, mu, m, deg, extraConfig)) / T;
    }

    double QuantumNumericalIntegrationScalarDensity(int statistics, double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0)           return BoltzmannScalarDensity(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return FermiZeroTScalarDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuScalarDensity(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        printf("**WARNING** QuantumNumericalIntegrationScalarDensity: Bose-Einstein condensation\n");
        calculationHadBECIssue = true;
        return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);

        // Use Eq. (5) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {

            // Numerical integration over pz
            double moverT = m / T;
            double muoverT = mu / T;
            double BoverT2 = extraConfig.MagneticField.B / T / T;
            for (int i = 0; i < 32; i++) {
              double tx = lagx32[i];
              double EoverT = sqrt(tx*tx + moverT*moverT + Qmod*BoverT2*(2.*l + 1. - 2.*sz));
              ret += lagw32[i] * T * moverT / EoverT / (exp(EoverT - muoverT) + statistics);
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        ret += lagw32[i] * T * tx * T * tx * T * moverT / sqrt(tx*tx + moverT * moverT) / (exp(sqrt(tx*tx + moverT * moverT) - muoverT) + statistics);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationT1dn1dmu1(int statistics, double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0)            return BoltzmannTdndmu(1, T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return 0.;
      if (statistics == 1 && mu > m)  return FermiNumericalIntegrationLargeMuT1dn1dmu1(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        printf("**WARNING** QuantumNumericalIntegrationT1dn1dmu1: Bose-Einstein condensation\n");
        calculationHadBECIssue = true;
        return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);

        // Use Eq. (5) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {

            // Numerical integration over pz
            double moverT = m / T;
            double muoverT = mu / T;
            double BoverT2 = extraConfig.MagneticField.B / T / T;
            for (int i = 0; i < 32; i++) {
              double tx = lagx32[i];
              double EoverT = sqrt(tx*tx + moverT*moverT + Qmod*BoverT2*(2.*l + 1. - 2.*sz));
              double Eexp = exp(EoverT - muoverT);
              ret += lagw32[i] * T * 1. / (1. + statistics / Eexp) / (Eexp + statistics);
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double Eexp = exp(sqrt(tx*tx + moverT * moverT) - muoverT);
        ret += lagw32[i] * T * tx * T * tx * T * 1. / (1. + statistics / Eexp) / (Eexp + statistics);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationT2dn2dmu2(int statistics, double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0)            return BoltzmannTdndmu(2, T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return 0.;
      if (statistics == 1 && mu > m)  return FermiNumericalIntegrationLargeMuT2dn2dmu2(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        printf("**WARNING** QuantumNumericalIntegrationT2dn2dmu2: Bose-Einstein condensation\n");
        calculationHadBECIssue = true;
        return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);

        // Use Eq. (5) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {

            // Numerical integration over pz
            double moverT = m / T;
            double muoverT = mu / T;
            double BoverT2 = extraConfig.MagneticField.B / T / T;
            for (int i = 0; i < 32; i++) {
              double tx = lagx32[i];
              double EoverT = sqrt(tx*tx + moverT*moverT + Qmod*BoverT2*(2.*l + 1. - 2.*sz));
              double Eexp = exp(EoverT - muoverT);
              ret += lagw32[i] * T * (1. - statistics / Eexp) / (1. + statistics / Eexp) / (1. + statistics / Eexp) / (Eexp + statistics);
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double Eexp = exp(sqrt(tx*tx + moverT * moverT) - muoverT);
        //ret += lagw32[i] * T * tx * T * tx * T * (Eexp*Eexp - statistics * Eexp) / (Eexp + statistics) / (Eexp + statistics) / (Eexp + statistics);
        ret += lagw32[i] * T * tx * T * tx * T * (1. - statistics / Eexp) / (1. + statistics / Eexp) / (1. + statistics / Eexp) / (Eexp + statistics);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationT3dn3dmu3(int statistics, double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0)            return BoltzmannTdndmu(3, T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return 0.;
      if (statistics == 1 && mu > m)  return FermiNumericalIntegrationLargeMuT3dn3dmu3(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        printf("**WARNING** QuantumNumericalIntegrationT3dn3dmu3: Bose-Einstein condensation\n");
        calculationHadBECIssue = true;
        return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;


      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);

        // Use Eq. (5) from https://arxiv.org/pdf/2104.06843.pdf
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {

            // Numerical integration over pz
            double moverT = m / T;
            double muoverT = mu / T;
            double BoverT2 = extraConfig.MagneticField.B / T / T;
            for (int i = 0; i < 32; i++) {
              double tx = lagx32[i];
              double EoverT = sqrt(tx*tx + moverT*moverT + Qmod*BoverT2*(2.*l + 1. - 2.*sz));
              double Eexp = exp(EoverT - muoverT);
              ret += lagw32[i] * T * (1. - 4.*statistics / Eexp + statistics * statistics / Eexp / Eexp) / (1. + statistics / Eexp) / (1. + statistics / Eexp) / (1. + statistics / Eexp) / (Eexp + statistics);
            }
          }
        }
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
      }

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double Eexp = exp(sqrt(tx*tx + moverT * moverT) - muoverT);
        //ret += lagw32[i] * T * tx * T * tx * T * (Eexp*Eexp*Eexp - 4.*statistics*Eexp*Eexp + statistics*statistics*Eexp) / (Eexp + statistics) / (Eexp + statistics) / (Eexp + statistics) / (Eexp + statistics);
        ret += lagw32[i] * T * tx * T * tx * T * (1. - 4.*statistics / Eexp + statistics * statistics / Eexp / Eexp) / (1. + statistics / Eexp) / (1. + statistics / Eexp) / (1. + statistics / Eexp) / (Eexp + statistics);
        //printf("%E %E  ", ret, Eexp);
      }
      //printf("\n");

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationTdndmu(int N, int statistics, double T, double mu, double m, double deg,
                                             const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (N < 0 || N>3) {
        printf("**WARNING** QuantumNumericalIntegrationTdndmu: N must be between 0 and 3!\n");
        calculationHadBECIssue = true;
        exit(1);
      }
      if (N == 0)
        return QuantumNumericalIntegrationDensity(statistics, T, mu, m, deg, extraConfig);

      if (N == 1)
        return QuantumNumericalIntegrationT1dn1dmu1(statistics, T, mu, m, deg, extraConfig);

      if (N == 2)
        return QuantumNumericalIntegrationT2dn2dmu2(statistics, T, mu, m, deg, extraConfig);

      return QuantumNumericalIntegrationT3dn3dmu3(statistics, T, mu, m, deg, extraConfig);
    }

    double QuantumNumericalIntegrationChiN(int N, int statistics, double T, double mu, double m, double deg,
                                           const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return QuantumNumericalIntegrationTdndmu(N - 1, statistics, T, mu, m, deg, extraConfig) / pow(T, 3) / xMath::GeVtoifm3();
    }

    double QuantumNumericalIntegrationChiNDimensionfull(int N, int statistics, double T, double mu, double m, double deg,
                                                        const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 1 && T == 0.0)
        return FermiZeroTChiNDimensionfull(N, mu, m, deg, extraConfig);
      if (statistics == -1 && T == 0.0) {
        if (mu >= m) {
          printf("**WARNING** QuantumNumericalIntegrationChiNDimensionfull: Bose-Einstein condensation\n");
          calculationHadBECIssue = true;
        }
        return 0.;
      }
      return QuantumNumericalIntegrationTdndmu(N - 1, statistics, T, mu, m, deg, extraConfig) / pow(T, N-1) / xMath::GeVtoifm3();
    }

    double psi(double x)
    {
      if (x == 0.0)
        return 1.;
      double x2 = x * x;
      double tmpsqrt = sqrt(1. + x2);
      return (1. + x2 / 2.) * tmpsqrt - x2 * x2 / 2. * log((1. + tmpsqrt) / x);
    }

    double psi2(double x)
    {
      if (x == 0.0)
        return 2.;
      double x2 = x * x;
      double tmpsqrt = sqrt(1. + x2);
      return 2. * tmpsqrt + 2. * x2 * log((1. + tmpsqrt) / x);
    }

    double FermiNumericalIntegrationLargeMuDensity(double T, double mu, double m, double deg,
                                                   const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationDensity(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf / (exp(-(sqrt(legx32[i] * legx32[i] * pf * pf + m * m) - mu) / T) + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        ret1 += lagw32[i] * T * tx * T * tx * T / (exp(sqrt(tx*tx + moverT * moverT) - muoverT) + 1.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      double ret2 = 0.;
      ret2 += deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() * pf * pf * pf / 3.;

      return ret1 + ret2;
    }

    double FermiNumericalIntegrationLargeMuPressure(double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationPressure(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double x2 = legx32[i] * pf * legx32[i] * pf;
        double E = sqrt(legx32[i] * legx32[i] * pf*pf + m * m);
        ret1 += -legw32[i] * pf * x2 * x2 / E / (exp(-(E - mu) / T) + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double x2 = tx * T * tx * T;
        double E = sqrt(tx*tx + moverT * moverT);
        ret1 += lagw32[i] * T * x2 * x2 / E / T / (exp(E - muoverT) + 1.);
      }

      ret1 *= deg / 6. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      double ret2 = 0.;
      ret2 += mu * pf * pf * pf;
      ret2 += -3. / 4. * pf * pf * pf * pf * psi(m / pf);
      ret2 *= deg / 6. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret1 + ret2;
    }

    double FermiNumericalIntegrationLargeMuEnergyDensity(double T, double mu, double m, double deg,
                                                         const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationEnergyDensity(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * sqrt(legx32[i] * legx32[i] * pf*pf + m * m) / (exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m * m) - mu) / T) + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        ret1 += lagw32[i] * T * tx * T * tx * T * sqrt(tx*tx + moverT * moverT) * T / (exp(sqrt(tx*tx + moverT * moverT) - muoverT) + 1.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      double ret2 = 0.;
      ret2 += deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() * pf * pf * pf * pf / 4. * psi(m / pf);

      return ret1 + ret2;
    }

    double FermiNumericalIntegrationLargeMuEntropyDensity(double T, double mu, double m, double deg,
                                                          const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return (FermiNumericalIntegrationLargeMuPressure(T, mu, m, deg, extraConfig)
      + FermiNumericalIntegrationLargeMuEnergyDensity(T, mu, m, deg, extraConfig)
      - mu * FermiNumericalIntegrationLargeMuDensity(T, mu, m, deg, extraConfig)) / T;
    }

    double FermiNumericalIntegrationLargeMuScalarDensity(double T, double mu, double m, double deg,
                                                         const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationScalarDensity(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * m / sqrt(legx32[i] * legx32[i] * pf*pf + m * m) / (exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m * m) - mu) / T) + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        ret1 += lagw32[i] * T * tx * T * tx * T * moverT / sqrt(tx*tx + moverT * moverT) / (exp(sqrt(tx*tx + moverT * moverT) - muoverT) + 1.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      double ret2 = 0.;
      ret2 += deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() * m * pf * (mu - pf / 4. * psi2(m / pf));

      return ret1 + ret2;
    }

    double FermiNumericalIntegrationLargeMuT1dn1dmu1(double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationT1dn1dmu1(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double Eexp = exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m * m) - mu) / T);
        ret1 += legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * 1. / (1. + 1./Eexp) / (Eexp + 1.);
        //ret1 += legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * Eexp / (Eexp + 1.) / (Eexp + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double Eexp = exp(sqrt(tx*tx + moverT * moverT) - muoverT);
        ret1 += lagw32[i] * T * tx * T * tx * T * 1. / (1. + 1. / Eexp) / (Eexp + 1.);
        //ret1 += lagw32[i] * T * tx * T * tx * T * Eexp / (Eexp + 1.) / (Eexp + 1.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret1;
    }

    double FermiNumericalIntegrationLargeMuT2dn2dmu2(double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationT2dn2dmu2(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double Eexp = exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m * m) - mu) / T);
        ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * (1. - 1. / Eexp) / (1. + 1. / Eexp) / (1. + 1. / Eexp) / (Eexp + 1.);
        //ret1 += -legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * (Eexp*Eexp - Eexp) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double Eexp = exp(sqrt(tx*tx + moverT * moverT) - muoverT);
        ret1 += lagw32[i] * T * tx * T * tx * T * (1. - 1. / Eexp) / (1. + 1. / Eexp) / (1. + 1. / Eexp) / (Eexp + 1.);
        //ret1 += lagw32[i] * T * tx * T * tx * T * (Eexp*Eexp - Eexp) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret1;
    }

    double FermiNumericalIntegrationLargeMuT3dn3dmu3(double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationT3dn3dmu3(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double Eexp = exp(-(sqrt(legx32[i] * legx32[i] * pf*pf + m * m) - mu) / T);
        ret1 += legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * (1. - 4./Eexp + 1./Eexp / Eexp) / (1. + 1. / Eexp) / (1. + 1. / Eexp) / (1. + 1. / Eexp) / (Eexp + 1.);
        //ret1 += legw32[i] * pf * legx32[i] * pf * legx32[i] * pf * (Eexp*Eexp*Eexp - 4.*Eexp*Eexp + Eexp) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double Eexp = exp(sqrt(tx*tx + moverT * moverT) - muoverT);
        ret1 += lagw32[i] * T * tx * T * tx * T * (1. - 4. / Eexp + 1. / Eexp / Eexp) / (1. + 1. / Eexp) / (1. + 1. / Eexp) / (1. + 1. / Eexp) / (Eexp + 1.);
        //ret1 += lagw32[i] * T * tx * T * tx * T * (Eexp*Eexp*Eexp - 4.*Eexp*Eexp + Eexp) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.) / (Eexp + 1.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret1;
    }

    double FermiNumericalIntegrationLargeMuTdndmu(int N, double T, double mu, double m, double deg,
                                                  const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (N < 0 || N>3) {
        printf("**ERROR** FermiNumericalIntegrationLargeMuTdndmu: N < 0 or N > 3\n");
        exit(1);
      }
      if (N == 0)
        return FermiNumericalIntegrationLargeMuDensity(T, mu, m, deg, extraConfig);

      if (N == 1)
        return FermiNumericalIntegrationLargeMuT1dn1dmu1(T, mu, m, deg, extraConfig);

      if (N == 2)
        return FermiNumericalIntegrationLargeMuT2dn2dmu2(T, mu, m, deg, extraConfig);

      return FermiNumericalIntegrationLargeMuT3dn3dmu3(T, mu, m, deg, extraConfig);
    }

    double FermiNumericalIntegrationLargeMuChiN(int N, double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return FermiNumericalIntegrationLargeMuTdndmu(N - 1, T, mu, m, deg, extraConfig) / pow(T, 3) / xMath::GeVtoifm3();
    }

    double FermiNumericalIntegrationLargeMuChiNDimensionfull(int N, double T, double mu, double m, double deg,
                                                             const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return FermiNumericalIntegrationLargeMuTdndmu(N - 1, T, mu, m, deg, extraConfig) / pow(T, N-1) / xMath::GeVtoifm3();
    }

    double FermiZeroTDensity(double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (m >= mu)
        return 0.0;
      double pf = sqrt(mu * mu - m * m);
      return deg / 6. / xMath::Pi() / xMath::Pi() * pf * pf * pf * xMath::GeVtoifm3();
    }

    double FermiZeroTPressure(double mu, double m, double deg,
                              const IdealGasFunctionsExtraConfig& extraConfig)
    {
      assert(extraConfig.MagneticField.B == 0.0);
      if (m >= mu)
        return 0.0;
      double pf = sqrt(mu * mu - m * m);
      if (m == 0.0) {
        return deg / 24. / xMath::Pi() / xMath::Pi() * pf * pf * pf * pf * xMath::GeVtoifm3();
      }
      double m2 = m * m;
      return deg / 48. / xMath::Pi() / xMath::Pi() *
        (
          mu * pf * (2. * mu * mu - 5. * m2) 
          - 3. * m2 * m2 * log(m / (mu + pf))
          ) * xMath::GeVtoifm3();
    }

    double FermiZeroTEnergyDensity(double mu, double m, double deg,
                                   const IdealGasFunctionsExtraConfig& extraConfig)
    {
      assert(extraConfig.MagneticField.B == 0.0);
      if (m >= mu)
        return 0.0;
      double pf = sqrt(mu * mu - m * m);
      if (m == 0.0) {
        return deg / 8. / xMath::Pi() / xMath::Pi() * pf * pf * pf * pf * xMath::GeVtoifm3();
      }
      double m2 = m * m;
      return deg / 16. / xMath::Pi() / xMath::Pi() *
        (
          mu * pf * (2. * mu * mu - m2) 
          + m2 * m2 * log(m / (mu + pf))
          ) * xMath::GeVtoifm3();
    }

    double FermiZeroTEntropyDensity(double mu, double m, double deg,
                                    const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return 0.0;
    }

    double FermiZeroTScalarDensity(double mu, double m, double deg,
                                   const IdealGasFunctionsExtraConfig& extraConfig)
    {
      assert(extraConfig.MagneticField.B == 0.0);
      if (m >= mu)
        return 0.0;
      double pf = sqrt(mu * mu - m * m);
      if (m == 0.0) {
        return 0.;
      }
      double m2 = m * m;
      return deg * m / 4. / xMath::Pi() / xMath::Pi() *
        (
          mu * pf 
          + m2 * log(m / (mu + pf))
          ) * xMath::GeVtoifm3();
    }

    double FermiZeroTdn1dmu1(double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig)
    {
      assert(extraConfig.MagneticField.B == 0.0);
      if (m >= mu)
        return 0.0;
      double pf = sqrt(mu * mu - m * m);
      return deg / 2. / xMath::Pi() / xMath::Pi() * mu * pf * xMath::GeVtoifm3();
    }

    double FermiZeroTdn2dmu2(double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig)
    {
      assert(extraConfig.MagneticField.B == 0.0);
      if (m >= mu)
        return 0.0;
      double pf = sqrt(mu * mu - m * m);
      return deg / 2. / xMath::Pi() / xMath::Pi() * (mu * mu + pf * pf) / pf * xMath::GeVtoifm3();
    }

    double FermiZeroTdn3dmu3(double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig)
    {
      assert(extraConfig.MagneticField.B == 0.0);
      if (m >= mu)
        return 0.0;
      double pf = sqrt(mu * mu - m * m);
      return deg / 2. / xMath::Pi() / xMath::Pi() * mu * (3. * pf * pf - mu * mu) / pf / pf / pf * xMath::GeVtoifm3();
    }

    double FermiZeroTdndmu(int N, double mu, double m, double deg,
                           const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (N < 0 || N>3) {
        printf("**ERROR** FermiNumericalIntegrationLargeMuTdndmu: N < 0 or N > 3\n");
        exit(1);
      }
      if (N == 0)
        return FermiZeroTDensity(mu, m, deg, extraConfig);

      if (N == 1)
        return FermiZeroTdn1dmu1(mu, m, deg, extraConfig);

      if (N == 2)
        return FermiZeroTdn2dmu2(mu, m, deg, extraConfig);

      return FermiZeroTdn3dmu3(mu, m, deg, extraConfig);
    }

    double FermiZeroTChiN(int N, double mu, double m, double deg,
                          const IdealGasFunctionsExtraConfig& extraConfig)
    {
      printf("**ERROR** FermiZeroTChiN: This quantity is infinite by definition at T = 0!\n");
      exit(1);
      return 0.0;
      //return FermiNumericalIntegrationLargeMuTdndmu(N - 1, mu, m, deg) / pow(T, 3) / xMath::GeVtoifm3();
    }

    double FermiZeroTChiNDimensionfull(int N, double mu, double m, double deg,
                                       const IdealGasFunctionsExtraConfig& extraConfig)
    {
      return FermiZeroTdndmu(N - 1, mu, m, deg, extraConfig) / xMath::GeVtoifm3();
    }

    double IdealGasQuantity(Quantity quantity, QStatsCalculationType calctype, int statistics, double T, double mu, double m, double deg, int order,
                            const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0) {
        if (quantity == ParticleDensity)
          return BoltzmannDensity(T, mu, m, deg, extraConfig);
        if (quantity == Pressure)
          return BoltzmannPressure(T, mu, m, deg, extraConfig);
        if (quantity == EnergyDensity)
          return BoltzmannEnergyDensity(T, mu, m, deg, extraConfig);
        if (quantity == EntropyDensity)
          return BoltzmannEntropyDensity(T, mu, m, deg, extraConfig);
        if (quantity == ScalarDensity)
          return BoltzmannScalarDensity(T, mu, m, deg, extraConfig);
        if (quantity == chi2)
          return BoltzmannChiN(2, T, mu, m, deg, extraConfig);
        if (quantity == chi3)
          return BoltzmannChiN(3, T, mu, m, deg, extraConfig);
        if (quantity == chi4)
          return BoltzmannChiN(4, T, mu, m, deg, extraConfig);
        if (quantity == chi2difull)
          return BoltzmannChiNDimensionfull(2, T, mu, m, deg, extraConfig);
        if (quantity == chi3difull)
          return BoltzmannChiNDimensionfull(3, T, mu, m, deg, extraConfig);
        if (quantity == chi4difull)
          return BoltzmannChiNDimensionfull(4, T, mu, m, deg, extraConfig);
      }
      else {
        if (calctype == ClusterExpansion) {
          if (quantity == ParticleDensity)
            return QuantumClusterExpansionDensity(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == Pressure)
            return QuantumClusterExpansionPressure(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == EnergyDensity)
            return QuantumClusterExpansionEnergyDensity(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == EntropyDensity)
            return QuantumClusterExpansionEntropyDensity(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == ScalarDensity)
            return QuantumClusterExpansionScalarDensity(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == chi2)
            return QuantumClusterExpansionChiN(2, statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == chi3)
            return QuantumClusterExpansionChiN(3, statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == chi4)
            return QuantumClusterExpansionChiN(4, statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == chi2difull)
            return QuantumClusterExpansionChiNDimensionfull(2, statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == chi3difull)
            return QuantumClusterExpansionChiNDimensionfull(3, statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == chi4difull)
            return QuantumClusterExpansionChiNDimensionfull(4, statistics, T, mu, m, deg, order, extraConfig);
        }
        else {
          if (quantity == ParticleDensity)
            return QuantumNumericalIntegrationDensity(statistics, T, mu, m, deg, extraConfig);
          if (quantity == Pressure)
            return QuantumNumericalIntegrationPressure(statistics, T, mu, m, deg, extraConfig);
          if (quantity == EnergyDensity)
            return QuantumNumericalIntegrationEnergyDensity(statistics, T, mu, m, deg, extraConfig);
          if (quantity == EntropyDensity)
            return QuantumNumericalIntegrationEntropyDensity(statistics, T, mu, m, deg, extraConfig);
          if (quantity == ScalarDensity)
            return QuantumNumericalIntegrationScalarDensity(statistics, T, mu, m, deg, extraConfig);
          if (quantity == chi2)
            return QuantumNumericalIntegrationChiN(2, statistics, T, mu, m, deg, extraConfig);
          if (quantity == chi3)
            return QuantumNumericalIntegrationChiN(3, statistics, T, mu, m, deg, extraConfig);
          if (quantity == chi4)
            return QuantumNumericalIntegrationChiN(4, statistics, T, mu, m, deg, extraConfig);
          if (quantity == chi2difull)
            return QuantumNumericalIntegrationChiNDimensionfull(2, statistics, T, mu, m, deg, extraConfig);
          if (quantity == chi3difull)
            return QuantumNumericalIntegrationChiNDimensionfull(3, statistics, T, mu, m, deg, extraConfig);
          if (quantity == chi4difull)
            return QuantumNumericalIntegrationChiNDimensionfull(4, statistics, T, mu, m, deg, extraConfig);
        }
      }
      printf("**WARNING** IdealGasFunctions::IdealGasQuantity: Unknown quantity\n");
      return 0.;
    }


  } // namespace IdealGasFunctions

} // namespace thermalfist
