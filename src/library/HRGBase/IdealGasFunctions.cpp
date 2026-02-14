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
#include <cassert>

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
        // Use e = -p + T (dp/dT) + mu (dp/mu) + B * (dp/dB)
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
        return ret + extraConfig.MagneticField.B * BoltzmannMagnetization(T, mu, m, deg, extraConfig) * xMath::GeVtoifm3();
      }

      // No magnetic field
      if (m == 0.)
        return 3 * T * BoltzmannDensity(T, mu, m, deg);
      return (3 * T + m * xMath::BesselK1exp(m / T) / xMath::BesselKexp(2, m / T)) * BoltzmannDensity(T, mu, m, deg);
    }

    double BoltzmannEntropyDensity(double T, double mu, double m, double deg,
                                   const IdealGasFunctionsExtraConfig& extraConfig) {
      double ret = (BoltzmannPressure(T, mu, m, deg, extraConfig) +
                    BoltzmannEnergyDensity(T, mu, m, deg, extraConfig) -
                    mu * BoltzmannDensity(T, mu, m, deg, extraConfig)) / T;

      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.)
        ret -= extraConfig.MagneticField.B * BoltzmannMagnetization(T, mu, m, deg, extraConfig) * xMath::GeVtoifm3() / T;

      return ret;
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
        if (m == 0.)
          ret += sign * cfug / static_cast<double>(i) / static_cast<double>(i) / static_cast<double>(i);
        else
          ret += sign * xMath::BesselKexp(2, i*moverT) * cfug / static_cast<double>(i);
        cfug *= tfug;
        if (signchange) sign = -sign;
      }

      double prefactor = 1.;
      if (m == 0.)
        prefactor = deg * T * T * T / 2. / xMath::Pi() / xMath::Pi() * 2.;
      else
        prefactor = deg * m * m * T / 2. / xMath::Pi() / xMath::Pi();

      return ret * prefactor * xMath::GeVtoifm3();
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
        if (m == 0.)
          ret += sign * cfug / static_cast<double>(i) / static_cast<double>(i) / static_cast<double>(i) / static_cast<double>(i);
        else
          ret += sign * xMath::BesselKexp(2, i*moverT) * cfug / static_cast<double>(i) / static_cast<double>(i);
        cfug *= tfug;
        if (signchange) sign = -sign;
      }

      double prefactor = 1.;
      if (m == 0.)
        prefactor = deg * T * T * T * T / 2. / xMath::Pi() / xMath::Pi() * 2.;
      else
        prefactor = deg * m * m * T * T / 2. / xMath::Pi() / xMath::Pi();

      return ret * prefactor * xMath::GeVtoifm3();
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
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3()
               + extraConfig.MagneticField.B * QuantumClusterExpansionMagnetization(statistics, T, mu, m, deg, order, extraConfig) * xMath::GeVtoifm3();
      }

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double cfug = tfug;
      double moverT = m / T;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        if (m == 0.)
          ret += sign * cfug / static_cast<double>(i) / static_cast<double>(i) / static_cast<double>(i) / static_cast<double>(i);
        else
          ret += sign * (xMath::BesselKexp(1, i*moverT) + 3. * xMath::BesselKexp(2, i*moverT) / moverT / static_cast<double>(i)) * cfug / static_cast<double>(i);
        cfug *= tfug;
        if (signchange) sign = -sign;
      }

      double prefactor = 1.;
      if (m == 0.)
        prefactor = deg * T * T * T * T / 2. / xMath::Pi() / xMath::Pi() * 6.;
      else
        prefactor = deg * m * m * m * T / 2. / xMath::Pi() / xMath::Pi();

      return ret * prefactor * xMath::GeVtoifm3();
    }

    double QuantumClusterExpansionEntropyDensity(int statistics, double T, double mu, double m, double deg, int order,
                                                 const IdealGasFunctionsExtraConfig& extraConfig)
    {
      double ret = (QuantumClusterExpansionPressure(statistics, T, mu, m, deg, order, extraConfig) +
                    QuantumClusterExpansionEnergyDensity(statistics, T, mu, m, deg, order, extraConfig) -
                    mu * QuantumClusterExpansionDensity(statistics, T, mu, m, deg, order, extraConfig)) / T;

      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.)
        ret -= extraConfig.MagneticField.B * QuantumClusterExpansionMagnetization(statistics, T, mu, m, deg, order, extraConfig) * xMath::GeVtoifm3() / T;

      return ret;
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
      
      if (m == 0.) 
        return ret;
      
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
        if (m == 0.)
          ret += sign * cfug * pow(static_cast<double>(i), N - 3);
        else
          ret += sign * xMath::BesselKexp(2, i*moverT) * cfug * pow(static_cast<double>(i), N - 1);
        cfug *= tfug;
        if (signchange) sign = -sign;
      } 

      double prefactor = 1.;
      if (m == 0.)
        prefactor = deg * T * T * T / 2. / xMath::Pi() / xMath::Pi() * 2.;
      else
        prefactor = deg * m * m * T / 2. / xMath::Pi() / xMath::Pi();

      return ret * prefactor * xMath::GeVtoifm3();
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
        std::cerr << "**WARNING** QuantumNumericalIntegrationDensity: Bose-Einstein condensation, mass = " << m << ", mu = " << mu << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
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
        std::cerr << "**WARNING** QuantumNumericalIntegrationPressure: Bose-Einstein condensation" << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
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
        std::cerr << "**WARNING** QuantumNumericalIntegrationEnergyDensity: Bose-Einstein condensation" << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
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
        return ret * Qmod * extraConfig.MagneticField.B / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3()
         + extraConfig.MagneticField.B * QuantumNumericalIntegrationMagnetization(statistics, T, mu, m, deg, extraConfig) * xMath::GeVtoifm3();
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

      if (statistics == -1 && mu > m) {
        std::cerr << "**WARNING** QuantumNumericalIntegrationEntropyDensity: Bose-Einstein condensation" << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
      }

      // Use direct entropy integrand for degenerate Fermi gas (avoids (P+e-mu*n)/T cancellation)
      if (statistics == 1 && mu > m)
        return FermiNumericalIntegrationLargeMuEntropyDensity(T, mu, m, deg, extraConfig);

      double ret = (QuantumNumericalIntegrationPressure(statistics, T, mu, m, deg, extraConfig)
                    + QuantumNumericalIntegrationEnergyDensity(statistics, T, mu, m, deg, extraConfig)
                    - mu * QuantumNumericalIntegrationDensity(statistics, T, mu, m, deg, extraConfig)) / T;

      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.)
        ret -= extraConfig.MagneticField.B * QuantumNumericalIntegrationMagnetization(statistics, T, mu, m, deg, extraConfig) * xMath::GeVtoifm3() / T;

      return ret;
    }

    double QuantumNumericalIntegrationScalarDensity(int statistics, double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (statistics == 0)           return BoltzmannScalarDensity(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return FermiZeroTScalarDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuScalarDensity(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        std::cerr << "**WARNING** QuantumNumericalIntegrationScalarDensity: Bose-Einstein condensation" << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
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
        std::cerr << "**WARNING** QuantumNumericalIntegrationT1dn1dmu1: Bose-Einstein condensation" << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
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
        std::cerr << "**WARNING** QuantumNumericalIntegrationT2dn2dmu2: Bose-Einstein condensation" << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
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
        std::cerr << "**WARNING** QuantumNumericalIntegrationT3dn3dmu3: Bose-Einstein condensation" << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
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
        std::cerr << "**WARNING** QuantumNumericalIntegrationTdndmu: N must be between 0 and 3!" << std::endl;
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
          std::cerr << "**WARNING** QuantumNumericalIntegrationChiNDimensionfull: Bose-Einstein condensation" << std::endl;
          calculationHadBECIssue = true;
          mu = m;
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

    /// \brief Sommerfeld-inspired variable substitution for the Legendre part
    ///        of the large-mu Fermi integrals.
    ///
    /// Maps s in [0,1] (Legendre node) to momentum p in [0, pF] such that
    /// nodes cluster near pF at low T (resolving the Fermi surface)
    /// while reducing to the uniform mapping p = s * pF at high T.
    ///
    /// Substitution:  u = -log(1 - beta * (1 - s)),   p = pF * (1 - u / alpha)
    /// where alpha = pF^2 / (mu * T),  beta = 1 - exp(-alpha).
    ///
    /// Returns p and dp/ds.
    inline void SommerfeldLegendreMap(double s, double pf, double mu, double T,
                                      double& p, double& dpds)
    {
      double alpha = pf * pf / (mu * T);
      if (alpha < 1.e-6) {
        // Small alpha: uniform mapping (avoid catastrophic cancellation in beta)
        p = s * pf;
        dpds = pf;
        return;
      }
      double expmAlpha = exp(-alpha);
      double beta = 1. - expmAlpha;
      double t = 1. - s;
      double oneMinusBetaT = 1. - beta * t;
      double u = -log(oneMinusBetaT);
      p = pf * (1. - u / alpha);
      dpds = pf / alpha * beta / oneMinusBetaT;
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
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        ret1 += -legw32[i] * dpds * p * p / (exp(-(sqrt(p * p + m * m) - mu) / T) + 1.);
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

      // Other terms from differentiating Eq. (55) in https://arxiv.org/pdf/1901.05249.pdf add up to zero

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
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double p2 = p * p;
        double E = sqrt(p2 + m * m);
        ret1 += -legw32[i] * dpds * p2 * p2 / E / (exp(-(E - mu) / T) + 1.);
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
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double E = sqrt(p * p + m * m);
        ret1 += -legw32[i] * dpds * p * p * E / (exp(-(E - mu) / T) + 1.);
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

    /// \brief Entropy integrand per mode: sigma(x) = -f*ln(f) - (1-f)*ln(1-f)
    ///        where f = 1/(e^x + 1) and x = (E - mu)/T.
    ///
    /// Numerically stable form: sigma = |x| / (e^{|x|} + 1) + ln(1 + e^{-|x|})
    inline double FermiEntropySigma(double x)
    {
      double ax = std::abs(x);
      return log(1. + exp(-ax)) + ax / (exp(ax) + 1.);
    }

    double FermiNumericalIntegrationLargeMuEntropyDensity(double T, double mu, double m, double deg,
                                                          const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationEntropyDensity(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu * mu - m * m);
      double ret1 = 0.;

      // Legendre on [0, pF] with Sommerfeld mapping: sigma is non-zero only near pF
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double E = sqrt(p * p + m * m);
        double x = (E - mu) / T;
        ret1 += legw32[i] * dpds * p * p * FermiEntropySigma(x);
      }

      // Shifted Laguerre above pF
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double E = sqrt(tx * tx + moverT * moverT);
        double x = E - muoverT;
        ret1 += lagw32[i] * T * tx * T * tx * T * FermiEntropySigma(x);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      // No ret2: at T=0 sigma=0 everywhere, so the analytic T=0 part is exactly zero
      return ret1;
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
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double E = sqrt(p * p + m * m);
        ret1 += -legw32[i] * dpds * p * p * m / E / (exp(-(E - mu) / T) + 1.);
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

      // Integration-by-parts (IBP) decomposition:
      //   int p^2 f(1-f) dp = T * int (2E - m^2/E) f dp
      //     (via df/dp = -f(1-f)/T * p/E)
      //   = T * [ int_0^pF (2E - m^2/E) dp  +  thermal correction ]
      //   = T * [ mu*pF + O(T^2) ]
      //
      // The analytic part mu*pF is extracted exactly, ensuring the correct
      // T->0 limit:  T * dn/dmu -> g/(2pi^2) * mu * pF  (no precision loss).
      // The thermal correction (ret1) vanishes as T->0.
      //
      // The Sommerfeld-Legendre quadrature splits the correction into:
      //   ret1 = -int_0^pF (2E - m^2/E)(1-f) dp  +  int_pF^inf (2E - m^2/E) f dp

      double m2 = m * m;
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double E = sqrt(p * p + m2);
        ret1 += -legw32[i] * dpds * (2. * E - m2 / E) / (exp(-(E - mu) / T) + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      double moverT2 = moverT * moverT;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double EoverT = sqrt(tx * tx + moverT2);
        ret1 += lagw32[i] * T * T * (2. * tx * tx + moverT2) / EoverT / (exp(EoverT - muoverT) + 1.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      double ret2 = deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() * mu * pf;

      return T * (ret1 + ret2);
    }

    double FermiNumericalIntegrationLargeMuT2dn2dmu2(double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationT2dn2dmu2(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);

      // Integration-by-parts (IBP) decomposition:
      //   This function returns T^2 * d^2n/dmu^2
      //     = g/(2pi^2) * int p^2 f(1-f)(1-2f) dp
      //     = g/(2pi^2) * T * int (2p^2+m^2)/E f(1-f) dp
      //
      // via f(1-f)(1-2f) = -T*(E/p) * d/dp[f(1-f)], integrate by parts,
      // boundary terms vanish.
      //
      // The integrand (2p^2+m^2)/E * f(1-f) is a smooth peaked function
      // ~1/(4 cosh^2((E-mu)/(2T))) with no explicit 1/T factor.
      // The Sommerfeld-Legendre + Laguerre quadrature concentrates points
      // near the Fermi surface.
      //
      // The T->0 limit of d^2n/dmu^2 is (mu^2+pF^2)/pF, recovered via
      //   d^2n/dmu^2 = (1/T) * int (2p^2+m^2)/E f(1-f) dp
      //              ~ (1/T) * T*(mu^2+pF^2)/pF = (mu^2+pF^2)/pF
      //
      // Analytic T=0 term is extracted exactly (ret2), thermal correction
      // computed via quadrature (ret1).  Since T^2 * (ret2 + ret1) / T^2
      // round-trips without precision loss in IEEE double, this ensures
      // accuracy even when the caller divides by T^2 to get d^2n/dmu^2.

      double m2 = m * m;
      double quad = 0.;
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double E = sqrt(p * p + m2);
        double x = (E - mu) / T;
        // f(1-f) = exp(x)/(1+exp(x))^2, stable for x < 0 (below pF)
        double ex = exp(x);
        double ff1mf = ex / ((1. + ex) * (1. + ex));
        quad += legw32[i] * dpds * (2. * p * p + m2) / E * ff1mf;
      }

      double moverT = m / T;
      double muoverT = mu / T;
      double moverT2 = moverT * moverT;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double EoverT = sqrt(tx * tx + moverT2);
        double x = EoverT - muoverT;
        // f(1-f) = exp(-x)/(1+exp(-x))^2, stable for x > 0 (above pF)
        double emx = exp(-x);
        double ff1mf = emx / ((1. + emx) * (1. + emx));
        quad += lagw32[i] * T * T * (2. * tx * tx + moverT2) / EoverT * ff1mf;
      }

      // quad = int (2p^2+m^2)/E f(1-f) dp  ~  T*(mu^2+pF^2)/pF as T->0
      // d^2n/dmu^2 = quad/T = (mu^2+pF^2)/pF + O(T^2)
      double ret2 = (mu * mu + pf * pf) / pf;   // analytic T=0 value of d^2n/dmu^2
      double ret1 = quad / T - ret2;             // thermal correction, -> 0 as T->0

      return T * T * (ret2 + ret1) * deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
    }

    double FermiNumericalIntegrationLargeMuT3dn3dmu3(double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (mu <= m)
        return QuantumNumericalIntegrationT3dn3dmu3(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);

      // Integration-by-parts (IBP) decomposition:
      //   This function returns T^3 * d^3n/dmu^3
      //     = g/(2pi^2) * int p^2 f(1-f)(1-6f(1-f)) dp
      //     = g/(2pi^2) * T * int (2p^2+m^2)/E f(1-f)(1-2f) dp
      //
      // via f(1-f)(1-6f(1-f)) = -T*(E/p) * d/dp[f(1-f)(1-2f)],
      // integrate by parts, boundary terms vanish.
      //
      // The T->0 value of d^3n/dmu^3 is mu*(3*pF^2-mu^2)/pF^3,
      // extracted analytically (ret2).  Thermal correction (ret1)
      // computed via quadrature, -> 0 as T->0.

      double m2 = m * m;
      double quad = 0.;
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double E = sqrt(p * p + m2);
        double x = (E - mu) / T;
        // f(1-f)(1-2f): use exp(x) form, stable for x < 0 (below pF)
        double ex = exp(x);
        double f = 1. / (ex + 1.);
        double ff1mf = ex / ((1. + ex) * (1. + ex));
        double ff1mf_1m2f = ff1mf * (1. - 2. * f);
        quad += legw32[i] * dpds * (2. * p * p + m2) / E * ff1mf_1m2f;
      }

      double moverT = m / T;
      double muoverT = mu / T;
      double moverT2 = moverT * moverT;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double EoverT = sqrt(tx * tx + moverT2);
        double x = EoverT - muoverT;
        // Above pF: x > 0, use exp(-x) forms for stability
        double emx = exp(-x);
        double f = emx / (1. + emx);
        double ff1mf = emx / ((1. + emx) * (1. + emx));
        double ff1mf_1m2f = ff1mf * (1. - 2. * f);
        quad += lagw32[i] * T * T * (2. * tx * tx + moverT2) / EoverT * ff1mf_1m2f;
      }

      // quad = int (2p^2+m^2)/E f(1-f)(1-2f) dp  ~  T^2 * mu*(3pF^2-mu^2)/pF^3 as T->0
      // d^3n/dmu^3 = quad/T^2 = mu*(3pF^2-mu^2)/pF^3 + O(T^2)
      double ret2 = mu * (3. * pf * pf - mu * mu) / (pf * pf * pf);  // analytic T=0
      double ret1 = quad / (T * T) - ret2;                            // thermal correction

      return T * T * T * (ret2 + ret1) * deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();
    }

    double FermiNumericalIntegrationLargeMuTdndmu(int N, double T, double mu, double m, double deg,
                                                  const IdealGasFunctionsExtraConfig& extraConfig)
    {
      if (N < 0 || N>3) {
        throw std::runtime_error("**ERROR** FermiNumericalIntegrationLargeMuTdndmu: N < 0 or N > 3");
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
        throw std::runtime_error("**ERROR** FermiNumericalIntegrationLargeMuTdndmu: N < 0 or N > 3");
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
      throw std::runtime_error("**ERROR** FermiZeroTChiN: This quantity is infinite by definition at T = 0!");
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
        // Temperature derivatives
        if (quantity == dndT)
          return BoltzmanndndT(T, mu, m, deg, extraConfig);
        if (quantity == d2ndT2)
          return Boltzmannd2ndT2(T, mu, m, deg, extraConfig);
        if (quantity == dedT)
          return BoltzmanndedT(T, mu, m, deg, extraConfig);
        if (quantity == dedmu)
          return Boltzmanndedmu(T, mu, m, deg, extraConfig);
        if (quantity == dchi2dT)
          return Boltzmanndchi2dT(T, mu, m, deg, extraConfig);
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
          // Temperature derivatives
          if (quantity == dndT)
            return QuantumClusterExpansiondndT(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == d2ndT2)
            return QuantumClusterExpansiond2ndT2(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == dedT)
            return QuantumClusterExpansiondedT(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == dedmu)
            return QuantumClusterExpansiondedmu(statistics, T, mu, m, deg, order, extraConfig);
          if (quantity == dchi2dT)
            return QuantumClusterExpansiondchi2dT(statistics, T, mu, m, deg, order, extraConfig);
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
          // Temperature derivatives
          if (quantity == dndT)
            return QuantumNumericalIntegrationdndT(statistics, T, mu, m, deg, extraConfig);
          if (quantity == d2ndT2)
            return QuantumNumericalIntegrationd2ndT2(statistics, T, mu, m, deg, extraConfig);
          if (quantity == dedT)
            return QuantumNumericalIntegrationdedT(statistics, T, mu, m, deg, extraConfig);
          if (quantity == dedmu)
            return QuantumNumericalIntegrationdedmu(statistics, T, mu, m, deg, extraConfig);
          if (quantity == dchi2dT)
            return QuantumNumericalIntegrationdchi2dT(statistics, T, mu, m, deg, extraConfig);
        }
      }
      std::cerr << "**WARNING** IdealGasFunctions::IdealGasQuantity: Unknown quantity" << std::endl;
      return 0.;
    }

    double BoltzmannMagnetization(double T, double mu, double m, double deg, const IdealGasFunctionsExtraConfig &extraConfig) {
      // Check for magnetic field effect for charged particles
      if (extraConfig.MagneticField.B != 0. && extraConfig.MagneticField.Q != 0.) {
        // TODO: Check that it's done correctly
        double Qmod = abs(extraConfig.MagneticField.Q);
        auto spins = GetSpins(extraConfig.MagneticField.degSpin);
        // Sum over spins
        double ret = 0.;
        for(double sz : spins) {
          // Sum over Landau levels
          for(int l = 0; l < extraConfig.MagneticField.lmax; ++l) {
            double e0 = sqrt(m*m + Qmod*extraConfig.MagneticField.B*(2.*l + 1. - 2.*sz));
            ret += (T * e0 * xMath::BesselKexp(1, e0 / T)
                    - Qmod * extraConfig.MagneticField.B * (l + 0.5 - sz) * xMath::BesselKexp(0, e0 / T)
                   ) * exp((mu - e0) / T);
          }
        }
        ret *= Qmod / 2. / xMath::Pi() / xMath::Pi();
        return ret;
      }

      // No magnetic field
      return 0.;
    }

    double QuantumClusterExpansionMagnetization(int statistics, double T, double mu, double m, double deg, int order,
                                                const IdealGasFunctionsExtraConfig &extraConfig) {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return BoltzmannMagnetization(T, mu, m, deg, extraConfig);

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
              ret += sign * (T / static_cast<double>(i) * e0 * xMath::BesselKexp(1, i*EoverT)
               - Qmod  * extraConfig.MagneticField.B * (l + 0.5 - sz) * xMath::BesselKexp(0, i*EoverT)
              ) * cfug;

              cfug *= tfug;
              if (signchange) sign = -sign;
            }
          }
        }
        return ret * Qmod / 2. / xMath::Pi() / xMath::Pi();
      }

      // No magnetic field
      return 0.;
    }

    double QuantumNumericalIntegrationMagnetization(int statistics, double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig &extraConfig) {
      if (statistics == 0)           return BoltzmannMagnetization(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return FermiZeroTMagnetization(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMuMagnetization(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        std::cerr << "**WARNING** QuantumNumericalIntegrationScalarDensity: Bose-Einstein condensation" << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
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
              ret += lagw32[i] * T * (tx * T * tx - Qmod * BoverT2 * T * (l + 0.5 - sz) ) / EoverT / (exp(EoverT - muoverT) + statistics);
            }
          }
        }
        return ret * Qmod / 2. / xMath::Pi() / xMath::Pi();
      }

      // No magnetic field
      return 0.;
    }

    double FermiNumericalIntegrationLargeMuMagnetization(double T, double mu, double m, double deg,
                                                         const IdealGasFunctionsExtraConfig &extraConfig) {
      if (mu <= m)
        return QuantumNumericalIntegrationMagnetization(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      return 0.;
    }

    double FermiZeroTMagnetization(double mu, double m, double deg, const IdealGasFunctionsExtraConfig &extraConfig) {
      assert(extraConfig.MagneticField.B == 0.0);
      return 0.;
    }

    /// Temperature derivatives

    double BoltzmanndndT(double T, double mu, double m, double deg, const IdealGasFunctionsExtraConfig &extraConfig) {
      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      if (m == 0.)
        return deg / xMath::Pi() / xMath::Pi() * T * exp(mu/ T) * (3. * T - mu) * xMath::GeVtoifm3();
      return deg * m * m / 2. / T / xMath::Pi() / xMath::Pi()
      * (m * xMath::BesselKexp(1, m / T) - (mu - 3. * T) * xMath::BesselKexp(2, m / T))
      * exp((mu - m) / T) * xMath::GeVtoifm3();
    }

    double Boltzmannd2ndT2(double T, double mu, double m, double deg, const IdealGasFunctionsExtraConfig &extraConfig) {
      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      if (m == 0.)
        return deg / xMath::Pi() / xMath::Pi() / T * exp(mu/ T) *
          (mu * mu - 4. * mu * T + 6. * T * T) * xMath::GeVtoifm3();
      return deg / 2. / T / T / T / xMath::Pi() / xMath::Pi()
             * m * (m * (m * m + mu * mu - 4. * mu * T + 6. * T * T) * xMath::BesselKexp(0, m / T)
             + (m * m * (3. * T - 2. * mu) + 2. * T * (mu * mu - 4. * mu * T + 6. * T * T)) * xMath::BesselKexp(1, m / T))
             * exp((mu - m) / T) * xMath::GeVtoifm3();
    }

    double BoltzmanndedT(double T, double mu, double m, double deg, const IdealGasFunctionsExtraConfig &extraConfig) {
      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      if (m == 0.)
        return 3. * deg / xMath::Pi() / xMath::Pi() * T * T * exp(mu/ T) *
               (4. * T - mu) * xMath::GeVtoifm3();

      return deg * m / 2. / T / xMath::Pi() / xMath::Pi()
      * (m * (m*m + 3*T*(4.*T-mu))*xMath::BesselKexp(0, m / T)
      + (6.*T*T*(4.*T-mu) - m*m*(mu-5.*T) )*xMath::BesselKexp(1, m / T))
      * exp((mu - m) / T) * xMath::GeVtoifm3();
    }

    double Boltzmanndedmu(double T, double mu, double m, double deg, const IdealGasFunctionsExtraConfig &extraConfig) {
      return BoltzmannEnergyDensity(T, mu, m, deg, extraConfig) / T;
    }

    double Boltzmanndchi2dT(double T, double mu, double m, double deg, const IdealGasFunctionsExtraConfig &extraConfig) {
      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      if (m == 0.)
        return -deg * mu / T / T / xMath::Pi() / xMath::Pi() * exp(mu/ T);
      return deg * m * m / 2. / T / T / T / T / xMath::Pi() / xMath::Pi()
        * (m * xMath::BesselKexp(1, m / T) - mu * xMath::BesselKexp(2, m / T))
        * exp((mu - m) / T);
    }

    double QuantumClusterExpansiondndT(int statistics, double T, double mu, double m, double deg, int order,
                                       const IdealGasFunctionsExtraConfig &extraConfig) {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return BoltzmanndndT(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double moverT = m / T;
      double cfug = tfug;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        if (m == 0.) {
          ret += sign * (
                - (mu - 3. * T/ static_cast<double>(i)) / static_cast<double>(i) / static_cast<double>(i)
                ) * cfug;
        }
        else {
          ret += sign * (
                m * xMath::BesselKexp(1, i*moverT)
                - (mu - 3. * T/ static_cast<double>(i)) * xMath::BesselKexp(2, i*moverT)
                ) * cfug;
        }
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      double prefactor = 1.;
      if (m == 0.) {
        prefactor = deg * T / xMath::Pi() / xMath::Pi();
      } else {
        prefactor = deg * m * m / T / 2. / xMath::Pi() / xMath::Pi();
      }
      return ret * prefactor * xMath::GeVtoifm3();
    }

    double QuantumClusterExpansiond2ndT2(int statistics, double T, double mu, double m, double deg, int order,
                                       const IdealGasFunctionsExtraConfig &extraConfig) {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return Boltzmannd2ndT2(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double moverT = m / T;
      double cfug = tfug;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        double k = static_cast<double>(i);
        if (m == 0.) {
          ret += sign * 2. * T * (mu * mu - 4. * mu * T / k + 6. * T / k * T / k) / k * cfug;
        } else {
          ret += sign * k * (
                (m * (m * m + mu * mu - 4. * mu * T / k + 6. * T * T / k / k) * xMath::BesselKexp(0, i*moverT)
                + (m * m * (3. * T / k - 2. * mu) + 2. * T / k * (mu * mu - 4. * mu * T / k + 6. * T / k * T / k)) * xMath::BesselKexp(1, i*moverT))
                ) * cfug;
        }
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      double prefactor = 1.;
      if (m == 0.) {
        prefactor = deg / T / T / 2. / xMath::Pi() / xMath::Pi();
      } else {
        prefactor = deg * m / T / T / T / 2. / xMath::Pi() / xMath::Pi();
      }
      return ret * prefactor * xMath::GeVtoifm3();
    }

    double QuantumClusterExpansiondedT(int statistics, double T, double mu, double m, double deg, int order,
                                       const IdealGasFunctionsExtraConfig &extraConfig) {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return BoltzmanndedT(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double moverT = m / T;
      double cfug = tfug;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        double k = static_cast<double>(i);
        if (m == 0.) {
          ret += sign * 6. * T / k * T / k * (4. * T / k - mu) / k * cfug;
        } else {
          ret += sign * (
                m * (m*m + 3*T/k*(4.*T/k-mu))*xMath::BesselKexp(0, i*moverT)
                + (6.*T/k*T/k*(4.*T/k-mu) - m*m*(mu-5.*T/k) )*xMath::BesselKexp(1, i*moverT)
                ) * cfug;
        }
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      double prefactor = 1.;
      if (m == 0.) {
        prefactor = deg / 2. / xMath::Pi() / xMath::Pi();
      } else {
        prefactor = deg * m / T / 2. / xMath::Pi() / xMath::Pi();
      }
      return ret * prefactor * xMath::GeVtoifm3();
    }

    double QuantumClusterExpansiondedmu(int statistics, double T, double mu, double m, double deg, int order,
                                       const IdealGasFunctionsExtraConfig &extraConfig) {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return Boltzmanndedmu(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double cfug = tfug;
      double moverT = m / T;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        if (m == 0.) {
          ret += sign * 3. / static_cast<double>(i) / static_cast<double>(i) / static_cast<double>(i) * cfug / T;
        } else {
          ret += sign * (xMath::BesselKexp(1, i*moverT) + 3. * xMath::BesselKexp(2, i*moverT) / moverT / static_cast<double>(i)) * cfug / T;
        }
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      double prefactor = 1.;
      if (m == 0.) {
        prefactor = deg * T * T * T * T / xMath::Pi() / xMath::Pi();
      } else {
        prefactor = deg * m * m * m * T / 2. / xMath::Pi() / xMath::Pi();
      }
      return ret * prefactor * xMath::GeVtoifm3();
    }


    double QuantumClusterExpansiondchi2dT(int statistics, double T, double mu, double m, double deg, int order,
                                       const IdealGasFunctionsExtraConfig &extraConfig) {
      bool signchange = true;
      if (statistics == 1) //Fermi
        signchange = true;
      else if (statistics == -1) //Bose
        signchange = false;
      else
        return Boltzmanndchi2dT(T, mu, m, deg, extraConfig);

      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      double tfug = exp((mu - m) / T);
      double moverT = m / T;
      double cfug = tfug;
      double sign = 1.;
      double ret = 0.;
      for (int i = 1; i <= order; ++i) {
        double k = static_cast<double>(i);
        if (m == 0.) {
          ret += sign * k * (-mu / k / k) * cfug;
        }
        else {
          ret += sign * k * (
                m * xMath::BesselKexp(1, i*moverT)
                - mu * xMath::BesselKexp(2, i*moverT)
                ) * cfug;
        }
        cfug *= tfug;
        if (signchange) sign = -sign;
      }
      double prefactor = 1.;
      if (m == 0.) {
        prefactor = deg / T / T / xMath::Pi() / xMath::Pi();
      } else {
        prefactor = deg * m * m / 2. / T / T / T / T / xMath::Pi() / xMath::Pi();
      }
      return ret * prefactor;
    }

    double QuantumNumericalIntegrationdndT(int statistics, double T, double mu, double m, double deg,
                                       const IdealGasFunctionsExtraConfig &extraConfig) {
      if (statistics == 0)            return BoltzmanndndT(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return 0.; //FermiZeroTDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m)  return FermiNumericalIntegrationLargeMudndT(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        std::cerr << "**WARNING** QuantumNumericalIntegrationDensity: Bose-Einstein condensation, mass = " << m << ", mu = " << mu << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double Eexp = exp(EoverT - muoverT);
        double f = 1. / (Eexp + statistics);
        ret += lagw32[i] * T * tx * T * tx * T * f * (EoverT - muoverT) * (1. - statistics * f);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() / T;

      return ret;
    }

    double QuantumNumericalIntegrationd2ndT2(int statistics, double T, double mu, double m, double deg,
                                       const IdealGasFunctionsExtraConfig &extraConfig) {
      if (statistics == 0)           return Boltzmannd2ndT2(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return 0.;//FermiZeroTDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMud2ndT2(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        std::cerr << "**WARNING** QuantumNumericalIntegrationDensity: Bose-Einstein condensation, mass = " << m << ", mu = " << mu << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double Eexp = exp(EoverT - muoverT);
        double f = 1. / (Eexp + statistics);
        ret += lagw32[i] * T * tx * T * tx * T * f
        * (EoverT - muoverT) * (1. - statistics * f)
        * ((EoverT - muoverT) * (1. - 2. * statistics * f) - 2.);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi()  * xMath::GeVtoifm3() / T / T;

      return ret;
    }

    double QuantumNumericalIntegrationdedT(int statistics, double T, double mu, double m, double deg,
                                        const IdealGasFunctionsExtraConfig &extraConfig) {
      if (statistics == 0)            return BoltzmanndedT(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return 0.; //FermiZeroTDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m)  return FermiNumericalIntegrationLargeMudedT(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        std::cerr << "**WARNING** QuantumNumericalIntegrationDensity: Bose-Einstein condensation, mass = " << m << ", mu = " << mu << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double Eexp = exp(EoverT - muoverT);
        double f = 1. / (Eexp + statistics);
        ret += lagw32[i] * T * tx * T * tx * T * EoverT * f * (EoverT - muoverT) * (1. - statistics * f);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationdedmu(int statistics, double T, double mu, double m, double deg,
                                        const IdealGasFunctionsExtraConfig &extraConfig) {
      if (statistics == 0)            return BoltzmanndedT(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return 0.; //FermiZeroTDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m)  return FermiNumericalIntegrationLargeMudedmu(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        std::cerr << "**WARNING** QuantumNumericalIntegrationDensity: Bose-Einstein condensation, mass = " << m << ", mu = " << mu << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double Eexp = exp(EoverT - muoverT);
        double f = 1. / (Eexp + statistics);
        ret += lagw32[i] * T * tx * T * tx * T * EoverT * f * (1. - statistics * f);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret;
    }

    double QuantumNumericalIntegrationdchi2dT(int statistics, double T, double mu, double m, double deg,
                                       const IdealGasFunctionsExtraConfig &extraConfig) {
      if (statistics == 0)           return Boltzmanndchi2dT(T, mu, m, deg, extraConfig);
      if (statistics == 1 && T == 0.) return 0.;//FermiZeroTDensity(mu, m, deg, extraConfig);
      if (statistics == 1 && mu > m) return FermiNumericalIntegrationLargeMudchi2dT(T, mu, m, deg, extraConfig);
      if (statistics == -1 && mu > m) {
        std::cerr << "**WARNING** QuantumNumericalIntegrationDensity: Bose-Einstein condensation, mass = " << m << ", mu = " << mu << std::endl;
        calculationHadBECIssue = true;
        mu = m;
        // return 0.;
      }
      if (statistics == -1 && T == 0.) return 0.;

      // Check for magnetic field effect for charged particles
      assert(extraConfig.MagneticField.B == 0);

      // No magnetic field
      double ret = 0.;
      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double Eexp = exp(EoverT - muoverT);
        double f = 1. / (Eexp + statistics);
        ret += lagw32[i] * T * tx * T * tx * T * f 
         * (1. - statistics * f)
        * ((EoverT - muoverT) * (1. - 2. * statistics * f) - 3.);
      }

      ret *= deg / 2. / xMath::Pi() / xMath::Pi() / T / T / T / T;

      return ret;
    }

    double FermiNumericalIntegrationLargeMudndT(double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig &extraConfig) {
      if (mu <= m)
        return QuantumNumericalIntegrationdndT(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double en = sqrt(p * p + m * m);
        double fbar = 1. / (exp(-(en - mu) / T) + 1.);
        ret1 += -legw32[i] * dpds * p * p * (mu - en) / T / T * fbar * (1. - fbar);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double f = 1. / (exp(EoverT - muoverT) + 1.);
        ret1 += lagw32[i] * T * tx * T * tx * T * (EoverT - muoverT) / T * f * (1. - f);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret1;
    }

    double FermiNumericalIntegrationLargeMud2ndT2(double T, double mu, double m, double deg,
                                                  const IdealGasFunctionsExtraConfig &extraConfig) {
      if (mu <= m)
        return QuantumNumericalIntegrationd2ndT2(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double en = sqrt(p * p + m * m);
        double fbar = 1. / (exp(-(en - mu) / T) + 1.);
        ret1 += -legw32[i] * dpds * p * p * (mu - en) / T / T / T
                * fbar * (1. - fbar) * ((mu - en) / T * (1.-2.*fbar) - 2.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double f = 1. / (exp(EoverT - muoverT) + 1.);
        ret1 += lagw32[i] * T * tx * T * tx * T * (EoverT - muoverT) / T  / T
                * f * (1. - f) * ((EoverT - muoverT) * (1. - 2. * f) - 2.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret1;
    }

    double FermiNumericalIntegrationLargeMudedT(double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig &extraConfig) {
      if (mu <= m)
        return QuantumNumericalIntegrationdedT(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double en = sqrt(p * p + m * m);
        double fbar = 1. / (exp(-(en - mu) / T) + 1.);
        ret1 += -legw32[i] * dpds * p * p * en * (mu - en) / T / T * fbar * (1. - fbar);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double f = 1. / (exp(EoverT - muoverT) + 1.);
        ret1 += lagw32[i] * T * tx * T * tx * T * EoverT * T * (EoverT - muoverT) / T * f * (1. - f);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3();

      return ret1;
    }

    double FermiNumericalIntegrationLargeMudedmu(double T, double mu, double m, double deg,
                                                 const IdealGasFunctionsExtraConfig &extraConfig) {
      if (mu <= m)
        return QuantumNumericalIntegrationdedmu(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double en = sqrt(p * p + m * m);
        double Eexp = exp(-(en - mu) / T);
        ret1 += legw32[i] * dpds * p * p * en * 1. / (1. + 1./Eexp) / (Eexp + 1.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double Eexp = exp(EoverT - muoverT);
        ret1 += lagw32[i] * T * tx * T * tx * T * EoverT * T * 1. / (1. + 1. / Eexp) / (Eexp + 1.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() * xMath::GeVtoifm3() / T;

      return ret1;
    }

    double FermiNumericalIntegrationLargeMudchi2dT(double T, double mu, double m, double deg,
                                                   const IdealGasFunctionsExtraConfig &extraConfig) {
      if (mu <= m)
        return QuantumNumericalIntegrationdchi2dT(1, T, mu, m, deg, extraConfig);

      assert(extraConfig.MagneticField.B == 0.0);

      double pf = sqrt(mu*mu - m * m);
      double ret1 = 0.;
      for (int i = 0; i < 32; i++) {
        double p, dpds;
        SommerfeldLegendreMap(legx32[i], pf, mu, T, p, dpds);
        double en = sqrt(p * p + m * m);
        double fbar = 1. / (exp(-(en - mu) / T) + 1.);
        ret1 += legw32[i] * dpds * p * p
                * fbar * (1. - fbar) * ((mu - en) / T * (1. - 2. * fbar) - 3.);
      }

      double moverT = m / T;
      double muoverT = mu / T;
      for (int i = 0; i < 32; i++) {
        double tx = pf / T + lagx32[i];
        double EoverT = sqrt(tx*tx + moverT * moverT);
        double f = 1. / (exp(EoverT - muoverT) + 1.);
        ret1 += lagw32[i] * T * tx * T * tx * T
                * f * (1. - f) * ((EoverT - muoverT) * (1. - 2. * f) - 3.);
      }

      ret1 *= deg / 2. / xMath::Pi() / xMath::Pi() / T / T / T / T;

      return ret1;
    }
  } // namespace IdealGasFunctions

} // namespace thermalfist
