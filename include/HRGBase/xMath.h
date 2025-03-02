/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALFIST_XMATH_H
#define THERMALFIST_XMATH_H
/**
 * \file xMath.h
 * 
 * \brief Contains some extra mathematical functions used in the code.
 * 
 */


namespace thermalfist {

  namespace xMath {

    /// Pi constant
    constexpr double Pi() { return 3.14159265358979323846; }
    /// A constant to transform GeV into fm\f$^{-1}\f$.
    constexpr double GeVtoifm() { return 5.06773; }

    /// A constant to transform GeV\f$^{2}\f$ into fm\f$^{-2}\f$.
    constexpr double GeVtoifm2() { return GeVtoifm() * GeVtoifm(); }

    /// A constant to transform GeV\f$^{3}\f$ into fm\f$^{-3}\f$.
    constexpr double GeVtoifm3() { return GeVtoifm() * GeVtoifm() * GeVtoifm(); }

    /// Nucleon's mass. Value as in UrQMD.
    constexpr double mnucleon() { return 0.938; }

    /// Pion's mass. Value as in UrQMD.
    constexpr double mpion() { return 0.138; }

    //@{
    /// Bessel and related special functions.
    /// Implementation of these special functions is adapted from
    /// the CERN-ROOT package: https://root.cern.ch/
    double BesselI(int n, double x);   ///< integer order modified Bessel function I_n(x)
    double BesselK(int n, double x);   ///< integer order modified Bessel function K_n(x)
    double BesselI0(double x);         ///< modified Bessel function I_0(x)
    double BesselK0(double x);         ///< modified Bessel function K_0(x)
    double BesselI1(double x);         ///< modified Bessel function I_1(x)
    double BesselK1(double x);         ///< modified Bessel function K_1(x)
    double BesselJ0(double x);         ///< Bessel function J0(x) for any real x
    double BesselJ1(double x);         ///< Bessel function J1(x) for any real x
    double BesselY0(double x);         ///< Bessel function Y0(x) for positive x
    double BesselY1(double x);         ///< Bessel function Y1(x) for positive x
    double StruveH0(double x);         ///< Struve functions of order 0
    double StruveH1(double x);         ///< Struve functions of order 1
    double StruveL0(double x);         ///< Modified Struve functions of order 0
    double StruveL1(double x);         ///< Modified Struve functions of order 1

    double BesselK0exp(double x);         ///< modified Bessel function K_0(x), divided by exponential factor
    double BesselK1exp(double x);         ///< modified Bessel function K_1(x), divided by exponential factor
    double BesselKexp(int n, double x);   ///< integer order modified Bessel function K_n(x), divided by exponential factor

    double BesselI0exp(double x);         ///< modified Bessel function I_0(x), divided by exponential factor
    double BesselI1exp(double x);         ///< modified Bessel function I_1(x), divided by exponential factor
    double BesselIexp(int n, double x);   ///< integer order modified Bessel function I_n(x), divided by exponential factor


    // Note that the functions Gamma and LogGamma are mutually dependent.
    double LogGamma(double);
    double Gamma(double);
    //@}

    /**
     * \brief Computes the Lambert W function (0-branch) using Halley's method.
     * 
     * The Lambert W function is the inverse function of \( f(W) = W e^W \).
     * This implementation uses Halley's method to achieve the desired accuracy within 10 * epsilon,
     * where epsilon is the machine precision. The initial guess for the iteration is z = 0.
     * 
     * \tparam T The type of the input value.
     * \param z The input value for which the Lambert W function is computed.
     * \return The computed value of the Lambert W function for the given input.
     */
    template<typename T> T LambertW0(T z);
  }

} // namespace thermalfist

// Implementation of template functions
#include "xMath.tcc"

#endif // THERMALFIST_XMATH_H
