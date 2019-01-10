/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELPARAMETERS_H
#define THERMALMODELPARAMETERS_H

namespace thermalfist {

  /**
   *   \brief Structure containing all thermal parameters of the model.
   * 
   *   Used in ThermalParticle class when calculating the ideal gas functions.
   */
  struct ThermalModelParameters {
    double    T;        /**< Temperature [GeV] */
    double    muB;      /**< Baryon chemical potential [GeV] */
    double    muS;      /**< Strangeness chemical potential [GeV] */
    double    muQ;      /**< Electric charge chemical potential [GeV] */
    double    muC;      /**< Charm chemical potential [GeV] */
    double    gammaq;   /**< Chemical non-equilibirum fugacity of light quarks */
    double    gammaS;   /**< Chemical non-equilibirum fugacity of strange quarks */
    double    gammaC;   /**< Chemical non-equilibirum fugacity of charm quarks */
    double    V;        /**< Total system volume [fm^3] */
    double    SVc;      /**< Canonical correlation volume [fm^3] */
    int       B;        /**< Total baryon charge (CE) */
    int       Q;        /**< Total electric charge (CE) */
    int       S;        /**< Total strangeness charge (CE) */
    int       C;        /**< Total charm charge (CE) */

    ThermalModelParameters(double pT = 0.155, double pmuB = 0.000, double pmuS = 0., double pmuQ = 0., double pgammaS = 1., double pV = 4000., double pSVc = 30., int pB = 2, int pQ = 2, int pS = 0, int pC = 0) :
      T(pT), muB(pmuB), muS(pmuS), muQ(pmuQ), muC(0.), gammaq(1.), gammaS(pgammaS), gammaC(1.), V(pV), SVc(pSVc), B(pB), Q(pQ), S(pS), C(pC) {
    }

    ThermalModelParameters(double pT, double pgammaS, double pV, int pB, int pQ, int pS, int pC = 0) :
      T(pT), gammaq(1.), gammaS(pgammaS), gammaC(1.), V(pV), SVc(pV), B(pB), Q(pQ), S(pS), C(pC) {
    }
  };

} // namespace thermalfist

#endif
