/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2016-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EXCLUDEDVOLUMEHELPER_H
#define EXCLUDEDVOLUMEHELPER_H

#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"
#include <cmath>

/**
 * \file ExcludedVolumeHelper.h
 * 
 * \brief Contains some functions do deal with excluded volumes. 
 * 
 */

namespace thermalfist {

  namespace CuteHRGHelper {
    /**
     * \brief Computes the excluded volume parameter from
     *        a given radius parameter value.
     * 
     * \param r       Radius parameter (fm)
     * \return double Excluded volume parameter (fm\f$^3\f$)
     */
    inline double vr(double r) { return (16. * xMath::Pi() / 3. * pow(r, 3)); }

    /**
     * \brief Computes the radius parameter from
     *        a given excluded volume parameter value.
     * 
     * \param v       Excluded volume parameter (fm\f$^3\f$)
     * \return double Radius parameter (fm)
     */
    inline double rv(double v) { return pow(v * 3. / (16. * xMath::Pi()), 1. / 3.); }

    /**
     * \brief Computes the symmetric 2nd virial coefficient
     *        \f$ b_{ij} \f$ of the classical hard spheres
     *        equation of state from the two radii.
     * 
     * \param r1       First radius (fm) 
     * \param r2       Second radius (fm) 
     * \return double  Virial coefficient (fm\f$^3\f$)
     */
    inline double brr(double r1, double r2) { return (2. * xMath::Pi() / 3.) * pow(r1 + r2, 3); }

    /**
     * \brief Returns the matrix of excluded volume
     *        coefficients \f$ \tilde{b}_{ij} \f$
     *        of an HRG model
     */
    std::vector< std::vector<double> > bijMatrix(const ThermalModelBase* model);

    /**
     * \brief Returns the matrix of van der Waals attraction
     *        coefficients \f$ a_{ij} \f$
     *        of an HRG model
     */
    std::vector< std::vector<double> > aijMatrix(const ThermalModelBase* model);

  }

} // namespace thermalfist

#endif