/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELFITPARAMETERS_H
#define THERMALMODELFITPARAMETERS_H

#include <string>

#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  /**
   *  \brief Structure for an arbitrary fit parameter.
   */
  struct FitParameter {
    /// Parameter value.  
    /// If parameter is fitted this value is used as a starting point. If not fitted this value is used throughout.
    double value;       

    double error;       ///< Parameter error (symmetric)
    double errp, errm;  ///< Asymmetric errors
    double xmin, xmax;  ///< Lower and uppers bounds on parameter value
    bool toFit;         ///< Whether the parameter is fitted or fixed
    std::string name;   ///< Name of the parameter

    /**
     * \brief Construct a new FitParameter object
     */
    FitParameter(std::string namep = "", bool toFitp = true, double val = 1., double err = 0., double vmin = -2., double vmax = 2.) :
      toFit(toFitp), value(val), error(err), errp(err), errm(err), xmin(vmin), xmax(vmax), name(namep) {
    }
  };

  /**
   *  \brief  Extended structure for calculating uncertainties in non-fit quantities resulting from
   *          uncertanties in fit parameters.
   */
  struct ThermalModelFitParametersExtended {
    FitParameter T;
    FitParameter muB;
    FitParameter muS;
    FitParameter muQ;
    FitParameter muC;
    FitParameter gammaq;
    FitParameter gammaS;
    FitParameter gammaC;
    FitParameter R;
    FitParameter Rc;

    FitParameter Tkin;

    /// Reduced \f$ \chi^2 \f$
    double chi2ndf;

    /// Auxiliary quantities
    FitParameter nH, rhoB, rhoQ, en, entropy, pressure;
    ThermalModelFitParametersExtended() { }
    ThermalModelFitParametersExtended(ThermalModelBase *model);
  };

  /**
   *  \brief Structure holding information about parameters of a thermal fit.
   */
  struct ThermalModelFitParameters {
    static const int ParameterCount = 11;

    bool GCE;  ///< 0 - CE, 1 - GCE

    /// All the fit parameters
    FitParameter T, muB, muS, muQ, muC, gammaq, gammaS, gammaC, R, Rc;

    /// Since version 1.3: Kinetic freeze-out temperature
    FitParameter Tkin;

    /// Vector of pointer to all the parameters
    std::vector<FitParameter*> ParameterList;

    /// Total charges (for CE)
    int B, S, Q, C;

    /// Value of the \f$ \chi^2 \f$ function
    double chi2; 
    
    /// Reduced \f$ \chi^2 \f$
    double chi2ndf;

    /// Number of degrees of freedom
    int ndf;

    /**
     * \brief Default constructor
     * 
     * Constructs fit parameters from given thermodynamic parameters
     * 
     * \param params The thermodynamic parameters of an HRG
     */
    ThermalModelFitParameters(const ThermalModelParameters &params = ThermalModelParameters());

    /// Copy constructor
    ThermalModelFitParameters(const ThermalModelFitParameters& op);

    /// Assignment operator
    ThermalModelFitParameters& operator=(const ThermalModelFitParameters& op);

    /// Fill ParameterList with pointers to all the parameters
    void FillParameterList();

    /// Find the index of the parameter with a given name in the ParameterList
    int IndexByName(const std::string& name) const;

    //@{
    /**
     * \brief Get the FitParameter by its name
     * 
     * \param name FitParameter name
     * \return FitParameter 
     */
    FitParameter GetParameter(const std::string& name) const;
    FitParameter& GetParameter(const std::string& name);
    //@}

    //@{
    /**
     * \brief Get the FitParameter by its index
     * 
     * \param index FitParameter index
     * \return FitParameter 
     */
    FitParameter GetParameter(const int index) const;
    FitParameter& GetParameter(const int index);
    //@}

    /// Set the FitParameter with a given name from a copy
    void SetParameter(const std::string& name, const FitParameter& param);

    /// Set the FitParameter with a given name
    /// \param val FitParameter value
    /// \param err FitParameter error
    /// \param xmin FitParameter value lower bound
    /// \param xmax FitParameter value upper bound
    void SetParameter(const std::string& name, double val, double err, double xmin, double xmax);

    /// Set the value of the FitParameter with a given name
    void SetParameterValue(const std::string& name, double value);

    /// Set whether the FitParameter with a given name should be fitted
    void SetParameterFitFlag(const std::string& name, bool toFit);

    /// Return ThermalModelParameters corresponding to the this ThermalModelFitParameters
    ThermalModelParameters GetThermalModelParameters();
  };

} // namespace thermalfist

#endif