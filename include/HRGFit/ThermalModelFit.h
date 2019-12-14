/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELFIT_H
#define THERMALMODELFIT_H

/**
 * \file ThermalModelFit.h
 * \brief ThermalModelFit class
 * 
 */

#include <string>

#include "HRGFit/ThermalModelFitParameters.h"
#include "HRGFit/ThermalModelFitQuantities.h"
#include "HRGBase/xMath.h"
#include "HRGPCE/ThermalModelPCE.h"

namespace thermalfist {

  /**
   * \brief Class implementing the thermal model fit procedure.
   * 
   * Performs a fit within a generic HRG model to
   * a set of measured multiplicities and/or yield ratios.
   * This is achieved by minimizing
   * \f[
   *   \chi^2 = \frac{\chi^2}{N_{\rm dof}}
       ~=~\frac{1}{N_{\rm dof}}\sum_{i=1}^N\frac{\left(N_i^{\rm exp}~-~N_i^{\rm HRG}\right)^2}{\sigma_i^2}
   * \f]
   * with MINUIT2. 
   * 
   * The actual HRG model to use is provided through a pointer to the
   * corresponding thermalfist::ThermalModelBase object in the constructor.
   * 
   * The fitted thermal parameters can be specified in a thermalfist::ThermalModelFitParameters
   * object and should be provided through the SetParameters() method.
   * 
   * The data to fit can be specified as vector of FittedQuantity objects.
   * It can be provided through the SetQuantities() method,
   * the data can be read from a file with loadExpDataFromFile(). 
   * 
   * Finally, the fit procedure is invoked through the PerformFit() method.
   * 
   */
  class ThermalModelFit
  {
  public:
    /**
     * \brief Construct a new ThermalModelFit object
     * 
     * \param model Pointer to the ThermalModelBase object
     *              which implements the HRG model to use in fits
     */
    ThermalModelFit(ThermalModelBase *model);

    /// \brief Destroy the Thermal Model Fit object
    ~ThermalModelFit(void);

    /**
     * \brief Set whether FitParameter with a given name should be fitted.
     * 
     * \param name FitParameter name
     * \param flag true -- fitted, false -- not fitted
     */
    void SetFitFlag(const std::string& name, bool flag) {
      m_Parameters.SetParameterFitFlag(name, flag);
    }



    /// \brief Sets the data to fit
    /// \param inData A vector of measurements to fit
    void SetData(const std::vector<FittedQuantity> & inData) {
      m_Quantities.resize(0);
      m_Ratios.resize(0);
      m_Multiplicities.resize(0);
      AddData(inData);
    }

    /// Same as SetData()
    void SetQuantities(const std::vector<FittedQuantity> & inQuantities) { SetData(inQuantities); }

    /// \brief Add more data to fit
    ///
    /// \param inData A vector of additional measurements to fit
    void AddData(const std::vector<FittedQuantity> & inData) {
      for (size_t i = 0; i < inData.size(); ++i)
        AddDataPoint(inData[i]);
    }

    /// Same as AddData()
    void AddQuantities(const std::vector<FittedQuantity> & inQuantities) { AddData(inQuantities); }

    /// \brief Add one more data point to fit
    ///
    /// \param inDataPoint Data point to fit
    void AddDataPoint(const FittedQuantity& inDataPoint) {
      m_Quantities.push_back(inDataPoint);
      if (inDataPoint.type == FittedQuantity::Ratio)
        m_Ratios.push_back(inDataPoint.ratio);
      else
        m_Multiplicities.push_back(inDataPoint.mult);
    }

    /// Same as AddDataPoint()
    void AddQuantity(const FittedQuantity& inQuantity) { AddDataPoint(inQuantity); }

    /// Clear the fitted data
    void ClearData() { m_Quantities.resize(0); }

    /// Same as ClearData()
    void ClearQuantities() { ClearData(); }

    //@{
    /// \deprecated Use SetQuantities(), AddQuantities(), and AddQuantity() instead
    void SetRatios(const std::vector<ExperimentRatio> & inRatios) {
      m_Ratios = inRatios;
    }

    void AddRatios(const std::vector<ExperimentRatio> & inRatios) {
      m_Ratios.insert(m_Ratios.end(), inRatios.begin(), inRatios.end());
      for (size_t i = 0; i < inRatios.size(); i++) {
        m_Quantities.push_back(FittedQuantity(inRatios[i]));
      }
    }

    void AddRatio(const ExperimentRatio& inRatio) {
      m_Ratios.push_back(inRatio);
      m_Quantities.push_back(FittedQuantity(inRatio));
    }

    void ClearRatios() { m_Ratios.resize(0); }

    void PrintRatios();

    void SetMultiplicities(const std::vector<ExperimentMultiplicity> & inMultiplicities) {
      m_Multiplicities = inMultiplicities;
    }

    void AddMultiplicities(const std::vector<ExperimentMultiplicity> & inMultiplicities) {
      m_Multiplicities.insert(m_Multiplicities.end(), inMultiplicities.begin(), inMultiplicities.end());
      for (size_t i = 0; i < inMultiplicities.size(); i++) {
        m_Quantities.push_back(FittedQuantity(inMultiplicities[i]));
      }
    }

    void AddMultiplicity(const ExperimentMultiplicity& inMultiplicity) {
      m_Multiplicities.push_back(inMultiplicity);
      m_Quantities.push_back(FittedQuantity(inMultiplicity));
    }

    void ClearMultiplicities() { m_Multiplicities.resize(0); }
    //@}

    //@{
    /// Print function 
    void PrintParameters();

    void PrintMultiplicities();

    void PrintYields();

    void PrintYieldsTable(std::string filename = "Yield.dat");

    void PrintYieldsTable2(std::string filename = "Yield.dat");

    void PrintYieldsLatex(std::string filename = "Yield.dat", std::string name = "A+A");

    void PrintYieldsLatexAll(std::string filename = "Yield.dat", std::string name = "A+A", bool asymm = false);
    //@}

    /// Prints the result of the fitting procedure to a file
    void PrintFitLog(std::string filename = "", std::string comment = "Thermal fit", bool asymm = false);

    //double chi2Ndf(double T, double muB);

    /**
     * \brief The thermal fit procedure.
     * 
     * Performs a thermal fit of thermal parameters, specified by SetParameters(),
     * to the experimental data, provided through SetQuantities(),
     * within a HRG model provided by a pointer to thermalfist::ThermalModelBase in constructor.
     * 
     * Returns a thermalfist::ThermalModelFitParameters parameters object containing
     * the values and errors of the fitted thermal parameters resulting from
     * \f$ \chi^2 \f$ minimization.
     * 
     * \param verbose     If true, additional output is shown on screen during the fitting
     * \param AsymmErrors If true, asymmetric error bars are computed
     * \return ThermalModelFitParameters The fitted parameters
     */
    ThermalModelFitParameters PerformFit(bool verbose = true, bool AsymmErrors = false);

    /// Number of degrees of freedom in the fit
    int GetNdf() const;

    /// Used by MINUIT
    void Increment() { m_Iters++; }

    /// Current fit parameters.
    /// Contains the fit result after PerformFit() was called. 
    const ThermalModelFitParameters& Parameters() const { return m_Parameters; }

    /// Sets the fit parameters
    void SetParameters(const ThermalModelFitParameters& params) { m_Parameters = params; }

    /// Sets the fit parameter with a given name
    void SetParameter(const std::string & name, const FitParameter & param) { m_Parameters.SetParameter(name, param); }
    
    /// Sets the fit parameter with a given name
    void SetParameter(const std::string & name, double val, double err, double xmin, double xmax) { m_Parameters.SetParameter(name, val, err, xmin, xmax); }
    
    /// Sets the (initial) value for the fit parameter with a given name
    void SetParameterValue(const std::string & name, double value) { m_Parameters.SetParameterValue(name, value); }
    
    /// Sets whether the fit parameter with a given name is fitted
    void SetParameterFitFlag(const std::string & name, bool toFit) { m_Parameters.SetParameterFitFlag(name, toFit); }
    

    const ThermalModelFitParametersExtended& ExtendedParameters() const { return m_ExtendedParameters; }
    
    /// Pointer to a ThermalModelBase object implementing
    /// the HRG model used
    ThermalModelBase* model() { return m_model; }

    /// Pointer to a ThermalModelBase object implementing
    /// the partial chemical equilibrium HRG model used
    ThermalModelPCE* modelpce() { return m_modelpce; }

    /// Set the entropy per baryon constraint for \f$\mu_B \f$
    /// \deprecated Entropy per baryon constraint should be set directly on a ThermalModelBase object
    void SetSBConstraint(double SB) {
      if (m_model != NULL)
        m_model->SetSoverB(SB);
      m_model->ConstrainMuB(true);
    }

    /// Set the electric-to-baryon charge ratio constraint for \f$\mu_Q \f$
    /// \deprecated Electric-to-baryon charge ratio constraint should be set directly on a ThermalModelBase object
    void SetQBConstraint(double QB) {
      if (m_model != NULL)
        m_model->SetQoverB(QB);
      m_model->ConstrainMuQ(true);
    }
    
    /// \deprecated
    double SoverB() const { return m_model->SoverB(); }
    /// \deprecated
    double QoverB() const { return m_model->QoverB(); }

    /// Return a vector of the fitted multiplicities (yields)
    /// Ignores the ratios
    const std::vector<ExperimentMultiplicity>&    Multiplicities()   const { return m_Multiplicities; }
    
    /// Return a vector of the fitted yield ratios
    /// Ignores the multiplicities
    const std::vector<ExperimentRatio>&           Ratios()           const { return m_Ratios; }

    /// Return a vector of the fitted quantities (data)
    const std::vector<FittedQuantity>&            FittedQuantities() const { return m_Quantities; }

    /**@{ 
     * Used for real-time monitoring in the GUI.
     */
    int& Iters() { return m_Iters; }

    double& Chi2() { return m_Chi2; }

    double& BT() { return m_BT; }

    double& QT() { return m_QT; }

    double& ST() { return m_ST; }

    double& CT() { return m_CT; }

    double& ModelData(int index) { return m_ModelData[index]; }

    int ModelDataSize() const { return m_ModelData.size(); }

    void ClearModelData() { m_ModelData.clear(); }

    int& Ndf() { return m_Ndf; }
    /**@} */

    //@{
    /// Whether the correlation volume is tied to the total volume
    /// or an independent parameter.
    void FixVcOverV(bool fix) { m_FixVcToV = fix; }
    bool FixVcOverV() const { return m_FixVcToV; }
    //@}

    //@{
    /// The value of the correlation volume over the total volume
    /// ratio, \f$ V_c / V \f$.
    /// If FixVcOverV() is set to true this value is
    /// used fix \f$ V_c \f$ by the total volume.
    void SetVcOverV(double VcOverV) { m_VcOverV = VcOverV; }
    double VcOverV() const { return m_VcOverV; }
    //@}

    /// Sets whether the yields should be evaluated at Tkin using partial chemical equilibrium
    void UseTkin(bool YieldsAtTkin) { m_YieldsAtTkin = YieldsAtTkin; }
    bool UseTkin() const { return m_YieldsAtTkin; }

    /// Sets whether the nuclear abundances are evaluated in PCE using the Saha equation
    void UseSahaForNuclei(bool UseSaha) { m_SahaForNuclei = UseSaha; }

    /// Sets whether the yields of long-lived resonance are frozen in the PCE
    void PCEFreezeLongLived(bool FreezeLongLived) { m_PCEFreezeLongLived = FreezeLongLived; }

    /// Sets the resonance width cut for freezeing the yields of long-lived resonances
    void SetPCEWidthCut(double WidthCut) { m_PCEWidthCut = WidthCut; }

    /// Returns a relative error of the data description (and its uncertainty estimate)
    std::pair< double, double > ModelDescriptionAccuracy() const;

    /// Load the experimental data from a file.
    static std::vector<FittedQuantity> loadExpDataFromFile(const std::string & filename);

    static void saveExpDataToFile(const std::vector<FittedQuantity> & outQuantities, const std::string & filename);

  private:
    static std::string GetCurrentTime();
    
    static std::vector<FittedQuantity> loadExpDataFromFile_OldFormat(std::fstream & fin);

    static std::vector<FittedQuantity> loadExpDataFromFile_NewFormat(std::fstream & fin);

    ThermalModelFitParameters m_Parameters;
    ThermalModelFitParametersExtended m_ExtendedParameters;
    ThermalModelBase *m_model;
    ThermalModelPCE  *m_modelpce;
    std::vector<ExperimentMultiplicity> m_Multiplicities;
    std::vector<ExperimentRatio> m_Ratios;
    std::vector<FittedQuantity> m_Quantities;
    int       m_Iters;
    double    m_Chi2;
    double    m_BT;
    double    m_QT;
    double    m_ST;
    double    m_CT;
    std::vector<double> m_ModelData;
    int       m_Ndf;
    bool      m_FixVcToV;
    double    m_VcOverV;

    bool      m_YieldsAtTkin;
    bool      m_SahaForNuclei;
    bool      m_PCEFreezeLongLived;
    double    m_PCEWidthCut;
  };

} // namespace thermalfist

#endif
