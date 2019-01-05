/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELFIT_H
#define THERMALMODELFIT_H

#include <string>

#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/xMath.h"

namespace thermalfist {

  /**
   *  Structure for an arbitrary fit parameter.
   */
  struct FitParameter {
    double value;
    double error;
    double errp, errm;
    double xmin, xmax;
    bool toFit;
    std::string name;
    FitParameter(std::string name_ = "", bool toFit_ = true, double val = 1., double err = 0., double vmin = -2., double vmax = 2.) :
      toFit(toFit_), value(val), error(err), errp(err), errm(-err), xmin(vmin), xmax(vmax), name(name_) {
    }
  };

  /**
   *  Extended structure for calculating uncertainties in non-fit quantities resulting from
   *  uncertanties in fit parameters.
   */
  struct ThermalModelFitParametersExtended {
    FitParameter T, muB, muS, muQ, muC, gammaq, gammaS, gammaC, R, Rc;
    double chi2ndf;
    FitParameter nH, rhoB, rhoQ, en, entropy, pressure;//, eta;
    ThermalModelFitParametersExtended() { }
    ThermalModelFitParametersExtended(ThermalModelBase *model);
  };

  /**
   *  Structure holding information about parameters of a thermal fit.
   */
  struct ThermalModelFitParameters {
    static const int ParameterCount = 10;
    bool GCE;  // 0 - CE, 1 - GCE
    FitParameter T, muB, muS, muQ, muC, gammaq, gammaS, gammaC, R, Rc;
    std::vector<FitParameter*> ParameterList;
    int B, S, Q, C;
    double chi2, chi2ndf;
    int ndf;

    // Default constructor
    ThermalModelFitParameters(const ThermalModelParameters &params = ThermalModelParameters());

    // Copy constructor
    ThermalModelFitParameters(const ThermalModelFitParameters& op);

    // Assignment operator
    ThermalModelFitParameters& operator=(const ThermalModelFitParameters& op);

    void FillParameterList();

    int IndexByName(const std::string& name) const;

    FitParameter GetParameter(const std::string& name) const;

    FitParameter& GetParameter(const std::string& name);

    FitParameter GetParameter(const int index) const;

    FitParameter& GetParameter(const int index);

    void SetParameter(const std::string& name, const FitParameter& param);

    void SetParameter(const std::string& name, double val, double err, double xmin, double xmax);

    void SetParameterValue(const std::string& name, double value);

    void SetParameterFitFlag(const std::string& name, bool toFit);

    ThermalModelParameters GetThermalModelParameters();
  };

  struct ExperimentMultiplicity {
    int fPDGID;
    double fValue, fError;
    Feeddown::Type fFeedDown; /// 0 - primordial, 1 - stability flags, 2 - strong + EM + weak, 3 - strong + EM, 4 - strong
    //int fFeedDown; /// 0 - primordial, 1 - +decays from unstable, 2 - +weak decays
    ExperimentMultiplicity(int PDGID = -211, double value = 300., double error = 20., Feeddown::Type fd = Feeddown::StabilityFlag) :
      fPDGID(PDGID), fValue(value), fError(error), fFeedDown(fd) { }
    void addSystematicError(double percent) {
      fError = sqrt(fError*fError + percent * percent*fValue*fValue);
    }
  };

  struct ExperimentRatio {
    int PDGID1, PDGID2;
    double fValue, fError;
    Feeddown::Type fFeedDown1, fFeedDown2; /// 0 - primordial, 1 - stability flags, 2 - strong + EM + weak, 3 - strong + EM, 4 - strong
    //int fFeedDown1, fFeedDown2; /// 0 - primordial, 1 - +decays from unstable, 2 - +weak decays
    ExperimentRatio(int PDGID1_ = 211, int PDGID2_ = -211, double value_ = 1., double error_ = 0.1, Feeddown::Type fd1 = Feeddown::StabilityFlag, Feeddown::Type fd2 = Feeddown::StabilityFlag) :
      PDGID1(PDGID1_), PDGID2(PDGID2_), fValue(value_), fError(error_), fFeedDown1(fd1), fFeedDown2(fd2) { }
    ExperimentRatio(int PDGID1_, int PDGID2_, double value1, double error1, double value2, double error2, Feeddown::Type fd1 = Feeddown::StabilityFlag, Feeddown::Type fd2 = Feeddown::StabilityFlag) :
      PDGID1(PDGID1_), PDGID2(PDGID2_), fFeedDown1(fd1), fFeedDown2(fd2) {
      fValue = value1 / value2;
      fError = sqrt(error1*error1 / value2 / value2 + value1 * value1 / value2 / value2 / value2 / value2 * error2 * error2);
    }
  };

  struct FittedQuantity {
    enum FittedQuantityType { Multiplicity = 0, Ratio = 1 };
    FittedQuantityType type; // 0 - Multiplicity, 1 - Ratio
    bool toFit;
    ExperimentMultiplicity mult;
    ExperimentRatio ratio;
    FittedQuantity() {
      toFit = true;
      type = FittedQuantity::Multiplicity;
      mult = ExperimentMultiplicity(-211, 10., 1.);
    }
    FittedQuantity(const ExperimentMultiplicity & op) {
      toFit = true;
      type = FittedQuantity::Multiplicity;
      mult = op;
    }
    FittedQuantity(const ExperimentRatio & op) {
      toFit = true;
      type = FittedQuantity::Ratio;
      ratio = op;
    }
    double Value() const {
      if (type == 0)
        return mult.fValue;
      else
        return ratio.fValue;
    }
    double ValueError() const {
      if (type == 0)
        return mult.fError;
      else
        return ratio.fError;
    }
  };

  class ThermalModelBase;

  class ThermalModelFit
  {
  public:
    ThermalModelFit(ThermalModelBase *model_);

    ~ThermalModelFit(void);

    void SetFitFlag(const std::string& name, bool flag) {
      m_Parameters.SetParameterFitFlag(name, flag);
    }

    void SetSBConstraint(double SB) {
      if (m_model != NULL)
        m_model->SetSoverB(SB);
    }

    void SetQBConstraint(double QB) {
      if (m_model != NULL)
        m_model->SetQoverB(QB);
    }

    void SetQuantities(const std::vector<FittedQuantity> & inQuantities) {
      m_Quantities.resize(0);
      m_Ratios.resize(0);
      m_Multiplicities.resize(0);
      AddQuantities(inQuantities);
    }

    void AddQuantities(const std::vector<FittedQuantity> & inQuantities) {
      for (int i = 0; i < inQuantities.size(); ++i)
        AddQuantity(inQuantities[i]);
    }

    void AddQuantity(const FittedQuantity& inQuantity) {
      m_Quantities.push_back(inQuantity);
      if (inQuantity.type == FittedQuantity::Ratio)
        m_Ratios.push_back(inQuantity.ratio);
      else
        m_Multiplicities.push_back(inQuantity.mult);
    }

    // To be deprecated
    void SetRatios(const std::vector<ExperimentRatio> & inRatios) {
      m_Ratios = inRatios;
    }

    void AddRatios(const std::vector<ExperimentRatio> & inRatios) {
      m_Ratios.insert(m_Ratios.end(), inRatios.begin(), inRatios.end());
      for (int i = 0; i < inRatios.size(); i++) {
        m_Quantities.push_back(FittedQuantity(inRatios[i]));
      }
    }

    void AddRatio(const ExperimentRatio& inRatio) {
      m_Ratios.push_back(inRatio);
      m_Quantities.push_back(FittedQuantity(inRatio));
    }

    void ClearRatios() { m_Ratios.resize(0); }

    void PrintRatios();

    // To be deprecated
    void SetMultiplicities(const std::vector<ExperimentMultiplicity> & inMultiplicities) {
      m_Multiplicities = inMultiplicities;
    }

    void AddMultiplicities(const std::vector<ExperimentMultiplicity> & inMultiplicities) {
      m_Multiplicities.insert(m_Multiplicities.end(), inMultiplicities.begin(), inMultiplicities.end());
      for (int i = 0; i < inMultiplicities.size(); i++) {
        m_Quantities.push_back(FittedQuantity(inMultiplicities[i]));
      }
    }

    void AddMultiplicity(const ExperimentMultiplicity& inMultiplicity) {
      m_Multiplicities.push_back(inMultiplicity);
      m_Quantities.push_back(FittedQuantity(inMultiplicity));
    }

    void ClearMultiplicities() { m_Multiplicities.resize(0); }

    void PrintParameters();

    void PrintMultiplicities();

    void PrintYields();

    void PrintYieldsTable(std::string filename = "Yield.dat");

    void PrintYieldsTable2(std::string filename = "Yield.dat");

    void PrintYieldsLatex(std::string filename = "Yield.dat", std::string name = "A+A");

    void PrintYieldsLatexAll(std::string filename = "Yield.dat", std::string name = "A+A", bool asymm = false);

    static std::string GetCurrentTime();
    void PrintFitLog(std::string filename = "", std::string comment = "Thermal fit", bool asymm = false);

    double chi2Ndf(double T, double muB);

    ThermalModelFitParameters PerformFit(bool verbose = true, bool AsymmErrors = false);

    int GetNdf() const;

    void Increment() { m_Iters++; }

    const ThermalModelFitParameters& Parameters() const { return m_Parameters; }
    void SetParameters(const ThermalModelFitParameters& params) { m_Parameters = params; }
    void SetParameter(const std::string & name, const FitParameter & param) { m_Parameters.SetParameter(name, param); }
    void SetParameter(const std::string & name, double val, double err, double xmin, double xmax) { m_Parameters.SetParameter(name, val, err, xmin, xmax); }
    void SetParameterValue(const std::string & name, double value) { m_Parameters.SetParameterValue(name, value); }
    void SetParameterFitFlag(const std::string & name, bool toFit) { m_Parameters.SetParameterFitFlag(name, toFit); }
    const ThermalModelFitParametersExtended& ExtendedParameters() const { return m_ExtendedParameters; }
    ThermalModelBase* model() { return m_model; }

    double SoverB() const { return m_model->SoverB(); }
    double QoverB() const { return m_model->QoverB(); }

    const std::vector<ExperimentMultiplicity>&    Multiplicities()   const { return m_Multiplicities; }
    const std::vector<ExperimentRatio>&           Ratios()           const { return m_Ratios; }
    const std::vector<FittedQuantity>&            FittedQuantities() const { return m_Quantities; }

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

    bool FixVcOverV() const { return m_FixVcToV; }
    void FixVcOverV(bool fix) { m_FixVcToV = fix; }
    double VcOverV() const { return m_VcOverV; }
    void SetVcOverV(double VcOverV) { m_VcOverV = VcOverV; }


    static std::vector<FittedQuantity> loadExpDataFromFile(const std::string & filename);
    static std::vector<FittedQuantity> loadExpDataFromFile_OldFormat(std::fstream & fin);
    static std::vector<FittedQuantity> loadExpDataFromFile_NewFormat(std::fstream & fin);
    static void saveExpDataToFile(const std::vector<FittedQuantity> & outQuantities, const std::string & filename);

  private:
    ThermalModelFitParameters m_Parameters;
    ThermalModelFitParametersExtended m_ExtendedParameters;
    ThermalModelBase *m_model;
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
  };

} // namespace thermalfist

#endif
