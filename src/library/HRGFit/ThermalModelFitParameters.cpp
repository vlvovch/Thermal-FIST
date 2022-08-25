/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGFit/ThermalModelFitParameters.h"

namespace thermalfist {

  ThermalModelFitParametersExtended::ThermalModelFitParametersExtended(ThermalModelBase * model)
  {
    T.value = model->Parameters().T;
    muB.value = model->Parameters().muB;
    muS.value = model->Parameters().muS;
    muQ.value = model->Parameters().muQ;
    muC.value = model->Parameters().muC;
    gammaq.value = model->Parameters().gammaq;
    gammaS.value = model->Parameters().gammaS;
    gammaC.value = model->Parameters().gammaC;
    //R.value = model->Parameters.R;
    R.value = 0.;
    Rc.value = 0.;
    Tkin.value = model->Parameters().T;
    nH.value = model->CalculateHadronDensity();
    rhoB.value = model->CalculateBaryonDensity();
    rhoQ.value = model->CalculateChargeDensity();
    en.value = model->CalculateEnergyDensity();
    entropy.value = model->CalculateEntropyDensity();
    pressure.value = model->CalculatePressure();
  }

  ThermalModelFitParameters::ThermalModelFitParameters(const ThermalModelParameters & params)
  {
    T = FitParameter("T", true, params.T, 0.05, 0.02, 0.300);
    muB = FitParameter("muB", true, params.muB, 0.05, -0.100, 0.900);
    muS = FitParameter("muS", false, params.muS, 0.05, -0.450, 0.450);
    muQ = FitParameter("muQ", false, params.muQ, 0.05, -0.130, 0.130);
    muC = FitParameter("muC", false, params.muC);
    gammaq = FitParameter("gammaq", false, params.gammaq, 0.5, 0.01, 3.);
    gammaS = FitParameter("gammaS", false, params.gammaS, 0.5, 0.01, 3.);
    gammaC = FitParameter("gammaC", false, params.gammaC, 0.5, 0.01, 50.);
    //V = FitParameter("V", true, V_, 2000., 1., 20000.);  // Volume no longer used
    R = FitParameter("R", true, pow(3. * params.V / 16. / xMath::Pi(), 1. / 3.), 1.0, 0., 25.0);
    Rc = FitParameter("Rc", true, pow(3. * params.SVc / 16. / xMath::Pi(), 1. / 3.), 1.0, 0., 10.0);
    B = params.B;
    Q = params.Q;
    S = params.S;
    C = params.C;

    Tkin = FitParameter("Tkin", false, params.T, 0.05, 0.050, 0.200);

    FillParameterList();
  }

  ThermalModelFitParameters::ThermalModelFitParameters(const ThermalModelFitParameters & op) :
    GCE(op.GCE), T(op.T), muB(op.muB), muS(op.muS), muQ(op.muQ), muC(op.muC),
    gammaq(op.gammaq), gammaS(op.gammaS), gammaC(op.gammaC), R(op.R), Rc(op.Rc),
    Tkin(op.Tkin),
    B(op.B), S(op.S), Q(op.Q), C(op.C), chi2(op.chi2), chi2ndf(op.chi2ndf),
    ndf(op.ndf)
  {
    FillParameterList();
  }

  ThermalModelFitParameters & ThermalModelFitParameters::operator=(const ThermalModelFitParameters & op)
  {
    GCE = op.GCE;
    T = op.T;
    muB = op.muB;
    muS = op.muS;
    muQ = op.muQ;
    muC = op.muC;
    gammaq = op.gammaq;
    gammaS = op.gammaS;
    gammaC = op.gammaC;
    R = op.R;
    Rc = op.Rc;
    B = op.B;
    S = op.S;
    Q = op.Q;
    C = op.C;
    chi2 = op.chi2;
    chi2ndf = op.chi2ndf;
    ndf = op.ndf;

    Tkin = op.Tkin;

    FillParameterList();
    return *this;
  }

  void ThermalModelFitParameters::FillParameterList()
  {
    ParameterList.clear();
    ParameterList.push_back(&T);
    ParameterList.push_back(&R);
    ParameterList.push_back(&Rc);
    ParameterList.push_back(&muB);
    ParameterList.push_back(&muQ);
    ParameterList.push_back(&muS);
    ParameterList.push_back(&muC);
    ParameterList.push_back(&gammaq);
    ParameterList.push_back(&gammaS);
    ParameterList.push_back(&gammaC);
    
    ParameterList.push_back(&Tkin);
  }

  int ThermalModelFitParameters::IndexByName(const std::string & name) const
  {
    int ret = -1;
    for (size_t i = 0; i < ParameterList.size(); ++i)
      if (ParameterList[i]->name == name)
        ret = i;
    return ret;
  }

  FitParameter ThermalModelFitParameters::GetParameter(const std::string & name) const
  {
    int ind = IndexByName(name);
    if (ind != -1)
      return *ParameterList[ind];
    // return T by default
    return T;
  }

  FitParameter & ThermalModelFitParameters::GetParameter(const std::string & name)
  {
    int ind = IndexByName(name);
    if (ind != -1)
      return *ParameterList[ind];
    // return T by default
    return T;
  }

  FitParameter ThermalModelFitParameters::GetParameter(const int index) const
  {
    if (index >= 0 && index < static_cast<int>(ParameterList.size()))
      return *ParameterList[index];
    // return T by default
    return T;
  }

  FitParameter & ThermalModelFitParameters::GetParameter(const int index)
  {
    if (index >= 0 && index < static_cast<int>(ParameterList.size()))
      return *ParameterList[index];
    // return T by default
    return T;
  }

  void ThermalModelFitParameters::SetParameter(const std::string & name, const FitParameter & param)
  {
    int ind = IndexByName(name);
    if (ind != -1)
      *ParameterList[ind] = param;
  }

  void ThermalModelFitParameters::SetParameter(const std::string & name, double val, double err, double xmin, double xmax)
  {
    int ind = IndexByName(name);
    if (ind != -1)
      *ParameterList[ind] = FitParameter(name, ParameterList[ind]->toFit, val, err, xmin, xmax);
  }

  void ThermalModelFitParameters::SetParameterValue(const std::string & name, double value)
  {
    int ind = IndexByName(name);
    if (ind != -1)
      ParameterList[ind]->value = value;
  }

  void ThermalModelFitParameters::SetParameterFitFlag(const std::string & name, bool toFit)
  {
    int ind = IndexByName(name);
    if (ind != -1)
      ParameterList[ind]->toFit = toFit;
  }

  ThermalModelParameters ThermalModelFitParameters::GetThermalModelParameters()
  {
    ThermalModelParameters ret(T.value, muB.value, muS.value, muQ.value, gammaS.value, 4. / 3. * xMath::Pi() * R.value * R.value * R.value);
    ret.SVc = 4. / 3. * xMath::Pi() * Rc.value * Rc.value * Rc.value;
    ret.gammaq = gammaq.value;
    ret.muC = muC.value;
    ret.gammaC = gammaC.value;
    ret.B = B;
    ret.Q = Q;
    ret.S = S;
    ret.C = C;
    return ret;
  }

} // namespace thermalfist