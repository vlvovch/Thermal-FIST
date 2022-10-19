/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ThermalParticle.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "HRGBase/Utility.h"
#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"
#include "HRGBase/ThermalParticleSystem.h"

using namespace std;

namespace thermalfist {

  ThermalParticle::ThermalParticle(bool Stable, std::string Name, long long PDGID, double Deg, int Stat, double Mass,
    int Strange, int Baryon, int Charge, double AbsS, double Width, double Threshold, int Charm, double AbsC, int Quark) :
    m_Stable(Stable), m_AntiParticle(false), m_Name(Name), m_PDGID(PDGID), m_Degeneracy(Deg), m_Statistics(Stat), m_StatisticsOrig(Stat), m_Mass(Mass),
    m_Strangeness(Strange), m_Baryon(Baryon), m_ElectricCharge(Charge), m_Charm(Charm), m_ArbitraryCharge(Baryon), m_AbsS(AbsS), m_AbsC(AbsC), m_Width(Width), m_Threshold(Threshold), m_Quark(Quark), m_Weight(1.),
    m_ResonanceWidthShape(ThermalParticle::RelativisticBreitWigner), m_ResonanceWidthIntegrationType(ThermalParticle::ZeroWidth)
  {
    if (!Disclaimer::DisclaimerPrinted) 
      Disclaimer::DisclaimerPrinted = Disclaimer::PrintDisclaimer();
    
    SetCalculationType(IdealGasFunctions::Quadratures);

    SetClusterExpansionOrder(3);
    if (m_Mass < 1.000) SetClusterExpansionOrder(5);
    if (m_Mass < 0.200) SetClusterExpansionOrder(10);

    SetResonanceWidthShape(RelativisticBreitWigner);
    //SetResonanceWidthIntegrationType(BWTwoGamma);
    SetResonanceWidthIntegrationType(ZeroWidth);

    m_DecayType = ParticleDecayType::Default;

    FillCoefficients();

    SetAbsoluteQuark(GetAbsQ());

  }


  ThermalParticle::~ThermalParticle(void)
  {
  }

  bool ThermalParticle::ZeroWidthEnforced() const
  {
    //if (PdgId() == 223) // omega(782)
    //  return true;
    return ((ResonanceWidth() / Mass()) < 0.01);
  }

  void ThermalParticle::SetResonanceWidth(double width)
  {
    m_Width = width;
    if (m_Width != 0.0) {
      FillCoefficients();
      FillCoefficientsDynamical();
    }
  }

  void ThermalParticle::SetDecayThresholdMass(double threshold)
  {
    m_Threshold = threshold;
    if (m_Threshold < 0.0) {
      printf("**WARNING** Trying to set negative decay threshold for %s, setting to zero instead", m_Name.c_str());
    }
    if (m_Width != 0.0) FillCoefficients();
  }

  void ThermalParticle::SetDecayThresholdMassDynamical(double threshold)
  {
    m_ThresholdDynamical = threshold;
    if (m_ThresholdDynamical < 0.0) {
      printf("**WARNING** Trying to set negative dynamical decay threshold for %s, setting to zero instead", m_Name.c_str());
    }
  }

  void ThermalParticle::CalculateAndSetDynamicalThreshold()
  {
    double Thr = m_Mass + m_Width;
    for (size_t i = 0; i < m_Decays.size(); ++i) {
      Thr = min(Thr, m_Decays[i].mM0);
    }
    m_ThresholdDynamical = Thr;
  }

  void ThermalParticle::SetResonanceWidthShape(ResonanceWidthShape shape)
  {
    if (shape != m_ResonanceWidthShape) {
      m_ResonanceWidthShape = shape;
      FillCoefficientsDynamical();
    }
  }

  void ThermalParticle::SetResonanceWidthIntegrationType(ResonanceWidthIntegration type)
  {
    m_ResonanceWidthIntegrationType = type;
    FillCoefficients();
    if (type == ThermalParticle::eBW || type == ThermalParticle::eBWconstBR)
      FillCoefficientsDynamical();
  }

  double ThermalParticle::MassDistribution(double m) const
  {
    return MassDistribution(m, m_Width);
  }

  double ThermalParticle::MassDistribution(double m, double width) const
  {
    if (width < 0.) width = m_Width;

    // Zero if outside the interval
    //if (m_ResonanceWidthIntegrationType == BWTwoGamma) {
    //  double a = max(m_Threshold, m_Mass - 2.*m_Width);
    //  double b = m_Mass + 2.*m_Width;
    //  if (m < a || m > b)
    //    return 0.;
    //}

    if (m_ResonanceWidthShape == RelativisticBreitWigner)
      return m_Mass * width * m / ((m * m - m_Mass * m_Mass)*(m * m - m_Mass * m_Mass) + m_Mass * m_Mass*width*width);
    else
      return width / ((m - m_Mass)*(m - m_Mass) + width * width / 4.);
  }

  void ThermalParticle::ReadDecays(string filename) {
    m_Decays.resize(0);
    ifstream fin(filename.c_str());
    if (fin.is_open()) {
      char cc[400];
      double tmpbr;
      while (fin >> tmpbr) {
        ParticleDecayChannel decay;
        decay.mBratio = tmpbr / 100.;
        fin.getline(cc, 350);
        stringstream ss;
        ss << cc;
        int tmpid;
        while (ss >> tmpid) {
          decay.mDaughters.push_back(tmpid);
        }
        m_Decays.push_back(decay);
      }
    }
  }

  void ThermalParticle::CalculateThermalBranchingRatios(const ThermalModelParameters & params, bool useWidth, double mu)
  {
    if (!useWidth || m_Width == 0.0 || m_Width / m_Mass < 1.e-2 || m_ResonanceWidthIntegrationType != eBW) {
      for (size_t j = 0; j < m_Decays.size(); ++j) {
        m_Decays[j].mBratioAverage = m_Decays[j].mBratio;
      }
    }
    else {
      if (!(params.gammaq == 1.))                  mu += log(params.gammaq) * GetAbsQ()  * params.T;
      if (!(params.gammaS == 1. || m_AbsS == 0.))  mu += log(params.gammaS) * m_AbsS     * params.T;
      if (!(params.gammaC == 1. || m_AbsC == 0.))  mu += log(params.gammaC) * m_AbsC     * params.T;

      for (size_t j = 0; j < m_Decays.size(); ++j) {
        m_Decays[j].mBratioAverage = 0.;
      }

      double ret1 = 0., ret2 = 0., tmp = 0.;
      for (size_t i = 0; i < m_xalldyn.size(); i++) {
        tmp = m_walldyn[i];
        double dens = IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ParticleDensity, m_QuantumStatisticsCalculationType, m_Statistics, params.T, mu, m_xalldyn[i], m_Degeneracy, m_ClusterExpansionOrder);
        ret1 += tmp * dens;
        ret2 += tmp;

        for (size_t j = 0; j < m_Decays.size(); ++j) {
          const_cast<double&>(m_Decays[j].mBratioAverage) += tmp * dens * m_Decays[j].mBratioVsM[i];
        }
      }

      for (size_t j = 0; j < m_Decays.size(); ++j) {
        if (ret1 != 0.0)
          m_Decays[j].mBratioAverage /= ret1;
        else
          m_Decays[j].mBratioAverage = m_Decays[j].mBratio;
      }
    }
  }

  std::vector<double> ThermalParticle::BranchingRatioWeights(const std::vector<double>& ms) const
  {
    std::vector<double> ret = ms;

    std::vector<double> mthr, br;
    double tsum = 0.;
    for (size_t i = 0; i < m_Decays.size(); ++i) {
      mthr.push_back(m_Decays[i].mM0);
      br.push_back(m_Decays[i].mBratio);
      tsum += m_Decays[i].mBratio;
    }
    if (tsum < 1.) {
      mthr.push_back(0.);
      br.push_back(1. - tsum);
    }
    else if (tsum > 1.) {
      for (size_t i = 0; i < br.size(); ++i)
        br[i] *= 1. / tsum;
    }

    for (size_t i = 0; i < ms.size(); ++i) {
      double tw = 0.;
      for (size_t j = 0; j < br.size(); ++j) {
        if (mthr[j] <= ms[i])
          tw += br[j];
      }
      ret[i] = tw;
    }
    return ret;
  }

  ThermalParticle ThermalParticle::GenerateAntiParticle(/*ThermalParticleSystem *TPS*/) const
  {
    string name = Name();
    if (BaryonCharge() == 0 && name[name.size() - 1] == '+')
      name[name.size() - 1] = '-';
    else if (BaryonCharge() == 0 && name[name.size() - 1] == '-')
      name[name.size() - 1] = '+';
    else
      name = "anti-" + name;

    ThermalParticle ret = *this;
    ret.SetName(name);
    ret.SetPdgId(-PdgId());
    ret.SetBaryonCharge(-BaryonCharge());
    ret.SetStrangenessCharge(-Strangeness());
    ret.SetElectricCharge(-ElectricCharge());
    ret.SetCharm(-Charm());
    ret.SetArbitraryCharge(-ArbitraryCharge());
    ret.SetAntiParticle(true);
    // Decays are to be done separately from ThermalParticleSystem instance
    // The code below is deprecated
    /*if (TPS != NULL) {
      ret.SetDecays(TPS->GetDecaysFromAntiParticle(Decays()));
      ret.SetDecaysOriginal(ret.Decays());
    }*/
    return ret;
  }

  bool ThermalParticle::operator==(const ThermalParticle & rhs) const
  {
    bool ret = true;

    /*ret &= m_xlag32 == rhs.m_xlag32;
    ret &= m_wlag32 == rhs.m_wlag32;
    ret &= m_xleg == rhs.m_xleg;
    ret &= m_wleg == rhs.m_wleg;
    ret &= m_xleg32 == rhs.m_xleg32;
    ret &= m_wleg32 == rhs.m_wleg32;*/

    ret &= m_Stable == rhs.m_Stable;
    ret &= m_AntiParticle == rhs.m_AntiParticle;
    ret &= m_Name == rhs.m_Name;
    ret &= m_PDGID == rhs.m_PDGID;
    ret &= m_Degeneracy == rhs.m_Degeneracy;
    ret &= m_Statistics == rhs.m_Statistics;
    ret &= m_StatisticsOrig == rhs.m_StatisticsOrig;
    ret &= m_Mass == rhs.m_Mass;
    ret &= m_QuantumStatisticsCalculationType == rhs.m_QuantumStatisticsCalculationType;
    ret &= m_ClusterExpansionOrder == rhs.m_ClusterExpansionOrder;
    ret &= m_Baryon == rhs.m_Baryon;
    ret &= m_ElectricCharge == rhs.m_ElectricCharge;
    ret &= m_Strangeness == rhs.m_Strangeness;
    ret &= m_Charm == rhs.m_Charm;
    ret &= m_Quark == rhs.m_Quark;
    ret &= m_ArbitraryCharge == rhs.m_ArbitraryCharge;
    ret &= m_AbsQuark == rhs.m_AbsQuark;
    ret &= m_AbsS == rhs.m_AbsS;
    ret &= m_AbsC == rhs.m_AbsC;
    ret &= m_Width == rhs.m_Width;
    ret &= m_Threshold == rhs.m_Threshold;
    ret &= m_ThresholdDynamical == rhs.m_ThresholdDynamical;
    ret &= m_ResonanceWidthShape == rhs.m_ResonanceWidthShape;
    ret &= m_ResonanceWidthIntegrationType == rhs.m_ResonanceWidthIntegrationType;
    ret &= m_Weight == rhs.m_Weight;
    ret &= m_DecayType == rhs.m_DecayType;
    ret &= m_Decays == rhs.m_Decays;
    ret &= m_DecaysOrig == rhs.m_DecaysOrig;
    
    return ret;
  }

  void ThermalParticle::NormalizeBranchingRatios() {
    double sum = 0.;
    for (size_t i = 0; i < m_Decays.size(); ++i) {
      sum += m_Decays[i].mBratio;
    }
    for (size_t i = 0; i < m_Decays.size(); ++i) {
      m_Decays[i].mBratio *= 1. / sum;
    }
    FillCoefficientsDynamical();
  }

  void ThermalParticle::RestoreBranchingRatios()
  {
    //m_Decays = m_DecaysOrig;
    for (size_t i = 0; i < m_Decays.size(); ++i) {
      m_Decays[i].mBratio = m_DecaysOrig[i].mBratio;
    }
    FillCoefficientsDynamical();
  }

  void ThermalParticle::FillCoefficients() {
    double a, b;
    if (m_ResonanceWidthIntegrationType != BWTwoGamma && m_Threshold >= 0.) {
      a = m_Threshold;
      b = m_Mass + 2.*m_Width;
      NumericalIntegration::GetCoefsIntegrateLegendre32(a, b, &m_xleg, &m_wleg);
      m_brweight = BranchingRatioWeights(m_xleg);
    }
    else {
      a = max(m_Threshold, m_Mass - 2.*m_Width);
      b = m_Mass + 2.*m_Width;
      NumericalIntegration::GetCoefsIntegrateLegendre10(a, b, &m_xleg, &m_wleg);
      m_brweight = BranchingRatioWeights(m_xleg);
    }

    // Old version
    //if (m_Width / m_Mass<1e-1) { NumericalIntegration::GetCoefsIntegrateLegendre10(a, b, &m_xleg, &m_wleg); }
    //else { NumericalIntegration::GetCoefsIntegrateLegendre32(a, b, &m_xleg, &m_wleg); }

    // New version
    NumericalIntegration::GetCoefsIntegrateLegendre32(0., 1., &m_xleg32, &m_wleg32);
    NumericalIntegration::GetCoefsIntegrateLaguerre32(&m_xlag32, &m_wlag32);
  }

  // Mass-dependent widths
  void ThermalParticle::FillCoefficientsDynamical() {
    if (m_Width == 0.0) return;

    double a, b;

    if (m_Decays.size() == 0)
      m_Threshold = m_ThresholdDynamical = m_Mass - 2.*m_Width + 1.e-6;

    //a = max(m_Threshold, m_ThresholdDynamical);
    a = m_ThresholdDynamical;

    b = m_Mass + 2.*m_Width;
    if (a >= m_Mass - 2.*m_Width) {
      m_xlegpdyn.resize(0);
      if (a >= m_Mass + 2.*m_Width)
        m_xlegdyn.resize(0);
      else
        NumericalIntegration::GetCoefsIntegrateLegendre32(a, b, &m_xlegdyn, &m_wlegdyn);
    }
    else {
      NumericalIntegration::GetCoefsIntegrateLegendre32(a, m_Mass - 2.*m_Width, &m_xlegpdyn, &m_wlegpdyn);
      NumericalIntegration::GetCoefsIntegrateLegendre32(m_Mass - 2.*m_Width, b, &m_xlegdyn, &m_wlegdyn);
    }
    NumericalIntegration::GetCoefsIntegrateLaguerre32(&m_xlagdyn, &m_wlagdyn);

    m_vallegpdyn = m_xlegpdyn;
    m_vallegdyn = m_xlegdyn;
    m_vallagdyn = m_xlagdyn;


    double tsumb = 0.;
    double tC = 0.;
    vector<double> tCP(m_Decays.size(), 0.);

    for (size_t i = 0; i < m_Decays.size(); ++i) {
      tsumb += m_Decays[i].mBratio;
      m_Decays[i].mBratioVsM.resize(0);
    }

    for (size_t j = 0; j < m_xlegpdyn.size(); ++j) {
      double twid = 0.;

      for (size_t i = 0; i < m_Decays.size(); ++i) {
        twid += m_Decays[i].ModifiedWidth(m_xlegpdyn[j]) * m_Width;
      }

      if (tsumb < 1.)
        twid += (1. - tsumb) * m_Width;

      if (twid == 0.0) {
        m_vallegpdyn[j] = 0.;
        for (size_t i = 0; i < m_Decays.size(); ++i)
          m_Decays[i].mBratioVsM.push_back(m_Decays[i].mBratio);
        continue;
      }

      for (size_t i = 0; i < m_Decays.size(); ++i) {
        double ttwid = m_Decays[i].ModifiedWidth(m_xlegpdyn[j]) * m_Width;
        m_Decays[i].mBratioVsM.push_back(ttwid / twid);
        tCP[i] += m_wlegpdyn[j] * ttwid / twid * MassDistribution(m_xlegpdyn[j], twid);
      }

      tC += m_wlegpdyn[j] * MassDistribution(m_xlegpdyn[j], twid);
      m_vallegpdyn[j] = MassDistribution(m_xlegpdyn[j], twid);
    }

    for (size_t j = 0; j < m_xlegdyn.size(); ++j) {
      double twid = 0.;

      for (size_t i = 0; i < m_Decays.size(); ++i) {
        twid += m_Decays[i].ModifiedWidth(m_xlegdyn[j]) * m_Width;
      }

      if (tsumb < 1.)
        twid += (1. - tsumb) * m_Width;

      if (twid == 0.0) {
        m_vallegdyn[j] = 0.;
        for (size_t i = 0; i < m_Decays.size(); ++i)
          m_Decays[i].mBratioVsM.push_back(m_Decays[i].mBratio);
        continue;
      }

      for (size_t i = 0; i < m_Decays.size(); ++i) {
        double ttwid = m_Decays[i].ModifiedWidth(m_xlegdyn[j]) * m_Width;
        m_Decays[i].mBratioVsM.push_back(ttwid / twid);
        tCP[i] += m_wlegdyn[j] * ttwid / twid * MassDistribution(m_xlegdyn[j], twid);
      }

      tC += m_wlegdyn[j] * MassDistribution(m_xlegdyn[j], twid);
      m_vallegdyn[j] = MassDistribution(m_xlegdyn[j], twid);
    }

    for (size_t j = 0; j < m_xlagdyn.size(); ++j) {
      double twid = 0.;
      double tx = m_Mass + 2.*m_Width + m_xlagdyn[j] * m_Width;

      for (size_t i = 0; i < m_Decays.size(); ++i) {
        twid += m_Decays[i].ModifiedWidth(tx) * m_Width;
      }

      if (tsumb < 1.)
        twid += (1. - tsumb) * m_Width;

      if (twid == 0.0) {
        m_vallagdyn[j] = 0.;
        for (size_t i = 0; i < m_Decays.size(); ++i)
          m_Decays[i].mBratioVsM.push_back(m_Decays[i].mBratio);
        continue;
      }

      for (size_t i = 0; i < m_Decays.size(); ++i) {
        double ttwid = m_Decays[i].ModifiedWidth(tx) * m_Width;
        m_Decays[i].mBratioVsM.push_back(ttwid / twid);
        tCP[i] += m_wlagdyn[j] * m_Width * ttwid / twid * MassDistribution(tx, twid);
      }

      tC += m_wlagdyn[j] * m_Width * MassDistribution(tx, twid);
      m_vallagdyn[j] = m_Width * MassDistribution(tx, twid);
    }

    m_xalldyn.resize(0);
    m_walldyn.resize(0);

    for (size_t j = 0; j < m_xlegpdyn.size(); ++j) {
      m_xalldyn.push_back(m_xlegpdyn[j]);
      m_walldyn.push_back(m_wlegpdyn[j] * m_vallegpdyn[j] / tC);
    }

    for (size_t j = 0; j < m_xlegdyn.size(); ++j) {
      m_xalldyn.push_back(m_xlegdyn[j]);
      m_walldyn.push_back(m_wlegdyn[j] * m_vallegdyn[j] / tC);
    }

    for (size_t j = 0; j < m_xlagdyn.size(); ++j) {
      m_xalldyn.push_back(m_Mass + 2.*m_Width + m_xlagdyn[j] * m_Width);
      m_walldyn.push_back(m_wlagdyn[j] * m_vallagdyn[j] / tC);
    }

    m_densalldyn.resize(m_xalldyn.size());

    double tsum = 0.;
    for (size_t j = 0; j < m_walldyn.size(); ++j) {
      tsum += m_walldyn[j];
    }


    // // Output to file for debugging
    // {
    //   char cc[300];
    //   sprintf(cc, "%d_width.dat", PdgId());
    //   FILE *f = fopen(cc, "w");

    //   for (int j = 0; j < m_xlegpdyn.size(); ++j)
    //     fprintf(f, "%15lf%15E\n", m_xlegpdyn[j], m_vallegpdyn[j] / tC);

    //   for (int j = 0; j < m_xlegdyn.size(); ++j)
    //     fprintf(f, "%15lf%15E\n", m_xlegdyn[j], m_vallegdyn[j] / tC);

    //   for (int j = 0; j < m_xlagdyn.size(); ++j)
    //     fprintf(f, "%15lf%15E\n", m_Mass + 2.*m_Width + m_xlagdyn[j] * m_Width, m_vallagdyn[j] / m_Width / tC);

    //   fclose(f);
    // }

  }

  double ThermalParticle::TotalWidtheBW(double M) const
  {
    //if (m_ResonanceWidthIntegrationType != eBW)
    //  return m_Width;

    double tsumb = 0.0;
    for (size_t i = 0; i < m_Decays.size(); ++i) {
      tsumb += m_Decays[i].mBratio;
    }

    double twid = 0.0;
    for (size_t i = 0; i < m_Decays.size(); ++i) {
      twid += m_Decays[i].ModifiedWidth(M) * m_Width;
    }

    if (tsumb < 1.)
      twid += (1. - tsumb) * m_Width;

    return twid;
  }

  std::vector<double> ThermalParticle::BranchingRatiosM(double M, bool eBW) const
  {
    std::vector<double> ret(m_Decays.size(), 0.);

    if (!eBW || m_Width / m_Mass < 0.01) {
      for (size_t i = 0; i < m_Decays.size(); ++i)
        ret[i] = m_Decays[i].mBratio;

      return ret;
    }

    double totwid = TotalWidtheBW(M);

    for (size_t i = 0; i < m_Decays.size(); ++i) {
      double partialwid = m_Decays[i].ModifiedWidth(M) * m_Width;
      ret[i] = partialwid / totwid;
    }

    return ret;
  }

  double ThermalParticle::ThermalMassDistribution(double M, double T, double Mu, double width)
  {
    return IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ParticleDensity, m_QuantumStatisticsCalculationType, m_Statistics, T, Mu, M, m_Degeneracy, m_ClusterExpansionOrder) * MassDistribution(M, width);
  }

  double ThermalParticle::ThermalMassDistribution(double M, double T, double Mu)
  {
    return ThermalMassDistribution(M, T, Mu, TotalWidtheBW(M));
  }

  void ThermalParticle::UseStatistics(bool enable) {
    if (!enable) m_Statistics = 0;
    else m_Statistics = m_StatisticsOrig;
  }

  void ThermalParticle::SetMass(double mass)
  {
    m_Mass = mass;
    if (m_Width != 0.0) {
      FillCoefficients();
      FillCoefficientsDynamical();
    }
  }

  int ThermalParticle::ConservedCharge(ConservedCharge::Name chg) const
  {
    if (chg == ConservedCharge::BaryonCharge)
      return BaryonCharge();
    if (chg == ConservedCharge::ElectricCharge)
      return ElectricCharge();
    if (chg == ConservedCharge::StrangenessCharge)
      return Strangeness();
    if (chg == ConservedCharge::CharmCharge)
      return Charm();
    return 0;
  }

  bool ThermalParticle::IsNeutral() const {
    return (m_Baryon == 0 && m_ElectricCharge == 0 && m_Strangeness == 0 && m_Charm == 0);
  }


  double ThermalParticle::Density(const ThermalModelParameters &params, IdealGasFunctions::Quantity type, bool useWidth, double mu) const {
    if (!(params.gammaq == 1.))                  mu += log(params.gammaq) * m_AbsQuark * params.T;
    if (!(params.gammaS == 1. || m_AbsS == 0.))  mu += log(params.gammaS) * m_AbsS     * params.T;
    if (!(params.gammaC == 1. || m_AbsC == 0.))  mu += log(params.gammaC) * m_AbsC     * params.T;

    //if (!useWidth || m_Mass == 0.0 || m_Width / m_Mass < 1.e-2 || m_ResonanceWidthIntegrationType == ZeroWidth) {
    if (!useWidth || m_Mass == 0.0 || ZeroWidthEnforced() || m_ResonanceWidthIntegrationType == ZeroWidth) {
      return IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, m_Statistics, params.T, mu, m_Mass, m_Degeneracy, m_ClusterExpansionOrder);
    }

    int ind = m_xleg.size();
    const vector<double> &x = m_xleg, &w = m_wleg;
    double ret1 = 0., ret2 = 0., tmp = 0.;

    // Integration from m0 or M-2*Gamma to M+2*Gamma
    if (m_ResonanceWidthIntegrationType != eBW && m_ResonanceWidthIntegrationType != eBWconstBR) {
      for (int i = 0; i < ind; i++) {

        tmp = w[i] * MassDistribution(x[i]);

        if (m_ResonanceWidthIntegrationType == FullIntervalWeighted)
          tmp *= m_brweight[i];

        double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, m_Statistics, params.T, mu, x[i], m_Degeneracy, m_ClusterExpansionOrder);

        ret1 += tmp * dens;
        ret2 += tmp;
      }
    }

    // Integration from M+2*Gamma to infinity
    if (m_ResonanceWidthIntegrationType == FullInterval || m_ResonanceWidthIntegrationType == FullIntervalWeighted) {
      int ind2 = m_xlag32.size();
      for (int i = 0; i < ind2; ++i) {
        double tmass = m_Mass + 2.*m_Width + m_xlag32[i] * m_Width;
        tmp = m_wlag32[i] * m_Width * MassDistribution(tmass);
        double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, m_Statistics, params.T, mu, tmass, m_Degeneracy, m_ClusterExpansionOrder);

        ret1 += tmp * dens;
        ret2 += tmp;
      }
    }

    if (m_ResonanceWidthIntegrationType == eBW || m_ResonanceWidthIntegrationType == eBWconstBR) {
      for (size_t i = 0; i < m_xalldyn.size(); i++) {
        tmp = m_walldyn[i];
        double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, m_Statistics, params.T, mu, m_xalldyn[i], m_Degeneracy, m_ClusterExpansionOrder);
        ret1 += tmp * dens;
        ret2 += tmp;
      }
    }

    return ret1 / ret2;
  }

  double ThermalParticle::DensityCluster(int n, const ThermalModelParameters & params, IdealGasFunctions::Quantity type, bool useWidth, double mu) const
  {
    double mn = 1.;
    if ((abs(BaryonCharge()) & 1) && !(n & 1))
      mn = -1.;

    if (!(params.gammaq == 1.))                  mu += log(params.gammaq) * m_AbsQuark /*GetAbsQ()*/  * params.T;
    if (!(params.gammaS == 1. || m_AbsS == 0.))  mu += log(params.gammaS) * m_AbsS     * params.T;
    if (!(params.gammaC == 1. || m_AbsC == 0.))  mu += log(params.gammaC) * m_AbsC     * params.T;

    //if (!useWidth || m_Width / m_Mass < 1.e-2 || m_ResonanceWidthIntegrationType == ZeroWidth) {
    if (!useWidth || ZeroWidthEnforced() || m_ResonanceWidthIntegrationType == ZeroWidth) {
      return mn * IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, 0, params.T / static_cast<double>(n), mu, m_Mass, m_Degeneracy);
    }

    int ind = m_xleg.size();
    const vector<double> &x = m_xleg, &w = m_wleg;
    double ret1 = 0., ret2 = 0., tmp = 0.;

    // Integration from m0 or M-2*Gamma to M+2*Gamma
    if (m_ResonanceWidthIntegrationType != eBW && m_ResonanceWidthIntegrationType != eBWconstBR) {
      for (int i = 0; i < ind; i++) {
        tmp = w[i] * MassDistribution(x[i]);

        if (m_ResonanceWidthIntegrationType == FullIntervalWeighted)
          tmp *= m_brweight[i];

        double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, 0, params.T / static_cast<double>(n), mu, x[i], m_Degeneracy);

        ret1 += tmp * dens;
        ret2 += tmp;
      }
    }

    // Integration from M+2*Gamma to infinity
    if (m_ResonanceWidthIntegrationType == FullInterval || m_ResonanceWidthIntegrationType == FullIntervalWeighted) {
      int ind2 = m_xlag32.size();
      for (int i = 0; i < ind2; ++i) {
        double tmass = m_Mass + 2.*m_Width + m_xlag32[i] * m_Width;
        tmp = m_wlag32[i] * m_Width * MassDistribution(tmass);
        double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, 0, params.T / static_cast<double>(n), mu, tmass, m_Degeneracy);

        ret1 += tmp * dens;
        ret2 += tmp;
      }
    }

    if (m_ResonanceWidthIntegrationType == eBW || m_ResonanceWidthIntegrationType == eBWconstBR) {
      for (size_t i = 0; i < m_xalldyn.size(); i++) {
        tmp = m_walldyn[i];
        double dens = IdealGasFunctions::IdealGasQuantity(type, m_QuantumStatisticsCalculationType, 0, params.T / static_cast<double>(n), mu, m_xalldyn[i], m_Degeneracy);
        ret1 += tmp * dens;
        ret2 += tmp;
      }
    }

    return mn * ret1 / ret2;
  }


  double ThermalParticle::chiDimensionfull(int index, const ThermalModelParameters& params, bool useWidth, double mu) const
  {
    if (index == 0) return Density(params, IdealGasFunctions::Pressure, useWidth, mu) / pow(xMath::GeVtoifm(), 3);
    if (index == 1) return Density(params, IdealGasFunctions::ParticleDensity, useWidth, mu) / pow(xMath::GeVtoifm(), 3);
    if (index == 2) return Density(params, IdealGasFunctions::chi2difull, useWidth, mu);
    if (index == 3) return Density(params, IdealGasFunctions::chi3difull, useWidth, mu);
    if (index == 4) return Density(params, IdealGasFunctions::chi4difull, useWidth, mu);
    return 1.;
  }

  double ThermalParticle::ScaledVariance(const ThermalModelParameters &params, bool useWidth, double mu) const {
    if (m_Degeneracy == 0.0) return 1.;
    if (m_Statistics == 0) return 1.;
    double dens = Density(params, IdealGasFunctions::ParticleDensity, useWidth, mu);
    if (dens == 0.) return 1.;
    double ret = params.T * chiDimensionfull(2, params, useWidth, mu) / chiDimensionfull(1, params, useWidth, mu);
    if (ret != ret) ret = 1.;
    return ret;
  }

  double ThermalParticle::Skewness(const ThermalModelParameters &params, bool useWidth, double mu) const
  {
    if (m_Degeneracy == 0) return 1.;
    if (m_Statistics == 0) return 1.;
    double dens = Density(params, IdealGasFunctions::ParticleDensity, useWidth, mu);
    if (dens == 0.) return 1.;
    double ret = params.T * chiDimensionfull(3, params, useWidth, mu) / chiDimensionfull(2, params, useWidth, mu);
    if (ret != ret) ret = 1.;
    return ret;
  }

  double ThermalParticle::Kurtosis(const ThermalModelParameters &params, bool useWidth, double mu) const
  {
    if (m_Degeneracy == 0) return 1.;
    if (m_Statistics == 0) return 1.;
    double dens = Density(params, IdealGasFunctions::ParticleDensity, useWidth, mu);
    if (dens == 0.) return 1.;
    double ret = params.T * params.T * chiDimensionfull(4, params, useWidth, mu) / chiDimensionfull(2, params, useWidth, mu);
    if (ret != ret) ret = 1.;
    return ret;
  }

  double ThermalParticle::FD(double k, double T, double mu, double m) const {
    double arg = (sqrt(k*k + m * m) - mu) / T;
    return 1. / (exp(arg) + 1.);
  }

  double ThermalParticle::GetAbsQ() const {
    if (m_Baryon == 0) return 2. - m_AbsC - m_AbsS;
    else return abs(m_Baryon) * 3. - m_AbsC - m_AbsS;
  }

  double ThermalParticle::GetCharge(int index) const {
    if (index == 0) return (double)m_Baryon;
    if (index == 1) return (double)m_ElectricCharge;
    if (index == 2) return (double)m_Strangeness;
    if (index == 3) return (double)m_Charm;
    return 0.;
  }

  double ThermalParticle::GetAbsCharge(int index) const {
    if (index == 0) return fabs((double)m_Baryon);
    if (index == 1) return fabs((double)m_ElectricCharge);
    if (index == 2) return m_AbsS;
    if (index == 3) return m_AbsC;
    return 0.;
  }

  double ThermalParticle::chi(int index, const ThermalModelParameters &params, bool useWidth, double mu) const {
    if (index == 0) return Density(params, IdealGasFunctions::Pressure, useWidth, mu) / pow(params.T, 4) / pow(xMath::GeVtoifm(), 3);
    if (index == 1) return Density(params, IdealGasFunctions::ParticleDensity, useWidth, mu) / pow(params.T, 3) / pow(xMath::GeVtoifm(), 3);
    if (index == 2) return Density(params, IdealGasFunctions::chi2, useWidth, mu);
    if (index == 3) return Density(params, IdealGasFunctions::chi3, useWidth, mu);
    if (index == 4) return Density(params, IdealGasFunctions::chi4, useWidth, mu);
    return 1.;
  }

} // namespace thermalfist