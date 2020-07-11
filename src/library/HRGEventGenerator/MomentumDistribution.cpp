/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEventGenerator/MomentumDistribution.h"

#include <vector>
#include <iostream>

#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"

namespace thermalfist {

  void SiemensRasmussenDistribution::Normalize() {
    m_Norm = 1.;
    double tmp = 0.;

    NumericalIntegration::GetCoefsIntegrateLaguerre32(&m_xlag, &m_wlag);
    for (int i = 0; i < 32; i++) {
      tmp += m_wlag[i] * dndp(m_xlag[i] * m_T);
    }
    tmp *= m_T;
    m_Norm = 1. / tmp;

    m_Normalized = true;
  }

  double SiemensRasmussenDistribution::dndp(double p) const {
    double alphap = alpha(p);
    double alphamn = 1.;
    if (alphap != 0.0)
      alphamn = sinh(alpha(p)) / alpha(p);
    return m_Norm * p * p * exp(-m_Gamma * w(p) / m_T) * ((1. + m_T / m_Gamma / w(p))*alphamn - m_T / m_Gamma / w(p)*cosh(alpha(p)));
  }

  double SiemensRasmussenDistribution::dndy(double y) const {
    std::vector<double> x, w;
    NumericalIntegration::GetCoefsIntegrateLaguerre32(&x, &w);
    double ret = 0.;
    for (int i = 0; i < 32; i++) {
      double tpt = x[i] * m_T;
      double tmt = sqrt(m_Mass*m_Mass + tpt * tpt);
      double tpz = tmt * sinh(y);
      double tp = sqrt(tpt*tpt + tpz * tpz);
      double ten = tmt * cosh(y);

      double alphap = alpha(tp);
      double alphamn = 1.;
      if (alphap != 0.0)
        alphamn = sinh(alpha(tp)) / alpha(tp);

      double tadd = w[i] * ten * tpt * exp(-m_Gamma * ten / m_T) * ((1. + m_T / m_Gamma / ten)*alphamn - m_T / m_Gamma / ten * cosh(alpha(tp)));
      if (tadd != tadd) break;
      ret += tadd;
    }
    return m_T * m_Norm / 2. * ret;
  }

  double SiemensRasmussenDistribution::dnmtdmt(double mt) const {
    std::vector<double> x, w;
    NumericalIntegration::GetCoefsIntegrateLegendre40(-4., 4., &x, &w);
    double ret = 0.;
    for (int i = 0; i < 32; i++) {
      double ty = x[i];
      double tmt = mt;
      double tpt = sqrt(tmt*tmt - m_Mass * m_Mass);
      double tpz = tmt * sinh(ty);
      double tp = sqrt(tpt*tpt + tpz * tpz);
      double ten = tmt * cosh(ty);

      double alphap = alpha(tp);
      double alphamn = 1.;
      if (alphap != 0.0)
        alphamn = sinh(alpha(tp)) / alpha(tp);

      double tadd = w[i] * ten * exp(-m_Gamma * ten / m_T) * ((1. + m_T / m_Gamma / ten)*alphamn - m_T / m_Gamma / ten * cosh(alpha(tp)));
      if (tadd != tadd) break;
      ret += tadd;
    }
    return m_Norm / 2. * ret;
  }

  double SiemensRasmussenDistribution::PAv() const {
    if (!m_useacc) {
      double tmp1 = 0., tmp2 = 0.;
      for (int i = 0; i < 32; i++) {
        double tmp = m_wlag[i] * dndp(m_xlag[i] * m_T);
        tmp1 += tmp * m_xlag[i] * m_T;
        tmp2 += tmp;
      }
      return tmp1 / tmp2;
    }
    else {
      std::vector<double> xy, wy;
      NumericalIntegration::GetCoefsIntegrateLegendre40(-4., 4., &xy, &wy);
      std::vector<double> xpt, wpt;
      NumericalIntegration::GetCoefsIntegrateLaguerre32(&xpt, &wpt);
      double tmp1 = 0., tmp2 = 0.;
      for (size_t iy = 0; iy < xy.size(); ++iy) {
        double ty = xy[iy];
        for (size_t ipt = 0; ipt < xpt.size(); ++ipt) {
          double tpt = xpt[ipt] * m_T;
          double tmp = wy[iy] * wpt[ipt] * d2ndptdy(tpt, ty);
          if (tmp <= 0. || tmp != tmp) continue;
          double tmt = sqrt(tpt*tpt + m_Mass * m_Mass);
          double tpz = tmt * sinh(ty);
          double tp = sqrt(tpt*tpt + tpz * tpz);
          tmp1 += tmp * tp;
          tmp2 += tmp;
        }
      }
      return tmp1 / tmp2;
    }
  }


  double SiemensRasmussenDistribution::d2ndptdy(double pt, double y) const {
    double tpt = pt, ty = y;
    double tmt = sqrt(m_Mass*m_Mass + tpt * tpt);
    double tpz = tmt * sinh(ty);
    double tp = sqrt(tpt*tpt + tpz * tpz);
    double ten = tmt * cosh(ty);

    double alphap = alpha(tp);
    double alphamn = 1.;
    if (alphap != 0.0)
      alphamn = sinh(alpha(tp)) / alpha(tp);

    double ret = m_Norm / 2. * ten * tpt * exp(-m_Gamma * ten / m_T) * ((1. + m_T / m_Gamma / ten)*alphamn - m_T / m_Gamma / ten * cosh(alpha(tp)));
    if (!m_useacc) return ret;
    else return ret * m_acc->getAcceptance(y + m_ycm, pt);
  }



  BoostInvariantMomentumDistribution::~BoostInvariantMomentumDistribution()
  {
    if (m_FreezeoutModel != NULL) {
      delete m_FreezeoutModel;
    }
  }

  void BoostInvariantMomentumDistribution::Normalize()
  {
    m_Norm = m_NormY = m_NormPt = 1.;
    double tmp = 0.;

    Initialize();

    tmp = 0.;
    for (size_t i = 0; i < m_xlag.size(); i++) {
      tmp += m_wlag[i] * dndpt(m_xlag[i] * m_T);
    }
    tmp *= m_T;
    m_NormPt = 1. / tmp;


    m_dndy.clearall();
    m_dndy.add_val(-100., 0.);
    m_dndyint.clearall();
    m_dndyint.add_val(-100., 0.);
    double dy = 0.05;
    double tmp2 = 0.;
    for (double yy = -4.; yy < 4.; yy += dy) {
      tmp = 0.;
      for (size_t j = 0; j < m_xlag.size(); j++) {
        tmp += m_wlag[j] * dndptsingle(m_xlag[j] * m_T, yy) * m_T;
      }
      tmp2 += tmp * dy;
      m_dndy.add_val(yy, tmp);
      m_dndyint.add_val(yy + 0.5 * dy, tmp2);
    }
    m_dndy.add_val(100., 0.);
    m_dndyint.add_val(100., tmp2);
    tmp = 0.;
    m_NormY = m_Norm = 1. / tmp2;

    if (m_xlegeta.size() > 1) {
      m_NormY *= 1. / 2. / m_EtaMax;
      m_Norm *= 1. / 2. / m_EtaMax;
    }

    m_Normalized = true;
  }

  double BoostInvariantMomentumDistribution::ZetaIntegrandpT(double zeta, double pt) const
  {
    double R = m_FreezeoutModel->Rfunc(zeta);
    double tau = m_FreezeoutModel->taufunc(zeta);
    double dRdZeta = m_FreezeoutModel->dRdZeta(zeta);
    double dTaudZeta = m_FreezeoutModel->dtaudZeta(zeta);
    double mT = sqrt(m_Mass * m_Mass + pt * pt);
    double coshetaperp = m_FreezeoutModel->coshetaperp(zeta);
    double sinhetaperp = m_FreezeoutModel->sinhetaperp(zeta);

    return R * tau * (mT * dRdZeta * xMath::BesselK1(mT * coshetaperp / m_T) * xMath::BesselI0(pt * sinhetaperp / m_T)
      - pt * dTaudZeta * xMath::BesselK0(mT * coshetaperp / m_T) * xMath::BesselI1(pt * sinhetaperp / m_T));
  }

  double BoostInvariantMomentumDistribution::dndy(double y) const
  {
    if (m_xlegeta.size() > 1)
      return (m_dndyint.f(m_EtaMax + y) - m_dndyint.f(-m_EtaMax + y)) * m_NormY;
    else
      return m_dndy.f(y) * m_NormY;
  }

  double BoostInvariantMomentumDistribution::dnmtdmt(double mt) const
  {
    double ret = 0.;
    double pt = sqrt(mt * mt - m_Mass * m_Mass);
    for (size_t i = 0; i < m_xlegT.size(); i++) {
      double tmp = m_wlegT[i] * ZetaIntegrandpT(m_xlegT[i], pt);
      if (tmp != tmp) break;
      ret += tmp;
    }
    return ret * m_NormPt;
  }

  double BoostInvariantMomentumDistribution::d2ndptdy(double pt, double y) const
  {
    double ret = 0.;
    double mt = sqrt(pt * pt + m_Mass * m_Mass);

    for (size_t i = 0; i < m_xlegeta.size(); i++) {
      for (size_t j = 0; j < m_xlegT.size(); j++) {
        double tmp = m_wlegeta[i] * m_wlegT[j] * ZetaIntegrandpTYSingleFireball(m_xlegeta[i], pt, y - m_xlegeta[i]);
        ret += tmp;
      }
    }

    if (!m_useacc) return ret * m_Norm * pt;
    else return ret * m_Norm * pt * m_acc->getAcceptance(y + m_ycm, pt);
  }

  double BoostInvariantMomentumDistribution::ZetaIntegrandpTYSingleFireball(double zeta, double pt, double y) const
  {
    double R = m_FreezeoutModel->Rfunc(zeta);
    double tau = m_FreezeoutModel->taufunc(zeta);
    double dRdZeta = m_FreezeoutModel->dRdZeta(zeta);
    double dTaudZeta = m_FreezeoutModel->dtaudZeta(zeta);
    double mT = sqrt(m_Mass * m_Mass + pt * pt);
    double coshetaperp = m_FreezeoutModel->coshetaperp(zeta);
    double sinhetaperp = m_FreezeoutModel->sinhetaperp(zeta);

    return R * tau * (mT * dRdZeta * cosh(y) * xMath::BesselI0(pt * sinhetaperp / m_T)
      - pt * dTaudZeta * xMath::BesselI1(pt * sinhetaperp / m_T)) * exp(-mT * coshetaperp * cosh(y) / m_T);
  }

  void BoostInvariantMomentumDistribution::Initialize()
  {
    NumericalIntegration::GetCoefsIntegrateLaguerre32(&m_xlag, &m_wlag);
    NumericalIntegration::GetCoefsIntegrateLegendre32(0., 1., &m_xlegT, &m_wlegT);
    NumericalIntegration::GetCoefsIntegrateLegendre32(-4., 4., &m_xlegY, &m_wlegY);

    if (m_EtaMax > 0.0)
      NumericalIntegration::GetCoefsIntegrateLegendre32(-m_EtaMax, m_EtaMax, &m_xlegeta, &m_wlegeta);
    else {
      m_xlegeta.resize(1);
      m_xlegeta[0] = 0.;
      m_wlegeta.resize(1);
      m_wlegeta[0] = 1.;
    }
  }

  double BoostInvariantMomentumDistribution::dndysingle(double y) const
  {
    return m_dndy.f(y);
  }

  double BoostInvariantMomentumDistribution::dndptsingle(double pt, double y) const
  {
    double ret = 0.;
    double mt = sqrt(pt * pt + m_Mass * m_Mass);

    for (size_t j = 0; j < m_xlegT.size(); j++) {
      double tmp = m_wlegT[j] * ZetaIntegrandpTYSingleFireball(m_xlegT[j], pt, y);
      ret += tmp;
    }

    return ret * pt;
  }

  double BoostInvariantMomentumDistribution::dndpt(double pt, double y) const
  {
    double ret = 0.;
    double mt = sqrt(pt * pt + m_Mass * m_Mass);

    for (size_t i = 0; i < m_xlegeta.size(); i++) {
      ret += m_wlegeta[i] * dndptsingle(pt, y - m_xlegeta[i]);
    }

    return ret * m_Norm * pt;
  }

} // namespace thermalfist


/////////////////////////////////////////////////////////////////////////////
// Deprecated code below

namespace thermalfist {
  void SSHDistribution::Initialize() {
    NumericalIntegration::GetCoefsIntegrateLaguerre32(&m_xlag, &m_wlag);
    NumericalIntegration::GetCoefsIntegrateLegendre32(0., 1., &m_xlegT, &m_wlegT);
    NumericalIntegration::GetCoefsIntegrateLegendre32(-4., 4., &m_xlegY, &m_wlegY);

    if (m_EtaMax > 0.0)
      NumericalIntegration::GetCoefsIntegrateLegendre32(-m_EtaMax, m_EtaMax, &m_xlegeta, &m_wlegeta);
    else {
      m_xlegeta.resize(1);
      m_xlegeta[0] = 0.;
      m_wlegeta.resize(1);
      m_wlegeta[0] = 1.;
    }
  }


  void SSHDistribution::Normalize() {
    m_Norm = m_NormY = m_NormPt = 1.;
    double tmp = 0.;

    NumericalIntegration::GetCoefsIntegrateLaguerre32(&m_xlag, &m_wlag);
    NumericalIntegration::GetCoefsIntegrateLegendre32(0., 1., &m_xlegT, &m_wlegT);
    NumericalIntegration::GetCoefsIntegrateLegendre32(-4., 4., &m_xlegY, &m_wlegY);

    if (m_EtaMax > 0.0)
      NumericalIntegration::GetCoefsIntegrateLegendre32(-m_EtaMax, m_EtaMax, &m_xlegeta, &m_wlegeta);
    else {
      m_xlegeta.resize(1);
      m_xlegeta[0] = 0.;
      m_wlegeta.resize(1);
      m_wlegeta[0] = 1.;
    }

    tmp = 0.;
    for (size_t i = 0; i < m_xlag.size(); i++) {
      tmp += m_wlag[i] * dndpt(m_xlag[i] * m_T);
    }
    tmp *= m_T;
    m_NormPt = 1. / tmp;


    m_dndy.clearall();
    m_dndy.add_val(-100., 0.);
    m_dndyint.clearall();
    m_dndyint.add_val(-100., 0.);
    double dy = 0.05;
    double tmp2 = 0.;
    for (double yy = -4.; yy < 4.; yy += dy) {
      tmp = 0.;
      for (size_t j = 0; j < m_xlag.size(); j++) {
        tmp += m_wlag[j] * dndptsingle(m_xlag[j] * m_T, yy) * m_T;
      }
      tmp2 += tmp * dy;
      m_dndy.add_val(yy, tmp);
      m_dndyint.add_val(yy + 0.5 * dy, tmp2);
    }
    m_dndy.add_val(100., 0.);
    m_dndyint.add_val(100., tmp2);
    tmp = 0.;
    m_NormY = m_Norm = 1. / tmp2;

    if (m_xlegeta.size() > 1) {
      m_NormY *= 1. / 2. / m_EtaMax;
      m_Norm *= 1. / 2. / m_EtaMax;
    }

    m_Normalized = true;
  }



  double SSHDistribution::dndy(double y) const {
    if (m_xlegeta.size() > 1)
      return (m_dndyint.f(m_EtaMax + y) - m_dndyint.f(-m_EtaMax + y)) * m_NormY;
    else
      return m_dndy.f(y) * m_NormY;
  }

  double SSHDistribution::dndysingle(double y) const {
    return m_dndy.f(y);
  }

  double SSHDistribution::dnmtdmt(double mt) const {
    double ret = 0.;
    double pt = sqrt(mt * mt - m_Mass * m_Mass);
    for (size_t i = 0; i < m_xlegT.size(); i++) {
      double tmp = m_wlegT[i] * m_xlegT[i] * mt * xMath::BesselI0(pt * sinh(rho(m_xlegT[i])) / m_T) * xMath::BesselK1(mt * cosh(rho(m_xlegT[i])) / m_T);
      if (tmp != tmp) break;
      ret += tmp;
    }
    return ret * m_NormPt;
  }

  double SSHDistribution::dndy(double y, double pt) const {
    double ret = 0.;
    double mt = sqrt(pt * pt + m_Mass * m_Mass);

    for (size_t i = 0; i < m_xlegeta.size(); i++) {
      for (size_t j = 0; j < m_xlegT.size(); j++) {
        double tmp = mt * m_wlegeta[i] * m_wlegT[j] * cosh(y - m_xlegeta[i]) * m_xlegT[j] * exp(-mt * cosh(rho(m_xlegT[j])) * cosh(y - m_xlegeta[i]) / m_T) *
          xMath::BesselI0(pt * sinh(rho(m_xlegT[j])) / m_T);
        ret += tmp;
      }
    }

    return ret * m_Norm;
  }

  double SSHDistribution::dndysingle(double y, double pt) const {
    double ret = 0.;
    double mt = sqrt(pt * pt + m_Mass * m_Mass);

    for (size_t j = 0; j < m_xlegT.size(); j++) {
      double tmp = mt * m_wlegT[j] * cosh(y) * m_xlegT[j] * exp(-mt * cosh(rho(m_xlegT[j])) * cosh(y) / m_T) *
        xMath::BesselI0(pt * sinh(rho(m_xlegT[j])) / m_T);
      ret += tmp;
    }

    return ret;
  }

  double SSHDistribution::dndpt(double pt) const {
    return pt * dnmtdmt(sqrt(pt * pt + m_Mass * m_Mass));
  }

  double SSHDistribution::dndpt(double pt, double y) const {
    double ret = 0.;
    double mt = sqrt(pt * pt + m_Mass * m_Mass);

    for (size_t i = 0; i < m_xlegeta.size(); i++) {
      for (size_t j = 0; j < m_xlegT.size(); j++) {
        double tmp = mt * m_wlegeta[i] * m_wlegT[j] * cosh(y - m_xlegeta[i]) * m_xlegT[j] * exp(-mt * cosh(rho(m_xlegT[j])) * cosh(y - m_xlegeta[i]) / m_T) *
          xMath::BesselI0(pt * sinh(rho(m_xlegT[j])) / m_T);
        ret += tmp;
      }
    }

    return ret * m_Norm * pt;
  }

  double SSHDistribution::dndptsingle(double pt, double y) const {
    double ret = 0.;
    double mt = sqrt(pt * pt + m_Mass * m_Mass);

    for (size_t j = 0; j < m_xlegT.size(); j++) {
      double tmp = mt * m_wlegT[j] * cosh(y) * m_xlegT[j] * exp(-mt * cosh(rho(m_xlegT[j])) * cosh(y) / m_T) *
        xMath::BesselI0(pt * sinh(rho(m_xlegT[j])) / m_T);
      ret += tmp;
    }

    return ret * pt;
  }

  double SSHDistribution::MtAv() const {
    if (0 && !m_useacc) {
      double tmp1 = 0., tmp2 = 0.;
      for (int i = 0; i < 32; i++) {
        double tmppt = m_xlag[i] * m_T;
        double tmpmt = sqrt(tmppt * tmppt + m_Mass * m_Mass);
        double tmp = m_wlag[i] * tmppt * dnmtdmt(tmpmt);
        tmp1 += tmp * tmpmt;
        tmp2 += tmp;
      }
      return tmp1 / tmp2;
    }
    else {
      std::vector<double> xy, wy;
      NumericalIntegration::GetCoefsIntegrateLegendre40(-4., 4., &xy, &wy);
      std::vector<double> xpt, wpt;
      NumericalIntegration::GetCoefsIntegrateLaguerre32(&xpt, &wpt);
      double tmp1 = 0., tmp2 = 0.;
      for (size_t iy = 0; iy < xy.size(); ++iy) {
        double ty = xy[iy];
        for (size_t ipt = 0; ipt < xpt.size(); ++ipt) {
          double tpt = xpt[ipt] * m_T;
          double tmp = wy[iy] * wpt[ipt] * d2ndptdy(tpt, ty);
          if (tmp <= 0. || tmp != tmp) continue;
          double tmt = sqrt(tpt * tpt + m_Mass * m_Mass);
          tmp1 += tmp * tmt;
          tmp2 += tmp;
        }
      }
      return tmp1 / tmp2;
    }
  }

  double SSHDistribution::y2Av() const {
    if (0 && !m_useacc) {
      double tmp1 = 0., tmp2 = 0.;
      double dy = 0.05;
      for (double yy = 0.; yy < 4.; yy += dy) {
        double tmp = dndysingle(yy);
        tmp1 += tmp * yy * yy;
        tmp2 += tmp;
      }
      return (tmp1 / tmp2 + m_EtaMax * m_EtaMax / 3.);
    }
    else {
      std::vector<double> xy, wy;
      NumericalIntegration::GetCoefsIntegrateLegendre40(-4., 4., &xy, &wy);
      std::vector<double> xpt, wpt;
      NumericalIntegration::GetCoefsIntegrateLaguerre32(&xpt, &wpt);
      double tmp1 = 0., tmp2 = 0.;
      for (size_t iy = 0; iy < xy.size(); ++iy) {
        double ty = xy[iy];
        for (size_t ipt = 0; ipt < xpt.size(); ++ipt) {
          double tpt = xpt[ipt] * m_T;
          double tmp = wy[iy] * wpt[ipt] * d2ndptdy(tpt, ty);
          if (tmp <= 0. || tmp != tmp) continue;
          //double tmt = sqrt(tpt*tpt + m_Mass * m_Mass);
          tmp1 += tmp * ty * ty;
          tmp2 += tmp;
        }
      }
      return tmp1 / tmp2;
    }
  }

  double SSHDistribution::d2ndptdy(double pt, double y) const {
    double ret = 0.;
    double mt = sqrt(pt * pt + m_Mass * m_Mass);

    for (size_t i = 0; i < m_xlegeta.size(); i++) {
      for (size_t j = 0; j < m_xlegT.size(); j++) {
        double tmp = mt * m_wlegeta[i] * m_wlegT[j] * cosh(y - m_xlegeta[i]) * m_xlegT[j] * exp(-mt * cosh(rho(m_xlegT[j])) * cosh(y - m_xlegeta[i]) / m_T) *
          xMath::BesselI0(pt * sinh(rho(m_xlegT[j])) / m_T);
        ret += tmp;
      }
    }

    if (!m_useacc) return ret * m_Norm * pt;
    else return ret * m_Norm * pt * m_acc->getAcceptance(y + m_ycm, pt);
  }
} // namespace thermalfist