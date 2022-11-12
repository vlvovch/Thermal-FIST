#include "HRGRealGas/MeanFieldModels.h"


namespace thermalfist {
  double MeanFieldModelVDF::v(double n) const
  {
    double ret = 0.;
    for (int k = 1; k < m_N; ++k) {
      ret += m_Ck[k] / m_bk[k] * pow(n / m_n0, m_bk[k] - 1);
    }
    ret *= n;
    return ret;
  }

  double MeanFieldModelVDF::dv(int order, double n) const
  {
    if (order == 0)
      return v(n);

    double ret = 0.;
    for (int k = 1; k < m_N; ++k) {
      double tret = 1.;
      for (int i = 1; i < order; ++i)
        tret *= (m_bk[k] - i + 1);
      ret += m_Ck[k] / m_bk[k] * tret * pow(n / m_n0, m_bk[k] - 1 - order) / pow(m_n0, order);
    }
    ret *= n;
    return ret;
  }

  double MeanFieldModelVDF::dvdT(double n) const
  {
    double ret = 0.;
    for (int k = 1; k < m_N; ++k) {
      ret += m_dCkdT[k] / m_bk[k] * pow(n / m_n0, m_bk[k] - 1);
      ret += -m_Ck[k] / m_bk[k] / m_bk[k] * m_dbkdT[k] * pow(n / m_n0, m_bk[k] - 1);
      ret += m_Ck[k] / m_bk[k] * pow(n / m_n0, m_bk[k] - 1) * log(n / m_n0) * m_dbkdT[k];
    }
    ret *= n;
    return ret;
  }
} // namespace thermalfist