/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "NumberStatistics.h"


#include <cmath>


void NumberStatistics::AddEvent(int tmpn, double weight) {
    double tmn = static_cast<double>(tmpn);
    AddEvent(tmn, weight);
}

void NumberStatistics::AddEvent(double tmpn, double weight) {
    double tmn = tmpn;
    n   += weight*tmn;
    n2  += weight*tmn*tmn;
    n3  += weight*tmn*tmn*tmn;
    n4  += weight*tmn*tmn*tmn*tmn;
    n5  += weight*tmn*tmn*tmn*tmn*tmn;
    n6  += weight*tmn*tmn*tmn*tmn*tmn*tmn;
    n7  += weight*tmn*tmn*tmn*tmn*tmn*tmn*tmn;
    n8  += weight*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn;
    n9  += weight*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn;
    n10 += weight*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn;
    n11 += weight*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn;
    n12 += weight*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn*tmn;
    wsum  += weight;
    w2sum += weight*weight;
    events++;
}

double NumberStatistics::GetMean() const {
    return means[1];
}

double NumberStatistics::GetMeanError() const {
    return sqrt((means[2]-means[1]*means[1])/(nE-1.));
}

double NumberStatistics::GetVariance() const {
    return (means[2]-means[1]*means[1]);
}

double NumberStatistics::GetStdDev() const {
    return sqrt(means[2]-means[1]*means[1]);
}

double NumberStatistics::GetScaledVariance() const {
    return GetVariance() / GetMean();
}

double NumberStatistics::GetSkewness() const {
    return (means[3]-3.*means[2]*means[1] + 2. * means[1] * means[1] * means[1]) / GetVariance();
}

double NumberStatistics::GetSkewnessError() const {
  double sig = GetStdDev();
  std::vector<double> sigv(7, 0.);
  sigv[0] = 1.;
  for(int i=1;i<sigv.size();++i)
    sigv[i] = sigv[i-1] * sig;
  return sqrt((9. - 6. * m[4] / sigv[4] + m[3] * m[3] / sigv[3] / sigv[3] * (6. + m[4] / sigv[4]) - 2. * m[3] / sigv[3] * m[5] / sigv[5] + m[6] / sigv[6]) * sigv[2]) / sqrt(nE);
}

double NumberStatistics::GetSkewnessError2() const {
  return sqrt(6.*events*(events-1.)/(events-2.)/(events+1.)/(events+3.) * GetVariance() + 0.*GetSkewness() * GetSkewness() * GetStdDevError() * GetStdDevError() / GetVariance());
}

double NumberStatistics::GetKurtosis() const {
    return ( means[4] - 4. * means[3] * means[1] + 6. * means[2]*means[1]*means[1] - 3. * means[1]*means[1]*means[1]*means[1] )/GetVariance() - 3. * GetVariance();
}

double NumberStatistics::GetKurtosisError() const {
    double sig = GetStdDev();
  std::vector<double> sigv(9, 0.);
  sigv[0] = 1.;
  for(int i=1;i<sigv.size();++i)
    sigv[i] = sigv[i-1] * sig;
  std::vector<double> mn(9, 0.);
  mn[0] = 1.;
  for(int i=1;i<sigv.size();++i)
    mn[i] = m[i] / sigv[i];
  return sqrt((-9.+6.*mn[4]*mn[4]+mn[4]*mn[4]*mn[4]+8.*mn[3]*mn[3]*(5.+mn[4])-8.*mn[3]*mn[5]+mn[4]*(9.-2.*mn[6])-6.*mn[6]+mn[8])*sigv[4]) / sqrt(nE);
}

double NumberStatistics::GetKurtosisError2() const {
  return sqrt(24.*events*(events-1.)*(events-1.)/(events-3.)/(events-2.)/(events+3.)/(events+5.)
    * GetVariance() * GetVariance() + 0.*GetKurtosis() * GetKurtosis() * GetVarianceError() * GetVarianceError() / GetVariance() / GetVariance());
}

double NumberStatistics::GetN2Error2() const {
    return 1. / (nE-1.) * (means[4] - means[2]*means[2]);
}

double NumberStatistics::GetVarianceError() const {

    double nav  = n  / events;
    double n2av = n2 / events;
    double n3av = n3 / events;
    double n4av = n4 / events;
    return sqrt(m[4] - m[2] * m[2]) / sqrt(nE);

    double tmp = n4av-4.*n3av*nav+8.*n2av*nav*nav-4.*nav*nav*nav*nav-n2av*n2av;

    return sqrt(n4av-4.*n3av*nav+8.*n2av*nav*nav-4.*nav*nav*nav*nav-n2av*n2av)/sqrt(events);
}

double NumberStatistics::GetScaledVarianceError() const {
    return GetVarianceError() / GetMean();
}

double NumberStatistics::GetStdDevError() const {
  return GetVarianceError() / 2. / GetStdDev();
}

void NumberStatistics::CalculateCentralMoments() {
    means.resize(13);
    means[0]  = 1.;
    means[1]  = n   / wsum;
    means[2]  = n2  / wsum;
    means[3]  = n3  / wsum;
    means[4]  = n4  / wsum;
    means[5]  = n5  / wsum;
    means[6]  = n6  / wsum;
    means[7]  = n7  / wsum;
    means[8]  = n8  / wsum;
    means[9]  = n9  / wsum;
    means[10] = n10 / wsum;
    means[11] = n11 / wsum;
    means[12] = n12 / wsum;
  m[1] = 0.;
  for(int mm=2;mm<=12;++mm) {
    m[mm] = 0.;
    double tmpn = 1;
    for(int i=0;i<=mm;++i) {
            m[mm] += Cnm(mm, i) * means[mm-i] * tmpn * ((i & 1) ? -1. : 1.);
            tmpn *= means[1];
    }
  }

    nE = events * wsum * wsum / w2sum / events;
}

double NumberStatistics::GetC5C2() const {
  return m[5] / m[2] - 10. * m[3];
}

double NumberStatistics::GetC5C2Error() const {
  double sig = GetStdDev();
  std::vector<double> sigv(11, 0.);
  sigv[0] = 1.;
  for(int i=1;i<sigv.size();++i)
    sigv[i] = sigv[i-1] * sig;
  std::vector<double> mn(11, 0.);
  mn[0] = 1.;
  for(int i=1;i<sigv.size();++i)
    mn[i] = m[i] / sigv[i];
    return sqrt((mn[10] - 100. * mn[3] * mn[3] + 10. * mn[3] * (-6. + mn[4]) * mn[5]
            + mn[4] * (125.*mn[4] + mn[5]*mn[5] - 10.*(90.-mn[6])) - 2.*mn[5]*mn[7]
            + 20.*(45.+mn[5]*mn[5]+8.*mn[6]-mn[8])) * sigv[6]) / sqrt(nE);
}

double NumberStatistics::GetC5C2Error2() const {
  double sig2 = GetVariance();
    return sqrt(120. * sig2 * sig2 * sig2 / nE);
}

double NumberStatistics::GetC6C2() const {
  return (m[6] - 15.*m[2]*m[4] - 10.*m[3]*m[3] + 30.*m[2]*m[2]*m[2])/m[2];
}

double NumberStatistics::GetC6C2Error() const {
  double sig = GetStdDev();
  std::vector<double> sigv(13, 0.);
  sigv[0] = 1.;
  for(int i=1;i<sigv.size();++i)
    sigv[i] = sigv[i-1] * sig;
  std::vector<double> mn(13, 0.);
  mn[0] = 1.;
  for(int i=1;i<sigv.size();++i)
    mn[i] = m[i] / sigv[i];
  
  return sqrt((10575.-30.*mn[10]+mn[12]+18300.*mn[3]*mn[3]
        +2600.*mn[3]*mn[3]*mn[3]*mn[3]-225.*(-3.+mn[4])*(-3.+mn[4])-7440*mn[3]*mn[5]
        -520*mn[3]*mn[3]*mn[3]*mn[5]+216.*mn[5]*mn[5]-2160.*mn[6]-200.*mn[3]*mn[3]*mn[6]+52.*mn[3]*mn[5]*mn[6]+33.*mn[6]*mn[6]
        + (-3.+mn[4])*(10.*(405.-390.*mn[3]*mn[3]+10.*mn[3]*mn[3]*mn[3]*mn[3]+24.*mn[3]*mn[5])-20.*(6.+mn[3]*mn[3])*mn[6]+mn[6]*mn[6])
            + 840.*mn[3]*mn[7]-12.*mn[5]*mn[7]+345.*mn[8]+20.*mn[3]*mn[3]*mn[8]-2.*mn[6]*mn[8]-40.*mn[3]*mn[9]) * sigv[8]) / sqrt(nE);
}

double NumberStatistics::GetC6C2Error2() const {
  double sig2 = GetVariance();
  return sqrt(720. * sig2 * sig2 * sig2 * sig2 / nE);
}

double NumberStatistics::GetC6C2alt() const {
  return m[6] / m[2] - 15.*m[4] - 10.*m[3]*m[3] / m[2] + 30.*m[2]*m[2];
}

double NumberStatistics::GetC6C2Erroralt() const {
  double sig = GetStdDev();
  
  double ret = 0.;
  ret += m[12]*m[2]*m[2];
  ret += -30. * m[10] * m[2] * m[2] * m[2];
  ret += -3600. * m[2] * m[2] * m[2] * m[2] * m[2] * m[2] * m[2] * m[2];
  ret += 30000. * m[2] * m[2] * m[2] * m[2] * m[2] * m[3] * m[3];
  ret += 2300. * m[2] * m[2] * m[3] * m[3] * m[3] * m[3];
  ret += 5400. * m[2] * m[2] * m[2] * m[2] * m[2] * m[2] * m[4];
  ret += -3900. * m[2] * m[2] * m[2] * m[3] * m[3] * m[4];
  ret += 100. * m[3] * m[3] * m[3] * m[3] * m[4];
  ret += -225. * m[2] * m[2] * m[2] * m[2] * m[4] * m[4];
  ret += -8160. * m[2] * m[2] * m[2] * m[2] * m[3] * m[5];
  ret += -520. * m[2] * m[3] * m[3] * m[3] * m[5];
  ret += 240. * m[2] * m[2] * m[3] * m[4] * m[5];
  ret += 216. * m[2] * m[2] * m[2] * m[5] * m[5];
  ret += -1800. * m[2] * m[2] * m[2] * m[2] * m[2] * m[6];
  ret += -140. * m[2] * m[2] * m[3] * m[3] * m[6];
  ret += -120. * m[2] * m[2] * m[2] * m[4] * m[6];
  ret += -20. * m[3] * m[3] * m[4] * m[6];
  ret += 52. * m[2] * m[3] * m[5] * m[6];
  ret += 30. * m[2] * m[2] * m[6] * m[6];
  ret += m[4] * m[6] * m[6];
  ret += 840. * m[2] * m[2] * m[2] * m[3] * m[7];
  ret += -12. * m[2] * m[2] * m[5] * m[7];
  ret += 345. * m[2] * m[2] * m[2] * m[2] * m[8];
  ret += 20. * m[2] * m[3] * m[3] * m[8];
  ret += -2. * m[2] * m[6] * m[8];
  ret += -40. * m[2] * m[2] * m[3] * m[9];
  ret /= m[2] * m[2] * m[2] * m[2];
  ret /= events;

  return sqrt(ret);
}
