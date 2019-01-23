/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <string>
#include <vector>
#define BINOMSZ 17

class NumberStatistics
{
  std::vector<double> m;
    long long events;
    std::string name;
    double mass;
    bool acc;
  std::vector< std::vector<int> > binom;
public:
  double n, n2, n3, n4;
    double n5, n6, n7, n8;
  double n9, n10, n11, n12;
    std::vector<double> means;
    double wsum;
    double w2sum;
    double nE;
    NumberStatistics(std::string name_ = "part"):n(0),n2(0),n3(0),n4(0),n5(0),n6(0),n7(0),n8(0),n9(0),n10(0),n11(0),n12(0),events(0),name(name_) {
        means.resize(13);
        wsum  = 0.;
        w2sum = 0.;
        nE    = 0.;
        acc      = false;
    binom.resize(BINOMSZ);
    for(int i=0;i<binom.size();++i) binom[i].resize(i+1);
    binom[0][0] = 1;
    binom[1][0] = 1;
    binom[1][1] = 1;
    for(int n=2;n<BINOMSZ;++n) {
      binom[n][0] = binom[n][n] = 1;
      for(int mm=1;mm<n;++mm)
        binom[n][mm] = binom[n-1][mm] + binom[n-1][mm-1];
    }
    m.resize(13);
    }
    ~NumberStatistics() {
    }
    void Reset() {
        n = n2 = n3 = n4 = 0;
        n5 = n6 = n7 = n8 = 0;
        events = 0;
        wsum  = 0.;
        w2sum = 0.;
        nE    = 0.;
    }

    void SetAcceptance(bool acc_) { acc = acc_; }

    bool GetAcceptance() const { return acc; }

  std::string GetName() const { return name; }

    void AddEvent(int tmpn, double weight=1.);
    void AddEvent(double tmpn, double weight=1.);

  double Cnm(int nn, int m) const {
    if (m<0 || m>nn) return -1;
    if (m==0 || m==nn) return 1;
    if (nn<BINOMSZ) return binom[nn][m];
    return Cnm(nn-1, m)  + Cnm(nn-1, m-1);
  }

    double GetMean() const;
    double GetMeanError() const;
    double GetVariance() const;
    double GetStdDev() const;
  
    double GetScaledVariance() const;
    double GetSkewness() const;
    double GetSkewnessError() const;
  double GetSkewnessError2() const;
    double GetKurtosis() const;
    double GetKurtosisError() const;
  double GetKurtosisError2() const;

    double GetN2Error2() const;
    double GetVarianceError() const;
    double GetScaledVarianceError() const;
  double GetStdDevError() const;

  double GetC5C2() const;
  double GetC5C2Error() const;
  double GetC5C2Error2() const;

  double GetC6C2() const;
  double GetC6C2Error() const;
  double GetC6C2Error2() const;

  double GetC6C2alt() const;
  double GetC6C2Erroralt() const;

  void CalculateCentralMoments();
};

