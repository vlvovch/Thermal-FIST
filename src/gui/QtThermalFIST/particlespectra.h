/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef PARTICLESPECTRA_H
#define PARTICLESPECTRA_H

#include <string>
#include <map>
#include <vector>
#include <cstdio>

#include "HRGBase/ThermalModelBase.h"
#include "HRGEventGenerator/SimpleParticle.h"
#include "HRGEventGenerator/SimpleEvent.h"
#include "HRGEventGenerator/MomentumDistribution.h"
#include "NumberStatistics.h"

struct stval {
    double sum;
    double sumsqr;
    int nev;
    double wsum, w2sum;
    stval() { sum = sumsqr = 0.; nev = 0; wsum = w2sum = 0.; }
    //double getAv() const { return sum / nev; }
    double getAv() const { return sum / wsum; }
    double getSigma() const { return sqrt((sumsqr / wsum - getAv() * getAv()) / (wsum * wsum / w2sum - 1.)); }
    //double getSigma() const { return sqrt((sumsqr - sum*sum/nev)/(nev-1.)/nev); }
    double getSigma2() const { return sqrt(sum) / nev; }
};

struct density
{
    double xst;
    double step;
    int counter;
    //double xend;
    int sz;
    std::vector<stval> data;
    std::vector<double> tsum;
    int valcount;
    int events;
    density()
    {
    }
    density(double xleft, double xright, int numb)
    {
        if (xright<xleft) printf("Right limit is less than left!\n");
        sz = numb;
        counter = 0;
        step = (xright-xleft)/sz;
        xst = xleft + 0.5*step;
        data.resize(sz);
        tsum.resize(sz);
        for(int i=0;i<tsum.size();++i) tsum[i] = 0.;
        valcount = 0;
        events = 0;
    }
    void insert(double value)
    {
        int tind = (int)((value-(xst-0.5*step))/step);
        if (tind>=0 && tind<sz)
        {
            tsum[tind]++;
        }
        counter++;
    }
    void updateEvent(double weight = 1.)
    {
        for(int i=0;i<tsum.size();++i) {
            data[i].sum    += weight * tsum[i];
            data[i].sumsqr += weight * tsum[i] * tsum[i];
            data[i].nev++;
            data[i].wsum  += weight;
            data[i].w2sum += weight*weight;
            tsum[i] = 0.;
        }
        events++;
    }
    double GetEntry(int ind) const {
        if (ind>=0 && ind<data.size()) {
            return data[ind].getAv() / step;
        }
        return 0.;
    }
    double GetEntryError(int ind) const {
        if (ind>=0 && ind<data.size()) {
            return data[ind].getSigma() / step;
        }
        return 0.;
    }
    double GetX(int ind) const {
        if (ind>=0 && ind<data.size()) {
            return xst + ind*step;
        }
        return 0.;
    }
    std::vector<double> GetXVector() const {
        std::vector<double> ret(data.size());
        for(int i=0;i<data.size();++i)
            ret[i] = xst + i*step;
        return ret;
    }
    std::vector<double> GetYVector() const {
        std::vector<double> ret(data.size(), 0.);
        if (events>0)
            for(int i=0;i<data.size();++i)
                ret[i] = data[i].getAv() / step;
        return ret;
    }
    std::vector<double> GetYErrorVector() const {
        std::vector<double> ret(data.size(), 0.);
        if (events>0)
            for(int i=0;i<data.size();++i)
                ret[i] = data[i].getSigma() / step;
        return ret;
    }
};

struct density2d
{
    double xst, yst;
    double stepx, stepy;
    int counter;
    int szx, szy;
    std::vector<stval> data;
    std::vector<double> tsum;
    int valcount;
    int events;
    density2d()
    {
    }
    density2d(double xleft, double xright, int numbx, double yleft, double yright, int numby)
    {
        if (xright<xleft || yright<yleft) printf("Right limit is less than left!\n");
        szx     = numbx;
        stepx   = (xright-xleft)/szx;
        xst     = xleft + 0.5*stepx;
        szy     = numby;
        stepy   = (yright-yleft)/szy;
        yst     = yleft + 0.5*stepy;
        counter = 0;
        data.resize(szx*szy);
        tsum.resize(szx*szy);
        for(int i=0;i<tsum.size();++i) tsum[i] = 0.;
        valcount = 0;
        events = 0;
    }
    void insert(double valuex, double valuey)
    {
        int tindx = (int)((valuex-(xst-0.5*stepx))/stepx);
        int tindy = (int)((valuey-(yst-0.5*stepy))/stepy);
        if (tindx>=0 && tindx<szx && tindy>=0 && tindy<szy)
        {
            int tind = tindx * szx + tindy;
            tsum[tind]++;
        }
        counter++;
    }
    void updateEvent(double weight = 1.)
    {
        for(int i=0;i<tsum.size();++i) {
            data[i].sum += weight * tsum[i];
            data[i].sumsqr += weight * tsum[i] * tsum[i];
            data[i].nev++;
            data[i].wsum += weight;
            data[i].w2sum += weight*weight;
            tsum[i] = 0.;
        }
        events++;
    }
    double GetEntry(int ind) const {
        if (ind>=0 && ind<data.size()) {
            return data[ind].getAv() / stepx / stepy;
        }
        return 0.;
    }
    double GetEntryError(int ind) const {
        if (ind>=0 && ind<data.size()) {
            return data[ind].getSigma() / stepx / stepy;
        }
        return 0.;
    }
    std::vector<double> GetXVector() const {
        std::vector<double> ret(data.size());
        for(int i=0;i<data.size();++i)
            ret[i] = xst + (i/szx)*stepx;
        return ret;
    }
    std::vector<double> GetYVector() const {
        std::vector<double> ret(data.size());
        for(int i=0;i<data.size();++i)
            ret[i] = yst + (i%szy)*stepy;
        return ret;
    }
    std::vector<double> GetZVector() const {
        std::vector<double> ret(data.size(), 0.);
        if (events>0)
            for(int i=0;i<data.size();++i)
                ret[i] = data[i].getAv() / stepx / stepy;
        return ret;
    }
    std::vector<double> GetZErrorVector() const {
        std::vector<double> ret(data.size(), 0.);
        if (events>0)
            for(int i=0;i<data.size();++i)
                ret[i] = data[i].getSigma() / stepx / stepy;
        return ret;
    }
};



class ParticleSpectrum
{
    double n, n2, n3, n4;
    double n5, n6, n7, n8;
    std::vector<double> means;
    int events;
    double wsum;
    double w2sum;
    bool isAveragesCalculated;
    double nE;
    int PDGID;
    double mass;
    int tmpn;
    thermalfist::MomentumDistributionBase *fDistribution;
    std::vector<double> xs;
    std::vector<double> ys;
    bool acc;
    double pTsum_event, pT2sum_event;
    double pTsum, pT2sum;
    double pTcnt, pTcnt2;
public:
    density   dndp;
    density   dndy;
    density   dndmt;
    density   dndpt;
    density2d d2ndptdy;
    ParticleSpectrum(int PDGID_ = 0, double mass_ = 1., double etamax = 0.,
      int bins = 500, int binsx = 40, int binsy = 40) :n(0), n2(0), n3(0), n4(0), n5(0), n6(0), n7(0), n8(0), wsum(0.), w2sum(0.), events(0), PDGID(PDGID_), mass(mass_), tmpn(0) {
        fDistribution = NULL;
        dndp     = density(0., 2. + mass, bins);
        dndy     = density(-3. - etamax, 3. + etamax, bins);
        dndmt    = density(mass, 2.+mass, bins);
        dndpt    = density(0., 3. + mass, bins);
        d2ndptdy = density2d(-3. - etamax, 3. + etamax, binsx, 0., 2., binsy);
        acc      = false;
        isAveragesCalculated = false;
        tmpn = 0;
        means.resize(9);
        pTsum = pT2sum = pTcnt = pTcnt2 = 0.;
        pTsum_event = pT2sum_event = 0.;
    }
    ~ParticleSpectrum() {
    }
    void Reset(int bins = 500, int binsx = 40, int binsy = 40) {
        n = n2 = n3 = n4 = 0;
        n5 = n6 = n7 = n8 = 0;
        events = 0;
        wsum = 0.;
        w2sum = 0.;
        tmpn = 0;
        dndp = density(0., 2. + mass, bins);
        dndy = density(-3., 3., bins);
        dndmt = density(mass, 2.+mass, bins);
        dndpt = density(0., 3. + mass, bins);
        d2ndptdy = density2d(-3., 3., binsx, 0., 2., binsy);
        pTsum = pT2sum = pTcnt = pTcnt2 = 0.;
        pTsum_event = pT2sum_event = 0.;
    }
    void SetDistribution(thermalfist::MomentumDistributionBase *distr) {
        fDistribution = distr;
    }
    void SetAcceptance(bool acc_) { acc = acc_; }

    bool GetAcceptance() const { return acc; }

    void CalculateAverages();

    int GetPDGID() const { return PDGID; }
    int GetEvents() const { return events; }
    void AddParticle(const thermalfist::SimpleParticle &part);
    void FinishEvent(double weight = 1.);

    double GetMean() const;
    double GetMeanError() const;
    double GetVariance() const;
    double GetStdDev() const;
    double GetScaledVariance() const;
    double GetSkewness() const;
    double GetSkewnessError() const;
    double GetKurtosis() const;
    double GetKurtosisError() const;

    double GetN2Error2() const;
    double GetVarianceError() const;
    double GetScaledVarianceError() const;

    double GetMeanPt() const;
    double GetMeanPtError() const;


    double GetModeldNdp(double p) const {
        if (fDistribution!=NULL && !fDistribution->isNormalized()) fDistribution->Normalize();
        if (fDistribution==NULL) return 0.;
        else return fDistribution->dndp(p);
    }
    double GetModeldNdy(double y) const {
        if (fDistribution!=NULL && !fDistribution->isNormalized()) fDistribution->Normalize();
        if (fDistribution==NULL) return 0.;
        else return fDistribution->dndy(y);
    }
    double GetModeldNmtdmt(double mt) const {
        if (fDistribution!=NULL && !fDistribution->isNormalized()) fDistribution->Normalize();
        if (fDistribution==NULL) return 0.;
        else return fDistribution->dnmtdmt(mt);
    }
    double GetModeldNdpt(double pt) const {
      if (fDistribution != NULL && !fDistribution->isNormalized()) fDistribution->Normalize();
      if (fDistribution == NULL) return 0.;
      else return fDistribution->dnmtdmt(sqrt(mass*mass + pt*pt)) * pt;
    }
    double GetModeld2Ndptdy(double pt, double y) const {
        if (fDistribution!=NULL && !fDistribution->isNormalized()) fDistribution->Normalize();
        if (fDistribution==NULL) return 0.;
        else return fDistribution->d2ndptdy(pt, y);
    }
    void FillModelDistribution(double tl = 0., double tr = 0., int iters = 0., int type = 0) {
        xs.resize(0);
        ys.resize(0);
        double dx = (tr - tl) / (iters);
        for(int i=0;i<iters;++i) {
            double tx = tl + (0.5+i) * dx;
            double ty = 0.;
            if (type==2) { tx += mass; }
            if (type==0) ty = GetModeldNdp(tx);
            else if (type==1) ty = GetModeldNdy(tx);
            else if (type == 4) ty = GetModeldNdpt(tx);
            else ty = GetModeldNmtdmt(tx);
            if (type!=2) xs.push_back(tx);
            else xs.push_back(tx-mass);
            ys.push_back(ty);
        }
    }
    std::vector<double> GetModelX() const {
        return xs;
    }
    std::vector<double> GetModelY() const {
        std::vector<double> ret = ys;
        double mean = GetMean();
        if (events>0)
            for(int i=0;i<ret.size();++i)
                ret[i] *= mean;
        return ret;
    }

    std::vector<double> GetXVector(int type=0) const {
        if (type==0) return dndp.GetXVector();
        else if (type==1) return dndy.GetXVector();
        else if (type==2) {
            std::vector<double> ret = dndmt.GetXVector();
            for(int i=0;i<ret.size();++i) ret[i] -= mass;
            return ret;
        }
        //else if (type == 4) {
        //  std::vector<double> ret = dndmt.GetXVector();
        //  for (int i = 0; i < ret.size(); ++i)
        //    ret[i] = sqrt(ret[i] * ret[i] - mass * mass);
        //  return ret;
        //}
        else if (type == 4) return dndpt.GetXVector();
        else return d2ndptdy.GetXVector();
    }
    std::vector<double> GetYVector(int type=0) const {
      if (type == 0) return dndp.GetYVector();
      else if (type == 1) return dndy.GetYVector();
      else if (type == 2) {
        std::vector<double> ret = dndmt.GetYVector();
        for (int i = 0; i < ret.size(); ++i) ret[i] /= dndmt.GetX(i);
        return ret;
      }
      else if (type == 4) return dndpt.GetYVector();
      //{
      //  std::vector<double> ret = dndmt.GetYVector();
      //  for (int i = 0; i < ret.size(); ++i) ret[i] *= sqrt(dndmt.GetX(i) * dndmt.GetX(i) - mass * mass) / dndmt.GetX(i);
      //  return ret;
      //}
      else return d2ndptdy.GetYVector();
    }
    std::vector<double> GetYErrorVector(int type=0) const {
        if (type==0) return dndp.GetYErrorVector();
        else if (type==1) return dndy.GetYErrorVector();
        else if (type == 3) {
            std::vector<double> ret = dndmt.GetYErrorVector();
            for(int i=0;i<ret.size();++i) ret[i] /= dndmt.GetX(i);
            return ret;
        }
        else return dndpt.GetYErrorVector();
        //{
        //  std::vector<double> ret = dndmt.GetYErrorVector();
        //  for (int i = 0; i < ret.size(); ++i) ret[i] *= sqrt(dndmt.GetX(i) * dndmt.GetX(i) - mass * mass) / dndmt.GetX(i);
        //  return ret;
        //}
    }
    std::vector<double> GetZVector(int type=0) const {
        return d2ndptdy.GetZVector();
    }
    std::vector<double> GetZErrorVector(int type=0) const {
        return d2ndptdy.GetZErrorVector();
    }
};

struct ParticleSpectraConfig {
  int fDistrType;
  double fEtaMax;
  int fStableOnly; // 0 - all, 1 - only stable, 2 - stable + list
  std::set<long long> fPdgCodes;
  double fT, fBeta, fNPow;
  int fBins, fBinsX, fBinsY;
  ParticleSpectraConfig(int distrtype = 0, double etamax = 0.5, int stableonly = 1,
    double T = 0.120, double beta = 0.5, double npow = 1.0,
    int bins = 500, int binsX = 40, int binsY = 40) :
    fDistrType(distrtype), fEtaMax(etamax), fStableOnly(stableonly),
    fT(T), fBeta(beta), fNPow(npow), fBins(bins), fBinsX(binsX), fBinsY(binsY) { }
};

class ParticlesSpectra {
public:
    std::vector<ParticleSpectrum> fParticles;
    std::vector<std::string> fNames;
    std::vector<NumberStatistics> fNetParticles;
    std::vector<NumberStatistics> fNetCharges;
    std::vector<NumberStatistics> fTotalCharges;
    std::vector<NumberStatistics> fPositiveCharges;
    std::vector<NumberStatistics> fNegativeCharges;
    std::vector< std::vector<int> > fParticleCharges;
    std::vector<double> fMasses;
    std::map<int, int> fPDGtoID;
    std::map<int, int> fPDGtoIDnet;
    std::map<int, int> fPDGtoIDall;
    std::vector<thermalfist::MomentumDistributionBase*> distrs;
    double fEtaMax;
    int fDistributionType;
    ParticleSpectraConfig fConfig;
    ParticlesSpectra(thermalfist::ThermalModelBase* model = NULL,
      const ParticleSpectraConfig& config = ParticleSpectraConfig());
      //double T = 0.120, double beta = 0.5, int distrtype = 0, double etamax = 0.5);
    ~ParticlesSpectra();
    void ProcessEvent(const thermalfist::SimpleEvent &evt);
    void Reset();
    void Reset(thermalfist::ThermalParticleSystem* TPS,
      const ParticleSpectraConfig& config = ParticleSpectraConfig());

      //double T = 0.120, double beta = 0.5, int distrtype = 0, double etamax = 0.5, double npow = 1.);
    void Reset(thermalfist::ThermalModelBase *model,
      const ParticleSpectraConfig& config = ParticleSpectraConfig());
      //double T = 0.120, double beta = 0.5, int distrtype = 0, double etamax = 0.5, double npow = 1.);
};



#endif // PARTICLESPECTRA_H
