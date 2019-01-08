/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ThermalParticleSystem.h"

#include <fstream>
#include <algorithm>
#include <queue>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>
#include <cmath>
#include <cstdlib>

using namespace std;

namespace thermalfist {

  namespace {
    /// For sorting particles by mass (default)
    bool cmpParticleMass(const ThermalParticle &a, const ThermalParticle &b) {
      return (a.Mass() < b.Mass());
    }
    
    /// For sorting particles by the quark content
    bool cmpParticlePDG(const ThermalParticle &a, const ThermalParticle &b) {
      if (abs(a.BaryonCharge()) != abs(b.BaryonCharge())) return (abs(a.BaryonCharge()) < abs(b.BaryonCharge()));
      if (abs(a.Charm()) != abs(b.Charm())) return (abs(a.Charm()) < abs(b.Charm()));
      if (abs(a.Strangeness()) != abs(b.Strangeness())) return (abs(a.Strangeness()) < abs(b.Strangeness()));
      return (a.Mass() < b.Mass());
    }
  }

  ThermalParticleSystem::ThermalParticleSystem(std::string InputFile, bool GenAntiP, double mcut)
  {
    m_NumberOfParticles = 0;
    m_Particles.resize(0);
    m_PDGtoID.clear();

    m_DecayContributionsByFeeddown.resize(Feeddown::NumberOfTypes);

    SetResonanceWidthIntegrationType(ThermalParticle::ZeroWidth);
    SetCalculationType(IdealGasFunctions::Quadratures);

    LoadTable(InputFile, GenAntiP, mcut);
  }


  ThermalParticleSystem::~ThermalParticleSystem(void)
  {
  }

  ThermalParticle::ParticleDecaysVector ThermalParticleSystem::GetDecaysFromAntiParticle(const ThermalParticle::ParticleDecaysVector& Decays) {
    ThermalParticle::ParticleDecaysVector ret = Decays;
    for (unsigned int i = 0; i < ret.size(); ++i) {
      for (unsigned int j = 0; j < ret[i].mDaughters.size(); ++j) {
        if (m_PDGtoID.count(-ret[i].mDaughters[j]) > 0) ret[i].mDaughters[j] = -ret[i].mDaughters[j];
      }
    }
    return ret;
  }

  void ThermalParticleSystem::ProcessDecays()
  {
    FillResonanceDecays(); 
    FillResonanceDecaysByFeeddown();
  }

  void ThermalParticleSystem::FillDecayProperties()
  {
    for (int i = 0; i < m_Particles.size(); ++i) {
      if (m_Particles[i].Decays().size() != 0) {
        double tsumb = 0.;

        for (int j = 0; j < m_Particles[i].Decays().size(); ++j) {

          m_Particles[i].Decays()[j].mPole = m_Particles[i].Mass();

          std::string tname = "";

          double M0 = 0.;
          double tS = 0.;
          for (int k = 0; k < m_Particles[i].Decays()[j].mDaughters.size(); ++k) {
            int tid = PdgToId(m_Particles[i].Decays()[j].mDaughters[k]);
            if (tid != -1) {
              M0 += m_Particles[tid].Mass();
              tS += max(0., (m_Particles[tid].Degeneracy() - 1.) / 2.);
              if (k != 0)
                tname += " - ";
              tname += m_Particles[tid].Name();
            }
          }
          m_Particles[i].Decays()[j].mM0 = M0;
          m_Particles[i].Decays()[j].mL = abs(max(0., (m_Particles[i].Degeneracy() - 1.) / 2.) - tS);

          tsumb += m_Particles[i].Decays()[j].mBratio;
          m_Particles[i].Decays()[j].mBratioAverage = m_Particles[i].Decays()[j].mBratio;

          if (tname.size() > 0)
            m_Particles[i].Decays()[j].mChannelName = tname;
        }
      }
      m_Particles[i].CalculateAndSetDynamicalThreshold();
      m_Particles[i].FillCoefficientsDynamical();
    }
  }

  void ThermalParticleSystem::FillDecayThresholds()
  {
    for (int i = 0; i < m_Particles.size(); ++i) {
      if (m_Particles[i].Decays().size() != 0) {
        for (int j = 0; j < m_Particles[i].Decays().size(); ++j) {
          double M0 = 0.;
          for (int k = 0; k < m_Particles[i].Decays()[j].mDaughters.size(); ++k) {
            if (PdgToId(m_Particles[i].Decays()[j].mDaughters[k]) != -1)
              M0 += m_Particles[PdgToId(m_Particles[i].Decays()[j].mDaughters[k])].Mass();
          }
          m_Particles[i].Decays()[j].mM0 = M0;
        }
        m_Particles[i].FillCoefficients();
      }
    }
  }

  void ThermalParticleSystem::FillResonanceDecays() {
    m_DecayContributionsByFeeddown[Feeddown::StabilityFlag].resize(m_Particles.size());
    m_DecayCumulants.resize(m_Particles.size());
    m_DecayProbabilities.resize(m_Particles.size());
    for (int i = 0; i < m_Particles.size(); ++i) {
      m_DecayContributionsByFeeddown[Feeddown::StabilityFlag][i].resize(0);
      m_DecayProbabilities[i].resize(0);
      m_DecayCumulants[i].resize(0);
    }
    for (int i = m_Particles.size() - 1; i >= 0; i--)
      if (!m_Particles[i].IsStable()) {
        GoResonance(i, i, 1.);
      }

    for (int i = 0; i < m_Particles.size(); ++i) {
      for (int j = 0; j < m_DecayContributionsByFeeddown[Feeddown::StabilityFlag][i].size(); ++j) {
        SingleDecayContribution &DecayContrib = m_DecayContributionsByFeeddown[Feeddown::StabilityFlag][i][j];
        vector<double> tmp = GoResonanceDecayProbs(DecayContrib.second, i, true);
        if (tmp.size() > 1) m_DecayProbabilities[i].push_back(make_pair(tmp, DecayContrib.second));
      }
      for (int j = 0; j < m_DecayProbabilities[i].size(); ++j) {
        double tmp = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
        for (int jj = 0; jj < m_DecayProbabilities[i][j].first.size(); ++jj) {
          tmp += m_DecayProbabilities[i][j].first[jj] * jj;
          tmp2 += m_DecayProbabilities[i][j].first[jj] * jj * jj;
          tmp3 += m_DecayProbabilities[i][j].first[jj] * jj * jj * jj;
          tmp4 += m_DecayProbabilities[i][j].first[jj] * jj * jj * jj * jj;
        }
        double n2 = 0., n3 = 0., n4 = 0.;
        n2 = tmp2 - tmp * tmp;
        n3 = tmp3 - 3. * tmp2 * tmp + 2. * tmp * tmp * tmp;
        n4 = tmp4 - 4. * tmp3 * tmp + 6. * tmp2 * tmp * tmp - 3. * tmp * tmp * tmp * tmp - 3. * n2 * n2;
        vector<double> moments(0);
        moments.push_back(tmp);
        moments.push_back(n2);
        moments.push_back(n3);
        moments.push_back(n4);
        m_DecayCumulants[i].push_back(make_pair(moments, m_DecayProbabilities[i][j].second));
      }
    }

    m_DecayDistributionsMap.resize(m_Particles.size());
    m_ResonanceFinalStatesDistributions.resize(m_Particles.size());
    for (int i = 0; i < m_Particles.size(); ++i) {
      m_ResonanceFinalStatesDistributions[i].resize(0);
      m_DecayDistributionsMap[i].resize(0);
    }
    for (int i = 0; i < m_Particles.size(); ++i) {
      m_ResonanceFinalStatesDistributions[i] = GoResonanceDecayDistributions(i, true);
    }
    // Clear m_DecayDistributionsMap and memory it occupies
    std::vector< std::vector< std::pair<double, std::vector<int> > > >().swap(m_DecayDistributionsMap);

    for (int i = 0; i < m_Particles.size(); ++i) {
      vector<int> nchtyp(0);
      nchtyp.push_back(0);
      nchtyp.push_back(1);
      nchtyp.push_back(-1);

      m_Particles[i].Nch().resize(0);
      m_Particles[i].DeltaNch().resize(0);

      for (int nti = 0; nti < 3; nti++) {
        vector<double> prob = GoResonanceDecayProbsCharge(i, nchtyp[nti], true);
        double tmp = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
        for (int jj = 0; jj < prob.size(); ++jj) {
          tmp += prob[jj] * jj;
          tmp2 += prob[jj] * jj * jj;
          tmp3 += prob[jj] * jj * jj * jj;
          tmp4 += prob[jj] * jj * jj * jj * jj;
        }
        double n2 = 0., n3 = 0., n4 = 0.;
        n2 = tmp2 - tmp * tmp;
        n3 = tmp3 - 3. * tmp2 * tmp + 2. * tmp * tmp * tmp;
        n4 = tmp4 - 4. * tmp3 * tmp + 6. * tmp2 * tmp * tmp - 3. * tmp * tmp * tmp * tmp - 3. * n2 * n2;
        m_Particles[i].Nch().push_back(tmp);
        m_Particles[i].DeltaNch().push_back(n2);
      }

    }
  }


  void ThermalParticleSystem::GoResonance(int ind, int startind, double BR) {
    DecayContributionsToParticle& DecayContrib = m_DecayContributionsByFeeddown[Feeddown::StabilityFlag][ind];
    if (ind != startind && DecayContrib.size() > 0 && DecayContrib[DecayContrib.size() - 1].second == startind)
    {
      DecayContrib[DecayContrib.size() - 1].first += BR;
    }
    else if (ind != startind) 
      DecayContrib.push_back(make_pair(BR, startind));

    if (!m_Particles[ind].IsStable()) {
      for (int i = 0; i < m_Particles[ind].Decays().size(); ++i) {
        const ParticleDecayChannel& decaychannel = m_Particles[ind].Decays()[i];
        double tbr = decaychannel.mBratio;

        if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && ind == startind)
          tbr = decaychannel.mBratioAverage;

        for (int j = 0; j < decaychannel.mDaughters.size(); ++j) {
          if (m_PDGtoID.count(decaychannel.mDaughters[j]) != 0)
            GoResonance(m_PDGtoID[decaychannel.mDaughters[j]], startind, BR*tbr);
        }
      }
    }
  }

  std::vector<double> ThermalParticleSystem::GoResonanceDecayProbs(int ind, int goalind, bool firstdecay) {
    std::vector<double> ret(1, 0.);
    if (m_Particles[ind].IsStable()) {
      if (ind == goalind) ret.push_back(1.);
      else ret[0] = 1.;
      return ret;
    }
    else if (ind == goalind) {
      ret.push_back(1.);
      return ret;
    }
    else {
      ret[0] = 0.;
      vector<double> tret;
      for (int i = 0; i < m_Particles[ind].Decays().size(); ++i) {
        double tbr = m_Particles[ind].Decays()[i].mBratio;
        if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && firstdecay)
          tbr = m_Particles[ind].Decays()[i].mBratioAverage;

        tret.resize(1);
        tret[0] = 1.;
        for (int j = 0; j < m_Particles[ind].Decays()[i].mDaughters.size(); ++j) {
          if (m_PDGtoID.count(m_Particles[ind].Decays()[i].mDaughters[j]) != 0) {
            vector<double> tmp = GoResonanceDecayProbs(m_PDGtoID[m_Particles[ind].Decays()[i].mDaughters[j]], goalind);
            vector<double> tmp2(tret.size() + tmp.size() - 1, 0.);
            for (int i1 = 0; i1 < tret.size(); ++i1)
              for (int i2 = 0; i2 < tmp.size(); ++i2)
                tmp2[i1 + i2] += tret[i1] * tmp[i2];
            tret = tmp2;
          }
        }
        if (ret.size() < tret.size()) ret.resize(tret.size(), 0.);
        for (int j = 0; j < tret.size(); ++j)
          ret[j] += tbr * tret[j];
      }
      double totprob = 0.;
      for (int i = 0; i < ret.size(); ++i)
        totprob += ret[i];
      if (totprob > 1.) {
        for (int i = 0; i < ret.size(); ++i)
          ret[i] *= 1. / totprob;
      }
      else {
        ret[0] += 1. - totprob;
      }
      return ret;
    }
    return ret;
  }

  std::vector<double> ThermalParticleSystem::GoResonanceDecayProbsCharge(int ind, int nch, bool firstdecay)
  {
    bool fl = false;
    int tQ = m_Particles[ind].ElectricCharge();
    if (nch == 0 && tQ != 0)
      fl = true;
    if (nch == 1 && tQ > 0)
      fl = true;
    if (nch == -1 && tQ < 0)
      fl = true;

    std::vector<double> ret(1, 0.);
    if (m_Particles[ind].IsStable()) {
      if (fl) ret.push_back(1.);
      else ret[0] = 1.;
      return ret;
    }
    else {
      ret[0] = 0.;
      vector<double> tret;
      for (int i = 0; i < m_Particles[ind].Decays().size(); ++i) {
        double tbr = m_Particles[ind].Decays()[i].mBratio;
        if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && firstdecay)
          tbr = m_Particles[ind].Decays()[i].mBratioAverage;

        tret.resize(1);
        tret[0] = 1.;
        for (int j = 0; j < m_Particles[ind].Decays()[i].mDaughters.size(); ++j) {
          if (m_PDGtoID.count(m_Particles[ind].Decays()[i].mDaughters[j]) != 0) {
            vector<double> tmp = GoResonanceDecayProbsCharge(m_PDGtoID[m_Particles[ind].Decays()[i].mDaughters[j]], nch);
            vector<double> tmp2(tret.size() + tmp.size() - 1, 0.);
            for (int i1 = 0; i1 < tret.size(); ++i1)
              for (int i2 = 0; i2 < tmp.size(); ++i2)
                tmp2[i1 + i2] += tret[i1] * tmp[i2];
            tret = tmp2;
          }
        }
        if (ret.size() < tret.size())
          ret.resize(tret.size(), 0.);
        for (int j = 0; j < tret.size(); ++j)
          ret[j] += tbr * tret[j];
      }
      double totprob = 0.;
      for (int i = 0; i < ret.size(); ++i)
        totprob += ret[i];
      if (totprob > 1.) {
        for (int i = 0; i < ret.size(); ++i)
          ret[i] *= 1. / totprob;
      }
      else {
        ret[0] += 1. - totprob;
      }
      return ret;
    }
    return ret;
  }

  ThermalParticleSystem::ResonanceFinalStatesDistribution ThermalParticleSystem::GoResonanceDecayDistributions(int ind, bool firstdecay)
  {
    if (!firstdecay && m_DecayDistributionsMap[ind].size() != 0)
      return m_DecayDistributionsMap[ind];


    std::vector<std::pair<double, std::vector<int>>> retorig(1);
    retorig[0].first = 1.;
    retorig[0].second = std::vector<int>(m_Particles.size(), 0);
    retorig[0].second[ind] = 1;

    std::vector<std::pair<double, std::vector<int>>> ret(0);

    ThermalParticle &tpart = m_Particles[ind];

    if (tpart.IsStable()) {
      m_ResonanceFinalStatesDistributions[ind] = retorig;
      return retorig;
    }

    for (int i = 0; i < tpart.Decays().size(); ++i) {
      double tbr = tpart.Decays()[i].mBratio;
      if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && firstdecay)
        tbr = m_Particles[ind].Decays()[i].mBratioAverage;

      std::vector<std::pair<double, std::vector<int> > > tret = retorig;

      for (int j = 0; j < tpart.Decays()[i].mDaughters.size(); ++j) {
        if (m_PDGtoID.count(tpart.Decays()[i].mDaughters[j]) != 0) {

          std::vector<std::pair<double, std::vector<int> > > tmp = GoResonanceDecayDistributions(m_PDGtoID[tpart.Decays()[i].mDaughters[j]]);
          std::vector<std::pair<double, std::vector<int> > > tmp2(tret.size() * tmp.size());
          for (int i1 = 0; i1 < tret.size(); ++i1) {
            for (int i2 = 0; i2 < tmp.size(); ++i2) {
              tmp2[i1*tmp.size() + i2].first = tret[i1].first * tmp[i2].first;
              tmp2[i1*tmp.size() + i2].second.resize(m_Particles.size());
              for (int jj = 0; jj < tmp2[i1*tmp.size() + i2].second.size(); ++jj)
                tmp2[i1*tmp.size() + i2].second[jj] = tret[i1].second[jj] + tmp[i2].second[jj];
            }
          }
          tret = tmp2;

          // Restrict maximum number of channels to 1500, otherwise memory is an issue, relevant for the THERMUS-3.0 table
          if (tret.size() > 1500) {
            printf("**WARNING** %s (%d) Decay Distributions: Too large array, cutting the number of channels to 1500!\n",
              m_Particles[ind].Name().c_str(),
              m_Particles[ind].PdgId());
            CuteHRGHelper::cutDecayDistributionsVector(tret);
          }
        }
      }

      for (int j = 0; j < tret.size(); ++j) {
        tret[j].first *= tbr;
        ret.push_back(tret[j]);
      }
    }

    // Restrict maximum number of channels to 1500, otherwise memory is an issue, relevant for the THERMUS-3.0 table
    if (ret.size() > 1500) {
      printf("**WARNING** %s (%d) Decay Distributions: Too large array, cutting the number of channels to 1500!\n",
        m_Particles[ind].Name().c_str(),
        m_Particles[ind].PdgId());
      CuteHRGHelper::cutDecayDistributionsVector(ret);
    }

    double totprob = 0.;
    for (int i = 0; i < ret.size(); ++i)
      totprob += ret[i].first;
    if (totprob > 1.) {
      for (int i = 0; i < ret.size(); ++i)
        ret[i].first *= 1. / totprob;
    }
    else if (totprob < 1.) {
      double emptyprob = 1. - totprob;
      ret.push_back(std::make_pair(emptyprob, retorig[0].second));
    }

    if (!firstdecay)
      m_DecayDistributionsMap[ind] = ret;

    // Debugging

    //if (tpart.BaryonCharge() == 1) {
    //  printf("Checking baryon number conservation in decays for %d\n", tpart.PdgId());
    //  double tBav = 0.;
    //  for (int i = 0; i < ret.size(); ++i) {
    //    double tbr = ret[i].first;
    //    for (int j = 0; j < ret[i].second.size(); ++j)
    //      if (m_Particles[j].IsStable())
    //        tBav += tbr * ret[i].second[j] * m_Particles[j].BaryonCharge();
    //  }
    //  printf("<B> = %lf\n", tBav);
    //  //printf("Decay distributions for %d\n", tpart.PdgId());
    //  //for (int i = 0; i < ret.size(); ++i) {
    //  //  printf("%lf\n", ret[i].first);
    //  //  for (int j = 0; j < ret[i].second.size(); ++j)
    //  //    if (ret[i].second[j] > 0)
    //  //      printf("%d: %d\n", m_Particles[j].PdgId(), ret[i].second[j]);
    //  //}
    //}

    return ret;
  }


  namespace CuteHRGHelper {
    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
      std::stringstream ss(s);
      std::string item;
      while (std::getline(ss, item, delim)) {
        elems.push_back(item);
      }
      return elems;
    }

    std::vector<std::string> split(const std::string &s, char delim) {
      std::vector<std::string> elems;
      split(s, delim, elems);
      return elems;
    }

    void cutDecayDistributionsVector(std::vector<std::pair<double, std::vector<int>>>& vect, int maxsize)
    {
      if (vect.size() > maxsize) {
        std::sort(vect.begin(), vect.end());
        std::reverse(vect.begin(), vect.end());
        vect.resize(1000);
      }
    }
  }

  void ThermalParticleSystem::LoadTable(std::string InputFile, bool GenAntiP, double mcut) {

    m_NumberOfParticles = 0;
    m_Particles.resize(0);
    m_PDGtoID.clear();

    m_NumBaryons = m_NumCharged = m_NumStrange = m_NumCharmed = 0;

    ifstream fin;
    fin.open(InputFile.c_str());
    if (fin.is_open()) {

      char tmpc[2000];
      fin.getline(tmpc, 2000);
      string tmp = string(tmpc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');

      int flnew = 0;
      if (elems.size() < 1 || CuteHRGHelper::split(elems[0], ';').size() < 4)
        flnew = 1;
      else
        flnew = 0;

      fin.clear();
      fin.seekg(0, ios::beg);

      if (flnew == 1)
        LoadTable_NewFormat(fin);
      else
        LoadTable_OldFormat(fin);

      fin.close();

    }

    string decayprefix = "";

    for (int i = InputFile.size() - 1; i >= 0; --i)
      if (InputFile[i] == '\\' || InputFile[i] == '/')
      {
        decayprefix = InputFile.substr(0, i + 1);
        break;
      }

    LoadDecays((decayprefix + "decays.dat").c_str(), GenAntiP);

    SetResonanceWidthIntegrationType(m_ResonanceWidthIntegrationType);
    SetCalculationType(m_QStatsCalculationType);
  }

  void ThermalParticleSystem::LoadTable_OldFormat(std::ifstream & fin, bool GenerateAntiParticles, double mcut)
  {
    if (fin.is_open()) {
      string tmp;
      char tmpc[500];
      fin.getline(tmpc, 500);
      tmp = string(tmpc);
      while (1) {
        vector<string> fields = CuteHRGHelper::split(tmp, ';');
        if (fields.size() < 14) break;
        int stable, pdgid, spin, stat, str, bary, chg, charm;
        double mass, width, threshold, abss, absc, radius = 0.5;
        string name, decayname = "";
        stable = atoi(fields[0].c_str());
        name = fields[1];
        pdgid = atoi(fields[2].c_str());
        spin = atoi(fields[3].c_str());
        stat = atoi(fields[4].c_str());
        mass = atof(fields[5].c_str());
        str = atoi(fields[6].c_str());
        bary = atoi(fields[7].c_str());
        chg = atoi(fields[8].c_str());
        charm = atoi(fields[9].c_str());
        abss = atof(fields[10].c_str());
        absc = atof(fields[11].c_str());
        width = atof(fields[12].c_str());
        threshold = atof(fields[13].c_str());
        if (fields.size() >= 15) radius = atof(fields[14].c_str());
        if (fields.size() == 16) decayname = fields[15];

        if (mass > mcut) {
          fin.getline(tmpc, 500);
          tmp = string(tmpc);
          continue;
        }

        if (bary != 0)  m_NumBaryons++;
        if (chg != 0)   m_NumCharged++;
        if (str != 0)   m_NumStrange++;
        if (charm != 0) m_NumCharmed++;

        m_Particles.push_back(ThermalParticle((bool)stable, name, pdgid, spin, stat, mass, str, bary, chg, abss, width, threshold, charm, absc, radius));
        m_NumberOfParticles++;

        if (GenerateAntiParticles && !(bary == 0 && chg == 0 && str == 0 && charm == 0)) {

          if (bary != 0)  m_NumBaryons++;
          if (chg != 0)   m_NumCharged++;
          if (str != 0)   m_NumStrange++;
          if (charm != 0) m_NumCharmed++;

          if (bary == 0 && name[name.size() - 1] == '+')
            name[name.size() - 1] = '-';
          else if (bary == 0 && name[name.size() - 1] == '-')
            name[name.size() - 1] = '+';
          else
            name = "anti-" + name;
          m_Particles.push_back(ThermalParticle((bool)stable, name, -pdgid, spin, stat, mass, -str, -bary, -chg, abss, width, threshold, -charm, absc, radius));
          m_Particles[m_Particles.size() - 1].SetAntiParticle(true);
        }

        fin.getline(tmpc, 500);
        tmp = string(tmpc);
      }

      FinalizeList();
    }
  }

  void ThermalParticleSystem::LoadTable_NewFormat(std::ifstream & fin, bool GenerateAntiParticles, double mcut)
  {
    if (fin.is_open()) {
      char cc[2000];
      while (!fin.eof()) {
        fin.getline(cc, 2000);
        string tmp = string(cc);
        vector<string> elems = CuteHRGHelper::split(tmp, '#');
        if (elems.size() < 1 || elems[0].size() == 0)
          continue;

        istringstream iss(elems[0]);

        int stable, pdgid, stat, str, bary, chg, charm;
        double mass, degeneracy, width, threshold, abss, absc;
        string name;

        if (iss >> pdgid
          >> name
          >> stable
          >> mass
          >> degeneracy
          >> stat
          >> bary
          >> chg
          >> str
          >> charm
          >> abss
          >> absc
          >> width
          >> threshold) {

          if (mass > mcut)
            continue;

          if (bary != 0)  m_NumBaryons++;
          if (chg != 0)   m_NumCharged++;
          if (str != 0)   m_NumStrange++;
          if (charm != 0) m_NumCharmed++;

          m_Particles.push_back(ThermalParticle((bool)stable, name, pdgid, degeneracy, stat, mass, str, bary, chg, abss, width, threshold, charm, absc));
          m_NumberOfParticles++;

          if (GenerateAntiParticles && !(bary == 0 && chg == 0 && str == 0 && charm == 0)) {

            if (bary != 0)  m_NumBaryons++;
            if (chg != 0)   m_NumCharged++;
            if (str != 0)   m_NumStrange++;
            if (charm != 0) m_NumCharmed++;

            if (bary == 0 && name[name.size() - 1] == '+')
              name[name.size() - 1] = '-';
            else if (bary == 0 && name[name.size() - 1] == '-')
              name[name.size() - 1] = '+';
            else
              name = "anti-" + name;
            m_Particles.push_back(ThermalParticle((bool)stable, name, -pdgid, degeneracy, stat, mass, -str, -bary, -chg, abss, width, threshold, -charm, absc));
            m_Particles[m_Particles.size() - 1].SetAntiParticle(true);
          }
        }
      }

      FinalizeList();
    }
  }

  void ThermalParticleSystem::SetTableFromVector(const std::vector<ThermalParticle>& part_in, bool GenerateAntiParticles)
  {
    m_Particles.resize(0);

    for (int i = 0; i < part_in.size(); ++i) {
      const ThermalParticle &part = part_in[i];
      if (!GenerateAntiParticles) {
        m_Particles.push_back(part);
      }
      else if (part.PdgId() > 0) {
        m_Particles.push_back(part);

        if (!part.IsNeutral()) {
          m_Particles.push_back(part.GenerateAntiParticle());
        }
      }
    }

    FinalizeList();

    if (GenerateAntiParticles) {
      for (int i = 0; i < m_Particles.size(); ++i) {
        if (m_Particles[i].IsAntiParticle() && PdgToId(-m_Particles[i].PdgId()) != -1) {
          m_Particles[i].SetDecays(GetDecaysFromAntiParticle(m_Particles[PdgToId(-m_Particles[i].PdgId())].Decays()));
        }
      }
    }

    FillDecayProperties();
    FillDecayThresholds();
    ProcessDecays();
  }

  void ThermalParticleSystem::WriteTableToFile(std::string OutputFile, bool WriteAntiParticles)
  {
    std::ofstream fout(OutputFile.c_str());
    if (fout.is_open()) {
      fout << "#"
        << std::setw(14) << "pdgid"
        << std::setw(20) << "name"
        << std::setw(15) << "stable"
        << std::setw(15) << "mass[GeV]"
        << std::setw(15) << "degeneracy"
        << std::setw(15) << "statistics"
        << std::setw(15) << "B"
        << std::setw(15) << "Q"
        << std::setw(15) << "S"
        << std::setw(15) << "C"
        << std::setw(15) << "|S|"
        << std::setw(15) << "|C|"
        << std::setw(15) << "width[GeV]"
        << std::setw(15) << "threshold[GeV]"
        << std::endl;

      for (int i = 0; i < m_Particles.size(); ++i) {
        const ThermalParticle& part = m_Particles[i];
        if (part.PdgId() < 0 && !WriteAntiParticles)
          continue;

        fout << std::setw(15) << part.PdgId()
          << std::setw(20) << part.Name()
          << std::setw(15) << static_cast<int>(part.IsStable())
          << std::setw(15) << part.Mass()
          << std::setw(15) << part.Degeneracy()
          << std::setw(15) << part.Statistics()
          << std::setw(15) << part.BaryonCharge()
          << std::setw(15) << part.ElectricCharge()
          << std::setw(15) << part.Strangeness()
          << std::setw(15) << part.Charm()
          << std::setw(15) << part.AbsoluteStrangeness()
          << std::setw(15) << part.AbsoluteCharm()
          << std::setw(15) << part.ResonanceWidth()
          << std::setw(15) << part.DecayThresholdMass()
          << std::endl;
      }
      fout.close();
    }
  }

  void ThermalParticleSystem::LoadDecays(std::string DecaysFile, bool GenerateAntiParticles)
  {
    for (int i = 0; i < m_Particles.size(); ++i)
      m_Particles[i].ClearDecays();

    ifstream fin(DecaysFile.c_str());

    if (fin.is_open()) {

      char tmpc[2000];
      fin.getline(tmpc, 2000);
      string tmp = string(tmpc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');

      int flnew = 0;
      if (tmp.size() == 0 || elems.size() >= 2)
        flnew = 1;
      else
        flnew = 0;

      fin.clear();
      fin.seekg(0, ios::beg);

      if (flnew == 1)
        ReadDecays_NewFormat(fin);
      else
        ReadDecays_OldFormat(fin);

      fin.close();

    }

    if (GenerateAntiParticles) {
      for (int i = 0; i < m_Particles.size(); ++i) {
        if (m_Particles[i].PdgId() < 0)
          m_Particles[i].SetDecays(GetDecaysFromAntiParticle(m_Particles[m_PDGtoID[-m_Particles[i].PdgId()]].Decays()));
      }
    }

    for (int i = 0; i < m_Particles.size(); ++i)
      m_Particles[i].SetDecaysOriginal(m_Particles[i].Decays());

    FillDecayProperties();
    FillDecayThresholds();
    ProcessDecays();
  }

  void ThermalParticleSystem::ReadDecays_NewFormat(std::ifstream & fin)
  {
    vector< ThermalParticle::ParticleDecaysVector > decays(0);
    map<int, int> decaymap;
    decaymap.clear();

    
    if (fin.is_open()) {
      char cc[2000];
      int index = 0;
      while (!fin.eof()) {
        fin.getline(cc, 2000);
        string tmp = string(cc);
        vector<string> elems = CuteHRGHelper::split(tmp, '#');
        if (elems.size() < 1 || elems[0].size() == 0)
          continue;

        int tpdgid, tdecaysnumber = 0;
        ThermalParticle::ParticleDecaysVector tdecays(0);

        istringstream iss(elems[0]);
        if (!(iss >> tpdgid)) continue;

        bool fl = false;
        while (!fl) {
          if (fin.eof()) break;
          fin.getline(cc, 500);
          tmp = string(cc);
          elems = CuteHRGHelper::split(tmp, '#');
          if (elems.size() < 1 || elems[0].size() == 0)
            continue;

          istringstream isstnum(elems[0]);
          if (!(isstnum >> tdecaysnumber)) {
            tdecaysnumber = 0;
            continue;
          }
          fl = true;
        }

        for (int i = 0; i < tdecaysnumber; ++i) {
          bool fl = false;
          while (!fl) {
            if (fin.eof()) break;
            fin.getline(cc, 500);
            tmp = string(cc);
            elems = CuteHRGHelper::split(tmp, '#');
            if (elems.size() < 1 || elems[0].size() == 0)
              continue;

            ParticleDecayChannel tdecay;
            istringstream issdec(elems[0]);
            if (!(issdec >> tdecay.mBratio)) continue;
            int tmpid;
            while (issdec >> tmpid) {
              tdecay.mDaughters.push_back(tmpid);
            }
            tdecays.push_back(tdecay);
            fl = true;
          }
        }

        if (tdecays.size() == tdecaysnumber && tdecays.size() != 0) {
          decays.push_back(tdecays);
          decaymap[tpdgid] = index;
          index++;
        }
      }

      for (int i = 0; i < m_Particles.size(); ++i) {
        if (decaymap.count(m_Particles[i].PdgId()) != 0)
          m_Particles[i].SetDecays(decays[decaymap[m_Particles[i].PdgId()]]);
      }
    }
  }

  void ThermalParticleSystem::WriteDecaysToFile(std::string OutputFile, bool WriteAntiParticles)
  {
    std::ofstream fout(OutputFile.c_str());
    if (fout.is_open()) {
      fout << "# the list of decays" << std::endl;
      fout << "# each entry consists of the following:" << std::endl;
      fout << "# a line with the pdgid of decaying particle" << std::endl;
      fout << "# a line with the number of decay channels" << std::endl;
      fout << "# for each channel a line containing whitespace-separated values of the channel branching ratio and pdg ids of the daughter products" << std::endl;
      fout << "# everything after the # symbol is treated as a comment and ignored" << std::endl;
      fout << "# decays of antiparticles are not listed but generated from the listed decays of particles" << std::endl;
      fout << std::endl;

      for (unsigned int i = 0; i < m_Particles.size(); ++i) {
        if (m_Particles[i].PdgId()>0 && m_Particles[i].Decays().size()>0) {
          fout << std::left << std::setw(36) << m_Particles[i].PdgId();
          fout << " # " << m_Particles[i].Name() << std::endl;

          fout << std::left << std::setw(36) << m_Particles[i].Decays().size();
          fout << " # " << m_Particles[i].Decays().size() << " decay channel";
          if (m_Particles[i].Decays().size() % 10 != 1 || m_Particles[i].Decays().size() % 100 == 11) fout << "s";
          fout << std::endl;

          for (unsigned int j = 0; j < m_Particles[i].Decays().size(); ++j) {
            fout << std::left << std::setw(15) << m_Particles[i].Decays()[j].mBratio << " ";
            std::ostringstream oss;
            for (unsigned int k = 0; k < m_Particles[i].Decays()[j].mDaughters.size(); ++k) {
              oss << m_Particles[i].Decays()[j].mDaughters[k];
              if (k != m_Particles[i].Decays()[j].mDaughters.size() - 1)
                oss << " ";
            }
            fout << std::left << std::setw(20) << oss.str();
            fout << " # " << m_Particles[i].Name() << " -> ";
            for (unsigned int k = 0; k < m_Particles[i].Decays()[j].mDaughters.size(); ++k) {
              if (m_PDGtoID.count(m_Particles[i].Decays()[j].mDaughters[k]) == 0) {
                if (m_Particles[i].Decays()[j].mDaughters[k] == 22) fout << "?gamma?";
                else fout << "???";
              }
              else
                fout << m_Particles[m_PDGtoID[m_Particles[i].Decays()[j].mDaughters[k]]].Name();
              if (k != m_Particles[i].Decays()[j].mDaughters.size() - 1)
                fout << " + ";
            }
            fout << std::endl;
          }
          fout << std::endl;
        }
      }

      fout.close();
    }
  }

  void ThermalParticleSystem::ReadDecays_OldFormat(std::ifstream & fin)
  {
    vector< ThermalParticle::ParticleDecaysVector > decays(0);
    vector<int> pdgids(0);
    map<int, int> decaymap;
    decaymap.clear();

    if (fin.is_open()) {
      int decaypartnumber = 0;
      fin >> decaypartnumber;
      decays.reserve(decaypartnumber);

      for (unsigned int i = 0; i < decaypartnumber; ++i) {
        int pdgid, decaysnumber, tmpid, daughters;
        double bratio;
        fin >> pdgid >> decaysnumber;
        decaymap[pdgid] = i;
        decays.push_back(ThermalParticle::ParticleDecaysVector(0));
        pdgids.push_back(pdgid);
        for (unsigned int j = 0; j < decaysnumber; ++j) {
          ParticleDecayChannel decay;
          fin >> bratio;
          decay.mBratio = bratio / 100.;
          fin >> daughters;
          decay.mDaughters.reserve(daughters);
          for (unsigned int k = 0; k < daughters; ++k) {
            fin >> tmpid;
            decay.mDaughters.push_back(tmpid);
          }
          decays[i].push_back(decay);
        }
      }
    }

    for (int i = 0; i < m_Particles.size(); ++i) {
      if (decaymap.count(m_Particles[i].PdgId()) != 0)
        m_Particles[i].SetDecays(decays[decaymap[m_Particles[i].PdgId()]]);
    }

  }

  std::string ThermalParticleSystem::GetNameFromPDG(int pdgid) {
    if (pdgid == 1) return string("Npart");
    if (pdgid == 33340) return string("Omega+Omegabar");
    if (pdgid == 22120) return string("p+n");
    if (m_PDGtoID.count(pdgid) == 0) return string("???");
    else return m_Particles[m_PDGtoID[pdgid]].Name();
  }

  void ThermalParticleSystem::NormalizeBranchingRatios() {
    for (int i = 0; i < m_Particles.size(); ++i) m_Particles[i].NormalizeBranchingRatios();
    ProcessDecays();
  }


  void ThermalParticleSystem::RestoreBranchingRatios() {
    for (int i = 0; i < m_Particles.size(); ++i) m_Particles[i].RestoreBranchingRatios();
    ProcessDecays();
  }

  void ThermalParticleSystem::SetCalculationType(IdealGasFunctions::QStatsCalculationType type)
  {
    m_QStatsCalculationType = type;
    for (int i = 0; i < m_Particles.size(); ++i)
      m_Particles[i].SetCalculationType(type);
  }

  void ThermalParticleSystem::SetClusterExpansionOrder(int order)
  {
    for (int i = 0; i < m_Particles.size(); ++i)
      m_Particles[i].SetClusterExpansionOrder(order);
  }

  void ThermalParticleSystem::SetResonanceWidthShape(ThermalParticle::ResonanceWidthShape shape)
  {
    for (int i = 0; i < m_Particles.size(); ++i)
      m_Particles[i].SetResonanceWidthShape(shape);
  }

  void ThermalParticleSystem::SetResonanceWidthIntegrationType(ThermalParticle::ResonanceWidthIntegration type)
  {
    bool dodecays = (type != m_ResonanceWidthIntegrationType);

    m_ResonanceWidthIntegrationType = type;

    for (int i = 0; i < m_Particles.size(); ++i)
      m_Particles[i].SetResonanceWidthIntegrationType(type);

    if (dodecays)
      ProcessDecays();
  }

  const ThermalParticle & ThermalParticleSystem::Particle(int id) const
  {
    if (id < 0 || id >= m_Particles.size()) {
      printf("**ERROR** ThermalParticleSystem::Particle(int id): id is out of bounds!");
      exit(1);
    }
    return m_Particles[id];
  }

  ThermalParticle & ThermalParticleSystem::Particle(int id)
  {
    if (id < 0 || id >= m_Particles.size()) {
      printf("**ERROR** ThermalParticleSystem::Particle(int id): id is out of bounds!\n");
      exit(1);
    }
    return m_Particles[id];
  }

  ThermalParticle & ThermalParticleSystem::ParticleByPDG(int pdgid)
  {
    if (m_PDGtoID.count(pdgid) == 0) {
      printf("**ERROR** ThermalParticleSystem::ParticleByPDG(int pdgid): pdgid %d is unknown\n", pdgid);
      exit(1);
    }
    return m_Particles[m_PDGtoID[pdgid]];
  }

  void ThermalParticleSystem::FillPdgMap()
  {
    m_NumBaryons = m_NumCharged = m_NumStrange = m_NumCharmed = 0;
    m_NumberOfParticles = 0;
    m_PDGtoID.clear();
    for (int i = 0; i < m_Particles.size(); ++i) {
      m_PDGtoID[m_Particles[i].PdgId()] = i;
      if (m_Particles[i].BaryonCharge() != 0)    m_NumBaryons++;
      if (m_Particles[i].ElectricCharge() != 0)  m_NumCharged++;
      if (m_Particles[i].Strangeness() != 0)     m_NumStrange++;
      if (m_Particles[i].Charm() != 0)           m_NumCharmed++;
      if (m_Particles[i].PdgId() > 0)            m_NumberOfParticles++;
    }

    for (int i = 0; i < m_DecayContributionsByFeeddown.size(); ++i)
      m_DecayContributionsByFeeddown[i].resize(m_Particles.size());
  }

  void ThermalParticleSystem::FinalizeList()
  {
    sort(m_Particles.begin(), m_Particles.end(), cmpParticleMass);
    //sort(m_Particles.begin(), m_Particles.end(), cmpParticlePDG);
    FillPdgMap();
    for (int i = 0; i < m_Particles.size(); ++i) {
      if (m_Particles[i].DecayType() == ParticleDecayType::Default)
        m_Particles[i].SetDecayType( DecayTypeByParticleType(m_Particles[i]) );
    }
  }

  void ThermalParticleSystem::AddParticle(const ThermalParticle & part)
  {
    m_Particles.push_back(part);
    FillPdgMap();
  }

  void ThermalParticleSystem::RemoveParticleAt(int ind)
  {
    if (ind >= 0 && ind < m_Particles.size()) {
      m_Particles.erase(m_Particles.begin() + ind);
      FillPdgMap();
    }
  }

  bool ThermalParticleSystem::CheckDecayChargesConservation(int ind) const
  {
    const ThermalParticle &part = Particles()[ind];
    int goalB = part.BaryonCharge();
    int goalQ = part.ElectricCharge();
    int goalS = part.Strangeness();
    int goalC = part.Charm();

    std::map<int, int> tPDGtoID = m_PDGtoID;

    for (int i = 0; i < part.Decays().size(); ++i) {
      int decB = 0, decQ = 0, decS = 0, decC = 0;
      for (int j = 0; j < part.Decays()[i].mDaughters.size(); ++j) {
        int tpdg = part.Decays()[i].mDaughters[j];
        if (tPDGtoID.count(tpdg) != 0) {
          int tid = tPDGtoID[tpdg];
          decB += Particles()[tid].BaryonCharge();
          decQ += Particles()[tid].ElectricCharge();
          decS += Particles()[tid].Strangeness();
          decC += Particles()[tid].Charm();
        }
      }
      if (goalB != decB || goalQ != decQ || goalS != decS || goalC != decC)
        return false;
    }

    return true;
  }

  bool ThermalParticleSystem::operator==(const ThermalParticleSystem & rhs) const
  {
    bool ret = true;
    ret &= m_PDGtoID == rhs.m_PDGtoID;
    ret &= m_NumBaryons == rhs.m_NumBaryons;
    ret &= m_NumCharged == rhs.m_NumCharged;
    ret &= m_NumStrange == rhs.m_NumStrange;
    ret &= m_NumCharmed == rhs.m_NumCharmed;
    ret &= m_NumberOfParticles == rhs.m_NumberOfParticles;
    ret &= m_ResonanceWidthIntegrationType == rhs.m_ResonanceWidthIntegrationType;
    ret &= m_DecayDistributionsMap == rhs.m_DecayDistributionsMap;

    ret &= m_Particles == rhs.m_Particles;

    return ret;
  }

  ParticleDecayType::DecayType ThermalParticleSystem::DecayTypeByParticleType(const ThermalParticle &part)
  {
    // Check if it's a known stable (wrt to any time scale of relevance) particle
    set<int> stablePDG;
    stablePDG.insert(2212); stablePDG.insert(-2212); // proton
    stablePDG.insert(2112); stablePDG.insert(-2112); // neutron
    stablePDG.insert(1000010020); stablePDG.insert(-1000010020); // deuteron
    stablePDG.insert(1000020030); stablePDG.insert(-1000020030); // He3
    stablePDG.insert(1000010030); stablePDG.insert(-1000010030); // triton
    stablePDG.insert(1000020040); stablePDG.insert(-1000020040); // He4

    if (stablePDG.count(part.PdgId()) > 0) {
      return ParticleDecayType::Stable;
    }
    
    // Check if it's a known weakly decaying particle
    set<int> weakPDG;
    weakPDG.insert(310); // K0S
    weakPDG.insert(130); // K0L
    weakPDG.insert(211); weakPDG.insert(-211);    // pi+-
    weakPDG.insert(321); weakPDG.insert(-321);    // K+-
    weakPDG.insert(3122); weakPDG.insert(-3122);  // Lambda
    weakPDG.insert(3222); weakPDG.insert(-3222);  // Sigma+
    weakPDG.insert(3112); weakPDG.insert(-3112);  // Sigma-
    weakPDG.insert(3322); weakPDG.insert(-3322);  // Ksi0
    weakPDG.insert(3312); weakPDG.insert(-3312);  // Ksi-
    weakPDG.insert(3334); weakPDG.insert(-3334);  // Omega
    weakPDG.insert(411); weakPDG.insert(-411);    // D+-
    weakPDG.insert(421); weakPDG.insert(-421);    // D0
    weakPDG.insert(431); weakPDG.insert(-431);    // Ds
    weakPDG.insert(4232); weakPDG.insert(-4232);  // Ksic+
    weakPDG.insert(4132); weakPDG.insert(-4132);  // Ksic0
    weakPDG.insert(4422); weakPDG.insert(-4422);  // Ksicc++
    weakPDG.insert(4412); weakPDG.insert(-4412);  // Ksicc+
    weakPDG.insert(4332); weakPDG.insert(-4332);  // Omegac

    if (weakPDG.count(part.PdgId()) > 0) {
      return ParticleDecayType::Weak;
    }

    // Check if it's a known electromagnetically decaying particle
    set<int> emPDG;
    emPDG.insert(111); // pi0
    emPDG.insert(221); // eta
    emPDG.insert(331); // eta'
    emPDG.insert(3212); emPDG.insert(-3212); // Sigma0

    if (emPDG.count(part.PdgId()) > 0) {
      return ParticleDecayType::Electromagnetic;
    }

    if (!part.IsStable()) // if particle not marked as stable assume it at least decays strongly
    {
      return ParticleDecayType::Strong;
    }
    else {
      // if contains strangeness (or charm), not stable under weak decays
      if (part.AbsoluteStrangeness() != 0 || part.AbsoluteCharm() != 0)
        return ParticleDecayType::Weak;

      return ParticleDecayType::Stable;
    }


    return ParticleDecayType::Default;
  }

  void ThermalParticleSystem::FillResonanceDecaysByFeeddown() {
    m_DecayContributionsByFeeddown[Feeddown::Weak].resize(m_Particles.size());
    m_DecayContributionsByFeeddown[Feeddown::Electromagnetic].resize(m_Particles.size());
    m_DecayContributionsByFeeddown[Feeddown::Strong].resize(m_Particles.size());
    for (int i = 0; i < m_Particles.size(); ++i) {
      m_DecayContributionsByFeeddown[Feeddown::Weak][i].resize(0);
      m_DecayContributionsByFeeddown[Feeddown::Electromagnetic][i].resize(0);
      m_DecayContributionsByFeeddown[Feeddown::Strong][i].resize(0);
    }
    for (int i = m_Particles.size() - 1; i >= 0; i--)
      if (m_Particles[i].DecayType() != ParticleDecayType::Stable && m_Particles[i].DecayType() != ParticleDecayType::Default) {
        GoResonanceByFeeddown(i, i, 1., Feeddown::Type(static_cast<int>(m_Particles[i].DecayType())));
      }
  }

  void ThermalParticleSystem::GoResonanceByFeeddown(int ind, int startind, double BR, Feeddown::Type feeddown) {
    for (int feed_index = static_cast<int>(Feeddown::Weak); feed_index <= static_cast<int>(Feeddown::Strong); ++feed_index) {
      if (static_cast<int>(feeddown) < feed_index) continue;
      //std::vector< std::pair<double, int> >& decayContributions = m_Particles[ind].DecayContributionsByFeeddown()[feed_index];
      DecayContributionsToParticle& decayContributions = m_DecayContributionsByFeeddown[feed_index][ind];
      if (ind != startind && decayContributions.size() > 0 && decayContributions[decayContributions.size() - 1].second == startind)
      {
        decayContributions[decayContributions.size() - 1].first += BR;
      }
      else if (ind != startind) {
        decayContributions.push_back(make_pair(BR, startind));
      }
    }

   
    if (m_Particles[ind].DecayType() != ParticleDecayType::Stable && m_Particles[ind].DecayType() != ParticleDecayType::Default) {
      for (int i = 0; i < m_Particles[ind].Decays().size(); ++i) {
        const ParticleDecayChannel& decaychannel = m_Particles[ind].Decays()[i];
        double tbr = decaychannel.mBratio;

        if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && ind == startind)
          tbr = decaychannel.mBratioAverage;

        for (int j = 0; j < decaychannel.mDaughters.size(); ++j) {
          if (m_PDGtoID.count(decaychannel.mDaughters[j]) != 0)
            GoResonanceByFeeddown(m_PDGtoID[decaychannel.mDaughters[j]], startind, BR*tbr, Feeddown::Type(static_cast<int>(m_Particles[ind].DecayType())));
        }
      }
    }
  }

} // namespace thermalfist