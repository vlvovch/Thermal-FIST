/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2018-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGVDW/ThermalModelVDWCanonicalStrangeness.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"

using namespace std;

namespace thermalfist {

  ThermalModelVDWCanonicalStrangeness::ThermalModelVDWCanonicalStrangeness(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
    ThermalModelCanonicalStrangeness(TPS_, params)
  {
    m_modelVDW = NULL;
    m_MuStar.resize(m_densities.size());
    m_Virial.resize(m_densities.size(), vector<double>(m_densities.size(), 0.));
    m_Attr = m_Virial;
    m_TAG = "ThermalModelVDWCanonicalStrangeness";

    m_Ensemble = SCE;
    m_InteractionModel = QvdW;
  }

  ThermalModelVDWCanonicalStrangeness::~ThermalModelVDWCanonicalStrangeness(void)
  {
    ClearModelVDW();
  }

  void ThermalModelVDWCanonicalStrangeness::CalculateDensitiesGCE() {
    if (m_modelVDW == NULL)
      PrepareModelVDW();

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_densitiesGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_MuStar[i]);
    }

    m_GCECalculated = true;
  }

  void ThermalModelVDWCanonicalStrangeness::CalculateEnergyDensitiesGCE() {
    if (m_modelVDW == NULL)
      PrepareModelVDW();

    m_energydensitiesGCE.resize(m_densitiesGCE.size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_energydensitiesGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_MuStar[i]);
    }
  }

  void ThermalModelVDWCanonicalStrangeness::CalculatePressuresGCE()
  {
    if (m_modelVDW == NULL)
      PrepareModelVDW();

    m_pressuresGCE.resize(m_densitiesGCE.size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_pressuresGCE[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
    }
  }

  void ThermalModelVDWCanonicalStrangeness::CalculateSums(const std::vector<double>& Vcs)
  {
    if (!m_GCECalculated)
      CalculateDensitiesGCE();

    m_Zsum.resize(m_StrVals.size());

    m_partialS.resize(m_StrVals.size());
    vector<double> xi(3, 0.), yi(3, 0.);

    for (size_t i = 0; i < m_StrVals.size(); ++i) {
      m_partialS[i] = 0.;
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j)
        if (m_StrVals[i] == m_TPS->Particles()[j].Strangeness())
          m_partialS[i] += m_densitiesGCE[j] * Vcs[j];
    }

    for (int i = 0; i < 3; ++i) {
      xi[i] = 2. * sqrt(m_partialS[m_StrMap[i + 1]] * m_partialS[m_StrMap[-(i + 1)]]);
      yi[i] = sqrt(m_partialS[m_StrMap[i + 1]] / m_partialS[m_StrMap[-(i + 1)]]);
    }

    // TODO: Choose iters dynamically based on total strangeness
    int iters = 20;

    for (unsigned int i = 0; i < m_StrVals.size(); ++i) {
      double res = 0.;

      for (int m = -iters; m <= iters; ++m)
        for (int n = -iters; n <= iters; ++n) {
          double tmp = xMath::BesselI(abs(3 * m + 2 * n - m_StrVals[i]), xi[0]) *
            xMath::BesselI(abs(n), xi[1]) *
            xMath::BesselI(abs(m), xi[2]) *
            pow(yi[0], m_StrVals[i] - 3 * m - 2 * n) *
            pow(yi[1], n) *
            pow(yi[2], m);
          if (tmp != tmp) continue;
          res += tmp;
        }
      m_Zsum[i] = res;
    }
  }


  void ThermalModelVDWCanonicalStrangeness::CalculatePrimordialDensities() {
    m_FluctuationsCalculated = false;

    m_energydensitiesGCE.resize(0);
    m_pressuresGCE.resize(0);

    PrepareModelVDW();

    CalculateDensitiesGCE();

    std::vector<double> Vcs(m_TPS->Particles().size(), 0.);
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      Vcs[i] = m_Parameters.SVc * m_Suppression[i];

    CalculateSums(Vcs);

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
        m_densities[i] = m_Suppression[i] * (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_densitiesGCE[i];
    }

    double tP = CalculatePressure();
    //printf("VDW-SCE calculation: The following parameter must be much smaller than unity. Otherwise we are in trouble...\n");
    printf("%s%lf\n", "PS/P = ", (tP - m_PNS) / tP);

    CalculateFeeddown();

    m_Calculated = true;
    m_LastCalculationSuccessFlag = true;
  }



  double ThermalModelVDWCanonicalStrangeness::CalculateEnergyDensity() {
    double ret = 0.;

    if (!m_Calculated)
      CalculateDensities();

    if (m_energydensitiesGCE.size() == 0)
      CalculateEnergyDensitiesGCE();

    ret += m_modelVDW->CalculateEnergyDensity();

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].Strangeness() != 0)
        if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
          ret += m_Suppression[i] * (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_energydensitiesGCE[i];
    return ret;
  }


  double ThermalModelVDWCanonicalStrangeness::CalculateEntropyDensity() {
    double ret = 0.;

    if (m_energydensitiesGCE.size() == 0)
      CalculateEnergyDensitiesGCE();

    ret += m_modelVDW->CalculateEntropyDensity();

    ret += log(m_Zsum[m_StrMap[0]]) / m_Parameters.SVc;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].Strangeness() != 0)
        if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
          ret += m_Suppression[i] * (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * ((m_energydensitiesGCE[i] - (m_MuStar[i]) * m_densitiesGCE[i]) / m_Parameters.T);

    return ret;
  }


  double ThermalModelVDWCanonicalStrangeness::CalculatePressure() {
    double ret = 0.;
    if (m_pressuresGCE.size() == 0)
      CalculatePressuresGCE();

    ret += m_modelVDW->CalculatePressure();

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      if (m_TPS->Particles()[i].Strangeness() != 0)
        if (m_StrMap.count(-m_TPS->Particles()[i].Strangeness()))
          ret += (m_Zsum[m_StrMap[-m_TPS->Particles()[i].Strangeness()]] / m_Zsum[m_StrMap[0]]) * m_pressuresGCE[i];
    return ret;
  }

  double ThermalModelVDWCanonicalStrangeness::MuShift(int id) const
  {
    if (id >= 0. && id < static_cast<int>(m_Virial.size()))
      return m_MuStar[id] - m_Chem[id];
    else
      return 0.0;
  }


  void ThermalModelVDWCanonicalStrangeness::PrepareModelVDW()
  {
    ClearModelVDW();

    ThermalParticleSystem *TPSnew = new ThermalParticleSystem(*m_TPS);

    // "Switch off" all strange particles
    for (size_t i = 0; i < TPSnew->Particles().size(); ++i) {
      ThermalParticle &part = TPSnew->Particle(i);
      if (part.Strangeness() != 0)
        part.SetDegeneracy(0.);
    }

    m_modelVDW = new ThermalModelVDWFull(TPSnew);
    m_modelVDW->FillVirialEV(m_Virial);
    m_modelVDW->FillAttraction(m_Attr);
    m_modelVDW->SetUseWidth(m_UseWidth);
    m_modelVDW->SetParameters(m_Parameters);
    m_modelVDW->SetVolume(m_Parameters.SVc);
    m_modelVDW->SetChemicalPotentials(m_Chem);

    m_modelVDW->CalculateDensities();

    std::vector<double> PidNS(m_modelVDW->Densities().size(), 0.);
    for (size_t j = 0; j < m_modelVDW->Densities().size(); ++j) {
      PidNS[j] = m_modelVDW->TPS()->Particles()[j].Density(m_modelVDW->Parameters(), IdealGasFunctions::Pressure, m_modelVDW->UseWidth(), m_modelVDW->MuStar(j));
    }

    m_PNS = m_modelVDW->CalculatePressure();

    m_Suppression.resize(TPS()->Particles().size());
    m_MuStar.resize(TPS()->Particles().size());
    for (size_t i = 0; i < m_Suppression.size(); ++i) {
      m_Suppression[i] = 1.0;
      for (size_t j = 0; j < m_modelVDW->Densities().size(); ++j) {
        m_Suppression[i] += -m_Virial[j][i] * m_modelVDW->Densities()[j];
      }

      m_MuStar[i] = m_Chem[i];
      for (size_t j = 0; j < m_modelVDW->Densities().size(); ++j) {
        m_MuStar[i] += -m_Virial[i][j] * PidNS[j];
        m_MuStar[i] += (m_Attr[i][j] + m_Attr[j][i]) * m_modelVDW->Densities()[j];
      }
    }
  }

  void ThermalModelVDWCanonicalStrangeness::ClearModelVDW()
  {
    if (m_modelVDW != NULL) {
      ThermalParticleSystem *TPSold = m_modelVDW->TPS();
      if (TPSold != NULL)
        delete TPSold;
      delete m_modelVDW;
      m_modelVDW = NULL;
    }
  }

  void ThermalModelVDWCanonicalStrangeness::FillVirial(const std::vector<double>& ri)
  {
    if (ri.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillVirial(const std::vector<double> & ri): size of ri does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_Virial.resize(m_TPS->Particles().size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_Virial[i].resize(m_TPS->Particles().size());
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j)
        m_Virial[i][j] = CuteHRGHelper::brr(ri[i], ri[j]);
    }

    // Correct R1=R2 and R2=0
    std::vector< std::vector<double> > fVirialTmp = m_Virial;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) {
        if (i == j) m_Virial[i][j] = fVirialTmp[i][j];
        else if ((fVirialTmp[i][i] + fVirialTmp[j][j]) > 0.0) m_Virial[i][j] = 2. * fVirialTmp[i][j] * fVirialTmp[i][i] / (fVirialTmp[i][i] + fVirialTmp[j][j]);
      }
  }

  void ThermalModelVDWCanonicalStrangeness::FillVirialEV(const std::vector< std::vector<double> > & bij)
  {
    if (bij.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillVirialEV(const std::vector<double> & bij): size of bij does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_Virial = bij;
  }

  void ThermalModelVDWCanonicalStrangeness::FillAttraction(const std::vector< std::vector<double> > & aij)
  {
    if (aij.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillAttraction(const std::vector<double> & aij): size of aij does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_Attr = aij;
  }

  void ThermalModelVDWCanonicalStrangeness::ReadInteractionParameters(const std::string & filename)
  {
    m_Virial = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));
    m_Attr = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));

    ifstream fin(filename.c_str());
    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1)
        continue;
      istringstream iss(elems[0]);
      int pdgid1, pdgid2;
      double b, a;
      if (iss >> pdgid1 >> pdgid2 >> b) {
        if (!(iss >> a))
          a = 0.;
        int ind1 = m_TPS->PdgToId(pdgid1);
        int ind2 = m_TPS->PdgToId(pdgid2);
        if (ind1 != -1 && ind2 != -1) {
          m_Virial[ind1][ind2] = b;
          m_Attr[ind1][ind2] = a;
        }
      }
    }
    fin.close();
  }

  void ThermalModelVDWCanonicalStrangeness::WriteInteractionParameters(const std::string & filename)
  {
    ofstream fout(filename.c_str());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) {
        fout << std::setw(15) << m_TPS->Particle(i).PdgId();
        fout << std::setw(15) << m_TPS->Particle(j).PdgId();
        fout << std::setw(15) << m_Virial[i][j];
        fout << std::setw(15) << m_Attr[i][j];
        fout << std::endl;
      }
    }
    fout.close();
  }

  double ThermalModelVDWCanonicalStrangeness::VirialCoefficient(int i, int j) const
  {
    if (i < 0 || i >= static_cast<int>(m_Virial.size()) || j < 0 || j >= static_cast<int>(m_Virial.size()))
      return 0.;
    return m_Virial[i][j];
  }

  double ThermalModelVDWCanonicalStrangeness::AttractionCoefficient(int i, int j) const
  {
    if (i < 0 || i >= static_cast<int>(m_Attr.size()) || j < 0 || j >= static_cast<int>(m_Attr.size()))
      return 0.;
    return m_Attr[i][j];
  }

  bool ThermalModelVDWCanonicalStrangeness::IsConservedChargeCanonical(ConservedCharge::Name charge) const {
    return (charge == ConservedCharge::StrangenessCharge);
  }

} // namespace thermalfist