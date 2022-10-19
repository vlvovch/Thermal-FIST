/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEV/ThermalModelEVCrosstermsLegacy.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"

#include <Eigen/Dense>

using namespace Eigen;

using namespace std;

namespace thermalfist {

  ThermalModelEVCrosstermsLegacy::ThermalModelEVCrosstermsLegacy(ThermalParticleSystem *TPS, const ThermalModelParameters& params) :
    ThermalModelBase(TPS, params)
  {
    m_densitiesid.resize(m_TPS->Particles().size());
    m_Volume = params.V;
    m_Ps.resize(m_TPS->Particles().size());
    FillVirial(std::vector<double>(m_TPS->Particles().size(), 0.));
    m_TAG = "ThermalModelEVCrosstermsLegacy";

    m_Ensemble = GCE;
    m_InteractionModel = CrosstermsEV;
  }

  ThermalModelEVCrosstermsLegacy::~ThermalModelEVCrosstermsLegacy(void)
  {
  }

  void ThermalModelEVCrosstermsLegacy::FillVirial(const std::vector<double> & ri) {
    if (ri.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillVirial(const std::vector<double> & ri): size %d of ri does not match number of hadrons %d in the list", m_TAG.c_str(), static_cast<int>(ri.size()), static_cast<int>(m_TPS->Particles().size()));
      return;
    }
    m_Virial.resize(m_TPS->Particles().size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_Virial[i].resize(m_TPS->Particles().size());
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j)
        m_Virial[i][j] = CuteHRGHelper::brr(ri[i], ri[j]);
    }

    // Correction for non-diagonal terms R1=R2 and R2=0
    std::vector< std::vector<double> > fVirialTmp = m_Virial;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) {
        if (i == j) m_Virial[i][j] = fVirialTmp[i][j];
        else if ((fVirialTmp[i][i] + fVirialTmp[j][j]) > 0.0) m_Virial[i][j] = 2. * fVirialTmp[i][j] * fVirialTmp[i][i] / (fVirialTmp[i][i] + fVirialTmp[j][j]);
      }
  }

  void ThermalModelEVCrosstermsLegacy::ReadInteractionParameters(const std::string & filename)
  {
    m_Virial = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));

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
      double b;
      if (iss >> pdgid1 >> pdgid2 >> b) {
        int ind1 = m_TPS->PdgToId(pdgid1);
        int ind2 = m_TPS->PdgToId(pdgid2);
        if (ind1 != -1 && ind2 != -1)
          m_Virial[ind1][ind2] = b;
      }
    }
    fin.close();
  }

  void ThermalModelEVCrosstermsLegacy::WriteInteractionParameters(const std::string & filename)
  {
    ofstream fout(filename.c_str());
    fout << "# List of crossterms parameters to be used in the Crossterms excluded-volume HRG model"
      << std::endl;
    fout << "# Only particle pairs with a non-zero eigenvolume parameter are listed here"
      << std::endl;
    /*fout << "#" << std::setw(14) << "pdg_i"
      << std::setw(15) << "pdg_j"
      << std::setw(15) << "b_{ij}[fm^3]"
      << std::endl;*/
    fout << "#" << " " << "pdg_i"
      << " " << "pdg_j"
      << " " << "b_{ij}[fm^3]"
      << std::endl;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) {
        if (m_Virial[i][j] != 0.) {
          //fout << std::setw(15) << m_TPS->Particle(i).PdgId();
          //fout << std::setw(15) << m_TPS->Particle(j).PdgId();
          //fout << std::setw(15) << m_Virial[i][j];
          fout << " " << m_TPS->Particle(i).PdgId();
          fout << " " << m_TPS->Particle(j).PdgId();
          fout << " " << m_Virial[i][j];
          fout << std::endl;
        }
      }
    }
    fout.close();
  }

  void ThermalModelEVCrosstermsLegacy::SetRadius(double rad) {
    FillVirial(vector<double>(m_TPS->Particles().size(), rad));
  }

  double ThermalModelEVCrosstermsLegacy::VirialCoefficient(int i, int j) const {
    if (i < 0 || i >= static_cast<int>(m_Virial.size()) || j < 0 || j > static_cast<int>(m_Virial.size()))
      return 0.;
    return m_Virial[i][j];
  }

  void ThermalModelEVCrosstermsLegacy::SetVirial(int i, int j, double b) {
    if (i >= 0 && i < static_cast<int>(m_Virial.size()) && j >= 0 && j < static_cast<int>(m_Virial.size())) 
      m_Virial[i][j] = b;
    else printf("**WARNING** Index overflow in ThermalModelEVCrosstermsLegacy::SetVirial\n");
  }

  void ThermalModelEVCrosstermsLegacy::ChangeTPS(ThermalParticleSystem *TPS) {
    ThermalModelBase::ChangeTPS(TPS);
    m_densitiesid.resize(m_TPS->Particles().size());
    FillVirial();
  }

  double ThermalModelEVCrosstermsLegacy::DensityId(int i, const std::vector<double>& pstars)
  {
    double dMu = 0.;
    if (pstars.size() == m_TPS->Particles().size())
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) 
        dMu += -m_Virial[i][j] * pstars[j];
    else
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) 
        dMu += -m_Virial[i][j] * m_Ps[j];

    return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i] + dMu);
  }

  double ThermalModelEVCrosstermsLegacy::Pressure(int i, const std::vector<double>& pstars)
  {
    double dMu = 0.;
    if (pstars.size() == m_TPS->Particles().size())
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) 
        dMu += -m_Virial[i][j] * pstars[j];
    else
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) 
        dMu += -m_Virial[i][j] * m_Ps[j];

    return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i] + dMu);
  }

  double ThermalModelEVCrosstermsLegacy::ScaledVarianceId(int i) {
    double dMu = 0.;
    for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];

    return m_TPS->Particles()[i].ScaledVariance(m_Parameters, m_UseWidth, m_Chem[i] + dMu);
  }

  double ThermalModelEVCrosstermsLegacy::PartialPressureDiagonal(int i, double P) {
    double dMu = 0.;
    dMu += -m_Virial[i][i] * P;

    return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i] + dMu);
  }


  double ThermalModelEVCrosstermsLegacy::PressureDiagonalTotal(double P) {
    double ret = 0.;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += PartialPressureDiagonal(i, P);
    return ret;
  }

  void ThermalModelEVCrosstermsLegacy::SolveDiagonal() {
    BroydenEquationsCRSDEV eqs(this);
    BroydenJacobian jac(&eqs);
    jac.SetDx(1.0E-8);
    Broyden broydn(&eqs, &jac);
    Broyden::BroydenSolutionCriterium crit(1.0E-8);

    m_Pressure = 0.;
    double x0 = m_Pressure;
    std::vector<double> x(1, x0);

    x = broydn.Solve(x, &crit);

    m_Pressure = x[0];
    for (size_t i = 0; i < m_Ps.size(); ++i)
      m_Ps[i] = PartialPressureDiagonal(i, m_Pressure);
  }


  void ThermalModelEVCrosstermsLegacy::SolvePressure(bool resetPartials) {
    if (resetPartials) {
      m_Ps.resize(m_TPS->Particles().size());
      for (size_t i = 0; i < m_Ps.size(); ++i) m_Ps[i] = 0.;
      SolveDiagonal();
    }

    BroydenEquationsCRS eqs(this);
    BroydenJacobianCRS  jac(this);
    Broyden broydn(&eqs, &jac);
    BroydenSolutionCriteriumCRS crit(this);

    m_Ps = broydn.Solve(m_Ps, &crit);
    m_Pressure = 0.;
    for (size_t i = 0; i < m_Ps.size(); ++i) 
      m_Pressure += m_Ps[i];

    if (broydn.Iterations() == broydn.MaxIterations())
      m_LastCalculationSuccessFlag = false;
    else m_LastCalculationSuccessFlag = true;

    m_MaxDiff = broydn.MaxDifference();
  }

  void ThermalModelEVCrosstermsLegacy::CalculatePrimordialDensities() {
    m_FluctuationsCalculated = false;

    map< vector<double>, int> m_MapEVcomponent;

    {
      int NN = m_densities.size();
      m_MapToEVComponent.resize(NN);
      m_MapFromEVComponent.clear();
      m_MapEVcomponent.clear();
      m_EVComponentIndices.clear();

      int tind = 0;
      for (int i = 0; i < NN; ++i) {
        vector<double> EVParam(0);
        for (int j = 0; j < NN; ++j) {
          EVParam.push_back(m_Virial[i][j]);
          EVParam.push_back(m_Virial[j][i]);
        }

        if (m_MapEVcomponent.count(EVParam) == 0) {
          m_MapEVcomponent[EVParam] = tind;
          m_MapToEVComponent[i] = tind;
          m_MapFromEVComponent.push_back(i);
          m_EVComponentIndices.push_back(vector<int>(1, i));
          tind++;
        }
        else {
          m_MapToEVComponent[i] = m_MapEVcomponent[EVParam];
          m_EVComponentIndices[m_MapEVcomponent[EVParam]].push_back(i);
        }
      }
    }

    // Pressure
    SolvePressure();
    vector<double> tN(m_densities.size());
    for (size_t i = 0; i < m_Ps.size(); ++i) 
      tN[i] = DensityId(i);

    // Densities

    int NN = m_densities.size();

    MatrixXd densMatrix(NN, NN);
    VectorXd solVector(NN), xVector(NN);

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = m_Virial[i][j] * tN[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    PartialPivLU<MatrixXd> decomp(densMatrix);

    for (int i = 0; i < NN; ++i) xVector[i] = 0.;
    for (int i = 0; i < NN; ++i) {
      xVector[i] = tN[i];
      solVector = decomp.solve(xVector);
      if (1) {
        m_densities[i] = 0.;
        for (int j = 0; j < NN; ++j)
          m_densities[i] += solVector[j];
      }
      else {
        cout << "Could not recover m_densities from partial pressures!\n";
      }
      xVector[i] = 0.;
    }

    std::vector<double> fEntropyP(m_densities.size());
    for (int i = 0; i < NN; ++i) {
      double dMu = 0.;
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];
      xVector[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i] + dMu);
    }

    solVector = decomp.solve(xVector);

    if (1) {
      m_TotalEntropyDensity = 0.;
      for (int i = 0; i < NN; ++i)
        m_TotalEntropyDensity += solVector[i];
    }
    else {
      cout << "**ERROR** Could not recover m_densities from partial pressures!\n";
      return;
    }

    m_Calculated = true;
    ValidateCalculation();
  }

  void ThermalModelEVCrosstermsLegacy::CalculatePrimordialDensitiesNoReset() {
    // Pressure
    SolvePressure(false);
    vector<double> tN(m_densities.size());
    for (size_t i = 0; i < m_Ps.size(); ++i) 
      tN[i] = DensityId(i);

    // Densities

    int NN = m_densities.size();

    MatrixXd densMatrix(NN, NN);
    VectorXd solVector(NN), xVector(NN);

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = m_Virial[i][j] * tN[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    PartialPivLU<MatrixXd> decomp(densMatrix);

    for (int i = 0; i < NN; ++i) xVector[i] = 0.;
    for (int i = 0; i < NN; ++i) {
      xVector[i] = tN[i];//m_Ps[i] / m_Parameters.T;
      //solVector = lu.Solve(xVector, ok);
      solVector = decomp.solve(xVector);
      //if (ok) {
      if (1) {
        //if (decomp.info()==Eigen::Success) {
        m_densities[i] = 0.;
        for (int j = 0; j < NN; ++j)
          m_densities[i] += solVector[j];
      }
      else {
        cout << "**ERROR** Could not recover m_densities from partial pressures!\n";
        return;
      }
      xVector[i] = 0.;
    }

    std::vector<double> fEntropyP(m_densities.size());
    for (int i = 0; i < NN; ++i) {
      double dMu = 0.;
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];
      xVector[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i] + dMu);
    }

    solVector = decomp.solve(xVector);

    //if (ok) {
    if (1) {
      //if (decomp.info()==Eigen::Success) {
      m_TotalEntropyDensity = 0.;
      for (int i = 0; i < NN; ++i)
        m_TotalEntropyDensity += solVector[i];
    }
    else {
      cout << "**ERROR** Could not recover m_densities from partial pressures!\n";
      return;
    }

    // Decays

    CalculateFeeddown();

    m_Calculated = true;
  }

  void ThermalModelEVCrosstermsLegacy::SolvePressureIter() {
    m_Ps.resize(m_TPS->Particles().size());
    for (size_t i = 0; i < m_Ps.size(); ++i) 
      m_Ps[i] = 0.;
    SolveDiagonal();
    vector<double> Pstmp = m_Ps;
    int iter = 0;
    double maxdiff = 0.;
    for (iter = 0; iter < 1000; ++iter)
    {
      maxdiff = 0.;
      for (size_t i = 0; i < m_Ps.size(); ++i) {
        Pstmp[i] = Pressure(i);
        maxdiff = max(maxdiff, fabs((Pstmp[i] - m_Ps[i]) / Pstmp[i]));
      }
      m_Ps = Pstmp;
      //cout << iter << "\t" << maxdiff << "\n";
      if (maxdiff < 1.e-10) break;
    }
    if (iter == 1000) cout << iter << "\t" << maxdiff << "\n";
    m_Pressure = 0.;
    for (size_t i = 0; i < m_Ps.size(); ++i) 
      m_Pressure += m_Ps[i];
  }

  void ThermalModelEVCrosstermsLegacy::CalculatePrimordialDensitiesIter() {
    // Pressure
    SolvePressureIter();

    int NN = m_densities.size();
    vector<double> tN(NN);
    for (int i = 0; i < NN; ++i) tN[i] = DensityId(i);

    MatrixXd densMatrix(NN, NN);
    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = m_Virial[j][i] * tN[i];
        if (i == j) densMatrix(i, j) += 1.;
      }
    //densMatrix.SetMatrixArray(matr.GetArray());

    VectorXd solVector(NN), xVector(NN);
    for (int i = 0; i < NN; ++i) xVector[i] = tN[i];

    PartialPivLU<MatrixXd> decomp(densMatrix);

    solVector = decomp.solve(xVector);

    if (1) {
      //if (decomp.info()==Eigen::Success) {
      for (int i = 0; i < NN; ++i)
        m_densities[i] = solVector[i];
    }
    else {
      cout << "**ERROR** Could not recover m_densities from partial pressures!\n";
      return;
    }

    // Entropy
    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = m_Virial[i][j] * tN[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    std::vector<double> fEntropyP(m_densities.size());
    for (int i = 0; i < NN; ++i) {
      double dMu = 0.;
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];
      xVector[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_Chem[i] + dMu);
    }

    decomp = PartialPivLU<MatrixXd>(densMatrix);
    solVector = decomp.solve(xVector);

    if (1) {
      //if (decomp.info()==Eigen::Success) {
      m_TotalEntropyDensity = 0.;
      for (int i = 0; i < NN; ++i)
        m_TotalEntropyDensity += solVector[i];
    }
    else {
      cout << "Could not recover entropy m_densities from partial pressures!\n";
    }

    m_Calculated = true;
  }

  void ThermalModelEVCrosstermsLegacy::CalculateTwoParticleCorrelations() {
    int NN = m_densities.size();
    vector<double> tN(NN);// , tW(NN);
    for (int i = 0; i < NN; ++i) tN[i] = DensityId(i);
    //for (int i = 0; i < NN; ++i) tW[i] = ScaledVarianceId(i);
    vector<double> chi2id(NN);
    for (int i = 0; i < NN; ++i) {
      double dMu = 0.;
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) dMu += -m_Virial[i][j] * m_Ps[j];
      chi2id[i] = m_densities[i] / tN[i] * m_TPS->Particles()[i].chiDimensionfull(2, m_Parameters, m_UseWidth, m_Chem[i] + dMu) * xMath::GeVtoifm3();
    }

    MatrixXd densMatrix(NN, NN);
    VectorXd solVector(NN), xVector(NN), xVector2(NN);

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = m_Virial[i][j] * tN[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    PartialPivLU<MatrixXd> decomp(densMatrix);

    vector< vector<double> > ders, coefs;

    ders.resize(NN);
    coefs.resize(NN);

    for (int i = 0; i < NN; ++i) {
      ders[i].resize(NN);
      coefs[i].resize(NN);
    }

    for (int i = 0; i < NN; ++i) xVector[i] = 0.;
    for (int i = 0; i < NN; ++i) {
      xVector[i] = tN[i];
      solVector = decomp.solve(xVector);
      if (1) {
        //if (decomp.info()==Eigen::Success) {
        for (int j = 0; j < NN; ++j) {
          ders[j][i] = solVector[j];
        }

        for (int l = 0; l < NN; ++l) {
          coefs[l][i] = 0.;
          for (int k = 0; k < NN; ++k) {
            coefs[l][i] += -m_Virial[l][k] * ders[k][i];
          }
          if (l == i) coefs[l][i] += 1.;
        }
      }
      else {
        cout << "**WARNING** Could not recover fluctuations!\n";
      }
      xVector[i] = 0.;
    }


    m_PrimCorrel.resize(NN);
    for (int i = 0; i < NN; ++i) m_PrimCorrel[i].resize(NN);
    m_TotalCorrel = m_PrimCorrel;

    for (int i = 0; i < NN; ++i)
      for (int j = i; j < NN; ++j) {
        for (int l = 0; l < NN; ++l)
          //xVector[l] = tN[l] / m_Parameters.T * tW[l] * coefs[l][i] * coefs[l][j];
          xVector[l] = chi2id[l] * coefs[l][i] * coefs[l][j];
        solVector = decomp.solve(xVector);
        if (1) {
          //if (decomp.info()==Eigen::Success) {
          m_PrimCorrel[i][j] = 0.;
          for (int k = 0; k < NN; ++k) {
            m_PrimCorrel[i][j] += solVector[k];
          }
          m_PrimCorrel[j][i] = m_PrimCorrel[i][j];
        }
        else {
          cout << "**WARNING** Could not recover fluctuations!\n";
        }
      }

    //cout << "Primaries solved!\n";

    for (int i = 0; i < NN; ++i) {
      m_wprim[i] = m_PrimCorrel[i][i];
      if (m_densities[i] > 0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
      else m_wprim[i] = 1.;
    }

    
  }

  // TODO include correlations
  void ThermalModelEVCrosstermsLegacy::CalculateFluctuations() {
    CalculateTwoParticleCorrelations();
    CalculateSusceptibilityMatrix();
    CalculateTwoParticleFluctuationsDecays();
    CalculateProxySusceptibilityMatrix();
    CalculateParticleChargeCorrelationMatrix();

    m_FluctuationsCalculated = true;

    for (size_t i = 0; i < m_wprim.size(); ++i) {
      m_skewprim[i] = 1.;
      m_kurtprim[i] = 1.;
      m_skewtot[i] = 1.;
      m_kurttot[i] = 1.;
    }
  }

  std::vector<double> ThermalModelEVCrosstermsLegacy::CalculateChargeFluctuations(const std::vector<double>& chgs, int order)
  {
    vector<double> ret(order + 1, 0.);

    // chi1
    for (size_t i = 0; i < m_densities.size(); ++i)
      ret[0] += chgs[i] * m_densities[i];

    ret[0] /= pow(m_Parameters.T * xMath::GeVtoifm(), 3);

    if (order < 2) return ret;
    // Preparing matrix for system of linear equations
    int NN = m_densities.size();

    vector<double> MuStar(NN, 0.);
    for (int i = 0; i < NN; ++i) {
      MuStar[i] = m_Chem[i] + MuShift(i);
    }

    MatrixXd densMatrix(2 * NN, 2 * NN);
    VectorXd solVector(2 * NN), xVector(2 * NN);

    vector<double> DensitiesId(m_densities.size()), chi2id(m_densities.size());
    for (int i = 0; i < NN; ++i) {
      DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, MuStar[i]);
      chi2id[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, MuStar[i]);
    }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = m_Virial[j][i] * DensitiesId[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j)
        densMatrix(i, NN + j) = 0.;

    for (int i = 0; i < NN; ++i) {
      densMatrix(i, NN + i) = 0.;
      for (int k = 0; k < NN; ++k) {
        densMatrix(i, NN + i) += m_Virial[k][i] * m_densities[k];
      }
      densMatrix(i, NN + i) = (densMatrix(i, NN + i) - 1.) * chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T;
    }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(NN + i, NN + j) = m_Virial[i][j] * DensitiesId[j];
        if (i == j) densMatrix(NN + i, NN + j) += 1.;
      }


    PartialPivLU<MatrixXd> decomp(densMatrix);

    // chi2
    vector<double> dni(NN, 0.), dmus(NN, 0.);

    for (int i = 0; i < NN; ++i) {
      xVector[i] = 0.;
      xVector[NN + i] = chgs[i];
    }

    solVector = decomp.solve(xVector);

    for (int i = 0; i < NN; ++i) {
      dni[i] = solVector[i];
      dmus[i] = solVector[NN + i];
    }

    for (int i = 0; i < NN; ++i)
      ret[1] += chgs[i] * dni[i];

    ret[1] /= pow(m_Parameters.T, 2) * pow(xMath::GeVtoifm(), 3);

    if (order < 3) return ret;
    // chi3
    vector<double> d2ni(NN, 0.), d2mus(NN, 0.);

    vector<double> chi3id(m_densities.size());
    for (int i = 0; i < NN; ++i)
      chi3id[i] = m_TPS->Particles()[i].chi(3, m_Parameters, m_UseWidth, MuStar[i]);

    for (int i = 0; i < NN; ++i) {
      xVector[i] = 0.;

      double tmp = 0.;
      for (int j = 0; j < NN; ++j) tmp += m_Virial[j][i] * dni[j];
      tmp = -2. * tmp * chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[i];
      xVector[i] += tmp;

      tmp = 0.;
      for (int j = 0; j < NN; ++j) tmp += m_Virial[j][i] * m_densities[j];
      tmp = -(tmp - 1.) * chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * dmus[i] * dmus[i];
      xVector[i] += tmp;
    }
    for (int i = 0; i < NN; ++i) {
      xVector[NN + i] = 0.;

      double tmp = 0.;
      for (int j = 0; j < NN; ++j) tmp += -m_Virial[i][j] * dmus[j] * chi2id[j] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[j];

      xVector[NN + i] = tmp;
    }

    solVector = decomp.solve(xVector);

    for (int i = 0; i < NN; ++i) {
      d2ni[i] = solVector[i];
      d2mus[i] = solVector[NN + i];
    }

    for (int i = 0; i < NN; ++i)
      ret[2] += chgs[i] * d2ni[i];

    ret[2] /= m_Parameters.T * pow(xMath::GeVtoifm(), 3);


    if (order < 4) return ret;

    // chi4
    vector<double> d3ni(NN, 0.), d3mus(NN, 0.);

    vector<double> chi4id(m_densities.size());
    for (int i = 0; i < NN; ++i)
      chi4id[i] = m_TPS->Particles()[i].chi(4, m_Parameters, m_UseWidth, MuStar[i]);

    vector<double> dnis(NN, 0.);
    for (int i = 0; i < NN; ++i) {
      dnis[i] = chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[i];
    }

    vector<double> d2nis(NN, 0.);
    for (int i = 0; i < NN; ++i) {
      d2nis[i] = chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * dmus[i] * dmus[i] +
        chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * d2mus[i];
    }

    for (int i = 0; i < NN; ++i) {
      xVector[i] = 0.;

      double tmp = 0.;
      for (int j = 0; j < NN; ++j) tmp += m_Virial[j][i] * dni[j];
      tmp = -3. * tmp * d2nis[i];
      xVector[i] += tmp;

      tmp = 0.;
      for (int j = 0; j < NN; ++j) tmp += m_Virial[j][i] * d2ni[j];
      tmp = -3. * tmp * dnis[i];
      xVector[i] += tmp;

      double tmps = 0.;
      for (int j = 0; j < NN; ++j) tmps += m_Virial[j][i] * m_densities[j];

      tmp = -(tmps - 1.) * chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * d2mus[i] * 3. * dmus[i];
      xVector[i] += tmp;

      tmp = -(tmps - 1.) * chi4id[i] * pow(xMath::GeVtoifm(), 3) * dmus[i] * dmus[i] * dmus[i];
      xVector[i] += tmp;
    }
    for (int i = 0; i < NN; ++i) {
      xVector[NN + i] = 0.;

      double tmp = 0.;
      for (int j = 0; j < NN; ++j) tmp += -2. * m_Virial[i][j] * d2mus[j] * dnis[j];
      xVector[NN + i] += tmp;

      tmp = 0.;
      for (int j = 0; j < NN; ++j) tmp += -m_Virial[i][j] * dmus[j] * d2nis[j];
      xVector[NN + i] += tmp;
    }

    solVector = decomp.solve(xVector);

    for (int i = 0; i < NN; ++i) {
      d3ni[i] = solVector[i];
      d3mus[i] = solVector[NN + i];
    }

    for (int i = 0; i < NN; ++i)
      ret[3] += chgs[i] * d3ni[i];

    ret[3] /= pow(xMath::GeVtoifm(), 3);

    return ret;
  }

  double ThermalModelEVCrosstermsLegacy::CalculateEnergyDensity() {
    double ret = 0.;
    ret += m_Parameters.T * CalculateEntropyDensity();
    ret += -CalculatePressure();
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      ret += m_Chem[i] * m_densities[i];
    return ret;
  }

  double ThermalModelEVCrosstermsLegacy::CalculateEntropyDensity() {
    if (!m_Calculated) CalculateDensities();
    return m_TotalEntropyDensity;
  }

  double ThermalModelEVCrosstermsLegacy::CalculatePressure() {
    if (!m_Calculated) CalculateDensities();
    return m_Pressure;
  }


  double ThermalModelEVCrosstermsLegacy::MuShift(int id) const
  {
    if (id >= 0. && id < static_cast<int>(m_Virial.size())) {
      double dMu = 0.;
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j)
        dMu += -m_Virial[id][j] * m_Ps[j];
      return dMu;
    }
    else
      return 0.0;
  }

  std::vector<double> ThermalModelEVCrosstermsLegacy::BroydenEquationsCRS::Equations(const std::vector<double>& x)
  {
    std::vector<double> ret(m_N);
    for (size_t i = 0; i < x.size(); ++i)
      ret[i] = x[i] - m_THM->Pressure(i, x);
    return ret;
  }

  std::vector<double> ThermalModelEVCrosstermsLegacy::BroydenJacobianCRS::Jacobian(const std::vector<double>& x)
  {
    int N = x.size();

    vector<double> tN(N);
    for (int i = 0; i < N; ++i) {
      tN[i] = m_THM->DensityId(i, x);
    }

    vector<double> ret(N*N, 0.);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        if (i == j) 
          ret[i*N + j] += 1.;
        ret[i*N + j] += m_THM->VirialCoefficient(i, j) * tN[i];
      }
    }

    return ret;
  }

  bool ThermalModelEVCrosstermsLegacy::BroydenSolutionCriteriumCRS::IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& /*xdelta*/) const
  {
    double maxdiff = 0.;
    for (size_t i = 0; i < x.size(); ++i) {
      maxdiff = std::max(maxdiff, fabs(f[i]) / x[i]);
    }
    return (maxdiff < m_MaximumError);
  }

  std::vector<double> ThermalModelEVCrosstermsLegacy::BroydenEquationsCRSDEV::Equations(const std::vector<double>& x)
  {
    std::vector<double> ret(1);
    ret[0] = x[0] - m_THM->PressureDiagonalTotal(x[0]);
    return ret;
  }
} // namespace thermalfist