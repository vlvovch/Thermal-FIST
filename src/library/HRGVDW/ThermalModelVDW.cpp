/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2016-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGVDW/ThermalModelVDW.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


#include <Eigen/Dense>

using namespace Eigen;

using namespace std;

namespace thermalfist {

  ThermalModelVDW::ThermalModelVDW(ThermalParticleSystem *TPS_, const ThermalModelParameters& params):
      ThermalModelBase(TPS_, params), m_SearchMultipleSolutions(false), m_TemperatureDependentAB(false), m_VDWComponentMapCalculated(false)
  {
    m_chi.resize(6);
    for(int i=0;i<6;++i) m_chi[i].resize(3);
    m_chiarb.resize(6);
    m_DensitiesId.resize(m_densities.size());
    m_MuStar.resize(m_densities.size());
    m_Virial.resize(m_densities.size(), vector<double>(m_densities.size(),0.));
    m_Attr = m_Virial;
    m_VirialdT = m_Virial;
    m_AttrdT   = m_AttrdT;
    m_Volume = params.V;
    m_TAG = "ThermalModelVDW";

    m_Ensemble = GCE;
    m_InteractionModel = QvdW;
  }


  ThermalModelVDW::~ThermalModelVDW(void)
  {
  }

  void ThermalModelVDW::FillChemicalPotentials() {
    ThermalModelBase::FillChemicalPotentials();
    for(size_t i = 0; i < m_MuStar.size(); ++i) 
      m_MuStar[i] = m_Chem[i];
  }

  void ThermalModelVDW::SetChemicalPotentials(const vector<double>& chem)
  {
    ThermalModelBase::SetChemicalPotentials(chem);
    for (size_t i = 0; i < m_MuStar.size(); ++i)
      m_MuStar[i] = m_Chem[i];
  }

  void ThermalModelVDW::FillVirial(const vector<double> & ri) {
    if (ri.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillVirial(const vector<double> & ri): size of ri does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_Virial.resize(m_TPS->Particles().size());
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_Virial[i].resize(m_TPS->Particles().size());
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j)
        m_Virial[i][j] = CuteHRGHelper::brr(ri[i], ri[j]);
    }

    // Correct R1=R2 and R2=0
    vector< vector<double> > fVirialTmp = m_Virial;
    for(int i=0;i<m_TPS->ComponentsNumber();++i)
      for(int j=0;j<m_TPS->ComponentsNumber();++j) {
        if (i==j) m_Virial[i][j] = fVirialTmp[i][j];
        else if ((fVirialTmp[i][i] + fVirialTmp[j][j]) > 0.0) m_Virial[i][j] = 2. * fVirialTmp[i][j] * fVirialTmp[i][i] / (fVirialTmp[i][i] + fVirialTmp[j][j]);
      }

    m_VirialdT = vector< vector<double> >(m_TPS->Particles().size(), vector<double>(m_TPS->Particles().size(),0.));
    m_AttrdT   = vector< vector<double> >(m_TPS->Particles().size(), vector<double>(m_TPS->Particles().size(), 0.));

    m_VDWComponentMapCalculated = false;
  }

  void ThermalModelVDW::FillVirialEV(const vector< vector<double> >& bij)
  {
    if (bij.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillVirialEV(const vector<double> & bij): size of bij does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_Virial = bij;

    m_VDWComponentMapCalculated = false;
  }

  void ThermalModelVDW::FillAttraction(const vector<vector<double> >& aij)
  {
    if (aij.size() != m_TPS->Particles().size()) {
      printf("**WARNING** %s::FillAttraction(const vector<double> & aij): size of aij does not match number of hadrons in the list", m_TAG.c_str());
      return;
    }
    m_Attr = aij;

    m_VDWComponentMapCalculated = false;
  }

  void ThermalModelVDW::ReadInteractionParameters(const string & filename)
  {
    m_Virial = vector< vector<double> >(m_TPS->Particles().size(), vector<double>(m_TPS->Particles().size(), 0.));
    m_Attr   = vector< vector<double> >(m_TPS->Particles().size(), vector<double>(m_TPS->Particles().size(), 0.));

    ifstream fin(filename.c_str());
    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1)
        continue;
      istringstream iss(elems[0]);
      long long pdgid1, pdgid2;
      double b, a;
      if (iss >> pdgid1 >> pdgid2 >> b) {
        if (!(iss >> a))
          a = 0.;
        int ind1 = m_TPS->PdgToId(pdgid1);
        int ind2 = m_TPS->PdgToId(pdgid2);
        if (ind1 != -1 && ind2 != -1) {
          m_Virial[ind1][ind2] = b;
          m_Attr[ind1][ind2]   = a;
        }
      }
    }
    fin.close();

    m_VDWComponentMapCalculated = false;
  }

  void ThermalModelVDW::WriteInteractionParameters(const string & filename)
  {
    ofstream fout(filename.c_str());
    fout << "# List of the van dar Waals interaction parameters to be used in the QvdW-HRG model"
      << std::endl;
    fout << "# Only particle pairs with a non-zero QvdW interaction are listed here"
      << std::endl;
    /*fout << "#" << std::setw(14) << "pdg_i"
      << std::setw(15) << "pdg_j"
      << std::setw(15) << "b_{ij}[fm^3]"
      << std::setw(20) << "a_{ij}[GeV*fm^3]"
      << std::endl;*/
    fout << "#" << " " << "pdg_i"
      << " " << "pdg_j"
      << " " << "b_{ij}[fm^3]"
      << " " << "a_{ij}[GeV*fm^3]"
      << std::endl;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) {
        if (m_Virial[i][j] != 0. || m_Attr[i][j] != 0.) {
          //fout << setw(15) << m_TPS->Particle(i).PdgId();
          //fout << setw(15) << m_TPS->Particle(j).PdgId();
          //fout << setw(15) << m_Virial[i][j];
          //fout << setw(20) << m_Attr[i][j];
          //fout << endl;
          fout << " " << m_TPS->Particle(i).PdgId();
          fout << " " << m_TPS->Particle(j).PdgId();
          fout << " " << m_Virial[i][j];
          fout << " " << m_Attr[i][j];
          fout << endl;
        }
      }
    }
    fout.close();
  }

  void ThermalModelVDW::ChangeTPS(ThermalParticleSystem *TPS_) {
      ThermalModelBase::ChangeTPS(TPS_);
      m_Virial = vector< vector<double> >(m_TPS->Particles().size(), vector<double>(m_TPS->Particles().size(), 0.));
      m_Attr = vector< vector<double> >(m_TPS->Particles().size(), vector<double>(m_TPS->Particles().size(), 0.));
      m_VirialdT = vector< vector<double> >(m_TPS->Particles().size(), vector<double>(m_TPS->Particles().size(), 0.));
      m_AttrdT = vector< vector<double> >(m_TPS->Particles().size(), vector<double>(m_TPS->Particles().size(), 0.));
      m_VDWComponentMapCalculated = false;
  }

  std::vector<double> ThermalModelVDW::ComputeNp(const std::vector<double>& dmustar)
  {
    vector<double> ns(m_densities.size(), 0.);
    for (size_t i = 0; i < m_densities.size(); ++i)
      ns[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i] + dmustar[m_MapTodMuStar[i]]);

    return ComputeNp(dmustar, ns);
  }

  std::vector<double> ThermalModelVDW::ComputeNp(const std::vector<double>& /*dmustar*/, const std::vector<double>& ns)
  {
    int NN = m_densities.size();
    int NNdmu = m_MapFromdMuStar.size();

    MatrixXd densMatrix(NNdmu, NNdmu);
    VectorXd solVector(NNdmu), xVector(NNdmu);

    for (int i = 0; i < NNdmu; ++i) {
      for (int j = 0; j < NNdmu; ++j) {
        densMatrix(i, j) = 0.;
        if (i == j)
          densMatrix(i, j) += 1.;

        for (size_t m = 0; m < m_dMuStarIndices[i].size(); ++m) {
          densMatrix(i, j) += m_Virial[m_MapFromdMuStar[j]][m_dMuStarIndices[i][m]] * ns[m_dMuStarIndices[i][m]];
        }
      }
    }

    PartialPivLU<MatrixXd> decomp(densMatrix);

    for (int kp = 0; kp < NNdmu; ++kp) {
      xVector[kp] = 0.;
      for (size_t m = 0; m < m_dMuStarIndices[kp].size(); ++m) {
        xVector[kp] += ns[m_dMuStarIndices[kp][m]];
      }
    }


    solVector = decomp.solve(xVector);

    vector<double> ntildek(NNdmu, 0.);
    for (int i = 0; i < NNdmu; ++i)
      ntildek[i] = solVector[i];

    vector<double> np(m_densities.size(), 0.);
    for (int i = 0; i < NN; ++i) {
      np[i] = 0.;
      for (int k = 0; k < NNdmu; ++k) {
        np[i] += m_Virial[m_MapFromdMuStar[k]][i] * solVector[k];
      }
      np[i] = (1. - np[i]) * ns[i];
    }

    return np;
  }

  void ThermalModelVDW::CalculateVDWComponentsMap()
  {
    map< vector<double>, int> MapVDW;

    int NN = m_densities.size();

    m_MapTodMuStar.resize(NN);
    m_MapFromdMuStar.clear();
    MapVDW.clear();
    m_dMuStarIndices.clear();

    int tind = 0;
    for (int i = 0; i < NN; ++i) {
      vector<double> VDWParam(0);
      for (int j = 0; j < NN; ++j) {
        VDWParam.push_back(m_Virial[i][j]);
      }
      for (int j = 0; j < NN; ++j) {
        VDWParam.push_back(m_Attr[i][j] + m_Attr[j][i]);
      }

      if (MapVDW.count(VDWParam) == 0) {
        MapVDW[VDWParam] = tind;
        m_MapTodMuStar[i] = tind;
        m_MapFromdMuStar.push_back(i);
        m_dMuStarIndices.push_back(vector<int>(1, i));
        tind++;
      }
      else {
        m_MapTodMuStar[i] = MapVDW[VDWParam];
        m_dMuStarIndices[MapVDW[VDWParam]].push_back(i);
      }
    }

    m_VDWComponentMapCalculated = true;
  }

  vector<double> ThermalModelVDW::SearchSingleSolution(const vector<double>& muStarInit)
  {
    int NN = m_densities.size();
    int NNdmu = m_MapFromdMuStar.size();

    vector<double> dmuscur(NNdmu, 0.);
    for (int i = 0; i < NNdmu; ++i)
      dmuscur[i] = muStarInit[m_MapFromdMuStar[i]] - m_Chem[m_MapFromdMuStar[i]];

    BroydenEquationsVDW eqs(this);
    BroydenJacobianVDW  jac(this);
    Broyden broydn(&eqs, &jac);
    BroydenSolutionCriteriumVDW crit(this);

    dmuscur = broydn.Solve(dmuscur, &crit);

    if (broydn.Iterations() == broydn.MaxIterations())
      m_LastBroydenSuccessFlag = false;
    else m_LastBroydenSuccessFlag = true;



    m_MaxDiff = broydn.MaxDifference();

    vector<double> ret(NN);
    for (int i = 0; i < NN; ++i)
      ret[i] = m_Chem[i] + dmuscur[m_MapTodMuStar[i]];

    return ret;
  }

  vector<double> ThermalModelVDW::SearchMultipleSolutions(int iters) {
    vector<double> csol(m_densities.size(), 0.);
    double Psol = 0.;
    bool solved = false;
    double muBmin = m_Parameters.muB - 0.5 * xMath::mnucleon();
    double muBmax = m_Parameters.muB + 0.5 * xMath::mnucleon();
    double dmu = (muBmax - muBmin) / iters;
    vector<double> curmust(m_densities.size(), 0.);
    double maxdif = 0.;
    for(int isol = 0; isol < iters; ++isol) {
      double tmu = muBmin + (0.5 + isol) * dmu;
      for(size_t j = 0; j < curmust.size(); ++j) {
        curmust[j] = m_Chem[j] + (tmu - m_Parameters.muB) * m_Chem[j] / m_Parameters.muB;
        if (m_TPS->Particles()[j].Statistics()==-1 && curmust[j] > m_TPS->Particles()[j].Mass()) 
          curmust[j] = 0.98 * m_TPS->Particles()[j].Mass();
      }

      vector<double> sol = SearchSingleSolution(curmust);

      bool fl = true;
      for(size_t i = 0; i < sol.size(); ++i)
        if (sol[i] != sol[i]) fl = false;
      fl &= m_LastBroydenSuccessFlag;
      if (!fl) continue;

      for(int i = 0; i < m_TPS->ComponentsNumber(); ++i) 
        m_DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, sol[i]);

      int NN = m_densities.size();

      MatrixXd densMatrix(NN, NN);
      VectorXd solVector(NN), xVector(NN);

      for(int i=0;i<NN;++i)
        for(int j=0;j<NN;++j) {
          densMatrix(i,j) = m_Virial[j][i] * m_DensitiesId[i];
          if (i==j) densMatrix(i,j) += 1.;
        }

      PartialPivLU<MatrixXd> decomp(densMatrix);

      for(int i=0;i<NN;++i) 
        xVector[i] = m_DensitiesId[i];
      solVector = decomp.solve(xVector);
      for(int i=0;i<NN;++i) 
        m_densities[i] = solVector[i];

      double tP = 0.;
      for(size_t i=0;i<m_densities.size();++i) 
        tP += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, sol[i]);
      for(size_t i=0;i<m_densities.size();++i)
        for(size_t j=0;j<m_densities.size();++j)
          tP += -m_Attr[i][j] * m_densities[i] * m_densities[j];

      if (!solved || tP > Psol) {
        solved = true;
        Psol = tP;
        csol = sol;
        maxdif = m_MaxDiff;
      }
    }
    m_LastBroydenSuccessFlag = solved;
    m_MaxDiff = maxdif;
    return csol;
  }

  void ThermalModelVDW::SolveEquations() {
    if (!m_SearchMultipleSolutions) {

      vector<double> muStarInit = m_MuStar;

      for(size_t i=0;i<muStarInit.size();++i) {
        if (m_TPS->Particles()[i].Statistics()==-1 && muStarInit[i] > m_TPS->Particles()[i].Mass()) 
          muStarInit[i] = 0.98 * m_TPS->Particles()[i].Mass();
      }
      m_MuStar = SearchSingleSolution(muStarInit);
    }
    else {
      m_MuStar = SearchMultipleSolutions(100);
    }
  }

  void ThermalModelVDW::CalculatePrimordialDensities() {
    CalculatePrimordialDensitiesNew();
    ValidateCalculation();
  }

  void ThermalModelVDW::CalculatePrimordialDensitiesOld() {
    m_FluctuationsCalculated = false;

    if (!m_VDWComponentMapCalculated)
      CalculateVDWComponentsMap();

    int NN = m_densities.size();
    
    //map< vector<double>, int> m_MapVDW;
    //{
    //  m_MapTodMuStar.resize(NN);
    //  m_MapFromdMuStar.clear();
    //  m_MapVDW.clear();
    //  m_dMuStarIndices.clear();

    //  int tind = 0;
    //  for (int i = 0; i < NN; ++i) {
    //    vector<double> VDWParam(0);
    //    for (int j = 0; j < NN; ++j) {
    //      VDWParam.push_back(m_Virial[i][j]);
    //    }
    //    for (int j = 0; j < NN; ++j) {
    //      VDWParam.push_back(m_Attr[i][j] + m_Attr[j][i]);
    //    }

    //    if (m_MapVDW.count(VDWParam) == 0) {
    //      m_MapVDW[VDWParam] = tind;
    //      m_MapTodMuStar[i]  = tind;
    //      m_MapFromdMuStar.push_back(i);
    //      m_dMuStarIndices.push_back(vector<int>(1, i));
    //      tind++;
    //    }
    //    else {
    //      m_MapTodMuStar[i] = m_MapVDW[VDWParam];
    //      m_dMuStarIndices[m_MapVDW[VDWParam]].push_back(i);
    //    }
    //  }

    //  printf("Optimization: %d --> %d\n", NN, static_cast<int>(m_MapFromdMuStar.size()));
    //}

    clock_t tbeg = clock();

    {

      vector<double> muStarInit = m_MuStar;

      for (size_t i = 0; i<muStarInit.size(); ++i) {
        if (m_TPS->Particles()[i].Statistics() == -1 && muStarInit[i] > m_TPS->Particles()[i].Mass())
          muStarInit[i] = 0.98 * m_TPS->Particles()[i].Mass();
      }


      m_MuStar = SearchSingleSolution(muStarInit);
    }
  

    printf("Solution time = %lf ms\n", (clock() - tbeg) / (double)(CLOCKS_PER_SEC) * 1.e3);

    tbeg = clock();

    for(int i=0;i<m_TPS->ComponentsNumber();++i) 
      m_DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_MuStar[i]);

    MatrixXd densMatrix(NN, NN);
    VectorXd solVector(NN), xVector(NN);

    for(int i=0;i<NN;++i)
      for(int j=0;j<NN;++j) {
        densMatrix(i,j) = m_Virial[j][i] * m_DensitiesId[i];
        if (i==j) densMatrix(i,j) += 1.;
      }

    PartialPivLU<MatrixXd> decomp(densMatrix);

    for(int i=0;i<NN;++i) xVector[i] = m_DensitiesId[i];
    solVector = decomp.solve(xVector);
    for(int i=0;i<NN;++i) m_densities[i] = solVector[i];

    // TODO: Scalar density properly
    for(int i=0;i<NN;++i) xVector[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ScalarDensity, m_UseWidth, m_MuStar[i]);
    solVector = decomp.solve(xVector);
    m_scaldens.resize(m_densities.size());
    for(int i=0;i<NN;++i) m_scaldens[i] = solVector[i];


    tbeg = clock();

    m_Calculated = true;
  }

  void ThermalModelVDW::CalculatePrimordialDensitiesNew() {
    m_FluctuationsCalculated = false;

    if (!m_VDWComponentMapCalculated)
      CalculateVDWComponentsMap();

    int NN = m_densities.size();

    clock_t tbeg = clock();

    SolveEquations();

    tbeg = clock();

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      m_DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_MuStar[i]);

  
    int NNdmu = m_MapFromdMuStar.size();

    MatrixXd densMatrix(NNdmu, NNdmu);
    VectorXd solVector(NNdmu), xVector(NNdmu);

    for (int i = 0; i < NNdmu; ++i) {
      for (int j = 0; j < NNdmu; ++j) {
        densMatrix(i, j) = 0.;
        if (i == j)
          densMatrix(i, j) += 1.;

        for (size_t m = 0; m < m_dMuStarIndices[i].size(); ++m) {
          densMatrix(i, j) += m_Virial[m_MapFromdMuStar[j]][m_dMuStarIndices[i][m]] * m_DensitiesId[m_dMuStarIndices[i][m]];
        }
      }
    }

    PartialPivLU<MatrixXd> decomp(densMatrix);

    for (int kp = 0; kp < NNdmu; ++kp) {
      xVector[kp] = 0.;
      for (size_t m = 0; m < m_dMuStarIndices[kp].size(); ++m) {
        xVector[kp] += m_DensitiesId[m_dMuStarIndices[kp][m]];
      }
    }

    solVector = decomp.solve(xVector);

    vector<double> ntildek(NNdmu, 0.);
    for (int i = 0; i < NNdmu; ++i)
      ntildek[i] = solVector[i];

    //vector<double> np(m_densities.size(), 0.);
    for (int i = 0; i < NN; ++i) {
      m_densities[i] = 0.;
      for (int k = 0; k < NNdmu; ++k) {
        m_densities[i] += m_Virial[m_MapFromdMuStar[k]][i] * solVector[k];
      }
      m_densities[i] = (1. - m_densities[i]) * m_DensitiesId[i];
    }

    // TODO: Scalar density properly
    m_scaldens = m_densities;

    m_Calculated = true;
  }

  vector<double> ThermalModelVDW::CalculateChargeFluctuations(const vector<double> &chgs, int order) {
    vector<double> ret(order + 1, 0.);
  
    // chi1
    for(size_t i=0;i<m_densities.size();++i)
      ret[0] += chgs[i] * m_densities[i];

    ret[0] /= pow(m_Parameters.T * xMath::GeVtoifm(), 3);

    if (order<2) return ret;
    // Preparing matrix for system of linear equations
    int NN = m_densities.size();
    MatrixXd densMatrix(2*NN, 2*NN);
    VectorXd solVector(2*NN), xVector(2*NN);

    vector<double> chi2id(m_densities.size());
    for(int i=0;i<NN;++i) 
      chi2id[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_MuStar[i]);

    for(int i=0;i<NN;++i)
      for(int j=0;j<NN;++j) {
        densMatrix(i,j) = m_Virial[j][i] * m_DensitiesId[i];
        if (i==j) densMatrix(i,j) += 1.;
      }

    for(int i=0;i<NN;++i)
      for(int j=0;j<NN;++j) 
        densMatrix(i,NN+j) = 0.;

    for(int i=0;i<NN;++i) {
      densMatrix(i,NN+i) = 0.;
      for(int k=0;k<NN;++k) {
        densMatrix(i,NN+i) += m_Virial[k][i] * m_densities[k];
      }
      densMatrix(i,NN+i) = (densMatrix(i,NN+i) - 1.) * chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T;
    }

    for(int i=0;i<NN;++i)
      for(int j=0;j<NN;++j) {
        densMatrix(NN+i,j) = -(m_Attr[i][j] + m_Attr[j][i]);
      }
  
    for(int i=0;i<NN;++i)
      for(int j=0;j<NN;++j) {
        densMatrix(NN+i,NN+j) = m_Virial[i][j] * m_DensitiesId[j];
        if (i==j) densMatrix(NN+i,NN+j) += 1.;
      }


    PartialPivLU<MatrixXd> decomp(densMatrix);

    // chi2
    vector<double> dni(NN, 0.), dmus(NN, 0.);

    for(int i=0;i<NN;++i) {
      xVector[i]    = 0.;
      xVector[NN+i] = chgs[i];
    }

    solVector = decomp.solve(xVector);

    for(int i=0;i<NN;++i) {
      dni[i]  = solVector[i];
      dmus[i] = solVector[NN+i];
    }

    for(int i=0;i<NN;++i)
      ret[1] += chgs[i] * dni[i];

    ret[1] /= pow(m_Parameters.T, 2) * pow(xMath::GeVtoifm(), 3);

    if (order<3) return ret;

    // chi3
    vector<double> d2ni(NN, 0.), d2mus(NN, 0.);

    vector<double> chi3id(m_densities.size());
    for(int i=0;i<NN;++i) 
      chi3id[i] = m_TPS->Particles()[i].chi(3, m_Parameters, m_UseWidth, m_MuStar[i]);

    for(int i=0;i<NN;++i) {
      xVector[i]    = 0.;

      double tmp = 0.;
      for(int j=0;j<NN;++j) tmp += m_Virial[j][i] * dni[j];
      tmp = -2. * tmp * chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[i];
      xVector[i] += tmp;

      tmp = 0.;
      for(int j=0;j<NN;++j) tmp += m_Virial[j][i] * m_densities[j];
      tmp = -(tmp - 1.) * chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * dmus[i] * dmus[i];
      xVector[i] += tmp;
    }
    for(int i=0;i<NN;++i) {
      xVector[NN+i]    = 0.;

      double tmp = 0.;
      for(int j=0;j<NN;++j) tmp += -m_Virial[i][j] * dmus[j] * chi2id[j] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[j];
    
      xVector[NN+i] = tmp;
    }

    solVector = decomp.solve(xVector);

    for(int i=0;i<NN;++i) {
      d2ni[i]  = solVector[i];
      d2mus[i] = solVector[NN+i];
    }

    for(int i=0;i<NN;++i)
      ret[2] += chgs[i] * d2ni[i];

    ret[2] /= m_Parameters.T * pow(xMath::GeVtoifm(), 3);


    if (order<4) return ret;

    // chi4
    vector<double> d3ni(NN, 0.), d3mus(NN, 0.);

    vector<double> chi4id(m_densities.size());
    for (int i = 0; i < NN; ++i)
      chi4id[i] = m_TPS->Particles()[i].chi(4, m_Parameters, m_UseWidth, m_MuStar[i]);

    vector<double> dnis(NN, 0.);
    for(int i=0;i<NN;++i) {
      dnis[i] = chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * dmus[i];
    }

    vector<double> d2nis(NN, 0.);
    for(int i=0;i<NN;++i) {
      d2nis[i] = chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * dmus[i] * dmus[i] + 
        chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T * d2mus[i];
    }

    for(int i=0;i<NN;++i) {
      xVector[i]    = 0.;

      double tmp = 0.;
      for(int j=0;j<NN;++j) tmp += m_Virial[j][i] * dni[j];
      tmp = -3. * tmp * d2nis[i];
      xVector[i] += tmp;

      tmp = 0.;
      for(int j=0;j<NN;++j) tmp += m_Virial[j][i] * d2ni[j];
      tmp = -3. * tmp * dnis[i];
      xVector[i] += tmp;

      double tmps = 0.;
      for(int j=0;j<NN;++j) tmps += m_Virial[j][i] * m_densities[j];

      tmp = -(tmps - 1.) * chi3id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * d2mus[i] * 3. * dmus[i];
      xVector[i] += tmp;

      tmp = -(tmps - 1.) * chi4id[i] * pow(xMath::GeVtoifm(), 3) * dmus[i] * dmus[i] * dmus[i];
      xVector[i] += tmp;
    }
    for(int i=0;i<NN;++i) {
      xVector[NN+i]    = 0.;

      double tmp = 0.;
      for(int j=0;j<NN;++j) tmp += -2. * m_Virial[i][j] * d2mus[j] * dnis[j];
      xVector[NN+i] += tmp;

      tmp = 0.;
      for(int j=0;j<NN;++j) tmp += -m_Virial[i][j] * dmus[j] * d2nis[j];
      xVector[NN+i] += tmp;
    }

    solVector = decomp.solve(xVector);

    for(int i=0;i<NN;++i) {
      d3ni[i]  = solVector[i];
      d3mus[i] = solVector[NN+i];
    }

    for(int i=0;i<NN;++i)
      ret[3] += chgs[i] * d3ni[i];

    ret[3] /= pow(xMath::GeVtoifm(), 3);

    return ret;
  }

  // TODO include correlations
  vector< vector<double> >  ThermalModelVDW::CalculateFluctuations(int order) {
    if (order<1) return m_chi;
  
    vector<double> chgs(m_densities.size());
    vector<double> chis;

    // Baryon charge
    for(size_t i=0;i<chgs.size();++i)
      chgs[i] = m_TPS->Particles()[i].BaryonCharge();
    chis = CalculateChargeFluctuations(chgs, order);
    for(int i=0;i<order;++i) m_chi[i][0] = chis[i];

    // Electric charge
    for(size_t i=0;i<chgs.size();++i)
      chgs[i] = m_TPS->Particles()[i].ElectricCharge();
    chis = CalculateChargeFluctuations(chgs, order);
    for(int i=0;i<order;++i) m_chi[i][1] = chis[i];

    // Strangeness charge
    for(size_t i=0;i<chgs.size();++i)
      chgs[i] = m_TPS->Particles()[i].Strangeness();
    chis = CalculateChargeFluctuations(chgs, order);
    for(int i=0;i<order;++i) m_chi[i][2] = chis[i];

    // Arbitrary charge
    for(size_t i=0;i<chgs.size();++i)
      chgs[i] = m_TPS->Particles()[i].ArbitraryCharge();
    chis = CalculateChargeFluctuations(chgs, order);
    for(int i=0;i<order;++i) m_chiarb[i] = chis[i];

    return m_chi;
  }

  void ThermalModelVDW::CalculateTwoParticleCorrelations()
  {
    int NN = m_densities.size();

    m_PrimCorrel.resize(NN);
    for (int i = 0; i < NN; ++i)
      m_PrimCorrel[i].resize(NN);
    m_TotalCorrel = m_PrimCorrel;

    MatrixXd densMatrix(2 * NN, 2 * NN);
    VectorXd solVector(2 * NN), xVector(2 * NN);

    vector<double> chi2id(m_densities.size());
    for (int i = 0; i<NN; ++i)
      //chi2id[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_MuStar[i]);
      chi2id[i] = m_TPS->Particles()[i].chiDimensionfull(2, m_Parameters, m_UseWidth, m_MuStar[i]);

    for (int i = 0; i<NN; ++i)
      for (int j = 0; j<NN; ++j) {
        densMatrix(i, j) = m_Virial[j][i] * m_DensitiesId[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    for (int i = 0; i<NN; ++i)
      for (int j = 0; j<NN; ++j)
        densMatrix(i, NN + j) = 0.;

    for (int i = 0; i<NN; ++i) {
      densMatrix(i, NN + i) = 0.;
      for (int k = 0; k<NN; ++k) {
        densMatrix(i, NN + i) += m_Virial[k][i] * m_densities[k];
      }
      //densMatrix(i, NN + i) = (densMatrix(i, NN + i) - 1.) * chi2id[i] * pow(xMath::GeVtoifm(), 3) * m_Parameters.T * m_Parameters.T;
      densMatrix(i, NN + i) = (densMatrix(i, NN + i) - 1.) * chi2id[i] * pow(xMath::GeVtoifm(), 3);
    }

    for (int i = 0; i<NN; ++i)
      for (int j = 0; j<NN; ++j) {
        densMatrix(NN + i, j) = -(m_Attr[i][j] + m_Attr[j][i]);
      }

    for (int i = 0; i<NN; ++i)
      for (int j = 0; j<NN; ++j) {
        densMatrix(NN + i, NN + j) = m_Virial[i][j] * m_DensitiesId[j];
        if (i == j) densMatrix(NN + i, NN + j) += 1.;
      }

    PartialPivLU<MatrixXd> decomp(densMatrix);

    for (int k = 0; k < NN; ++k) {
      vector<double> dni(NN, 0.), dmus(NN, 0.);

      for (int i = 0; i < NN; ++i) {
        xVector[i] = 0.;
        xVector[NN + i] = static_cast<int>(i == k);
      }

      solVector = decomp.solve(xVector);

      for (int i = 0; i < NN; ++i) {
        dni[i]  = solVector[i];
        dmus[i] = solVector[NN + i];
      }

      for (int j = 0; j < NN; ++j) {
        m_PrimCorrel[j][k] = dni[j];
      }
    }

    for (int i = 0; i < NN; ++i) {
      m_wprim[i] = m_PrimCorrel[i][i];
      if (m_densities[i] > 0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
      else m_wprim[i] = 1.;
    }
    
  }

  void ThermalModelVDW::CalculateFluctuations()
  {
    CalculateTwoParticleCorrelations();
    CalculateSusceptibilityMatrix();
    CalculateTwoParticleFluctuationsDecays();
    CalculateProxySusceptibilityMatrix();
    CalculateParticleChargeCorrelationMatrix();

    for (size_t i = 0; i < m_wprim.size(); ++i) {
      m_skewprim[i] = 1.;
      m_kurtprim[i] = 1.;
      m_skewtot[i] = 1.;
      m_kurttot[i] = 1.;
    }

    m_FluctuationsCalculated = true;
  }


  double ThermalModelVDW::CalculateEnergyDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for(size_t i=0;i<m_densities.size();++i) 
      if (m_densities[i]>0.) 
        ret += m_densities[i] / m_DensitiesId[i] * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_MuStar[i]);
    for(size_t i=0;i<m_densities.size();++i)
      for(size_t j=0;j<m_densities.size();++j)
        ret += -m_Attr[i][j] * m_densities[i] * m_densities[j];

    if (m_TemperatureDependentAB) {
      for (size_t i = 0; i < m_densities.size(); ++i) {
        double tPid = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
        for (size_t j = 0; j < m_densities.size(); ++j) {
          ret += -tPid * m_densities[j] * m_Parameters.T * m_VirialdT[j][i];
          ret += m_Parameters.T * m_AttrdT[i][j] * m_densities[i] * m_densities[j];
        }
      }
    }

    return ret;
  }

  double ThermalModelVDW::CalculateEntropyDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for(size_t i=0;i<m_densities.size();++i) 
      if (m_densities[i]>0.) 
        ret += m_densities[i] / m_DensitiesId[i] * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_MuStar[i]);
  
    if (m_TemperatureDependentAB) {
      for (size_t i = 0; i < m_densities.size(); ++i) {
        double tPid = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
        for (size_t j = 0; j < m_densities.size(); ++j) {
          ret += -tPid * m_densities[j] * m_VirialdT[j][i];
          ret += m_AttrdT[i][j] * m_densities[i] * m_densities[j];
        }
      }
    }
  
    return ret;
  }

  double ThermalModelVDW::CalculatePressure() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for(size_t i=0;i<m_densities.size();++i) 
      ret += m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
    for(size_t i=0;i<m_densities.size();++i)
      for(size_t j=0;j<m_densities.size();++j)
        ret += -m_Attr[i][j] * m_densities[i] * m_densities[j];
    return ret;
  }


  double ThermalModelVDW::ParticleScalarDensity(int part) {
    if (!m_Calculated) CalculateDensities();

    return m_scaldens[part];
  }

  double ThermalModelVDW::MuShift(int id) const
  {
    if (id >= 0. && id < static_cast<int>(m_Virial.size()))
      return m_MuStar[id] - m_Chem[id];
    else
      return 0.0;
  }

  double ThermalModelVDW::VirialCoefficient(int i, int j) const
  {
    if (i<0 || i >= static_cast<int>(m_Virial.size()) || j < 0 || j >= static_cast<int>(m_Virial.size()))
      return 0.;
    return m_Virial[i][j];
  }

  double ThermalModelVDW::AttractionCoefficient(int i, int j) const
  {
    if (i<0 || i >= static_cast<int>(m_Attr.size()) || j < 0 || j >= static_cast<int>(m_Attr.size()))
      return 0.;
    return m_Attr[i][j];
  }

  double ThermalModelVDW::VirialCoefficientdT(int i, int j) const
  {
    if (i<0 || i >= static_cast<int>(m_VirialdT.size()) || j < 0 || j >= static_cast<int>(m_VirialdT.size()))
      return 0.;
    return m_VirialdT[i][j];
  }

  double ThermalModelVDW::AttractionCoefficientdT(int i, int j) const
  {
    if (i<0 || i >= static_cast<int>(m_AttrdT.size()) || j < 0 || j >= static_cast<int>(m_AttrdT.size()))
      return 0.;
    return m_AttrdT[i][j];
  }

  std::vector<double> ThermalModelVDW::BroydenEquationsVDW::Equations(const std::vector<double>& x)
  {
    int NN = m_THM->Densities().size();
    vector<double> Ps(NN, 0.);
    for (int i = 0; i < NN; ++i) {
      Ps[i] = m_THM->TPS()->Particles()[i].Density(m_THM->Parameters(),
        IdealGasFunctions::Pressure,
        m_THM->UseWidth(),
        m_THM->ChemicalPotential(i) + x[m_THM->m_MapTodMuStar[i]]
      );
    }

    vector<double> ns(NN, 0.);
    for (int i = 0; i < NN; ++i) {
      ns[i] = m_THM->TPS()->Particles()[i].Density(m_THM->Parameters(),
        IdealGasFunctions::ParticleDensity,
        m_THM->UseWidth(),
        m_THM->ChemicalPotential(i) + x[m_THM->m_MapTodMuStar[i]]
      );
    }

    vector<double> np = m_THM->ComputeNp(x, ns);


    vector<double> ret(m_N, 0.);
    for (size_t i = 0; i < ret.size(); ++i) {
      ret[i] = x[i];
      for (int j = 0; j < NN; ++j)
        ret[i] += m_THM->VirialCoefficient(m_THM->m_MapFromdMuStar[i], j) * Ps[j]
        - (m_THM->AttractionCoefficient(m_THM->m_MapFromdMuStar[i], j)
          + m_THM->AttractionCoefficient(j, m_THM->m_MapFromdMuStar[i])) * np[j];
    }
    return ret;
  }

  std::vector<double> ThermalModelVDW::BroydenJacobianVDW::Jacobian(const std::vector<double>& x)
  {
    int NN = m_THM->m_densities.size();
    int NNdmu = m_THM->m_MapFromdMuStar.size();

    bool attrfl = false;
    for (int i = 0; i < NN; ++i) {
      for (int j = 0; j < NN; ++j) {
        if (m_THM->AttractionCoefficient(i, j) != 0.0) {
          attrfl = true;
          break;
        }
      }
      if (attrfl) break;
    }

    MatrixXd densMatrix(NNdmu, NNdmu);
    VectorXd solVector(NNdmu), xVector(NNdmu);

    std::vector<double> ret(NNdmu*NNdmu, 0.);
    {
      vector<double> Ps(NN, 0.);
      for (int i = 0; i<NN; ++i)
        Ps[i] = m_THM->TPS()->Particles()[i].Density(m_THM->Parameters(),
          IdealGasFunctions::Pressure,
          m_THM->UseWidth(),
          m_THM->ChemicalPotential(i) + x[m_THM->m_MapTodMuStar[i]]
        );

      vector<double> ns(NN, 0.);
      for (int i = 0; i<NN; ++i)
        ns[i] = m_THM->TPS()->Particles()[i].Density(m_THM->Parameters(),
          IdealGasFunctions::ParticleDensity,
          m_THM->UseWidth(),
          m_THM->ChemicalPotential(i) + x[m_THM->m_MapTodMuStar[i]]
        );

      vector<double> chi2s(NN, 0.);
      for (int i = 0; i<NN; ++i)
        chi2s[i] = m_THM->TPS()->Particles()[i].chiDimensionfull(2, m_THM->Parameters(),
          m_THM->UseWidth(),
          m_THM->ChemicalPotential(i) + x[m_THM->m_MapTodMuStar[i]]
        );

      for (int i = 0; i < NNdmu; ++i) {
        for (int j = 0; j < NNdmu; ++j) {
          densMatrix(i, j) = 0.;
          if (i == j)
            densMatrix(i, j) += 1.;

          for (size_t m = 0; m < m_THM->m_dMuStarIndices[i].size(); ++m) {
            densMatrix(i, j) += m_THM->m_Virial[m_THM->m_MapFromdMuStar[j]][m_THM->m_dMuStarIndices[i][m]] * ns[m_THM->m_dMuStarIndices[i][m]];
          }
        }
      }

      PartialPivLU<MatrixXd> decomp(densMatrix);


      for (int kp = 0; kp < NNdmu; ++kp) {
        xVector[kp] = 0.;
        for (size_t m = 0; m < m_THM->m_dMuStarIndices[kp].size(); ++m) {
          xVector[kp] += ns[m_THM->m_dMuStarIndices[kp][m]];
        }
      }


      solVector = decomp.solve(xVector);

      vector<double> ntildek(NNdmu, 0.);
      for (int i = 0; i < NNdmu; ++i)
        ntildek[i] = solVector[i];

      vector<double> np(NN, 0.);
      for (int i = 0; i < NN; ++i) {
        np[i] = 0.;
        for (int k = 0; k < NNdmu; ++k) {
          np[i] += m_THM->m_Virial[m_THM->m_MapFromdMuStar[k]][i] * solVector[k];
        }
        np[i] = (1. - np[i]) * ns[i];
      }

      for (int kp = 0; kp < NNdmu; ++kp) {

        if (attrfl) {
          for (int l = 0; l < NNdmu; ++l) {
            xVector[l] = 0.;
            for (size_t m = 0; m < m_THM->m_dMuStarIndices[l].size(); ++m) {
              int ti = m_THM->m_dMuStarIndices[l][m];
              if (m_THM->m_MapTodMuStar[ti] != kp)
                continue;

              double tmps = 1.;
              if (ns[ti] != 0.)
                tmps = np[ti] / ns[ti];
              xVector[l] += chi2s[ti] * pow(xMath::GeVtoifm(), 3) * tmps;
            }
          }

          solVector = decomp.solve(xVector);
          for (int i = 0; i < NNdmu; ++i)
            if (solVector[i] > 1.) solVector[i] = 1.;  // Stabilizer
        }

        vector<double> dnjdmukp(NN, 0.);
        if (attrfl) {
          for (int j = 0; j < NN; ++j) {
            for (int kk = 0; kk < NNdmu; ++kk) {
              dnjdmukp[j] += -m_THM->m_Virial[m_THM->m_MapFromdMuStar[kk]][j] * solVector[kk] * ns[j];
            }

            if (m_THM->m_MapTodMuStar[j] == kp) {
              double tmps = 1.;
              if (ns[j] != 0.)
                tmps = np[j] / ns[j];
              dnjdmukp[j] += tmps * chi2s[j] * pow(xMath::GeVtoifm(), 3);
            }
          }
        }


        for (int k = 0; k < NNdmu; ++k) {
          if (k == kp)
            ret[k*NNdmu + kp] += 1.;
          for (size_t m = 0; m < m_THM->m_dMuStarIndices[kp].size(); ++m) {
            int tj = m_THM->m_dMuStarIndices[kp][m];
            ret[k*NNdmu + kp] += m_THM->m_Virial[m_THM->m_MapFromdMuStar[k]][tj] * ns[tj];
          }

          if (attrfl) {
            for (int j = 0; j < NN; ++j) {
              ret[k*NNdmu + kp] += -(m_THM->m_Attr[m_THM->m_MapFromdMuStar[k]][j] + m_THM->m_Attr[j][m_THM->m_MapFromdMuStar[k]]) * dnjdmukp[j];
            }
          }
        }

      }
    }

    return ret;
  }

  bool ThermalModelVDW::BroydenSolutionCriteriumVDW::IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta) const
  {
    if (xdelta.size() == x.size()) {
      double maxdiff = 0.;
      for (size_t i = 0; i < xdelta.size(); ++i) {
        maxdiff = std::max(maxdiff, fabs(xdelta[i]));
        maxdiff = std::max(maxdiff, fabs(f[i]));
      }
      return (maxdiff < m_MaximumError);
    }
    else {
      return Broyden::BroydenSolutionCriterium::IsSolved(x, f, xdelta);
    }
  }

} // namespace thermalfist