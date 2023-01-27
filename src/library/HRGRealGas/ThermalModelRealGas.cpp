#include "HRGRealGas/ThermalModelRealGas.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <ctime>
#include <iostream>

#include <Eigen/Dense>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


using namespace Eigen;
using namespace std;

namespace thermalfist {

  ThermalModelRealGas::ThermalModelRealGas(ThermalParticleSystem* TPS_, const ThermalModelParameters& params) :
    ThermalModelBase(TPS_, params), m_SearchMultipleSolutions(false)
  {
    m_chi.resize(6);
    for (int i = 0; i < 6; ++i) m_chi[i].resize(3);
    m_chiarb.resize(6);
    m_DensitiesId.resize(m_densities.size());
    m_MuStar.resize(m_densities.size());
    m_Volume = params.V;
    m_TAG = "ThermalModelRealGas";

    m_exvolmodideal = new ExcludedVolumeModelMultiBase(m_densities.size());
    m_mfmodideal = new MeanFieldModelMultiBase(m_densities.size());			
    m_exvolmod = m_exvolmodideal;
    m_mfmod = m_mfmodideal;

    m_Ensemble = GCE;
    m_InteractionModel = RealGas;
  }

  ThermalModelRealGas::~ThermalModelRealGas(void)
  {
    if (m_exvolmodideal != NULL) {
      delete m_exvolmodideal;
      m_exvolmodideal = NULL;
    }
    if (m_mfmodideal != NULL) {
      delete m_mfmodideal;
      m_mfmodideal = NULL;
    }
    if (m_exvolmod != NULL) {
      delete m_exvolmod;
      m_exvolmod = NULL;
    }
    if (m_mfmod != NULL) {
      delete m_mfmod;
      m_mfmod = NULL;
    }
  }

  void ThermalModelRealGas::FillChemicalPotentials() {
    ThermalModelBase::FillChemicalPotentials();
    for (size_t i = 0; i < m_MuStar.size(); ++i)
      m_MuStar[i] = m_Chem[i];
  }

  void ThermalModelRealGas::SetChemicalPotentials(const vector<double>& chem)
  {
    ThermalModelBase::SetChemicalPotentials(chem);
    for (size_t i = 0; i < m_MuStar.size(); ++i)
      m_MuStar[i] = m_Chem[i];
  }

  vector<double> ThermalModelRealGas::SearchSingleSolution(const vector<double>& muStarInit)
  {
    int NN = m_densities.size();

    vector<double> dmuscur(NN, 0.);
    for (int i = 0; i < NN; ++i)
      dmuscur[i] = muStarInit[i] - m_Chem[i];

    BroydenEquationsRealGas eqs(this);
    BroydenJacobianRealGas  jac(this);
    Broyden broydn(&eqs, &jac);
    BroydenSolutionCriteriumRealGas crit(this);

    dmuscur = broydn.Solve(dmuscur, &crit);

    if (broydn.Iterations() == broydn.MaxIterations())
      m_LastBroydenSuccessFlag = false;
    else m_LastBroydenSuccessFlag = true;

    //printf("Iters: %d\n", broydn.Iterations());

    m_MaxDiff = broydn.MaxDifference();

    vector<double> ret(NN);
    for (int i = 0; i < NN; ++i)
      ret[i] = m_Chem[i] + dmuscur[i];

    return ret;
  }

  vector<double> ThermalModelRealGas::SearchMultipleSolutions(int iters) {
    vector<double> csol(m_densities.size(), 0.);
    double Psol = 0.;
    bool solved = false;
    double muBmin = m_Parameters.muB - 0.5 * xMath::mnucleon();
    double muBmax = m_Parameters.muB + 0.5 * xMath::mnucleon();
    double dmu = (muBmax - muBmin) / iters;
    vector<double> curmust(m_densities.size(), 0.);
    double maxdif = 0.;
    for (int isol = 0; isol < iters; ++isol) {
      double tmu = muBmin + (0.5 + isol) * dmu;
      for (size_t j = 0; j < curmust.size(); ++j) {
        curmust[j] = m_Chem[j] + (tmu - m_Parameters.muB) * m_Chem[j] / m_Parameters.muB;
        if (m_TPS->Particles()[j].Statistics() == -1 && curmust[j] > m_TPS->Particles()[j].Mass())
          curmust[j] = 0.98 * m_TPS->Particles()[j].Mass();
      }

      vector<double> sol = SearchSingleSolution(curmust);

      bool fl = true;
      for (size_t i = 0; i < sol.size(); ++i)
        if (sol[i] != sol[i]) fl = false;
      fl &= m_LastBroydenSuccessFlag;
      if (!fl) continue;

      for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
        m_DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, sol[i]);

      m_densities = m_exvolmod->nsol(m_DensitiesId);
      m_exvolmod->SetDensities(m_densities);
      m_mfmod->SetDensities(m_densities);

      double tP = 0.;
      for (size_t i = 0; i < m_densities.size(); ++i) {
        double tfsum = 0.;
        for (size_t j = 0; j < m_densities.size(); ++j) {
          tfsum += m_densities[j] * m_exvolmod->df(i, j);
        }
        tfsum = m_exvolmod->f(i) - tfsum;
        tP += tfsum * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, sol[i]);
      }

      tP += -m_mfmod->v();
      for (size_t i = 0; i < m_densities.size(); ++i)
        tP += m_mfmod->dv(i) * m_densities[i];

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

  void ThermalModelRealGas::SolveEquations() {
    if (!m_SearchMultipleSolutions) {

      vector<double> muStarInit = m_MuStar;

      for (size_t i = 0; i < muStarInit.size(); ++i) {
        if (m_TPS->Particles()[i].Statistics() == -1 && muStarInit[i] > m_TPS->Particles()[i].Mass())
          muStarInit[i] = 0.98 * m_TPS->Particles()[i].Mass();
      }
      m_MuStar = SearchSingleSolution(muStarInit);
    }
    else {
      m_MuStar = SearchMultipleSolutions(100);
    }
  }

  void ThermalModelRealGas::CalculatePrimordialDensities() {
    m_FluctuationsCalculated = false;

    int NN = m_densities.size();

    SolveEquations();

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i)
      m_DensitiesId[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_MuStar[i]);

    m_densities = m_exvolmod->nsol(m_DensitiesId);
    m_exvolmod->SetDensities(m_densities);
    m_mfmod->SetDensities(m_densities);

    // TODO: Scalar density properly
    m_scaldens = m_densities;

    m_Calculated = true;

    ValidateCalculation();
  }

  vector<double> ThermalModelRealGas::CalculateChargeFluctuations(const vector<double>& chgs, int order) {
    vector<double> ret(order + 1, 0.);

    // chi1
    for (size_t i = 0; i < m_densities.size(); ++i)
      ret[0] += chgs[i] * m_densities[i];

    ret[0] /= pow(m_Parameters.T * xMath::GeVtoifm(), 3);

    if (order < 2) return ret;
    // Preparing matrix for system of linear equations
    int NN = m_densities.size();
    int Nevcomp = m_exvolmod->ComponentsNumber();
    const vector<int>& evinds = m_exvolmod->ComponentIndices();
    const vector<int>& evindsfrom = m_exvolmod->ComponentIndicesFrom();

    int Nmfcomp = m_mfmod->ComponentsNumber();
    const vector<int>& mfinds = m_mfmod->ComponentIndices();
    const vector<int>& mfindsfrom = m_mfmod->ComponentIndicesFrom();

    MatrixXd densMatrix(2 * NN, 2 * NN);
    VectorXd solVector(2 * NN), xVector(2 * NN);

    vector<double> chi2id(m_densities.size()), Ps(m_densities.size());
    for (int i = 0; i < NN; ++i) {
      //chi2id[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_MuStar[i]);
      Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
      chi2id[i] = m_TPS->Particles()[i].chiDimensionfull(2, m_Parameters, m_UseWidth, m_MuStar[i]) * pow(xMath::GeVtoifm(), 3);
    }

    vector<double> evc_chi2id(Nevcomp, 0.), evc_Ps(Nevcomp, 0.), evc_ns(Nevcomp, 0.);
    for (int i = 0; i < NN; ++i) {
      evc_Ps[evinds[i]]     += Ps[i];
      evc_ns[evinds[i]]     += m_DensitiesId[i];
      evc_chi2id[evinds[i]] += chi2id[i];
    }

    vector<vector<double>> pkd2fkij(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      //int indi = evinds[i];
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        //int indj = evinds[j];
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  pkd2fkij[indi][indj] += Ps[k] * m_exvolmod->d2f(k, i, j);
        //}
        for (int k = 0; k < Nevcomp; ++k) {
          pkd2fkij[indi][indj] += evc_Ps[k] *  m_exvolmod->d2f(evindsfrom[k], i, j);
        }
      }
    }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = -m_exvolmod->df(i, j) * m_DensitiesId[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j)
        densMatrix(i, NN + j) = 0.;

    for (int i = 0; i < NN; ++i) {
      densMatrix(i, NN + i) = -m_exvolmod->f(i) * chi2id[i];
    }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(NN + i, j) = m_mfmod->d2v(i, j);
        //for (int k = 0; k < NN; ++k)
        //  densMatrix(NN + i, j) += -m_exvolmod->d2f(k, i, j) * Ps[k];
        densMatrix(NN + i, j) += -pkd2fkij[evinds[i]][evinds[j]];
      }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(NN + i, NN + j) = -m_exvolmod->df(j, i) * m_DensitiesId[j];
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

    vector<double> evc_dn(Nevcomp, 0.), evc_dmus(Nevcomp, 0.), evc_nsdmus(Nevcomp, 0.);
    for (int i = 0; i < NN; ++i) {
      evc_dn[evinds[i]] += dni[i];
      evc_dmus[evinds[i]] += dmus[i];
      evc_nsdmus[evinds[i]] += m_DensitiesId[i] * dmus[i];
    }

    vector<double> mfc_dn(Nmfcomp, 0.);
    for (int i = 0; i < NN; ++i) {
      mfc_dn[mfinds[i]] += dni[i];
    }

    // chi3
    vector<double> d2ni(NN, 0.), d2mus(NN, 0.);

    vector<double> chi3id(m_densities.size());
    for (int i = 0; i < NN; ++i)
      chi3id[i] = m_TPS->Particles()[i].chiDimensionfull(3, m_Parameters, m_UseWidth, m_MuStar[i]) * pow(xMath::GeVtoifm(), 3);


    vector<vector<double>> d2fijkdnk(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  d2fijkdnk[indi][indj] += dni[k] * m_exvolmod->d2f(i, j, k);
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          d2fijkdnk[indi][indj] += evc_dn[indk] * m_exvolmod->d2f(i, j, k);
        }
      }
    }

    vector<double> dfikdnk(Nevcomp, 0.);
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      //for (int k = 0; k < NN; ++k) {
      //  dfikdnk[indi] += dni[k] * m_exvolmod->df(i, k);
      //}
      for (int indk = 0; indk < Nevcomp; ++indk) {
        int k = evindsfrom[indk];
        dfikdnk[indi] += evc_dn[indk] * m_exvolmod->df(i, k);
      }
    }

    vector<vector<double>> d2fkijnskmusk(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  d2fkijnskmusk[indi][indj] += m_exvolmod->d2f(k, i, j) * m_DensitiesId[k] * dmus[k];
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          d2fkijnskmusk[indi][indj] += m_exvolmod->d2f(k, i, j) * evc_nsdmus[indk];
        }
      }
    }

    vector<vector<double>> pkd3fkijmdnm(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  for (int m = 0; m < NN; ++m) {
        //    pkd3fkijmdnm[indi][indj] += Ps[k] * m_exvolmod->d3f(k, i, j, m) * dni[m];
        //  }
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk]; 
          for (int indm = 0; indm < Nevcomp; ++indm) {
            int m = evindsfrom[indm];
            pkd3fkijmdnm[indi][indj] += evc_Ps[indk] * m_exvolmod->d3f(k, i, j, m) * evc_dn[indm];
          }
        }
      }
    }

    vector<vector<double>> d3vijkdnk(Nmfcomp, vector<double>(Nmfcomp, 0.));
    for (int indi = 0; indi < Nmfcomp; ++indi) {
      int i = mfindsfrom[indi];
      for (int indj = 0; indj < Nmfcomp; ++indj) {
        int j = mfindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  d3vijkdnk[indi][indj] += dni[k] * m_mfmod->d3v(i, j, k);
        //}
        for (int indk = 0; indk < Nmfcomp; ++indk) {
          int k = mfindsfrom[indk];
          d3vijkdnk[indi][indj] += mfc_dn[indk] * m_mfmod->d3v(i, j, k);
        }
      }
    }

    vector< vector<double> > daij11, daij12, daij21, daij22;

    daij11.resize(NN);
    daij12.resize(NN);
    daij21.resize(NN);
    daij22.resize(NN);
    for (int i = 0; i < NN; ++i) {
      //cout << "chi3 iter: " << i << "\n";
      daij11[i].resize(NN);
      daij12[i].resize(NN);
      daij21[i].resize(NN);
      daij22[i].resize(NN);
      for (int j = 0; j < NN; ++j) {
        daij11[i][j] = 0.;
        //for (int k = 0; k < NN; ++k)
        //  daij11[i][j] += -m_exvolmod->d2f(i, j, k) * dni[k] * m_DensitiesId[i];
        daij11[i][j] += -d2fijkdnk[evinds[i]][evinds[j]] * m_DensitiesId[i];
        daij11[i][j] += -m_exvolmod->df(i, j) * chi2id[i] * dmus[i];

        daij12[i][j] = 0.;
        if (i == j) {
          //for (int k = 0; k < NN; ++k)
          //  daij12[i][j] += -m_exvolmod->df(i, k) * chi2id[i] * dni[k];
          daij12[i][j] += -dfikdnk[evinds[i]] * chi2id[i];
          daij12[i][j] += -m_exvolmod->f(i) * chi3id[i] * dmus[i];
        }


        daij21[i][j] = 0.;
        daij21[i][j] += d3vijkdnk[mfinds[i]][mfinds[j]];
        //for (int k = 0; k < NN; ++k) {
        //  daij21[i][j] += m_mfmod->d3v(i, j, k) * dni[k];
        //  //daij21[i][j] += -m_exvolmod->d2f(k, i, j) * m_DensitiesId[k] * dmus[k];
        //  //for (int m = 0; m < NN; ++m)
        //  //  daij21[i][j] += -Ps[k] * m_exvolmod->d3f(k, i, j, m) * dni[m];
        //}
        daij21[i][j] += -d2fkijnskmusk[evinds[i]][evinds[j]];
        daij21[i][j] += -pkd3fkijmdnm[evinds[i]][evinds[j]];

        daij22[i][j] = 0.;
        daij22[i][j] += -m_exvolmod->df(j, i) * chi2id[j] * dmus[j];
        //for (int k = 0; k < NN; ++k)
        //  daij22[i][j] += -m_exvolmod->d2f(j, i, k) * m_DensitiesId[j] * dni[k];
        daij22[i][j] += -m_DensitiesId[j] * d2fijkdnk[evinds[j]][evinds[i]];
      }
    }


    for (int i = 0; i < NN; ++i) {
      xVector[i] = 0.;

      for (int j = 0; j < NN; ++j)
        xVector[i] += -daij11[i][j] * dni[j];

      for (int j = 0; j < NN; ++j)
        xVector[i] += -daij12[i][j] * dmus[j];
    }
    for (int i = 0; i < NN; ++i) {
      xVector[NN + i] = 0.;

      for (int j = 0; j < NN; ++j)
        xVector[NN + i] += -daij21[i][j] * dni[j];

      for (int j = 0; j < NN; ++j)
        xVector[NN + i] += -daij22[i][j] * dmus[j];
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

    vector<double> evc_d2n(Nevcomp, 0.), evc_d2mus(Nevcomp, 0.), evc_nsd2mus(Nevcomp, 0.), evc_chi2iddmus2(Nevcomp, 0.);
    for (int i = 0; i < NN; ++i) {
      evc_d2n[evinds[i]] += d2ni[i];
      evc_d2mus[evinds[i]] += d2mus[i];
      evc_nsd2mus[evinds[i]] += m_DensitiesId[i] * d2mus[i];
      evc_chi2iddmus2[evinds[i]] += chi2id[i] * dmus[i] * dmus[i];
    }

    vector<double> mfc_d2n(Nmfcomp, 0.);
    for (int i = 0; i < NN; ++i) {
      mfc_d2n[mfinds[i]] += d2ni[i];
    }

    // chi4
    vector<double> d3ni(NN, 0.), d3mus(NN, 0.);

    vector<double> chi4id(m_densities.size());
    for (int i = 0; i < NN; ++i)
      chi4id[i] = m_TPS->Particles()[i].chiDimensionfull(4, m_Parameters, m_UseWidth, m_MuStar[i]) * pow(xMath::GeVtoifm(), 3);

    vector<vector<double>> d2fijkd2nk(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  d2fijkd2nk[indi][indj] += d2ni[k] * m_exvolmod->d2f(i, j, k);
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          d2fijkd2nk[indi][indj] += evc_d2n[indk] * m_exvolmod->d2f(i, j, k);
        }
      }
    }

    vector<vector<double>> d3fijkmdnkdnm(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  for (int m = 0; m < NN; ++m) {
        //    d3fijkmdnkdnm[indi][indj] += m_exvolmod->d3f(i, j, k, m) * dni[k] * dni[m];
        //  }
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          for (int indm = 0; indm < Nevcomp; ++indm) {
            int m = evindsfrom[indm];
            d3fijkmdnkdnm[indi][indj] += m_exvolmod->d3f(i, j, k, m) * evc_dn[indk] * evc_dn[indm];
          }
        }
      }
    }

    vector<double> dfikd2nk(Nevcomp, 0.);
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      //for (int k = 0; k < NN; ++k) {
      //  dfikd2nk[indi] += d2ni[k] * m_exvolmod->df(i, k);
      //}
      for (int indk = 0; indk < Nevcomp; ++indk) {
        int k = evindsfrom[indk];
        dfikd2nk[indi] += evc_d2n[indk] * m_exvolmod->df(i, k);
      }
    }

    vector<double> d2fikmdnkdnm(Nevcomp, 0.);
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      //for (int k = 0; k < NN; ++k) {
      //  for (int m = 0; m < NN; ++m) {
      //    d2fikmdnkdnm[indi] += m_exvolmod->d2f(i, k, m) * dni[k] * dni[m];
      //  }
      //}
      for (int indk = 0; indk < Nevcomp; ++indk) {
        int k = evindsfrom[indk];
        for (int indm = 0; indm < Nevcomp; ++indm) {
          int m = evindsfrom[indm];
          d2fikmdnkdnm[indi] += m_exvolmod->d2f(i, k, m) * evc_dn[indk] * evc_dn[indm];
        }
      }
    }

    vector<vector<double>> d2fkijnskd2musk(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  d2fkijnskd2musk[indi][indj] += m_exvolmod->d2f(k, i, j) * m_DensitiesId[k] * d2mus[k];
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          d2fkijnskd2musk[indi][indj] += m_exvolmod->d2f(k, i, j) * evc_nsd2mus[indk];
        }
      }
    }

    vector<vector<double>> d2fkijc2kdmuskdmusk(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  d2fkijc2kdmuskdmusk[indi][indj] += m_exvolmod->d2f(k, i, j) * chi2id[k] * dmus[k] * dmus[k];
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          d2fkijc2kdmuskdmusk[indi][indj] += m_exvolmod->d2f(k, i, j) * evc_chi2iddmus2[indk];
        }
      }
    }

    vector<vector<double>> nskd3fkijmdmuskdnm(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  for (int m = 0; m < NN; ++m) {
        //    nskd3fkijmdmuskdnm[indi][indj] += m_DensitiesId[k] * m_exvolmod->d3f(k, i, j, m) * dmus[k] * dni[m];
        //  }
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          for (int indm = 0; indm < Nevcomp; ++indm) {
            int m = evindsfrom[indm];
            nskd3fkijmdmuskdnm[indi][indj] += evc_nsdmus[indk] * m_exvolmod->d3f(k, i, j, m) * evc_dn[indm];
          }
        }
      }
    }

    vector<vector<double>> pkd3fkijmd2nm(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  for (int m = 0; m < NN; ++m) {
        //    pkd3fkijmd2nm[indi][indj] += Ps[k] * m_exvolmod->d3f(k, i, j, m) * d2ni[m];
        //  }
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          for (int indm = 0; indm < Nevcomp; ++indm) {
            int m = evindsfrom[indm];
            pkd3fkijmd2nm[indi][indj] += evc_Ps[indk] * m_exvolmod->d3f(k, i, j, m) * evc_d2n[indm];
          }
        }
      }
    }

    vector<vector<double>> pkd4fkijmldnmdnl(Nevcomp, vector<double>(Nevcomp, 0.));
    for (int indi = 0; indi < Nevcomp; ++indi) {
      int i = evindsfrom[indi];
      for (int indj = 0; indj < Nevcomp; ++indj) {
        int j = evindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  for (int m = 0; m < NN; ++m) {
        //    for (int l = 0; m < NN; ++m) {
        //      pkd4fkijmldnmdnl[indi][indj] += Ps[k] * m_exvolmod->d4f(k, i, j, m, l) * dni[m] * dni[l];
        //    }
        //  }
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          for (int indm = 0; indm < Nevcomp; ++indm) {
            int m = evindsfrom[indm];
            for (int indl = 0; indl < Nevcomp; ++indl) {
              int l = evindsfrom[indl];
              pkd4fkijmldnmdnl[indi][indj] += evc_Ps[indk] * m_exvolmod->d4f(k, i, j, m, l) * evc_dn[indm] * evc_dn[indl];
            }
          }
        }
      }
    }

    vector<vector<double>> d3vijkd2nk(Nmfcomp, vector<double>(Nmfcomp, 0.));
    for (int indi = 0; indi < Nmfcomp; ++indi) {
      int i = mfindsfrom[indi];
      for (int indj = 0; indj < Nmfcomp; ++indj) {
        int j = mfindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  d3vijkd2nk[indi][indj] += d2ni[k] * m_mfmod->d3v(i, j, k);
        //}
        for (int indk = 0; indk < Nmfcomp; ++indk) {
          int k = mfindsfrom[indk];
          d3vijkd2nk[indi][indj] += mfc_d2n[indk] * m_mfmod->d3v(i, j, k);
        }
      }
    }

    vector<vector<double>> d4vijkmdnkdnm(Nmfcomp, vector<double>(Nmfcomp, 0.));
    for (int indi = 0; indi < Nmfcomp; ++indi) {
      int i = mfindsfrom[indi];
      for (int indj = 0; indj < Nmfcomp; ++indj) {
        int j = mfindsfrom[indj];
        //for (int k = 0; k < NN; ++k) {
        //  for (int m = 0; m < NN; ++m) {
        //    d4vijkmdnkdnm[indi][indj] += dni[k] * dni[m] * m_mfmod->d4v(i, j, k, m);
        //  }
        //}
        for (int indk = 0; indk < Nmfcomp; ++indk) {
          int k = mfindsfrom[indk];
          for (int indm = 0; indm < Nmfcomp; ++indm) {
            int m = mfindsfrom[indm];
            d4vijkmdnkdnm[indi][indj] += mfc_dn[indk] * mfc_dn[indm] * m_mfmod->d4v(i, j, k, m);
          }
        }
      }
    }

    vector< vector<double> > d2aij11, d2aij12, d2aij21, d2aij22;

    d2aij11.resize(NN);
    d2aij12.resize(NN);
    d2aij21.resize(NN);
    d2aij22.resize(NN);
    for (int i = 0; i < NN; ++i) {
      //cout << "chi4 iter: " << i << "\n";
      d2aij11[i].resize(NN);
      d2aij12[i].resize(NN);
      d2aij21[i].resize(NN);
      d2aij22[i].resize(NN);
      for (int j = 0; j < NN; ++j) {
        d2aij11[i][j] = 0.;
        d2aij11[i][j] += -m_exvolmod->df(i, j) * chi3id[i] * dmus[i] * dmus[i];
        d2aij11[i][j] += -m_exvolmod->df(i, j) * chi2id[i] * d2mus[i];

        d2aij11[i][j] += -2. * d2fijkdnk[evinds[i]][evinds[j]] * chi2id[i] * dmus[i];
        d2aij11[i][j] += -d2fijkd2nk[evinds[i]][evinds[j]] * m_DensitiesId[i];
        d2aij11[i][j] += -d3fijkmdnkdnm[evinds[i]][evinds[j]] * m_DensitiesId[i];

        //for (int k = 0; k < NN; ++k) {
        //  //d2aij11[i][j] += -m_exvolmod->d2f(i, j, k) * chi2id[i] * dmus[i] * dni[k];
        //  //d2aij11[i][j] += -m_exvolmod->d2f(i, j, k) * m_DensitiesId[i] * d2ni[k];
        //  //d2aij11[i][j] += -m_exvolmod->d2f(i, j, k) * chi2id[i] * dmus[i] * dni[k];
        //  //for (int m = 0; m < NN; ++m) {
        //  //  d2aij11[i][j] += -m_exvolmod->d3f(i, j, k, m) * m_DensitiesId[i] * dni[k] * dni[m];
        //  //}
        //}

        d2aij12[i][j] = 0.;
        if (i == j) {
          d2aij12[i][j] += -m_exvolmod->f(i) * chi3id[i] * d2mus[i];
          d2aij12[i][j] += -m_exvolmod->f(i) * chi4id[i] * dmus[i] * dmus[i];

          d2aij12[i][j] += -2. * dfikdnk[evinds[i]] * chi3id[i] * dmus[i];
          d2aij12[i][j] += -dfikd2nk[evinds[i]] * chi2id[i];
          d2aij12[i][j] += -d2fikmdnkdnm[evinds[i]] * chi2id[i];

          //for (int k = 0; k < NN; ++k) {
          //  d2aij12[i][j] += -2. * m_exvolmod->df(i, k) * chi3id[i] * dmus[i] * dni[k];
          //  d2aij12[i][j] += -m_exvolmod->df(i, k) * chi2id[i] * d2ni[k];

          //  //for (int m = 0; m < NN; ++m) {
          //  //  d2aij12[i][j] += -m_exvolmod->d2f(i, k, m) * chi2id[i] * dni[k] * dni[m];
          //  //}
          //}
        }

        d2aij21[i][j] = 0.;

        d2aij21[i][j] += -d2fkijnskd2musk[evinds[i]][evinds[j]];
        d2aij21[i][j] += -d2fkijc2kdmuskdmusk[evinds[i]][evinds[j]];
        d2aij21[i][j] += -2. * nskd3fkijmdmuskdnm[evinds[i]][evinds[j]];
        d2aij21[i][j] += -pkd3fkijmd2nm[evinds[i]][evinds[j]];
        d2aij21[i][j] += -pkd4fkijmldnmdnl[evinds[i]][evinds[j]];

        d2aij21[i][j] += d3vijkd2nk[mfinds[i]][mfinds[j]];
        d2aij21[i][j] += d4vijkmdnkdnm[mfinds[i]][mfinds[j]];

        //for (int k = 0; k < NN; ++k) {
        //  ////daij21[i][j] += m_mfmod->d3v(i, j, k) * dni[k];
        //  d2aij21[i][j] += m_mfmod->d3v(i, j, k) * d2ni[k];
        //  for (int m = 0; m < NN; ++m)
        //    d2aij21[i][j] += m_mfmod->d4v(i, j, k, m) * dni[k] * dni[m];
        //  ////daij21[i][j] += -m_exvolmod->d2f(k, i, j) * m_DensitiesId[k] * dmus[k];
        //  //d2aij21[i][j] += -m_exvolmod->d2f(k, i, j) * m_DensitiesId[k] * d2mus[k];
        //  //d2aij21[i][j] += -m_exvolmod->d2f(k, i, j) * chi2id[k] * dmus[k] * dmus[k];
        //  //for (int m = 0; m < NN; ++m)
        //  //  d2aij21[i][j] += -m_exvolmod->d3f(k, i, j, m) * m_DensitiesId[k] * dmus[k] * dni[m];

        //  //for (int m = 0; m < NN; ++m) {
        //  //  ////daij21[i][j] += -Ps[k] * m_exvolmod->d3f(k, i, j, m) * dni[m];
        //  //  //d2aij21[i][j] += -m_DensitiesId[k] * m_exvolmod->d3f(k, i, j, m) * dni[m] * dmus[k];
        //  //  //d2aij21[i][j] += -Ps[k] * m_exvolmod->d3f(k, i, j, m) * d2ni[m];
        //  //  //for (int l = 0; l < NN; ++l)
        //  //  //  d2aij21[i][j] += -Ps[k] * m_exvolmod->d4f(k, i, j, m, l) * dni[m] * dni[l];
        //  //}
        //}

        d2aij22[i][j] = 0.;
        //daij22[i][j] += -m_exvolmod->df(j, i) * chi2id[j] * dmus[j];
        d2aij22[i][j] += -m_exvolmod->df(j, i) * chi3id[j] * dmus[j] * dmus[j];
        d2aij22[i][j] += -m_exvolmod->df(j, i) * chi2id[j] * d2mus[j];

        d2aij22[i][j] += -2. * d2fijkdnk[evinds[j]][evinds[i]] * chi2id[j] * dmus[j];

        //for (int k = 0; k < NN; ++k)
        //  d2aij22[i][j] += -m_exvolmod->d2f(j, i, k) * chi2id[j] * dmus[j] * dni[k];

        d2aij22[i][j] += -d2fijkd2nk[evinds[j]][evinds[i]] * m_DensitiesId[j];
        d2aij22[i][j] += -d3fijkmdnkdnm[evinds[j]][evinds[i]] * m_DensitiesId[j];

        //for (int k = 0; k < NN; ++k) {
        //  ////daij22[i][j] += -m_exvolmod->d2f(j, i, k) * m_DensitiesId[j] * dni[k];
        //  //d2aij22[i][j] += -m_exvolmod->d2f(j, i, k) * chi2id[j] * dni[k] * dmus[j];
        //  //d2aij22[i][j] += -m_exvolmod->d2f(j, i, k) * m_DensitiesId[j] * d2ni[k];
        //  //for (int m = 0; m < NN; ++m) {
        //  //  d2aij22[i][j] += -m_exvolmod->d3f(j, i, k, m) * m_DensitiesId[j] * dni[k] * dni[m];
        //  //}
        //}
      }
    }


    for (int i = 0; i < NN; ++i) {
      xVector[i] = 0.;

      for (int j = 0; j < NN; ++j)
        xVector[i] += -2. * daij11[i][j] * d2ni[j] - d2aij11[i][j] * dni[j];

      for (int j = 0; j < NN; ++j)
        xVector[i] += -2. * daij12[i][j] * d2mus[j] - d2aij12[i][j] * dmus[j];
    }
    for (int i = 0; i < NN; ++i) {
      xVector[NN + i] = 0.;

      for (int j = 0; j < NN; ++j)
        xVector[NN + i] += -2. * daij21[i][j] * d2ni[j] - d2aij21[i][j] * dni[j];

      for (int j = 0; j < NN; ++j)
        xVector[NN + i] += -2. * daij22[i][j] * d2mus[j] - d2aij22[i][j] * dmus[j];
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

  vector<double> ThermalModelRealGas::CalculateChargeFluctuationsOld(const vector<double>& chgs, int order) {
    vector<double> ret(order + 1, 0.);

    // chi1
    for (size_t i = 0; i < m_densities.size(); ++i)
      ret[0] += chgs[i] * m_densities[i];

    ret[0] /= pow(m_Parameters.T * xMath::GeVtoifm(), 3);

    if (order < 2) return ret;
    // Preparing matrix for system of linear equations
    int NN = m_densities.size();

    MatrixXd densMatrix(2 * NN, 2 * NN);
    VectorXd solVector(2 * NN), xVector(2 * NN);

    vector<double> chi2id(m_densities.size()), Ps(m_densities.size());
    for (int i = 0; i < NN; ++i) {
      //chi2id[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_MuStar[i]);
      Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
      chi2id[i] = m_TPS->Particles()[i].chiDimensionfull(2, m_Parameters, m_UseWidth, m_MuStar[i]) * pow(xMath::GeVtoifm(), 3);
    }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = -m_exvolmod->df(i, j) * m_DensitiesId[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j)
        densMatrix(i, NN + j) = 0.;

    for (int i = 0; i < NN; ++i) {
      densMatrix(i, NN + i) = -m_exvolmod->f(i) * chi2id[i];
    }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(NN + i, j) = m_mfmod->d2v(i, j);
        for (int k = 0; k < NN; ++k) {
          densMatrix(NN + i, j) += -m_exvolmod->d2f(k, i, j) * Ps[k];
        }
      }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(NN + i, NN + j) = -m_exvolmod->df(j, i) * m_DensitiesId[j];
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
      chi3id[i] = m_TPS->Particles()[i].chiDimensionfull(3, m_Parameters, m_UseWidth, m_MuStar[i]) * pow(xMath::GeVtoifm(), 3);

    vector< vector<double> > daij11, daij12, daij21, daij22;

    daij11.resize(NN);
    daij12.resize(NN);
    daij21.resize(NN);
    daij22.resize(NN);
    for (int i = 0; i < NN; ++i) {
      cout << "chi3 iter: " << i << "\n";
      daij11[i].resize(NN);
      daij12[i].resize(NN);
      daij21[i].resize(NN);
      daij22[i].resize(NN);
      for (int j = 0; j < NN; ++j) {
        daij11[i][j] = 0.;
        for (int k = 0; k < NN; ++k)
          daij11[i][j] += -m_exvolmod->d2f(i, j, k) * dni[k] * m_DensitiesId[i];
        daij11[i][j] += -m_exvolmod->df(i, j) * chi2id[i] * dmus[i];

        daij12[i][j] = 0.;
        if (i == j) {
          for (int k = 0; k < NN; ++k)
            daij12[i][j] += -m_exvolmod->df(i, k) * chi2id[i] * dni[k];
          daij12[i][j] += -m_exvolmod->f(i) * chi3id[i] * dmus[i];
        }


        daij21[i][j] = 0.;
        for (int k = 0; k < NN; ++k) {
          daij21[i][j] += m_mfmod->d3v(i, j, k) * dni[k];
          daij21[i][j] += -m_exvolmod->d2f(k, i, j) * m_DensitiesId[k] * dmus[k];
          for (int m = 0; m < NN; ++m)
            daij21[i][j] += -Ps[k] * m_exvolmod->d3f(k, i, j, m) * dni[m];
        }

        daij22[i][j] = 0.;
        daij22[i][j] += -m_exvolmod->df(j, i) * chi2id[j] * dmus[j];
        for (int k = 0; k < NN; ++k)
          daij22[i][j] += -m_exvolmod->d2f(j, i, k) * m_DensitiesId[j] * dni[k];
      }
    }


    for (int i = 0; i < NN; ++i) {
      xVector[i] = 0.;

      for (int j = 0; j < NN; ++j)
        xVector[i] += -daij11[i][j] * dni[j];

      for (int j = 0; j < NN; ++j)
        xVector[i] += -daij12[i][j] * dmus[j];
    }
    for (int i = 0; i < NN; ++i) {
      xVector[NN + i] = 0.;

      for (int j = 0; j < NN; ++j)
        xVector[NN + i] += -daij21[i][j] * dni[j];

      for (int j = 0; j < NN; ++j)
        xVector[NN + i] += -daij22[i][j] * dmus[j];
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
      chi4id[i] = m_TPS->Particles()[i].chiDimensionfull(4, m_Parameters, m_UseWidth, m_MuStar[i]) * pow(xMath::GeVtoifm(), 3);

    vector< vector<double> > d2aij11, d2aij12, d2aij21, d2aij22;

    d2aij11.resize(NN);
    d2aij12.resize(NN);
    d2aij21.resize(NN);
    d2aij22.resize(NN);
    for (int i = 0; i < NN; ++i) {
      cout << "chi4 iter: " << i << "\n";
      d2aij11[i].resize(NN);
      d2aij12[i].resize(NN);
      d2aij21[i].resize(NN);
      d2aij22[i].resize(NN);
      for (int j = 0; j < NN; ++j) {
        d2aij11[i][j] = 0.;
        d2aij11[i][j] += -m_exvolmod->df(i, j) * chi3id[i] * dmus[i] * dmus[i];
        d2aij11[i][j] += -m_exvolmod->df(i, j) * chi2id[i] * d2mus[i];
        for (int k = 0; k < NN; ++k) {
          d2aij11[i][j] += -m_exvolmod->d2f(i, j, k) * chi2id[i] * dmus[i] * dni[k];
          d2aij11[i][j] += -m_exvolmod->d2f(i, j, k) * m_DensitiesId[i] * d2ni[k];
          d2aij11[i][j] += -m_exvolmod->d2f(i, j, k) * chi2id[i] * dmus[i] * dni[k];
          for (int m = 0; m < NN; ++m) {
            d2aij11[i][j] += -m_exvolmod->d3f(i, j, k, m) * m_DensitiesId[i] * dni[k] * dni[m];
          }
        }

        d2aij12[i][j] = 0.;
        if (i == j) {
          d2aij12[i][j] += -m_exvolmod->f(i) * chi3id[i] * d2mus[i];
          d2aij12[i][j] += -m_exvolmod->f(i) * chi4id[i] * dmus[i] * dmus[i];
          for (int k = 0; k < NN; ++k) {
            d2aij12[i][j] += -2. * m_exvolmod->df(i, k) * chi3id[i] * dmus[i] * dni[k];
            d2aij12[i][j] += -m_exvolmod->df(i, k) * chi2id[i] * d2ni[k];
            for (int m = 0; m < NN; ++m) {
              d2aij12[i][j] += -m_exvolmod->d2f(i, k, m) * chi2id[i] * dni[k] * dni[m];
            }
          }
        }

        d2aij21[i][j] = 0.;
        for (int k = 0; k < NN; ++k) {
          //daij21[i][j] += m_mfmod->d3v(i, j, k) * dni[k];
          d2aij21[i][j] += m_mfmod->d3v(i, j, k) * d2ni[k];
          for (int m = 0; m < NN; ++m)
            d2aij21[i][j] += m_mfmod->d4v(i, j, k, m) * dni[k] * dni[m];
          //daij21[i][j] += -m_exvolmod->d2f(k, i, j) * m_DensitiesId[k] * dmus[k];
          d2aij21[i][j] += -m_exvolmod->d2f(k, i, j) * m_DensitiesId[k] * d2mus[k];
          d2aij21[i][j] += -m_exvolmod->d2f(k, i, j) * chi2id[k] * dmus[k] * dmus[k];
          for (int m = 0; m < NN; ++m)
            d2aij21[i][j] += -m_exvolmod->d3f(k, i, j, m) * m_DensitiesId[k] * dmus[k] * dni[m];

          for (int m = 0; m < NN; ++m) {
            //daij21[i][j] += -Ps[k] * m_exvolmod->d3f(k, i, j, m) * dni[m];
            d2aij21[i][j] += -m_DensitiesId[k] * m_exvolmod->d3f(k, i, j, m) * dni[m] * dmus[k];
            d2aij21[i][j] += -Ps[k] * m_exvolmod->d3f(k, i, j, m) * d2ni[m];
            for (int l = 0; l < NN; ++l)
              d2aij21[i][j] += -Ps[k] * m_exvolmod->d4f(k, i, j, m, l) * dni[m] * dni[l];
          }
        }

        d2aij22[i][j] = 0.;
        //daij22[i][j] += -m_exvolmod->df(j, i) * chi2id[j] * dmus[j];
        d2aij22[i][j] += -m_exvolmod->df(j, i) * chi3id[j] * dmus[j] * dmus[j];
        d2aij22[i][j] += -m_exvolmod->df(j, i) * chi2id[j] * d2mus[j];
        for (int k = 0; k < NN; ++k)
          d2aij22[i][j] += -m_exvolmod->d2f(j, i, k) * chi2id[j] * dmus[j] * dni[k];
        for (int k = 0; k < NN; ++k) {
          //daij22[i][j] += -m_exvolmod->d2f(j, i, k) * m_DensitiesId[j] * dni[k];
          d2aij22[i][j] += -m_exvolmod->d2f(j, i, k) * chi2id[j] * dni[k] * dmus[j];
          d2aij22[i][j] += -m_exvolmod->d2f(j, i, k) * m_DensitiesId[j] * d2ni[k];
          for (int m = 0; m < NN; ++m) {
            d2aij22[i][j] += -m_exvolmod->d3f(j, i, k, m) * m_DensitiesId[j] * dni[k] * dni[m];
          }
        }
      }
    }


    for (int i = 0; i < NN; ++i) {
      xVector[i] = 0.;

      for (int j = 0; j < NN; ++j)
        xVector[i] += -2. * daij11[i][j] * d2ni[j] - d2aij11[i][j] * dni[j];

      for (int j = 0; j < NN; ++j)
        xVector[i] += -2. * daij12[i][j] * d2mus[j] - d2aij12[i][j] * dmus[j];
    }
    for (int i = 0; i < NN; ++i) {
      xVector[NN + i] = 0.;

      for (int j = 0; j < NN; ++j)
        xVector[NN + i] += -2. * daij21[i][j] * d2ni[j] - d2aij21[i][j] * dni[j];

      for (int j = 0; j < NN; ++j)
        xVector[NN + i] += -2. * daij22[i][j] * d2mus[j] - d2aij22[i][j] * dmus[j];
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


  // TODO include correlations
  vector< vector<double> > ThermalModelRealGas::CalculateFluctuations(int order) {
    if (order < 1) return m_chi;

    vector<double> chgs(m_densities.size());
    vector<double> chis;

    // Baryon charge
    for (size_t i = 0; i < chgs.size(); ++i)
      chgs[i] = m_TPS->Particles()[i].BaryonCharge();
    chis = CalculateChargeFluctuations(chgs, order);
    for (int i = 0; i < order; ++i) m_chi[i][0] = chis[i];

    // Electric charge
    for (size_t i = 0; i < chgs.size(); ++i)
      chgs[i] = m_TPS->Particles()[i].ElectricCharge();
    chis = CalculateChargeFluctuations(chgs, order);
    for (int i = 0; i < order; ++i) m_chi[i][1] = chis[i];

    // Strangeness charge
    for (size_t i = 0; i < chgs.size(); ++i)
      chgs[i] = m_TPS->Particles()[i].Strangeness();
    chis = CalculateChargeFluctuations(chgs, order);
    for (int i = 0; i < order; ++i) m_chi[i][2] = chis[i];

    // Arbitrary charge
    for (size_t i = 0; i < chgs.size(); ++i)
      chgs[i] = m_TPS->Particles()[i].ArbitraryCharge();
    chis = CalculateChargeFluctuations(chgs, order);
    for (int i = 0; i < order; ++i) m_chiarb[i] = chis[i];

    return m_chi;
  }

  void ThermalModelRealGas::CalculateTwoParticleCorrelations()
  {
    int NN = m_densities.size();
    int Nevcomp = m_exvolmod->ComponentsNumber();
    const vector<int>& evinds = m_exvolmod->ComponentIndices();
    const vector<int>& evindsfrom = m_exvolmod->ComponentIndicesFrom();


    m_PrimCorrel.resize(NN);
    for (int i = 0; i < NN; ++i)
      m_PrimCorrel[i].resize(NN);
    m_TotalCorrel = m_PrimCorrel;

    MatrixXd densMatrix(2 * NN, 2 * NN);
    VectorXd solVector(2 * NN), xVector(2 * NN);

    vector<double> chi2id(m_densities.size()), Ps(m_densities.size());
    for (int i = 0; i < NN; ++i) {
      //chi2id[i] = m_TPS->Particles()[i].chi(2, m_Parameters, m_UseWidth, m_MuStar[i]);
      Ps[i] = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
      chi2id[i] = m_TPS->Particles()[i].chiDimensionfull(2, m_Parameters, m_UseWidth, m_MuStar[i]);
    }

    vector<double>  evc_Ps(Nevcomp, 0.);
    for (int i = 0; i < NN; ++i) {
      evc_Ps[evinds[i]] += Ps[i];
    }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(i, j) = -m_exvolmod->df(i, j) * m_DensitiesId[i];
        if (i == j) densMatrix(i, j) += 1.;
      }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j)
        densMatrix(i, NN + j) = 0.;

    for (int i = 0; i < NN; ++i) {
      densMatrix(i, NN + i) = -m_exvolmod->f(i) * chi2id[i] * pow(xMath::GeVtoifm(), 3);
    }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(NN + i, j) = m_mfmod->d2v(i, j);
        //for (int k = 0; k < NN; ++k) {
        //  densMatrix(NN + i, j) += -m_exvolmod->d2f(k, i, j) * Ps[k];
        //}
        for (int indk = 0; indk < Nevcomp; ++indk) {
          int k = evindsfrom[indk];
          densMatrix(NN + i, j) += -m_exvolmod->d2f(k, i, j) * evc_Ps[indk];
        }
      }

    for (int i = 0; i < NN; ++i)
      for (int j = 0; j < NN; ++j) {
        densMatrix(NN + i, NN + j) = -m_exvolmod->df(j,i) * m_DensitiesId[j];
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
        dni[i] = solVector[i];
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

  void ThermalModelRealGas::CalculateFluctuations()
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

  double ThermalModelRealGas::CalculateEnergyDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (size_t i = 0; i < m_densities.size(); ++i)
      ret += m_exvolmod->f(i) * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_MuStar[i]);
    ret += m_mfmod->v();

    if (1 /*&& m_TemperatureDependentAB*/) {
      for (size_t i = 0; i < m_densities.size(); ++i) {
        double tPid = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
        ret += -m_Parameters.T * tPid * m_exvolmod->dfdT(i);
      }
      ret += m_Parameters.T * m_mfmod->dvdT();
    }

    return ret;
  }

  double ThermalModelRealGas::CalculateEntropyDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (size_t i = 0; i < m_densities.size(); ++i)
      ret += m_exvolmod->f(i) * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::EntropyDensity, m_UseWidth, m_MuStar[i]);

    if (1 /*&& m_TemperatureDependentAB*/) {
      for (size_t i = 0; i < m_densities.size(); ++i) {
        double tPid = m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
        ret += -m_Parameters.T * tPid * m_exvolmod->dfdT(i);
      }
      ret += m_Parameters.T * m_mfmod->dvdT();
    }

    return ret;
  }

  double ThermalModelRealGas::CalculatePressure() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;
    for (size_t i = 0; i < m_densities.size(); ++i) {
      double tfsum = 0.;
      for (size_t j = 0; j < m_densities.size(); ++j) {
        tfsum += m_densities[j] * m_exvolmod->df(i, j);
      }
      tfsum = m_exvolmod->f(i) - tfsum;
      ret += tfsum * m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_MuStar[i]);
    }

    ret += -m_mfmod->v();
    for (size_t i = 0; i < m_densities.size(); ++i)
      ret += m_mfmod->dv(i) * m_densities[i];

    return ret;
  }

  double ThermalModelRealGas::ParticleScalarDensity(int part) {
    if (!m_Calculated) CalculateDensities();

    return m_scaldens[part];
  }

  double ThermalModelRealGas::MuShift(int id) const
  {
    if (id >= 0. && id < static_cast<int>(ComponentsNumber()))
      return m_MuStar[id] - m_Chem[id];
    else
      return 0.0;
  }

  std::vector<double> ThermalModelRealGas::BroydenEquationsRealGas::Equations(const std::vector<double>& x)
  {
    int NN = m_THM->Densities().size();
    vector<double> Ps(NN, 0.);
    for (int i = 0; i < NN; ++i) {
      Ps[i] = m_THM->TPS()->Particles()[i].Density(m_THM->Parameters(),
        IdealGasFunctions::Pressure,
        m_THM->UseWidth(),
        m_THM->ChemicalPotential(i) + x[i]
      );
    }

    vector<double> ns(NN, 0.);
    for (int i = 0; i < NN; ++i) {
      ns[i] = m_THM->TPS()->Particles()[i].Density(m_THM->Parameters(),
        IdealGasFunctions::ParticleDensity,
        m_THM->UseWidth(),
        m_THM->ChemicalPotential(i) + x[i]
      );
    }

    vector<double> np = m_THM->m_exvolmod->nsol(ns);
    m_THM->m_exvolmod->SetDensities(np);
    m_THM->m_mfmod->SetDensities(np);

    vector<double> ret(m_N, 0.);
    for (size_t i = 0; i < ret.size(); ++i) {
      ret[i] = x[i];
      for (int j = 0; j < NN; ++j)
        ret[i] += -m_THM->m_exvolmod->df(j, i) * Ps[j];
      ret[i] += m_THM->m_mfmod->dv(i);
    }
    return ret;
  }

  std::vector<double> ThermalModelRealGas::BroydenJacobianRealGas::Jacobian(const std::vector<double>& x)
  {
    int NN = m_THM->m_densities.size();

    MatrixXd densMatrix(NN, NN);
    VectorXd solVector(NN), xVector(NN);

    std::vector<double> ret(NN * NN, 0.);
    {
      vector<double> Ps(NN, 0.);
      for (int i = 0; i < NN; ++i)
        Ps[i] = m_THM->TPS()->Particles()[i].Density(m_THM->Parameters(),
          IdealGasFunctions::Pressure,
          m_THM->UseWidth(),
          m_THM->ChemicalPotential(i) + x[i]
        );

      vector<double> ns(NN, 0.);
      for (int i = 0; i < NN; ++i)
        ns[i] = m_THM->TPS()->Particles()[i].Density(m_THM->Parameters(),
          IdealGasFunctions::ParticleDensity,
          m_THM->UseWidth(),
          m_THM->ChemicalPotential(i) + x[i]
        );

      vector<double> chi2s(NN, 0.);
      for (int i = 0; i < NN; ++i)
        chi2s[i] = m_THM->TPS()->Particles()[i].chiDimensionfull(2, m_THM->Parameters(),
          m_THM->UseWidth(),
          m_THM->ChemicalPotential(i) + x[i]
        );

      vector<double> np = m_THM->m_exvolmod->nsol(ns);
      m_THM->m_exvolmod->SetDensities(np);
      m_THM->m_mfmod->SetDensities(np);

      for (int i = 0; i < NN; ++i) {
        for (int j = 0; j < NN; ++j) {
          densMatrix(i, j) = 0.;
          if (i == j)
            densMatrix(i, j) += 1.;

          densMatrix(i, j) += -m_THM->m_exvolmod->df(i, j) * ns[i];
        }
      }

      int Nevcomp = m_THM->m_exvolmod->ComponentsNumber();
      const vector<int>& evinds = m_THM->m_exvolmod->ComponentIndices();
      const vector<int>& evindsfrom = m_THM->m_exvolmod->ComponentIndicesFrom();

      int Nmfcomp = m_THM->m_mfmod->ComponentsNumber();
      const vector<int>& mfinds = m_THM->m_mfmod->ComponentIndices();
      const vector<int>& mfindsfrom = m_THM->m_mfmod->ComponentIndicesFrom();

      vector<double> evc_Ps(Nevcomp, 0.);
      for (int i = 0; i < NN; ++i)
        evc_Ps[m_THM->m_exvolmod->ComponentIndices()[i]] += Ps[i];

      PartialPivLU<MatrixXd> decomp(densMatrix);

      for (int kp = 0; kp < NN; ++kp) {

        if (1 /*attrfl*/) {
          for (int l = 0; l < NN; ++l) {
            xVector[l] = 0.;
            if (l == kp) {
              xVector[l] = chi2s[l] * pow(xMath::GeVtoifm(), 3) * m_THM->m_exvolmod->f(l);
            }
          }

          solVector = decomp.solve(xVector);
          for (int i = 0; i < NN; ++i)
            if (solVector[i] > 1.) solVector[i] = 1.;  // Stabilizer
        }

        vector<double> dnjdmukp(NN, 0.);
        vector<double> evc_dnjdmukp(Nevcomp, 0.);
        vector<double> evc_d2fmul(Nevcomp, 0.);
        vector<double> mfc_dnjdmukp(Nmfcomp, 0.);
        vector<double> mfc_d2vmul(Nmfcomp, 0.);
        if (1 /*attrfl*/) {
          for (int j = 0; j < NN; ++j) {
            dnjdmukp[j] = solVector[j];
            evc_dnjdmukp[evinds[j]] += solVector[j];
            mfc_dnjdmukp[mfinds[j]] += solVector[j];
          }
          for (int nn = 0; nn < Nevcomp; ++nn) {
            for (int in = 0; in < Nevcomp; ++in) {
              for (int inp = 0; inp < Nevcomp; ++inp) {
                int n  = evindsfrom[in];
                int np = evindsfrom[inp];
                int k  = evindsfrom[nn];
                evc_d2fmul[nn] += -evc_Ps[in] * m_THM->m_exvolmod->d2f(n, k, np) * evc_dnjdmukp[inp];
              }
            }
          }

          for (int nn = 0; nn < Nmfcomp; ++nn) {
            for (int in = 0; in < Nmfcomp; ++in) {
              int n = mfindsfrom[in];
              int k = mfindsfrom[nn];
              mfc_d2vmul[nn] += m_THM->m_mfmod->d2v(k, n) * mfc_dnjdmukp[in];
            }
          }
        }


        for (int k = 0; k < NN; ++k) {
          if (k == kp)
            ret[k * NN + kp] += 1.;
          ret[k * NN + kp] += -m_THM->m_exvolmod->df(kp, k) * ns[kp];

          if (1 /*attrfl*/) {
            //for (int n = 0; n < NN; ++n) {
            //  for (int np = 0; np < NN; ++np) {
            //    ret[k * NN + kp] += -Ps[n] * m_THM->m_exvolmod->d2f(n, k, np) * dnjdmukp[np];
            //  }
            //}

            ret[k * NN + kp] += evc_d2fmul[evinds[k]];
            //for (int in = 0; in < Nevcomp; ++in) {
            //  for (int inp = 0; inp < Nevcomp; ++inp) {
            //    int n  = m_THM->m_exvolmod->ComponentIndicesFrom()[in];
            //    int np = m_THM->m_exvolmod->ComponentIndicesFrom()[inp];
            //    ret[k * NN + kp] += -evc_Ps[in] * m_THM->m_exvolmod->d2f(n, k, np) * evc_dnjdmukp[inp];
            //  }
            //}

            //for (int n = 0; n < NN; ++n) {
            //  ret[k * NN + kp] += m_THM->m_mfmod->d2v(k, n) * dnjdmukp[n];
            //}

            ret[k * NN + kp] += mfc_d2vmul[mfinds[k]];

            //for (int in = 0; in < Nmfcomp; ++in) {
            //  int n = m_THM->m_mfmod->ComponentIndicesFrom()[in];
            //  ret[k * NN + kp] += m_THM->m_mfmod->d2v(k, n) * mfc_dnjdmukp[in];
            //}
          }
        }
      }
    }

    return ret;
  }

  bool ThermalModelRealGas::BroydenSolutionCriteriumRealGas::IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta) const
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

}