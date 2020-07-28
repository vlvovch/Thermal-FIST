/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ThermalModelCanonical.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"

using namespace std;

namespace thermalfist {

  ThermalModelCanonical::ThermalModelCanonical(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
    ThermalModelBase(TPS_, params), m_BCE(1), m_QCE(1), m_SCE(1), m_CCE(1), m_IntegrationIterationsMultiplier(1)
  {

    m_TAG = "ThermalModelCanonical";

    m_Ensemble = CE;
    m_InteractionModel = Ideal;

    m_modelgce = NULL;

    m_Banalyt = false;
  }


  ThermalModelCanonical::~ThermalModelCanonical(void)
  {
    CleanModelGCE();
  }

  void ThermalModelCanonical::ChangeTPS(ThermalParticleSystem *TPS_) {
    ThermalModelBase::ChangeTPS(TPS_);
  }

  void ThermalModelCanonical::CalculateQuantumNumbersRange(bool computeFluctuations)
  {
    m_BMAX = 0;
    m_QMAX = 0;
    m_SMAX = 0;
    m_CMAX = 0;


    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ThermalParticle &part = m_TPS->Particle(i);

      if (part.Statistics() != 0 && IsParticleCanonical(part)) {
        part.SetCalculationType(IdealGasFunctions::ClusterExpansion);

        m_BMAX = max(m_BMAX, abs(part.BaryonCharge() * part.ClusterExpansionOrder()));
        m_QMAX = max(m_QMAX, abs(part.ElectricCharge() * part.ClusterExpansionOrder()));
        m_SMAX = max(m_SMAX, abs(part.Strangeness() * part.ClusterExpansionOrder()));
        m_CMAX = max(m_CMAX, abs(part.Charm() * part.ClusterExpansionOrder()));
      }
      else {
        m_BMAX = max(m_BMAX, abs(part.BaryonCharge()));
        m_QMAX = max(m_QMAX, abs(part.ElectricCharge()));
        m_SMAX = max(m_SMAX, abs(part.Strangeness()));
        m_CMAX = max(m_CMAX, abs(part.Charm()));
      }
    }

    m_BMAX_list = m_BMAX;
    m_QMAX_list = m_QMAX;
    m_SMAX_list = m_SMAX;
    m_CMAX_list = m_CMAX;

    if (computeFluctuations) {
      m_BMAX *= 2;
      m_QMAX *= 2;
      m_SMAX *= 2;
      m_CMAX *= 2;
    }

    // Some charges may be treated grand-canonically
    m_BMAX *= m_BCE;
    m_QMAX *= m_QCE;
    m_SMAX *= m_SCE;
    m_CMAX *= m_CCE;

    printf("BMAX = %d\tQMAX = %d\tSMAX = %d\tCMAX = %d\n", m_BMAX, m_QMAX, m_SMAX, m_CMAX);

    m_QNMap.clear();
    m_QNvec.resize(0);

    m_Corr.resize(0);
    m_PartialZ.resize(0);

    int ind = 0;
    for (int iB = -m_BMAX; iB <= m_BMAX; ++iB)
      for (int iQ = -m_QMAX; iQ <= m_QMAX; ++iQ)
        for (int iS = -m_SMAX; iS <= m_SMAX; ++iS)
          for (int iC = -m_CMAX; iC <= m_CMAX; ++iC) {

            QuantumNumbers qn(iB, iQ, iS, iC);
            m_QNMap[qn] = ind;
            m_QNvec.push_back(qn);

            m_PartialZ.push_back(0.);
            m_Corr.push_back(1.);


            ind++;
          }

  }

  void ThermalModelCanonical::SetStatistics(bool stats) {
    m_QuantumStats = stats;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      m_TPS->Particle(i).UseStatistics(stats);

      // Only cluster expansion method supported for particles with canonically conserved charges
      if (stats && IsParticleCanonical(m_TPS->Particle(i)))
        m_TPS->Particle(i).SetCalculationType(IdealGasFunctions::ClusterExpansion);
    }
    m_PartialZ.clear();
    //CalculateQuantumNumbersRange();
  }

  void ThermalModelCanonical::FixParameters()
  {
    m_ConstrainMuB &= !m_BCE;
    m_ConstrainMuQ &= !m_QCE;
    m_ConstrainMuS &= !m_SCE;
    m_ConstrainMuC &= !m_CCE;
    ThermalModelBase::FixParameters();
  }

  void ThermalModelCanonical::FixParametersNoReset()
  {
    m_ConstrainMuB &= !m_BCE;
    m_ConstrainMuQ &= !m_QCE;
    m_ConstrainMuS &= !m_SCE;
    m_ConstrainMuC &= !m_CCE;
    ThermalModelBase::FixParametersNoReset();
  }


  void ThermalModelCanonical::CalculatePrimordialDensities() {
    m_FluctuationsCalculated = false;

    if (m_PartialZ.size() == 0)
      CalculateQuantumNumbersRange();

    if (m_BMAX_list == 1 && m_BCE && m_QCE && m_SCE && m_CCE && !UsePartialChemicalEquilibrium()) {
      m_Banalyt = true;
      m_Parameters.muB = 0.0;
      m_Parameters.muQ = 0.0;
      m_Parameters.muS = 0.0;
      m_Parameters.muC = 0.0;
    }
    else {
      m_Banalyt = false;
      if (m_BCE)
        m_Parameters.muB = 0.0;
      if (m_QCE)
        m_Parameters.muQ = 0.0;
      if (m_SCE)
        m_Parameters.muS = 0.0;
      if (m_CCE)
        m_Parameters.muC = 0.0;

      //PrepareModelGCE(); // Plan B, may work better when quantum numbers are large
    }

    CalculatePartitionFunctions();

    for (size_t i = 0; i < m_densities.size(); ++i) {
      ThermalParticle &tpart = m_TPS->Particle(i);
      m_densities[i] = 0.;

      if (!IsParticleCanonical(tpart)) {
        m_densities[i] = tpart.Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
      }
      else if (tpart.Statistics() == 0
        || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion)
      {
        int ind = m_QNMap[QuantumNumbers(m_BCE * tpart.BaryonCharge(), m_QCE * tpart.ElectricCharge(), m_SCE * tpart.Strangeness(), m_CCE * tpart.Charm())];

        if (ind < static_cast<int>(m_Corr.size()))
          m_densities[i] = m_Corr[ind] * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
      }
      else {
        for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
          int ind = m_QNMap[QuantumNumbers(m_BCE*n*tpart.BaryonCharge(), m_QCE*n*tpart.ElectricCharge(), m_SCE*n*tpart.Strangeness(), m_CCE*n*tpart.Charm())];
          if (ind < static_cast<int>(m_Corr.size()))
            m_densities[i] += m_Corr[ind] * tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
        }
      }
    }

    m_Calculated = true;
    ValidateCalculation();
  }

  void ThermalModelCanonical::ValidateCalculation()
  {
    ThermalModelBase::ValidateCalculation();

    char cc[1000];

    double TOL = 1.e-4;

    // Checking that the CE calculation is valid
    if (m_BCE && m_BMAX != 0) {
      double totB = CalculateBaryonDensity() * m_Parameters.SVc;
      if (fabs(m_Parameters.B - totB) > TOL) {

        sprintf(cc, "**WARNING** ThermalModelCanonical: Inaccurate calculation of total baryon number.\
\n\
Expected: %d\n\
Obtained: %lf\n\
\n", m_Parameters.B, totB);

        printf("%s", cc);

        m_ValidityLog.append(cc);

        m_LastCalculationSuccessFlag = false;
      }
    }


    if (m_QCE && m_QMAX != 0) {
      double totQ = CalculateChargeDensity() * m_Parameters.SVc;
      if (fabs(m_Parameters.Q - totQ) > TOL) {
        sprintf(cc, "**WARNING** ThermalModelCanonical: Inaccurate calculation of total electric charge.\
\n\
Expected: %d\n\
Obtained: %lf\n\
\n", m_Parameters.Q, totQ);

        printf("%s", cc);

        m_ValidityLog.append(cc);

        m_LastCalculationSuccessFlag = false;
      }
    }


    if (m_SCE && m_SMAX != 0) {
      double totS = CalculateStrangenessDensity() * m_Parameters.SVc;
      if (fabs(m_Parameters.S - totS) > TOL) {
        sprintf(cc, "**WARNING** ThermalModelCanonical: Inaccurate calculation of total strangeness.\
\n\
Expected: %d\n\
Obtained: %lf\n\
\n", m_Parameters.S, totS);

        printf("%s", cc);

        m_ValidityLog.append(cc);

        m_LastCalculationSuccessFlag = false;
      }
    }


    if (m_CCE && m_CMAX != 0) {
      double totC = CalculateCharmDensity() * m_Parameters.SVc;
      if (fabs(m_Parameters.C - totC) > TOL) {
        sprintf(cc, "**WARNING** ThermalModelCanonical: Inaccurate calculation of total charm.\
\n\
Expected: %d\n\
Obtained: %lf\n\
\n", m_Parameters.C, totC);

        printf("%s", cc);

        m_ValidityLog.append(cc);

        m_LastCalculationSuccessFlag = false;
      }
    }
  }

  void ThermalModelCanonical::CalculatePartitionFunctions(double Vc)
  {
    if (Vc < 0.0)
      Vc = m_Parameters.SVc;

    if (!UsePartialChemicalEquilibrium()) 
      FillChemicalPotentials();
    else {
      // Partial chemical equilibrium canonical ensemble currently works only if there is particle-antiparticle symmetry
      for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
        int i2 = m_TPS->PdgToId(-m_TPS->Particle(i).PdgId());
        if (i2 != -1) {
          if (fabs(m_Chem[i] - m_Chem[i2]) > 1.e-8) {
            printf("**ERROR** ThermalModelCanonical::CalculatePartitionFunctions: Partial chemical equilibrium canonical ensemble only supported if particle-antiparticle fugacities are symmetric!\n");
            exit(1);
          }
        }
      }
    }

    bool AllMuZero = true;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ThermalParticle &tpart = m_TPS->Particle(i);
      if (IsParticleCanonical(tpart) && m_Chem[i] != 0.0)
      {
        AllMuZero = false;
        break;
      }
    }

    vector<double> Nsx(m_PartialZ.size(), 0.);
    vector<double> Nsy(m_PartialZ.size(), 0.);

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ThermalParticle &tpart = m_TPS->Particle(i);

      if (!IsParticleCanonical(tpart)) {
        int ind = m_QNMap[QuantumNumbers(m_BCE * tpart.BaryonCharge(), m_QCE * tpart.ElectricCharge(), m_SCE * tpart.Strangeness(), m_CCE * tpart.Charm())];
        if (ind != m_QNMap[QuantumNumbers(0, 0, 0, 0)]) {
          printf("**ERROR** ThermalModelCanonical: neutral particle cannot have non-zero ce charges\n");
          exit(1);
        }
        if (ind < static_cast<int>(Nsx.size()))
          Nsx[ind] += tpart.Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
      }
      else if (tpart.Statistics() == 0
        || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion) {
        int ind = m_QNMap[QuantumNumbers(m_BCE * tpart.BaryonCharge(), m_QCE * tpart.ElectricCharge(), m_SCE * tpart.Strangeness(), m_CCE * tpart.Charm())];
        if (ind < static_cast<int>(Nsx.size())) {
          if (!UsePartialChemicalEquilibrium()) {
            double tdens = tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
            Nsx[ind] += tdens * cosh(m_Chem[i] / m_Parameters.T);
            Nsy[ind] += tdens * sinh(m_Chem[i] / m_Parameters.T);
          }
          // Currently only works at mu = 0!!
          else {
            double tdens = tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
            Nsx[ind] += tdens;
          }
        }
      }
      else {
        for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
          int ind = m_QNMap[QuantumNumbers(m_BCE*n*tpart.BaryonCharge(), m_QCE*n*tpart.ElectricCharge(), m_SCE*n*tpart.Strangeness(), m_CCE*n*tpart.Charm())];
          if (ind < static_cast<int>(Nsx.size())) {
            if (!UsePartialChemicalEquilibrium()) {
              double tdens = tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.) / static_cast<double>(n); // TODO: Check
              Nsx[ind] += tdens * cosh(n * m_Chem[i] / m_Parameters.T);
              Nsy[ind] += tdens * sinh(n * m_Chem[i] / m_Parameters.T);
            }
            // Currently only works at mu = 0!!
            else {
              double tdens = tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]) / static_cast<double>(n); // TODO: Check
              Nsx[ind] += tdens;
            }
          }
        }
      }
    }

    for (int i = 0; i < static_cast<int>(Nsx.size()); ++i) {
      Nsx[i] *= Vc;
      Nsy[i] *= Vc;
      m_PartialZ[i] = 0.;
    }

    int nmax = max(3, (int)sqrt(m_Parameters.B*m_Parameters.B + m_Parameters.Q*m_Parameters.Q + m_Parameters.S*m_Parameters.S + m_Parameters.C*m_Parameters.C));
    if (m_Parameters.B == 0 && m_Parameters.Q == 0 && m_Parameters.S == 0 && m_Parameters.C == 0)
      nmax = 4;


    
    int nmaxB = max(4, m_Parameters.B);
    int nmaxQ = max(4, m_Parameters.Q);
    int nmaxS = max(4, m_Parameters.S);
    int nmaxC = max(4, m_Parameters.C);

    // UPDATE: allow to increase the number of interations externally
    nmax *= m_IntegrationIterationsMultiplier;
    nmaxB = nmaxQ = nmaxS = nmaxC = nmax;





    m_MultExp = 0.;
    m_MultExpBanalyt = 0.;
    for (size_t i = 0; i < m_PartialZ.size(); ++i) {
      if (!m_Banalyt || m_QNvec[i].B == 0)
        m_MultExp += Nsx[i];
      if (m_Banalyt && (m_QNvec[i].B == 1 || m_QNvec[i].B == -1))
        m_MultExpBanalyt += Nsx[i];
    }

    double dphiB = xMath::Pi() / nmaxB;
    int maxB = 2 * nmaxB;
    if (m_BMAX == 0 || m_Banalyt)
      maxB = 1;

    for (int iB = 0; iB < maxB; ++iB) {

      vector<double> xlegB, wlegB;

      if (m_BMAX != 0 && !m_Banalyt) {
        double aB = iB * dphiB;
        if (iB >= nmaxB) aB = xMath::Pi() - (iB + 1) * dphiB;
        double bB = aB + dphiB;
        NumericalIntegration::GetCoefsIntegrateLegendre10(aB, bB, &xlegB, &wlegB);
      }
      else {
        xlegB.resize(1);
        xlegB[0] = 0.;
        wlegB.resize(1);
        wlegB[0] = 1.;
      }


      double dphiS = xMath::Pi() / nmaxS;
      int maxS = 2 * nmaxS;
      if (m_SMAX == 0)
        maxS = 1;

      for (int iS = 0; iS < maxS; ++iS) {
        vector<double> xlegS, wlegS;

        if (m_SMAX != 0) {
          double aS = iS * dphiS;
          if (iS >= nmaxS) aS = xMath::Pi() - (iS + 1) * dphiS;
          double bS = aS + dphiS;
          NumericalIntegration::GetCoefsIntegrateLegendre10(aS, bS, &xlegS, &wlegS);
        }
        else {
          xlegS.resize(1);
          xlegS[0] = 0.;
          wlegS.resize(1);
          wlegS[0] = 1.;
        }

        double dphiQ = xMath::Pi() / nmaxQ;
        int maxQ = nmaxQ;
        if (m_QMAX == 0)
          maxQ = 1;

        for (int iQ = 0; iQ < maxQ; ++iQ) {
          vector<double> xlegQ, wlegQ;

          if (m_QMAX != 0) {
            double aQ = iQ * dphiQ;
            double bQ = aQ + dphiQ;
            NumericalIntegration::GetCoefsIntegrateLegendre10(aQ, bQ, &xlegQ, &wlegQ);
          }
          else {
            xlegQ.resize(1);
            xlegQ[0] = 0.;
            wlegQ.resize(1);
            wlegQ[0] = 1.;
          }

          double dphiC = xMath::Pi() / nmaxC;
          int maxC = 2 * nmaxC;
          if (m_CMAX == 0)
            maxC = 1;

          for (int iC = 0; iC < maxC; ++iC) {
            vector<double> xlegC, wlegC;
            if (m_CMAX != 0) {
              double aC = iC * dphiC;
              if (iC >= nmaxC) aC = xMath::Pi() - (iC + 1) * dphiC;
              double bC = aC + dphiC;
              NumericalIntegration::GetCoefsIntegrateLegendre10(aC, bC, &xlegC, &wlegC);
            }
            else {
              xlegC.resize(1);
              xlegC[0] = 0.;
              wlegC.resize(1);
              wlegC[0] = 1.;
            }

            for (size_t iBt = 0; iBt < xlegB.size(); ++iBt) {
              for (size_t iSt = 0; iSt < xlegS.size(); ++iSt) {
                for (size_t iQt = 0; iQt < xlegQ.size(); ++iQt) {
                  for (size_t iCt = 0; iCt < xlegC.size(); ++iCt) {
                    vector<double> cosph(m_PartialZ.size(), 0.), sinph(m_PartialZ.size(), 0.);
                    double wx = 0., wy = 0., mx = 0., my = 0.;
                    for (size_t i = 0; i < m_PartialZ.size(); ++i) {
                      int tB = m_QNvec[i].B;
                      int tQ = m_QNvec[i].Q;
                      int tS = m_QNvec[i].S;
                      int tC = m_QNvec[i].C;

                      if (m_Banalyt) {
                        cosph[i] = cos(tS*xlegS[iSt] + tQ * xlegQ[iQt] + tC * xlegC[iCt]);
                        sinph[i] = sin(tS*xlegS[iSt] + tQ * xlegQ[iQt] + tC * xlegC[iCt]);

                        if (m_QNvec[i].B == 1) {
                          wx += Nsx[i] * cosph[i];
                          wy += Nsx[i] * sinph[i];
                        }
                        else if (m_QNvec[i].B == 0) {
                          mx += Nsx[i] * (cosph[i] - 1.);
                        }
                      }
                      else {
                        cosph[i] = cos(tB*xlegB[iBt] + tS * xlegS[iSt] + tQ * xlegQ[iQt] + tC * xlegC[iCt]);
                        mx += Nsx[i] * (cosph[i] - 1.);

                        if (!AllMuZero) {
                          sinph[i] = sin(tB*xlegB[iBt] + tS * xlegS[iSt] + tQ * xlegQ[iQt] + tC * xlegC[iCt]);
                          my += Nsy[i] * sinph[i];
                        }
                      }
                    }

                    double wmod = 0.;
                    double warg = 0.;

                    if (m_Banalyt) {
                      wmod = sqrt(wx*wx + wy * wy);
                      warg = atan2(wy, wx);
                    }

                    for (size_t iN = 0; iN < m_PartialZ.size(); ++iN) {
                      int tBg = m_Parameters.B - m_QNvec[iN].B;
                      int tQg = m_Parameters.Q - m_QNvec[iN].Q;
                      int tSg = m_Parameters.S - m_QNvec[iN].S;
                      int tCg = m_Parameters.C - m_QNvec[iN].C;

                      if (m_Banalyt) {
                        //m_PartialZ[iN] += wlegB[iBt] * wlegS[iSt] * wlegQ[iQt] * wlegC[iCt] * cos(tBg*xlegB[iBt] + tSg * xlegS[iSt] + tQg * xlegQ[iQt] + tCg * xlegC[iCt] - tBg * warg) * exp(mx) * xMath::BesselI(tBg, 2. * wmod);
                        m_PartialZ[iN] += wlegB[iBt] * wlegS[iSt] * wlegQ[iQt] * wlegC[iCt] * 
                          cos(tBg * xlegB[iBt] + tSg * xlegS[iSt] + tQg * xlegQ[iQt] + tCg * xlegC[iCt] - tBg * warg) * 
                          exp(mx + 2. * wmod - m_MultExpBanalyt) * 
                          xMath::BesselIexp(tBg, 2. * wmod);
                      }
                      else {
                        if (AllMuZero)
                          m_PartialZ[iN] += wlegB[iBt] * wlegS[iSt] * wlegQ[iQt] * wlegC[iCt] * cos(tBg*xlegB[iBt] + tSg * xlegS[iSt] + tQg * xlegQ[iQt] + tCg * xlegC[iCt]) * exp(mx);
                        else
                          m_PartialZ[iN] += wlegB[iBt] * wlegS[iSt] * wlegQ[iQt] * wlegC[iCt] * exp(mx)
                          * (cos(tBg*xlegB[iBt] + tSg * xlegS[iSt] + tQg * xlegQ[iQt] + tCg * xlegC[iCt]) * cos(my)
                            + sin(tBg*xlegB[iBt] + tSg * xlegS[iSt] + tQg * xlegQ[iQt] + tCg * xlegC[iCt]) * sin(my));
                      }
                    }
                  }
                }
              }
            }

            //int cind = iB * maxS * maxQ * maxC + iS * maxQ * maxC + iQ * maxC + iC;
            //int tot = maxB * maxS * maxQ * maxC;
          }
        }
      }
    }

    for (size_t iN = 0; iN < m_PartialZ.size(); ++iN) {
      if (m_BMAX != 0 && m_BMAX != 1 && !m_Banalyt) // TODO: cross-check
        m_PartialZ[iN] /= 2. * xMath::Pi();
      if (m_QMAX != 0) {
        m_PartialZ[iN] /= 2. * xMath::Pi();
        m_PartialZ[iN] *= 2.; /// TODO: Extra cross-check the factor 2, can be important for entropy density, irrelevant for everything else
      }
      if (m_SMAX != 0)
        m_PartialZ[iN] /= 2. * xMath::Pi();
      if (m_CMAX != 0)
        m_PartialZ[iN] /= 2. * xMath::Pi();
    }


    m_Corr.resize(m_PartialZ.size());
    for (size_t iN = 0; iN < m_PartialZ.size(); ++iN) {
      m_Corr[iN] = m_PartialZ[iN] / m_PartialZ[m_QNMap[QuantumNumbers(0, 0, 0, 0)]];
    }
  }

  double ThermalModelCanonical::ParticleScaledVariance(int part)
  {
    ThermalParticle &tpart = m_TPS->Particle(part);
    double ret1 = 0., ret2 = 0., ret3 = 0.;

    if (!IsParticleCanonical(tpart)) {
      return tpart.ScaledVariance(m_Parameters, m_UseWidth, m_Chem[part]);
    }
    else if (tpart.Statistics() == 0
      || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion)
    {
      int ind = m_QNMap[QuantumNumbers(m_BCE*tpart.BaryonCharge(), m_QCE*tpart.ElectricCharge(), m_SCE*tpart.Strangeness(), m_CCE*tpart.Charm())];
      int ind2 = m_QNMap[QuantumNumbers(m_BCE * 2 * tpart.BaryonCharge(), m_QCE * 2 * tpart.ElectricCharge(), m_SCE * 2 * tpart.Strangeness(), m_CCE * 2 * tpart.Charm())];

      ret1 = 1.;
      if (ind < static_cast<int>(m_Corr.size()) && ind2 < static_cast<int>(m_Corr.size()))
        ret2 = m_Corr[ind2] / m_Corr[ind] * m_Parameters.SVc * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[part]);

      if (ind < static_cast<int>(m_Corr.size()))
        ret3 = -m_Corr[ind] * m_Parameters.SVc * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[part]);
    }
    else {
      double ret1num = 0., ret1zn = 0.;
      for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
        int ind = m_QNMap[QuantumNumbers(m_BCE*n*tpart.BaryonCharge(), m_QCE*n*tpart.ElectricCharge(), m_SCE*n*tpart.Strangeness(), m_CCE*n*tpart.Charm())];

        double densityClusterN = tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[part]);

        if (ind < static_cast<int>(m_Corr.size())) {
          ret1num += m_Corr[ind] * n * densityClusterN;
          ret1zn += m_Corr[ind] * densityClusterN;
        }

        for (int n2 = 1; n2 <= tpart.ClusterExpansionOrder(); ++n2) {
          if (m_QNMap.count(QuantumNumbers(m_BCE*(n + n2)*tpart.BaryonCharge(), m_QCE*(n + n2)*tpart.ElectricCharge(), m_SCE*(n + n2)*tpart.Strangeness(), m_CCE*(n + n2)*tpart.Charm())) != 0) {
            int ind2 = m_QNMap[QuantumNumbers(m_BCE*(n + n2)*tpart.BaryonCharge(), m_QCE*(n + n2)*tpart.ElectricCharge(), m_SCE*(n + n2)*tpart.Strangeness(), m_CCE*(n + n2)*tpart.Charm())];
            if (ind < static_cast<int>(m_Corr.size()) && ind2 < static_cast<int>(m_Corr.size()))
              ret2 += densityClusterN * m_Corr[ind2] * m_Parameters.SVc * tpart.DensityCluster(n2, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[part]);
          }
        }
      }

      if (ret1zn == 0.0)
        return 1.;

      ret1 = ret1num / ret1zn;
      ret2 = ret2 / ret1zn;
      ret3 = -ret1zn * m_Parameters.SVc;
    }
    return ret1 + ret2 + ret3;
  }

  void ThermalModelCanonical::CalculateTwoParticleCorrelations() {
    int NN = m_densities.size();

    vector<double> yld(NN, 0);
    vector<double> ret1num(NN, 0);
    vector< vector<double> > ret2num(NN, vector<double>(NN, 0.));

    for (int i = 0; i < NN; ++i)
      yld[i] = m_densities[i] * m_Parameters.SVc;

    for (int i = 0; i < NN; ++i) {
      ThermalParticle &tpart = m_TPS->Particle(i);
      if (!IsParticleCanonical(tpart)) {
        ret1num[i] = tpart.ScaledVariance(m_Parameters, m_UseWidth, m_Chem[i]) * yld[i];
      }
      else if (tpart.Statistics() == 0
        || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion)
      {
        ret1num[i] = yld[i];
      }
      else {
        for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
          int ind = m_QNMap[QuantumNumbers(m_BCE*n*tpart.BaryonCharge(), m_QCE*n*tpart.ElectricCharge(), m_SCE*n*tpart.Strangeness(), m_CCE*n*tpart.Charm())];

          double densityClusterN = tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);

          if (ind < static_cast<int>(m_Corr.size()))
            ret1num[i] += m_Corr[ind] * n * densityClusterN * m_Parameters.SVc;
        }
      }
    }

    for (int i = 0; i < NN; ++i) {
      for (int j = 0; j < NN; ++j) {
        ThermalParticle &tpart1 = m_TPS->Particle(i);
        ThermalParticle &tpart2 = m_TPS->Particle(j);

        int n1max = tpart1.ClusterExpansionOrder();
        int n2max = tpart2.ClusterExpansionOrder();

        if (tpart1.Statistics() == 0 || tpart1.CalculationType() != IdealGasFunctions::ClusterExpansion)
          n1max = 1;
        if (tpart2.Statistics() == 0 || tpart2.CalculationType() != IdealGasFunctions::ClusterExpansion)
          n2max = 1;

        if (!IsParticleCanonical(tpart1) || !IsParticleCanonical(tpart2)) {
          ret2num[i][j] = yld[i] * yld[j];
        }
        else {
          for (int n1 = 1; n1 <= n1max; ++n1) {
            for (int n2 = 1; n2 <= n2max; ++n2) {
              int ind = m_QNMap[QuantumNumbers(
                m_BCE*(n1*tpart1.BaryonCharge() + n2 * tpart2.BaryonCharge()),
                m_QCE*(n1*tpart1.ElectricCharge() + n2 * tpart2.ElectricCharge()),
                m_SCE*(n1*tpart1.Strangeness() + n2 * tpart2.Strangeness()),
                m_CCE*(n1*tpart1.Charm() + n2 * tpart2.Charm()))];

              double densityClusterN1 = tpart1.DensityCluster(n1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
              double densityClusterN2 = tpart2.DensityCluster(n2, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[j]);

              if (ind < static_cast<int>(m_Corr.size()))
                ret2num[i][j] += m_Corr[ind] * densityClusterN1 * densityClusterN2 * m_Parameters.SVc * m_Parameters.SVc;
            }
          }
        }
      }
    }


    m_PrimCorrel.resize(NN);
    for (int i = 0; i < NN; ++i)
      m_PrimCorrel[i].resize(NN);
    m_TotalCorrel = m_PrimCorrel;

    for (int i = 0; i < NN; ++i) {
      for (int j = 0; j < NN; ++j) {
        m_PrimCorrel[i][j] = 0.;
        if (i == j)
          m_PrimCorrel[i][j] += ret1num[i] / yld[i];
        m_PrimCorrel[i][j] += ret2num[i][j] / yld[i];
        m_PrimCorrel[i][j] += -yld[j];
        m_PrimCorrel[i][j] *= yld[i] / m_Parameters.SVc / m_Parameters.T;
        if (yld[i] == 0.0)
          m_PrimCorrel[i][j] = 0.;
      }
    }

    for (int i = 0; i < NN; ++i) {
      m_wprim[i] = m_PrimCorrel[i][i];
      if (m_densities[i] > 0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
      else m_wprim[i] = 1.;
    }

  }

  void ThermalModelCanonical::CalculateFluctuations() {
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

  double ThermalModelCanonical::CalculateEnergyDensity() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ThermalParticle &tpart = m_TPS->Particle(i);
      {
        if (!IsParticleCanonical(tpart)) {
          ret += tpart.Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);
        }
        else if (tpart.Statistics() == 0
          || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion) {
          int ind = m_QNMap[QuantumNumbers(tpart.BaryonCharge(), tpart.ElectricCharge(), tpart.Strangeness(), tpart.Charm())];

          if (ind < static_cast<int>(m_Corr.size()))
            ret += m_Corr[ind] * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);
        }
        else {
          for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
            int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];
            if (ind < static_cast<int>(m_Corr.size()))
              ret += m_Corr[ind] * tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);
          }
        }
      }
    }

    return ret;
  }

  double ThermalModelCanonical::CalculatePressure() {
    if (!m_Calculated) CalculateDensities();
    double ret = 0.;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ThermalParticle &tpart = m_TPS->Particle(i);
      {
        if (!IsParticleCanonical(tpart)) {
          ret += tpart.Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i]);
        }
        else if (tpart.Statistics() == 0
          || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion) {
          int ind = m_QNMap[QuantumNumbers(tpart.BaryonCharge(), tpart.ElectricCharge(), tpart.Strangeness(), tpart.Charm())];

          if (ind < static_cast<int>(m_Corr.size()))
            ret += m_Corr[ind] * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i]);
        }
        else {
          for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
            int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];

            if (ind < static_cast<int>(m_Corr.size()))
              ret += m_Corr[ind] * tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i]);
          }
        }
      }
    }

    return ret;
  }

  double ThermalModelCanonical::CalculateEntropyDensity()
  {
    double ret = (CalculateEnergyDensity() / m_Parameters.T) + (m_MultExp + m_MultExpBanalyt + log(m_PartialZ[m_QNMap[QuantumNumbers(0, 0, 0, 0)]])) / m_Parameters.SVc;

    if (m_BCE)
      ret += -m_Parameters.muB / m_Parameters.T * m_Parameters.B / m_Parameters.SVc;
    else
      ret += -m_Parameters.muB * CalculateBaryonDensity() / m_Parameters.T;

    if (m_QCE)
      ret += -m_Parameters.muQ / m_Parameters.T * m_Parameters.Q / m_Parameters.SVc;
    else
      ret += -m_Parameters.muQ * CalculateChargeDensity() / m_Parameters.T;

    if (m_SCE)
      ret += -m_Parameters.muS / m_Parameters.T * m_Parameters.S / m_Parameters.SVc;
    else
      ret += -m_Parameters.muS * CalculateStrangenessDensity() / m_Parameters.T;

    if (m_CCE)
      ret += -m_Parameters.muC / m_Parameters.T * m_Parameters.C / m_Parameters.SVc;
    else
      ret += -m_Parameters.muC * CalculateCharmDensity() / m_Parameters.T;

    return ret;
  }

  double ThermalModelCanonical::GetGCEDensity(int i) const
  {
    return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_Chem[i]);
  }

  bool ThermalModelCanonical::IsParticleCanonical(const ThermalParticle & part)
  {
    return !(
      (part.BaryonCharge() == 0 || m_BCE == 0)
      && (part.ElectricCharge() == 0 || m_QCE == 0)
      && (part.Strangeness() == 0 || m_SCE == 0)
      && (part.Charm() == 0 || m_CCE == 0)
      );
  }

  bool ThermalModelCanonical::IsConservedChargeCanonical(ConservedCharge::Name charge) const
  {
    if (charge == ConservedCharge::BaryonCharge)
      return (m_BCE != 0);
    else if (charge == ConservedCharge::ElectricCharge)
      return (m_QCE != 0);
    else if (charge == ConservedCharge::StrangenessCharge)
      return (m_SCE != 0);
    else if (charge == ConservedCharge::CharmCharge)
      return (m_CCE != 0);
    return 0;
  }

  void ThermalModelCanonical::PrepareModelGCE()
  {
    CleanModelGCE();

    m_modelgce = new ThermalModelIdeal(m_TPS, m_Parameters);
    m_modelgce->SetUseWidth(m_UseWidth);
    m_modelgce->SetChemicalPotentials(m_Chem);

    if (m_BCE)
      m_Parameters.muB = 0.0;

    if (!m_BCE && m_SCE)
      m_Parameters.muS = m_Parameters.muB / 3.;

    if (!m_BCE && m_QCE)
      m_Parameters.muQ = -m_Parameters.muB / 30.;

    m_Parameters.muC = 0.;

    m_modelgce->SolveChemicalPotentials(m_Parameters.B, m_Parameters.Q, m_Parameters.S, m_Parameters.C,
      m_Parameters.muB, m_Parameters.muQ, m_Parameters.muS, m_Parameters.muC,
      static_cast<bool>(m_BCE), 
      static_cast<bool>(m_QCE), 
      static_cast<bool>(m_SCE), 
      static_cast<bool>(m_CCE));

    m_Parameters.muB = m_modelgce->Parameters().muB;
    m_Parameters.muQ = m_modelgce->Parameters().muQ;
    m_Parameters.muS = m_modelgce->Parameters().muS;
    m_Parameters.muC = m_modelgce->Parameters().muC;


    // Possible alternative below
    //double tdens = 0.;
    //for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
    //  ThermalParticle &part = m_TPS->Particle(i);
    //  if (part.BaryonCharge() == 1)
    //    tdens += part.Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
    //}
    //m_Parameters.muB = m_Parameters.T * asinh(m_Parameters.B / m_Parameters.SVc / 2. / tdens);
    //m_Parameters.muS = m_Parameters.muB / 3.;
    //m_Parameters.muQ = -m_Parameters.muB / 30.;
  }

  void ThermalModelCanonical::CleanModelGCE()
  {
    if (m_modelgce != NULL) {
      delete m_modelgce;
      m_modelgce = NULL;
    }
  }

} // namespace thermalfist


