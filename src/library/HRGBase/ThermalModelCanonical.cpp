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
#include <cassert>

#include "HRGBase/ThermalModelIdeal.h"
#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"

#include <Eigen/Dense>

using namespace std;

namespace thermalfist {

  ThermalModelCanonical::ThermalModelCanonical(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
    ThermalModelBase(TPS_, params), 
    m_BCE(1), 
    m_QCE(1), 
    m_SCE(1), 
    m_CCE(1),
    m_IntegrationIterationsMultiplier(1),
    m_Method(GaussLegendre)
  {

    m_TAG = "ThermalModelCanonical";

    m_Ensemble = CE;
    m_InteractionModel = Ideal;

    m_modelgce = NULL;
    m_modelgce_default = new ThermalModelIdeal(m_TPS, m_Parameters);
    m_modelgce_ext = NULL;

    m_Banalyt = false;
    m_PartialZCalculated = false;
    m_PartialZCalculatedFlucts = false;

    m_SaddlePointDim = 0;
    m_SaddlePointLogDetSigma = 0.0;
    m_NLOEnergyDensity = 0.0;
    m_NLOPressure = 0.0;
    m_NLOEntropyDensity = 0.0;
  }


  
  ThermalModelCanonical::~ThermalModelCanonical(void)
  {
    delete m_modelgce_default;
  }

  void ThermalModelCanonical::ChangeTPS(ThermalParticleSystem *TPS_) {
    ThermalModelBase::ChangeTPS(TPS_);
    delete m_modelgce_default;
    m_modelgce_default = new ThermalModelIdeal(m_TPS, m_Parameters);
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

    m_PartialZCalculatedFlucts = computeFluctuations;
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

  // Solve non-canonical chemical potentials using a temporary GCE model.
  // This avoids running the expensive (and potentially unstable at large muB)
  // canonical partition function integration during each Broyden iteration.
  void ThermalModelCanonical::SolveChemicalPotentialsGCE(bool resetInitialValues)
  {
    // Build a temporary GCE model and solve for the non-canonical chemical potentials
    ThermalModelIdeal tempGCE(m_TPS, m_Parameters);
    tempGCE.SetUseWidth(m_UseWidth);

    // Only constrain charges that are NOT canonical (canonical mus will be set to 0 later)
    bool cB = m_ConstrainMuB && !m_BCE;
    bool cQ = m_ConstrainMuQ && !m_QCE;
    bool cS = m_ConstrainMuS && !m_SCE;
    bool cC = m_ConstrainMuC && !m_CCE;

    if (!cB && !cQ && !cS && !cC) {
      // Nothing to solve for
      return;
    }

    if (resetInitialValues) {
      // Use reasonable initial guesses
      if (cS)
        tempGCE.SetStrangenessChemicalPotential(m_Parameters.muB / 3.);
      if (cQ)
        tempGCE.SetElectricChemicalPotential(-m_Parameters.muB / 30.);
    }

    // Copy constraint settings
    tempGCE.ConstrainMuB(cB);
    tempGCE.ConstrainMuQ(cQ);
    tempGCE.ConstrainMuS(cS);
    tempGCE.ConstrainMuC(cC);
    tempGCE.SetQoverB(QoverB());
    tempGCE.SetSoverB(SoverB());

    if (resetInitialValues)
      tempGCE.FixParameters();
    else
      tempGCE.FixParametersNoReset();

    // Copy solved chemical potentials back
    m_Parameters.muB = tempGCE.Parameters().muB;
    m_Parameters.muQ = tempGCE.Parameters().muQ;
    m_Parameters.muS = tempGCE.Parameters().muS;
    m_Parameters.muC = tempGCE.Parameters().muC;
  }

  void ThermalModelCanonical::FixParameters()
  {
    // Initial guess for gammaC if constrained
    if (m_ConstrainGammaC) {
      m_Parameters.gammaC = 1.0;
    }

    m_ConstrainMuB &= !m_BCE;
    m_ConstrainMuQ &= !m_QCE;
    m_ConstrainMuS &= !m_SCE;
    m_ConstrainMuC &= !m_CCE;

    // Step 1: GCE solve for stable initial chemical potentials.
    // This avoids the CE Broyden solver wandering into extreme mu values
    // where the GL integration may break down.
    SolveChemicalPotentialsGCE(true);

    // Step 2: CE refinement using the base class Broyden solver,
    // starting from the GCE-solved values.
    // Each Broyden step calls CalculatePrimordialDensities(),
    // which accounts for canonical corrections in the constraint equations.
    // Since we start near the solution, the Broyden converges quickly.
    // Note: For GL at finite muB, the integration may still be inaccurate.
    //       Use the SaddlePoint method for best results at finite muB.
    ThermalModelBase::FixParametersNoReset();
  }

  void ThermalModelCanonical::FixParametersNoReset()
  {
    m_ConstrainMuB &= !m_BCE;
    m_ConstrainMuQ &= !m_QCE;
    m_ConstrainMuS &= !m_SCE;
    m_ConstrainMuC &= !m_CCE;

    // Step 1: GCE solve for stable initial chemical potentials.
    SolveChemicalPotentialsGCE(false);

    // Step 2: CE refinement using the base class Broyden solver.
    ThermalModelBase::FixParametersNoReset();
  }


  void ThermalModelCanonical::CalculatePrimordialDensities() {
    assert(m_IGFExtraConfig.MagneticField.B == 0.); // No magnetic field supported currently

    m_FluctuationsCalculated = false;

    if (m_PartialZ.size() == 0)
      CalculateQuantumNumbersRange();

    // Non-ideal interaction models are only supported with SaddlePointLO; any
    // other method (incl. GaussLegendre and plain SaddlePoint) falls back to it.
    EnforceSaddlePointLOForNonIdeal();

    // SaddlePointNLO/LO: compute means directly from CGF expansion, bypass partition functions
    if (m_Method == SaddlePointNLO || m_Method == SaddlePointLO) {
      m_Banalyt = false;
      ComputeAnalyticCumulants(false);
      m_Calculated = true;
      ValidateCalculation();
      return;
    }

    if (m_Method == SaddlePoint) {
      // Saddle-point method: chemical potentials will be handled inside
      // ComputeSaddlePointPartitionFunctions() — zeroed after saddle point is found
      m_Banalyt = false;
    }
    else if (m_BMAX_list == 1 && m_BCE && m_QCE && m_SCE && m_CCE && !UsePartialChemicalEquilibrium()) {
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
        else {
          cout << "ThermalModelCanonical::CalculatePrimordialDensities: Warning! No canonical partition function for this particle!" << endl;
          cout << "B = " << m_BCE * tpart.BaryonCharge() << "\tQ = " << m_QCE * tpart.ElectricCharge() << "\tS = " << m_SCE * tpart.Strangeness() << "\tC = " << m_CCE * tpart.Charm() << endl;
        }
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
    // This method is public, so a direct external call must not bypass the
    // non-ideal -> SaddlePointLO fallback that guards CalculateDensities() and
    // CalculateFluctuations().  SaddlePointLO/NLO do not use precomputed
    // partition functions (means and fluctuations come from
    // ComputeAnalyticCumulants()), so there is nothing to compute here for them.
    EnforceSaddlePointLOForNonIdeal();
    if (m_Method == SaddlePointLO || m_Method == SaddlePointNLO)
      return;

    if (Vc < 0.0)
      Vc = m_Parameters.SVc;

    if (!UsePartialChemicalEquilibrium())
      FillChemicalPotentials();
    else {
      // Partial chemical equilibrium canonical ensemble currently works only if there is particle-antiparticle symmetry
      for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
        int i2 = m_TPS->PdgToId(-m_TPS->Particle(i).PdgId());
        if (i2 != -1) {
          if (std::abs(m_Chem[i] - m_Chem[i2]) > 1.e-8) {
            throw std::runtime_error("ThermalModelCanonical::CalculatePartitionFunctions: Partial chemical equilibrium canonical ensemble only supported if particle-antiparticle fugacities are symmetric!");
          }
        }
      }
    }

    // Saddle-point: self-contained, no need for Fourier coefficients
    if (m_Method == SaddlePoint) {
      ComputeSaddlePointPartitionFunctions();

      m_Corr.resize(m_PartialZ.size());
      for (size_t iN = 0; iN < m_PartialZ.size(); ++iN) {
        m_Corr[iN] = m_PartialZ[iN] / m_PartialZ[m_QNMap[QuantumNumbers(0, 0, 0, 0)]];
      }
      m_PartialZCalculated = true;
      return;
    }

    // Contour shift: set canonical chemical potentials to their
    // GCE saddle-point values for dramatically improved GL convergence.
    // This is mathematically exact (equivalent to shifting the Fourier
    // contour through the saddle point). The fugacity factors from the
    // shifted densities and the Z-tilde ratios cancel, giving the
    // correct canonical result for all observables.
    //
    // Not used with m_Banalyt (analytical baryon fugacity via Bessel
    // functions), which assumes N_+ = N_- symmetry that breaks at mu != 0.
    // m_Banalyt is only active when ALL four charges are canonical, where
    // the GL integrand is already well-behaved at mu = 0.
    if (!m_Banalyt) {
      PrepareModelGCE();
      CleanModelGCE();
      if (!UsePartialChemicalEquilibrium())
        FillChemicalPotentials();
    }

    // --- Gauss-Legendre quadrature ---

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
          throw std::invalid_argument("ThermalModelCanonical: neutral particle cannot have non-zero ce charges");
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
          int ind = m_QNMap[QuantumNumbers(m_BCE * n * tpart.BaryonCharge(), m_QCE * n * tpart.ElectricCharge(), m_SCE * n * tpart.Strangeness(), m_CCE * n * tpart.Charm())];
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

    {
    // --- Gauss-Legendre quadrature ---

    int nmax = max(3, (int)sqrt((double)(m_Parameters.B*m_Parameters.B + m_Parameters.Q*m_Parameters.Q + m_Parameters.S*m_Parameters.S + m_Parameters.C*m_Parameters.C)));
    if (m_Parameters.B == 0 && m_Parameters.Q == 0 && m_Parameters.S == 0 && m_Parameters.C == 0)
      nmax = 4;

    nmax *= m_IntegrationIterationsMultiplier;


    m_MultExp = 0.;
    m_MultExpBanalyt = 0.;
    for (size_t i = 0; i < m_PartialZ.size(); ++i) {
      if (!m_Banalyt || m_QNvec[i].B == 0)
        m_MultExp += Nsx[i];
      if (m_Banalyt && (m_QNvec[i].B == 1 || m_QNvec[i].B == -1))
        m_MultExpBanalyt += Nsx[i];
    }

    double dphiB = xMath::Pi() / nmax;
    int maxB = 2 * nmax;
    if (m_BMAX == 0 || m_Banalyt)
      maxB = 1;

    for (int iB = 0; iB < maxB; ++iB) {

      vector<double> xlegB, wlegB;

      if (m_BMAX != 0 && !m_Banalyt) {
        double aB = iB * dphiB;
        if (iB >= nmax) aB = xMath::Pi() - (iB + 1) * dphiB;
        double bB = aB + dphiB;
        NumericalIntegration::GetCoefsIntegrateLegendre10(aB, bB, &xlegB, &wlegB);
      }
      else {
        xlegB.resize(1);
        xlegB[0] = 0.;
        wlegB.resize(1);
        wlegB[0] = 1.;
      }


      double dphiS = xMath::Pi() / nmax;
      int maxS = 2 * nmax;
      if (m_SMAX == 0)
        maxS = 1;

      for (int iS = 0; iS < maxS; ++iS) {
        vector<double> xlegS, wlegS;

        if (m_SMAX != 0) {
          double aS = iS * dphiS;
          if (iS >= nmax) aS = xMath::Pi() - (iS + 1) * dphiS;
          double bS = aS + dphiS;
          NumericalIntegration::GetCoefsIntegrateLegendre10(aS, bS, &xlegS, &wlegS);
        }
        else {
          xlegS.resize(1);
          xlegS[0] = 0.;
          wlegS.resize(1);
          wlegS[0] = 1.;
        }

        double dphiQ = xMath::Pi() / nmax;
        int maxQ = nmax;
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

          double dphiC = xMath::Pi() / nmax;
          int maxC = 2 * nmax;
          if (m_CMAX == 0)
            maxC = 1;

          for (int iC = 0; iC < maxC; ++iC) {
            vector<double> xlegC, wlegC;
            if (m_CMAX != 0) {
              double aC = iC * dphiC;
              if (iC >= nmax) aC = xMath::Pi() - (iC + 1) * dphiC;
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
      if (m_BMAX != 0 && !m_Banalyt)
        m_PartialZ[iN] /= 2. * xMath::Pi();
      if (m_QMAX != 0) {
        m_PartialZ[iN] /= 2. * xMath::Pi();
        m_PartialZ[iN] *= 2.; // Factor 2 from half-range [0,pi] integration for phi_Q
      }
      if (m_SMAX != 0)
        m_PartialZ[iN] /= 2. * xMath::Pi();
      if (m_CMAX != 0)
        m_PartialZ[iN] /= 2. * xMath::Pi();
    }

    } // end GL quadrature block

    m_Corr.resize(m_PartialZ.size());
    for (size_t iN = 0; iN < m_PartialZ.size(); ++iN) {
      m_Corr[iN] = m_PartialZ[iN] / m_PartialZ[m_QNMap[QuantumNumbers(0, 0, 0, 0)]];
    }

    m_PartialZCalculated = true;
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
      int ind = m_QNMap[QuantumNumbers(m_BCE * tpart.BaryonCharge(), m_QCE * tpart.ElectricCharge(), m_SCE * tpart.Strangeness(), m_CCE * tpart.Charm())];
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
          if (m_QNMap.count(QuantumNumbers(m_BCE * (n + n2)*tpart.BaryonCharge(), m_QCE * (n + n2)*tpart.ElectricCharge(), m_SCE * (n + n2)*tpart.Strangeness(), m_CCE * (n + n2)*tpart.Charm())) != 0) {
            int ind2 = m_QNMap[QuantumNumbers(m_BCE * (n + n2)*tpart.BaryonCharge(), m_QCE * (n + n2)*tpart.ElectricCharge(), m_SCE * (n + n2)*tpart.Strangeness(), m_CCE * (n + n2)*tpart.Charm())];
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
          int ind = m_QNMap[QuantumNumbers(m_BCE * n * tpart.BaryonCharge(), m_QCE * n * tpart.ElectricCharge(), m_SCE * n * tpart.Strangeness(), m_CCE * n * tpart.Charm())];

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
      if (m_wprim[i] != m_wprim[i]) m_wprim[i] = 1.;
    }

    m_TwoParticleCorrelationsCalculated = true;
  }

  void ThermalModelCanonical::CalculateFluctuations() {
    EnforceSaddlePointLOForNonIdeal();
    if (m_Method == SaddlePointNLO || m_Method == SaddlePointLO) {
      // SaddlePointNLO/LO: compute means + covariance in one shot
      if (m_PartialZ.size() == 0)
        CalculateQuantumNumbersRange();
      ComputeAnalyticCumulants(true);
      m_Calculated = true;
      // Skip CalculateTwoParticleCorrelations — already filled by ComputeAnalyticCumulants
    }
    else {
      if (!m_PartialZCalculatedFlucts) {
        CalculateQuantumNumbersRange(true);
        CalculatePrimordialDensities();
      }
      CalculateTwoParticleCorrelations();
    }

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

    // SaddlePointNLO/LO: energy density from saddle-point expansion.
    //   NLO: ε = ε^GCE(μ*) - (1/2) Σ_j E_j G_jj
    //   LO:  ε = ε^GCE(μ*)
    // Computed and stored by ComputeAnalyticCumulants().
    if (m_Method == SaddlePointNLO || m_Method == SaddlePointLO)
      return m_NLOEnergyDensity;

    double ret = 0.;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ThermalParticle &tpart = m_TPS->Particle(i);
      {
        if (!IsParticleCanonical(tpart)) {
          ret += tpart.Density(m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);
        }
        else if (tpart.Statistics() == 0
          || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion) {
          int ind = m_QNMap[QuantumNumbers(m_BCE * tpart.BaryonCharge(), m_QCE * tpart.ElectricCharge(), m_SCE * tpart.Strangeness(), m_CCE * tpart.Charm())];

          if (ind < static_cast<int>(m_Corr.size()))
            ret += m_Corr[ind] * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_Chem[i]);
        }
        else {
          for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
            int ind = m_QNMap[QuantumNumbers(m_BCE * n * tpart.BaryonCharge(), m_QCE * n * tpart.ElectricCharge(), m_SCE * n * tpart.Strangeness(), m_CCE * n * tpart.Charm())];
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

    // SaddlePointNLO/LO: pressure from saddle-point expansion.
    //   NLO: P = p^GCE(μ*) - d·T/(2·Vc)
    //   LO:  P = p^GCE(μ*)
    // Computed and stored by ComputeAnalyticCumulants().
    if (m_Method == SaddlePointNLO || m_Method == SaddlePointLO)
      return m_NLOPressure;

    double ret = 0.;

    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      ThermalParticle &tpart = m_TPS->Particle(i);
      {
        if (!IsParticleCanonical(tpart)) {
          ret += tpart.Density(m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i]);
        }
        else if (tpart.Statistics() == 0
          || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion) {
          int ind = m_QNMap[QuantumNumbers(m_BCE * tpart.BaryonCharge(), m_QCE * tpart.ElectricCharge(), m_SCE * tpart.Strangeness(), m_CCE * tpart.Charm())];

          if (ind < static_cast<int>(m_Corr.size()))
            ret += m_Corr[ind] * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, m_Chem[i]);
        }
        else {
          for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
            int ind = m_QNMap[QuantumNumbers(m_BCE * n * tpart.BaryonCharge(), m_QCE * n * tpart.ElectricCharge(), m_SCE * n * tpart.Strangeness(), m_CCE * n * tpart.Charm())];

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
    if (!m_Calculated) CalculateDensities();

    // SaddlePointLO: entropy is the GCE entropy at the saddle point μ*,
    // consistent with the LO pressure and energy density (the ln det Sigma
    // Gaussian term is an NLO effect).  Stored by ComputeAnalyticCumulants().
    // NLO and the other methods use the ln Z_SP form below.
    if (m_Method == SaddlePointLO)
      return m_NLOEntropyDensity;

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

  // PrepareModelGCE() selects and configures the GCE model used to evaluate
  // densities, susceptibilities, and two-particle correlations at the
  // saddle-point chemical potentials.
  //
  // Uses the external model (set via SetModelGCE) if available,
  // otherwise uses the persistent default ThermalModelIdeal.
  // Both models share the same TPS, so statistics and width settings
  // are already in sync.  Only thermal parameters need to be updated.
  void ThermalModelCanonical::PrepareModelGCE()
  {
    m_modelgce = (m_modelgce_ext != NULL) ? m_modelgce_ext : m_modelgce_default;

    // Sync thermal parameters (T, V, mu, gammas, etc.)
    m_modelgce->SetParameters(m_Parameters);

    if (m_BCE)
      m_Parameters.muB = 0.0;

    if (!m_BCE && m_SCE)
      m_Parameters.muS = m_Parameters.muB / 3.;

    if (!m_BCE && m_QCE)
      m_Parameters.muQ = -m_Parameters.muB / 30.;

    // Only reset muC when charm is canonical; otherwise muC is a fixed
    // grand-canonical input and must be preserved through the GCE solve.
    if (m_CCE)
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
  }

  void ThermalModelCanonical::CleanModelGCE()
  {
    m_modelgce = NULL;
  }

  void ThermalModelCanonical::SetModelGCE(ThermalModelBase* model)
  {
    m_modelgce_ext = model;
    if (model != NULL)
      m_InteractionModel = model->InteractionModel();
    else
      m_InteractionModel = Ideal;
    m_PartialZCalculated = false;
  }

  void ThermalModelCanonical::EnforceSaddlePointLOForNonIdeal()
  {
    // Non-ideal equations of state are only supported by the SaddlePointLO
    // method.  Every other method (GaussLegendre, SaddlePoint, SaddlePointNLO)
    // assumes an ideal gas, so silently mixing a non-ideal GCE model with them
    // would be inconsistent (ratios/thermodynamics from the interacting model,
    // densities from ideal-gas cluster sums).  Fall back to SaddlePointLO.
    if (m_InteractionModel != Ideal && m_Method != SaddlePointLO) {
      std::cerr << "**WARNING** " << m_TAG
                << ": Non-ideal interaction model is only supported with the SaddlePointLO method. "
                << "Falling back to SaddlePointLO." << std::endl;
      m_Method = SaddlePointLO;
    }
  }

  // --- Saddle-point approximation methods ---

  void ThermalModelCanonical::SolveSaddlePointEquations()
  {
    // Use PrepareModelGCE to solve the GCE constraint equations:
    //   <Q_a>_GCE(mu*/T) = N_a / V_c
    // This IS the saddle-point condition.
    PrepareModelGCE();

    // Store saddle-point chemical potentials mu*_a/T
    m_SaddlePointMu.resize(4);
    m_SaddlePointMu[0] = m_Parameters.muB / m_Parameters.T;
    m_SaddlePointMu[1] = m_Parameters.muQ / m_Parameters.T;
    m_SaddlePointMu[2] = m_Parameters.muS / m_Parameters.T;
    m_SaddlePointMu[3] = m_Parameters.muC / m_Parameters.T;

    // Store per-particle chemical potentials at the saddle point
    int N = m_TPS->ComponentsNumber();
    m_SaddlePointMuStar.resize(N);
    for (int j = 0; j < N; ++j) {
      m_SaddlePointMuStar[j] = m_modelgce->ChemicalPotential(j);
    }

    // NOTE: m_modelgce is intentionally NOT cleaned up here.
    // The caller manages cleanup so that ComputeAnalyticCumulants can use
    // m_modelgce to compute GCE two-particle correlations (cross-susceptibilities
    // for non-ideal models).

    // Zero the canonical chemical potentials
    if (m_BCE)
      m_Parameters.muB = 0.0;
    if (m_QCE)
      m_Parameters.muQ = 0.0;
    if (m_SCE)
      m_Parameters.muS = 0.0;
    if (m_CCE)
      m_Parameters.muC = 0.0;

    FillChemicalPotentials();
  }

  void ThermalModelCanonical::SetupSaddlePointChargeIndices()
  {
    m_SaddlePointChargeIndices.clear();
    if (m_BCE && m_BMAX_list > 0) m_SaddlePointChargeIndices.push_back(0);  // B
    if (m_QCE && m_QMAX_list > 0) m_SaddlePointChargeIndices.push_back(1);  // Q
    if (m_SCE && m_SMAX_list > 0) m_SaddlePointChargeIndices.push_back(2);  // S
    if (m_CCE && m_CMAX_list > 0) m_SaddlePointChargeIndices.push_back(3);  // C
    m_SaddlePointDim = static_cast<int>(m_SaddlePointChargeIndices.size());
  }

  double ThermalModelCanonical::ComputeSigmaLogDetFromModel()
  {
    int d = m_SaddlePointDim;
    if (d == 0)
      return 0.0;

    int Nparts = m_TPS->ComponentsNumber();
    double Vc = m_Parameters.SVc;
    double T  = m_Parameters.T;
    const auto& pc = m_modelgce->PrimCorrel();

    // Build Sigma_{ab} = sum_{j,k} q_{a,j} q_{b,k} chi^GCE_{jk}
    // where chi^GCE_{jk} = PrimCorrel[j][k] * Vc * T
    Eigen::MatrixXd Sigma = Eigen::MatrixXd::Zero(d, d);
    for (int j = 0; j < Nparts; ++j) {
      for (int k = 0; k < Nparts; ++k) {
        if (pc[j][k] == 0.0) continue;
        double chi_jk = pc[j][k] * Vc * T;
        for (int ia = 0; ia < d; ++ia) {
          int a = m_SaddlePointChargeIndices[ia];
          double q_aj = m_TPS->Particle(j).GetCharge(a);
          if (q_aj == 0.0) continue;
          for (int ib = 0; ib < d; ++ib) {
            int b = m_SaddlePointChargeIndices[ib];
            Sigma(ia, ib) += q_aj * m_TPS->Particle(k).GetCharge(b) * chi_jk;
          }
        }
      }
    }

    // log(det(Sigma)) via Cholesky
    double logDet;
    Eigen::LLT<Eigen::MatrixXd> llt(Sigma);
    if (llt.info() == Eigen::Success) {
      logDet = 0.0;
      for (int i = 0; i < d; ++i)
        logDet += 2.0 * log(llt.matrixL()(i, i));
    } else {
      double det = Sigma.determinant();
      logDet = (det > 0.0) ? log(det) : -1.e100;
    }

    return logDet;
  }

  void ThermalModelCanonical::ComputeSaddlePointPartitionFunctions()
  {
    // Step 1: Find saddle-point chemical potentials.
    // m_modelgce is kept alive (no CleanModelGCE) so we can use it
    // for re-solving at shifted charges, pressure, and Sigma.
    SolveSaddlePointEquations();

    // Step 2: Build susceptibility matrix from m_modelgce
    SetupSaddlePointChargeIndices();
    m_modelgce->CalculateTwoParticleCorrelations();
    m_SaddlePointLogDetSigma = ComputeSigmaLogDetFromModel();

    int d = m_SaddlePointDim;
    m_MultExp = 0.0;
    m_MultExpBanalyt = 0.0;

    // Reference charges
    double Nref[4] = {
      static_cast<double>(m_Parameters.B),
      static_cast<double>(m_Parameters.Q),
      static_cast<double>(m_Parameters.S),
      static_cast<double>(m_Parameters.C)
    };

    int Nparts = m_TPS->ComponentsNumber();
    double Vc = m_Parameters.SVc;
    double T  = m_Parameters.T;

    // Reference pressure from the GCE model at the saddle point
    double p0 = m_modelgce->CalculatePressure();

    // Compute ln Z_SP(N) = Vc*p0/T - sum_a mu*_a N_a/T - d/2 ln(2pi) - 1/2 ln(det Sigma)
    {
      double term1 = p0 * Vc / T;
      double term2 = 0.0;
      for (int ia = 0; ia < d; ++ia) {
        int a = m_SaddlePointChargeIndices[ia];
        term2 -= m_SaddlePointMu[a] * Nref[a];
      }
      double term3 = -0.5 * d * log(2.0 * xMath::Pi()) - 0.5 * m_SaddlePointLogDetSigma;
      m_MultExp = term1 + term2 + term3;
    }

    double muBguess = m_SaddlePointMu[0] * T;
    double muQguess = m_SaddlePointMu[1] * T;
    double muSguess = m_SaddlePointMu[2] * T;
    double muCguess = m_SaddlePointMu[3] * T;

    // Compute reference open charm density at the saddle point (for analytical muC initial guess).
    // n_charm_ref = (1/2) * sum_{|C_j|>0} |C_j| * n_j(mu*)  ≈  half the absolute charm density
    double nCharmRef = 0.0;
    {
      const auto& dens = m_modelgce->Densities();
      for (int j = 0; j < Nparts; ++j) {
        int absC = abs(m_TPS->Particle(j).Charm());
        if (absC > 0)
          nCharmRef += absC * dens[j];
      }
    }
    nCharmRef *= 0.5;  // half = positive-charm density ≈ negative-charm density at mu*

    // Reference mu . N term
    double muDotN_0 = 0.0;
    for (int ia = 0; ia < d; ++ia) {
      int a = m_SaddlePointChargeIndices[ia];
      muDotN_0 += m_SaddlePointMu[a] * Nref[a];
    }

    double logDetSigma_0 = m_SaddlePointLogDetSigma;

    int zeroIdx = m_QNMap[QuantumNumbers(0, 0, 0, 0)];

    for (size_t iN = 0; iN < m_PartialZ.size(); ++iN) {
      if (static_cast<int>(iN) == zeroIdx) {
        m_PartialZ[iN] = 1.0;
        continue;
      }

      int charges[4] = { m_QNvec[iN].B, m_QNvec[iN].Q, m_QNvec[iN].S, m_QNvec[iN].C };

      // Shifted target charges
      double Nshifted[4] = { Nref[0], Nref[1], Nref[2], Nref[3] };
      for (int ia = 0; ia < d; ++ia) {
        int a = m_SaddlePointChargeIndices[ia];
        Nshifted[a] -= charges[a];
      }

      // Analytical muC estimate for the shifted charm target.
      // Avoids overshooting in the first Newton step when charm is rare.
      // <C> ≈ 2 * nCharmRef * V * sinh(muC/T)  →  muC = T * arcsinh(Cshifted / (2 * nCharmRef * V))
      double muCguess_shifted = muCguess;
      if (m_CCE && nCharmRef > 0.0 && Nshifted[3] != 0.0) {
        double arg = Nshifted[3] / (2.0 * nCharmRef * m_Parameters.V);
        muCguess_shifted = m_Parameters.T * asinh(arg);
      }

      // For zero-shifted charges (charges[a]==0), use the saddle-point mu as the initial guess
      // rather than the previous sector's solved mu.  This prevents garbage propagation from
      // sectors where the solver converged to an unphysical fixed point due to ill-conditioning.
      // Shifted charges use the running muGuess from nearby sectors for faster convergence.
      double muBinit = (charges[0] != 0) ? muBguess : m_SaddlePointMu[0] * m_Parameters.T;
      double muQinit = (charges[1] != 0) ? muQguess : m_SaddlePointMu[1] * m_Parameters.T;
      double muSinit = (charges[2] != 0) ? muSguess : m_SaddlePointMu[2] * m_Parameters.T;
      double muCinit = (charges[3] != 0) ? muCguess_shifted : m_SaddlePointMu[3] * m_Parameters.T;

      // Solve GCE constraint equations for shifted charges using m_modelgce.
      // This uses the correct equation of state (ideal or interacting).
      // All canonical charges are constrained to capture cross-correlations
      // (e.g. mu_C shifts when mu_B shifts, through charmed baryons like Lambda_c).
      bool converged = m_modelgce->SolveChemicalPotentials(
        Nshifted[0], Nshifted[1], Nshifted[2], Nshifted[3],
        muBinit, muQinit, muSinit, muCinit,
        static_cast<bool>(m_BCE), static_cast<bool>(m_QCE),
        static_cast<bool>(m_SCE), static_cast<bool>(m_CCE));

      if (!converged) {
        // Solver did not converge for this QN sector: treat as inaccessible
        printf("**WARNING** ThermalModelCanonical::ComputeSaddlePointPartitionFunctions: "
               "Solver did not converge for QN sector (B,Q,S,C) = (%.0f,%.0f,%.0f,%.0f). Skipping.\n",
               Nshifted[0], Nshifted[1], Nshifted[2], Nshifted[3]);
        m_PartialZ[iN] = 0.0;
        continue;
      }

      double muOverT_c[4];
      muOverT_c[0] = m_modelgce->Parameters().muB / T;
      muOverT_c[1] = m_modelgce->Parameters().muQ / T;
      muOverT_c[2] = m_modelgce->Parameters().muS / T;
      muOverT_c[3] = m_modelgce->Parameters().muC / T;

      // Check for NaN/inf in solved chemical potentials
      bool muValid = true;
      for (int ia = 0; ia < 4; ++ia) {
        if (!std::isfinite(muOverT_c[ia])) {
          muValid = false;
          break;
        }
      }
      if (!muValid) {
        printf("**WARNING** ThermalModelCanonical::ComputeSaddlePointPartitionFunctions: "
               "NaN/inf chemical potentials for QN sector (B,Q,S,C) = (%.0f,%.0f,%.0f,%.0f). Skipping.\n",
               Nshifted[0], Nshifted[1], Nshifted[2], Nshifted[3]);
        m_PartialZ[iN] = 0.0;
        continue;
      }

      // Update muGuess only for charges that are shifted in this sector.
      // For zero-shifted charges, the solved mu captures cross-correlations but
      // is not a useful initial guess for other sectors with different shift patterns.
      if (charges[0] != 0) muBguess = m_modelgce->Parameters().muB;
      if (charges[1] != 0) muQguess = m_modelgce->Parameters().muQ;
      if (charges[2] != 0) muSguess = m_modelgce->Parameters().muS;
      if (charges[3] != 0) muCguess = m_modelgce->Parameters().muC;

      std::vector<double> muStar_c(Nparts);
      for (int j = 0; j < Nparts; ++j) {
        muStar_c[j] = m_modelgce->ChemicalPotential(j);
      }

      // Term 1: Vc/T * [p_c - p_0]
      double p_c = m_modelgce->CalculatePressure();
      double deltaP = Vc * (p_c - p0) / T;

      // Term 2: -[mu_c . N' - mu_0 . N]/T
      double muDotN_c = 0.0;
      for (int ia = 0; ia < d; ++ia) {
        int a = m_SaddlePointChargeIndices[ia];
        muDotN_c += muOverT_c[a] * Nshifted[a];
      }
      double deltaMuN = muDotN_c - muDotN_0;

      // Term 3: -1/2 * [ln det Sigma_c - ln det Sigma_0]
      // Use the interacting model's PrimCorrel at the shifted mu
      m_modelgce->CalculateTwoParticleCorrelations();
      double logDetSigma_c = ComputeSigmaLogDetFromModel();
      double deltaLogDet = logDetSigma_c - logDetSigma_0;

      // ln R(c) = deltaP - deltaMuN - 1/2 * deltaLogDet
      double lnR = deltaP - deltaMuN - 0.5 * deltaLogDet;

      // Guard against NaN propagation
      if (!std::isfinite(lnR)) {
        m_PartialZ[iN] = 0.0;
        continue;
      }

      m_PartialZ[iN] = exp(lnR);
    }

    CleanModelGCE();
  }

  void ThermalModelCanonical::ComputeAnalyticCumulants(bool computeCovariance)
  {
    // Analytic LO (or LO+NLO) from the CGF expansion.
    // Single saddle-point solve, then direct computation of means and (optionally) covariances.
    // For SaddlePointLO, only the leading-order terms are kept.
    //
    // All GCE quantities (densities, thermodynamics, susceptibilities) are read
    // from m_modelgce, which is always valid (either the external interacting model
    // or the default ThermalModelIdeal).  This ensures correct results for both
    // ideal and non-ideal equations of state.

    bool isNLO = (m_Method == SaddlePointNLO);

    // Step 1: Solve saddle-point equations for mu*/T,
    //         then compute GCE densities and two-particle correlations at the saddle point.
    SolveSaddlePointEquations();

    m_modelgce->CalculateTwoParticleCorrelations();

    // Set up charge indices (B, Q, S, C) that are treated canonically.
    // We do NOT call BuildSusceptibilityMatrix() here because that builds
    // Sigma from ideal-gas chi2.  Instead, we build Sigma directly from
    // m_modelgce->PrimCorrel() below, which is correct for both ideal and
    // non-ideal equations of state.
    SetupSaddlePointChargeIndices();

    int Nparts = m_TPS->ComponentsNumber();
    double Vc = m_Parameters.SVc;
    double T  = m_Parameters.T;

    // Step 2: Read GCE densities, thermodynamic quantities,
    //         and susceptibilities from the GCE model.
    const auto& modelDensities = m_modelgce->Densities();
    const auto& pc = m_modelgce->PrimCorrel();

    double epsGCE = m_modelgce->CalculateEnergyDensity();
    double pGCE   = m_modelgce->CalculatePressure();

    // W1[j] = <N_j>_GCE = density_j * Vc  (canonical species only)
    // W2[j] = var_GCE(N_j) = chi^GCE_{jj} = PrimCorrel[j][j] * Vc * T
    std::vector<double> W1(Nparts, 0.0);
    std::vector<double> W2(Nparts, 0.0);
    for (int j = 0; j < Nparts; ++j) {
      if (!IsParticleCanonical(m_TPS->Particle(j)))
        continue;
      W1[j] = modelDensities[j] * Vc;
      W2[j] = pc[j][j] * Vc * T;
    }

    // Build the charge susceptibility matrix over a given set of charge axes
    // (indices into {B=0,Q=1,S=2,C=3}):
    //   Sigma_{ab} = sum_{j,k} q_{a,j} q_{b,k} chi^GCE_{jk}
    // For ideal gas (diagonal PrimCorrel), this reduces to sum_j q_aj q_bj W2[j].
    // For interacting models, off-diagonal chi^GCE_{jk} are included.
    auto buildSigma = [&](const std::vector<int>& idx) -> Eigen::MatrixXd {
      int dd = static_cast<int>(idx.size());
      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(dd, dd);
      for (int j = 0; j < Nparts; ++j) {
        for (int k = 0; k < Nparts; ++k) {
          if (pc[j][k] == 0.0) continue;
          double chi_jk = pc[j][k] * Vc * T;
          for (int ia = 0; ia < dd; ++ia) {
            double q_aj = m_TPS->Particle(j).GetCharge(idx[ia]);
            if (q_aj == 0.0) continue;
            for (int ib = 0; ib < dd; ++ib)
              S(ia, ib) += q_aj * m_TPS->Particle(k).GetCharge(idx[ib]) * chi_jk;
          }
        }
      }
      return S;
    };

    // Drop linearly dependent (redundant) charge axes so that Sigma is full
    // rank before it is inverted.  Example: a proton-only gas with canonical B
    // and Q has q_B = q_Q for every species, so Sigma is rank 1 and the second
    // axis must be removed.  Greedily keep an axis only if it strictly raises
    // the rank (keeps Sigma positive definite).  The surviving count is then
    // d = rank(Sigma), the effective dimensionality used in the NLO pressure.
    {
      std::vector<int> cand = m_SaddlePointChargeIndices;
      std::vector<int> kept;
      for (size_t c = 0; c < cand.size(); ++c) {
        std::vector<int> trial = kept;
        trial.push_back(cand[c]);
        Eigen::MatrixXd St = buildSigma(trial);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(St, Eigen::EigenvaluesOnly);
        double maxEv = es.eigenvalues()(trial.size() - 1);
        double minEv = es.eigenvalues()(0);
        if (maxEv > 0.0 && minEv > 1.e-10 * maxEv)
          kept = trial;
        else
          std::cerr << "**WARNING** " << m_TAG
                    << ": dropping linearly dependent conserved-charge axis "
                    << cand[c] << " (B=0,Q=1,S=2,C=3) from the saddle-point basis." << std::endl;
      }
      m_SaddlePointChargeIndices = kept;
      m_SaddlePointDim = static_cast<int>(kept.size());
    }

    int d = m_SaddlePointDim;
    Eigen::MatrixXd Sigma = buildSigma(m_SaddlePointChargeIndices);

    // Compute log(det(Sigma)) for entropy
    double logDetSigma = 0.0;
    if (d > 0) {
      Eigen::LLT<Eigen::MatrixXd> llt(Sigma);
      if (llt.info() == Eigen::Success) {
        for (int i = 0; i < d; ++i)
          logDetSigma += 2.0 * log(llt.matrixL()(i, i));
      } else {
        double det = Sigma.determinant();
        logDetSigma = (det > 0.0) ? log(det) : -1.e100;
      }
    }
    m_SaddlePointLogDetSigma = logDetSigma;

    // Compute m_MultExp = ln Z_SP for entropy calculation.
    // ln Z_SP = Vc * p_GCE / T - sum_a mu*_a N_a / T - d/2 ln(2pi) - 1/2 ln(det Sigma)
    {
      double term1 = pGCE * Vc / T;
      double term2 = 0.0;
      for (int ia = 0; ia < d; ++ia) {
        int a = m_SaddlePointChargeIndices[ia];
        double Na = (a == 0) ? m_Parameters.B : (a == 1) ? m_Parameters.Q : (a == 2) ? m_Parameters.S : m_Parameters.C;
        term2 -= m_SaddlePointMu[a] * Na;
      }
      double term3 = -0.5 * d * log(2.0 * xMath::Pi()) - 0.5 * logDetSigma;
      m_MultExp = term1 + term2 + term3;
    }
    m_MultExpBanalyt = 0.0;

    // Set reference partition function to 1 so that log(m_PartialZ[0]) = 0
    // in CalculateEntropyDensity(), which uses m_MultExp + log(m_PartialZ[0]).
    m_PartialZ[m_QNMap[QuantumNumbers(0, 0, 0, 0)]] = 1.0;

    // Step 2b (NLO only): Compute higher cluster moments W3, W4 and
    // energy-weighted cluster sums Ej from cluster integrals.
    // These are specific to the ideal-gas NLO expansion and cannot be
    // read from m_modelgce.
    std::vector<double> W3(Nparts, 0.0);
    std::vector<double> W4(Nparts, 0.0);
    std::vector<double> Ej(Nparts, 0.0);

    if (isNLO) {
      for (int j = 0; j < Nparts; ++j) {
        ThermalParticle& tpart = m_TPS->Particle(j);

        if (!IsParticleCanonical(tpart))
          continue;

        if (tpart.Statistics() == 0
          || tpart.CalculationType() != IdealGasFunctions::ClusterExpansion)
        {
          // Boltzmann: all W^(k) equal
          double omega = W2[j];
          W3[j] = omega;
          W4[j] = omega;
          Ej[j] = tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_SaddlePointMuStar[j]);
        }
        else {
          // Quantum statistics with cluster expansion
          for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
            double dens_n = tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, m_SaddlePointMuStar[j]);
            double omega_n = dens_n * Vc / static_cast<double>(n);
            W3[j] += n * n * n * omega_n;
            W4[j] += n * n * n * n * omega_n;

            double edens_n = tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, m_SaddlePointMuStar[j]);
            Ej[j] += n * n * edens_n;
          }
        }
      }
    }

    // Step 3: Precompute Sigma^{-1}, P, Gjj, computeG.
    //
    // P_{j,a} = sum_b Sigma^{-1}_{ab} q_{b,j}
    // G_{jk}  = q_j^T Sigma^{-1} q_k = sum_a q_{a,j} P_{k,a}
    //
    // These are used by:
    //   - Ideal-gas fast-path covariance: cov = delta_lm W2[l] - W2[l] W2[m] G_lm
    //   - NLO corrections to thermodynamics, means, and covariances
    // For non-ideal LO, they are not needed (the general covariance path
    // uses the full chi^GCE matrix and Sigma^{-1} directly).
    bool isIdealGas = (m_InteractionModel == Ideal);

    Eigen::MatrixXd SigmaInv = Sigma.inverse();

    std::vector<double> P(Nparts * d, 0.0);
    std::vector<double> Gjj(Nparts, 0.0);

    if (isIdealGas) {
      for (int j = 0; j < Nparts; ++j) {
        for (int ia = 0; ia < d; ++ia) {
          double val = 0.0;
          for (int ib = 0; ib < d; ++ib) {
            int b = m_SaddlePointChargeIndices[ib];
            val += SigmaInv(ia, ib) * m_TPS->Particle(j).GetCharge(b);
          }
          P[j * d + ia] = val;
        }
      }

      for (int j = 0; j < Nparts; ++j) {
        double g = 0.0;
        for (int ia = 0; ia < d; ++ia) {
          int a = m_SaddlePointChargeIndices[ia];
          g += m_TPS->Particle(j).GetCharge(a) * P[j * d + ia];
        }
        Gjj[j] = g;
      }
    }

    // Helper lambda: compute G_jk from precomputed P
    auto computeG = [&](int j, int k) -> double {
      double g = 0.0;
      for (int ia = 0; ia < d; ++ia) {
        int a = m_SaddlePointChargeIndices[ia];
        g += m_TPS->Particle(j).GetCharge(a) * P[k * d + ia];
      }
      return g;
    };

    // Step 3b: Compute thermodynamic quantities.
    //
    // NLO:
    //   Energy:  ε = ε^GCE(μ*) - (1/2) Σ_j E_j G_jj
    //   Pressure: P = p^GCE(μ*) - d·T/(2Vc)
    //   Entropy: s = ln Z_SP / Vc + ε/T  (computed in CalculateEntropyDensity)
    // LO:
    //   Energy:  ε = ε^GCE(μ*)
    //   Pressure: P = p^GCE(μ*)
    //   Entropy: s = s^GCE(μ*)  (pure GCE, consistent with the LO P and ε;
    //            the ln det Sigma Gaussian term is an NLO effect and is dropped)
    if (isNLO) {
      double deltaEps = 0.0;
      for (int j = 0; j < Nparts; ++j) {
        if (Ej[j] == 0.0 || Gjj[j] == 0.0)
          continue;
        deltaEps += Ej[j] * Gjj[j];
      }
      m_NLOEnergyDensity = epsGCE - 0.5 * deltaEps;
      m_NLOPressure = pGCE - d * T / (2.0 * Vc);
      m_NLOEntropyDensity = 0.; // unused for NLO; entropy comes from ln Z_SP path
    }
    else {
      m_NLOEnergyDensity = epsGCE;
      m_NLOPressure = pGCE;
      m_NLOEntropyDensity = m_modelgce->CalculateEntropyDensity(); // s^GCE(μ*)
    }

    // Step 4: Compute means (LO, or LO + NLO if isNLO)
    for (int l = 0; l < Nparts; ++l) {
      if (!IsParticleCanonical(m_TPS->Particle(l))) {
        m_densities[l] = modelDensities[l];
        continue;
      }

      // LO mean
      double mean = W1[l];

      if (isNLO) {
        // NLO correction: -1/2 * sum_j G_jj * W_j^(3) * D_jl
        double nloTerm1 = -0.5 * W3[l] * Gjj[l];

        double sumGjjW3Gjl = 0.0;
        for (int j = 0; j < Nparts; ++j) {
          if (W3[j] == 0.0 || Gjj[j] == 0.0)
            continue;
          double G_jl = computeG(j, l);
          sumGjjW3Gjl += W3[j] * Gjj[j] * G_jl;
        }
        double nloTerm2 = 0.5 * W2[l] * sumGjjW3Gjl;

        mean += nloTerm1 + nloTerm2;
      }

      m_densities[l] = mean / Vc;
    }

    if (!computeCovariance) {
      CleanModelGCE();
      return;
    }

    // Step 5: Compute LO + NLO covariance matrix
    int NN = Nparts;

    m_PrimCorrel.resize(NN);
    for (int i = 0; i < NN; ++i)
      m_PrimCorrel[i].resize(NN, 0.0);
    m_TotalCorrel = m_PrimCorrel;

    // For ideal gas, non-canonical species have no canonical correction
    // (chi^GCE is diagonal, so v_l = 0 when particle l carries no canonical charges).
    // Their GCE susceptibility is the final answer.
    // For non-ideal models, the full SAM formula in the general covariance path
    // below handles all species (including non-canonical) correctly.
    if (isIdealGas) {
      for (int l = 0; l < NN; ++l) {
        if (!IsParticleCanonical(m_TPS->Particle(l)))
          m_PrimCorrel[l][l] = pc[l][l];
      }
    }

    // =====================================================================
    // Step 5a: LO covariance
    //
    // General formula (SAM / Schur complement of Eq. (2.86)):
    //   cov(N_l, N_m) = chi^GCE_{lm} - v_l^T  Sigma_chi^{-1}  v_m
    //
    // For the ideal gas, chi_{jk} = delta_{jk} W2[j], which simplifies to:
    //   cov = delta_lm W2[l] - W2[l] W2[m] G_lm
    // where G_lm = q_l^T Sigma^{-1} q_m.
    // This avoids the O(N^2) chiGCE allocation and double sums.
    // =====================================================================

    if (isIdealGas) {
      // ---- Ideal gas fast path ----
      // cov(l,m) = delta_lm * W2[l] - W2[l] * W2[m] * G_lm
      // NLO terms also only use W2, computeG, P -- no chiGCE needed.
      for (int l = 0; l < NN; ++l) {
        if (!IsParticleCanonical(m_TPS->Particle(l)))
          continue;

        for (int m = l; m < NN; ++m) {
          if (!IsParticleCanonical(m_TPS->Particle(m)))
            continue;

          double cov = -W2[l] * W2[m] * computeG(l, m);
          if (l == m)
            cov += W2[l];

          if (isNLO) {
            // NLO: -1/2 B_lm + 1/2 C_lm  (identical to general formula)

            // --- C_lm ---
            std::vector<double> Ml(d * d, 0.0);
            std::vector<double> Mm(d * d, 0.0);

            for (int j = 0; j < Nparts; ++j) {
              if (W3[j] == 0.0)
                continue;

              double D_jl = (j == l ? 1.0 : 0.0) - W2[l] * computeG(j, l);
              double D_jm = (j == m ? 1.0 : 0.0) - W2[m] * computeG(j, m);

              double w3Djl = W3[j] * D_jl;
              double w3Djm = W3[j] * D_jm;

              for (int ic = 0; ic < d; ++ic) {
                double Pjc = P[j * d + ic];
                for (int ia = 0; ia < d; ++ia) {
                  int a = m_SaddlePointChargeIndices[ia];
                  double q_aj = m_TPS->Particle(j).GetCharge(a);
                  Ml[ic * d + ia] += w3Djl * Pjc * q_aj;
                  Mm[ic * d + ia] += w3Djm * Pjc * q_aj;
                }
              }
            }

            double C_lm = 0.0;
            for (int ic = 0; ic < d; ++ic)
              for (int ia = 0; ia < d; ++ia)
                C_lm += Ml[ic * d + ia] * Mm[ia * d + ic];

            // --- B_lm ---
            std::vector<double> S_blm(d, 0.0);
            for (int k = 0; k < Nparts; ++k) {
              if (W3[k] == 0.0)
                continue;
              double D_kl = (k == l ? 1.0 : 0.0) - W2[l] * computeG(k, l);
              double D_km = (k == m ? 1.0 : 0.0) - W2[m] * computeG(k, m);
              double w3DklDkm = W3[k] * D_kl * D_km;
              for (int ib = 0; ib < d; ++ib) {
                int b = m_SaddlePointChargeIndices[ib];
                S_blm[ib] += m_TPS->Particle(k).GetCharge(b) * w3DklDkm;
              }
            }

            std::vector<double> SinvS(d, 0.0);
            for (int ia = 0; ia < d; ++ia) {
              for (int ib = 0; ib < d; ++ib)
                SinvS[ia] += SigmaInv(ia, ib) * S_blm[ib];
            }

            double B_lm = 0.0;
            for (int j = 0; j < Nparts; ++j) {
              if (Gjj[j] == 0.0)
                continue;

              double D_jl = (j == l ? 1.0 : 0.0) - W2[l] * computeG(j, l);
              double D_jm = (j == m ? 1.0 : 0.0) - W2[m] * computeG(j, m);

              double term1 = W4[j] * D_jl * D_jm;

              double d2alpha_j = 0.0;
              for (int ia = 0; ia < d; ++ia) {
                int a = m_SaddlePointChargeIndices[ia];
                d2alpha_j -= m_TPS->Particle(j).GetCharge(a) * SinvS[ia];
              }
              double term2 = W3[j] * d2alpha_j;

              B_lm += Gjj[j] * (term1 + term2);
            }

            cov += -0.5 * B_lm + 0.5 * C_lm;
          }

          double susc = cov / (Vc * T);
          m_PrimCorrel[l][m] = susc;
          if (l != m)
            m_PrimCorrel[m][l] = susc;
        }
      }
    } else {
      // ---- Non-ideal (general) path ----
      // Use the full N x N GCE susceptibility matrix from m_modelgce.
      //   cov(l,m) = chi^GCE_{lm} - v_l^T Sigma^{-1} v_m

      // Precompute v[l][a] = sum_j q_{a,j} * chi^GCE_{jl}
      // and u[l][a] = sum_b SigmaInv(a,b) * v[l][b]
      //
      // For non-canonical species l (all canonical charges zero), v_l is NOT
      // necessarily zero in interacting models: chi^GCE_{jl} can have off-diagonal
      // terms from interactions, so sum_j q_{a,j} chi^GCE_{jl} != 0 in general.
      // Therefore we compute v_l and u_l for ALL species, not just canonical ones.
      std::vector<std::vector<double>> vl(Nparts, std::vector<double>(d, 0.0));
      std::vector<std::vector<double>> ul(Nparts, std::vector<double>(d, 0.0));
      for (int l = 0; l < Nparts; ++l) {
        for (int j = 0; j < Nparts; ++j) {
          if (pc[j][l] == 0.0)
            continue;
          double chi_jl = pc[j][l] * Vc * T;
          for (int ia = 0; ia < d; ++ia) {
            int a = m_SaddlePointChargeIndices[ia];
            vl[l][ia] += m_TPS->Particle(j).GetCharge(a) * chi_jl;
          }
        }
        for (int ia = 0; ia < d; ++ia) {
          for (int ib = 0; ib < d; ++ib)
            ul[l][ia] += SigmaInv(ia, ib) * vl[l][ib];
        }
      }

      // LO covariance for ALL species pairs (canonical and non-canonical).
      // The SAM formula cov(l,m) = chi^GCE_{lm} - v_l^T Sigma^{-1} v_m
      // applies to all species: interactions can induce canonical corrections
      // even for non-canonical species through off-diagonal chi^GCE_{jl}.
      for (int l = 0; l < NN; ++l) {
        for (int m = l; m < NN; ++m) {
          double cov = pc[l][m] * Vc * T;
          for (int ia = 0; ia < d; ++ia)
            cov -= vl[l][ia] * ul[m][ia];

          double susc = cov / (Vc * T);
          m_PrimCorrel[l][m] = susc;
          if (l != m)
            m_PrimCorrel[m][l] = susc;
        }
      }
    }

    // Set scaled variances
    for (int i = 0; i < NN; ++i) {
      m_wprim[i] = m_PrimCorrel[i][i];
      if (m_densities[i] > 0.)
        m_wprim[i] *= T / m_densities[i];
      else
        m_wprim[i] = 1.;
      if (m_wprim[i] != m_wprim[i])
        m_wprim[i] = 1.;
    }

    m_TwoParticleCorrelationsCalculated = true;

    CleanModelGCE();
  }

} // namespace thermalfist
