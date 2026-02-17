#include "CosmicEos/EffectiveMassModel.h"

#include <iostream>
#include <vector>

using namespace std;

namespace thermalfist {
  EffectiveMassModel::EffectiveMassModel(const ThermalParticle& particle, EMMFieldPressure* FieldPressure) :
    GeneralizedDensity(),
    m_Particle(particle), m_FieldPressure(FieldPressure)
  {
    m_T  = -1.;
    m_Mu = 0.;
    m_isSolved = false;
    m_isBEC = false;
    m_TBEC = 0.050;
    m_Meff = m_Particle.Mass();
    m_Pressure = m_Density = m_DensityScalar = m_EntropyDensity = m_Chi2 = 0.;
  }

  EffectiveMassModel::~EffectiveMassModel(void)
  {
    if (m_FieldPressure != NULL)
      delete m_FieldPressure;
    m_FieldPressure = NULL;
  }

  EffectiveMassModel::EffectiveMassModel(const EffectiveMassModel& other) :
    GeneralizedDensity(),
    m_Particle(other.m_Particle), m_FieldPressure(NULL)
  {
    m_T = other.m_T;
    m_Mu = other.m_Mu;
    m_isSolved = other.m_isSolved;
    m_isBEC = other.m_isBEC;
    m_TBEC = other.m_TBEC;
    m_Meff = other.m_Meff;
    m_Pressure = other.m_Pressure;
    m_Density = other.m_Density;
    m_DensityScalar = other.m_DensityScalar;
    m_EntropyDensity = other.m_EntropyDensity;
    m_Chi2 = other.m_Chi2;
    
    // Use previous value of TBEC as starting guess
      m_FieldPressure = other.m_FieldPressure->clone();
  }

  double EffectiveMassModel::ComputeTBEC(double Mu, double TBECinit) const
  {
    // Assume effective mass always larger than vacuum mass (repulsive interaction only)
    if (m_Particle.Statistics() != -1 || Mu <= m_Particle.Mass())
      return 0.0;

    if ((Mu - m_Particle.Mass()) < 1.e-9)
      return 0.;

    // Use previous value of TBEC as strating guess
    if (TBECinit < 0.)
      TBECinit = m_TBEC;

    // Use small non-zero initial value of TBEC
    if (TBECinit == 0.0)
      TBECinit = 0.001;

    BroydenEquationsEMMTBEC eqs(this, Mu);
    Broyden broydn(&eqs);

    Broyden::BroydenSolutionCriterium criterium;

    vector<double> vars(1, log(TBECinit / m_Particle.Mass()));
    vars = broydn.Solve(vars, &criterium);

    double TBEC = m_Particle.Mass() * exp(vars[0]);

    // Fit did not converge, assume there is no condensation for given chemical potential
    if (broydn.Iterations() == broydn.MaxIterations() || vars[0] != vars[0] || std::isinf(vars[0]))
      TBEC = 0.;

    return TBEC;
  }

  void EffectiveMassModel::SetParameters(double T, double Mu)
  {
    if (m_T == 0.) {
      m_TBEC = 0.;
      m_isBEC = (Mu > m_Particle.Mass());
    }
    else if (m_Particle.Statistics() == -1 && (m_T < 0. || Mu != m_Mu)) {
      m_TBEC = ComputeTBEC(Mu);
      m_isBEC = (T <= m_TBEC);
    }

    if (T != m_T || Mu != m_Mu)
      m_isSolved = false;

    m_T = T;
    m_Mu = Mu;
  }

  void EffectiveMassModel::SolveMeff(double meffinit)
  {
    if (!m_isSolved) {
      if (m_T == 0.0) {
        if (!IsBECPhase())
          m_Meff = m_Particle.Mass();
        else
          m_Meff = m_Mu;
        m_isSolved = true;
        return;
      }
      if (!IsBECPhase()) {
        if (meffinit < 0.)
          meffinit = m_Meff;

        if (std::isnan(meffinit))
          meffinit = m_Particle.Mass();

        BroydenEquations* eqs;
        if (m_Particle.Statistics() == -1 && m_Mu > m_Particle.Mass())
          eqs = new BroydenEquationsEMMMeffConstrained(this);
        else
          eqs = new BroydenEquationsEMMMeff(this);
        Broyden broydn(eqs);

        Broyden::BroydenSolutionCriterium criterium;

        vector<double> vars(1, meffinit);
        vars = broydn.Solve(vars, &criterium);

        if (m_Particle.Statistics() == -1 && m_Mu > m_Particle.Mass())
          m_Meff = m_Mu * (1. + exp(vars[0]));
        else
          m_Meff = vars[0];

        delete eqs;
      }
      else {
        m_Meff = m_Mu;
      }
      m_isSolved = true;
    }
  }

  double EffectiveMassModel::Pressure() const
  {
    if (!IsSolved())
      return 0.0;

    if (m_T == 0.0) {
      if (IsBECPhase())
        return m_FieldPressure->pf(m_Meff);
      else
        return 0.;
    }

    ThermalModelParameters params;
    params.T = m_T;

    return IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::Pressure, IdealGasFunctions::Quadratures, m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy())
      + m_FieldPressure->pf(m_Meff);
  }

  double EffectiveMassModel::Density() const
  {
    if (!IsSolved())
      return 0.0;

    if (m_T == 0.0) {
      if (IsBECPhase())
        return m_FieldPressure->Dpf(m_Meff);
      else
        return 0.;
    }

    if (!IsBECPhase()) {
      return IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ParticleDensity, IdealGasFunctions::Quadratures, m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy());
    }
    else {
      return IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ParticleDensity, IdealGasFunctions::Quadratures, m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy())
        - IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures, m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy())
        + m_FieldPressure->Dpf(m_Meff);
    }
  }

  double EffectiveMassModel::EntropyDensity() const
  {
    if (!IsSolved())
      return 0.0;

    if (m_T == 0.0)
      return 0.;

    return IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::EntropyDensity, IdealGasFunctions::Quadratures, m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy());
  }

  double EffectiveMassModel::EnergyDensity() const
  {
    if (!IsSolved())
      return 0.0;

    return -Pressure() + m_T * EntropyDensity() + m_Mu * Density();
  }

  double EffectiveMassModel::Quantity(IdealGasFunctions::Quantity quantity, double T, double mu)
  {
    SetParameters(T, mu);
    SolveMeff();

    if (quantity == IdealGasFunctions::ParticleDensity)
      return Density();

    if (quantity == IdealGasFunctions::EnergyDensity)
      return EnergyDensity();

    if (quantity == IdealGasFunctions::EntropyDensity)
      return EntropyDensity();

    if (quantity == IdealGasFunctions::Pressure)
      return Pressure();

    // Chi2: evaluate using numerical derivative
    if (quantity == IdealGasFunctions::chi2difull) {
      double dmu = 0.001 * m_Particle.Mass();
      double n1  = Quantity(IdealGasFunctions::ParticleDensity, T, mu - 0.5 * dmu);
      double n2 =  Quantity(IdealGasFunctions::ParticleDensity, T, mu + 0.5 * dmu);
      SetParameters(T, mu);
      double ret = (n2 - n1) / dmu;
      ret /= xMath::GeVtoifm3();
      return ret;
    }

    // Chi3: evaluate using numerical derivative
    if (quantity == IdealGasFunctions::chi3difull) {
      double dmu = 0.001 * m_Particle.Mass();
      double n0  = Quantity(IdealGasFunctions::ParticleDensity, T, mu);
      double nm1 = Quantity(IdealGasFunctions::ParticleDensity, T, mu - dmu);
      double np1 = Quantity(IdealGasFunctions::ParticleDensity, T, mu + dmu);
      double ret = (np1 - 2.*n0 + nm1) / dmu / dmu;
      SetParameters(T, mu);
      ret /= xMath::GeVtoifm3();
      return ret;
    }

    // Chi4: evaluate using numerical derivative
    if (quantity == IdealGasFunctions::chi4difull) {
      double dmu = 0.001 * m_Particle.Mass();
      //double n0  = Quantity(IdealGasFunctions::ParticleDensity, T, mu);
      double nm2 = Quantity(IdealGasFunctions::ParticleDensity, T, mu - 2. * dmu);
      double nm1 = Quantity(IdealGasFunctions::ParticleDensity, T, mu - dmu);
      double np1 = Quantity(IdealGasFunctions::ParticleDensity, T, mu + dmu);
      double np2 = Quantity(IdealGasFunctions::ParticleDensity, T, mu + 2. * dmu);
      double ret = (0.5 * np2 - np1 + nm1 - 0.5 * nm2) / dmu / dmu / dmu;
      SetParameters(T, mu);
      ret /= xMath::GeVtoifm3();
      return ret;
    }

    if (quantity == IdealGasFunctions::chi2) {
      return Quantity(IdealGasFunctions::chi2difull, T, mu) / T / T;
    }

    if (quantity == IdealGasFunctions::chi3) {
      return Quantity(IdealGasFunctions::chi3difull, T, mu) / T;
    }

    if (quantity == IdealGasFunctions::chi4) {
      return Quantity(IdealGasFunctions::chi4difull, T, mu);
    }

    // Temperature derivatives: dQ/dT = (dQ^id/dT)|_{m*} + (dQ^id/dm*)|_T * dm*/dT
    // In BEC phase m* = mu (fixed), so dm*/dT = 0

    if (quantity == IdealGasFunctions::dndT) {
      if (m_T == 0.0)
        return 0.;

      if (!IsBECPhase()) {
        // Normal phase: dn/dT = (dn^id/dT)|_{m*} + (dn^id/dm*) * dm*/dT
        double dndT_fixedm = IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::dndT, IdealGasFunctions::Quadratures,
          m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy());
        double dndm = IdealGasQuantityMassDerivative(IdealGasFunctions::ParticleDensity);
        return dndT_fixedm + dndm * DmeffDT();
      }
      else {
        // BEC phase: n = n^id - rho_scalar^id + p_f'(m*), with m* = mu fixed
        // dn/dT = (dn^id/dT)|_{m*} - (drho_scalar^id/dT)|_{m*}
        double dndT_fixedm = IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::dndT, IdealGasFunctions::Quadratures,
          m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy());
        // d(rho_scalar)/dT at fixed m* via finite differences
        double dT = 1.e-4 * m_T;
        double rhos_p = IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures,
          m_Particle.Statistics(), m_T + 0.5 * dT, m_Mu, m_Meff, m_Particle.Degeneracy());
        double rhos_m = IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures,
          m_Particle.Statistics(), m_T - 0.5 * dT, m_Mu, m_Meff, m_Particle.Degeneracy());
        double drhosdT = (rhos_p - rhos_m) / dT;
        return dndT_fixedm - drhosdT;
      }
    }

    if (quantity == IdealGasFunctions::d2ndT2) {
      // Use numerical finite differences of dndT
      if (m_T == 0.0)
        return 0.;
      double dT = 1.e-4 * m_T;
      double dndT_p = Quantity(IdealGasFunctions::dndT, T + 0.5 * dT, mu);
      double dndT_m = Quantity(IdealGasFunctions::dndT, T - 0.5 * dT, mu);
      SetParameters(T, mu);
      SolveMeff();
      return (dndT_p - dndT_m) / dT;
    }

    if (quantity == IdealGasFunctions::dedT) {
      if (m_T == 0.0)
        return 0.;
      // de/dT = T * ds/dT + s + mu * dn/dT  (from e = -P + Ts + mu*n)
      // Or use numerical finite differences for simplicity and consistency
      double dT = 1.e-4 * m_T;
      double e_p = Quantity(IdealGasFunctions::EnergyDensity, T + 0.5 * dT, mu);
      double e_m = Quantity(IdealGasFunctions::EnergyDensity, T - 0.5 * dT, mu);
      SetParameters(T, mu);
      SolveMeff();
      return (e_p - e_m) / dT;
    }

    if (quantity == IdealGasFunctions::dedmu) {
      // de/dmu at fixed T: use finite differences
      double dmu = 0.001 * m_Particle.Mass();
      double e_p = Quantity(IdealGasFunctions::EnergyDensity, T, mu + 0.5 * dmu);
      double e_m = Quantity(IdealGasFunctions::EnergyDensity, T, mu - 0.5 * dmu);
      SetParameters(T, mu);
      SolveMeff();
      return (e_p - e_m) / dmu;
    }

    if (quantity == IdealGasFunctions::dchi2dT) {
      // dchi2/dT where chi2 = T^{-2} * dn/dmu (dimensionless)
      // Use numerical finite differences of chi2
      if (m_T == 0.0)
        return 0.;
      double dT = 1.e-4 * m_T;
      double chi2_p = Quantity(IdealGasFunctions::chi2, T + 0.5 * dT, mu);
      double chi2_m = Quantity(IdealGasFunctions::chi2, T - 0.5 * dT, mu);
      SetParameters(T, mu);
      SolveMeff();
      return (chi2_p - chi2_m) / dT;
    }

    if (quantity == IdealGasFunctions::dsdT) {
      if (m_T == 0.0)
        return 0.;

      // ds/dT = (ds^id/dT)|_{m*} + (ds^id/dm*)|_T * dm*/dT
      // The entropy is s = s^id(T, mu, m*) regardless of BEC phase
      // In BEC phase m* = mu (fixed), dm*/dT = 0
      double dsdT_fixedm = IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::dsdT, IdealGasFunctions::Quadratures,
        m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy());
      if (!IsBECPhase()) {
        double dsdm = IdealGasQuantityMassDerivative(IdealGasFunctions::EntropyDensity);
        return dsdT_fixedm + dsdm * DmeffDT();
      }
      else {
        return dsdT_fixedm;
      }
    }

    throw std::runtime_error("**ERROR** EffectiveMassModel::Quantity(): Calculate " + std::to_string(static_cast<int>(quantity)) + " quantity is not implemented!");

    return 0.0;
  }

  double EffectiveMassModel::DmeffDT() const
  {
    if (!IsSolved() || m_T == 0.0)
      return 0.;

    // In BEC phase, m* = mu which is independent of T
    if (IsBECPhase())
      return 0.;

    // Gap equation: p_f'(m*) = rho_scalar(T, mu, m*)
    // Implicit differentiation:
    // p_f''(m*) * dm*/dT = (d rho_scalar/dT)|_{m*} + (d rho_scalar/dm*)|_T * dm*/dT
    // => dm*/dT = (d rho_scalar/dT)|_{m*} / [p_f''(m*) - (d rho_scalar/dm*)|_T]

    // d rho_scalar / dT at fixed m*
    double dT = 1.e-4 * m_T;
    double rhos_p = IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures,
      m_Particle.Statistics(), m_T + 0.5 * dT, m_Mu, m_Meff, m_Particle.Degeneracy());
    double rhos_m = IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures,
      m_Particle.Statistics(), m_T - 0.5 * dT, m_Mu, m_Meff, m_Particle.Degeneracy());
    double drhosdT = (rhos_p - rhos_m) / dT;

    // d rho_scalar / dm* at fixed T
    double drhosdm = IdealGasQuantityMassDerivative(IdealGasFunctions::ScalarDensity);

    double D2pf = m_FieldPressure->D2pf(m_Meff);
    double denom = D2pf - drhosdm;

    // If denominator is too small, the gap equation is degenerate
    if (std::abs(denom) < 1.e-20)
      return 0.;

    return drhosdT / denom;
  }

  double EffectiveMassModel::IdealGasQuantityMassDerivative(IdealGasFunctions::Quantity quantity) const
  {
    if (!IsSolved() || m_T == 0.0)
      return 0.;

    double dm = 1.e-4 * m_Meff;
    if (dm < 1.e-8)
      dm = 1.e-8;

    double Q_p = IdealGasFunctions::IdealGasQuantity(quantity, IdealGasFunctions::Quadratures,
      m_Particle.Statistics(), m_T, m_Mu, m_Meff + 0.5 * dm, m_Particle.Degeneracy());
    double Q_m = IdealGasFunctions::IdealGasQuantity(quantity, IdealGasFunctions::Quadratures,
      m_Particle.Statistics(), m_T, m_Mu, m_Meff - 0.5 * dm, m_Particle.Degeneracy());

    return (Q_p - Q_m) / dm;
  }

  std::vector<double> EffectiveMassModel::BroydenEquationsEMMTBEC::Equations(const std::vector<double> &x)
  {
    double TBEC = m_EMM->m_Particle.Mass() * exp(x[0]);
    double ret = m_EMM->m_FieldPressure->Dpf(m_MuBEC) /
      IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures, m_EMM->m_Particle.Statistics(), TBEC, m_MuBEC, m_MuBEC, m_EMM->m_Particle.Degeneracy())
      - 1.;
    return std::vector<double>(1, ret);
  }

  std::vector<double> EffectiveMassModel::BroydenEquationsEMMMeff::Equations(const std::vector<double>& x)
  {
    double ret = m_EMM->m_FieldPressure->Dpf(x[0]) /
      IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures, m_EMM->m_Particle.Statistics(), m_EMM->m_T, m_EMM->m_Mu, x[0], m_EMM->m_Particle.Degeneracy())
      - 1.;
    return std::vector<double>(1, ret);
  }

  std::vector<double> EffectiveMassModel::BroydenEquationsEMMMeffConstrained::Equations(const std::vector<double>& x)
  {
    double meff = m_EMM->m_Mu * (1. + exp(x[0]));
    double ret = m_EMM->m_FieldPressure->Dpf(meff) /
      IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures, m_EMM->m_Particle.Statistics(), m_EMM->m_T, m_EMM->m_Mu, meff, m_EMM->m_Particle.Degeneracy())
      - 1.;
    return std::vector<double>(1, ret);
  }


  double EffectiveMassModel::BECFraction() const
  {
    if (!IsSolved() || !IsBECPhase())
      return 0.0;

    double n0 = m_FieldPressure->Dpf(m_Meff);
    if (m_T != 0.0)
      n0 -= IdealGasFunctions::IdealGasQuantity(IdealGasFunctions::ScalarDensity, IdealGasFunctions::Quadratures, m_Particle.Statistics(), m_T, m_Mu, m_Meff, m_Particle.Degeneracy());

    return n0 / Density();
  }

}