#include "CosmicEos/EffectiveMassModel.h"

#include <iostream>

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

  EffectiveMassModel::EffectiveMassModel(const EffectiveMassModel& other)
  {
    *this = other;
    if (other.m_FieldPressure != NULL)
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

        if (meffinit != meffinit)
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
//      cout << m_TBEC << " " << m_Meff << " " << m_Mu << "\n";
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

    printf("**ERROR** EffectiveMassModel::Quantity(): Calculate %d quantity is not implemented!", static_cast<int>(quantity));
    exit(1);

    return 0.0;
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
    //cout << x[0] << " ";
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