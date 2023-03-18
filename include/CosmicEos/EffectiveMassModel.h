#ifndef EFFECTIVEMASSMODEL_H
#define EFFECTIVEMASSMODEL_H
#include <map>

#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/Broyden.h"

namespace thermalfist {

    /** \file EffectiveMassModel.h
    *   \brief Header with effective mass model implementation, including the Bose-condensed phase
    *
    */

    /**
    *   \brief Base class implementing field pressure contribution function in the effective mass model.
    *          Default is linear scalar interaction, as in https://arxiv.org/abs/2004.09004 [Eq. (26)]
    */
    class EMMFieldPressure {
    protected:
      double m_mass;
      double m_c;
    public:
      EMMFieldPressure(double mass = 0.135, double c = 0) : m_mass(mass), m_c(c) {}

      virtual ~EMMFieldPressure() {}

      virtual EMMFieldPressure *clone() const { return new EMMFieldPressure(*this); }

      /// Field pressure as a function of effective mass x
      virtual double pf(double x) const { return (m_mass - x) * (m_mass - x) / 2. / m_c; }

      /// Derivative of the field pressure with respect to effective mass x
      virtual double Dpf(double x) const { return -(m_mass - x) / m_c; }

      /// Second derivative of the field pressure with respect to effective mass x
      virtual double D2pf(double x) const { return 1. / m_c; }
    };

    /**
    *   \brief Effective mass model matched to chiral perturbation theory for pions at T = 0.
    *          See supplemental material section I in https://arxiv.org/abs/2009.02309
    */
    class EMMFieldPressureChPT : public EMMFieldPressure {
    protected:
      double m_fpi;
    public:
      EMMFieldPressureChPT(double mass = 0.135, double fpi = 0.133) : EMMFieldPressure(mass), m_fpi(fpi) {}

      virtual ~EMMFieldPressureChPT() {}

      virtual EMMFieldPressure *clone() const { return new EMMFieldPressureChPT(*this); }

      virtual double pf(double x) const {
        return x * x * m_fpi * m_fpi / 4. * (1 - m_mass * m_mass / x / x) * (1 - m_mass * m_mass / x / x) *
               pow(xMath::GeVtoifm(), 3);
      }

      virtual double Dpf(double x) const {
        return x * m_fpi * m_fpi / 2. * (1 - pow(m_mass / x, 4)) * pow(xMath::GeVtoifm(), 3);
      }

      virtual double D2pf(double x) const {
        return m_fpi * m_fpi / 2. * (1. + 3. * pow(m_mass / x, 4)) * pow(xMath::GeVtoifm(), 3);
      }
    };

    /**
     * \brief Class implementing an effective mass model for single particle species.
     *
     */
    class EffectiveMassModel : public GeneralizedDensity {
      ThermalParticle m_Particle;
      EMMFieldPressure *m_FieldPressure;

      double m_T, m_Mu;
      bool m_isSolved;
      bool m_isBEC;

      double m_TBEC;

      double m_Meff;
      double m_Pressure;
      double m_Density;
      double m_DensityScalar;
      double m_EntropyDensity;
      double m_Chi2;
    public:
      /// Constructor.
      /// Default: pi+
      EffectiveMassModel(const ThermalParticle &particle = ThermalParticle(true, "pi", 211, 1.0, -1, 0.135, 0, 0, 1),
                         EMMFieldPressure *FieldPressure = NULL);

      virtual ~EffectiveMassModel(void);

      // copy constructor
      // TODO: proper way using swap function
      EffectiveMassModel(const EffectiveMassModel &other);

      /// Whether EMM equations are solved
      bool IsSolved() const { return m_isSolved; }

      /// Whether the particle is in a BEC phase
      virtual bool IsBECPhase() const { return m_isBEC; }

      /// Calculate temperature of BEC formation for given Mu
      double ComputeTBEC(double Mu, double TBECinit = -1.) const;

      /// Set the temperature and chemical potential of given particle
      void SetParameters(double T, double Mu);

      /// The present temperature
      double Temperature() const { return m_T; }

      /// The present chemical potential
      double ChemicalPotential() const { return m_Mu; }

      /// Find the effective mass for given parameters
      void SolveMeff(double meffinit = -1.);

      /// Set the effective mass manually
      void SetMeff(double meff) { m_Meff = meff; }

      /// The present effective mass
      double Meff() const { return m_Meff; }

      /// Pressure [GeV/fm^3]
      double Pressure() const;

      /// Number density [1/fm^3]
      double Density() const;

      /// Entropy density [1/fm^3]
      double EntropyDensity() const;

      /// Energy density [GeV/fm^3]
      double EnergyDensity() const;

      virtual double Quantity(IdealGasFunctions::Quantity quantity, double T, double mu);

      virtual double EffectiveMass() const { return Meff(); }

      virtual double BECFraction() const;

      /// Effective mass model equation to determine BEC onset at given chemical potential MuBEC
      /// To be solved using Broyden's method
      class BroydenEquationsEMMTBEC : public BroydenEquations {
      public:
        BroydenEquationsEMMTBEC(const EffectiveMassModel *model, double MuBEC)
                : BroydenEquations(), m_EMM(model), m_MuBEC(MuBEC) { m_N = 1; }

        std::vector<double> Equations(const std::vector<double> &x);

      private:
        const EffectiveMassModel *m_EMM;
        double m_MuBEC;
      };

      /// Effective mass model gap equation
      /// To be solved using Broyden's method
      class BroydenEquationsEMMMeff : public BroydenEquations {
      public:
        BroydenEquationsEMMMeff(const EffectiveMassModel *model) : BroydenEquations(), m_EMM(model) { m_N = 1; }

        std::vector<double> Equations(const std::vector<double> &x);

      private:
        const EffectiveMassModel *m_EMM;
      };

      /// Effective mass model gap equation
      /// Uses variable transformation to ensure m* >= mu making the method stable wrt BEC formation
      /// To be solved using Broyden's method
      class BroydenEquationsEMMMeffConstrained : public BroydenEquations {
      public:
        BroydenEquationsEMMMeffConstrained(const EffectiveMassModel *model)
                : BroydenEquations(), m_EMM(model) { m_N = 1; }

        std::vector<double> Equations(const std::vector<double> &x);

      private:
        const EffectiveMassModel *m_EMM;
      };

    };

} // namespace thermalfist

#endif

