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
      /**
       * \brief Constructor for the EMMFieldPressure class.
       * 
       * \param mass The mass parameter (default: 0.135 GeV).
       * \param c The coupling parameter (default: 0).
       */
      EMMFieldPressure(double mass = 0.135, double c = 0)
        : m_mass(mass),
          m_c(c) {}

      /**
       * \brief Virtual destructor for the EMMFieldPressure class.
       */
      virtual ~EMMFieldPressure() {}

      /**
       * \brief Creates a clone of this EMMFieldPressure object.
       * 
       * \return A pointer to a new EMMFieldPressure object.
       */
      virtual EMMFieldPressure *clone() const { return new EMMFieldPressure(*this); }

      /**
       * \brief Field pressure as a function of effective mass x.
       * 
       * \param x The effective mass value.
       * \return The field pressure value.
       */
      virtual double pf(double x) const { return (m_mass - x) * (m_mass - x) / 2. / m_c; }

      /**
       * \brief Derivative of the field pressure with respect to effective mass x.
       * 
       * \param x The effective mass value.
       * \return The derivative of the field pressure.
       */
      virtual double Dpf(double x) const { return -(m_mass - x) / m_c; }

      /**
       * \brief Second derivative of the field pressure with respect to effective mass x.
       * 
       * \param x The effective mass value.
       * \return The second derivative of the field pressure.
       */
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
      /**
       * \brief Constructor for the EMMFieldPressureChPT class.
       * 
       * \param mass The mass parameter (default: 0.135 GeV).
       * \param fpi The pion decay constant (default: 0.133 GeV).
       */
      EMMFieldPressureChPT(double mass = 0.135, double fpi = 0.133) : EMMFieldPressure(mass), m_fpi(fpi) {}

      /**
       * \brief Virtual destructor for the EMMFieldPressureChPT class.
       */
      virtual ~EMMFieldPressureChPT() {}

      /**
       * \brief Creates a clone of this EMMFieldPressureChPT object.
       * 
       * \return A pointer to a new EMMFieldPressureChPT object.
       */
      virtual EMMFieldPressure *clone() const { return new EMMFieldPressureChPT(*this); }

      /**
       * \brief Field pressure as a function of effective mass x according to ChPT.
       * 
       * \param x The effective mass value.
       * \return The field pressure value.
       */
      virtual double pf(double x) const {
        return x * x * m_fpi * m_fpi / 4. * (1 - m_mass * m_mass / x / x) * (1 - m_mass * m_mass / x / x) *
               pow(xMath::GeVtoifm(), 3);
      }

      /**
       * \brief Derivative of the field pressure with respect to effective mass x according to ChPT.
       * 
       * \param x The effective mass value.
       * \return The derivative of the field pressure.
       */
      virtual double Dpf(double x) const {
        return x * m_fpi * m_fpi / 2. * (1 - pow(m_mass / x, 4)) * pow(xMath::GeVtoifm(), 3);
      }

      /**
       * \brief Second derivative of the field pressure with respect to effective mass x according to ChPT.
       * 
       * \param x The effective mass value.
       * \return The second derivative of the field pressure.
       */
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
      /**
       * \brief Constructor for the EffectiveMassModel class.
       * 
       * \param particle The thermal particle to model (default: pi+).
       * \param FieldPressure Pointer to the field pressure model (default: NULL).
       */
      EffectiveMassModel(
          const ThermalParticle &particle = ThermalParticle(true, "pi", 211, 1.0, -1, 0.135, 0, 0, 1),
          EMMFieldPressure *FieldPressure = NULL);

      /**
       * \brief Destructor for the EffectiveMassModel class.
       */
      virtual ~EffectiveMassModel(void);

      /**
       * \brief Copy constructor for the EffectiveMassModel class.
       * 
       * \param other The EffectiveMassModel to copy.
       */
      EffectiveMassModel(const EffectiveMassModel &other);

      /**
       * \brief Checks if the EMM equations are solved.
       * 
       * \return True if the equations are solved, false otherwise.
       */
      bool IsSolved() const { return m_isSolved; }

      /**
       * \brief Checks if the particle is in a Bose-Einstein condensed phase.
       * 
       * \return True if the particle is in a BEC phase, false otherwise.
       */
      virtual bool IsBECPhase() const { return m_isBEC; }

      /**
       * \brief Calculates the temperature of BEC formation for a given chemical potential.
       * 
       * \param Mu The chemical potential value.
       * \param TBECinit Initial guess for the BEC temperature (default: -1).
       * \return The BEC formation temperature.
       */
      double ComputeTBEC(double Mu, double TBECinit = -1.) const;

      /**
       * \brief Sets the temperature and chemical potential of the particle.
       * 
       * \param T The temperature value.
       * \param Mu The chemical potential value.
       */
      void SetParameters(double T, double Mu);

      /**
       * \brief Gets the current temperature.
       * 
       * \return The temperature value.
       */
      double Temperature() const { return m_T; }

      /**
       * \brief Gets the current chemical potential.
       * 
       * \return The chemical potential value.
       */
      double ChemicalPotential() const { return m_Mu; }

      /**
       * \brief Solves for the effective mass using the given parameters.
       * 
       * \param meffinit Initial guess for the effective mass (default: -1).
       */
      void SolveMeff(double meffinit = -1.);

      /**
       * \brief Sets the effective mass manually.
       * 
       * \param meff The effective mass value to set.
       */
      void SetMeff(double meff) { m_Meff = meff; }

      /**
       * \brief Gets the current effective mass.
       * 
       * \return The effective mass value.
       */
      double Meff() const { return m_Meff; }

      /**
       * \brief Calculates the pressure in GeV/fm^3.
       * 
       * \return The pressure value.
       */
      double Pressure() const;

      /**
       * \brief Calculates the number density in 1/fm^3.
       * 
       * \return The number density value.
       */
      double Density() const;

      /**
       * \brief Calculates the entropy density in 1/fm^3.
       * 
       * \return The entropy density value.
       */
      double EntropyDensity() const;

      /**
       * \brief Calculates the energy density in GeV/fm^3.
       * 
       * \return The energy density value.
       */
      double EnergyDensity() const;

      /**
       * \brief Calculates a thermodynamic quantity.
       * 
       * \param quantity The type of quantity to calculate.
       * \param T The temperature value.
       * \param mu The chemical potential value.
       * \return The calculated quantity value.
       */
      virtual double Quantity(IdealGasFunctions::Quantity quantity, double T, double mu);

      /**
       * \brief Computes the temperature derivative of the effective mass dm_eff/dT at fixed mu.
       *
       * From implicit differentiation of the gap equation p_f'(m_eff) = rho_scalar(T, mu, m_eff):
       * dm_eff/dT = (d rho_scalar / dT) / [p_f''(m_eff) - (d rho_scalar / dm_eff)]
       *
       * \return dm_eff/dT in GeV/GeV (dimensionless).
       */
      double DmeffDT() const;

      /**
       * \brief Computes the derivative of an ideal gas quantity with respect to mass,
       *        (dQ^id/dm)|_{T,mu}, evaluated at the current effective mass m*.
       *        Uses central finite differences.
       *
       * \param quantity The ideal gas quantity.
       * \return The mass derivative.
       */
      double IdealGasQuantityMassDerivative(IdealGasFunctions::Quantity quantity) const;

      /**
       * \brief Gets the effective mass.
       *
       * \return The effective mass value.
       */
      virtual double EffectiveMass() const { return Meff(); }

      /**
       * \brief Calculates the fraction of particles in the BEC phase.
       * 
       * \return The BEC fraction value.
       */
      virtual double BECFraction() const;

      /**
       * \brief Class implementing the effective mass model equation to determine BEC onset.
       * To be solved using Broyden's method.
       */
      class BroydenEquationsEMMTBEC : public BroydenEquations {
      public:
        /**
         * \brief Constructor for the BroydenEquationsEMMTBEC class.
         * 
         * \param model Pointer to the EffectiveMassModel.
         * \param MuBEC The chemical potential at BEC onset.
         */
        BroydenEquationsEMMTBEC(const EffectiveMassModel *model, double MuBEC)
                : BroydenEquations(), m_EMM(model), m_MuBEC(MuBEC) { m_N = 1; }

        /**
         * \brief Implements the equations to be solved.
         * 
         * \param x The vector of variables.
         * \return The vector of equation values.
         */
        std::vector<double> Equations(const std::vector<double> &x);

      private:
        const EffectiveMassModel *m_EMM;
        double m_MuBEC;
      };

      /**
       * \brief Class implementing the effective mass model gap equation.
       * To be solved using Broyden's method.
       */
      class BroydenEquationsEMMMeff : public BroydenEquations {
      public:
        /**
         * \brief Constructor for the BroydenEquationsEMMMeff class.
         * 
         * \param model Pointer to the EffectiveMassModel.
         */
        BroydenEquationsEMMMeff(const EffectiveMassModel *model) : BroydenEquations(), m_EMM(model) { m_N = 1; }

        /**
         * \brief Implements the equations to be solved.
         * 
         * \param x The vector of variables.
         * \return The vector of equation values.
         */
        std::vector<double> Equations(const std::vector<double> &x);

      private:
        const EffectiveMassModel *m_EMM;
      };

      /**
       * \brief Class implementing the effective mass model gap equation with constraints.
       * Uses variable transformation to ensure m* >= mu making the method stable wrt BEC formation.
       * To be solved using Broyden's method.
       */
      class BroydenEquationsEMMMeffConstrained : public BroydenEquations {
      public:
        /**
         * \brief Constructor for the BroydenEquationsEMMMeffConstrained class.
         * 
         * \param model Pointer to the EffectiveMassModel.
         */
        BroydenEquationsEMMMeffConstrained(const EffectiveMassModel *model)
                : BroydenEquations(), m_EMM(model) { m_N = 1; }

        /**
         * \brief Implements the equations to be solved.
         * 
         * \param x The vector of variables.
         * \return The vector of equation values.
         */
        std::vector<double> Equations(const std::vector<double> &x);

      private:
        const EffectiveMassModel *m_EMM;
      };

    };

} // namespace thermalfist

#endif
