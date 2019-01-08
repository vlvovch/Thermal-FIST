/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef RANDOMGENERATORS_H
#define RANDOMGENERATORS_H

#include <cmath>

#include "MersenneTwister/MersenneTwister.h"
#include "HRGEventGenerator/MomentumDistribution.h"
#include "HRGBase/ThermalParticle.h"

namespace thermalfist {

  /// \brief Contains random generator functions used in
  ///        the Monte Carlo Thermal Event Generator 
  namespace RandomGenerators {

    /// \brief The Mersenne Twister random number generator
    extern MTRand randgenMT;

    /// \brief Generates random integer distributed by Poisson with specified mean
    /// Uses randgenMT
    /// \param mean Mean of the Poisson distribution
    int RandomPoisson(double mean);

    /// \brief Same as randgenMT(double) but uses the provided instance
    ///        of the  Mersenne Twister random number generator
    /// \param mean \copydoc mean
    /// \param rangen A Mersenne Twister random number generator to use
    int RandomPoisson(double mean, MTRand &rangen);


    /// \brief Base class for Monte Carlo sampling of particle momenta
    class ParticleMomentumGenerator
    {
    public:
      /// Default constructor
      ParticleMomentumGenerator() {}

      /// Destructor
      virtual ~ParticleMomentumGenerator() { }

      /// Samples the 3-momentum of a particle
      /// \return std::vector<double> A vector containing the sampled
      ///         \f$p_x\f$, \f$p_y\f$, \f$p_z\f$ components of the three-momentum
      virtual std::vector<double> GetMomentum() const = 0;
    };

    /**
     * \brief Class for generating the momentum of a particle 
     *              in accordance with the Siemens-Rasmussen formula.
     * 
     * Used in the spherically symmetric blast-wave event generator.
     * 
     */
    class SiemensRasmussenMomentumGenerator
      : public ParticleMomentumGenerator
    {
    public:

      SiemensRasmussenMomentumGenerator() { }

      /**
       * \brief Construct a new SiemensRasmussenGenerator object
       * 
       * \param T    The kinetic temperature (in GeV)
       * \param beta Transverse flow velocity
       * \param mass Particle mass (in GeV)
       */
      SiemensRasmussenMomentumGenerator(double T, double beta, double mass) :m_T(T), m_Beta(beta), m_Mass(mass) {
        m_Gamma = 1. / sqrt(1. - m_Beta * m_Beta);
        FixParameters();
      }

      ~SiemensRasmussenMomentumGenerator() { }

      /**
       * \brief Sets the parameters of the Siemens-Rasmussen distribution
       * 
       * \param T    The kinetic temperature (in GeV)
       * \param beta Transverse flow velocity
       * \param mass Particle mass (in GeV)
       */
      void SetParameters(double T, double beta, double mass) {
        m_T = T;
        m_Beta = beta;
        m_Mass = mass;
        m_Gamma = 1. / sqrt(1. - m_Beta * m_Beta);
        FixParameters();
      }

      // Override functions begin

      std::vector<double> GetMomentum() const;

      // Override functions end

    private:
      double w(double p) const {
        return sqrt(p*p + m_Mass * m_Mass);
      }

      double alpha(double p) const {
        return m_Gamma * m_Beta * p / m_T;
      }

      /// Unnormalized probability density of x = exp(-p)
      double g(double x) const;
      double g2(double x, double tp) const;

      /// Finds the maximum of g(x) using the ternary search
      void FixParameters();

      /// Generates random momentum p from Siemens-Rasmussen distribution
      /// Initially x = exp(-p) is generated in [0,1] where p is given in GeV
      /// Then p is recovered as p = -log(x)
      double GetRandom() const;

      double m_T;
      double m_Beta;
      double m_Mass;
      double m_Gamma;
      double m_Norm;
      double m_Max;
    };


    // Class for generating momentum of particle with mass m in accordance with Schnedermann-Sollfrank-Heinz formula with temperature T, transverse flow beta, and longitudinal flow etamax
    /**
     * \brief Class for generating the momentum of a particle 
     *              in accordance with the longitudinally symmetric
     *              blast-wave model.
     * 
     * Used in the longitudinally symmetric blast-wave event generator.
     * 
     */
    class SSHMomentumGenerator
      : public ParticleMomentumGenerator
    {
    public:
      SSHMomentumGenerator() { }

      /**
       * \brief Construct a new SSHGenerator object
       * \param T      The kinetic temperature (in GeV)
       * \param beta   The transverse flow velocity
       * \param etamax The longitudinal space-time rapidity cut-off
       * \param npow   The power in the transverse flow profile function
       * \param mass   Particle mass (in GeV) 
       */
      SSHMomentumGenerator(double T, double beta, double etamax, double npow, double mass) :m_T(T), m_Beta(beta), m_EtaMax(etamax), m_n(npow), m_Mass(mass) {
        m_distr = SSHDistribution(0, m_Mass, m_T, m_Beta, m_EtaMax, m_n, false);
        m_dPt = 0.02;
        m_dy = 0.05;
        FixParameters2();

      }
      
      ~SSHMomentumGenerator() { }

      /**
       * \brief Sets the parameters of the distribution
       * 
       * \param T      The kinetic temperature (in GeV)
       * \param beta   The transverse flow velocity
       * \param etamax The longitudinal space-time rapidity cut-off
       * \param npow   The power in the transverse flow profile function
       * \param mass   Particle mass (in GeV) 
       */
      void SetParameters(double T, double beta, double etamax, double npow, double mass) {
        m_T = T;
        m_Beta = beta;
        m_EtaMax = etamax;
        m_Mass = mass;
        m_distr = SSHDistribution(0, m_Mass, m_T, m_Beta, m_EtaMax, m_n);
        m_dPt = 0.02;
        m_dy = 0.05;
        FixParameters2();
      }

      // Override functions begin

      std::vector<double> GetMomentum() const;

      // Override functions end

    private:
      double w(double p) const {
        return sqrt(p*p + m_Mass * m_Mass);
      }

      double g(double x) const {
        return m_distr.dndpt(-log(x)) / x;
      }

      double g2(double x) const {
        return m_dndpt.f(x);
      }

      void FixParameters();

      void FixParameters2();

      // Finds the maximum of g(x) using the ternary search
      void FindMaximumPt();

      void FindMaximumY(double pt);

      void FindMaximumY2(double pt);

      // Generates random pt and y
      std::pair<double, double> GetRandom();

      std::pair<double, double> GetRandom2() const;

      double m_T, m_Beta, m_EtaMax, m_n, m_Mass;
      double m_MaxY, m_MaxPt;
      SSHDistribution m_distr;
      SplineFunction m_dndpt;
      std::vector<SplineFunction> m_dndy;
      std::vector<double> m_MaxYs;
      double m_dPt;
      double m_dy;
    };


    // Class for generating mass of resonance in accordance with the relativistic Breit-Wigner distribution
    class BreitWignerGenerator
    {
    public:
      BreitWignerGenerator() { }
      BreitWignerGenerator(double M, double gamma, double mthr) : m_M(M), m_Gamma(gamma), m_Mthr(mthr) { FixParameters(); }
      ~BreitWignerGenerator() { }

      void SetParameters(double M, double gamma, double mthr);

      // Unnormalized probability density of x
      double f(double x) const;

      // Finds the maximum of f(x) using the ternary search
      void FixParameters();

      // Generates random resonance mass m from relativistic Breit-Wigner distribution
      double GetRandom() const;

    private:
      double m_M;
      double m_Gamma;
      double m_Mthr;
      double m_Norm;
      double m_Max;
    };


    // Class for generating mass of resonance in accordance with its fixed Breit-Wigner distribution multiplied by the Boltzmann factor
    class ThermalBreitWignerGenerator
    {
    public:
      ThermalBreitWignerGenerator() { }
      ThermalBreitWignerGenerator(ThermalParticle *part, double T, double Mu) : m_part(part), m_T(T), m_Mu(Mu) { FixParameters(); }
      virtual ~ThermalBreitWignerGenerator() { }

      void SetParameters(ThermalParticle *part, double T, double Mu);


      virtual void FixParameters();

      // Unnormalized probability density of thermal resonance mass M
      virtual double f(double M) const;

      // Generates random resonance mass m from relativistic Breit-Wigner distribution
      double GetRandom() const;

    protected:
      ThermalParticle *m_part;
      double m_T, m_Mu;
      double m_Xmin, m_Xmax;
      double m_Max;

      /*double m_M;
      double m_Gamma;
      double m_Mthr;
      double m_Norm;
      double m_Max;*/
    };

    // Class for generating mass of resonance in accordance with its fixed Breit-Wigner distribution multiplied by the Boltzmann factor
    class ThermalEnergyBreitWignerGenerator
      : public ThermalBreitWignerGenerator
    {
    public:
      ThermalEnergyBreitWignerGenerator() : ThermalBreitWignerGenerator() { }
      ThermalEnergyBreitWignerGenerator(ThermalParticle *part, double T, double Mu) :
        ThermalBreitWignerGenerator(part, T, Mu) {
        FixParameters();
      }
      ~ThermalEnergyBreitWignerGenerator() { }

      virtual void FixParameters();

      // Unnormalized probability density of thermal resonance mass M
      virtual double f(double M) const;
    };
  }

} // namespace thermalfist


#endif
