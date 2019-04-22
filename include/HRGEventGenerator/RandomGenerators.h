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

#include "MersenneTwister.h"
#include "HRGEventGenerator/MomentumDistribution.h"
#include "HRGBase/ThermalParticle.h"

namespace thermalfist {

  /// \brief Contains random generator functions used in
  ///        the Monte Carlo Thermal Event Generator 
  namespace RandomGenerators {

    /// \brief The Mersenne Twister random number generator
    extern MTRand randgenMT;

    /// \brief Set the seed of the random number generator randgenMT
    void SetSeed(const unsigned int seed);

    /// \brief Generates random integer distributed by Poisson with specified mean
    /// Uses randgenMT
    /// \param mean Mean of the Poisson distribution
    int RandomPoisson(double mean);

    /// \brief Same as randgenMT(double) but uses the provided instance
    ///        of the  Mersenne Twister random number generator
    /// \param mean \copydoc mean
    /// \param rangen A Mersenne Twister random number generator to use
    int RandomPoisson(double mean, MTRand &rangen);

    /// \brief Probability of a Skellam distributed random variable with Poisson means
    ///        mu1 and mu2 to have the value of k.
    double SkellamProbability(int k, double mu1, double mu2);


    /// \brief Generator of a random number from the Bessel distribution (a, nu), nu is integer
    ///        Uses methods from https://www.sciencedirect.com/science/article/pii/S016771520200055X
    ///        Used in event generator with exact conservation of charges to
    ///        generate two Poisson numbers with fixed difference, as described in https://arxiv.org/pdf/1609.01087.pdf
    struct BesselDistributionGenerator
    {
      static double pn(int n, double a, int nu);

      //static double R(double x, int nu) { return xMath::BesselI(nu + 1, x) / xMath::BesselI(nu, x); }
      static double R(double x, int nu);

      static double mu(double a, int nu) { return a * R(a, nu) / 2; }

      static double chi2(double a, int nu);

      static int    m(double a, int nu) { return static_cast<int>((sqrt(a*a+static_cast<double>(nu)*nu) - nu)/2.); }

      static double sig2(double a, int nu);

      static double Q2(double a, int nu);

      static int RandomBesselPoisson(double a, int nu, MTRand &rangen);

      static int RandomBesselPoisson(double a, int nu) { return RandomBesselPoisson(a, nu, randgenMT); }

      static double pmXmOverpm(int X, int tm, double a, int nu);

      static int RandomBesselDevroye1(double a, int nu, MTRand &rangen);

      static int RandomBesselDevroye1(double a, int nu) { return RandomBesselDevroye1(a, nu, randgenMT); }

      static int RandomBesselDevroye2(double a, int nu, MTRand &rangen);

      static int RandomBesselDevroye2(double a, int nu) { return RandomBesselDevroye2(a, nu, randgenMT); }

      static int RandomBesselDevroye3(double a, int nu, MTRand &rangen);

      static int RandomBesselDevroye3(double a, int nu) { return RandomBesselDevroye3(a, nu, randgenMT); }

      static int RandomBesselNormal(double a, int nu, MTRand &rangen);

      static int RandomBesselNormal(double a, int nu) { return RandomBesselNormal(a, nu, randgenMT); }

      static int RandomBesselCombined(double a, int nu, MTRand &rangen);

      static int RandomBesselCombined(double a, int nu) { return RandomBesselCombined(a, nu, randgenMT); }
    };


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
        m_n = npow;
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


    /// \brief Class for generating mass of resonance in accordance 
    ///        with the relativistic Breit-Wigner distribution
    class BreitWignerGenerator
    {
    public:

      BreitWignerGenerator() { }

      /**
       * \brief Construct a new BreitWignerGenerator object
       * 
       * Currently not used.
       * 
       * \param M     Pole mass of the resonance (in GeV)
       * \param gamma Width of the resonance (in GeV)
       * \param mthr  The threshold mass of the resonance (in GeV)
       */
      BreitWignerGenerator(double M, double gamma, double mthr) : m_M(M), m_Gamma(gamma), m_Mthr(mthr) { FixParameters(); }
      
      ~BreitWignerGenerator() { }

      /**
       * \brief Set the Breit-Wigner spectral function parameters
       * 
       * \param M     Pole mass of the resonance (in GeV)
       * \param gamma Width of the resonance (in GeV)
       * \param mthr  The threshold mass of the resonance (in GeV)
       */
      void SetParameters(double M, double gamma, double mthr);

      // Generates random resonance mass m from relativistic Breit-Wigner distribution
      /**
       * \brief Samples the resonance mass from a
       *        relativistic Breit-Wigner distribution
       *        with a constant width
       * 
       * \return The sampled mass (in GeV)
       */
      double GetRandom() const;

    private:
      /// Unnormalized probability density of x
      double f(double x) const;

      /// Finds the maximum of f(x) using the ternary search
      void FixParameters();

      double m_M;
      double m_Gamma;
      double m_Mthr;
      double m_Norm;
      double m_Max;
    };


    /// \brief Class for generating mass of resonance 
    ///        in accordance with constant width Breit-Wigner distribution 
    ///        multiplied by the thermal density.
    ///
    /// Samples

    /**
     * \brief Class for generating mass of resonance
     *        in accordance with the constant width Breit-Wigner distribution
     *        multiplied by the thermal density.
     * 
     * Sample the mass from the distribution
     * \f[
     *   \rho(M) \sim \rho_{\rm BW} (M) \, n_{\rm th}^{\rm id} (T,\mu;M)~.
     * \f]
     * 
     */
    class ThermalBreitWignerGenerator
    {
    public:
      ThermalBreitWignerGenerator() { }

      /**
       * \brief Construct a new ThermalBreitWignerGenerator object
       * 
       * \param part A pointer to the ThermalParticle object representing
       *             the species sampled.
       * \param T    Tempeature (in GeV)
       * \param Mu   Chemical potential of the sampled particle (in GeV)
       */
      ThermalBreitWignerGenerator(ThermalParticle *part, double T, double Mu) : m_part(part), m_T(T), m_Mu(Mu) { FixParameters(); }
      
      virtual ~ThermalBreitWignerGenerator() { }

      /**
       * \brief Sets the parameters of the distribution.
       * 
       * \param part A pointer to the ThermalParticle object representing
       *             the species sampled.
       * \param T    Tempeature (in GeV)
       * \param Mu   Chemical potential of the sampled particle (in GeV)
       */
      void SetParameters(ThermalParticle *part, double T, double Mu);

      /**
       * \brief Samples the mass.
       * 
       * \return double The sampled mass (in GeV)
       */
      double GetRandom() const;

    protected:
      /// Computes some auxiliary stuff needed for sampling
      virtual void FixParameters();

      /// Unnormalized resonance mass probability density
      virtual double f(double M) const;

      ThermalParticle *m_part;
      double m_T, m_Mu;
      double m_Xmin, m_Xmax;
      double m_Max;
    };

    // Class for generating mass of resonance in accordance with its fixed Breit-Wigner distribution multiplied by the Boltzmann factor
    /**
     * \brief Class for generating mass of resonance
     *        in accordance with the energy-dependent Breit-Wigner distribution
     *        multiplied by the thermal density.
     * 
     * Sample the mass from the distribution
     * \f[
     *   \rho(M) \sim \rho_{\rm eBW} (M) \, n_{\rm th}^{\rm id} (T,\mu;M)~.
     * \f]
     * 
     */
    class ThermalEnergyBreitWignerGenerator
      : public ThermalBreitWignerGenerator
    {
    public:

      ThermalEnergyBreitWignerGenerator() : ThermalBreitWignerGenerator() { }

      /**
       * \brief Construct a new ThermalEnergyBreitWignerGenerator object
       * 
       * \copydetails ThermalBreitWignerGenerator(ThermalParticle*,double,double)
       */
      ThermalEnergyBreitWignerGenerator(ThermalParticle *part, double T, double Mu) :
        ThermalBreitWignerGenerator(part, T, Mu) {
        FixParameters();
      }

      ~ThermalEnergyBreitWignerGenerator() { }

    protected:
      virtual void FixParameters();
      
      virtual double f(double M) const;
    };
  }

} // namespace thermalfist


#endif
