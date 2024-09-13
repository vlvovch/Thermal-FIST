/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2017-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef IDEALGASFUNCTIONS_H
#define IDEALGASFUNCTIONS_H

/**
 * \file IdealGasFunctions.h
 * 
 * A collection of thermodynamic functions
 * corresponding to an ideal Maxwell-Boltzmann,
 * Fermi-Dirac, or Bose-Einstein gas in the
 * grand canonical ensemble.
 * 
 */

namespace thermalfist {

  /// \brief Contains implementation of the thermodynamic functions
  ///        corresponding to an ideal gas of particles in the grand-canonical
  ///        ensemble.
  namespace IdealGasFunctions {
    /**
     * \brief Identifies the thermodynamic function.
     * 
     */
    enum Quantity { ParticleDensity, EnergyDensity, EntropyDensity, Pressure, chi2, chi3, chi4, ScalarDensity, 
      chi2difull, chi3difull, chi4difull,
      dndT, d2ndT2, dedT, dedmu, dchi2dT
    };
    /**
     * \brief Identifies whether quantum statistics
     *        are to be computed using the cluster expansion
     *        or numerical integration using 32-point Gauss-Legendre quadratures.
     * 
     */
    enum QStatsCalculationType { ClusterExpansion, Quadratures };

    /// \brief Whether \mu > m Bose-Einstein condensation issue was encountered for a Bose gas
    extern bool calculationHadBECIssue;

    /// Magnetic field configuration
    struct MagneticFieldConfiguration {
      double B;        ///< Magnetic field strength [GeV^2].
      int lmax;        ///< Number of the Landau levels.
      double degSpin;  ///< Particle's spin degeneracy.
      double Q;        ///< Particle's electric charge [e].
      MagneticFieldConfiguration(double in_B = 0.,
                                 int in_lmax = 100,
                                 double in_degSpin = 1.,
                                 double in_Q = 0.) :
                                 B(in_B),lmax(in_lmax),degSpin(in_degSpin),Q(in_Q)
                                 { }
    };



    struct IdealGasFunctionsExtraConfig {
      MagneticFieldConfiguration MagneticField;  ///< Magnetic field configuration.
      IdealGasFunctionsExtraConfig(const MagneticFieldConfiguration& in_MagneticField = MagneticFieldConfiguration()) :
                                   MagneticField(in_MagneticField)
                                   { }
    };

    /**
     * \brief Computes the particle number density of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double BoltzmannDensity(double T, double mu, double m, double deg,
                            const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the pressure of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double BoltzmannPressure(double T, double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the energy density of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double BoltzmannEnergyDensity(double T, double mu, double m, double deg,
                                  const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the entropy density of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double BoltzmannEntropyDensity(double T, double mu, double m, double deg,
                                   const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the scalar density of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double BoltzmannScalarDensity(double T, double mu, double m, double deg,
                                  const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());  // TODO: Check for correctness
    
    /**
     * \brief Computes the chemical potential derivative of density for a Maxwell-Boltzmann gas.
     * 
     * Computes \f$ T \frac{\partial n}{\partial mu} \f$ for a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ T \frac{\partial n}{\partial mu} \f$ [GeV fm-3].
     */
    double BoltzmannTdndmu(int N, double T, double mu, double m, double deg,
                           const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the n-th order dimensionless susceptibility for a Maxwell-Boltzmann gas.
     * 
     * Computes \f$ \chi_n \equiv \frac{\partial^n p/T^4}{\partial (mu/T)^n} \f$ 
     * for a Maxwell-Boltzmann gas.
     * 
     * \param N Susceptibility order.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$.
     */
    double BoltzmannChiN(int N, double T, double mu, double m, double deg,
                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the n-th order dimensionfull susceptibility for a Maxwell-Boltzmann gas.
     *
     * Computes \f$ \chi_n \equiv \frac{\partial^n p}{\partial mu^n} \f$
     * for a Maxwell-Boltzmann gas.
     *
     * \param N Susceptibility order.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$.
     */
    double BoltzmannChiNDimensionfull(int N, double T, double mu, double m, double deg,
                                      const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the thermal part of the magnetization of a Maxwell-Boltzmann gas, m_B = dP/dB.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \param extraConfig Extra parameters such as magnetic field.
     * \return Magnetization [GeV^2].
     */
    double BoltzmannMagnetization(double T, double mu, double m, double deg,
                            const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the particle number density of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double QuantumClusterExpansionDensity(int statistics, double T, double mu, double m, double deg, int order = 1,
                                          const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the pressure of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double QuantumClusterExpansionPressure(int statistics, double T, double mu, double m, double deg, int order = 1,
                                           const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the energy density of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double QuantumClusterExpansionEnergyDensity(int statistics, double T, double mu, double m, double deg, int order = 1,
                                                const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());  // TODO: Check for correctness
    
    /**
     * \brief Computes the entropy density of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double QuantumClusterExpansionEntropyDensity(int statistics, double T, double mu, double m, double deg, int order = 1,
                                                 const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the scalar density of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double QuantumClusterExpansionScalarDensity(int statistics, double T, double mu, double m, double deg, int order = 1,
                                                const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());   // TODO: Check for correctness
    
    /**
     * \brief Computes the chemical potential derivative of density for a quantum ideal gas using cluster expansion.
     * 
     * Computes \f$ T \frac{\partial n}{\partial mu} \f$ for a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ T \frac{\partial n}{\partial mu} \f$ [GeV fm-3].
     */
    double QuantumClusterExpansionTdndmu(int N, int statistics, double T, double mu, double m, double deg, int order = 1,
                                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the n-th order dimensionless susceptibility for a quantum ideal gas using cluster expansion.
     * 
     * Computes \f$ \chi_n \equiv \frac{\partial^n p/T^4}{\partial (mu/T)^n} \f$ 
     * for a quantum ideal gas using cluster expansion.
     * 
     * \param N Susceptibility order.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$.
     */
    double QuantumClusterExpansionChiN(int N, int statistics, double T, double mu, double m, double deg, int order = 1,
                                       const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the n-th order dimensionfull susceptibility for a quantum ideal gas using cluster expansion.
     *
     * Computes \f$ \chi_n \equiv \frac{\partial^n p}{\partial mu^n} \f$
     * for a quantum ideal gas using cluster expansion.
     *
     * \param N Susceptibility order.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$ [GeV^{4-N}].
     */
    double QuantumClusterExpansionChiNDimensionfull(int N, int statistics, double T, double mu, double m, double deg, int order = 1,
                                                    const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());



    /**
     * \brief Computes the thermal part of the magnetization of a Maxwell-Boltzmann gas, m_B = dP/dB.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Magnetization [GeV^2].
     */
    double QuantumClusterExpansionMagnetization(int statistics, double T, double mu, double m, double deg, int order = 1,
                                                const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the particle number density of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double QuantumNumericalIntegrationDensity(int statistics, double T, double mu, double m, double deg,
                                              const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the pressure of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double QuantumNumericalIntegrationPressure(int statistics, double T, double mu, double m, double deg,
                                               const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the energy density of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double QuantumNumericalIntegrationEnergyDensity(int statistics, double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the entropy density of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double QuantumNumericalIntegrationEntropyDensity(int statistics, double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the scalar density of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double QuantumNumericalIntegrationScalarDensity(int statistics, double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());  // TODO: Check for correctness

    /**
     * \brief Computes the magnetization of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Magnetization [GeV^2].
     */
    double QuantumNumericalIntegrationMagnetization(int statistics, double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    double QuantumNumericalIntegrationT1dn1dmu1(int statistics, double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double QuantumNumericalIntegrationT2dn2dmu2(int statistics, double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double QuantumNumericalIntegrationT3dn3dmu3(int statistics, double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double QuantumNumericalIntegrationTdndmu(int N, int statistics, double T, double mu, double m, double deg,
                                             const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the n-th order dimensionless susceptibility for a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * Computes \f$ \chi_n \equiv \frac{\partial^n p/T^4}{\partial (mu/T)^n} \f$ 
     * for a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param N Susceptibility order.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$.
     */
    double QuantumNumericalIntegrationChiN(int N, int statistics, double T, double mu, double m, double deg,
                                           const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the n-th order dimensionfull susceptibility for a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     *
     * Computes \f$ \chi_n \equiv \frac{\partial^n p}{\partial mu^n} \f$
     * for a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     *
     * \param N Susceptibility order.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$ [GeV^{n-4}.
     */
    double QuantumNumericalIntegrationChiNDimensionfull(int N, int statistics, double T, double mu, double m, double deg,
                                                        const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Auxiliary function.
     * 
     * Auxiliary function used for Fermi-Dirac integrals at m > mu.
     * See Eq. (A.2) in [https://arxiv.org/pdf/0901.1430.pdf](https://arxiv.org/pdf/0901.1430.pdf)
     * 
     * \param x Variable.
     * \return Function value.
     */
    double psi(double x);

    /**
     * \brief Auxiliary function.
     * 
     * Auxiliary function used for Fermi-Dirac integrals at m > mu.
     * 
     * \param x Variable.
     * \return Function value.
     */
    double psi2(double x);

    /**
     * \brief Computes the particle number density of a Fermi-Dirac ideal gas 
     *        at mu > m.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double FermiNumericalIntegrationLargeMuDensity(double T, double mu, double m, double deg,
                                                   const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the pressure of a Fermi-Dirac ideal gas 
     *        at mu > m.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double FermiNumericalIntegrationLargeMuPressure(double T, double mu, double m, double deg,
                                                    const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the energy density of a Fermi-Dirac ideal gas 
     *        at mu > m.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double FermiNumericalIntegrationLargeMuEnergyDensity(double T, double mu, double m, double deg,
                                                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the entropy density of a Fermi-Dirac ideal gas 
     *        at mu > m.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double FermiNumericalIntegrationLargeMuEntropyDensity(double T, double mu, double m, double deg,
                                                          const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the scalar density of a Fermi-Dirac ideal gas 
     *        at mu > m.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double FermiNumericalIntegrationLargeMuScalarDensity(double T, double mu, double m, double deg,
                                                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());  // TODO: Check for correctness


    double FermiNumericalIntegrationLargeMuMagnetization(double T, double mu, double m, double deg,
                                                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());  // TODO: Check for correctness

    double FermiNumericalIntegrationLargeMuT1dn1dmu1(double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double FermiNumericalIntegrationLargeMuT2dn2dmu2(double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double FermiNumericalIntegrationLargeMuT3dn3dmu3(double T, double mu, double m, double deg,
                                                     const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double FermiNumericalIntegrationLargeMuTdndmu(int N, double T, double mu, double m, double deg,
                                                  const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    
    /**
     * \brief Computes the n-th order dimensionless susceptibility for a Fermi-Dirac ideal gas 
     *        at mu > m.
     * 
     * Computes \f$ \chi_n \equiv \frac{\partial^n p/T^4}{\partial (mu/T)^n} \f$ 
     * for a Fermi-Dirac ideal gas at mu > m.
     * 
     * \param N Susceptibility order.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$.
     */
    double FermiNumericalIntegrationLargeMuChiN(int N, double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the n-th order dimensionfull susceptibility for a Fermi-Dirac ideal gas
     *        at mu > m.
     *
     * Computes \f$ \chi_n \equiv \frac{\partial^n p}{\partial mu^n} \f$
     * for a Fermi-Dirac ideal gas at mu > m.
     *
     * \param N Susceptibility order.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$.
     */
    double FermiNumericalIntegrationLargeMuChiNDimensionfull(int N, double T, double mu, double m, double deg,
                                                             const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());


    /**
     * \brief Computes the particle number density of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double FermiZeroTDensity(double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the pressure of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double FermiZeroTPressure(double mu, double m, double deg,
                              const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the energy density of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double FermiZeroTEnergyDensity(double mu, double m, double deg,
                                   const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the entropy density of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double FermiZeroTEntropyDensity(double mu, double m, double deg,
                                    const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the scalar density of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double FermiZeroTScalarDensity(double mu, double m, double deg,
                                   const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());  // TODO: Check for correctness

    double FermiZeroTMagnetization(double mu, double m, double deg,
                                   const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());  // TODO: Check for correctness


    double FermiZeroTdn1dmu1(double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double FermiZeroTdn2dmu2(double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double FermiZeroTdn3dmu3(double mu, double m, double deg,
                             const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
    double FermiZeroTdndmu(int N, double mu, double m, double deg,
                           const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the n-th order susceptibility for a Fermi-Dirac ideal gas at zero temperature.
     *
     * Computes \f$ \chi_n \equiv \frac{\partial^n p/T^4}{\partial (mu/T)^n} \f$
     * for a Fermi-Dirac ideal gas at mu > m.
     *
     * \param N Susceptibility order.
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$.
     */
    double FermiZeroTChiN(int N, double mu, double m, double deg,
                          const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the n-th order susceptibility scaled by T^n for a Fermi-Dirac ideal gas at zero temperature.
     *
     * Computes \f$ \chi_n \equiv \frac{\partial^n p}{\partial (mu)^n} \f$
     * for a Fermi-Dirac ideal gas at mu > m.
     * Useful for zero temperature where the dimensionless susceptibility would be singular.
     *
     * \param N Susceptibility order.
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return \f$ \chi_n \f$ [GeV^{4-N}].
     */
    double FermiZeroTChiNDimensionfull(int N, double mu, double m, double deg,
                                       const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Calculation of a generic ideal gas function.
     * 
     * Calculations of the ideal gas thermodynamic function
     * specified by \param quantity
     * 
     * \param quantity Identifies the thermodynamic function to calculate.
     * \param calctype Method used to perform the calculation if quantum statistics used.
     * \param statistics 0 -- Maxwell-Boltzmann, +1 -- Fermi-Dirac, -1 -- Bose-Einstein.
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \param order Number of terms in the cluster expansion if this method is used.
     * \return Computed thermodynamic function.
     */
    double IdealGasQuantity(Quantity quantity, QStatsCalculationType calctype, int statistics, double T, double mu, double m, double deg, int order = 1,
                            const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /// Temperature derivatives
    /**
     * \brief Computes the derivative of particle number density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dndT [fm-3 * GeV^-1].
     */
    double BoltzmanndndT(double T, double mu, double m, double deg,
                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
    * \brief Computes the 2nd derivative of particle number density wrt T at constant mu for a Maxwell-Boltzmann gas.
    *
    * \param T Temperature [GeV].
    * \param mu Chemical potential [GeV].
    * \param m  Particle's mass [GeV].
    * \param deg Internal degeneracy factor.
    * \return dndT [fm-3 * GeV^-2].
    */
    double Boltzmannd2ndT2(double T, double mu, double m, double deg,
                           const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of energy density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dedT [fm-3].
     */
    double BoltzmanndedT(double T, double mu, double m, double deg,
                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of energy density wrt mu at constant T for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dedmu [fm-3].
     */
    double Boltzmanndedmu(double T, double mu, double m, double deg,
                          const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of the susceptibility wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dchi2dT [fm-3].
     */
    double Boltzmanndchi2dT(double T, double mu, double m, double deg,
                            const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    // Cluster expansion
    /**
     * \brief Computes the derivative of particle number density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dndT [fm-3 * GeV^-1].
     */
    double QuantumClusterExpansiondndT(int statistics, double T, double mu, double m, double deg, int order = 1,
                                       const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
    * \brief Computes the 2nd derivative of particle number density wrt T at constant mu for a Maxwell-Boltzmann gas.
    *
    * \param T Temperature [GeV].
    * \param mu Chemical potential [GeV].
    * \param m  Particle's mass [GeV].
    * \param deg Internal degeneracy factor.
    * \return dndT [fm-3 * GeV^-2].
    */
    double QuantumClusterExpansiond2ndT2(int statistics, double T, double mu, double m, double deg, int order = 1,
                                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of energy density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dedT [fm-3].
     */
    double QuantumClusterExpansiondedT(int statistics, double T, double mu, double m, double deg, int order = 1,
                                       const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of energy density wrt mu at constant T for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dedmu [fm-3].
     */
    double QuantumClusterExpansiondedmu(int statistics, double T, double mu, double m, double deg, int order = 1,
                                        const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of the susceptibility wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dchi2dT [fm-3].
     */
    double QuantumClusterExpansiondchi2dT(int statistics, double T, double mu, double m, double deg, int order = 1,
                                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());


    /**
     * \brief Computes the derivative of particle number density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dndT [fm-3 * GeV^-1].
     */
    double QuantumNumericalIntegrationdndT(int statistics, double T, double mu, double m, double deg,
                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
    * \brief Computes the 2nd derivative of particle number density wrt T at constant mu for a Maxwell-Boltzmann gas.
    *
    * \param T Temperature [GeV].
    * \param mu Chemical potential [GeV].
    * \param m  Particle's mass [GeV].
    * \param deg Internal degeneracy factor.
    * \return dndT [fm-3 * GeV^-2].
    */
    double QuantumNumericalIntegrationd2ndT2(int statistics, double T, double mu, double m, double deg,
                           const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of energy density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dedT [fm-3].
     */
    double QuantumNumericalIntegrationdedT(int statistics, double T, double mu, double m, double deg,
                         const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of energy density wrt mu at constant T for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dedmu [fm-3].
     */
    double QuantumNumericalIntegrationdedmu(int statistics, double T, double mu, double m, double deg,
                          const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of the susceptibility wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dchi2dT [fm-3].
     */
    double QuantumNumericalIntegrationdchi2dT(int statistics, double T, double mu, double m, double deg,
                            const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of particle number density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dndT [fm-3 * GeV^-1].
     */
    double FermiNumericalIntegrationLargeMudndT(double T, double mu, double m, double deg,
                                           const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the second derivative of particle number density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dndT [fm-3 * GeV^-2].
     */
    double FermiNumericalIntegrationLargeMud2ndT2(double T, double mu, double m, double deg,
                                                const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of energy density wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dedT [fm-3].
     */
    double FermiNumericalIntegrationLargeMudedT(double T, double mu, double m, double deg,
                                           const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of energy density wrt mu at constant T for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dedmu [fm-3].
     */
    double FermiNumericalIntegrationLargeMudedmu(double T, double mu, double m, double deg,
                                            const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());

    /**
     * \brief Computes the derivative of the susceptibility wrt T at constant mu for a Maxwell-Boltzmann gas.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return dchi2dT [fm-3].
     */
    double FermiNumericalIntegrationLargeMudchi2dT(double T, double mu, double m, double deg,
                                              const IdealGasFunctionsExtraConfig& extraConfig = IdealGasFunctionsExtraConfig());
  }

  /// \brief Implements the possibility of a generalized calculation of the densities.
  ///        For example, effective mass model.
  ///        Abstract class.
  class GeneralizedDensity {
  public:
    GeneralizedDensity() {}
    virtual ~GeneralizedDensity() {}

    virtual double Quantity(IdealGasFunctions::Quantity quantity, double T, double mu) = 0;
    virtual double EffectiveMass() const { return -1.; }
    virtual bool   IsBECPhase() const { return false; }
    virtual double BECFraction() const { return 0.; }
  };

} // namespace thermalfist

#endif // IDEALGASFUNCTIONS_H
