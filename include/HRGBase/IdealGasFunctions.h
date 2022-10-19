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
      chi2difull, chi3difull, chi4difull
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

    /**
     * \brief Computes the particle number density of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double BoltzmannDensity(double T, double mu, double m, double deg);

    /**
     * \brief Computes the pressure of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double BoltzmannPressure(double T, double mu, double m, double deg);

    /**
     * \brief Computes the energy density of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double BoltzmannEnergyDensity(double T, double mu, double m, double deg);

    /**
     * \brief Computes the entropy density of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double BoltzmannEntropyDensity(double T, double mu, double m, double deg);

    /**
     * \brief Computes the scalar density of a Maxwell-Boltzmann gas.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double BoltzmannScalarDensity(double T, double mu, double m, double deg);  // TODO: Check for correctness
    
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
    double BoltzmannTdndmu(int N, double T, double mu, double m, double deg);

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
    double BoltzmannChiN(int N, double T, double mu, double m, double deg);

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
    double BoltzmannChiNDimensionfull(int N, double T, double mu, double m, double deg);

    /**
     * \brief Computes the particle number density of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double QuantumClusterExpansionDensity(int statistics, double T, double mu, double m, double deg, int order = 1);
    
    /**
     * \brief Computes the pressure of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double QuantumClusterExpansionPressure(int statistics, double T, double mu, double m, double deg, int order = 1);
    
    /**
     * \brief Computes the energy density of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double QuantumClusterExpansionEnergyDensity(int statistics, double T, double mu, double m, double deg, int order = 1);  // TODO: Check for correctness
    
    /**
     * \brief Computes the entropy density of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double QuantumClusterExpansionEntropyDensity(int statistics, double T, double mu, double m, double deg, int order = 1);
    
    /**
     * \brief Computes the scalar density of a quantum ideal gas using cluster expansion.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double QuantumClusterExpansionScalarDensity(int statistics, double T, double mu, double m, double deg, int order = 1);   // TODO: Check for correctness
    
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
    double QuantumClusterExpansionTdndmu(int N, int statistics, double T, double mu, double m, double deg, int order = 1);
    
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
    double QuantumClusterExpansionChiN(int N, int statistics, double T, double mu, double m, double deg, int order = 1);

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
    double QuantumClusterExpansionChiNDimensionfull(int N, int statistics, double T, double mu, double m, double deg, int order = 1);

    /**
     * \brief Computes the particle number density of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double QuantumNumericalIntegrationDensity(int statistics, double T, double mu, double m, double deg);
    
    /**
     * \brief Computes the pressure of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double QuantumNumericalIntegrationPressure(int statistics, double T, double mu, double m, double deg);
    
    /**
     * \brief Computes the energy density of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double QuantumNumericalIntegrationEnergyDensity(int statistics, double T, double mu, double m, double deg);
    
    /**
     * \brief Computes the entropy density of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double QuantumNumericalIntegrationEntropyDensity(int statistics, double T, double mu, double m, double deg);
    
    /**
     * \brief Computes the scalar density of a quantum ideal gas using 32-point Gauss-Laguerre quadratures.
     * 
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double QuantumNumericalIntegrationScalarDensity(int statistics, double T, double mu, double m, double deg);  // TODO: Check for correctness
    
    double QuantumNumericalIntegrationT1dn1dmu1(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationT2dn2dmu2(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationT3dn3dmu3(int statistics, double T, double mu, double m, double deg);
    double QuantumNumericalIntegrationTdndmu(int N, int statistics, double T, double mu, double m, double deg);
    
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
    double QuantumNumericalIntegrationChiN(int N, int statistics, double T, double mu, double m, double deg);

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
    double QuantumNumericalIntegrationChiNDimensionfull(int N, int statistics, double T, double mu, double m, double deg);

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
    double FermiNumericalIntegrationLargeMuDensity(double T, double mu, double m, double deg);
    
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
    double FermiNumericalIntegrationLargeMuPressure(double T, double mu, double m, double deg);
    
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
    double FermiNumericalIntegrationLargeMuEnergyDensity(double T, double mu, double m, double deg);
    
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
    double FermiNumericalIntegrationLargeMuEntropyDensity(double T, double mu, double m, double deg);
    
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
    double FermiNumericalIntegrationLargeMuScalarDensity(double T, double mu, double m, double deg);  // TODO: Check for correctness
    
    double FermiNumericalIntegrationLargeMuT1dn1dmu1(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuT2dn2dmu2(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuT3dn3dmu3(double T, double mu, double m, double deg);
    double FermiNumericalIntegrationLargeMuTdndmu(int N, double T, double mu, double m, double deg);
    
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
    double FermiNumericalIntegrationLargeMuChiN(int N, double T, double mu, double m, double deg);

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
    double FermiNumericalIntegrationLargeMuChiNDimensionfull(int N, double T, double mu, double m, double deg);


    /**
     * \brief Computes the particle number density of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Particle number density [fm-3].
     */
    double FermiZeroTDensity(double mu, double m, double deg);

    /**
     * \brief Computes the pressure of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Pressure [GeV fm-3].
     */
    double FermiZeroTPressure(double mu, double m, double deg);

    /**
     * \brief Computes the energy density of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Energy density [GeV fm-3].
     */
    double FermiZeroTEnergyDensity(double mu, double m, double deg);

    /**
     * \brief Computes the entropy density of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Entropy density [GeV fm-3].
     */
    double FermiZeroTEntropyDensity(double mu, double m, double deg);

    /**
     * \brief Computes the scalar density of a Fermi-Dirac ideal gas at zero temperature.
     *
     * \param T Temperature [GeV].
     * \param mu Chemical potential [GeV].
     * \param m  Particle's mass [GeV].
     * \param deg Internal degeneracy factor.
     * \return Scalar density [fm-3].
     */
    double FermiZeroTScalarDensity(double mu, double m, double deg);  // TODO: Check for correctness

    double FermiZeroTdn1dmu1(double mu, double m, double deg);
    double FermiZeroTdn2dmu2(double mu, double m, double deg);
    double FermiZeroTdn3dmu3(double mu, double m, double deg);
    double FermiZeroTdndmu(int N, double mu, double m, double deg);

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
    double FermiZeroTChiN(int N, double mu, double m, double deg);

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
    double FermiZeroTChiNDimensionfull(int N, double mu, double m, double deg);

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
    double IdealGasQuantity(Quantity quantity, QStatsCalculationType calctype, int statistics, double T, double mu, double m, double deg, int order = 1);
  }

} // namespace thermalfist

#endif // IDEALGASFUNCTIONS_H
