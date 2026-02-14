#ifndef THERMALMODELREALGAS_H
#define THERMALMODELREALGAS_H

#include "HRGBase/ThermalModelBase.h"
#include "ExcludedVolumeModelsMulti.h"
#include "MeanFieldModelsMulti.h"

namespace thermalfist {

	/**
 * \brief Class implementing the quantum real gas HRG model.
 *
 * The model formulation can be found in
 *
 * V. Vovchenko, 
 * Phys. Rev. C **96**, 015206 (2017),
 * [https://arxiv.org/pdf/1701.06524](https://arxiv.org/pdf/1701.06524)
 *
 * The system of transcendental equations for
 * the "shifted" chemical potentials of hadrons
 * is solved using the Broyden's method.
 *
 */
	class ThermalModelRealGas : public ThermalModelBase
	{
	public:
		/**
		 * \brief Construct a new ThermalModelRealGas object
		 *
		 * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
		 * \param params ThermalModelParameters object with current thermal parameters
		 */
		ThermalModelRealGas(ThermalParticleSystem* TPS_, const ThermalModelParameters& params = ThermalModelParameters());

		/**
		 * \brief Destroy the ThermalModelRealGas object
		 *
		 */
		virtual ~ThermalModelRealGas(void);

		/**
		 * \brief Set the excluded volume model for the real gas HRG model.
		 *
		 * \param exvolmod Pointer to the excluded volume model to be used.
		 */
		void SetExcludedVolumeModel(ExcludedVolumeModelMultiBase* exvolmod) { m_exvolmod = exvolmod; m_ComponentMapCalculated = false; }

		/**
		 * \brief Set the mean field model for the real gas HRG model.
		 *
		 * \param mfmod Pointer to the mean field model to be used.
		 */
		void SetMeanFieldModel(MeanFieldModelMultiBase* mfmod) { m_mfmod = mfmod; m_ComponentMapCalculated = false; }

		/**
		 * \brief Get the excluded volume model used in the real gas HRG model.
		 *
		 * \return Pointer to the excluded volume model used.
		 */
		ExcludedVolumeModelMultiBase* ExcludedVolumeModel() const { return m_exvolmod; }

		/**
		 * \brief Get the mean field model used in the real gas HRG model.
		 *
		 * \return Pointer to the mean field model used.
		 */
		MeanFieldModelMultiBase* MeanFieldModel() const { return m_mfmod; }

		//void SetExcludedVolumeModelCopy(const ExcludedVolumeModelMultiBase& exvolmod) { *m_exvolmod = exvolmod; }
		//void SetMeanFieldModelCopy(const MeanFieldModelMultiBase& mfmod) { *m_mfmod = mfmod; }

		/**
		 * \brief Whether to search for multiple solutions of the real gas model equations
		 * by considering different initial guesses in the Broyden's method.
		 *
		 * Multiple solutions in the real gas model appear e.g. below the
		 * critical temperature of the liquid-gas phase transition.
		 *
		 * \param search Whether multiple solutions of the QvdW equations
		 *               should be considered. False by default.
		 */
		virtual void SetMultipleSolutionsMode(bool search) { m_SearchMultipleSolutions = search; }

		/**
		 * \brief Whether to search for multiple solutions of the real gas model equations
		 * by considering different initial guesses in the Broyden's method.
		 *
		 * \return true  Multiple solutions considered
		 * \return false Multiple solutions not considered
		 */
		virtual bool UseMultipleSolutionsMode() const { return m_SearchMultipleSolutions; }

		/**
		 * \brief Get the shifted chemical potential of a particle species.
		 *
		 * \param i 0-based particle species index.
		 * \return The shifted chemical potential of the particle species.
		 */
		double MuStar(int i) const { return m_MuStar[i]; }

		/**
		 * \brief Get the vector of shifted chemical potentials for all particle species.
		 *
		 * \return Vector of shifted chemical potentials, one element per particle species.
		 */
		std::vector<double> GetMuStar() const { return m_MuStar; }

		/**
		 * \brief Set the vector of shifted chemical potentials for all particle species.
		 *
		 * \param MuStar Vector of shifted chemical potentials, one element per particle species.
		 */
		void SetMuStar(const std::vector<double>& MuStar) { m_MuStar = MuStar; }

		// Override functions begin

		/**
		 * \brief Fill the chemical potentials of the particle species.
		 */
		void FillChemicalPotentials();

		/**
		 * \brief Set the chemical potentials of the particle species.
		 *
		 * \param chem Vector of chemical potentials, one element per particle species.
		 *             If empty, the chemical potentials are set to zero.
		 */
		virtual void SetChemicalPotentials(const std::vector<double>& chem = std::vector<double>(0));

		/**
		 * \brief Calculate the primordial densities of the particle species.
		 */
		virtual void CalculatePrimordialDensities();

		/**
		 * \brief Calculate the charge fluctuations of the particle species.
		 *
		 * \param chgs Vector of charges, one element per particle species.
		 * \param order Order of the fluctuations (default is 4).
		 * \return Vector of charge fluctuations, one element per particle species.
		 */
		virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double>& chgs, int order = 4, bool dimensionfull = false);

		/**
		 * \brief Calculate the charge fluctuations of the particle species (old method).
		 *
		 * \param chgs Vector of charges, one element per particle species.
		 * \param order Order of the fluctuations (default is 4).
		 * \return Vector of charge fluctuations, one element per particle species.
		 */
		virtual std::vector<double> CalculateChargeFluctuationsOld(const std::vector<double>& chgs, int order = 4);

		/**
		 * \brief Calculate the fluctuations of the particle species.
		 *
		 * \param order Order of the fluctuations (default is not specified).
		 * \return Vector of fluctuations, one element per particle species.
		 */
		virtual std::vector< std::vector<double> >  CalculateFluctuations(int order);

		/**
		 * \brief Calculate the two-particle correlations of the particle species.
		 */
		void CalculateTwoParticleCorrelations();

		/**
		 * \brief Calculate the fluctuations.
		 */
		void CalculateFluctuations();

		/**
		 * \brief Calculate the pressure of the system.
		 *
		 * \return The pressure of the system.
		 */
		virtual double CalculatePressure();

		/**
		 * \brief Calculate the energy density of the system.
		 *
		 * \return The energy density of the system.
		 */
		virtual double CalculateEnergyDensity();

		/**
		 * \brief Calculate the entropy density of the system.
		 *
		 * \return The entropy density of the system.
		 */
		virtual double CalculateEntropyDensity();

    /**
     * \brief Calculate the derivative of the energy density with respect to temperature.
     *
     * \return The derivative of the energy density with respect to temperature.
     */
    virtual double CalculateEnergyDensityDerivativeT();

    virtual double CalculateEntropyDensityDerivativeTZeroTemperature();

    /**
     * \brief Calculate the temperature derivatives of the system.
     */
    virtual void CalculateTemperatureDerivatives();

		// Dummy
		/**
		 * \brief Calculate the baryon matter entropy density of the system.
		 *
		 * \return The baryon matter entropy density of the system (always returns 0).
		 */
		virtual double CalculateBaryonMatterEntropyDensity() { return 0.; }

		/**
		 * \brief Calculate the meson matter entropy density of the system.
		 *
		 * \return The meson matter entropy density of the system (always returns 0).
		 */
		virtual double CalculateMesonMatterEntropyDensity() { return 0.; }

		/**
		 * \brief Calculate the scalar density of a particle species.
		 *
		 * \param part 0-based particle species index.
		 * \return The scalar density of the particle species.
		 */
		virtual double ParticleScalarDensity(int part);

		/**
		 * \brief Check if the last solution was successful.
		 *
		 * \return True if the last solution was successful, false otherwise.
		 */
		bool   IsLastSolutionOK() const { return m_LastBroydenSuccessFlag && m_LastCalculationSuccessFlag; }

		/**
		 * \brief Get the ideal gas density of a particle species.
		 *
		 * \param index 0-based particle species index.
		 * \return The ideal gas density of the particle species.
		 */
		double DensityId(int index) { return m_DensitiesId[index]; }

		// Override functions end

		/**
		 * \brief Get the component indices of the particle species.
		 *
		 * \return Vector of component indices, one element per particle species.
		 */
		const std::vector< std::vector<int> >& ComponentIndices() const { return m_dMuStarIndices; }

		/**
		 * \brief Get the delta mu (chemical potential shift due to interactions) of a particle species.
		 *
		 * \param i 0-based particle species index.
		 * \return The delta mu of the particle species.
		 */
		virtual double DeltaMu(int i) const { return MuShift(i); }

		/**
		 * \brief Vector of computed susceptibilities values.
		 */
		std::vector< std::vector<double> > m_chi;

		/**
		 * \brief Vector of computed susceptibilities for a specified arbitraty charge.
		 */
		std::vector<double> m_chiarb;
		double ChiArb(int charge);

	protected:

	    /**
	     * \brief Partitions particles species into sets that have identical pair interactions.
	     *
	     * This function is used to optimize the calculation of the real gas model equations.
	     */
    	void CalculateComponentsMap();

		/**
		 * \brief Uses the Broyden method with a provided initial guess
		 *        to determine the shifted chemical potentials
		 *        by solving the transcendental equations with the
		 *        Broyden's method.
		 *
		 * \param muStarInit Initial guess for the shifted chemical potentials
		 * \return std::vector<double> The solved shifted chemical potentials
		 */
		virtual std::vector<double> SearchSingleSolution(const std::vector<double>& muStarInit);

		/**
		 * \brief Uses the Broyden method with a provided initial guess
		 *        to determine the shifted chemical potentials
		 *        by solving the transcendental equations with the
		 *        Broyden's method, using the component mapping.
		 *
		 * \param muStarInit Initial guess for the shifted chemical potentials
		 * \return std::vector<double> The solved shifted chemical potentials
		 */
		virtual std::vector<double> SearchSingleSolutionUsingComponents(const std::vector<double>& muStarInit);

		/**
		 * \brief Uses the Broyden method with a provided initial guess
		 *        to determine the shifted chemical potentials
		 *        by solving the transcendental equations for all particles with the
		 *        Broyden's method.
		 *
		 * \param muStarInit Initial guess for the shifted chemical potentials
		 * \return std::vector<double> The solved shifted chemical potentials
		 */
		virtual std::vector<double> SearchSingleSolutionDirect(const std::vector<double>& muStarInit);

		/**
		 * \brief Uses the Broyden method with different initial guesses
		 *        to look for different possible solutions
		 *        of the transcendental equations for
		 *        shifted chemical potentials
		 *
		 * Looks for the solution with the largest pressure.
		 *
		 * \param iters Number of different initial guesses to try
		 * \return std::vector<double> The solution with the largest pressure among those which were found
		 */
		std::vector<double> SearchMultipleSolutions(int iters = 300);

		/**
		 * \brief Try different initial guesses and return the first valid solution found.
		 *
		 * Used as a fallback when the default single-solution search fails.
		 *
		 * \param iters Number of different initial guesses to try
		 * \return std::vector<double> The first valid solution found
		 */
		std::vector<double> SearchFirstSolution(int iters = 50);

		/// Solve the transcedental equations for the
		/// shifted chemical potentials
		void SolveEquations();

		/**
		 * \brief The shift in the chemical potential
		 *        of particle species i due to the
		 *        QvdW interactions.
		 *
		 * \param i 0-based particle specie index
		 * \return  The shift in the chemical potential
		 */
		virtual double MuShift(int id) const;

		/// Vector of ideal gas densities with shifted chemical potentials
		std::vector<double> m_DensitiesId;

		/// Vector of scalar densities. Not used.
		std::vector<double> m_scaldens;

		/// Whether multiple solutions are considered
		bool   m_SearchMultipleSolutions;

		/// Whether Broyden's method was successfull
		bool   m_LastBroydenSuccessFlag;

		/// Whether the mapping to components with the same VDW parameters has been calculated
		bool   m_ComponentMapCalculated;

		/// Vector of the shifted chemical potentials
		std::vector<double> m_MuStar;

		/// Mapping from particle species to dMuStar indices
		std::vector<int> m_MapTodMuStar;

		/// Mapping from dMuStar indices to particle species
		std::vector<int> m_MapFromdMuStar;

		/// Vector of component indices for each particle species
		std::vector< std::vector<int> > m_dMuStarIndices;


		/// Excluded volume model used in the real gas HRG model
		ExcludedVolumeModelMultiBase* m_exvolmod;

		/// Mean field model used in the real gas HRG model
		MeanFieldModelMultiBase* m_mfmod;

		/// Excluded volume model object in the ideal gas limit
		ExcludedVolumeModelMultiBase* m_exvolmodideal;

		/// Mean field model object in the ideal gas limit
		MeanFieldModelMultiBase* m_mfmodideal;
	private:
		/**
		 * \brief Class implementing the Broyden equations for the real gas model.
		 *
		 * This class is used to solve the transcendental equations for the shifted chemical potentials.
		 */
		class BroydenEquationsRealGas : public BroydenEquations
		{
		public:
			/**
			 * \brief Construct a new BroydenEquationsRealGas object
			 *
			 * \param model Pointer to the ThermalModelRealGas object
			 */
			BroydenEquationsRealGas(ThermalModelRealGas* model) : BroydenEquations(), m_THM(model) { m_N = model->ComponentsNumber(); }

			/**
			 * \brief Evaluate the Broyden equations for the real gas model.
			 *
			 * \param x Vector of shifted chemical potentials
			 * \return Vector of Broyden equations
			 */
			std::vector<double> Equations(const std::vector<double>& x);
		private:
			/// Pointer to the ThermalModelRealGas object
			ThermalModelRealGas* m_THM;
		};

		/**
		 * \brief Class implementing the Broyden Jacobian for the real gas model.
		 *
		 * This class is used to calculate the Jacobian of the Broyden equations.
		 */
		class BroydenJacobianRealGas : public BroydenJacobian
		{
		public:
			/**
			 * \brief Construct a new BroydenJacobianRealGas object
			 *
			 * \param model Pointer to the ThermalModelRealGas object
			 */
			BroydenJacobianRealGas(ThermalModelRealGas* model) : BroydenJacobian(), m_THM(model) { }

			/**
			 * \brief Evaluate the Broyden Jacobian for the real gas model.
			 *
			 * \param x Vector of shifted chemical potentials
			 * \return Vector of Broyden Jacobian
			 */
			std::vector<double> Jacobian(const std::vector<double>& x);
		private:
			/// Pointer to the ThermalModelRealGas object
			ThermalModelRealGas* m_THM;
		};

		/**
		 * \brief Class implementing the Broyden solution criterium for the real gas model.
		 *
		 * This class is used to determine whether the Broyden method has converged.
		 */
		class BroydenSolutionCriteriumRealGas : public Broyden::BroydenSolutionCriterium
		{
		public:
			/**
			 * \brief Construct a new BroydenSolutionCriteriumRealGas object
			 *
			 * \param model Pointer to the ThermalModelRealGas object
			 * \param relative_error Relative error tolerance (default is Broyden::TOL)
			 */
			BroydenSolutionCriteriumRealGas(ThermalModelRealGas* model, double relative_error = Broyden::TOL) : Broyden::BroydenSolutionCriterium(relative_error), m_THM(model) { }

			/**
			 * \brief Check whether the Broyden method has converged.
			 *
			 * \param x Vector of shifted chemical potentials
			 * \param f Vector of Broyden equations
			 * \param xdelta Vector of changes in the shifted chemical potentials (optional)
			 * \return True if the Broyden method has converged, false otherwise
			 */
			virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
		protected:
			/// Pointer to the ThermalModelRealGas object
			ThermalModelRealGas* m_THM;
		};

		/**
		 * \brief Class implementing the Broyden equations for the real gas model using components.
		 *
		 * This class is used to solve the transcendental equations for the shifted chemical potentials using the component mapping.
		 */
		class BroydenEquationsRealGasComponents : public BroydenEquations
		{
		public:
			/**
			 * \brief Construct a new BroydenEquationsRealGasComponents object
			 *
			 * \param model Pointer to the ThermalModelRealGas object
			 */
			BroydenEquationsRealGasComponents(ThermalModelRealGas* model) : BroydenEquations(), m_THM(model) { m_N = model->m_MapFromdMuStar.size(); }

			/**
			 * \brief Evaluate the Broyden equations for the real gas model using components.
			 *
			 * \param x Vector of shifted chemical potentials
			 * \return Vector of Broyden equations
			 */
			std::vector<double> Equations(const std::vector<double>& x);
		private:
			/// Pointer to the ThermalModelRealGas object
			ThermalModelRealGas* m_THM;
		};

		/**
		 * \brief Class implementing the Broyden Jacobian for the real gas model using components.
		 *
		 * This class is used to calculate the Jacobian of the Broyden equations using the component mapping.
		 */
		class BroydenJacobianRealGasComponents : public BroydenJacobian
		{
		public:
			/**
			 * \brief Construct a new BroydenJacobianRealGasComponents object
			 *
			 * \param model Pointer to the ThermalModelRealGas object
			 */
			BroydenJacobianRealGasComponents(ThermalModelRealGas* model) : BroydenJacobian(), m_THM(model) { }

			/**
			 * \brief Evaluate the Broyden Jacobian for the real gas model using components.
			 *
			 * \param x Vector of shifted chemical potentials
			 * \return Vector of Broyden Jacobian
			 */
			std::vector<double> Jacobian(const std::vector<double>& x);
		private:
			/// Pointer to the ThermalModelRealGas object
			ThermalModelRealGas* m_THM;
		};

	};

} // namespace thermalfist

#endif
