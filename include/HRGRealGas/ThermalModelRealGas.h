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

		void SetExcludedVolumeModel(ExcludedVolumeModelMultiBase* exvolmod) { m_exvolmod = exvolmod; }
		void SetMeanFieldModel(MeanFieldModelMultiBase* mfmod) { m_mfmod = mfmod; }

		ExcludedVolumeModelMultiBase* ExcludedVolumeModel() const { return m_exvolmod; }
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
		bool UseMultipleSolutionsMode() const { return m_SearchMultipleSolutions; }

		/// The shifted chemical potential of particle species i.
		double MuStar(int i) const { return m_MuStar[i]; }

		/// Returns vector of shifted chemical potentials,
		/// one element per each species
		std::vector<double> GetMuStar() const { return m_MuStar; }

		/// Set the vector of shifted chemical potentials
		void SetMuStar(const std::vector<double>& MuStar) { m_MuStar = MuStar; }

		// Override functions begin

		void FillChemicalPotentials();

		virtual void SetChemicalPotentials(const std::vector<double>& chem = std::vector<double>(0));

		virtual void CalculatePrimordialDensities();

		virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double>& chgs, int order = 4);
		virtual std::vector<double> CalculateChargeFluctuationsOld(const std::vector<double>& chgs, int order = 4);

		virtual std::vector< std::vector<double> >  CalculateFluctuations(int order);

		void CalculateTwoParticleCorrelations();

		void CalculateFluctuations();

		virtual double CalculatePressure();

		virtual double CalculateEnergyDensity();

		virtual double CalculateEntropyDensity();

		// Dummy
		virtual double CalculateBaryonMatterEntropyDensity() { return 0.; }

		virtual double CalculateMesonMatterEntropyDensity() { return 0.; }

		virtual double ParticleScalarDensity(int part);

		bool   IsLastSolutionOK() const { return m_LastBroydenSuccessFlag && m_LastCalculationSuccessFlag; }

		double DensityId(int index) { return m_DensitiesId[index]; }

		// Override functions end

		virtual double DeltaMu(int i) const { return MuShift(i); }

		std::vector< std::vector<double> > m_chi;

		std::vector<double> m_chiarb;

	protected:

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

		/// Vector of the shifted chemical potentials
		std::vector<double> m_MuStar;


		ExcludedVolumeModelMultiBase* m_exvolmod;
		MeanFieldModelMultiBase* m_mfmod;

		ExcludedVolumeModelMultiBase* m_exvolmodideal;
		MeanFieldModelMultiBase* m_mfmodideal;
	private:
		class BroydenEquationsRealGas : public BroydenEquations
		{
		public:
			BroydenEquationsRealGas(ThermalModelRealGas* model) : BroydenEquations(), m_THM(model) { m_N = model->ComponentsNumber(); }
			std::vector<double> Equations(const std::vector<double>& x);
		private:
			ThermalModelRealGas* m_THM;
		};

		class BroydenJacobianRealGas : public BroydenJacobian
		{
		public:
			BroydenJacobianRealGas(ThermalModelRealGas* model) : BroydenJacobian(), m_THM(model) { }
			std::vector<double> Jacobian(const std::vector<double>& x);
		private:
			ThermalModelRealGas* m_THM;
		};

		class BroydenSolutionCriteriumRealGas : public Broyden::BroydenSolutionCriterium
		{
		public:
			BroydenSolutionCriteriumRealGas(ThermalModelRealGas* model, double relative_error = Broyden::TOL) : Broyden::BroydenSolutionCriterium(relative_error), m_THM(model) { }
			virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
		protected:
			ThermalModelRealGas* m_THM;
		};
	};

} // namespace thermalfist

#endif

