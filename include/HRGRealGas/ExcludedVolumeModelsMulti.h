/*
 * Thermal-FIST package
 *
 * Copyright (c) 2017-2022 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EXCLUDEVOLUMEMODELSMULTI_H
#define EXCLUDEVOLUMEMODELSMULTI_H
#include <cmath>
#include <vector>

#include "ExcludedVolumeModels.h"
#include "HRGBase/Broyden.h"

namespace thermalfist {

	/**
	 * \brief Base class for multi-component excluded volume models.
	 * 
	 * This class serves as the base for all multi-component excluded volume models,
	 * providing the interface for calculating the excluded volume effects in a system
	 * with multiple particle species.
	 */
	class ExcludedVolumeModelMultiBase {
	public:
		/**
		 * \brief Constructor for the ExcludedVolumeModelMultiBase class.
		 * 
		 * \param N Number of particle species in the system.
		 */
		ExcludedVolumeModelMultiBase(int N) : 
			m_N(N), 
			m_densities(std::vector<double>(m_N,0.)), 
			m_components(std::vector<int>(m_N, 0)) 
		{
			ComputeComponents();
		}
		virtual ~ExcludedVolumeModelMultiBase() { }
		
		/**
		 * \brief Calculates the suppression factor for species i.
		 * 
		 * \param i Index of the particle species.
		 * \return Suppression factor for species i.
		 */
		virtual double f(int i) const { return 1.; }
		
		/**
		 * \brief Calculates the first derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return First derivative of the suppression factor.
		 */
		virtual double df(int i, int j) const { return 0.; }
		
		/**
		 * \brief Calculates the second derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Second derivative of the suppression factor.
		 */
		virtual double d2f(int i, int j, int k) const { return 0.; }
		
		/**
		 * \brief Calculates the third derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Third derivative of the suppression factor.
		 */
		virtual double d3f(int i, int j, int k, int l) const { return 0.; }
		
		/**
		 * \brief Calculates the fourth derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \param m Fifth species index.
		 * \return Fourth derivative of the suppression factor.
		 */
		virtual double d4f(int i, int j, int k, int l, int m) const { return 0.; }
		
		/**
		 * \brief Calculates the temperature derivative of the suppression factor.
		 * 
		 * \param i Species index.
		 * \return Temperature derivative of the suppression factor.
		 */
		virtual double dfdT(int i) const { return 0.; }
		
		/**
		 * \brief Solves for the actual densities given the ideal gas densities.
		 * 
		 * \param ntil Vector of ideal gas densities.
		 * \return Vector of actual densities.
		 */
		virtual std::vector<double> nsol(const std::vector<double>& ntil) { return ntil; }
		
		/**
		 * \brief Solves for the actual densities using Broyden's method.
		 * 
		 * \param ntil Vector of ideal gas densities.
		 * \return Vector of actual densities.
		 */
		virtual std::vector<double> nsolBroyden(const std::vector<double>& ntil);
		
		/**
		 * \brief Solves for the actual densities using Broyden's method, considering components.
		 * 
		 * \param ntil Vector of ideal gas densities.
		 * \return Vector of actual densities.
		 */
		virtual std::vector<double> nsolBroydenComponents(const std::vector<double>& ntil);
		
		/**
		 * \brief Sets the densities of particle species.
		 * 
		 * \param n Vector of densities.
		 */
		virtual void SetDensities(const std::vector<double>& n) { m_densities = n; }
		
		/**
		 * \brief Gets the component indices.
		 * 
		 * \return Vector of component indices.
		 */
		virtual const std::vector<int>& ComponentIndices() const { return m_components; }
		
		/**
		 * \brief Gets the component indices from.
		 * 
		 * \return Vector of component indices from.
		 */
		virtual const std::vector<int>& ComponentIndicesFrom() const { return m_componentsFrom; }
		
		/**
		 * \brief Gets the number of components.
		 * 
		 * \return Number of components.
		 */
		virtual const int ComponentsNumber() const { return m_componentsNumber; }
	protected:
		/**
		 * \brief Computes the components based on the excluded volume parameters.
		 */
		virtual void ComputeComponents();
		int m_N;
		std::vector<double> m_densities;
		std::vector<int> m_components;
		std::vector<int> m_componentsFrom;
		int m_componentsNumber;

	  class BroydenEquationsEVMulti : public BroydenEquations
		{
		public:
			/**
			 * \brief Constructor for the BroydenEquationsEVMulti class.
			 * 
			 * \param evmodel Pointer to the excluded volume model.
			 * \param ntil Pointer to the ideal gas densities.
			 * \param componentsMode Flag indicating whether to consider components.
			 */
			BroydenEquationsEVMulti(ExcludedVolumeModelMultiBase* evmodel, const std::vector<double>* ntil, bool componentsMode = false) : 
				BroydenEquations(),
				m_EVM(evmodel),
				m_ntil(ntil),
				m_componentsMode(componentsMode)
			{ 
				m_N = evmodel->m_N; 
				if (m_componentsMode)
					m_N = evmodel->m_componentsNumber;
			}

			/**
			 * \brief Calculates the equations for Broyden's method.
			 * 
			 * \param x Vector of variables.
			 * \return Vector of equations.
			 */
			std::vector<double> Equations(const std::vector<double>& x);
		private:
			ExcludedVolumeModelMultiBase* m_EVM;
			const std::vector<double>* m_ntil;
			bool m_componentsMode;
		};

		class BroydenJacobianEVMulti : public BroydenJacobian
		{
		public:
			/**
			 * \brief Constructor for the BroydenJacobianEVMulti class.
			 * 
			 * \param evmodel Pointer to the excluded volume model.
			 * \param ntil Pointer to the ideal gas densities.
			 * \param componentsMode Flag indicating whether to consider components.
			 */
			BroydenJacobianEVMulti(ExcludedVolumeModelMultiBase* evmodel, 
				const std::vector<double>* ntil,
				bool componentsMode = false) :
				BroydenJacobian(), 
				m_EVM(evmodel),
				m_ntil(ntil),
				m_componentsMode(componentsMode)
			{
			}

			/**
			 * \brief Calculates the Jacobian matrix for Broyden's method.
			 * 
			 * \param x Vector of variables.
			 * \return Jacobian matrix.
			 */
			std::vector<double> Jacobian(const std::vector<double>& x);
		private:
			ExcludedVolumeModelMultiBase* m_EVM;
			const std::vector<double>* m_ntil;
			bool m_componentsMode;
		};


	};

	/**
	 * \brief Implementation of the diagonal van der Waals excluded volume model.
	 * 
	 * This class implements the diagonal van der Waals excluded volume model,
	 * where only the self-interactions of particles are considered.
	 */
	class ExcludedVolumeModelDiagonalVDW :
		public ExcludedVolumeModelMultiBase {
	public:
		/**
		 * \brief Constructor for the ExcludedVolumeModelDiagonalVDW class.
		 * 
		 * \param b Vector of excluded volumes for each species.
		 * \param dbdT Vector of temperature derivatives of excluded volumes (optional).
		 */
		ExcludedVolumeModelDiagonalVDW(const std::vector<double>& b, const std::vector<double>& dbdT = std::vector<double>())
			: ExcludedVolumeModelMultiBase(b.size()), m_b(b), m_dbdT(dbdT) {
			if (m_dbdT.size() == 0)
				m_dbdT = std::vector<double>(m_b.size(), 0.);
			ComputeComponents();
		}
		/**
		 * \brief Calculates the suppression factor for species i.
		 * 
		 * \param i Index of the particle species.
		 * \return Suppression factor for species i.
		 */
		virtual double f(int i) const;
		/**
		 * \brief Calculates the first derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return First derivative of the suppression factor.
		 */
		virtual double df(int i, int j)											const;
		/**
		 * \brief Calculates the second derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Second derivative of the suppression factor.
		 */
		virtual double d2f(int i, int j, int k)							const { return 0.; }
		/**
		 * \brief Calculates the third derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Third derivative of the suppression factor.
		 */
		virtual double d3f(int i, int j, int k, int l)				const { return 0.; }
		/**
		 * \brief Calculates the fourth derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \param m Fifth species index.
		 * \return Fourth derivative of the suppression factor.
		 */
		virtual double d4f(int i, int j, int k, int l, int m) const { return 0.; }
		/**
		 * \brief Calculates the temperature derivative of the suppression factor.
		 * 
		 * \param i Species index.
		 * \return Temperature derivative of the suppression factor.
		 */
		virtual double dfdT(int i)											      const;
		/**
		 * \brief Solves for the actual densities given the ideal gas densities.
		 * 
		 * \param nid Vector of ideal gas densities.
		 * \return Vector of actual densities.
		 */
		virtual std::vector<double> nsol(const std::vector<double>& nid);
	protected:
		/**
		 * \brief Computes the components based on the excluded volume parameters.
		 */
		virtual void ComputeComponents();
		std::vector<double> m_b;
		std::vector<double> m_dbdT;
	};

	/**
	 * \brief Implementation of a diagonal generalized excluded volume model.
	 * 
	 * This class implements a diagonal generalized excluded volume model,
	 * where the excluded volume effects are described by a general model
	 * specified by the ExcludedVolumeModelBase object.
	 * The object pointed by ExcludedVolumeModelBase will be deleted on destruction.
	 */
	class ExcludedVolumeModelDiagonalGeneralized :
		public ExcludedVolumeModelMultiBase {
	public:
		/**
		 * \brief Constructor for the ExcludedVolumeModelDiagonalGeneralized class.
		 * 
		 * \param evmodelsingle Pointer to the excluded volume model for a single component.
		 * \param b Vector of excluded volumes for each species.
		 * \param dbdT Vector of temperature derivatives of excluded volumes (optional).
		 */
		ExcludedVolumeModelDiagonalGeneralized(
			ExcludedVolumeModelBase *evmodelsingle,
			const std::vector<double>& b, 
			const std::vector<double>& dbdT = std::vector<double>())
			: ExcludedVolumeModelMultiBase(b.size()), m_evmodelsingle(evmodelsingle), m_b(b), m_dbdT(dbdT), m_eta(0.) {
			if (m_dbdT.size() == 0)
				m_dbdT = std::vector<double>(m_b.size(), 0.);
			ComputeComponents();
		}
		/**
		 * \brief Destructor for the ExcludedVolumeModelDiagonalGeneralized class.
		 */
		virtual ~ExcludedVolumeModelDiagonalGeneralized();
		/**
		 * \brief Calculates the suppression factor for species i.
		 * 
		 * \param i Index of the particle species.
		 * \return Suppression factor for species i.
		 */
		virtual double f(int i) const;
		/**
		 * \brief Calculates the first derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return First derivative of the suppression factor.
		 */
		virtual double df(int i, int j)											const;
		/**
		 * \brief Calculates the second derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Second derivative of the suppression factor.
		 */
		virtual double d2f(int i, int j, int k)							const;
		/**
		 * \brief Calculates the third derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Third derivative of the suppression factor.
		 */
		virtual double d3f(int i, int j, int k, int l)				const;
		/**
		 * \brief Calculates the fourth derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \param m Fifth species index.
		 * \return Fourth derivative of the suppression factor.
		 */
		virtual double d4f(int i, int j, int k, int l, int m) const;
		/**
		 * \brief Calculates the temperature derivative of the suppression factor.
		 * 
		 * \param i Species index.
		 * \return Temperature derivative of the suppression factor.
		 */
		virtual double dfdT(int i)											      const;
		/**
		 * \brief Solves for the actual densities given the ideal gas densities.
		 * 
		 * \param nid Vector of ideal gas densities.
		 * \return Vector of actual densities.
		 */
		virtual std::vector<double> nsol(const std::vector<double>& nid);
		/**
		 * \brief Sets the densities of particle species.
		 * 
		 * \param n Vector of densities.
		 */
		virtual void SetDensities(const std::vector<double>& n);
	protected:
		/**
		 * \brief Computes the components based on the excluded volume parameters.
		 */
		virtual void ComputeComponents();
		/**
		 * \brief Calculates the eta parameter for the given densities.
		 * 
		 * \param n Vector of densities.
		 * \return Eta parameter.
		 */
		double GetEta(const std::vector<double>& n) const;
		ExcludedVolumeModelBase* m_evmodelsingle;
		std::vector<double> m_b;
		std::vector<double> m_dbdT;
		double m_eta;
	};


	/**
	 * \brief Implementation of the crossterms van der Waals excluded volume model.
	 * 
	 * This class implements the crossterms van der Waals excluded volume model,
	 * where interactions between different particle species are considered.
	 */
	class ExcludedVolumeModelCrosstermsVDW :
		public ExcludedVolumeModelMultiBase {
	public:
		/**
		 * \brief Constructor for the ExcludedVolumeModelCrosstermsVDW class.
		 * 
		 * \param b Matrix of excluded volumes between species.
		 * \param dbdT Matrix of temperature derivatives of excluded volumes (optional).
		 */
		ExcludedVolumeModelCrosstermsVDW(const std::vector< std::vector<double> >& b, const std::vector< std::vector<double> >& dbdT = std::vector< std::vector<double> >())
			: ExcludedVolumeModelMultiBase(b.size()), m_b(b), m_dbdT(dbdT) {
			if (m_dbdT.size() == 0)
				m_dbdT = std::vector< std::vector<double> >(m_b.size(), std::vector<double>(m_b.size(), 0.));
			ComputeComponents();
		}
		/**
		 * \brief Calculates the suppression factor for species i.
		 * 
		 * \param i Index of the particle species.
		 * \return Suppression factor for species i.
		 */
		virtual double f(int i) const;
		/**
		 * \brief Calculates the first derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return First derivative of the suppression factor.
		 */
		virtual double df(int i, int j)											  const;
		/**
		 * \brief Calculates the second derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Second derivative of the suppression factor.
		 */
		virtual double d2f(int i, int j, int k)								const { return 0.; }
		/**
		 * \brief Calculates the third derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Third derivative of the suppression factor.
		 */
		virtual double d3f(int i, int j, int k, int l)				const { return 0.; }
		/**
		 * \brief Calculates the fourth derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \param m Fifth species index.
		 * \return Fourth derivative of the suppression factor.
		 */
		virtual double d4f(int i, int j, int k, int l, int m) const { return 0.; }
		/**
		 * \brief Calculates the temperature derivative of the suppression factor.
		 * 
		 * \param i Species index.
		 * \return Temperature derivative of the suppression factor.
		 */
		virtual double dfdT(int i)											      const;
		/**
		 * \brief Solves for the actual densities given the ideal gas densities.
		 * 
		 * \param nid Vector of ideal gas densities.
		 * \return Vector of actual densities.
		 */
		virtual std::vector<double> nsol(const std::vector<double>& nid);
	protected:
		/**
		 * \brief Computes the components based on the excluded volume parameters.
		 */
		virtual void ComputeComponents();
		std::vector< std::vector<double> > m_b;
		std::vector< std::vector<double> > m_dbdT;
	};

	/**
	 * \brief Implementation of a crossterms generalized excluded volume model.
	 * 
	 * This class implements a crossterms generalized excluded volume model,
	 * where the excluded volume effects between different particle species
	 * are described by a general model specified by the ExcludedVolumeModelBase object.
	 * The object pointed by ExcludedVolumeModelBase will be deleted on destruction.
	 */
	class ExcludedVolumeModelCrosstermsGeneralized :
		public ExcludedVolumeModelMultiBase {
	public:
		/**
		 * \brief Constructor for the ExcludedVolumeModelCrosstermsGeneralized class.
		 * 
		 * \param evmodelsingle Pointer to the excluded volume model for a single component.
		 * \param b Matrix of excluded volumes between species.
		 * \param dbdT Matrix of temperature derivatives of excluded volumes (optional).
		 */
		ExcludedVolumeModelCrosstermsGeneralized(
			ExcludedVolumeModelBase* evmodelsingle, 
			const std::vector< std::vector<double> >& b, 
			const std::vector< std::vector<double> >& dbdT = std::vector< std::vector<double> >()
		)
			: ExcludedVolumeModelMultiBase(b.size()), m_evmodelsingle(evmodelsingle), m_b(b), m_dbdT(dbdT), m_componentsDisconnected(false) {
			if (m_dbdT.size() == 0)
				m_dbdT = std::vector< std::vector<double> >(m_b.size(), std::vector<double>(m_b.size(), 0.));
			ComputeComponents();
			m_etas = std::vector<double>(m_componentsNumber, 0.);
		}
		/**
		 * \brief Destructor for the ExcludedVolumeModelCrosstermsGeneralized class.
		 */
		virtual ~ExcludedVolumeModelCrosstermsGeneralized();
		/**
		 * \brief Calculates the suppression factor for species i.
		 * 
		 * \param i Index of the particle species.
		 * \return Suppression factor for species i.
		 */
		virtual double f(int i) const;
		/**
		 * \brief Calculates the first derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return First derivative of the suppression factor.
		 */
		virtual double df(int i, int j)											  const;
		/**
		 * \brief Calculates the second derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Second derivative of the suppression factor.
		 */
		virtual double d2f(int i, int j, int k)								const;
		/**
		 * \brief Calculates the third derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Third derivative of the suppression factor.
		 */
		virtual double d3f(int i, int j, int k, int l)				const;
		/**
		 * \brief Calculates the fourth derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \param m Fifth species index.
		 * \return Fourth derivative of the suppression factor.
		 */
		virtual double d4f(int i, int j, int k, int l, int m) const;
		/**
		 * \brief Calculates the temperature derivative of the suppression factor.
		 * 
		 * \param i Species index.
		 * \return Temperature derivative of the suppression factor.
		 */
		virtual double dfdT(int i)											      const;
		/**
		 * \brief Solves for the actual densities given the ideal gas densities.
		 * 
		 * \param nid Vector of ideal gas densities.
		 * \return Vector of actual densities.
		 */
		virtual std::vector<double> nsol(const std::vector<double>& nid);
		/**
		 * \brief Sets the densities of particle species.
		 * 
		 * \param n Vector of densities.
		 */
		virtual void SetDensities(const std::vector<double>& n);
	protected:
		/**
		 * \brief Computes the components based on the excluded volume parameters.
		 */
		virtual void ComputeComponents();
		/**
		 * \brief Calculates the eta parameter for the given densities and species index.
		 * 
		 * \param i Species index.
		 * \param n Vector of densities.
		 * \return Eta parameter.
		 */
		double GetEta(int i, const std::vector<double>& n) const;
		ExcludedVolumeModelBase* m_evmodelsingle;
		std::vector< std::vector<double> > m_b;
		std::vector< std::vector<double> > m_dbdT;
		std::vector<double> m_etas;
		bool m_componentsDisconnected;
	};

	/**
	 * \brief Implementation of an excluded volume model with components.
	 * 
	 * This class implements an excluded volume model where particles are
	 * grouped into components with similar excluded volume properties.
	 */
	class ExcludedVolumeModelComponents :
		public ExcludedVolumeModelMultiBase {
	public:
		/**
		 * \brief Constructor for the ExcludedVolumeModelComponents class.
		 * 
		 * \param components Number of components.
		 * \param evmods Vector of excluded volume models for each component.
		 * \param ind Vector of component indices for each particle species.
		 * \param b Vector of excluded volumes for each species.
		 * \param dbdT Vector of temperature derivatives of excluded volumes (optional).
		 */
		ExcludedVolumeModelComponents(
			int components, 
			const std::vector<ExcludedVolumeModelBase*>& evmods, 
			const std::vector<int>& ind, 
			const std::vector<double>& b, 
			const std::vector<double>& dbdT = std::vector<double>()
		)
			: ExcludedVolumeModelMultiBase(ind.size()), 
			//m_N_components(components), 
			m_evmodels(evmods),
			//m_indices(ind), 
			m_b(b), 
			m_dbdT(dbdT), 
			m_densities_components(std::vector<double>(b.size(), 0.))
		{
			if (m_dbdT.size() == 0)
				m_dbdT.resize(m_b.size(), 0.);
			m_components = ind;
			m_componentsNumber = components;
			ComputeComponents();
		}
		/**
		 * \brief Destructor for the ExcludedVolumeModelComponents class.
		 */
		virtual ~ExcludedVolumeModelComponents();
		/**
		 * \brief Calculates the suppression factor for species i.
		 * 
		 * \param i Index of the particle species.
		 * \return Suppression factor for species i.
		 */
		virtual double f(int i) const;
		/**
		 * \brief Calculates the first derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return First derivative of the suppression factor.
		 */
		virtual double df(int i, int j) const;
		/**
		 * \brief Calculates the second derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Second derivative of the suppression factor.
		 */
		virtual double d2f(int i, int j, int k) const;
		/**
		 * \brief Calculates the third derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Third derivative of the suppression factor.
		 */
		virtual double d3f(int i, int j, int k, int l) const;
		/**
		 * \brief Calculates the fourth derivative of the suppression factor.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \param m Fifth species index.
		 * \return Fourth derivative of the suppression factor.
		 */
		virtual double d4f(int i, int j, int k, int l, int m) const;
		/**
		 * \brief Calculates the temperature derivative of the suppression factor.
		 * 
		 * \param i Species index.
		 * \return Temperature derivative of the suppression factor.
		 */
		virtual double dfdT(int i)											      const;
		/**
		 * \brief Solves for the actual densities given the ideal gas densities.
		 * 
		 * \param nid Vector of ideal gas densities.
		 * \return Vector of actual densities.
		 */
		virtual std::vector<double> nsol(const std::vector<double>& nid);
		/**
		 * \brief Sets the densities of particle species.
		 * 
		 * \param n Vector of densities.
		 */
		virtual void SetDensities(const std::vector<double>& n);
	protected:
		/**
		 * \brief Computes the components based on the excluded volume parameters.
		 */
		virtual void ComputeComponents();
		//int m_N_components;
		std::vector<ExcludedVolumeModelBase*> m_evmodels;
		//std::vector<int> m_indices;
		std::vector<double> m_b;
		std::vector<double> m_dbdT;
		std::vector<double> m_densities_components;
	};

} // namespace thermalfist

#endif // REALGASMODELS_H
