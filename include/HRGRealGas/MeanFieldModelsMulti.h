/*
 * Thermal-FIST package
 *
 * Copyright (c) 2017-2022 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef MEANFIELDMODELSMULTI_H
#define MEANFIELDMODELSMULTI_H
#include <cmath>
#include <vector>

#include "MeanFieldModels.h"

namespace thermalfist {

	/**
	 * \brief Base class for multi-component mean field models.
	 * 
	 * This class serves as the base for all multi-component mean field models,
	 * providing the interface for calculating the mean field effects in a system
	 * with multiple particle species. By default, it implements the ideal gas case
	 * with no mean field.
	 */
	class MeanFieldModelMultiBase {
	public:
		/**
		 * \brief Constructor for the MeanFieldModelMultiBase class.
		 * 
		 * \param N Number of particle species in the system.
		 */
		MeanFieldModelMultiBase(int N) : m_N(N) { MeanFieldModelMultiBase::ComputeComponents(); }
		
		/**
		 * \brief Destructor for the MeanFieldModelMultiBase class.
		 */
		virtual ~MeanFieldModelMultiBase() { }
		
		/**
		 * \brief Calculates the mean field value.
		 * 
		 * \return Mean field value in units of GeV/fm^3.
		 */
		virtual double v()	const { return 0.; }
		
		/**
		 * \brief Calculates the first derivative of the mean field.
		 * 
		 * \param i Species index.
		 * \return First derivative of the mean field.
		 */
		virtual double dv(int i)											  const { return 0.; }
		
		/**
		 * \brief Calculates the second derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return Second derivative of the mean field.
		 */
		virtual double d2v(int i, int j)								const { return 0.; }
		
		/**
		 * \brief Calculates the third derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Third derivative of the mean field.
		 */
		virtual double d3v(int i, int j, int k)				  const { return 0.; }
		
		/**
		 * \brief Calculates the fourth derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Fourth derivative of the mean field.
		 */
		virtual double d4v(int i, int j, int k, int l)  const { return 0.; }
		
		/**
		 * \brief Calculates the temperature derivative of the mean field.
		 * 
		 * \return Temperature derivative of the mean field.
		 */
		virtual double dvdT()										        const { return 0.; }
		
		/**
		 * \brief Sets the densities of particle species.
		 * 
		 * \param n Vector of densities.
		 */
		virtual void SetDensities(const std::vector<double>& n);
		
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
		 * \brief Computes the components based on the mean field parameters.
		 */
		virtual void ComputeComponents();
		int m_N;
		std::vector<double> m_densities;
		std::vector<int> m_components;
		std::vector<int> m_componentsFrom;
		int m_componentsNumber;
		std::vector<double> m_densities_components;
	};

	/**
	 * \brief Implementation of the van der Waals mean field model for multiple components.
	 * 
	 * This class implements the van der Waals mean field model for a system
	 * with multiple particle species, where interactions between different
	 * particle species are considered.
	 */
	class MeanFieldModelMultiVDW :
		public MeanFieldModelMultiBase {
	public:
		/**
		 * \brief Constructor for the MeanFieldModelMultiVDW class.
		 * 
		 * \param a Matrix of attraction parameters between species.
		 * \param dadT Matrix of temperature derivatives of attraction parameters (optional).
		 */
		MeanFieldModelMultiVDW(
			const std::vector< std::vector<double> >& a, 
			const std::vector< std::vector<double> >& dadT = std::vector< std::vector<double> >()
		)
			: MeanFieldModelMultiBase(a.size()), 
			m_a(a),
			m_dadT(dadT)
		{
			if (m_dadT.size() == 0)
				m_dadT = std::vector< std::vector<double> >(m_a.size(), std::vector<double>(m_a.size(), 0.));
			ComputeComponents();
		}
		
		/**
		 * \brief Destructor for the MeanFieldModelMultiVDW class.
		 */
		virtual ~MeanFieldModelMultiVDW() { }

		/**
		 * \brief Sets the attraction parameters between species.
		 * 
		 * \param a Matrix of attraction parameters between species.
		 */
		void setAij(const std::vector< std::vector<double> >& a) { m_a = a; }
		
		/**
		 * \brief Sets the temperature derivatives of the attraction parameters.
		 * 
		 * \param dadT Matrix of temperature derivatives of attraction parameters (optional).
		 */
		void setdAijdT(const std::vector< std::vector<double> >& dadT) { m_dadT = dadT; }
		
		/**
		 * \brief Calculates the mean field value.
		 * 
		 * \return Mean field value in units of GeV/fm^3.
		 */
		virtual double v() const;
		
		/**
		 * \brief Calculates the first derivative of the mean field.
		 * 
		 * \param i Species index.
		 * \return First derivative of the mean field.
		 */
		virtual double dv(int i) const;
		
		/**
		 * \brief Calculates the second derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return Second derivative of the mean field.
		 */
		virtual double d2v(int i, int j) const;
		
		/**
		 * \brief Calculates the third derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Third derivative of the mean field.
		 */
		virtual double d3v(int i, int j, int k) const { return 0.; }
		
		/**
		 * \brief Calculates the fourth derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Fourth derivative of the mean field.
		 */
		virtual double d4v(int i, int j, int k, int l) const { return 0.; }
		
		/**
		 * \brief Calculates the temperature derivative of the mean field.
		 * 
		 * \return Temperature derivative of the mean field.
		 */
		virtual double dvdT() const;
	protected:
		/**
		 * \brief Computes the components based on the mean field parameters.
		 */
		virtual void ComputeComponents();
		std::vector< std::vector<double> > m_a;
		std::vector< std::vector<double> > m_dadT;
	};

	/**
	 * \brief Implementation of a mean field model with components.
	 * 
	 * This class implements a mean field model where particles are
	 * grouped into components with similar mean field properties.
	 */
	class MeanFieldModelComponents :
		public MeanFieldModelMultiBase {
	public:
		/**
		 * \brief Constructor for the MeanFieldModelComponents class.
		 * 
		 * \param components Number of components.
		 * \param mfmods Vector of mean field models for each component.
		 * \param ind Vector of component indices for each particle species.
		 */
		MeanFieldModelComponents(
					int components,
					const std::vector<MeanFieldModelBase*>& mfmods,
					const std::vector<int>& ind
				) 
					: MeanFieldModelMultiBase(ind.size()), 
					m_mfmodels(mfmods)
				{
					m_components = ind;
					m_componentsNumber = components;
					ComputeComponents();
				}
		
		/**
		 * \brief Destructor for the MeanFieldModelComponents class.
		 */
		virtual ~MeanFieldModelComponents();
		
		/**
		 * \brief Calculates the mean field value.
		 * 
		 * \return Mean field value in units of GeV/fm^3.
		 */
		virtual double v() const;
		
		/**
		 * \brief Calculates the first derivative of the mean field.
		 * 
		 * \param i Species index.
		 * \return First derivative of the mean field.
		 */
		virtual double dv(int i) const;
		
		/**
		 * \brief Calculates the second derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return Second derivative of the mean field.
		 */
		virtual double d2v(int i, int j) const;
		
		/**
		 * \brief Calculates the third derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Third derivative of the mean field.
		 */
		virtual double d3v(int i, int j, int k) const;
		
		/**
		 * \brief Calculates the fourth derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Fourth derivative of the mean field.
		 */
		virtual double d4v(int i, int j, int k, int l) const;
		
		/**
		 * \brief Calculates the temperature derivative of the mean field.
		 * 
		 * \return Temperature derivative of the mean field.
		 */
		virtual double dvdT() const;
		//virtual void SetDensities(const std::vector<double>& n);
		protected:
		/**
		 * \brief Computes the components based on the mean field parameters.
		 */
		virtual void ComputeComponents();
		//int m_N_components;
		std::vector<MeanFieldModelBase*> m_mfmodels;
		//std::vector<int> m_indices;
	};

	/**
	 * \brief Implementation of a charge density dependent mean field model.
	 * 
	 * This class implements a mean field model where the mean field
	 * depends on the charge density of the system rather than the
	 * individual particle densities.
	 */
	class MeanFieldModelChargeDensityDependent :
		public MeanFieldModelMultiBase {
	public:
		/**
		 * \brief Constructor for the MeanFieldModelChargeDensityDependent class.
		 * 
		 * \param mfmod Pointer to the mean field model for a single component.
		 * \param chgs Vector of charge indices for each particle species.
		 */
		MeanFieldModelChargeDensityDependent(
				  MeanFieldModelBase* mfmod,
					const std::vector<double>& chgs
				)
					: MeanFieldModelMultiBase(chgs.size()),
					m_mfmodel(mfmod),
					m_charges(chgs),
					m_nB(0.)
				{
					ComputeComponents();
				}
		
		/**
		 * \brief Destructor for the MeanFieldModelChargeDensityDependent class.
		 */
		virtual ~MeanFieldModelChargeDensityDependent();
		
		/**
		 * \brief Calculates the mean field value.
		 * 
		 * \return Mean field value in units of GeV/fm^3.
		 */
		virtual double v() const;
		
		/**
		 * \brief Calculates the first derivative of the mean field.
		 * 
		 * \param i Species index.
		 * \return First derivative of the mean field.
		 */
		virtual double dv(int i) const;
		
		/**
		 * \brief Calculates the second derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \return Second derivative of the mean field.
		 */
		virtual double d2v(int i, int j) const;
		
		/**
		 * \brief Calculates the third derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \return Third derivative of the mean field.
		 */
		virtual double d3v(int i, int j, int k) const;
		
		/**
		 * \brief Calculates the fourth derivative of the mean field.
		 * 
		 * \param i First species index.
		 * \param j Second species index.
		 * \param k Third species index.
		 * \param l Fourth species index.
		 * \return Fourth derivative of the mean field.
		 */
		virtual double d4v(int i, int j, int k, int l) const;
		
		/**
		 * \brief Calculates the temperature derivative of the mean field.
		 * 
		 * \return Temperature derivative of the mean field.
		 */
		virtual double dvdT() const;
		
		/**
		 * \brief Sets the densities of particle species.
		 * 
		 * \param n Vector of densities.
		 */
		virtual void SetDensities(const std::vector<double>& n);
		protected:
		/**
		 * \brief Computes the components based on the mean field parameters.
		 */
		virtual void ComputeComponents();
		MeanFieldModelBase* m_mfmodel;
		std::vector<double> m_charges;
		double m_nB;
	};

} // namespace thermalfist

#endif // REALGASMODELS_H
