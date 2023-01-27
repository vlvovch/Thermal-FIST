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

	// Default: point-particle
	class ExcludedVolumeModelMultiBase {
	public:
		ExcludedVolumeModelMultiBase(int N) : 
			m_N(N), 
			m_densities(std::vector<double>(m_N,0.)), 
			m_components(std::vector<int>(m_N, 0)) 
		{
			ComputeComponents();
		}
		virtual ~ExcludedVolumeModelMultiBase() { }
		virtual double f(int i) const { return 1.; }
		virtual double df(int i, int j) const { return 0.; }
		virtual double d2f(int i, int j, int k) const { return 0.; }
		virtual double d3f(int i, int j, int k, int l) const { return 0.; }
		virtual double d4f(int i, int j, int k, int l, int m) const { return 0.; }
		virtual double dfdT(int i) const { return 0.; }
		virtual std::vector<double> nsol(const std::vector<double>& ntil) { return ntil; }
		virtual std::vector<double> nsolBroyden(const std::vector<double>& ntil);
		virtual std::vector<double> nsolBroydenComponents(const std::vector<double>& ntil);
		virtual void SetDensities(const std::vector<double>& n) { m_densities = n; }
		virtual const std::vector<int>& ComponentIndices() const { return m_components; }
		virtual const std::vector<int>& ComponentIndicesFrom() const { return m_componentsFrom; }
		virtual const int ComponentsNumber() const { return m_componentsNumber; }
	protected:
		virtual void ComputeComponents();
		int m_N;
		std::vector<double> m_densities;
		std::vector<int> m_components;
		std::vector<int> m_componentsFrom;
		int m_componentsNumber;

	  class BroydenEquationsEVMulti : public BroydenEquations
		{
		public:
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

			std::vector<double> Equations(const std::vector<double>& x);
		private:
			ExcludedVolumeModelMultiBase* m_EVM;
			const std::vector<double>* m_ntil;
			bool m_componentsMode;
		};

		class BroydenJacobianEVMulti : public BroydenJacobian
		{
		public:
			BroydenJacobianEVMulti(ExcludedVolumeModelMultiBase* evmodel, 
				const std::vector<double>* ntil,
				bool componentsMode = false) :
				BroydenJacobian(), 
				m_EVM(evmodel),
				m_ntil(ntil),
				m_componentsMode(componentsMode)
			{
			}

			std::vector<double> Jacobian(const std::vector<double>& x);
		private:
			ExcludedVolumeModelMultiBase* m_EVM;
			const std::vector<double>* m_ntil;
			bool m_componentsMode;
		};


	};

	// Diagonal VDW
	class ExcludedVolumeModelDiagonalVDW :
		public ExcludedVolumeModelMultiBase {
	public:
		ExcludedVolumeModelDiagonalVDW(const std::vector<double>& b, const std::vector<double>& dbdT = std::vector<double>())
			: ExcludedVolumeModelMultiBase(b.size()), m_b(b), m_dbdT(dbdT) {
			if (m_dbdT.size() == 0)
				m_dbdT = std::vector<double>(m_b.size(), 0.);
			ComputeComponents();
		}
		virtual double f(int i) const;
		virtual double df(int i, int j)											const;
		virtual double d2f(int i, int j, int k)							const { return 0.; }
		virtual double d3f(int i, int j, int k, int l)				const { return 0.; }
		virtual double d4f(int i, int j, int k, int l, int m) const { return 0.; }
		virtual double dfdT(int i)											      const;
		virtual std::vector<double> nsol(const std::vector<double>& nid);
	protected:
		virtual void ComputeComponents();
		std::vector<double> m_b;
		std::vector<double> m_dbdT;
	};

	// Diagonal generalized
	// The objected pointed by ExcludedVolumeModelBase will be deleted on destruction
	class ExcludedVolumeModelDiagonalGeneralized :
		public ExcludedVolumeModelMultiBase {
	public:
		ExcludedVolumeModelDiagonalGeneralized(
			ExcludedVolumeModelBase *evmodelsingle,
			const std::vector<double>& b, 
			const std::vector<double>& dbdT = std::vector<double>())
			: ExcludedVolumeModelMultiBase(b.size()), m_evmodelsingle(evmodelsingle), m_b(b), m_dbdT(dbdT), m_eta(0.) {
			if (m_dbdT.size() == 0)
				m_dbdT = std::vector<double>(m_b.size(), 0.);
			ComputeComponents();
		}
		virtual ~ExcludedVolumeModelDiagonalGeneralized();
		virtual double f(int i) const;
		virtual double df(int i, int j)											const;
		virtual double d2f(int i, int j, int k)							const;
		virtual double d3f(int i, int j, int k, int l)				const;
		virtual double d4f(int i, int j, int k, int l, int m) const;
		virtual double dfdT(int i)											      const;
		virtual std::vector<double> nsol(const std::vector<double>& nid);
		virtual void SetDensities(const std::vector<double>& n);
	protected:
		virtual void ComputeComponents();
		double GetEta(const std::vector<double>& n) const;
		ExcludedVolumeModelBase* m_evmodelsingle;
		std::vector<double> m_b;
		std::vector<double> m_dbdT;
		double m_eta;
	};


	// Crossterms VDW
	class ExcludedVolumeModelCrosstermsVDW :
		public ExcludedVolumeModelMultiBase {
	public:
		ExcludedVolumeModelCrosstermsVDW(const std::vector< std::vector<double> >& b, const std::vector< std::vector<double> >& dbdT = std::vector< std::vector<double> >())
			: ExcludedVolumeModelMultiBase(b.size()), m_b(b), m_dbdT(dbdT) {
			if (m_dbdT.size() == 0)
				m_dbdT = std::vector< std::vector<double> >(m_b.size(), std::vector<double>(m_b.size(), 0.));
			ComputeComponents();
		}
		virtual double f(int i) const;
		virtual double df(int i, int j)											  const;
		virtual double d2f(int i, int j, int k)								const { return 0.; }
		virtual double d3f(int i, int j, int k, int l)				const { return 0.; }
		virtual double d4f(int i, int j, int k, int l, int m) const { return 0.; }
		virtual double dfdT(int i)											      const;
		virtual std::vector<double> nsol(const std::vector<double>& nid);
	protected:
		virtual void ComputeComponents();
		std::vector< std::vector<double> > m_b;
		std::vector< std::vector<double> > m_dbdT;
	};

	// Crossterms generalized
	// The objected pointed by ExcludedVolumeModelBase will be deleted on destruction
	class ExcludedVolumeModelCrosstermsGeneralized :
		public ExcludedVolumeModelMultiBase {
	public:
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
		virtual ~ExcludedVolumeModelCrosstermsGeneralized();
		virtual double f(int i) const;
		virtual double df(int i, int j)											  const;
		virtual double d2f(int i, int j, int k)								const;
		virtual double d3f(int i, int j, int k, int l)				const;
		virtual double d4f(int i, int j, int k, int l, int m) const;
		virtual double dfdT(int i)											      const;
		virtual std::vector<double> nsol(const std::vector<double>& nid);
		virtual void SetDensities(const std::vector<double>& n);
	protected:
		virtual void ComputeComponents();
		double GetEta(int i, const std::vector<double>& n) const;
		ExcludedVolumeModelBase* m_evmodelsingle;
		std::vector< std::vector<double> > m_b;
		std::vector< std::vector<double> > m_dbdT;
		std::vector<double> m_etas;
		bool m_componentsDisconnected;
	};

	// Components
	class ExcludedVolumeModelComponents :
		public ExcludedVolumeModelMultiBase {
	public:
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
		virtual ~ExcludedVolumeModelComponents();

		virtual double f(int i) const;
		virtual double df(int i, int j) const;
		virtual double d2f(int i, int j, int k) const;
		virtual double d3f(int i, int j, int k, int l) const;
		virtual double d4f(int i, int j, int k, int l, int m) const;
		virtual double dfdT(int i)											      const;
		virtual std::vector<double> nsol(const std::vector<double>& nid);
		virtual void SetDensities(const std::vector<double>& n);
	protected:
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
