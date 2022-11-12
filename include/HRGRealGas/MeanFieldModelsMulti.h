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

	// Default: no mean field
	class MeanFieldModelMultiBase {
	public:
		MeanFieldModelMultiBase(int N) : m_N(N) { }
		virtual double v()	const { return 0.; }
		virtual double dv(int i)											  const { return 0.; }
		virtual double d2v(int i, int j)								const { return 0.; }
		virtual double d3v(int i, int j, int k)				  const { return 0.; }
		virtual double d4v(int i, int j, int k, int l)  const { return 0.; }
		virtual double dvdT()										        const { return 0.; }
		virtual void SetDensities(const std::vector<double>& n) { m_densities = n; }
	protected:
		int m_N;
		std::vector<double> m_densities;
	};

	class MeanFieldModelMultiVDW :
		public MeanFieldModelMultiBase {
	public:
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
		}
		void setAij(const std::vector< std::vector<double> >& a) { m_a = a; }
		void setdAijdT(const std::vector< std::vector<double> >& dadT) { m_dadT = dadT; }
		virtual double v()                              const;
		virtual double dv(int i)											  const;
		virtual double d2v(int i, int j)								const;
		virtual double d3v(int i, int j, int k)				  const { return 0.; }
		virtual double d4v(int i, int j, int k, int l)  const { return 0.; }
		virtual double dvdT()										        const;
	protected:
		std::vector< std::vector<double> > m_a;
		std::vector< std::vector<double> > m_dadT;
	};

	class MeanFieldModelComponents :
		public MeanFieldModelMultiBase {
	public:
		MeanFieldModelComponents(
			int components,
			const std::vector<MeanFieldModelBase*>& mfmods,
			const std::vector<int>& ind
		) 
			: MeanFieldModelMultiBase(ind.size()), 
			m_N_components(components),
			m_mfmodels(mfmods),
			m_indices(ind),
			m_densities_components(std::vector<double>(components, 0.))
		{
		}
		virtual ~MeanFieldModelComponents();

		virtual double v()                              const;
		virtual double dv(int i)											  const;
		virtual double d2v(int i, int j)								const;
		virtual double d3v(int i, int j, int k)				  const;
		virtual double d4v(int i, int j, int k, int l)  const;
		virtual double dvdT()										        const;
		virtual void SetDensities(const std::vector<double>& n);
	protected:
		int m_N_components;
		std::vector<MeanFieldModelBase*> m_mfmodels;
		std::vector<int> m_indices;
		std::vector<double> m_densities_components;
	};

	class MeanFieldModelChargeDensityDependent :
		public MeanFieldModelMultiBase {
	public:
		MeanFieldModelChargeDensityDependent(
		  MeanFieldModelBase* mfmod,
			const std::vector<double>& ind
		)
			: MeanFieldModelMultiBase(ind.size()),
			m_mfmodel(mfmod),
			m_charges(ind),
			m_nB(0.)
		{
		}
		virtual ~MeanFieldModelChargeDensityDependent();

		virtual double v()                              const;
		virtual double dv(int i)											  const;
		virtual double d2v(int i, int j)								const;
		virtual double d3v(int i, int j, int k)				  const;
		virtual double d4v(int i, int j, int k, int l)  const;
		virtual double dvdT()										        const;
		virtual void SetDensities(const std::vector<double>& n);
	protected:
		MeanFieldModelBase* m_mfmodel;
		std::vector<double> m_charges;
		double m_nB;
	};

} // namespace thermalfist

#endif // REALGASMODELS_H
