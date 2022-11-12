#include "HRGRealGas/ExcludedVolumeModelsMulti.h"
#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;

namespace thermalfist {

	double ExcludedVolumeModelDiagonalVDW::f(int i) const {
		const std::vector<double>& n = m_densities;
		double ret = 1.;
		for (int ii = 0; ii < n.size(); ++ii)
			ret -= m_b[ii] * n[ii];
		return ret;
	}

	double ExcludedVolumeModelDiagonalVDW::df(int i, int j) const {
		const std::vector<double>& n = m_densities;
		return -m_b[j];
	}

	double ExcludedVolumeModelDiagonalVDW::dfdT(int i) const
	{
		const std::vector<double>& n = m_densities;
		double ret = 0.;
		for (int ii = 0; ii < n.size(); ++ii)
			ret -= m_dbdT[ii] * n[ii];
		return ret;
	}

	std::vector<double> ExcludedVolumeModelDiagonalVDW::nsol(const std::vector<double>& nid) {
		std::vector<double> ret = nid;
		double denom = 1.;
		for (int ii = 0; ii < nid.size(); ++ii)
			denom += m_b[ii] * nid[ii];
		for (int ii = 0; ii < nid.size(); ++ii)
			ret[ii] = nid[ii] / denom;
		return ret;
	}
	
	double ExcludedVolumeModelCrosstermsVDW::f(int i) const {
		const std::vector<double>& n = m_densities;
		double ret = 1.;
		for (int j = 0; j < n.size(); ++j)
			ret -= m_b[j][i] * n[j];
		return ret;
	}

	double ExcludedVolumeModelCrosstermsVDW::df(int i, int j) const {
		const std::vector<double>& n = m_densities;
		return -m_b[j][i];
	}

	double ExcludedVolumeModelCrosstermsVDW::dfdT(int i) const {
		const std::vector<double>& n = m_densities;
		double ret = 0.;
		for (int j = 0; j < n.size(); ++j)
			ret -= m_dbdT[j][i] * n[j];
		return ret;
	}

	std::vector<double> ExcludedVolumeModelCrosstermsVDW::nsol(const std::vector<double>& nid) {
		std::vector<double> ret = nid;
		int NN = nid.size();

		MatrixXd densMatrix(NN, NN);
		VectorXd solVector(NN), xVector(NN);

		for (int i = 0; i < NN; ++i)
			for (int j = 0; j < NN; ++j) {
				densMatrix(i, j) = m_b[j][i] * nid[i];
				if (i == j) densMatrix(i, j) += 1.;
			}

		PartialPivLU<MatrixXd> decomp(densMatrix);

		for (int i = 0; i < NN; ++i) xVector[i] = nid[i];
		solVector = decomp.solve(xVector);

		for (int i = 0; i < NN; ++i) ret[i] = solVector[i];

		return ret;
	}

	ExcludedVolumeModelDiagonalGeneralized::~ExcludedVolumeModelDiagonalGeneralized()
	{
		if (m_evmodelsingle != NULL) {
			delete m_evmodelsingle;
			m_evmodelsingle = NULL;
		}
	}

	double ExcludedVolumeModelDiagonalGeneralized::f(int i) const
	{
		return m_evmodelsingle->f(m_eta);
	}

	double ExcludedVolumeModelDiagonalGeneralized::df(int i, int j) const
	{
		return m_evmodelsingle->df(1, m_eta) * (m_b[j] / 4.);
	}

	double ExcludedVolumeModelDiagonalGeneralized::d2f(int i, int j, int k) const
	{
		return m_evmodelsingle->df(2, m_eta) * (m_b[j] / 4.) * (m_b[k] / 4.);
	}

	double ExcludedVolumeModelDiagonalGeneralized::d3f(int i, int j, int k, int l) const
	{
		return m_evmodelsingle->df(3, m_eta) * (m_b[j] / 4.) * (m_b[k] / 4.) * (m_b[l] / 4.);
	}

	double ExcludedVolumeModelDiagonalGeneralized::d4f(int i, int j, int k, int l, int m) const
	{
		return m_evmodelsingle->df(4, m_eta) * (m_b[j] / 4.) * (m_b[k] / 4.) * (m_b[l] / 4.) * (m_b[m] / 4.);
	}

	double ExcludedVolumeModelDiagonalGeneralized::dfdT(int i) const
	{
		double ret = 0.;
		for (int i = 0; i < m_N; ++i)
			ret += m_dbdT[i] * m_densities[i] / 4.;
		ret *= m_evmodelsingle->df(1, m_eta);
		return ret;
	}

	std::vector<double> ExcludedVolumeModelDiagonalGeneralized::nsol(const std::vector<double>& nid)
	{
		double etaid = GetEta(nid);
		double eta = m_evmodelsingle->etasol(etaid);
		double fval = m_evmodelsingle->f(eta);
		std::vector<double> ret(m_N, 0.);
		for (int i = 0; i < m_N; ++i)
			ret[i] = fval * nid[i];
		return ret;
	}

	void ExcludedVolumeModelDiagonalGeneralized::SetDensities(const std::vector<double>& n)
	{
		ExcludedVolumeModelMultiBase::SetDensities(n);
		m_eta = GetEta(m_densities);
	}

	double ExcludedVolumeModelDiagonalGeneralized::GetEta(const std::vector<double>& n) const
	{
		double ret = 0.0;
		for (int i = 0; i < n.size(); ++i)
			ret += m_b[i] * n[i];
		ret /= 4.;
		return ret;
	}

	std::vector<double> ExcludedVolumeModelMultiBase::BroydenEquationsEVMulti::Equations(const std::vector<double>& x)
	{
		std::vector<double> ret(m_N, 0.);

		std::vector<double> dens(m_N, 0.);
		for (int i = 0; i < m_N; ++i) {
			dens[i] = x[i] * m_ntil->operator[](i);
		}
		m_EVM->SetDensities(dens);

		for (int i = 0; i < m_N; ++i) {
			ret[i] = x[i] - m_EVM->f(i);
		}

		return ret;
	}

	std::vector<double> ExcludedVolumeModelMultiBase::BroydenJacobianEVMulti::Jacobian(const std::vector<double>& x)
	{
		int N = x.size();
		std::vector<double> ret(N * N);

		std::vector<double> dens(N, 0.);
		for (int i = 0; i < N; ++i) {
			dens[i] = x[i] * m_ntil->operator[](i);
		}
		m_EVM->SetDensities(dens);

		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				ret[i * N + j] = 0.;
				if (i == j)
					ret[i * N + j] += 1.;
				ret[i * N + j] += -m_EVM->df(i, j) * m_ntil->operator[](j);
			}
		}

		return ret;;
	}

	std::vector<double> ExcludedVolumeModelMultiBase::nsolBroyden(const std::vector<double>& ntil)
	{
		BroydenEquationsEVMulti eqs(this, &ntil);
		BroydenJacobianEVMulti jac(this, &ntil);
		Broyden broydn(&eqs, &jac);
		Broyden::BroydenSolutionCriterium crit(1.0E-10);

		std::vector<double> ret(m_N, 0.);

		ret = broydn.Solve(ret, &crit);
		for (int i = 0; i < m_N; ++i)
			ret[i] *= ntil[i];

		return ret;
	}


	void ExcludedVolumeModelCrosstermsGeneralized::SetDensities(const std::vector<double>& n)
	{
		ExcludedVolumeModelMultiBase::SetDensities(n);
		for (int i = 0; i < m_N; ++i) {
			m_etas[i] = GetEta(i, m_densities);
		}
	}

	double ExcludedVolumeModelCrosstermsGeneralized::GetEta(int i, const std::vector<double>& n) const
	{
		double ret = 0.0;
		for (int j = 0; j < n.size(); ++j)
			ret += m_b[j][i] * n[j];
		ret /= 4.;
		return ret;
	}

	ExcludedVolumeModelCrosstermsGeneralized::~ExcludedVolumeModelCrosstermsGeneralized()
	{
		if (m_evmodelsingle != NULL) {
			delete m_evmodelsingle;
			m_evmodelsingle = NULL;
		}
	}

	double ExcludedVolumeModelCrosstermsGeneralized::f(int i) const
	{
		return m_evmodelsingle->f(m_etas[i]);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::df(int i, int j) const
	{
		return m_evmodelsingle->df(1, m_etas[i]) * (m_b[j][i] / 4.);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::d2f(int i, int j, int k) const
	{
		return m_evmodelsingle->df(2, m_etas[i]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::d3f(int i, int j, int k, int l) const
	{
		return m_evmodelsingle->df(3, m_etas[i]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.) * (m_b[l][i] / 4.);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::d4f(int i, int j, int k, int l, int m) const
	{
		return m_evmodelsingle->df(4, m_etas[i]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.) * (m_b[l][i] / 4.) * (m_b[m][i] / 4.);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::dfdT(int i) const
	{
		double ret = 0.;
		for (int j = 0; j < m_N; ++j)
			ret += m_dbdT[j][i] * m_densities[j] / 4.;
		ret *= m_evmodelsingle->df(1, m_etas[i]);
		return ret;
	}

	std::vector<double> ExcludedVolumeModelCrosstermsGeneralized::nsol(const std::vector<double>& nid)
	{
		return nsolBroyden(nid);
		//std::vector<double> etasid(m_N);
		//for (int i = 0; i < m_N; ++i) 
		//	etasid[i] = GetEta(i, nid);
		//for (int i = 0; i < m_N; ++i) {
		//	double etaid = etasid[i];
		//	double eta = m_evmodelsingle->etasol(etaid);
		//	double fval = m_evmodelsingle->f(eta);
		//	ret[i] = fval * nid[i];
		//}
		//return ret;
	}

	void ExcludedVolumeModelComponents::SetDensities(const std::vector<double>& n)
	{
		ExcludedVolumeModelMultiBase::SetDensities(n);

		for (int i = 0; i < m_densities_components.size(); ++i)
			m_densities_components[i] = 0.;

		for (int i = 0; i < n.size(); ++i)
			m_densities_components[m_indices[i]] += n[i];
	}

	ExcludedVolumeModelComponents::~ExcludedVolumeModelComponents()
	{
		for (int i = 0; i < m_N_components; ++i) {
			if (m_evmodels[i] != NULL) {
				delete m_evmodels[i];
				m_evmodels[i] = NULL;
			}
		}
	}

	double ExcludedVolumeModelComponents::f(int i) const
	{
		int tind = m_indices[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->f(eta);
	}

	double ExcludedVolumeModelComponents::df(int i, int j) const
	{
		bool fl = 1;
		fl &= (m_indices[i] == m_indices[j]);
		if (!fl) return 0.;

		int tind = m_indices[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(1, eta) * (m_b[tind] / 4.);
	}

	double ExcludedVolumeModelComponents::d2f(int i, int j, int k) const
	{
		bool fl = 1;
		fl &= (m_indices[i] == m_indices[j]);
		fl &= (m_indices[i] == m_indices[k]);
		if (!fl) return 0.;

		int tind = m_indices[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(2, eta) * (m_b[tind] / 4.) * (m_b[tind] / 4.);
	}

	double ExcludedVolumeModelComponents::d3f(int i, int j, int k, int l) const
	{
		bool fl = 1;
		fl &= (m_indices[i] == m_indices[j]);
		fl &= (m_indices[i] == m_indices[k]);
		fl &= (m_indices[i] == m_indices[l]);
		if (!fl) return 0.;

		int tind = m_indices[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(3, eta) * (m_b[tind] / 4.) * (m_b[tind] / 4.) * (m_b[tind] / 4.);
	}

	double ExcludedVolumeModelComponents::d4f(int i, int j, int k, int l, int m) const
	{
		bool fl = 1;
		fl &= (m_indices[i] == m_indices[j]);
		fl &= (m_indices[i] == m_indices[k]);
		fl &= (m_indices[i] == m_indices[l]);
		fl &= (m_indices[i] == m_indices[m]);
		if (!fl) return 0.;

		int tind = m_indices[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(4, eta) * (m_b[tind] / 4.) * (m_b[tind] / 4.) * (m_b[tind] / 4.) * (m_b[tind] / 4.);
	}

	double ExcludedVolumeModelComponents::dfdT(int i) const
	{
		int tind = m_indices[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(1, eta) * (m_dbdT[tind] * m_densities_components[tind] / 4.);
	}

	std::vector<double> ExcludedVolumeModelComponents::nsol(const std::vector<double>& nid)
	{
		std::vector<double> nids_comp(m_N_components, 0.);
		for (int i = 0; i < nid.size(); ++i)
			nids_comp[m_indices[i]] += nid[i];
		
		std::vector<double> fvals(m_N_components, 0.);
		for (int ic = 0; ic < m_N_components; ++ic) {
			double etaid = m_b[ic] * nids_comp[ic] / 4.;
			double eta = m_evmodels[ic]->etasol(etaid);
			fvals[ic] = m_evmodels[ic]->f(eta);
		}
		std::vector<double> ret(m_N, 0.);
		for (int i = 0; i < nid.size(); ++i) {
			ret[i] = fvals[m_indices[i]] * nid[i];
		}
		return ret;
	}

} // namespace thermalfist