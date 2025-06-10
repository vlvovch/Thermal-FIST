#include "HRGRealGas/ExcludedVolumeModelsMulti.h"
#include <cmath>
#include <map>

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

	void ExcludedVolumeModelDiagonalVDW::ComputeComponents()
	{
		m_components = std::vector<int>(m_N, 0);
		m_componentsFrom = std::vector<int>();
		std::map<std::pair<double,double>, int> MapEV;
		m_componentsNumber = 0;
		for (int i = 0; i < m_N; ++i) {
			std::pair<double, double> param = std::make_pair(m_b[i], m_dbdT[i]);
			if (MapEV.count(param) == 0) {
				MapEV[param] = m_componentsNumber;
				m_componentsFrom.push_back(i);
				m_componentsNumber++;
			}
			m_components[i] = MapEV[param];
		}
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

	void ExcludedVolumeModelCrosstermsVDW::ComputeComponents()
	{
		m_components = std::vector<int>(m_N, 0);
		m_componentsFrom = std::vector<int>();
		std::map<std::vector<double>, int> MapEV;
		m_componentsNumber = 0;
		for (int i = 0; i < m_N; ++i) {
			std::vector<double> param(0);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_b[i][j]);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_b[j][i]);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_dbdT[i][j]);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_dbdT[j][i]);

			if (MapEV.count(param) == 0) {
				MapEV[param] = m_componentsNumber;
				m_componentsFrom.push_back(i);
				m_componentsNumber++;
			}
			m_components[i] = MapEV[param];
		}
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

	void ExcludedVolumeModelDiagonalGeneralized::ComputeComponents()
	{
		m_components = std::vector<int>(m_N, 0);
		m_componentsFrom = std::vector<int>();
		std::map<std::pair<double, double>, int> MapEV;
		m_componentsNumber = 0;
		for (int i = 0; i < m_N; ++i) {
			std::pair<double, double> param = std::make_pair(m_b[i], m_dbdT[i]);
			if (MapEV.count(param) == 0) {
				MapEV[param] = m_componentsNumber;
				m_componentsFrom.push_back(i);
				m_componentsNumber++;
			}
			m_components[i] = MapEV[param];
		}
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

		std::vector<double> dens(m_EVM->m_N, 0.);
		for (int i = 0; i < m_EVM->m_N; ++i) {
			int ti = i;
			if (m_componentsMode)
				ti = m_EVM->ComponentIndices()[i];
			dens[i] = x[ti] * m_ntil->operator[](i);
		}
		m_EVM->SetDensities(dens);

		for (int i = 0; i < m_N; ++i) {
			int fi = i;
			if (m_componentsMode)
				fi = m_EVM->ComponentIndicesFrom()[i];
			ret[i] = x[i] - m_EVM->f(fi);
		}

		return ret;
	}

	std::vector<double> ExcludedVolumeModelMultiBase::BroydenJacobianEVMulti::Jacobian(const std::vector<double>& x)
	{
		int N = x.size();
		std::vector<double> ret(N * N);

		std::vector<double> dens(m_EVM->m_N, 0.);
		for (int i = 0; i < m_EVM->m_N; ++i) {
			int ti = i;
			if (m_componentsMode)
				ti = m_EVM->ComponentIndices()[i];
			dens[i] = x[ti] * m_ntil->operator[](i);
		}
		m_EVM->SetDensities(dens);

		std::vector<double> ntilcomp(N, 0.);
		for (int i = 0; i < m_EVM->m_N; ++i) {
			int ti = i;
			if (m_componentsMode)
				ti = m_EVM->ComponentIndices()[i];
			ntilcomp[ti] = m_ntil->operator[](i);
		}

		for (int i = 0; i < N; ++i) {
			int fi = i;
			if (m_componentsMode)
				fi = m_EVM->ComponentIndicesFrom()[i];
			for (int j = 0; j < N; ++j) {
				int fj = j;
				if (m_componentsMode)
					fj = m_EVM->ComponentIndicesFrom()[j];
				ret[i * N + j] = 0.;
				if (i == j)
					ret[i * N + j] += 1.;
				ret[i * N + j] += -m_EVM->df(fi, fj) * ntilcomp[j];
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

	std::vector<double> ExcludedVolumeModelMultiBase::nsolBroydenComponents(const std::vector<double>& ntil)
	{
		BroydenEquationsEVMulti eqs(this, &ntil, true);
		BroydenJacobianEVMulti jac(this, &ntil, true);
		Broyden broydn(&eqs, &jac);
		Broyden::BroydenSolutionCriterium crit(1.0E-10);

		std::vector<double> ret(m_N, 0.);
		std::vector<double> sol(m_componentsNumber, 0.);

	  sol = broydn.Solve(sol, &crit);
		for (int i = 0; i < m_N; ++i)
			ret[i] = ntil[i] * sol[m_components[i]];

		return ret;
	}

	void ExcludedVolumeModelMultiBase::ComputeComponents()
	{
		m_components = std::vector<int>(m_N, 0);
		m_componentsNumber = 1;
		m_componentsFrom = std::vector<int>(1, 0);
	}


	void ExcludedVolumeModelCrosstermsGeneralized::SetDensities(const std::vector<double>& n)
	{
		ExcludedVolumeModelMultiBase::SetDensities(n);
		for (int i = 0; i < m_componentsNumber; ++i) {
			m_etas[i] = GetEta(m_componentsFrom[i], m_densities);
		}
		//for (int i = 0; i < m_N; ++i) {
		//	m_etas[i] = GetEta(i, m_densities);
		//}
	}

	void ExcludedVolumeModelCrosstermsGeneralized::ComputeComponents()
	{
		m_components = std::vector<int>(m_N, 0);
		m_componentsFrom = std::vector<int>();
		std::map<std::vector<double>, int> MapEV;
		m_componentsNumber = 0;
		for (int i = 0; i < m_N; ++i) {
			std::vector<double> param(0);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_b[i][j]);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_b[j][i]);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_dbdT[i][j]);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_dbdT[j][i]);

			if (MapEV.count(param) == 0) {
				MapEV[param] = m_componentsNumber;
				m_componentsFrom.push_back(i);
				m_componentsNumber++;
			}
			m_components[i] = MapEV[param];
		}

		m_componentsDisconnected = true;
		for (int i = 0; i < m_componentsNumber; ++i) {
			for (int j = 0; j < m_componentsNumber; ++j) {
				if (i != j) {
					m_componentsDisconnected &= (   m_b[m_componentsFrom[i]][m_componentsFrom[j]] == 0.0);
					m_componentsDisconnected &= (m_dbdT[m_componentsFrom[i]][m_componentsFrom[j]] == 0.0);
				}
			}
		}

		printf("EV Components: %d  Are Disconnected: %d\n", m_componentsNumber, (int)m_componentsDisconnected);
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
		//return m_evmodelsingle->f(m_etas[i]);
		return m_evmodelsingle->f(m_etas[ComponentIndices()[i]]);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::df(int i, int j) const
	{
		//return m_evmodelsingle->df(1, m_etas[i]) * (m_b[j][i] / 4.);
		return m_evmodelsingle->df(1, m_etas[ComponentIndices()[i]]) * (m_b[j][i] / 4.);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::d2f(int i, int j, int k) const
	{
		//return m_evmodelsingle->df(2, m_etas[i]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.);
		return m_evmodelsingle->df(2, m_etas[ComponentIndices()[i]]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::d3f(int i, int j, int k, int l) const
	{
		//return m_evmodelsingle->df(3, m_etas[i]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.) * (m_b[l][i] / 4.);
		return m_evmodelsingle->df(3, m_etas[ComponentIndices()[i]]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.) * (m_b[l][i] / 4.);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::d4f(int i, int j, int k, int l, int m) const
	{
		//return m_evmodelsingle->df(4, m_etas[i]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.) * (m_b[l][i] / 4.) * (m_b[m][i] / 4.);
		return m_evmodelsingle->df(4, m_etas[ComponentIndices()[i]]) * (m_b[j][i] / 4.) * (m_b[k][i] / 4.) * (m_b[l][i] / 4.) * (m_b[m][i] / 4.);
	}

	double ExcludedVolumeModelCrosstermsGeneralized::dfdT(int i) const
	{
		double ret = 0.;
		for (int j = 0; j < m_N; ++j)
			ret += m_dbdT[j][i] * m_densities[j] / 4.;
		//ret *= m_evmodelsingle->df(1, m_etas[i]);
		ret *= m_evmodelsingle->df(1, m_etas[ComponentIndices()[i]]);
		return ret;
	}

	std::vector<double> ExcludedVolumeModelCrosstermsGeneralized::nsol(const std::vector<double>& nid)
	{
		if (!m_componentsDisconnected) {
			return nsolBroydenComponents(nid);
			return nsolBroyden(nid);
		}

		// Disconnected components, use binary search separately for each component
		std::vector<double> nids_comp(m_componentsNumber, 0.);
		for (int i = 0; i < nid.size(); ++i)
			nids_comp[m_components[i]] += nid[i];

		std::vector<double> fvals(m_componentsNumber, 0.);
		for (int ic = 0; ic < m_componentsNumber; ++ic) {
			int ti = m_componentsFrom[ic];
			double etaid = m_b[ti][ti] * nids_comp[ic] / 4.;
			double eta = m_evmodelsingle->etasol(etaid);
			fvals[ic] = m_evmodelsingle->f(eta);
		}
		std::vector<double> ret(m_N, 0.);
		for (int i = 0; i < nid.size(); ++i) {
			ret[i] = fvals[m_components[i]] * nid[i];
		}
		return ret;



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
			m_densities_components[m_components[i]] += n[i];
	}

	void ExcludedVolumeModelComponents::ComputeComponents()
	{
		m_componentsFrom.resize(m_componentsNumber);
		for (int i = m_N - 1; i >= 0; --i) {
			m_componentsFrom[m_components[i]] = i;
		}
	}

	ExcludedVolumeModelComponents::~ExcludedVolumeModelComponents()
	{
		for (int i = 0; i < m_componentsNumber; ++i) {
			if (m_evmodels[i] != NULL) {
				delete m_evmodels[i];
				m_evmodels[i] = NULL;
			}
		}
	}

	double ExcludedVolumeModelComponents::f(int i) const
	{
		int tind = m_components[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->f(eta);
	}

	double ExcludedVolumeModelComponents::df(int i, int j) const
	{
		bool fl = 1;
		fl &= (m_components[i] == m_components[j]);
		if (!fl) return 0.;

		int tind = m_components[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(1, eta) * (m_b[tind] / 4.);
	}

	double ExcludedVolumeModelComponents::d2f(int i, int j, int k) const
	{
		bool fl = 1;
		fl &= (m_components[i] == m_components[j]);
		fl &= (m_components[i] == m_components[k]);
		if (!fl) return 0.;

		int tind = m_components[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(2, eta) * (m_b[tind] / 4.) * (m_b[tind] / 4.);
	}

	double ExcludedVolumeModelComponents::d3f(int i, int j, int k, int l) const
	{
		bool fl = 1;
		fl &= (m_components[i] == m_components[j]);
		fl &= (m_components[i] == m_components[k]);
		fl &= (m_components[i] == m_components[l]);
		if (!fl) return 0.;

		int tind = m_components[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(3, eta) * (m_b[tind] / 4.) * (m_b[tind] / 4.) * (m_b[tind] / 4.);
	}

	double ExcludedVolumeModelComponents::d4f(int i, int j, int k, int l, int m) const
	{
		bool fl = 1;
		fl &= (m_components[i] == m_components[j]);
		fl &= (m_components[i] == m_components[k]);
		fl &= (m_components[i] == m_components[l]);
		fl &= (m_components[i] == m_components[m]);
		if (!fl) return 0.;

		int tind = m_components[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(4, eta) * (m_b[tind] / 4.) * (m_b[tind] / 4.) * (m_b[tind] / 4.) * (m_b[tind] / 4.);
	}

	double ExcludedVolumeModelComponents::dfdT(int i) const
	{
		int tind = m_components[i];
		double eta = m_b[tind] * m_densities_components[tind] / 4.;
		return m_evmodels[tind]->df(1, eta) * (m_dbdT[tind] * m_densities_components[tind] / 4.);
	}

	std::vector<double> ExcludedVolumeModelComponents::nsol(const std::vector<double>& nid)
	{
		std::vector<double> nids_comp(m_componentsNumber, 0.);
		for (int i = 0; i < nid.size(); ++i)
			nids_comp[m_components[i]] += nid[i];
		
		std::vector<double> fvals(m_componentsNumber, 0.);
		for (int ic = 0; ic < m_componentsNumber; ++ic) {
			double etaid = m_b[ic] * nids_comp[ic] / 4.;
			double eta = m_evmodels[ic]->etasol(etaid);
			fvals[ic] = m_evmodels[ic]->f(eta);
		}
		std::vector<double> ret(m_N, 0.);
		for (int i = 0; i < nid.size(); ++i) {
			ret[i] = fvals[m_components[i]] * nid[i];
		}
		return ret;
	}

} // namespace thermalfist