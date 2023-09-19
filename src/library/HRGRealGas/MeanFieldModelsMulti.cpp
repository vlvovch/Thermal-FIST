#include "HRGRealGas/MeanFieldModelsMulti.h"
#include <cmath>
#include <map>

namespace thermalfist {

	double MeanFieldModelMultiVDW::v() const {
		double ret = 0.;
		//for (int i = 0; i < m_N; ++i) {
		//	for (int j = 0; j < m_N; ++j) {
		//		ret += -m_a[i][j] * m_densities[i] * m_densities[j];
		//	}
		//}
		for (int i = 0; i < m_componentsNumber; ++i) {
			for (int j = 0; j < m_componentsNumber; ++j) {
				int ti = m_componentsFrom[i];
				int tj = m_componentsFrom[j];
				ret += -m_a[ti][tj] * m_densities_components[i] * m_densities_components[j];
			}
		}
		return ret;
	}

	double MeanFieldModelMultiVDW::dv(int i) const
	{
		double ret = 0.;
		//for (int j = 0; j < m_N; ++j) {
		//	ret += -(m_a[i][j] + m_a[j][i]) * m_densities[j];
		//}
		for (int j = 0; j < m_componentsNumber; ++j) {
			int tj = m_componentsFrom[j];
			ret += -(m_a[i][tj] + m_a[tj][i]) * m_densities_components[j];
		}
		return ret;
	}

	double MeanFieldModelMultiVDW::d2v(int i, int j) const
	{
    return -(m_a[i][j] + m_a[j][i]);
	}

	double MeanFieldModelMultiVDW::dvdT() const
	{
		double ret = 0.;
		//for (int i = 0; i < m_N; ++i) {
		//	for (int j = 0; j < m_N; ++j) {
		//		ret += -m_dadT[i][j] * m_densities[i] * m_densities[j];
		//	}
		//}
		for (int i = 0; i < m_componentsNumber; ++i) {
			for (int j = 0; j < m_componentsNumber; ++j) {
				int ti = m_componentsFrom[i];
				int tj = m_componentsFrom[j];
				ret += -m_dadT[ti][tj] * m_densities_components[i] * m_densities_components[j];
			}
		}
		return ret;
	}

	void MeanFieldModelMultiVDW::ComputeComponents()
	{
		m_components = std::vector<int>(m_N, 0);
		m_componentsFrom = std::vector<int>();
		std::map<std::vector<double>, int> MapVDW;
		m_componentsNumber = 0;
		for (int i = 0; i < m_N; ++i) {
			std::vector<double> param(0);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_a[i][j] + m_a[j][i]);
			for (int j = 0; j < m_N; ++j)
				param.push_back(m_dadT[i][j] + m_dadT[j][i]);

			if (MapVDW.count(param) == 0) {
				MapVDW[param] = m_componentsNumber;
				m_componentsFrom.push_back(i);
				m_componentsNumber++;
			}
			m_components[i] = MapVDW[param];
		}
		printf("VDW Components: %d\n", m_componentsNumber);
	}

	MeanFieldModelComponents::~MeanFieldModelComponents()
	{
		for (int i = 0; i < m_componentsNumber; ++i) {
			if (m_mfmodels[i] != NULL) {
				delete m_mfmodels[i];
				m_mfmodels[i] = NULL;
			}
		}
	}

	double MeanFieldModelComponents::v() const
	{
		double ret = 0.;
		for (int i = 0; i < m_componentsNumber; ++i) {
			ret += m_mfmodels[i]->v(m_densities_components[i]);
		}
		return ret;
	}

	double MeanFieldModelComponents::dv(int i) const
	{
		int tind = m_components[i];
		return m_mfmodels[tind]->dv(1, m_densities_components[tind]);
	}

	double MeanFieldModelComponents::d2v(int i, int j) const
	{
		bool fl = 1;
		fl &= (m_components[i] == m_components[j]);
		if (!fl) return 0.;

		int tind = m_components[i];
		return m_mfmodels[tind]->dv(2, m_densities_components[tind]);
	}

	double MeanFieldModelComponents::d3v(int i, int j, int k) const
	{
		bool fl = 1;
		fl &= (m_components[i] == m_components[j]);
		fl &= (m_components[i] == m_components[k]);
		if (!fl) return 0.;

		int tind = m_components[i];
		return m_mfmodels[tind]->dv(3, m_densities_components[tind]);
	}

	double MeanFieldModelComponents::d4v(int i, int j, int k, int l) const
	{
		bool fl = 1;
		fl &= (m_components[i] == m_components[j]);
		fl &= (m_components[i] == m_components[k]);
		fl &= (m_components[i] == m_components[l]);
		if (!fl) return 0.;

		int tind = m_components[i];
		return m_mfmodels[tind]->dv(4, m_densities_components[tind]);
	}

	double MeanFieldModelComponents::dvdT() const
	{
		double ret = 0.;
		for (int i = 0; i < m_componentsNumber; ++i) {
			ret += m_mfmodels[i]->dvdT(m_densities_components[i]);
		}
		return ret;
	}

	//void MeanFieldModelComponents::SetDensities(const std::vector<double>& n)
	//{
	//	MeanFieldModelMultiBase::SetDensities(n);

	//	m_densities_components = std::vector<double>(m_componentsNumber, 0.);
	//	for (int i = 0; i < n.size(); ++i)
	//		m_densities_components[m_components[i]] += n[i];
	//}

	void MeanFieldModelComponents::ComputeComponents()
	{
		m_componentsFrom.resize(m_componentsNumber);
		for (int i = m_N - 1; i >= 0; --i) {
			m_componentsFrom[m_components[i]] = i;
		}
	}

	MeanFieldModelChargeDensityDependent::~MeanFieldModelChargeDensityDependent()
	{
		if (m_mfmodel != NULL) {
			delete m_mfmodel;
			m_mfmodel = NULL;
		}
	}

	double MeanFieldModelChargeDensityDependent::v() const
	{
		return m_mfmodel->v(m_nB);
	}

	double MeanFieldModelChargeDensityDependent::dv(int i) const
	{
		return m_mfmodel->dv(1, m_nB) * m_charges[i];
	}

	double MeanFieldModelChargeDensityDependent::d2v(int i, int j) const
	{
		return m_mfmodel->dv(2, m_nB) * m_charges[i];
	}

	double MeanFieldModelChargeDensityDependent::d3v(int i, int j, int k) const
	{
		return m_mfmodel->dv(3, m_nB) * m_charges[i];
	}

	double MeanFieldModelChargeDensityDependent::d4v(int i, int j, int k, int l) const
	{
		return m_mfmodel->dv(4, m_nB) * m_charges[i];
	}

	double MeanFieldModelChargeDensityDependent::dvdT() const
	{
		return m_mfmodel->dvdT(m_nB);
	}

	void MeanFieldModelChargeDensityDependent::SetDensities(const std::vector<double>& n)
	{
		MeanFieldModelMultiBase::SetDensities(n);
		m_nB = 0.0;
		for (int i = 0; i < m_N; ++i)
			m_nB += m_charges[i] * m_densities[i];
	}

	void MeanFieldModelChargeDensityDependent::ComputeComponents()
	{
		m_components = std::vector<int>(m_N, 0);
		m_componentsFrom = std::vector<int>();
		std::map<double, int> MapMF;
		m_componentsNumber = 0;
		for (int i = 0; i < m_N; ++i) {
			double param = m_charges[i];

			if (MapMF.count(param) == 0) {
				MapMF[param] = m_componentsNumber;
				m_componentsFrom.push_back(i);
				m_componentsNumber++;
			}
			m_components[i] = MapMF[param];
		}
		printf("Mean-field Components: %d\n", m_componentsNumber);
	}

	void MeanFieldModelMultiBase::SetDensities(const std::vector<double>& n)
	{
		m_densities = n;

		m_densities_components = std::vector<double>(m_componentsNumber, 0.);
		for (int i = 0; i < n.size(); ++i)
			m_densities_components[m_components[i]] += n[i];
	}

	void MeanFieldModelMultiBase::ComputeComponents()
	{
		m_components = std::vector<int>(m_N, 0);
		for (int i = 0; i < m_N; ++i)
			m_components[i] = i;
		m_componentsNumber = m_N;
		m_componentsFrom = m_components;
	}

} // namespace thermalfist