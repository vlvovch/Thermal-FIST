#include "HRGRealGas/MeanFieldModelsMulti.h"
#include <cmath>

namespace thermalfist {

	double MeanFieldModelMultiVDW::v() const {
		double ret = 0.;
		for (int i = 0; i < m_N; ++i) {
			for (int j = 0; j < m_N; ++j) {
				ret += -m_a[i][j] * m_densities[i] * m_densities[j];
			}
		}
		return ret;
	}

	double MeanFieldModelMultiVDW::dv(int i) const
	{
		double ret = 0.;
		for (int j = 0; j < m_N; ++j) {
			ret += -(m_a[i][j] + m_a[j][i]) * m_densities[j];
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
		for (int i = 0; i < m_N; ++i) {
			for (int j = 0; j < m_N; ++j) {
				ret += -m_dadT[i][j] * m_densities[i] * m_densities[j];
			}
		}
		return ret;
	}

	MeanFieldModelComponents::~MeanFieldModelComponents()
	{
		for (int i = 0; i < m_N_components; ++i) {
			if (m_mfmodels[i] != NULL) {
				delete m_mfmodels[i];
				m_mfmodels[i] = NULL;
			}
		}
	}

	double MeanFieldModelComponents::v() const
	{
		double ret = 0.;
		for (int i = 0; i < m_N_components; ++i) {
			ret += m_mfmodels[i]->v(m_densities_components[i]);
		}
		return ret;
	}

	double MeanFieldModelComponents::dv(int i) const
	{
		int tind = m_indices[i];
		return m_mfmodels[tind]->dv(1, m_densities_components[tind]);
	}

	double MeanFieldModelComponents::d2v(int i, int j) const
	{
		bool fl = 1;
		fl &= (m_indices[i] == m_indices[j]);
		if (!fl) return 0.;

		int tind = m_indices[i];
		return m_mfmodels[tind]->dv(2, m_densities_components[tind]);
	}

	double MeanFieldModelComponents::d3v(int i, int j, int k) const
	{
		bool fl = 1;
		fl &= (m_indices[i] == m_indices[j]);
		fl &= (m_indices[i] == m_indices[k]);
		if (!fl) return 0.;

		int tind = m_indices[i];
		return m_mfmodels[tind]->dv(3, m_densities_components[tind]);
	}

	double MeanFieldModelComponents::d4v(int i, int j, int k, int l) const
	{
		bool fl = 1;
		fl &= (m_indices[i] == m_indices[j]);
		fl &= (m_indices[i] == m_indices[k]);
		fl &= (m_indices[i] == m_indices[l]);
		if (!fl) return 0.;

		int tind = m_indices[i];
		return m_mfmodels[tind]->dv(4, m_densities_components[tind]);
	}

	double MeanFieldModelComponents::dvdT() const
	{
		double ret = 0.;
		for (int i = 0; i < m_N_components; ++i) {
			ret += m_mfmodels[i]->dvdT(m_densities_components[i]);
		}
		return ret;
	}

	void MeanFieldModelComponents::SetDensities(const std::vector<double>& n)
	{
		MeanFieldModelMultiBase::SetDensities(n);

		m_densities_components = std::vector<double>(m_N_components, 0.);

		for (int i = 0; i < n.size(); ++i)
			m_densities_components[m_indices[i]] += n[i];
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

} // namespace thermalfist