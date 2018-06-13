#pragma once
#ifndef THERMALMODELVDWFULL_H
#define THERMALMODELVDWFULL_H

#include "HRGBase/ThermalModelBase.h"


class ThermalModelVDWFull : public ThermalModelBase
{
public:
	ThermalModelVDWFull(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters(), int mode = 0);

	virtual ~ThermalModelVDWFull(void);

	void FillChemicalPotentials();
	virtual void SetChemicalPotentials(const std::vector<double> & chem = std::vector<double>(0));

	void FillVirial(const std::vector<double> & ri = std::vector<double>(0));
	void FillVirialEV(const std::vector< std::vector<double> > & bij = std::vector< std::vector<double> > (0));
	void FillAttraction(const std::vector< std::vector<double> > & aij = std::vector< std::vector<double> >(0));

	virtual void ReadInteractionParameters(const std::string &filename);
	virtual void WriteInteractionParameters(const std::string &filename);

	void SetVirial(int i, int j, double b) { if (i >= 0 && i < m_Virial.size() && j >= 0 && j < m_Virial[i].size()) m_Virial[i][j] = b; }
	void SetAttraction(int i, int j, double a) { if (i >= 0 && i < m_Attr.size() && j >= 0 && j < m_Attr[i].size())     m_Attr[i][j] = a; }

	double VirialCoefficient(int i, int j) const;
	double AttractionCoefficient(int i, int j) const;

	void SetVirialdT(int i, int j, double dbdT) { if (i >= 0 && i < m_VirialdT.size() && j >= 0 && j < m_VirialdT[i].size()) m_VirialdT[i][j] = dbdT; }
	void SetAttractiondT(int i, int j, double dadT) { if (i >= 0 && i < m_AttrdT.size() && j >= 0 && j < m_AttrdT[i].size())     m_AttrdT[i][j] = dadT; }

	double VirialCoefficientdT(int i, int j) const;
	double AttractionCoefficientdT(int i, int j) const;

	void SetTemperatureDependentAB(bool Tdep) { m_TemperatureDependentAB = Tdep; }
	bool TemperatureDependentAB() const { return m_TemperatureDependentAB; }

	virtual void SetParameters(double T, double muB, double muS, double muQ, double gammaS, double V, double R);

	virtual void SetParameters(const ThermalModelParameters& params);

	void UpdateParameters();

	virtual void ChangeTPS(ThermalParticleSystem *TPS_);

	std::vector<double> SearchSolutionsSingle(const std::vector<double> & muStarInit);
	std::vector<double> SearchSolutionsSingleBroyden(const std::vector<double> & muStarInit);
	std::vector<double> SearchSolutionsSingleBroydenOptimized(const std::vector<double> & muStarInit);
	std::vector<double> SearchSolutionsSingleBroydenOptimized2(const std::vector<double> & muStarInit);

	std::vector<double> SearchSolutionsMulti(int iters = 300);
	void SearchSolutions();

	virtual void CalculateDensities();

	virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);
	virtual std::vector< std::vector<double> >  CalculateFluctuations(int order);


	void CalculateTwoParticleCorrelations();
	// TODO higher orders
	void CalculateFluctuations();


	virtual double CalculateEnergyDensity();

	virtual double CalculateEntropyDensity();

	// Dummy
	virtual double CalculateBaryonMatterEntropyDensity();

	virtual double CalculateMesonMatterEntropyDensity();

	virtual double CalculatePressure();

	virtual double ParticleScalarDensity(int part);

	void SetMultipleSolutionsMode(bool search) { m_SearchMultipleSolutions = search; }
	bool UseMultipleSolutionsMode() const { return m_SearchMultipleSolutions; }

	double GetMaxDiff() const { return m_MaxDiff; }
	bool   IsLastSolutionOK() const { return m_LastBroydenSuccessFlag; }

	double DensityId(int index) { return m_DensitiesId[index]; }

	double MuStar(int index) const { return m_MuStar[index]; }
	std::vector<double> GetMuStar() const { return m_MuStar; }
	void SetMuStar(const std::vector<double> & MuStar) { m_MuStar = MuStar; }

protected:
	std::vector<double> m_DensitiesId;
	std::vector<double> m_scaldens;

	std::vector< std::vector<double> > m_Virial;
	std::vector< std::vector<double> > m_Attr;
	std::vector< std::vector<double> > m_VirialdT;
	std::vector< std::vector<double> > m_AttrdT;

	std::vector< std::vector<double> > m_chi;
	std::vector<double> m_chiarb;

	bool   m_TemperatureDependentAB;


	virtual void CalculateDensitiesOld();
	virtual void CalculateDensitiesNew();
private:
	bool   m_SearchMultipleSolutions;
	bool   m_LastBroydenSuccessFlag;
	double m_MaxDiff;

	std::vector<double> m_MuStar;

	std::vector<int> m_MapTodMuStar;
	std::vector<int> m_MapFromdMuStar;
	std::vector< std::vector<int> > m_dMuStarIndices;
};

#endif

