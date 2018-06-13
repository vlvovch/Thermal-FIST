#ifndef EVENTGENERATORBASE_H
#define EVENTGENERATORBASE_H


#include <sstream>

#include "HRGEventGenerator/SimpleEvent.h"
#include "HRGEventGenerator/Acceptance.h"
#include "HRGEventGenerator/RandomGenerators.h"
#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"


template <typename T>
std::string to_string_fix(T value)
{
  //create an output string stream
  std::ostringstream os ;

  //throw the value into the string stream
  os << value ;

  //convert the string stream into a string and return
  return os.str() ;
}


//class ExcludedVolumeModel;

struct EventGeneratorConfiguration {
	enum Ensemble { GCE, SCE, CE };
	enum ModelType { PointParticle, DiagonalEV, CrosstermsEV, MeanFieldEV, QvdW };
	Ensemble fEnsemble;
	ModelType fModelType;
	ThermalModelParameters Parameters;
	double T, muB, muS, muQ, gammaS, gammaq, R;
	int B, Q, S;
};

// Base class for generating events from Thermal Model
class EventGeneratorBase
{
public:
	EventGeneratorBase() { m_THM = NULL; fCEAccepted = fCETotal = 0; }
	virtual ~EventGeneratorBase();// { }

	void ClearMomentumGenerators();

	void SetCollisionKineticEnergy(double ekin_) {
		SetCollisionCMSEnergy(sqrt(2.*xMath::mnucleon()*(ekin_ + 2. * xMath::mnucleon())));
	}
	void SetCollisionLabEnergy(double elab_) {
		SetCollisionCMSEnergy(sqrt(2.*xMath::mnucleon()*(elab_ + xMath::mnucleon())));
	}
	void SetCollisionCMSEnergy(double ssqrt_) {
		m_ssqrt = ssqrt_;
		m_ekin = m_ssqrt * m_ssqrt / 2. / xMath::mnucleon() - 2. * xMath::mnucleon();
		m_elab = xMath::mnucleon() + m_ekin;
		double plab = sqrt(m_elab*m_elab - xMath::mnucleon() * xMath::mnucleon());
		m_ycm = 0.5 * log((m_elab + xMath::mnucleon() + plab) / (m_elab + xMath::mnucleon() - plab));
	}

	//void SetConfiguration(const ThermalModelParameters& params, EventGeneratorConfiguration::Ensemble ensemble, EventGeneratorConfiguration::ModelType modeltype, ThermalParticleSystem *TPS, ThermalModelBase *original, ExcludedVolumeModel *exmod);
	void SetConfiguration(const ThermalModelParameters& params, EventGeneratorConfiguration::Ensemble ensemble, EventGeneratorConfiguration::ModelType modeltype, ThermalParticleSystem *TPS, ThermalModelBase *original, ThermalModelBase *THMEVVDW);

	std::vector<Acceptance::AcceptanceFunction>& GetAcceptance() { return m_acc; }
	virtual void ReadAcceptance(std::string accfolder);

	double getYcm() const { return m_ycm; }

	void PrepareMultinomials();
	std::vector<int> GenerateTotals() const;// { return vector<int>(); }
	std::vector<int> GenerateTotalsGCE() const;
	std::vector<int> GenerateTotalsSCE() const;
	std::vector<int> GenerateTotalsSCEnew() const;
	std::vector<int> GenerateTotalsSCESubVolume(double VolumeSC) const;
	std::vector<int> GenerateTotalsCE() const;

	virtual SimpleEvent GetEvent(bool PerformDecays = true) const;

	static int fCEAccepted, fCETotal;

protected:
	double m_ekin, m_ycm, m_ssqrt, m_elab;
	std::vector<Acceptance::AcceptanceFunction> m_acc;
	bool m_OnlyStable;
	EventGeneratorConfiguration m_Config;
	ThermalModelBase *m_THM;

	// Ideal gas densities for an interacting HRG
	std::vector<double> m_DensitiesIdeal;

	// Holds indexes and multinomial probabilities for efficient CE sampling
	std::vector< std::pair<double, int> > m_Baryons;
	std::vector< std::pair<double, int> > m_AntiBaryons;
	std::vector< std::pair<double, int> > m_StrangeMesons;
	std::vector< std::pair<double, int> > m_AntiStrangeMesons;
	std::vector< std::pair<double, int> > m_ChargeMesons;
	std::vector< std::pair<double, int> > m_AntiChargeMesons;

	std::vector<double> m_BaryonsProbs;
	std::vector<double> m_AntiBaryonsProbs;
	std::vector<double> m_StrangeMesonsProbs;
	std::vector<double> m_AntiStrangeMesonsProbs;
	std::vector<double> m_ChargeMesonsProbs;
	std::vector<double> m_AntiChargeMesonsProbs;

	double m_MeanB, m_MeanAB;
	double m_MeanSM, m_MeanASM;
	double m_MeanCM, m_MeanACM;

	static double m_LastWeight;
	static double m_LastLogWeight;
	
	std::vector<RandomGenerators::ParticleMomentumGenerator*>  	m_MomentumGens;
	std::vector<RandomGenerators::BreitWignerGenerator>				  m_BWGens;
};

#endif
