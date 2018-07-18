#ifndef THERMALPARTICLESYSTEM_H
#define THERMALPARTICLESYSTEM_H

#include <map>
#include <vector>
#include <fstream>

#include "HRGBase/ThermalParticle.h"


/**
* Class which contains the whole particle table
*/
class ThermalParticleSystem
{
	public:
		enum ConservedCharge { BaryonCharge = 0, ElectricCharge = 1, StrangenessCharge = 2, CharmCharge = 3 };

		ThermalParticleSystem(std::string InputFile="", bool GenAntiP = true, double mcut = 1.e9);

		~ThermalParticleSystem(void);

		std::vector<ParticleDecay> GetDecaysFromAntiParticle(const std::vector<ParticleDecay> &Decays);

		void ProcessDecays() { FillResonanceDecays(); SeparateDecaysIntoWeakAndStrong(); FillResonanceWeakDecays(); }

		void FillDecayProperties();

		void FillDecayThresholds();

		void FillResonanceDecays();

		void FillResonanceWeakDecays();

		void GoResonance(int ind, int startind, double BR);

		void GoResonanceWeak(int ind, int startind, double BR);

		std::vector<double> GoResonanceDecayProbs(int ind, int goalind, bool firstdecay = false);

		std::vector<double> GoResonanceDecayProbsCharge(int ind, int nch, bool firstdecay = false);

		void LoadTable(std::string InputFile="", bool GenerateAntiParticles = true, double mcut = 1.e9);
		void LoadTable_OldFormat(std::ifstream &fin, bool GenerateAntiParticles = true, double mcut = 1.e9);
		void LoadTable_NewFormat(std::ifstream &fin, bool GenerateAntiParticles = true, double mcut = 1.e9);
		void WriteTableToFile(std::string OutputFile = "", bool WriteAntiParticles = false);

		void LoadDecays(std::string DecaysFile = "", bool GenerateAntiParticles = true);

		void NormalizeBranchingRatios();

		void RestoreBranchingRatios();

		void SetCalculationType(IdealGasFunctions::QStatsCalculationType type);
		void SetClusterExpansionOrder(int order);
		void SetResonanceWidthShape(ThermalParticle::ResonanceWidthShape shape);
		void SetResonanceWidthIntegrationType(ThermalParticle::ResonanceWidthIntegration type);
		ThermalParticle::ResonanceWidthIntegration ResonanceWidthIntegrationType() const { return m_ResonanceWidthIntegrationType; }

		std::string GetNameFromPDG(int pdgid);

		void SeparateDecaysIntoWeakAndStrong();

		bool hasBaryons() const { return (m_NumBaryons > 0); }
		bool hasCharged() const { return (m_NumCharged > 0); }
		bool hasStrange() const { return (m_NumStrange > 0); }
		bool hasCharmed() const { return (m_NumCharmed > 0); }

		const std::vector<ThermalParticle>& Particles() const { return m_Particles; }
		
		const ThermalParticle& Particle(int id)					const;// { return (id >= 0 && id < m_Particles.size()) ? m_Particles[id] : ThermalParticle(); }
					ThermalParticle& Particle(int id)							 ;// { return (id >= 0 && id < m_Particles.size()) ? m_Particles[id] : ThermalParticle(); }
		
					ThermalParticle& ParticleByPDG(int pdgid);// { return (m_PDGtoID.count(pdgid) > 0) ? m_Particles[m_PDGtoID[pdgid]] : ThermalParticle(); }
		
		int		PdgToId(int pdgid)		/*const*/ { return (m_PDGtoID.count(pdgid) > 0) ? m_PDGtoID[pdgid] : -1; }
		int		IdToPdg(int id)			const { return (id >= 0 && id < m_Particles.size()) ? m_Particles[id].PdgId() : 0; }

	private:
		std::vector<ThermalParticle>		m_Particles;
		std::map<int, int>							m_PDGtoID;
		std::vector<int>								m_IDtoPDG;
		int m_NumBaryons;
		int m_NumCharged;
		int m_NumStrange;
		int m_NumCharmed;

		int m_NumberOfParticles;

		ThermalParticle::ResonanceWidthIntegration m_ResonanceWidthIntegrationType;
};

namespace CuteHRGHelper {
	std::vector<std::string> split(const std::string &s, char delim);
}

#endif
