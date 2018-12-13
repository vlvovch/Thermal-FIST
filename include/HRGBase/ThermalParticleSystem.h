/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALPARTICLESYSTEM_H
#define THERMALPARTICLESYSTEM_H

#include <map>
#include <vector>
#include <fstream>

#include "HRGBase/ThermalParticle.h"

namespace thermalfist {

  /**
  * Class which contains the whole particle list
  */
  class ThermalParticleSystem
  {
  public:
    enum ConservedCharge { BaryonCharge = 0, ElectricCharge = 1, StrangenessCharge = 2, CharmCharge = 3 };

    ThermalParticleSystem(std::string InputFile = "", bool GenAntiP = true, double mcut = 1.e9);

    ~ThermalParticleSystem(void);

    std::vector<ParticleDecay> GetDecaysFromAntiParticle(const std::vector<ParticleDecay> &Decays);

    void ProcessDecays();// { FillResonanceDecays(); SeparateDecaysIntoWeakAndStrong(); FillResonanceWeakDecays(); }

    void FillDecayProperties();

    void FillDecayThresholds();

    void FillResonanceDecays();

    //void FillResonanceWeakDecays();
    void FillResonanceDecaysByFeeddown();

    void GoResonance(int ind, int startind, double BR);

    //void GoResonanceWeak(int ind, int startind, double BR);
    void GoResonanceByFeeddown(int ind, int startind, double BR, Feeddown::Type feeddown);

    std::vector<double> GoResonanceDecayProbs(int ind, int goalind, bool firstdecay = false);

    std::vector<double> GoResonanceDecayProbsCharge(int ind, int nch, bool firstdecay = false);

    std::vector< std::pair<double, std::vector<int> > > GoResonanceDecayDistributions(int ind, bool firstdecay = false);

    void LoadTable(std::string InputFile = "", bool GenerateAntiParticles = true, double mcut = 1.e9);
    void LoadTable_OldFormat(std::ifstream &fin, bool GenerateAntiParticles = true, double mcut = 1.e9);
    void LoadTable_NewFormat(std::ifstream &fin, bool GenerateAntiParticles = true, double mcut = 1.e9);
    void SetTableFromVector(const std::vector<ThermalParticle> &part_in, bool GenerateAntiParticles = true);
    void WriteTableToFile(std::string OutputFile = "", bool WriteAntiParticles = false);

    void LoadDecays(std::string DecaysFile = "", bool GenerateAntiParticles = true);
    void ReadDecays_OldFormat(std::ifstream &fin);
    void ReadDecays_NewFormat(std::ifstream &fin);
    void WriteDecaysToFile(std::string OutputFile = "", bool WriteAntiParticles = false);

    void NormalizeBranchingRatios();

    void RestoreBranchingRatios();

    void SetCalculationType(IdealGasFunctions::QStatsCalculationType type);
    void SetClusterExpansionOrder(int order);
    void SetResonanceWidthShape(ThermalParticle::ResonanceWidthShape shape);
    void SetResonanceWidthIntegrationType(ThermalParticle::ResonanceWidthIntegration type);
    ThermalParticle::ResonanceWidthIntegration ResonanceWidthIntegrationType() const { return m_ResonanceWidthIntegrationType; }

    std::string GetNameFromPDG(int pdgid);

    //void SeparateDecaysIntoWeakAndStrong();

    bool hasBaryons() const { return (m_NumBaryons > 0); }
    bool hasCharged() const { return (m_NumCharged > 0); }
    bool hasStrange() const { return (m_NumStrange > 0); }
    bool hasCharmed() const { return (m_NumCharmed > 0); }

    const std::vector<ThermalParticle>& Particles() const { return m_Particles; }

    const ThermalParticle& Particle(int id)          const;// { return (id >= 0 && id < m_Particles.size()) ? m_Particles[id] : ThermalParticle(); }
    ThermalParticle& Particle(int id);// { return (id >= 0 && id < m_Particles.size()) ? m_Particles[id] : ThermalParticle(); }

    ThermalParticle& ParticleByPDG(int pdgid);// { return (m_PDGtoID.count(pdgid) > 0) ? m_Particles[m_PDGtoID[pdgid]] : ThermalParticle(); }

    int    PdgToId(int pdgid)    /*const*/ { return (m_PDGtoID.count(pdgid) > 0) ? m_PDGtoID[pdgid] : -1; }
    int    IdToPdg(int id)      const { return (id >= 0 && id < m_Particles.size()) ? m_Particles[id].PdgId() : 0; }

    void FillPdgMap();
    void FinalizeList();

    void AddParticle(const ThermalParticle & part);
    void RemoveParticleAt(int ind);

    /**
    * Checks whether cumulative charges (B, Q, S, C) of decay products match those of decaying particle with index ind
    */
    bool CheckDecayChargesConservation(int ind) const;

    bool operator==(const ThermalParticleSystem &rhs) const;
    bool operator!=(const ThermalParticleSystem &rhs) const { return !(*this == rhs); }

    static ParticleDecay::DecayType DecayTypeByParticleType(const ThermalParticle &part);
  private:
    std::vector<ThermalParticle>    m_Particles;
    std::map<int, int>              m_PDGtoID;
    int m_NumBaryons;
    int m_NumCharged;
    int m_NumStrange;
    int m_NumCharmed;

    int m_NumberOfParticles;

    ThermalParticle::ResonanceWidthIntegration m_ResonanceWidthIntegrationType;


    // Map for DP-based calculations of decay distributions
    std::vector< std::vector< std::pair<double, std::vector<int> > > > m_DecayDistributionsMap;
  };

  namespace CuteHRGHelper {
    std::vector<std::string> split(const std::string &s, char delim);
    void cutDecayDistributionsVector(std::vector<std::pair<double, std::vector<int> > > &vect, int maxsize = 1000);
  }

} // namespace thermalfist

#endif
