/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2025 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  // Strip charge suffix from a particle name to get its isospin multiplet base.
  // "Sigma(1670)+" -> "Sigma(1670)", "Delta(1232)++" -> "Delta(1232)",
  // "pi0" -> "pi", "p" -> "p", "n" -> "n"
  std::string GetMultipletBase(const std::string& name) {
    // Order matters: try "++" before "+"
    if (name.size() >= 2 && name.substr(name.size()-2) == "++")
      return name.substr(0, name.size()-2);
    if (!name.empty() && name.back() == '+')
      return name.substr(0, name.size()-1);
    if (!name.empty() && name.back() == '-')
      return name.substr(0, name.size()-1);
    // Trailing "0" is a charge label only if the char before it is not a digit
    if (name.size() >= 2 && name.back() == '0' && !std::isdigit(name[name.size()-2]))
      return name.substr(0, name.size()-1);
    return name;
  }

  // Apply isospin grouping rules to a base name.
  // Groups particles that mix due to isospin CG coefficients.
  std::string ApplyIsospinGrouping(const std::string& base) {
    if (base == "p" || base == "n") return "N";
    // Ground-state Lambda (I=0) and Sigma (I=1) mix in Sigma-resonance
    // decays via isospin CG coefficients; group them together
    if (base == "Lambda" || base == "Sigma") return "LambdaSigma";
    return base;
  }

  // Map a daughter particle to its "channel type" base name.
  // Handles antiparticle detection, nucleon grouping (p/n -> "N"),
  // Lambda/Sigma grouping, and self-conjugate meson identification.
  std::string GetDaughterBase(long long pdgid, ThermalParticleSystem* TPS) {
    int idx = TPS->PdgToId(pdgid);
    if (idx >= 0) {
      const ThermalParticle& part = TPS->Particle(idx);
      // If this is an auto-generated antiparticle, use the particle's base
      if (part.IsAntiParticle()) {
        int idxP = TPS->PdgToId(-pdgid);
        if (idxP >= 0) {
          const ThermalParticle& p = TPS->Particle(idxP);
          std::string base = ApplyIsospinGrouping(GetMultipletBase(p.Name()));
          // Self-conjugate mesons (B=0, S=0, C=0): same multiplet
          if (p.BaryonCharge() == 0 && p.Strangeness() == 0 && p.Charm() == 0)
            return base;
          return "anti-" + base;
        }
      }
      return ApplyIsospinGrouping(GetMultipletBase(part.Name()));
    }
    // Antiparticle not in list: look up the particle by negating PID
    int idxP = TPS->PdgToId(-pdgid);
    if (idxP >= 0) {
      const ThermalParticle& p = TPS->Particle(idxP);
      std::string base = ApplyIsospinGrouping(GetMultipletBase(p.Name()));
      if (p.BaryonCharge() == 0 && p.Strangeness() == 0 && p.Charm() == 0)
        return base;
      return "anti-" + base;
    }
    // Elementary particles (gamma, leptons): use PID as name
    return std::to_string(pdgid);
  }

  // Build a canonical channel type string from a decay channel's daughters.
  std::string GetChannelType(const ParticleDecayChannel& ch,
                             ThermalParticleSystem* TPS) {
    std::vector<std::string> bases;
    for (size_t i = 0; i < ch.mDaughters.size(); ++i) {
      bases.push_back(GetDaughterBase(ch.mDaughters[i], TPS));
    }
    std::sort(bases.begin(), bases.end());
    std::string result;
    for (size_t i = 0; i < bases.size(); ++i) {
      if (i > 0) result += " + ";
      result += bases[i];
    }
    return result;
  }

  //============================================================================
  // Particle List Validation Tests (PDG2025)
  //============================================================================

  class ParticleListValidation : public ::testing::Test {
  protected:
    void SetUp() override {
      inputDir = std::string(ThermalFIST_INPUT_FOLDER);
      TPS = new ThermalParticleSystem(
        inputDir + "/list/PDG2025/list.dat",
        inputDir + "/list/PDG2025/decays.dat");
    }

    void TearDown() override {
      delete TPS;
    }

    std::string inputDir;
    ThermalParticleSystem* TPS;
  };

  // Verify branching ratios sum to exactly 1.0 for all unstable particles.
  // The generator normalizes all BRs to sum to 1.0 in decays.dat,
  // so any deviation beyond floating-point precision indicates a bug.
  TEST_F(ParticleListValidation, BranchingRatioSums) {
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      const ThermalParticle& part = TPS->Particle(i);
      if (part.IsStable() || part.Decays().size() == 0)
        continue;

      double brSum = 0.0;
      for (size_t j = 0; j < part.Decays().size(); ++j) {
        brSum += part.Decays()[j].mBratio;
      }

      EXPECT_NEAR(brSum, 1.0, 1.e-6)
        << "BR sum deviates from 1.0 (" << brSum << ") for "
        << part.Name() << " (PDG " << part.PdgId() << ")";
    }
  }

  // Verify charge conservation (B, Q, S, C) in every decay channel
  TEST_F(ParticleListValidation, DecayChargeConservation) {
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      const ThermalParticle& part = TPS->Particle(i);
      if (part.IsStable() || part.Decays().size() == 0)
        continue;

      EXPECT_TRUE(TPS->CheckDecayChargesConservation(i))
        << "Charge conservation violated for "
        << part.Name() << " (PDG " << part.PdgId() << ")";
    }
  }

  // Verify all unstable particles have decay channels specified
  TEST_F(ParticleListValidation, DecayChannelsSpecified) {
    EXPECT_TRUE(TPS->CheckDecayChannelsAreSpecified())
      << "Some unstable particles are missing decay channels";
  }

  // Verify expected particle counts and list properties
  TEST_F(ParticleListValidation, ParticleCountConsistency) {
    // PDG2025 list.dat has 246 entries, generating ~440 with antiparticles
    EXPECT_EQ(TPS->ComponentsNumber(), 440)
      << "Unexpected particle count for PDG2025 list.dat";

    EXPECT_TRUE(TPS->hasBaryons())
      << "PDG2025 list should contain baryons";
    EXPECT_TRUE(TPS->hasStrange())
      << "PDG2025 list should contain strange particles";
  }

  // Verify modular list files produce the same particle system
  TEST_F(ParticleListValidation, ModularListConsistency) {
    // Load modular hadrons + charm
    std::vector<std::string> listFiles = {
      inputDir + "/list/PDG2025/modular/list-hadrons.dat",
      inputDir + "/list/PDG2025/modular/list-charm.dat"
    };
    std::vector<std::string> decayFiles = {
      inputDir + "/list/PDG2025/modular/decays-hadrons.dat",
      inputDir + "/list/PDG2025/modular/decays-charm.dat"
    };
    std::set<std::string> flags;
    ThermalParticleSystem TPSmod(listFiles, decayFiles, flags, -1.);

    // Load monolithic list-withcharm.dat for comparison
    ThermalParticleSystem TPSmono(
      inputDir + "/list/PDG2025/list-withcharm.dat",
      inputDir + "/list/PDG2025/decays.dat");

    // Same total number of particles
    EXPECT_EQ(TPSmod.ComponentsNumber(), TPSmono.ComponentsNumber())
      << "Modular and monolithic particle counts differ";

    // Spot-check: pion (111), proton (2212), Lambda (3122)
    long long checkPdgs[] = { 111, 211, 2212, 3122, 421 };
    for (long long pdg : checkPdgs) {
      int idMod = TPSmod.PdgToId(pdg);
      int idMono = TPSmono.PdgToId(pdg);

      ASSERT_GE(idMod, 0) << "PDG " << pdg << " not found in modular list";
      ASSERT_GE(idMono, 0) << "PDG " << pdg << " not found in monolithic list";

      EXPECT_DOUBLE_EQ(TPSmod.Particle(idMod).Mass(),
                        TPSmono.Particle(idMono).Mass())
        << "Mass mismatch for PDG " << pdg;

      EXPECT_EQ(TPSmod.Particle(idMod).Decays().size(),
                TPSmono.Particle(idMono).Decays().size())
        << "Decay channel count mismatch for PDG " << pdg;
    }
  }

  // Verify PDG2020 list still loads correctly (backward compatibility)
  TEST_F(ParticleListValidation, PDG2020Compatibility) {
    EXPECT_NO_THROW({
      ThermalParticleSystem TPS2020(inputDir + "/list/PDG2020/list.dat");
      EXPECT_GT(TPS2020.ComponentsNumber(), 0)
        << "PDG2020 list loaded but has no particles";
      EXPECT_EQ(TPS2020.ComponentsNumber(), 434)
        << "Unexpected particle count for PDG2020 list.dat";
    }) << "Failed to load PDG2020 particle list";
  }

  // Verify isospin multiplet decay BR consistency for baryon resonances.
  // For each baryon isospin multiplet (e.g. Sigma(1670)+/0/-), the total
  // BR per "channel type" should be the same across all members.  Channel
  // types are formed by mapping daughter particles to their multiplet base
  // names (e.g. "Sigma + f(0)(500)", "Lambda + N + pi"), so that different
  // charge sub-channels are aggregated together.
  //
  // We restrict to baryons because meson multiplets have legitimate
  // asymmetries (EM decays for neutrals, KKbar phase space differences).
  //
  // This test would have caught the Sigma(1670)/Sigma(1915) bugs where
  // wrong charge states of daughters were assigned to some multiplet members.
  TEST_F(ParticleListValidation, MultipletBRConsistency) {
    // Group baryon particles into isospin multiplets by stripping charge labels
    std::map<std::string, std::vector<int>> multiplets;  // base -> list of indices
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      const ThermalParticle& p = TPS->Particle(i);
      // Only include unstable baryons (B > 0) with decay channels
      if (p.IsStable() || p.Decays().empty())
        continue;
      if (p.BaryonCharge() <= 0)
        continue;
      std::string base = GetMultipletBase(p.Name());
      multiplets[base].push_back(i);
    }

    int nChecked = 0;
    for (auto& kv : multiplets) {
      const std::string& base = kv.first;
      const std::vector<int>& members = kv.second;
      if (members.size() < 2)
        continue;

      // Build BR-by-channel-type map for each member
      std::vector<std::map<std::string, double>> memberBRs(members.size());
      for (size_t m = 0; m < members.size(); ++m) {
        const ThermalParticle& part = TPS->Particle(members[m]);
        for (size_t j = 0; j < part.Decays().size(); ++j) {
          std::string chType = GetChannelType(part.Decays()[j], TPS);
          memberBRs[m][chType] += part.Decays()[j].mBratio;
        }
      }

      // Collect all channel types
      std::set<std::string> allTypes;
      for (size_t m = 0; m < members.size(); ++m)
        for (auto& ct : memberBRs[m])
          allTypes.insert(ct.first);

      // Compare BR per channel type across members
      for (const std::string& chType : allTypes) {
        double maxBR = 0.0, minBR = 1.0;
        for (size_t m = 0; m < members.size(); ++m) {
          double br = 0.0;
          auto it = memberBRs[m].find(chType);
          if (it != memberBRs[m].end())
            br = it->second;
          maxBR = std::max(maxBR, br);
          minBR = std::min(minBR, br);
        }

        // 2% absolute tolerance — allows for CG coefficient rounding
        // and mass-difference effects within multiplet.
        // Only flag if the max BR is significant (>1%)
        if (maxBR > 0.01) {
          EXPECT_NEAR(maxBR, minBR, 0.02)
            << "BR mismatch in " << base << " multiplet for channel type ["
            << chType << "]: ";
        }
      }
      ++nChecked;
    }

    // Sanity check: we should have checked a reasonable number of multiplets
    // (N, Delta, Sigma, Xi families — expect ~30+)
    EXPECT_GT(nChecked, 20)
      << "Too few baryon multiplets checked — possible issue with name parsing";
  }

}
