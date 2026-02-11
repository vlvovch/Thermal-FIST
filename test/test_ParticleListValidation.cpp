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
#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

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

  // Verify branching ratios sum to ~1.0 for all unstable particles.
  // Raw BR sums in decays.dat may be slightly below 1.0 for poorly-measured
  // resonances; FIST normalizes them at runtime. We check that sums are
  // within [0.8, 1.001] to catch gross errors without flagging known gaps.
  TEST_F(ParticleListValidation, BranchingRatioSums) {
    for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
      const ThermalParticle& part = TPS->Particle(i);
      if (part.IsStable() || part.Decays().size() == 0)
        continue;

      double brSum = 0.0;
      for (size_t j = 0; j < part.Decays().size(); ++j) {
        brSum += part.Decays()[j].mBratio;
      }

      EXPECT_GE(brSum, 0.8)
        << "BR sum suspiciously low (" << brSum << ") for "
        << part.Name() << " (PDG " << part.PdgId() << ")";
      EXPECT_LE(brSum, 1.001)
        << "BR sum exceeds 1.0 (" << brSum << ") for "
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

}
