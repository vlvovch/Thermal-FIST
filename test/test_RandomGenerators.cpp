/*
 * Thermal-FIST package
 *
 * Copyright (c) 2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

// Tests of the random number generation in the event generator:
//  - RandomGenerators::RandomShuffle(), the random permutation used to
//    shuffle the order in which the momenta of the particles are sampled;
//  - reproducibility of the generated events for a fixed seed set with
//    RandomGenerators::SetSeed(). All the random numbers of the event generators
//    come from RandomGenerators::randgenMT, so the seed fixes the events.

#include <cmath>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalParticleSystem.h"
#include "HRGEventGenerator/RandomGenerators.h"
#include "HRGEventGenerator/SphericalBlastWaveEventGenerator.h"
#include "HRGEventGenerator/CylindricalBlastWaveEventGenerator.h"
#include "HRGEventGenerator/HypersurfaceSampler.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  // ---- RandomShuffle -----------------------------------------------------------

  // Empty and single-element ranges are left as they are and use no random numbers
  TEST(RandomShuffle, SmallSizes)
  {
    MTRand rangen(1), reference(1);

    std::vector<int> empty;
    RandomGenerators::RandomShuffle(empty.begin(), empty.end(), rangen);
    EXPECT_TRUE(empty.empty());

    std::vector<int> single(1, 42);
    RandomGenerators::RandomShuffle(single.begin(), single.end(), rangen);
    ASSERT_EQ(single.size(), 1u);
    EXPECT_EQ(single[0], 42);

    EXPECT_EQ(rangen.randInt(), reference.randInt());
  }

  // The result is a permutation of the input
  TEST(RandomShuffle, Permutation)
  {
    MTRand rangen(2);
    const int n = 1000;
    std::vector<int> v(n);
    for (int i = 0; i < n; ++i)
      v[i] = i;
    RandomGenerators::RandomShuffle(v.begin(), v.end(), rangen);

    int fixedPoints = 0;
    for (int i = 0; i < n; ++i)
      if (v[i] == i)
        fixedPoints++;
    // The mean number of fixed points of a random permutation is one
    EXPECT_LT(fixedPoints, 10);

    std::sort(v.begin(), v.end());
    for (int i = 0; i < n; ++i)
      ASSERT_EQ(v[i], i);
  }

  // All the 4! = 24 permutations of four elements are equally likely.
  // A Fisher-Yates shuffle with an off-by-one error in the range of the swap index
  // (which yields only the 6 cyclic permutations) or the naive shuffle
  // (swap index drawn from the whole range at each step) fail this test by far.
  TEST(RandomShuffle, UniformOverPermutations)
  {
    MTRand rangen(3);
    const int n = 4, nperm = 24, samples = 240000;
    std::vector<int> counts(nperm, 0);
    // Position of each element, to check the one-dimensional marginals
    std::vector< std::vector<int> > positions(n, std::vector<int>(n, 0));

    std::vector<int> identity(n);
    for (int i = 0; i < n; ++i)
      identity[i] = i;
    std::vector< std::vector<int> > permutations;
    {
      std::vector<int> perm = identity;
      do {
        permutations.push_back(perm);
      } while (std::next_permutation(perm.begin(), perm.end()));
    }
    ASSERT_EQ(static_cast<int>(permutations.size()), nperm);

    for (int is = 0; is < samples; ++is) {
      std::vector<int> v = identity;
      RandomGenerators::RandomShuffle(v.begin(), v.end(), rangen);
      int index = static_cast<int>(std::find(permutations.begin(), permutations.end(), v) - permutations.begin());
      ASSERT_LT(index, nperm);
      counts[index]++;
      for (int i = 0; i < n; ++i)
        positions[v[i]][i]++;
    }

    // Pearson chi^2 with 23 degrees of freedom; P(chi^2 > 70) ~ 1e-6
    double expected = static_cast<double>(samples) / nperm;
    double chi2 = 0.;
    for (int i = 0; i < nperm; ++i)
      chi2 += (counts[i] - expected) * (counts[i] - expected) / expected;
    EXPECT_LT(chi2, 70.);

    // Each element is equally likely at each position, 5 standard deviations
    double p = 1. / n;
    double sigma = std::sqrt(samples * p * (1. - p));
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        EXPECT_NEAR(positions[i][j], samples * p, 5. * sigma) << "element " << i << " at position " << j;
  }

  // By default the main generator randgenMT is used
  TEST(RandomShuffle, UsesMainGeneratorByDefault)
  {
    const unsigned int seed = 20261115;
    std::vector<int> v1(50), v2, v3;
    for (int i = 0; i < 50; ++i)
      v1[i] = i;
    v2 = v1;
    v3 = v1;

    RandomGenerators::SetSeed(seed);
    RandomGenerators::RandomShuffle(v1.begin(), v1.end());
    RandomGenerators::SetSeed(seed);
    RandomGenerators::RandomShuffle(v2.begin(), v2.end());
    MTRand rangen(seed);
    RandomGenerators::RandomShuffle(v3.begin(), v3.end(), rangen);

    EXPECT_EQ(v1, v2);
    EXPECT_EQ(v1, v3);
  }

  // ---- Reproducibility of the generated events ----------------------------------

  /// Pions, kaons, and nucleons with their antiparticles, all stable
  void FillLightHadrons(ThermalParticleSystem& TPS)
  {
    std::vector<ThermalParticle> parts;
    parts.push_back(ThermalParticle(true, "pi+", 211, 1., -1, 0.13957, 0, 0, 1));
    parts.push_back(ThermalParticle(true, "pi0", 111, 1., -1, 0.13498, 0, 0, 0));
    parts.push_back(ThermalParticle(true, "K+", 321, 1., -1, 0.49368, 1, 0, 1));
    parts.push_back(ThermalParticle(true, "p", 2212, 2., 1, 0.938272, 0, 1, 1));
    parts.push_back(ThermalParticle(true, "n", 2112, 2., 1, 0.939565, 0, 1, 0));
    TPS.SetTableFromVector(parts, true);
  }

  /// Grand-canonical ensemble of point particles at zero chemical potentials
  EventGeneratorConfiguration MakeConfig(double V)
  {
    EventGeneratorConfiguration config;
    config.fEnsemble = EventGeneratorConfiguration::GCE;
    config.fModelType = EventGeneratorConfiguration::PointParticle;
    config.CFOParameters = ThermalModelParameters(0.150, 0., 0., 0., 1., V, V, 0, 0, 0, 0);
    return config;
  }

  bool SameParticle(const SimpleParticle& p1, const SimpleParticle& p2)
  {
    return p1.PDGID == p2.PDGID
      && p1.p0 == p2.p0 && p1.px == p2.px && p1.py == p2.py && p1.pz == p2.pz
      && p1.r0 == p2.r0 && p1.rx == p2.rx && p1.ry == p2.ry && p1.rz == p2.rz;
  }

  /// Compares the events of two runs particle by particle (exact equality).
  /// Requires at least minParticles particles of at least two species,
  /// otherwise the order of the sampling would not matter.
  void ExpectIdenticalRuns(const std::vector<SimpleEvent>& run1, const std::vector<SimpleEvent>& run2,
    size_t minParticles)
  {
    ASSERT_EQ(run1.size(), run2.size());
    size_t compared = 0, differences = 0;
    std::string firstDifference;
    std::set<long long> species;
    for (size_t iev = 0; iev < run1.size(); ++iev) {
      ASSERT_EQ(run1[iev].Particles.size(), run2[iev].Particles.size()) << "event " << iev;
      for (size_t i = 0; i < run1[iev].Particles.size(); ++i) {
        const SimpleParticle& p1 = run1[iev].Particles[i];
        const SimpleParticle& p2 = run2[iev].Particles[i];
        species.insert(p1.PDGID);
        compared++;
        if (!SameParticle(p1, p2)) {
          if (differences == 0) {
            std::ostringstream os;
            os << "event " << iev << ", particle " << i << ": PDG " << p1.PDGID << " vs " << p2.PDGID
              << ", px " << p1.px << " vs " << p2.px;
            firstDifference = os.str();
          }
          differences++;
        }
      }
    }
    EXPECT_EQ(differences, 0u) << differences << " of " << compared << " particles differ, first: " << firstDifference;
    EXPECT_GE(compared, minParticles);
    EXPECT_GE(species.size(), 2u);
  }

  /// Generates nevents events with the generator after setting the seed
  std::vector<SimpleEvent> GenerateRun(const EventGeneratorBase& gen, unsigned int seed, int nevents, bool decays)
  {
    RandomGenerators::SetSeed(seed);
    std::vector<SimpleEvent> ret;
    for (int iev = 0; iev < nevents; ++iev)
      ret.push_back(gen.GetEvent(decays));
    return ret;
  }

  // Both runs are in the same process, so a second random number generator
  // that is not reset by SetSeed() would make them differ
  TEST(EventGeneratorFixedSeed, SphericalBlastWave)
  {
    ThermalParticleSystem TPS;
    FillLightHadrons(TPS);
    SphericalBlastWaveEventGenerator gen(&TPS, MakeConfig(2000.), 0.120, 0.5);

    std::vector<SimpleEvent> run1 = GenerateRun(gen, 11, 10, false);
    std::vector<SimpleEvent> run2 = GenerateRun(gen, 11, 10, false);
    ExpectIdenticalRuns(run1, run2, 1000);
  }

  // Full hadron list with resonance decays, new generator object for the second run
  TEST(EventGeneratorFixedSeed, CylindricalBlastWaveWithDecays)
  {
    ThermalParticleSystem TPS(std::string(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list.dat");

    std::vector<SimpleEvent> runs[2];
    for (int irun = 0; irun < 2; ++irun) {
      CylindricalBlastWaveEventGenerator gen(&TPS, MakeConfig(500.), 0.120, 0.5, 1.0, 1., 6.5);
      runs[irun] = GenerateRun(gen, 12, 10, true);
    }
    ExpectIdenticalRuns(runs[0], runs[1], 1000);
  }

  /// Static hypersurface: n x n volume elements at rest at a fixed proper time,
  /// on a grid in the transverse plane with the given spacing (in fm).
  /// The particles from an element all start at its position.
  /// The small space-like part of dsigma keeps the rejection weight of the momentum
  /// sampling below unity, without it rounding errors trigger warnings.
  ParticlizationHypersurface MakeHypersurface(int n, double spacing, double dV)
  {
    ParticlizationHypersurface ret;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        ParticlizationHypersurfaceElement elem;
        elem.tau = 10.;
        elem.x = spacing * (i - 0.5 * (n - 1));
        elem.y = spacing * (j - 0.5 * (n - 1));
        elem.eta = 0.;
        elem.dsigma[0] = dV;
        elem.u[0] = 1.;
        for (int mu = 1; mu < 4; ++mu) {
          elem.dsigma[mu] = 0.;
          elem.u[mu] = 0.;
        }
        elem.dsigma[1] = 0.05 * dV;
        elem.T = 0.150;
        elem.muB = elem.muQ = elem.muS = 0.;
        elem.edens = 0.3;
        elem.rhoB = 0.;
        for (int k = 0; k < 10; ++k)
          elem.pi[k] = 0.;
        elem.Pi = 0.;
        elem.P = 0.;
        ret.push_back(elem);
      }
    }
    return ret;
  }

  TEST(EventGeneratorFixedSeed, HypersurfaceSampler)
  {
    ThermalParticleSystem TPS;
    FillLightHadrons(TPS);
    ParticlizationHypersurface hypersurface = MakeHypersurface(8, 2., 2000. / 64.);
    HypersurfaceEventGenerator gen(&TPS, MakeConfig(2000.), &hypersurface);

    std::vector<SimpleEvent> run1 = GenerateRun(gen, 13, 10, false);
    std::vector<SimpleEvent> run2 = GenerateRun(gen, 13, 10, false);
    ExpectIdenticalRuns(run1, run2, 1000);
  }

  // Mesons, baryons, and antibaryons are shuffled separately here.
  // The hard-core radius of (anti)baryons is well below the grid spacing,
  // so that about nine (anti)baryons per event can always be placed.
  TEST(EventGeneratorFixedSeed, HypersurfaceSamplerEVHRG)
  {
    ThermalParticleSystem TPS;
    FillLightHadrons(TPS);
    ParticlizationHypersurface hypersurface = MakeHypersurface(8, 2., 2000. / 64.);
    HypersurfaceEventGeneratorEVHRG gen(&TPS, MakeConfig(2000.), &hypersurface);
    gen.SetExcludedVolume(1.);
    gen.SetBaryonRadius(0.3);

    std::vector<SimpleEvent> run1 = GenerateRun(gen, 14, 10, false);
    std::vector<SimpleEvent> run2 = GenerateRun(gen, 14, 10, false);
    ExpectIdenticalRuns(run1, run2, 1000);
  }

}
