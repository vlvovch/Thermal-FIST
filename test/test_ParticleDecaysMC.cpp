/*
 * Thermal-FIST package
 *
 * Copyright (c) 2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

// Tests of the Monte Carlo decay kinematics of ParticleDecaysMC
// for decays at, or within rounding errors of, the threshold.
//
// ManyBodyDecay() moves the mass of a mother below the threshold of the decay channel
// to the threshold plus 1e-7 GeV, and in a decay into four or more particles it draws
// the invariant mass of all daughters but the last one uniformly between its threshold
// and the mother mass minus the mass of the last daughter. The three-body decay
// that follows can then have no phase space, up to rounding errors. The rejection
// sampling of the invariant mass m12 in GetRandomThreeBodym12() never terminated in
// this case (seen in a production run when a uniform draw returned exactly zero).
// The tests below hang, up to the timeout of ctest, if this problem comes back.

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "HRGEventGenerator/ParticleDecaysMC.h"
#include "HRGEventGenerator/RandomGenerators.h"
#include "HRGEventGenerator/SimpleParticle.h"
#include "gtest/gtest.h"

using namespace thermalfist;

namespace {

  const double kPionMass = 0.13957;
  const double kKaonMass = 0.49368;

  /// Checks that a decay gives the expected number of products with finite momenta,
  /// and that it conserves the four-momentum of the mother (relative tolerance tol)
  void CheckDecayProducts(const SimpleParticle& mother, const std::vector<SimpleParticle>& products,
    size_t nproducts, double tol, const std::string& label)
  {
    ASSERT_EQ(products.size(), nproducts) << label;
    double E = 0., px = 0., py = 0., pz = 0.;
    for (size_t i = 0; i < products.size(); ++i) {
      const SimpleParticle& part = products[i];
      ASSERT_TRUE(std::isfinite(part.p0) && std::isfinite(part.px) && std::isfinite(part.py) && std::isfinite(part.pz))
        << label << ": product " << i << " has p = (" << part.p0 << ", " << part.px << ", " << part.py << ", " << part.pz << ")";
      E += part.p0;
      px += part.px;
      py += part.py;
      pz += part.pz;
    }
    EXPECT_NEAR(E, mother.p0, tol * mother.p0) << label;
    EXPECT_NEAR(px, mother.px, tol * mother.p0) << label;
    EXPECT_NEAR(py, mother.py, tol * mother.p0) << label;
    EXPECT_NEAR(pz, mother.pz, tol * mother.p0) << label;
  }

  /// Invariant mass of the products with indices [first, last)
  double InvariantMass(const std::vector<SimpleParticle>& products, size_t first, size_t last)
  {
    double E = 0., px = 0., py = 0., pz = 0.;
    for (size_t i = first; i < last; ++i) {
      E += products[i].p0;
      px += products[i].px;
      py += products[i].py;
      pz += products[i].pz;
    }
    return std::sqrt(E * E - px * px - py * py - pz * pz);
  }

  std::string Label(const std::string& what, double M)
  {
    std::ostringstream os;
    os.precision(17);
    os << what << ", M = " << M;
    return os.str();
  }

  /// Inverse of the tempering transformation applied by MTRand::randInt()
  /// to the words of the state of the Mersenne Twister
  MTRand::uint32 Untemper(MTRand::uint32 y)
  {
    const MTRand::uint32 mask32 = 0xffffffffUL;
    y &= mask32;
    y ^= y >> 18;
    y ^= (y << 15) & 0xefc60000UL;
    MTRand::uint32 t = y;
    for (int i = 0; i < 4; ++i)
      t = y ^ ((t << 7) & 0x9d2c5680UL);
    y = t & mask32;
    t = y;
    for (int i = 0; i < 2; ++i)
      t = y ^ (t >> 11);
    return t & mask32;
  }

  /// Makes the k-th upcoming 32-bit output of randgenMT equal to the given value.
  /// Used right after SetSeed(), when the state holds MTRand::N fresh words.
  void ForceUpcomingOutput(int k, MTRand::uint32 output)
  {
    MTRand::uint32 state[MTRand::SAVE];
    RandomGenerators::randgenMT.save(state);
    int left = static_cast<int>(state[MTRand::N]);
    ASSERT_LT(k, left);
    state[MTRand::N - left + k] = Untemper(output);
    RandomGenerators::randgenMT.load(state);
  }

  TEST(ParticleDecaysMC, ForcedGeneratorOutput)
  {
    const MTRand::uint32 outputs[] = { 0x0UL, 0xffffffffUL, 0x12345678UL, 0x80000001UL };
    for (int i = 0; i < 4; ++i) {
      RandomGenerators::SetSeed(1);
      MTRand reference(1);
      ForceUpcomingOutput(2, outputs[i]);
      EXPECT_EQ(RandomGenerators::randgenMT.randInt(), reference.randInt());
      EXPECT_EQ(RandomGenerators::randgenMT.randInt(), reference.randInt());
      EXPECT_EQ(RandomGenerators::randgenMT.randInt(), outputs[i]);
      reference.randInt();
      EXPECT_EQ(RandomGenerators::randgenMT.randInt(), reference.randInt());
    }
    RandomGenerators::SetSeed(1);
    ForceUpcomingOutput(0, 0x0UL);
    EXPECT_EQ(RandomGenerators::randgenMT.rand(), 0.);
    RandomGenerators::SetSeed(1);
    ForceUpcomingOutput(0, 0xffffffffUL);
    EXPECT_EQ(RandomGenerators::randgenMT.rand(), 1.);
    RandomGenerators::SetSeed(1);
    ForceUpcomingOutput(0, 0x0UL);
    EXPECT_GT(RandomGenerators::randgenMT.randDblExc(), 0.);
    RandomGenerators::SetSeed(1);
    ForceUpcomingOutput(0, 0xffffffffUL);
    EXPECT_LT(RandomGenerators::randgenMT.randDblExc(), 1.);
  }

  // A decay into three pions of a mother at the three-pion threshold, as computed
  // in ManyBodyDecay(), and slightly above it.
  // M = 0.41871000000000003 GeV is the mass of the three-pion state in the
  // stuck four-body decay: the uniform draw of its mass returned exactly zero.
  TEST(ParticleDecaysMC, ThreeBodyDecayAtThreshold)
  {
    const double m = kPionMass;
    std::vector<double> masses(3, m);
    std::vector<long long> pdgs;
    pdgs.push_back(211);
    pdgs.push_back(-211);
    pdgs.push_back(111);

    const double threshold = m + m + m;
    EXPECT_EQ(threshold, 0.41871000000000003);

    std::vector<double> motherMasses;
    double M = threshold;
    for (int i = 0; i <= 8; ++i) {
      motherMasses.push_back(M);
      M = std::nextafter(M, 1.);
    }
    const double offsets[] = { 1.e-15, 1.e-14, 1.e-13, 1.e-12, 1.e-11, 1.e-10, 1.e-9, 1.e-7, 1.e-3 };
    for (size_t i = 0; i < sizeof(offsets) / sizeof(offsets[0]); ++i)
      motherMasses.push_back(threshold + offsets[i]);

    RandomGenerators::SetSeed(1);
    for (size_t iM = 0; iM < motherMasses.size(); ++iM) {
      // At rest and moving
      SimpleParticle atRest(0., 0., 0., motherMasses[iM], 9999);
      SimpleParticle moving(0.3, -0.2, 1.5, motherMasses[iM], 9999);
      for (int rep = 0; rep < 20; ++rep) {
        CheckDecayProducts(atRest, ParticleDecaysMC::ManyBodyDecay(atRest, masses, pdgs), 3, 1.e-12, Label("at rest", motherMasses[iM]));
        CheckDecayProducts(moving, ParticleDecaysMC::ManyBodyDecay(moving, masses, pdgs), 3, 1.e-12, Label("moving", motherMasses[iM]));
      }
    }
  }

  // m12 sampling for degenerate or invalid inputs
  TEST(ParticleDecaysMC, ThreeBodym12DegenerateInputs)
  {
    const double m = kPionMass;
    const double threshold = m + m + m;
    const double eps = std::numeric_limits<double>::epsilon();

    // No phase space, at and below the threshold
    const double motherMasses[] = { threshold, std::nextafter(threshold, 1.), std::nextafter(threshold, 0.), threshold - 1.e-9 };
    for (int iM = 0; iM < 4; ++iM) {
      double M = motherMasses[iM];
      double fmax = ParticleDecaysMC::TernaryThreeBodym12Maximum(M, m, m, m);
      EXPECT_TRUE(std::isfinite(fmax)) << Label("maximum of the m12 density", M);
      EXPECT_GE(fmax, 0.) << Label("maximum of the m12 density", M);
      double m12 = ParticleDecaysMC::GetRandomThreeBodym12(M, m, m, m, 1.01 * fmax);
      EXPECT_TRUE(std::isfinite(m12)) << Label("m12", M);
      EXPECT_GE(m12, (m + m) * (1. - 4. * eps)) << Label("m12", M);
      EXPECT_LE(m12, std::max(m + m, M - m) * (1. + 4. * eps)) << Label("m12", M);
    }

    // Invalid maximum of the density, with a regular phase space
    const double M = 1.;
    const double maxima[] = { std::numeric_limits<double>::quiet_NaN(), 0., -1.,
      std::numeric_limits<double>::infinity(), 1.e300 };
    for (int i = 0; i < 5; ++i) {
      double m12 = ParticleDecaysMC::GetRandomThreeBodym12(M, m, m, m, maxima[i]);
      EXPECT_TRUE(std::isfinite(m12)) << "fm12max = " << maxima[i];
      EXPECT_GE(m12, m + m) << "fm12max = " << maxima[i];
      EXPECT_LE(m12, M - m) << "fm12max = " << maxima[i];
    }
  }

  // Four-pion decay with the uniform draw of the three-pion mass forced to the
  // smallest (zero) and the largest 32-bit output of the generator.
  // For the mother below the threshold, ManyBodyDecay() sets its mass to the threshold plus 1e-7 GeV.
  TEST(ParticleDecaysMC, FourBodyDecayWithExtremeDraws)
  {
    const double m = kPionMass;
    std::vector<double> masses(4, m);
    std::vector<long long> pdgs;
    pdgs.push_back(211);
    pdgs.push_back(-211);
    pdgs.push_back(211);
    pdgs.push_back(-211);
    const double threshold = m + m + m + m;

    const double motherMasses[] = { 1.0, threshold - 0.01, threshold };
    const MTRand::uint32 outputs[] = { 0x0UL, 0xffffffffUL };
    for (int iM = 0; iM < 3; ++iM) {
      for (int io = 0; io < 2; ++io) {
        SimpleParticle mother(0.2, -0.1, 0.4, motherMasses[iM], 9999);
        SimpleParticle adjusted = mother;
        if (adjusted.m < threshold) {
          adjusted.m = threshold + 1.e-7;
          adjusted.p0 = std::sqrt(adjusted.px * adjusted.px + adjusted.py * adjusted.py + adjusted.pz * adjusted.pz + adjusted.m * adjusted.m);
        }
        std::string label = Label(io == 0 ? "zero draw" : "largest draw", motherMasses[iM]);

        RandomGenerators::SetSeed(2);
        // Draws 0-2 shuffle the daughters, draw 3 gives the mass of the first three
        ForceUpcomingOutput(3, outputs[io]);
        std::vector<SimpleParticle> products = ParticleDecaysMC::ManyBodyDecay(mother, masses, pdgs);
        // Energy is not conserved exactly at the threshold, see below
        CheckDecayProducts(adjusted, products, 4, iM == 2 ? 1.e-6 : 1.e-12, label);
        if (products.size() != 4)
          continue;

        // The first product is the last daughter, the other three come from the three-pion state
        double m3pi = InvariantMass(products, 1, 4);
        if (io == 0)
          EXPECT_LT(m3pi - 3. * m, 1.e-6) << label << ": the three-pion state should be at its threshold";
        else
          EXPECT_LT(adjusted.m - m - m3pi, 1.e-6) << label << ": the three-pion state should be at its maximum mass";
      }
    }
  }

  // A decay into two particles at the threshold: the daughters are at rest
  // in the frame of the mother, without NaNs from rounding errors
  TEST(ParticleDecaysMC, TwoBodyDecayAtThreshold)
  {
    const double daughterMasses[][2] = { { kPionMass, kPionMass }, { kPionMass, kKaonMass }, { kKaonMass, kPionMass }, { 0.938272, kPionMass } };
    RandomGenerators::SetSeed(3);
    for (int ic = 0; ic < 4; ++ic) {
      std::vector<double> masses(daughterMasses[ic], daughterMasses[ic] + 2);
      std::vector<long long> pdgs(2, 211);
      double M = masses[0] + masses[1];
      for (int i = 0; i <= 4; ++i) {
        SimpleParticle atRest(0., 0., 0., M, 9999);
        SimpleParticle moving(-0.4, 0.1, 2.0, M, 9999);
        for (int rep = 0; rep < 20; ++rep) {
          CheckDecayProducts(atRest, ParticleDecaysMC::ManyBodyDecay(atRest, masses, pdgs), 2, 1.e-12, Label("at rest", M));
          CheckDecayProducts(moving, ParticleDecaysMC::ManyBodyDecay(moving, masses, pdgs), 2, 1.e-12, Label("moving", M));
          CheckDecayProducts(atRest, ParticleDecaysMC::TwoBodyDecay(atRest, masses[0], pdgs[0], masses[1], pdgs[1]), 2, 1.e-12, Label("TwoBodyDecay at rest", M));
        }
        M = std::nextafter(M, 10.);
      }
    }
  }

}
