/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGPhaseShifts/PhaseShiftChannel.h"
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>

#include "HRGBase/ThermalParticleSystem.h"
#include "HRGBase/ThermalParticle.h"

namespace thermalfist {

  namespace PhaseShifts {

    namespace {

      double fact(int n) {
        double r = 1.;
        for (int i = 2; i <= n; ++i) r *= i;
        return r;
      }
      int ri(double x) { return (int)std::lround(x); }

      // Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M> (Racah formula).
      // Arguments may be half-integer.
      double clebschGordan(double j1, double m1, double j2, double m2, double J, double M) {
        if (std::fabs(m1 + m2 - M) > 1e-9) return 0.;
        if (J > j1 + j2 + 1e-9 || J < std::fabs(j1 - j2) - 1e-9) return 0.;
        if (std::fabs(m1) > j1 + 1e-9 || std::fabs(m2) > j2 + 1e-9 || std::fabs(M) > J + 1e-9) return 0.;

        double pref = std::sqrt((2. * J + 1.)
                       * fact(ri(j1 + j2 - J)) * fact(ri(j1 - j2 + J)) * fact(ri(-j1 + j2 + J))
                       / fact(ri(j1 + j2 + J + 1)))
                    * std::sqrt(fact(ri(j1 + m1)) * fact(ri(j1 - m1))
                       * fact(ri(j2 + m2)) * fact(ri(j2 - m2))
                       * fact(ri(J + M)) * fact(ri(J - M)));

        int klo = std::max(0, std::max(ri(j2 - J - m1), ri(j1 - J + m2)));
        int khi = std::min(ri(j1 + j2 - J), std::min(ri(j1 - m1), ri(j2 + m2)));
        double sum = 0.;
        for (int k = klo; k <= khi; ++k) {
          double denom = fact(k) * fact(ri(j1 + j2 - J) - k) * fact(ri(j1 - m1) - k)
                       * fact(ri(j2 + m2) - k) * fact(ri(J - j2 + m1) + k) * fact(ri(J - j1 - m2) + k);
          sum += ((k % 2) ? -1. : 1.) / denom;
        }
        return pref * sum;
      }

      // 2*Iz of the lowest non-negative member (0 for integer I, 1 for half-integer).
      int twoIzStart(int twoI) { return (twoI % 2 == 0) ? 0 : 1; }

      bool selfConjugate(const PhaseShiftChannel& ch, int twoIz) {
        return twoIz == 0 && ch.B == 0 && ch.S == 0 && ch.C == 0;
      }

      std::string IzLabel(int twoIz) {
        if (twoIz == 0) return "Iz=0";
        const std::string sign = (twoIz > 0) ? "+" : "-";
        const int a = std::abs(twoIz);
        return (a % 2 == 0) ? ("Iz=" + sign + std::to_string(a / 2))
                            : ("Iz=" + sign + std::to_string(a) + "/2");
      }

    } // anonymous namespace

    long long PhaseShiftPdgId(const PhaseShiftChannel& ch, int twoIz, int twoJplus1) {
      if (twoJplus1 < 1 || twoJplus1 > 9)
        throw std::invalid_argument("PhaseShiftPdgId: 2J+1 must be a single digit (J <= 4)");
      const long long z = std::llabs((long long)twoIz);  // 2|Iz|
      // 9 9 | F | (2I) | (2|Iz|) | n n | (2J+1) ; excitation slot n n = 00.
      return 99LL * 1000000LL
           + (long long)ch.family * 100000LL
           + (long long)ch.twoI   * 10000LL
           +              z        * 1000LL
           + (long long)twoJplus1;
    }

    std::vector<std::pair<double, std::pair<long long, long long> > >
      ChannelDecays(const PhaseShiftChannel& ch, int twoIz) {
      const double I  = ch.twoI / 2.0,        Iz = twoIz / 2.0;
      const double Ia = ch.a.twoIsospin / 2.0, Ib = ch.b.twoIsospin / 2.0;

      std::map<std::pair<long long, long long>, double> acc;
      for (std::map<int, long long>::const_iterator pa = ch.a.chargeStates.begin();
           pa != ch.a.chargeStates.end(); ++pa) {
        for (std::map<int, long long>::const_iterator pb = ch.b.chargeStates.begin();
             pb != ch.b.chargeStates.end(); ++pb) {
          if (pa->first + pb->first != twoIz) continue;       // 2m_a + 2m_b = 2Iz
          double cg = clebschGordan(Ia, pa->first / 2.0, Ib, pb->first / 2.0, I, Iz);
          double w = cg * cg;
          if (w < 1e-12) continue;
          long long d1 = pa->second, d2 = pb->second;
          if (d1 > d2) std::swap(d1, d2);                     // merge orderings (identical constituents)
          acc[std::make_pair(d1, d2)] += w;
        }
      }

      std::vector<std::pair<double, std::pair<long long, long long> > > out;
      for (std::map<std::pair<long long, long long>, double>::const_iterator it = acc.begin();
           it != acc.end(); ++it)
        out.push_back(std::make_pair(it->second, it->first));
      return out;
    }

    void WritePhaseShiftListFile(const PhaseShiftChannel& ch, const PhaseShiftModel& model,
                                 const std::string& listFile) {
      std::ofstream fout(listFile.c_str());
      if (!fout.is_open())
        throw std::runtime_error("WritePhaseShiftListFile: cannot open " + listFile);
      const double mass = ch.m1 + ch.m2;     // nominal; density override supplies thermodynamics
      fout << "# S-matrix phase-shift channel '" << ch.name << "' (model: " << model.name << ")\n";
      fout << "# pdgid name stable mass deg stat B Q S C |S| |C| width threshold\n";
      for (int twoIz = twoIzStart(ch.twoI); twoIz <= ch.twoI; twoIz += 2) {
        const int Q = (twoIz + ch.B + ch.S + ch.C) / 2;
        for (size_t w = 0; w < model.waves.size(); ++w) {
          const int twoJp1 = model.waves[w].twoJplus1;
          const long long pdg = PhaseShiftPdgId(ch, twoIz, twoJp1);
          const std::string name = ch.name + "[" + IzLabel(twoIz) + ",2J+1=" + std::to_string(twoJp1) + "]";
          fout << std::setw(12) << pdg << "  " << std::setw(28) << std::left << name << std::right
               << "  " << 0                      // stable=0 (decays to constituents)
               << "  " << std::setprecision(8) << mass
               << "  " << 1                      // degeneracy (2J+1 lives in the spectral weight)
               << "  " << ch.statistics
               << "  " << ch.B << "  " << Q << "  " << ch.S << "  " << ch.C
               << "  " << std::abs(ch.S) << "  " << std::abs(ch.C)
               << "  " << 0.0                    // width
               << "  " << mass                   // threshold
               << "\n";
        }
      }
    }

    void WritePhaseShiftDecaysFile(const PhaseShiftChannel& ch, const PhaseShiftModel& model,
                                   const std::string& decaysFile) {
      std::ofstream fout(decaysFile.c_str());
      if (!fout.is_open())
        throw std::runtime_error("WritePhaseShiftDecaysFile: cannot open " + decaysFile);
      fout << "# S-matrix phase-shift channel '" << ch.name << "' decays (isospin Clebsch-Gordan)\n";
      for (int twoIz = twoIzStart(ch.twoI); twoIz <= ch.twoI; twoIz += 2) {
        std::vector<std::pair<double, std::pair<long long, long long> > > dec = ChannelDecays(ch, twoIz);
        for (size_t w = 0; w < model.waves.size(); ++w) {
          const int twoJp1 = model.waves[w].twoJplus1;
          const long long pdg = PhaseShiftPdgId(ch, twoIz, twoJp1);
          fout << pdg << "  # " << ch.name << "[" << IzLabel(twoIz) << ",2J+1=" << twoJp1 << "]\n";
          fout << dec.size() << "\n";
          for (size_t i = 0; i < dec.size(); ++i)
            fout << std::setprecision(10) << dec[i].first << "  "
                 << dec[i].second.first << " " << dec[i].second.second << "\n";
        }
      }
    }

    void WritePhaseShiftFiles(const PhaseShiftChannel& ch, const PhaseShiftModel& model,
                              const std::string& dir) {
      std::string base = dir;
      if (!base.empty() && base[base.size() - 1] != '/' && base[base.size() - 1] != '\\')
        base += "/";
      WritePhaseShiftListFile(ch, model, base + "list-" + ch.name + ".dat");
      WritePhaseShiftDecaysFile(ch, model, base + "decays-" + ch.name + ".dat");
    }

    std::vector<long long> AttachDensities(ThermalParticleSystem& TPS,
                                           const PhaseShiftChannel& ch,
                                           const PhaseShiftModel& model) {
      std::vector<long long> touched;
      for (int twoIz = twoIzStart(ch.twoI); twoIz <= ch.twoI; twoIz += 2) {
        const bool selfConj = selfConjugate(ch, twoIz);
        for (size_t w = 0; w < model.waves.size(); ++w) {
          const long long mag = PhaseShiftPdgId(ch, twoIz, model.waves[w].twoJplus1);
          long long pdgs[2] = { mag, -mag };
          const int npdg = selfConj ? 1 : 2;
          for (int s = 0; s < npdg; ++s) {
            if (TPS.PdgToId(pdgs[s]) < 0) continue;
            std::vector<PhaseShiftPartialWave> oneWave(1, model.waves[w]);
            TPS.ParticleByPDG(pdgs[s]).SetGeneralizedDensity(
              new PhaseShiftDensity(oneWave, ch.m1, ch.m2, ch.Mmax, ch.statistics, ch.quadratureNodes));
            touched.push_back(pdgs[s]);
          }
        }
      }
      return touched;
    }

    std::vector<long long> SubsumeResonances(ThermalParticleSystem& TPS,
                                             const PhaseShiftChannel& ch) {
      std::vector<long long> removed;
      for (size_t k = 0; k < ch.subsumedPdg.size(); ++k) {
        const long long code = ch.subsumedPdg[k];
        if (code == 0) continue;
        long long signs[2] = { code, -code };
        const int nsign = (code == -code) ? 1 : 2;
        for (int s = 0; s < nsign; ++s) {
          const int id = TPS.PdgToId(signs[s]);
          if (id >= 0) { TPS.RemoveParticleAt(id); removed.push_back(signs[s]); }
        }
      }
      // Rebuild decay bookkeeping if anything was removed (removed resonances may
      // have been decay daughters of other species).
      if (!removed.empty()) {
        TPS.FillDecayProperties();
        TPS.FillDecayThresholds();
        TPS.ProcessDecays();
      }
      return removed;
    }

    std::vector<long long> AddPhaseShiftChannel(ThermalParticleSystem& TPS,
                                                const PhaseShiftChannel& ch,
                                                const PhaseShiftModel& model) {
      const double mass = ch.m1 + ch.m2;   // nominal; density override supplies thermodynamics
      std::vector<long long> added;

      // ---- add members (with isospin-CG decays) and their antiparticles ----
      for (int twoIz = twoIzStart(ch.twoI); twoIz <= ch.twoI; twoIz += 2) {
        const int Q = (twoIz + ch.B + ch.S + ch.C) / 2;

        std::vector<std::pair<double, std::pair<long long, long long> > > dec = ChannelDecays(ch, twoIz);
        ThermalParticle::ParticleDecaysVector pdecays;
        for (size_t i = 0; i < dec.size(); ++i) {
          ParticleDecayChannel c;
          c.mBratio = dec[i].first;
          c.mDaughters.push_back(dec[i].second.first);
          c.mDaughters.push_back(dec[i].second.second);
          pdecays.push_back(c);
        }
        const bool selfConj = selfConjugate(ch, twoIz);

        for (size_t w = 0; w < model.waves.size(); ++w) {
          const int twoJp1 = model.waves[w].twoJplus1;
          const long long mag = PhaseShiftPdgId(ch, twoIz, twoJp1);
          const std::string name = ch.name + "[" + IzLabel(twoIz) + ",2J+1=" + std::to_string(twoJp1) + "]";

          ThermalParticle part(/*Stable*/ false, name, mag, /*Deg*/ 1., ch.statistics, mass,
                               ch.S, ch.B, Q, std::abs(ch.S), 0., mass, ch.C, std::abs(ch.C), 0);
          part.SetDecays(pdecays);
          part.SetDecaysOriginal(pdecays);
          TPS.AddParticle(part);
          added.push_back(mag);

          if (!selfConj) {
            TPS.AddParticle(part.GenerateAntiParticle());
            // charge-conjugate the antiparticle decays (pi+ pi+ -> pi- pi-, ...)
            ThermalParticle::ParticleDecaysVector adec = TPS.GetDecaysFromAntiParticle(pdecays);
            TPS.ParticleByPDG(-mag).SetDecays(adec);
            TPS.ParticleByPDG(-mag).SetDecaysOriginal(adec);
            added.push_back(-mag);
          }
        }
      }

      TPS.FinalizeList();                  // sort + PDG map + decay types
      AttachDensities(TPS, ch, model);     // bind the chosen model's delta(M)
      SubsumeResonances(TPS, ch);          // graceful resonance removal (no-op if none)

      // rebuild decay bookkeeping so feeddown includes the new clusters
      TPS.FillDecayProperties();
      TPS.FillDecayThresholds();
      TPS.ProcessDecays();
      return added;
    }

    PhaseShiftChannel ChannelByName(const std::string& name) {
      if (name == "pipi_I2") return PiPi_I2_Channel();
      throw std::invalid_argument("ChannelByName: unknown channel '" + name + "'");
    }

    PhaseShiftModel ParsePhaseShiftModel(const std::string& spec) {
      if (spec.substr(0, 4) == "tab:") {
        std::vector<std::pair<int, std::string> > waveFiles;
        std::string rest = spec.substr(4);
        std::stringstream ss(rest);
        std::string item;
        while (std::getline(ss, item, ',')) {
          std::string::size_type eq = item.find('=');
          if (eq == std::string::npos)
            throw std::invalid_argument("ParsePhaseShiftModel: expected '<2J+1>=<file>' in '" + item + "'");
          int twoJp1 = std::atoi(item.substr(0, eq).c_str());
          std::string file = item.substr(eq + 1);
          waveFiles.push_back(std::make_pair(twoJp1, file));
        }
        return TabulatedModel(spec, waveFiles);
      }
      return AnalyticModel(spec);
    }

    std::vector<long long> AddPhaseShiftChannelsFromFile(ThermalParticleSystem& TPS,
                                                         const std::string& configFile) {
      std::ifstream fin(configFile.c_str());
      if (!fin.is_open())
        throw std::runtime_error("AddPhaseShiftChannelsFromFile: cannot open " + configFile);
      std::vector<long long> all;
      std::string line;
      while (std::getline(fin, line)) {
        std::string::size_type hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        std::istringstream iss(line);
        std::string chName, modelSpec;
        if (!(iss >> chName >> modelSpec)) continue;       // skip blank/comment lines
        PhaseShiftChannel ch = ChannelByName(chName);
        PhaseShiftModel  model = ParsePhaseShiftModel(modelSpec);
        std::vector<long long> a = AddPhaseShiftChannel(TPS, ch, model);
        all.insert(all.end(), a.begin(), a.end());
      }
      return all;
    }

    PhaseShiftChannel PiPi_I2_Channel() {
      PhaseShiftChannel ch;
      ch.name       = "pipi_I2";
      ch.family     = FamilyPiPi;
      ch.m1         = PionMass();
      ch.m2         = PionMass();
      ch.Mmax       = PiPiI2_Mmax();
      ch.statistics = -1;          // bosonic cluster
      ch.twoI       = 4;           // I = 2
      ch.B = ch.S = ch.C = 0;
      // both constituents are the pion isospin triplet: 2Iz -> pdg
      ch.a.twoIsospin = 2;
      ch.a.chargeStates[+2] = 211;   // pi+
      ch.a.chargeStates[ 0] = 111;   // pi0
      ch.a.chargeStates[-2] = -211;  // pi-
      ch.b = ch.a;
      // Purely repulsive, non-resonant: subsumes no resonances.
      ch.quadratureNodes = 64;
      return ch;
    }

  } // namespace PhaseShifts

} // namespace thermalfist
