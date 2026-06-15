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
#include <cctype>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <set>

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

      // Spectroscopic letter of a partial wave from its degeneracy 2J+1.
      std::string WaveLabel(int twoJplus1) {
        switch (twoJplus1) {
          case 1: return "S"; case 3: return "P"; case 5: return "D";
          case 7: return "F"; case 9: return "G";
        }
        return "J" + std::to_string((twoJplus1 - 1) / 2);
      }

      // Inverse: partial-wave token ("S".."G", or a numeric 2J+1) -> 2J+1.
      int WaveTokenToTwoJplus1(const std::string& tok) {
        if (tok.size() == 1) {
          switch (std::toupper(tok[0])) {
            case 'S': return 1; case 'P': return 3; case 'D': return 5;
            case 'F': return 7; case 'G': return 9;
          }
        }
        return std::atoi(tok.c_str());   // numeric 2J+1
      }

      // Directory part of a path (with trailing separator), or "" if none.
      std::string dirOf(const std::string& path) {
        std::string::size_type p = path.find_last_of("/\\");
        return (p == std::string::npos) ? std::string("") : path.substr(0, p + 1);
      }

      bool isAbsPath(const std::string& p) {
        return !p.empty() && (p[0] == '/' || p[0] == '\\'
                              || (p.size() > 1 && p[1] == ':'));   // C:\ ...
      }

      // Resolve a (possibly relative) path against a base directory.
      std::string resolvePath(const std::string& dir, const std::string& file) {
        return isAbsPath(file) ? file : dir + file;
      }

      // Read a phase-shift decays file (the new free format) into pdg -> decays.
      // Does not touch any ThermalParticleSystem, so it is additive by design.
      std::map<long long, ThermalParticle::ParticleDecaysVector>
        readClusterDecays(const std::string& file) {
        std::ifstream f(file.c_str());
        if (!f.is_open())
          throw std::runtime_error("AddPhaseShiftChannelsFromFile: cannot open decay file " + file);
        std::map<long long, ThermalParticle::ParticleDecaysVector> out;
        std::string raw, s;
        // next non-blank, comment-stripped line
        struct Next { std::ifstream& f; std::string& raw;
          bool operator()(std::string& s) {
            while (std::getline(f, raw)) {
              std::string::size_type h = raw.find('#');
              if (h != std::string::npos) raw = raw.substr(0, h);
              std::string::size_type a = raw.find_first_not_of(" \t\r\n");
              if (a == std::string::npos) continue;
              s = raw.substr(a);
              return true;
            }
            return false;
          }
        } next{f, raw};
        while (next(s)) {
          std::istringstream hs(s);
          long long pdg;
          if (!(hs >> pdg)) continue;
          std::string cs;
          if (!next(cs)) break;
          const int ndec = std::atoi(cs.c_str());
          ThermalParticle::ParticleDecaysVector decs;
          for (int i = 0; i < ndec; ++i) {
            std::string ds;
            if (!next(ds)) break;
            std::istringstream is(ds);
            ParticleDecayChannel c;
            if (!(is >> c.mBratio)) continue;
            long long d;
            while (is >> d) c.mDaughters.push_back(d);
            decs.push_back(c);
          }
          out[pdg] = decs;
        }
        return out;
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

    void WritePhaseShiftListFile(const PhaseShiftChannel& ch, const std::vector<PhaseShiftPartialWave>& waves,
                                 const std::string& listFile) {
      std::ofstream fout(listFile.c_str());
      if (!fout.is_open())
        throw std::runtime_error("WritePhaseShiftListFile: cannot open " + listFile);
      const double mass = ch.m1 + ch.m2;     // nominal; density override supplies thermodynamics
      fout << "# S-matrix phase-shift channel '" << ch.name << "', wave" << (waves.size() == 1 ? "" : "s") << ":";
      for (size_t w = 0; w < waves.size(); ++w) fout << " " << WaveLabel(waves[w].twoJplus1);
      fout << "\n";
      fout << "# pdgid name stable mass deg stat B Q S C |S| |C| width threshold\n";
      for (int twoIz = twoIzStart(ch.twoI); twoIz <= ch.twoI; twoIz += 2) {
        const int Q = (twoIz + ch.B + ch.S + ch.C) / 2;
        for (size_t w = 0; w < waves.size(); ++w) {
          const int twoJp1 = waves[w].twoJplus1;
          const long long pdg = PhaseShiftPdgId(ch, twoIz, twoJp1);
          const std::string name = ch.name + "[" + IzLabel(twoIz) + "," + WaveLabel(twoJp1) + "]";
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

    void WritePhaseShiftDecaysFile(const PhaseShiftChannel& ch, const std::vector<PhaseShiftPartialWave>& waves,
                                   const std::string& decaysFile) {
      std::ofstream fout(decaysFile.c_str());
      if (!fout.is_open())
        throw std::runtime_error("WritePhaseShiftDecaysFile: cannot open " + decaysFile);
      fout << "# S-matrix phase-shift channel '" << ch.name << "' decays (isospin Clebsch-Gordan), wave"
           << (waves.size() == 1 ? "" : "s") << ":";
      for (size_t w = 0; w < waves.size(); ++w) fout << " " << WaveLabel(waves[w].twoJplus1);
      fout << "\n";
      for (int twoIz = twoIzStart(ch.twoI); twoIz <= ch.twoI; twoIz += 2) {
        std::vector<std::pair<double, std::pair<long long, long long> > > dec = ChannelDecays(ch, twoIz);
        for (size_t w = 0; w < waves.size(); ++w) {
          const int twoJp1 = waves[w].twoJplus1;
          const long long pdg = PhaseShiftPdgId(ch, twoIz, twoJp1);
          fout << pdg << "  # " << ch.name << "[" << IzLabel(twoIz) << "," << WaveLabel(twoJp1) << "]\n";
          fout << dec.size() << "\n";
          for (size_t i = 0; i < dec.size(); ++i)
            fout << std::setprecision(10) << dec[i].first << "  "
                 << dec[i].second.first << " " << dec[i].second.second << "\n";
        }
      }
    }

    void WritePhaseShiftFiles(const PhaseShiftChannel& ch, const std::vector<PhaseShiftPartialWave>& waves,
                              const std::string& dir) {
      std::string base = dir;
      if (!base.empty() && base[base.size() - 1] != '/' && base[base.size() - 1] != '\\')
        base += "/";
      // One list + one decay file per partial wave, named list-<channel>_<wave>.dat
      // and decays-<channel>_<wave>.dat (wave = S/P/D/F/G).
      for (size_t w = 0; w < waves.size(); ++w) {
        const std::vector<PhaseShiftPartialWave> one(1, waves[w]);
        const std::string suffix = "_" + WaveLabel(waves[w].twoJplus1) + ".dat";
        WritePhaseShiftListFile(ch, one, base + "list-" + ch.name + suffix);
        WritePhaseShiftDecaysFile(ch, one, base + "decays-" + ch.name + suffix);
      }
    }

    std::vector<long long> AttachDensities(ThermalParticleSystem& TPS,
                                           const PhaseShiftChannel& ch,
                                           const std::vector<PhaseShiftPartialWave>& waves) {
      std::vector<long long> touched;
      for (int twoIz = twoIzStart(ch.twoI); twoIz <= ch.twoI; twoIz += 2) {
        const bool selfConj = selfConjugate(ch, twoIz);
        for (size_t w = 0; w < waves.size(); ++w) {
          const long long mag = PhaseShiftPdgId(ch, twoIz, waves[w].twoJplus1);
          long long pdgs[2] = { mag, -mag };
          const int npdg = selfConj ? 1 : 2;
          for (int s = 0; s < npdg; ++s) {
            if (TPS.PdgToId(pdgs[s]) < 0) continue;
            std::vector<PhaseShiftPartialWave> oneWave(1, waves[w]);
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
                                                const std::vector<PhaseShiftPartialWave>& waves) {
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

        for (size_t w = 0; w < waves.size(); ++w) {
          const int twoJp1 = waves[w].twoJplus1;
          const long long mag = PhaseShiftPdgId(ch, twoIz, twoJp1);
          const std::string name = ch.name + "[" + IzLabel(twoIz) + "," + WaveLabel(twoJp1) + "]";

          // Idempotent: if the cluster is already in the list (e.g. loaded from a
          // list-<name>.dat file), don't re-add it - AttachDensities() below will
          // bind the density to the existing entry.
          if (TPS.PdgToId(mag) < 0) {
            ThermalParticle part(/*Stable*/ false, name, mag, /*Deg*/ 1., ch.statistics, mass,
                                 ch.S, ch.B, Q, std::abs(ch.S), 0., mass, ch.C, std::abs(ch.C), 0);
            part.SetDecays(pdecays);
            part.SetDecaysOriginal(pdecays);
            TPS.AddParticle(part);
          }
          added.push_back(mag);

          if (!selfConj) {
            if (TPS.PdgToId(-mag) < 0) {
              ThermalParticle anti = TPS.ParticleByPDG(mag).GenerateAntiParticle();
              TPS.AddParticle(anti);
              // charge-conjugate the antiparticle decays (pi+ pi+ -> pi- pi-, ...)
              ThermalParticle::ParticleDecaysVector adec = TPS.GetDecaysFromAntiParticle(pdecays);
              TPS.ParticleByPDG(-mag).SetDecays(adec);
              TPS.ParticleByPDG(-mag).SetDecaysOriginal(adec);
            }
            added.push_back(-mag);
          }
        }
      }

      TPS.FinalizeList();                  // sort + PDG map + decay types
      AttachDensities(TPS, ch, waves);     // bind each wave's delta(M)
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

    std::vector<long long> AddPhaseShiftChannelsFromFile(ThermalParticleSystem& TPS,
                                                         const std::string& configFile) {
      std::ifstream fin(configFile.c_str());
      if (!fin.is_open())
        throw std::runtime_error("AddPhaseShiftChannelsFromFile: cannot open " + configFile);

      const std::string dir = dirOf(configFile);

      // One line per wave: "<channel>:<wave>  <listFile>  <decayFile>  <model>".
      // The list and decay files (per wave) carry the cluster members and their
      // isospin-CG decays; the model column gives the delta(M):
      //   <model> = <param>        analytic parametrization name (-> AnalyticWave)
      //           = tab:<file>     a tabulated (M, delta) table
      // File paths are resolved relative to the config file's directory.
      struct Entry {
        std::string channel; int twoJplus1;
        std::string listFile, decayFile, model;
      };
      std::vector<Entry> entries;
      std::string line;
      while (std::getline(fin, line)) {
        std::string::size_type hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        std::istringstream iss(line);
        std::string key, listF, decF, model;
        if (!(iss >> key >> listF >> decF >> model)) continue;   // skip blank/comment lines
        std::string::size_type colon = key.find(':');
        if (colon == std::string::npos)
          throw std::invalid_argument("AddPhaseShiftChannelsFromFile: column 1 must be "
                                      "'<channel>:<wave>', got '" + key + "'");
        Entry e;
        e.channel    = key.substr(0, colon);
        e.twoJplus1  = WaveTokenToTwoJplus1(key.substr(colon + 1));
        e.listFile   = resolvePath(dir, listF);
        e.decayFile  = resolvePath(dir, decF);
        e.model      = model;
        entries.push_back(e);
      }

      // 1. Add the per-wave cluster members from each (unique) list file. Thermal-
      //    FIST generates the antiparticles; existing species are left untouched.
      std::set<std::string> addedLists;
      for (size_t i = 0; i < entries.size(); ++i)
        if (addedLists.insert(entries[i].listFile).second)
          TPS.AddParticlesToListFromFile(entries[i].listFile);
      TPS.FinalizeList();

      // 2. Attach each wave's isospin-CG decays from its (unique) decay file.
      //    Additive: only the listed cluster pdgs get decays; antiparticle decays
      //    are charge-conjugated. The rest of the list keeps its decays.
      std::set<std::string> loadedDecays;
      for (size_t i = 0; i < entries.size(); ++i) {
        if (!loadedDecays.insert(entries[i].decayFile).second) continue;
        std::map<long long, ThermalParticle::ParticleDecaysVector> dm =
          readClusterDecays(entries[i].decayFile);
        for (std::map<long long, ThermalParticle::ParticleDecaysVector>::const_iterator it =
               dm.begin(); it != dm.end(); ++it) {
          if (TPS.PdgToId(it->first) >= 0) {
            TPS.ParticleByPDG(it->first).SetDecays(it->second);
            TPS.ParticleByPDG(it->first).SetDecaysOriginal(it->second);
          }
          if (TPS.PdgToId(-it->first) >= 0) {
            ThermalParticle::ParticleDecaysVector adec = TPS.GetDecaysFromAntiParticle(it->second);
            TPS.ParticleByPDG(-it->first).SetDecays(adec);
            TPS.ParticleByPDG(-it->first).SetDecaysOriginal(adec);
          }
        }
      }

      // 3. Attach the per-wave delta(M) density to each wave's clusters.
      std::vector<long long> all;
      std::vector<std::string> channelOrder;
      std::set<std::string> channelSeen;
      for (size_t i = 0; i < entries.size(); ++i) {
        const Entry& e = entries[i];
        PhaseShiftChannel ch = ChannelByName(e.channel);
        // A tabulated wave takes its 2J+1 from the wave key (col 1); an analytic
        // wave's name is wave-specific, so its 2J+1 must match the wave key.
        PhaseShiftPartialWave w = (e.model.substr(0, 4) == "tab:")
          ? TabulatedWave(e.twoJplus1, resolvePath(dir, e.model.substr(4)))
          : AnalyticWave(e.channel, e.model);
        if (w.twoJplus1 != e.twoJplus1)
          throw std::invalid_argument("AddPhaseShiftChannelsFromFile: wave key '"
            + e.channel + ":" + WaveLabel(e.twoJplus1) + "' does not match model '"
            + e.model + "' (a " + WaveLabel(w.twoJplus1) + "-wave)");
        std::vector<long long> t = AttachDensities(TPS, ch, std::vector<PhaseShiftPartialWave>(1, w));
        all.insert(all.end(), t.begin(), t.end());
        if (channelSeen.insert(e.channel).second) channelOrder.push_back(e.channel);
      }

      // 4. Remove any resonances subsumed by these channels (no-op if absent).
      for (size_t i = 0; i < channelOrder.size(); ++i)
        SubsumeResonances(TPS, ChannelByName(channelOrder[i]));

      // 5. Rebuild decay bookkeeping so feeddown includes the new clusters.
      TPS.FillDecayProperties();
      TPS.FillDecayThresholds();
      TPS.ProcessDecays();
      return all;
    }

    int SetPhaseShiftsEnabled(ThermalParticleSystem& TPS, bool enabled) {
      int n = 0;
      for (int i = 0; i < TPS.ComponentsNumber(); ++i) {
        GeneralizedDensity* gd = TPS.Particle(i).GetGeneralizedDensity();
        PhaseShiftDensity* psd = dynamic_cast<PhaseShiftDensity*>(gd);
        if (psd) { psd->SetEnabled(enabled); ++n; }
      }
      return n;
    }

    int CountPhaseShiftDensities(const ThermalParticleSystem& TPS) {
      int n = 0;
      for (int i = 0; i < TPS.ComponentsNumber(); ++i)
        if (dynamic_cast<PhaseShiftDensity*>(TPS.Particles()[i].GetGeneralizedDensity())) ++n;
      return n;
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
