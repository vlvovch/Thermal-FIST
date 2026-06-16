/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGPhaseShifts/PhaseShiftChannel.h"
#include "HRGPhaseShifts/LightMesonPhaseShifts.h"
#include "HRGPhaseShifts/MesonBaryonPhaseShifts.h"

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

      // Sum the absolute light-quark / strangeness / charm content of the
      // channel's two scattering constituents (looked up in the list), so the
      // cluster scales with the correct non-equilibrium fugacity powers: e.g.
      // gammaq^4 for a pi-pi cluster, gammaq^3 gammaS for pi-K. Without this a
      // cluster would inherit single-meson powers from GetAbsQ() = 2 - |S| - |C|.
      // Returns false (leaving the defaults) if a constituent is not in the list.
      bool ConstituentAbsCharges(ThermalParticleSystem& TPS, const PhaseShiftChannel& ch,
                                 double& absQ, double& absS, double& absC) {
        absQ = absS = absC = 0.;
        const PhaseShiftConstituent* cons[2] = { &ch.a, &ch.b };
        for (int c = 0; c < 2; ++c) {
          int id = -1;
          for (std::map<int, long long>::const_iterator it = cons[c]->chargeStates.begin();
               it != cons[c]->chargeStates.end(); ++it) {
            int tid = TPS.PdgToId(it->second);
            if (tid >= 0) { id = tid; break; }
          }
          if (id < 0) return false;
          const ThermalParticle& p = TPS.Particle(id);
          absQ += p.GetAbsQ();
          absS += p.AbsoluteStrangeness();
          absC += p.AbsoluteCharm();
        }
        return true;
      }

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

      // Does the antiparticle of a member live in the SAME Iz-multiplet? True for
      // S=B=C=0 (e.g. pi-pi): then Iz<0 members are the antiparticles of Iz>0.
      // False for strangeness/baryon/charm-carrying channels (e.g. pi-K, S=+1):
      // every Iz is a distinct member and the antiparticles form a separate
      // (conjugate) multiplet.
      bool selfConjMultiplet(const PhaseShiftChannel& ch) {
        return ch.B == 0 && ch.S == 0 && ch.C == 0;
      }

      // 2Iz values of the members to create explicitly. Self-conjugate multiplet:
      // only Iz >= 0 (Iz < 0 are the antiparticles). Otherwise: every Iz in
      // [-I, +I] (the antiparticles are the conjugate multiplet).
      std::vector<int> memberTwoIz(const PhaseShiftChannel& ch) {
        std::vector<int> out;
        const int lo = selfConjMultiplet(ch) ? twoIzStart(ch.twoI) : -ch.twoI;
        for (int t = lo; t <= ch.twoI; t += 2) out.push_back(t);
        return out;
      }

      std::string IzLabel(int twoIz) {
        if (twoIz == 0) return "Iz=0";
        const std::string sign = (twoIz > 0) ? "+" : "-";
        const int a = std::abs(twoIz);
        return (a % 2 == 0) ? ("Iz=" + sign + std::to_string(a / 2))
                            : ("Iz=" + sign + std::to_string(a) + "/2");
      }

      // Spectroscopic ORBITAL letter (S/P/D/F/G) of a partial wave. Meson waves
      // have odd 2J+1 with J=l, so the orbital follows from 2J+1. Baryon waves
      // have even 2J+1 with J=l+-1/2, so 2J+1 does NOT fix l - the orbital is then
      // taken from the channel's excitation slot (which carries l, e.g. S31->0,
      // P31->1), giving correct labels for J=3/2 waves (P33/P13 are "P", not "J1").
      std::string WaveLabel(int twoJplus1, int excitation = 0) {
        const int l = (twoJplus1 % 2 == 1) ? (twoJplus1 - 1) / 2 : excitation;
        switch (l) {
          case 0: return "S"; case 1: return "P"; case 2: return "D";
          case 3: return "F"; case 4: return "G";
        }
        return "L" + std::to_string(l);
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
      // A channel may reuse real resonance PDG codes (subsumption by coincidence):
      // each member 2Iz maps to a real code, which the phase-shift density then
      // overrides (if present) or is created under (if absent, e.g. the kappa).
      if (!ch.memberPdg.empty()) {
        std::map<int, long long>::const_iterator it = ch.memberPdg.find(twoIz);
        if (it != ch.memberPdg.end()) return it->second;
        throw std::invalid_argument("PhaseShiftPdgId: memberPdg has no entry for 2Iz="
                                    + std::to_string(twoIz) + " in channel '" + ch.name + "'");
      }
      if (twoJplus1 < 1 || twoJplus1 > 9)
        throw std::invalid_argument("PhaseShiftPdgId: 2J+1 must be a single digit (J <= 4)");
      // Iz field: for a self-conjugate multiplet use 2|Iz| (Iz<0 is the
      // antiparticle = the negative code). For a charge/strangeness-carrying
      // multiplet every Iz is a distinct member, so encode the multiplet index
      // I+Iz = (2I+2Iz)/2 in 0..2I (the antiparticle stays the negative code).
      const long long z = selfConjMultiplet(ch)
        ? std::llabs((long long)twoIz)
        : (long long)((ch.twoI + twoIz) / 2);
      if (ch.excitation < 0 || ch.excitation > 99)
        throw std::invalid_argument("PhaseShiftPdgId: excitation (n n slot) must be 0..99");
      // 9 9 | F | (2I) | (Iz field) | n n | (2J+1) ; n n = excitation (orbital l
      // for baryon waves, 0 for mesons) occupies the tens+hundreds digits.
      return 99LL * 1000000LL
           + (long long)ch.family    * 100000LL
           + (long long)ch.twoI      * 10000LL
           +              z           * 1000LL
           + (long long)ch.excitation * 10LL
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
      for (size_t w = 0; w < waves.size(); ++w) fout << " " << WaveLabel(waves[w].twoJplus1, ch.excitation);
      fout << "\n";
      fout << "# pdgid name stable mass deg stat B Q S C |S| |C| width threshold\n";
      const std::vector<int> izs = memberTwoIz(ch);
      for (size_t ii = 0; ii < izs.size(); ++ii) {
        const int twoIz = izs[ii];
        const int Q = (twoIz + ch.B + ch.S + ch.C) / 2;
        for (size_t w = 0; w < waves.size(); ++w) {
          const int twoJp1 = waves[w].twoJplus1;
          const long long pdg = PhaseShiftPdgId(ch, twoIz, twoJp1);
          const std::string name = ch.name + "[" + IzLabel(twoIz) + "," + WaveLabel(twoJp1, ch.excitation) + "]";
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
      for (size_t w = 0; w < waves.size(); ++w) fout << " " << WaveLabel(waves[w].twoJplus1, ch.excitation);
      fout << "\n";
      const std::vector<int> izs = memberTwoIz(ch);
      for (size_t ii = 0; ii < izs.size(); ++ii) {
        const int twoIz = izs[ii];
        std::vector<std::pair<double, std::pair<long long, long long> > > dec = ChannelDecays(ch, twoIz);
        for (size_t w = 0; w < waves.size(); ++w) {
          const int twoJp1 = waves[w].twoJplus1;
          const long long pdg = PhaseShiftPdgId(ch, twoIz, twoJp1);
          fout << pdg << "  # " << ch.name << "[" << IzLabel(twoIz) << "," << WaveLabel(twoJp1, ch.excitation) << "]\n";
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
        const std::string suffix = "_" + WaveLabel(waves[w].twoJplus1, ch.excitation) + ".dat";
        WritePhaseShiftListFile(ch, one, base + "list-" + ch.name + suffix);
        WritePhaseShiftDecaysFile(ch, one, base + "decays-" + ch.name + suffix);
      }
    }

    std::vector<long long> AttachDensities(ThermalParticleSystem& TPS,
                                           const PhaseShiftChannel& ch,
                                           const std::vector<PhaseShiftPartialWave>& waves) {
      std::vector<long long> touched;
      // Constituent fugacity powers (gammaq^Nq gammaS^Ns gammaC^Nc), shared by
      // every member of the channel. If a constituent is missing from the list
      // the cluster keeps its (single-meson) defaults.
      double absQ = 0., absS = 0., absC = 0.;
      const bool haveAbs = ConstituentAbsCharges(TPS, ch, absQ, absS, absC);
      const std::vector<int> izs = memberTwoIz(ch);
      for (size_t ii = 0; ii < izs.size(); ++ii) {
        const int twoIz = izs[ii];
        const bool selfConj = selfConjugate(ch, twoIz);
        for (size_t w = 0; w < waves.size(); ++w) {
          const long long mag = PhaseShiftPdgId(ch, twoIz, waves[w].twoJplus1);
          long long pdgs[2] = { mag, -mag };
          const int npdg = selfConj ? 1 : 2;
          for (int s = 0; s < npdg; ++s) {
            if (TPS.PdgToId(pdgs[s]) < 0) continue;
            std::vector<PhaseShiftPartialWave> oneWave(1, waves[w]);
            ThermalParticle& cl = TPS.ParticleByPDG(pdgs[s]);
            cl.SetGeneralizedDensity(
              new PhaseShiftDensity(oneWave, ch.m1, ch.m2, ch.Mmax, ch.statistics, ch.quadratureNodes));
            if (haveAbs) {
              // Set strangeness/charm first (their setters recompute AbsQuark via
              // the single-hadron formula), then override AbsQuark last.
              cl.SetAbsoluteStrangeness(absS);
              cl.SetAbsoluteCharm(absC);
              cl.SetAbsoluteQuark(absQ);
            }
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
      const std::vector<int> izs = memberTwoIz(ch);
      for (size_t ii = 0; ii < izs.size(); ++ii) {
        const int twoIz = izs[ii];
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
          const std::string name = ch.name + "[" + IzLabel(twoIz) + "," + WaveLabel(twoJp1, ch.excitation) + "]";

          // Idempotent create-or-override: if the member is already in the list
          // (a synthetic cluster from a list file, or a real resonance reused via
          // memberPdg, e.g. the sigma / f0(980)), don't re-add it or touch its
          // decays - AttachDensities() below overrides only its thermal
          // contribution. If absent (e.g. the kappa, whose code is excluded from
          // the list) it is created here with the channel's isospin-CG decays.
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
      if (name == "pipi_I0") return PiPi_I0_Channel();
      if (name == "pipi_I0_f0980") return PiPi_I0_f0980_Channel();
      if (name == "pipi_I1") return PiPi_I1_Channel();
      if (name == "pipi_I1_F") return PiPi_I1_F_Channel();
      if (name == "piK_I32") return PiK_I32_Channel();
      if (name == "piK_I12") return PiK_I12_Channel();
      if (name == "piK_K892") return PiK_K892_Channel();
      if (name == "piN_Delta") return PiN_Delta_Channel();
      if (name == "piN_S31") return PiN_S31_Channel();
      if (name == "piN_S11") return PiN_S11_Channel();
      if (name == "piN_P31") return PiN_P31_Channel();
      if (name == "piN_P11") return PiN_P11_Channel();
      if (name == "piN_P13") return PiN_P13_Channel();
      if (name == "KN_S01") return KN_S01_Channel();
      if (name == "KN_P01") return KN_P01_Channel();
      if (name == "KN_P03") return KN_P03_Channel();
      if (name == "KN_S11") return KN_S11_Channel();
      if (name == "KN_P11") return KN_P11_Channel();
      if (name == "KN_P13") return KN_P13_Channel();
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
      // File paths are resolved relative to the config file's directory. A "-" in
      // the list and/or decay column means "no file": the channel reuses an
      // existing resonance (memberPdg), so nothing is added/loaded and the model
      // is just attached to it (subsumption by PDG coincidence).
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
        e.listFile   = (listF == "-") ? std::string() : resolvePath(dir, listF);
        e.decayFile  = (decF  == "-") ? std::string() : resolvePath(dir, decF);
        e.model      = model;
        entries.push_back(e);
      }

      // 1. Add the per-wave cluster members from each (unique) list file. Thermal-
      //    FIST generates the antiparticles; existing species are left untouched.
      //    Entries with no list file ("-") reuse an existing resonance: nothing added.
      std::set<std::string> addedLists;
      for (size_t i = 0; i < entries.size(); ++i)
        if (!entries[i].listFile.empty() && addedLists.insert(entries[i].listFile).second)
          TPS.AddParticlesToListFromFile(entries[i].listFile);
      TPS.FinalizeList();

      // 2. Attach each wave's isospin-CG decays from its (unique) decay file.
      //    Additive: only the listed cluster pdgs get decays; antiparticle decays
      //    are charge-conjugated. The rest of the list keeps its decays.
      std::set<std::string> loadedDecays;
      for (size_t i = 0; i < entries.size(); ++i) {
        if (entries[i].decayFile.empty()) continue;          // reused resonance keeps its decays
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

      // 3. Bind each wave's delta(M). File-based entries already have their cluster
      //    members in the list (added in step 1), so we just attach the density.
      //    "-" entries reuse real resonance codes (memberPdg): route them through
      //    the catalog AddPhaseShiftChannel, which create-or-overrides those codes
      //    (creates the kappa, overrides the sigma / f0(980)).
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
            + e.channel + ":" + WaveLabel(e.twoJplus1, ch.excitation) + "' does not match model '"
            + e.model + "' (a " + WaveLabel(w.twoJplus1, ch.excitation) + "-wave)");
        const std::vector<PhaseShiftPartialWave> one(1, w);
        std::vector<long long> t = e.listFile.empty()
          ? AddPhaseShiftChannel(TPS, ch, one)   // reuse real codes (create-or-override)
          : AttachDensities(TPS, ch, one);        // synthetic clusters from the list file
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

    int CountOverriddenResonances(const ThermalParticleSystem& TPS) {
      // A phase-shift density on a REUSED real resonance (memberPdg) rather than a
      // synthetic cluster: its PDG code is not in the synthetic 99-prefixed block
      // [99000000, 99999999]. For these the cheap enable/disable toggle is not
      // exact (a disabled override gives 0, not the pole-mass term), so the list
      // must be rebuilt to toggle them.
      int n = 0;
      for (int i = 0; i < TPS.ComponentsNumber(); ++i) {
        if (!dynamic_cast<PhaseShiftDensity*>(TPS.Particles()[i].GetGeneralizedDensity())) continue;
        const long long pdg = std::llabs((long long)TPS.Particles()[i].PdgId());
        if (pdg < 99000000LL || pdg >= 100000000LL) ++n;   // not a synthetic id
      }
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

    PhaseShiftChannel PiPi_I0_Channel() {
      PhaseShiftChannel ch;
      ch.name       = "pipi_I0";
      ch.family     = FamilyPiPi;
      ch.m1         = PionMass();
      ch.m2         = PionMass();
      ch.Mmax       = PiPiI0_Mmax();   // K-Kbar threshold (2 M_K): elastic sigma region
      ch.statistics = -1;              // bosonic cluster
      ch.twoI       = 0;               // I = 0 (isoscalar; single neutral member)
      ch.B = ch.S = ch.C = 0;
      ch.a.twoIsospin = 2;             // pion isospin triplet
      ch.a.chargeStates[+2] = 211;     // pi+
      ch.a.chargeStates[ 0] = 111;     // pi0
      ch.a.chargeStates[-2] = -211;    // pi-
      ch.b = ch.a;
      // The I=0 S-wave IS the sigma/f0(500): reuse its real PDG code (9000221).
      // The sigma is in the list with degeneracy 0, so the phase-shift density
      // overrides a zero baseline (no double counting), and it stays a decay
      // product. Toggling phase shifts off restores its (deg=0) list entry.
      ch.memberPdg[0] = 9000221;       // f0(500)/sigma
      ch.quadratureNodes = 64;
      return ch;
    }

    PhaseShiftChannel PiPi_I0_f0980_Channel() {
      // Same I=0 structure as the sigma channel, but covering the part of
      // delta_0^0 ABOVE the K-Kbar threshold (up to 1.42 GeV) - the f0(980).
      // It REUSES the real f0(980) PDG code (9010221): the phase-shift density is
      // attached to that existing resonance, overriding its thermal contribution
      // (subsumption by PDG coincidence). The f0(980) stays in the
      // list, so it is still a decay product of heavier resonances, and keeps its
      // own decays. (Above K-Kbar the I=0 wave is inelastic, so the elastic
      // Beth-Uhlenbeck term is an approximation.)
      PhaseShiftChannel ch = PiPi_I0_Channel();
      ch.name      = "pipi_I0_f0980";
      ch.Mmax      = PiPiI2_Mmax();      // 1.42 GeV (full Garcia-Martin upper range)
      ch.memberPdg.clear();
      ch.memberPdg[0] = 9010221;         // reuse the real f0(980) -> override it
      return ch;
    }

    PhaseShiftChannel PiPi_I1_Channel() {
      PhaseShiftChannel ch;
      ch.name       = "pipi_I1";
      ch.family     = FamilyPiPi;
      ch.m1         = PionMass();
      ch.m2         = PionMass();
      ch.Mmax       = PiPiI1_Mmax();   // K-Kbar threshold (2 M_K): elastic rho region
      ch.statistics = -1;              // bosonic cluster
      ch.twoI       = 2;               // I = 1 (the rho isospin triplet)
      ch.B = ch.S = ch.C = 0;
      ch.a.twoIsospin = 2;             // pion isospin triplet
      ch.a.chargeStates[+2] = 211;     // pi+
      ch.a.chargeStates[ 0] = 111;     // pi0
      ch.a.chargeStates[-2] = -211;    // pi-
      ch.b = ch.a;
      // The I=1 P-wave IS the rho(770): reuse its real PDG codes (the rho is in
      // the list, so the density overrides its contribution; it stays a decay
      // product). Self-conjugate multiplet: rho0 (Iz=0) + rho+ (Iz=+1, antiparticle
      // rho-). The P-wave's 2J+1=3 (in the spectral weight) carries the rho spin.
      ch.memberPdg[0]  = 113;          // rho0  (Iz=0,  Q=0)
      ch.memberPdg[+2] = 213;          // rho+  (Iz=+1, Q=1; antiparticle -213 = rho-)
      ch.quadratureNodes = 64;
      return ch;
    }

    PhaseShiftChannel PiPi_I1_F_Channel() {
      // pi-pi I=1 F-wave: non-resonant below 1.42 GeV (rho3(1690) is higher and
      // out of range), so unlike the rho there is no resonance to reuse - it is a
      // synthetic cluster. Different Mmax (1.42, treated as elastic) than the rho
      // P-wave (2 M_K), hence a separate channel.
      PhaseShiftChannel ch = PiPi_I1_Channel();
      ch.name = "pipi_I1_F";
      ch.Mmax = PiPiI2_Mmax();         // 1.42 GeV (elastic per the parametrization)
      ch.memberPdg.clear();            // synthetic cluster (no resonance to reuse)
      return ch;
    }

    // Shared structure of a pi-K channel (S = +1): pion triplet x kaon doublet.
    static PhaseShiftChannel PiK_ChannelBase() {
      PhaseShiftChannel ch;
      ch.family     = FamilyPiK;
      ch.m1         = PionMass();
      ch.m2         = KaonMass();
      ch.statistics = -1;            // bosonic cluster
      ch.B = ch.C = 0;
      ch.S = +1;                     // strangeness of the kaon
      ch.a.twoIsospin = 2;           // pion isospin triplet
      ch.a.chargeStates[+2] = 211;   // pi+
      ch.a.chargeStates[ 0] = 111;   // pi0
      ch.a.chargeStates[-2] = -211;  // pi-
      ch.b.twoIsospin = 1;           // kaon isospin doublet (S=+1)
      ch.b.chargeStates[+1] = 321;   // K+
      ch.b.chargeStates[-1] = 311;   // K0
      ch.quadratureNodes = 64;
      return ch;
    }

    PhaseShiftChannel PiK_I32_Channel() {
      PhaseShiftChannel ch = PiK_ChannelBase();
      ch.name = "piK_I32";
      ch.twoI = 3;                   // I = 3/2 (repulsive, non-resonant)
      ch.Mmax = PiK_I32_Mmax();      // elastic up to ~1.74 GeV
      // Non-resonant: subsumes no resonances.
      return ch;
    }

    PhaseShiftChannel PiK_I12_Channel() {
      PhaseShiftChannel ch = PiK_ChannelBase();
      ch.name = "piK_I12";
      ch.twoI = 1;                   // I = 1/2 (attractive, the kappa/K0*(700))
      ch.Mmax = PiK_I12_Mmax();      // elastic only below the K-eta threshold
      // The I=1/2 S-wave IS the kappa/K0*(700): reuse its real PDG codes. The
      // kappa is usually excluded from the HRG list (its contribution being the
      // I=1/2 phase shift), so these are created with the channel's isospin-CG
      // decays; if present they are overridden instead. (Subsumption by PDG
      // coincidence: 2Iz -> code.)
      ch.memberPdg[+1] = 9000321;    // K(0)*(700)+  (Iz=+1/2, Q=+1)
      ch.memberPdg[-1] = 9000311;    // K(0)*(700)0  (Iz=-1/2, Q=0)
      return ch;
    }

    PhaseShiftChannel PiK_K892_Channel() {
      // I=1/2 P-wave = the K*(892): elastic below the K-eta threshold, resonant.
      // Reuses the real K*(892) codes (in the list), overriding their contribution.
      // Separate channel from the kappa (piK_I12): same isospin but a different
      // wave reusing a different resonance.
      PhaseShiftChannel ch = PiK_ChannelBase();
      ch.name = "piK_K892";
      ch.twoI = 1;                   // I = 1/2
      ch.Mmax = PiK_I12_Mmax();      // elastic only below the K-eta threshold
      ch.memberPdg[+1] = 323;        // K*(892)+  (Iz=+1/2, Q=+1)
      ch.memberPdg[-1] = 313;        // K*(892)0  (Iz=-1/2, Q=0)
      return ch;
    }

    // ---- pi-N (meson-baryon, B = +1; not a self-conjugate multiplet) --------

    static PhaseShiftChannel PiN_ChannelBase() {
      PhaseShiftChannel ch;
      ch.family     = FamilyPiN;
      ch.m1         = PionMass();
      ch.m2         = NucleonMass();
      ch.statistics = +1;            // fermionic cluster (a baryon)
      ch.B = +1;                     // baryon number of the nucleon
      ch.S = ch.C = 0;
      ch.a.twoIsospin = 2;           // pion isospin triplet
      ch.a.chargeStates[+2] = 211;   // pi+
      ch.a.chargeStates[ 0] = 111;   // pi0
      ch.a.chargeStates[-2] = -211;  // pi-
      ch.b.twoIsospin = 1;           // nucleon isospin doublet (B=+1)
      ch.b.chargeStates[+1] = 2212;  // proton   (Iz=+1/2, Q=+1)
      ch.b.chargeStates[-1] = 2112;  // neutron  (Iz=-1/2, Q=0)
      ch.quadratureNodes = 64;
      return ch;
    }

    PhaseShiftChannel PiN_Delta_Channel() {
      // I=3/2 P-wave = the Delta(1232) (P33): resonant (branch-tracked through
      // 90 deg), elastic across the resonance up to sqrt(s) = 1.38 GeV. REUSES the
      // real Delta codes, overriding their thermal contribution while they stay in
      // the list as decay products (subsumption by PDG coincidence). A baryon
      // channel (B=+1), so every Iz is a distinct member and the antiparticles
      // (anti-Delta) form the conjugate multiplet.
      PhaseShiftChannel ch = PiN_ChannelBase();
      ch.name = "piN_Delta";
      ch.twoI = 3;                   // I = 3/2
      ch.Mmax = PiN_Delta_Mmax();    // elastic up to the matching point ~1.38 GeV
      ch.memberPdg[+3] = 2224;       // Delta(1232)++  (Iz=+3/2, Q=+2)
      ch.memberPdg[+1] = 2214;       // Delta(1232)+   (Iz=+1/2, Q=+1)
      ch.memberPdg[-1] = 2114;       // Delta(1232)0   (Iz=-1/2, Q=0)
      ch.memberPdg[-3] = 1114;       // Delta(1232)-   (Iz=-3/2, Q=-1)
      return ch;
    }

    // Non-resonant pi-N background waves (Roy-Steiner). Synthetic clusters (no
    // resonance to reuse in the elastic region); each is its own single-wave
    // channel so it gets its own elastic cutoff Mmax. The orbital l goes in the
    // synthetic-id excitation slot so same-(I,J) waves (e.g. S31/P31, S11/P11,
    // which share 2J+1=2) get distinct codes.
    static PhaseShiftChannel PiN_BackgroundChannel(const std::string& name,
                                                   int twoI, int orbitalL, double Mmax) {
      PhaseShiftChannel ch = PiN_ChannelBase();
      ch.name       = name;
      ch.twoI       = twoI;
      ch.excitation = orbitalL;      // n n slot = orbital l (S=0, P=1)
      ch.Mmax       = Mmax;
      return ch;                     // synthetic (no memberPdg)
    }

    PhaseShiftChannel PiN_S31_Channel() {  // I=3/2 S-wave (repulsive), elastic to 1.38
      return PiN_BackgroundChannel("piN_S31", 3, 0, PiN_Mmax_elastic());
    }
    PhaseShiftChannel PiN_S11_Channel() {  // I=1/2 S-wave (attractive), elastic to ~1.22
      return PiN_BackgroundChannel("piN_S11", 1, 0, PiN_Mmax_inelastic());
    }
    PhaseShiftChannel PiN_P31_Channel() {  // I=3/2 P-wave (J=1/2), elastic to ~1.22
      return PiN_BackgroundChannel("piN_P31", 3, 1, PiN_Mmax_inelastic());
    }
    PhaseShiftChannel PiN_P11_Channel() {  // I=1/2 P-wave (J=1/2, Roper tail), elastic to ~1.22
      return PiN_BackgroundChannel("piN_P11", 1, 1, PiN_Mmax_inelastic());
    }
    PhaseShiftChannel PiN_P13_Channel() {  // I=1/2 P-wave (J=3/2), elastic to 1.38
      return PiN_BackgroundChannel("piN_P13", 1, 1, PiN_Mmax_elastic());
    }

    // ---- K-N (kaon-nucleon, exotic S=+1; B=+1, fermionic; Gibbs-Arceo) -------

    static PhaseShiftChannel KN_ChannelBase() {
      PhaseShiftChannel ch;
      ch.family     = FamilyKN;
      ch.m1         = KaonMass();
      ch.m2         = NucleonMass();
      ch.statistics = +1;            // fermionic cluster (a baryon)
      ch.B = +1;                     // baryon number of the nucleon
      ch.S = +1;                     // strangeness of the kaon (exotic S=+1)
      ch.C = 0;
      ch.a.twoIsospin = 1;           // kaon isospin doublet (S=+1)
      ch.a.chargeStates[+1] = 321;   // K+  (Iz=+1/2)
      ch.a.chargeStates[-1] = 311;   // K0  (Iz=-1/2)
      ch.b.twoIsospin = 1;           // nucleon isospin doublet (B=+1)
      ch.b.chargeStates[+1] = 2212;  // proton  (Iz=+1/2)
      ch.b.chargeStates[-1] = 2112;  // neutron (Iz=-1/2)
      ch.quadratureNodes = 64;
      return ch;
    }

    // Non-resonant K-N waves: synthetic clusters, each its own single-wave channel,
    // all elastic to the K-pi-N threshold. The orbital l goes in the excitation slot
    // so same-(I,J) waves (S/P1/2 share 2J+1=2) get distinct codes.
    static PhaseShiftChannel KN_WaveChannel(const std::string& name, int twoI, int orbitalL) {
      PhaseShiftChannel ch = KN_ChannelBase();
      ch.name = name; ch.twoI = twoI; ch.excitation = orbitalL; ch.Mmax = KN_Mmax();
      return ch;
    }
    PhaseShiftChannel KN_S01_Channel() { return KN_WaveChannel("KN_S01", 0, 0); }
    PhaseShiftChannel KN_P01_Channel() { return KN_WaveChannel("KN_P01", 0, 1); }
    PhaseShiftChannel KN_P03_Channel() { return KN_WaveChannel("KN_P03", 0, 1); }
    PhaseShiftChannel KN_S11_Channel() { return KN_WaveChannel("KN_S11", 2, 0); }
    PhaseShiftChannel KN_P11_Channel() { return KN_WaveChannel("KN_P11", 2, 1); }
    PhaseShiftChannel KN_P13_Channel() { return KN_WaveChannel("KN_P13", 2, 1); }

  } // namespace PhaseShifts

} // namespace thermalfist
