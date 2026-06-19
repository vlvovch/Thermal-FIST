/*
 * Thermal-FIST package
 *
 * Copyright (c) 2024-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef PHASESHIFTCHANNEL_H
#define PHASESHIFTCHANNEL_H

/**
 * \file PhaseShiftChannel.h
 *
 * \brief Physical description of an S-matrix two-body scattering channel (an
 *        isospin multiplet) plus the helpers that turn it into a Thermal-FIST
 *        list module: a particle table (list-<name>.dat), a decays file
 *        (decays-<name>.dat) and the routine that attaches a (selectable)
 *        per-wave delta(M) to each cluster member after loading.
 *
 * The channel is model-independent: it fixes the constituents, quantum numbers
 * and the decay (isospin-Clebsch-Gordan) structure. The dynamics, delta(M), come
 * from a set of per-wave models (std::vector<PhaseShiftPartialWave>, one entry
 * per partial wave), so each wave can use a different parametrization or a
 * tabulated phase shift.
 */

#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include "HRGPhaseShifts/PhaseShiftModel.h"

namespace thermalfist {

  class ThermalParticleSystem; // forward declaration

  namespace PhaseShifts {

    /**
     * \brief Flavour-family digit of the synthetic PDG id (see PhaseShiftPdgId).
     */
    enum ChannelFamily {
      FamilyPiPi  = 1, ///< pi pi
      FamilyPiK   = 2, ///< pi K
      FamilyKK    = 3, ///< K Kbar
      FamilyPiN   = 4, ///< pi N
      FamilyPiEta = 5, ///< pi eta
      FamilyKN    = 6  ///< K N (exotic, S=+1)
    };

    /**
     * \brief One scattering constituent: its isospin and the PDG codes of its
     *        charge (Iz) states, used to build the isospin-CG decay channels.
     */
    struct PhaseShiftConstituent {
      int twoIsospin;                          ///< 2 * isospin of the constituent
      std::map<int, long long> chargeStates;   ///< 2*Iz -> PDG code of that charge state
    };

    /**
     * \brief Physical (model-independent) definition of a scattering channel.
     *
     * The electric charge of each member follows from Q = Iz + (B + S + C)/2.
     */
    struct PhaseShiftChannel {
      std::string name;                ///< label, e.g. "pipi_I2"
      int    family;                   ///< ChannelFamily flavour digit
      double m1, m2;                   ///< constituent masses [GeV]
      double Mmax;                     ///< integration cutoff [GeV]
      int    statistics;               ///< cluster stat (-1 Bose, +1 Fermi, 0 Boltzmann)
      int    twoI;                     ///< 2 * total isospin
      int    B, S, C;                  ///< cluster baryon / strangeness / charm
      int    excitation;               ///< the synthetic-id "n n" slot (0 by default).
                                       ///< Distinguishes channels with the same
                                       ///< (family, 2I, Iz, 2J+1): for baryon waves
                                       ///< (J = l +- 1/2) it carries the orbital l, so
                                       ///< e.g. S31 (l=0) and P31 (l=1) - which share
                                       ///< 2J+1=2 - get distinct codes.
      PhaseShiftConstituent a, b;      ///< the two scattering constituents (for decays)
      std::vector<long long> subsumedPdg; ///< PDG codes of resonances this channel replaces
      int    quadratureNodes;          ///< Gauss-Legendre nodes for the q-integral
      std::map<int, long long> memberPdg; ///< if non-empty, each member REUSES a real PDG
                                       ///< code (2Iz -> code) instead of a synthetic id:
                                       ///< "subsumption by PDG coincidence". If that code
                                       ///< is already in the list the phase-shift density
                                       ///< overrides its thermal contribution (and it
                                       ///< keeps its decays / stays a decay product); if
                                       ///< absent it is created (e.g. the kappa) with the
                                       ///< channel's isospin-CG decays. Must cover every
                                       ///< member 2Iz of the channel.
      double resMass;                  ///< for memberPdg (subsumed-resonance) channels: the
                                       ///< real resonance mass [GeV] written to the generated
                                       ///< list file and used when the code is absent and must
                                       ///< be created (so the created entry, and the pole-mass
                                       ///< fallback if phase shifts are toggled off, are
                                       ///< physical). <= 0 -> fall back to the cluster mass
                                       ///< m1+m2. Ignored for synthetic clusters.
      int    resDeg;                   ///< real spin degeneracy 2J+1 of the subsumed resonance
                                       ///< (may legitimately be 0, e.g. the broad sigma kept
                                       ///< as a phase-shift-only placeholder). Ignored for
                                       ///< synthetic clusters (which use deg 1).
      double resWidth;                 ///< real width [GeV] of the subsumed resonance for the
                                       ///< created/fallback entry (0 default; ignored when the
                                       ///< phase-shift density overrides the thermodynamics).

      PhaseShiftChannel()
        : family(FamilyPiPi), m1(0.), m2(0.), Mmax(0.), statistics(-1),
          twoI(0), B(0), S(0), C(0), excitation(0), quadratureNodes(64),
          resMass(0.), resDeg(0), resWidth(0.) {}
    };

    /**
     * \brief Synthetic PDG id (magnitude) of one partial wave of a channel member,
     *        following the spirit of the PDG MC numbering.
     *
     * Layout (8 digits): "9 9 F (2I) (2|Iz|) n n (2J+1)". The LAST digit is the
     * wave's spin degeneracy 2J+1 (as in PDG); the two "n n" digits are the
     * channel's excitation index (0 by default), which distinguishes two channels
     * with the same (family, 2I, Iz, 2J+1) - e.g. sigma (00) vs f0(980) (01) in
     * pi-pi I=0; electric charge, baryon number and strangeness are carried by the
     * ThermalParticle fields, not packed into the code; the 8-digit 99-prefix
     * stays clear of all real PDG codes. The antiparticle uses the negative.
     *
     * \param ch          Channel (provides family and twoI).
     * \param twoIz        2*Iz of the member (magnitude uses |2Iz|).
     * \param twoJplus1    2J+1 of the wave (1..9); becomes the last digit.
     */
    long long PhaseShiftPdgId(const PhaseShiftChannel& ch, int twoIz, int twoJplus1);

    /**
     * \brief Isospin-Clebsch-Gordan decay channels of the (I, Iz) member into its
     *        constituent charge states. Independent of the partial wave.
     *
     * \return List of (branching ratio, {daughter PDG, daughter PDG}). For
     *         identical constituents the two orderings are merged; the branchings
     *         sum to 1.
     */
    std::vector<std::pair<double, std::pair<long long, long long> > >
      ChannelDecays(const PhaseShiftChannel& ch, int twoIz);

    /**
     * \brief Write the particle table (list-<name>.dat) for the channel/model.
     *
     * One row per (Iz >= 0, partial wave) member, marked unstable; antiparticles
     * are generated by Thermal-FIST at load time.
     */
    void WritePhaseShiftListFile(const PhaseShiftChannel& ch, const std::vector<PhaseShiftPartialWave>& waves,
                                 const std::string& listFile);

    /**
     * \brief Write the decays file (decays-<name>.dat) for the channel/model.
     *
     * Each (Iz >= 0, partial wave) member gets the isospin-CG decay channels into
     * its constituents; antiparticle decays are generated by Thermal-FIST.
     */
    void WritePhaseShiftDecaysFile(const PhaseShiftChannel& ch, const std::vector<PhaseShiftPartialWave>& waves,
                                   const std::string& decaysFile);

    /**
     * \brief Convenience: write both files into \p dir as list-<name>.dat and
     *        decays-<name>.dat.
     */
    void WritePhaseShiftFiles(const PhaseShiftChannel& ch, const std::vector<PhaseShiftPartialWave>& waves,
                              const std::string& dir);

    /**
     * \brief After a list build, attach a single-wave PhaseShiftDensity (from the
     *        chosen model) to every cluster member of the channel present in the
     *        system, particles and antiparticles alike.
     *
     * The model's waves are matched to members by 2J+1. Call this after LoadList()
     * (and before model calculations).
     *
     * \return The PDG ids that received a density (signed).
     */
    std::vector<long long> AttachDensities(ThermalParticleSystem& TPS,
                                           const PhaseShiftChannel& ch,
                                           const std::vector<PhaseShiftPartialWave>& waves);

    /**
     * \brief Graceful resonance subsumption: remove ch.subsumedPdg (both charge
     *        signs) from the system if present, to avoid double counting. A
     *        no-op for codes not in the list. Returns the codes actually removed.
     */
    std::vector<long long> SubsumeResonances(ThermalParticleSystem& TPS,
                                             const PhaseShiftChannel& ch);

    /// Catalog: the repulsive pi-pi I=2 channel (constituents + structure only).
    PhaseShiftChannel PiPi_I2_Channel();

    /// Catalog: the attractive pi-pi I=0 channel (the sigma/f0(500)); a single
    /// neutral isoscalar member. The sigma is kept (deg=0) in the list, not
    /// subsumed. Integrated up to the K-Kbar threshold (the elastic limit).
    PhaseShiftChannel PiPi_I0_Channel();

    /// Catalog: the pi-pi I=0 f0(980) channel - the part of delta_0^0 above the
    /// K-Kbar threshold (up to 1.42 GeV). It REUSES the real f0(980) PDG code
    /// (9010221, memberPdg): the phase-shift density is attached to that existing
    /// resonance, overriding its thermal contribution (subsumption by PDG
    /// coincidence). The f0(980) stays in the list (still a decay product).
    PhaseShiftChannel PiPi_I0_f0980_Channel();

    /// Catalog: the pi-pi I=1 P-wave channel - the rho(770). It REUSES the real
    /// rho codes (rho0=113, rho+=213, memberPdg): the density overrides the rho's
    /// contribution; the rho stays in the list (still a decay product).
    PhaseShiftChannel PiPi_I1_Channel();

    /// Catalog: the pi-pi I=1 F-wave channel - non-resonant below 1.42 GeV
    /// (rho3(1690) is higher), so a synthetic cluster (no resonance to reuse).
    PhaseShiftChannel PiPi_I1_F_Channel();

    /// Catalog: the pi-K I=3/2 (repulsive) and I=1/2 (attractive, kappa/K0*(700))
    /// channels (constituents + structure only; S = +1). The members span all Iz
    /// of the multiplet (a strange channel is not self-conjugate).
    PhaseShiftChannel PiK_I32_Channel();
    PhaseShiftChannel PiK_I12_Channel();

    /// Catalog: the pi-K I=1/2 P-wave channel - the K*(892). REUSES the real K*
    /// codes (K*(892)+=323, K*(892)0=313), overriding their contribution; elastic
    /// below the K-eta threshold. Separate from the kappa (different wave/resonance).
    PhaseShiftChannel PiK_K892_Channel();

    /// Catalog: the pi-N I=3/2 P-wave channel - the Delta(1232) (P33). A baryon
    /// channel (B=+1, fermionic), elastic across the resonance up to ~1.38 GeV.
    /// REUSES the real Delta codes (Delta++=2224, Delta+=2214, Delta0=2114,
    /// Delta-=1114), overriding their contribution; they stay as decay products.
    /// (Roy-Steiner conformal parametrization, arXiv:1510.06039.)
    PhaseShiftChannel PiN_Delta_Channel();

    /// Catalog: the non-resonant pi-N Roy-Steiner background waves (same reference,
    /// Schenk parametrization). Synthetic clusters (no resonance to reuse): each is
    /// its own single-wave channel with its own elastic cutoff. S31, P13 are
    /// elastic to ~1.38 GeV; S11, P31, P11 only below the pi-pi-N threshold ~1.22
    /// GeV. The orbital l goes in the synthetic-id excitation slot so same-(I,J)
    /// waves (S31/P31, S11/P11) get distinct codes.
    PhaseShiftChannel PiN_S31_Channel();
    PhaseShiftChannel PiN_S11_Channel();
    PhaseShiftChannel PiN_P31_Channel();
    PhaseShiftChannel PiN_P11_Channel();
    PhaseShiftChannel PiN_P13_Channel();

    /// Catalog: the K-N (kaon-nucleon, exotic S=+1) S- and P-waves (Gibbs, Arceo,
    /// Phys. Rev. C75 (2007) 054005). Meson-baryon channels (B=+1, S=+1, fermionic),
    /// non-resonant (no S=+1 resonances) -> synthetic clusters, each its own
    /// single-wave channel; elastic below the K-pi-N threshold ~1.57 GeV. Naming
    /// L_{2I,2J}: S01/P01/P03 are I=0, S11/P11/P13 are I=1.
    PhaseShiftChannel KN_S01_Channel();
    PhaseShiftChannel KN_P01_Channel();
    PhaseShiftChannel KN_P03_Channel();
    PhaseShiftChannel KN_S11_Channel();
    PhaseShiftChannel KN_P11_Channel();
    PhaseShiftChannel KN_P13_Channel();

    // ------------------------------------------------------------------------
    // High-level convenience API (no file juggling).
    // ------------------------------------------------------------------------

    /**
     * \brief One-call: add a fully-wired phase-shift channel to an already-loaded
     *        ThermalParticleSystem.
     *
     * Adds the cluster members and their antiparticles, sets the isospin-CG decay
     * channels (charge-conjugated for antiparticles), attaches the chosen model's
     * delta(M) density to each member, removes any subsumed resonances, and
     * reprocesses the decays so feeddown includes the new clusters. The channel's
     * constituent species (e.g. the pions) must already be present in the list.
     *
     * Typical use (a per-wave model set; each wave is explicit):
     * \code
     *   ThermalParticleSystem TPS("input/list/PDG2020/list.dat");
     *   std::vector<PhaseShiftPartialWave> waves;
     *   waves.push_back(PhaseShifts::AnalyticWave("pipi_I2", "GarciaMartin2011_S")); // S
     *   waves.push_back(PhaseShifts::AnalyticWave("pipi_I2", "GarciaMartin2011_D")); // D
     *   PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I2_Channel(), waves);
     *   ThermalModelIdeal model(&TPS);
     * \endcode
     *
     * \return The PDG ids of the members added (signed).
     */
    std::vector<long long> AddPhaseShiftChannel(ThermalParticleSystem& TPS,
                                                const PhaseShiftChannel& ch,
                                                const std::vector<PhaseShiftPartialWave>& waves);

    /// Look up a catalog channel by name (e.g. "pipi_I2").
    /// \throws std::invalid_argument for an unknown name.
    PhaseShiftChannel ChannelByName(const std::string& name);

    /**
     * \brief One-call: add every channel listed in a single config file.
     *
     * One line per wave, "<channel>:<wave>  <listFile>  <decayFile>  <model>".
     * Column 1 is a unique wave key (channel before the ':', wave after it); the
     * per-wave list and decay files carry the cluster members and their isospin-CG
     * decays; the model column gives delta(M):
     * \verbatim
     *   # wave       list                  decays                 model
     *   pipi_I2:S    list-pipi_I2_S.dat    decays-pipi_I2_S.dat   GarciaMartin2011_S
     *   pipi_I2:D    list-pipi_I2_D.dat    decays-pipi_I2_D.dat   GarciaMartin2011_D
     *   pipi_I2:D    list-pipi_I2_D.dat    decays-pipi_I2_D.dat   tab:delta_pipi_I2_D.dat
     *   pipi_I0_f0980:S   -    -    GarciaMartin2011_S     # reuse existing f0(980)
     * \endverbatim
     * <wave> is S/P/D/F/G or a numeric 2J+1; <model> is a wave-specific analytic
     * parametrization name or "tab:<file>". The analytic name is per-wave (two
     * waves of one paper get two names) and must agree with the wave key, else it
     * throws. File paths are resolved relative to the config file's directory. The
     * referenced .dat files are loaded into \p TPS, so they are the (hand-editable)
     * source of truth, decoupled from the in-code catalog. A "-" list/decay column
     * means "no file": the channel reuses an existing resonance (memberPdg) and the
     * model is attached to it (subsumption by PDG coincidence), nothing is added.
     *
     * \return The PDG ids of all members added across all channels.
     */
    std::vector<long long> AddPhaseShiftChannelsFromFile(ThermalParticleSystem& TPS,
                                                         const std::string& configFile);

    /**
     * \brief As AddPhaseShiftChannelsFromFile, but skips every config entry whose
     *        channel name is in \p skipChannels.
     *
     * Lets a caller (e.g. a GUI) apply a config with individual channels turned off:
     * a skipped channel is not loaded/overridden at all, so a reused resonance stays
     * the plain (pole-mass) particle from the base list. The no-skip overload above
     * delegates here with an empty set.
     */
    std::vector<long long> AddPhaseShiftChannelsFromFile(ThermalParticleSystem& TPS,
                                                         const std::string& configFile,
                                                         const std::set<std::string>& skipChannels);

    /**
     * \brief Summary of one channel referenced by a phase-shift config file: its
     *        name, the partial-wave labels present, and whether it reuses a real
     *        resonance's PDG codes (subsumption) or is a synthetic cluster.
     */
    struct PhaseShiftConfigChannel {
      std::string name;                    ///< channel name (config column 1 before ':')
      std::vector<std::string> waves;      ///< partial-wave labels (S/P/D/F/G)
      bool reusesResonance;                ///< true if it reuses real PDG codes
      std::vector<long long> reusedCodes;  ///< those codes (empty for a synthetic cluster)
      PhaseShiftConfigChannel() : reusesResonance(false) {}
    };

    /**
     * \brief Enumerate the channels referenced by a phase-shift config file - one
     *        entry per distinct channel, in file order - WITHOUT modifying any
     *        system. Useful to present the channels (e.g. a toggle list in a UI)
     *        before applying them. Commented and "-" entries are handled as in
     *        AddPhaseShiftChannelsFromFile.
     */
    std::vector<PhaseShiftConfigChannel> ListPhaseShiftConfigChannels(const std::string& configFile);

    /**
     * \brief Enable or disable every phase-shift channel already in the system,
     *        without rebuilding the list.
     *
     * Flips the enabled flag on each PhaseShiftDensity attached to a species.
     * A disabled channel contributes nothing (its density returns 0), so its
     * effect on densities, susceptibilities and feeddown vanishes while the
     * clusters stay in the list (and keep their EV/vdW parameters). Useful for a
     * cheap on/off toggle. Returns the number of channels affected.
     */
    int SetPhaseShiftsEnabled(ThermalParticleSystem& TPS, bool enabled);

    /// Number of species in the system carrying a PhaseShiftDensity.
    int CountPhaseShiftDensities(const ThermalParticleSystem& TPS);

    /// Number of phase-shift densities attached to a REUSED real resonance (a real
    /// PDG code, not a synthetic 99-prefixed cluster). When > 0 the cheap
    /// enable/disable toggle is not exact - a disabled override yields 0 instead of
    /// the resonance's own (pole-mass) contribution - so rebuild the list to toggle.
    int CountOverriddenResonances(const ThermalParticleSystem& TPS);

  } // namespace PhaseShifts

} // namespace thermalfist

#endif
