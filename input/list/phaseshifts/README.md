# S-matrix / phase-shift list modules

Everything here is **per partial wave**. Each wave of a scattering channel is a
drop-in Thermal-FIST list module — a particle table (`list-<channel>_<wave>.dat`)
plus a decays file (`decays-<channel>_<wave>.dat`) — that adds the wave's isospin
multiplet as effective "cluster" degrees of freedom. The clusters carry the
channel's conserved charges and decay (with isospin Clebsch-Gordan branchings)
into their constituent hadrons.

The synthetic PDG ids follow `99 F (2I) (2|Iz|) nn (2J+1)`: the last digit is
the wave's spin degeneracy `2J+1`, the `99` prefix keeps them clear of real PDG
codes, and the sign selects particle/antiparticle (generated automatically at
load). Charges/baryon/strangeness live in the table columns, not in the id.

These files only define the *structure*. The actual phase shift `delta(M)` — the
dynamics — is a per-wave model: a wave-specific analytic parametrization name
(`PhaseShifts::AnalyticWave(channel, param)`, e.g. `param = "GarciaMartin2011_S"`)
or a tabulated phase shift (`PhaseShifts::TabulatedWave(2J+1, file)`). The model
name is per-wave (two waves of one paper get two names — they have different phase
shifts); there is no "all waves" model.

## Usage

Single config file (recommended) — one line **per wave**, referencing that wave's
list/decay files and its model. `AddPhaseShiftChannelsFromFile` loads the
referenced `.dat` files (members + antiparticles + decays), attaches each wave's
`delta(M)`, and removes any subsumed resonances. A ready-made `pipi.conf` ships
here:

```cpp
#include "HRGPhaseShifts/PhaseShiftChannel.h"
using namespace thermalfist;

ThermalParticleSystem TPS("input/list/PDG2020/list.dat");   // your usual base list
PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, "input/list/phaseshifts/pipi.conf");
ThermalModelIdeal model(&TPS);
```
```
# pipi.conf -- one line per wave: <channel>:<wave> <list> <decays> <model>
pipi_I2:S   list-pipi_I2_S.dat   decays-pipi_I2_S.dat   GarciaMartin2011_S   # S (J=0)
pipi_I2:D   list-pipi_I2_D.dat   decays-pipi_I2_D.dat   GarciaMartin2011_D   # D (J=2)
# swap one wave for a tabulated phase shift (table is M[GeV] delta[rad]):
# pipi_I2:D list-pipi_I2_D.dat   decays-pipi_I2_D.dat   tab:delta_pipi_I2_D.dat
```

Column 1 is a unique wave key: the channel before the `:`, the wave (S/P/D/F/G or
a numeric 2J+1) after it. Columns 2-3 are the per-wave list/decay files (paths
resolved relative to the config file). Column 4 is the model: an analytic
`<param>` name, or `tab:<file>` for a tabulated wave. The `.dat` files are the
hand-editable source of truth — edit them (or point at your own) and the loader
picks them up. Cluster names show the wave letter, e.g. `pipi_I2[Iz=+1,D]`.

Programmatic (no files) — build a wave set in code and add it in one call:

```cpp
ThermalParticleSystem TPS("input/list/PDG2020/list.dat");
std::vector<PhaseShiftPartialWave> waves;
waves.push_back(PhaseShifts::AnalyticWave("pipi_I2", "GarciaMartin2011_S"));  // S
waves.push_back(PhaseShifts::AnalyticWave("pipi_I2", "GarciaMartin2011_D"));  // D
PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I2_Channel(), waves);
```

Low-level — load the per-wave `.dat` files through the standard list build, then
attach each wave's model:

```cpp
ThermalParticleSystem TPS(
  { "input/list/PDG2020/list.dat",   "input/list/phaseshifts/list-pipi_I2_S.dat",
                                      "input/list/phaseshifts/list-pipi_I2_D.dat" },
  { "input/list/PDG2020/decays.dat", "input/list/phaseshifts/decays-pipi_I2_S.dat",
                                      "input/list/phaseshifts/decays-pipi_I2_D.dat" });
PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I2_Channel();
PhaseShifts::AttachDensities(TPS, ch, { PhaseShifts::AnalyticWave("pipi_I2", "GarciaMartin2011_S") });
PhaseShifts::AttachDensities(TPS, ch, { PhaseShifts::AnalyticWave("pipi_I2", "GarciaMartin2011_D") });
PhaseShifts::SubsumeResonances(TPS, ch);
```

The `.dat` files here are generated from the channel definitions via
`PhaseShifts::WritePhaseShiftFiles(channel, waves, dir)` (one list + one decay
file per wave) and committed for inspection; regenerate them if a channel's
structure changes.

## Excluded volume / van der Waals

The cluster members are ordinary hadrons (a pi-pi cluster has B = 0, so it is a
meson), so they get excluded-volume / vdW parameters from the same
quantum-number-based assignment as every other species — no special-casing.

## Channels

- `pipi_I2` — repulsive pi-pi I=2 (S- and D-waves; Garcia-Martin et al.,
  Phys. Rev. D83 (2011) 074004). Non-resonant: subsumes no resonances.
  Config: `pipi.conf`.
- `pipi_I0` — attractive pi-pi I=0 (S-wave; same reference, Appendix A), i.e. the
  sigma/f0(500). The phase passes through 90 deg, so it is branch-tracked
  (`atan2`). Isoscalar and neutral, so it contributes to the EoS and pion feeddown
  but not to charge fluctuations. Integrated up to the K-Kbar threshold (2 M_K),
  the elastic limit. It **reuses the real f0(500) code (9000221)**: the density
  overrides the sigma's (deg=0) list entry. Config: `pipi.conf`.
- `pipi_I0_f0980` — the part of delta_0^0 above the K-Kbar threshold (the
  f0(980)). It **reuses the real f0(980) code (9010221)**, overriding its thermal
  contribution. **OFF by default** (commented out in `pipi.conf`): that region is
  strongly inelastic (eta_0^0 ~ 0.35, the f0(980) is mostly K-Kbar), which the
  phase-shift-only Beth-Uhlenbeck ignores - a knowingly-crude ~2% effect. Doing it
  right needs the coupled-channel (pi-pi/K-Kbar) treatment. Uncomment to enable.
- `pipi_I1` — the rho(770): pi-pi I=1 P-wave (same reference, Eq. A7),
  branch-tracked through 90 deg at the rho mass; elastic up to the K-Kbar
  threshold. It **reuses the real rho codes (rho0=113, rho+=213)**, overriding
  the rho's contribution (the P-wave 2J+1=3 carries the rho spin). Config:
  `pipi.conf`.
- `pipi_I1_F` — pi-pi I=1 F-wave (same reference, Eq. A14). Non-resonant below
  1.42 GeV (rho3(1690) is higher) and small (k^7-suppressed); elastic to 1.42 GeV.
  A synthetic cluster (no resonance to reuse). Config: `pipi.conf`.
- `piK_I32` — repulsive pi-K I=3/2 (S, P, D waves; Pelaez, Rodas, Phys. Rev. D93
  (2016) 074025). Non-resonant, elastic up to ~1.74 GeV; the P and D waves are
  tiny. Synthetic clusters. Config: `piK.conf`.
- `piK_I12` — attractive pi-K I=1/2 S-wave (same reference), i.e. the
  kappa/K0*(700). Elastic only below the K-eta threshold. It **reuses the real
  kappa codes (9000321, 9000311)**; since the kappa is usually excluded from the
  list, these are created with the isospin-CG decays (if present, overridden).
  Config: `piK.conf`.
- `piK_K892` — pi-K I=1/2 P-wave = the K*(892) (same reference, Eq. 28). Resonant
  (branch-tracked through 90 deg), elastic below the K-eta threshold. It **reuses
  the real K*(892) codes (323, 313)**, overriding their contribution. Config:
  `piK.conf`.
- `piN_Delta` — pi-N I=3/2 P-wave = the **Delta(1232)** (P33; Hoferichter et al.,
  Phys. Rept. 625 (2016) 1, Roy-Steiner conformal parametrization). A
  **meson-baryon** channel (B=+1, **fermionic**), resonant (branch-tracked through
  90 deg), elastic across the resonance up to ~1.38 GeV. It **reuses the real Delta
  codes (2224, 2214, 2114, 1114)**, overriding their contribution. Config:
  `piN.conf`.
- `piN_S31`, `piN_S11`, `piN_P31`, `piN_P11`, `piN_P13` — the non-resonant pi-N
  Roy-Steiner background waves (same reference, Schenk parametrization). **Synthetic
  clusters** (no resonance to reuse), each its own single-wave channel with its own
  elastic cutoff: S31 (repulsive) and P13 are elastic to ~1.38 GeV; S11
  (attractive), P31 and P11 only below the pi-pi-N threshold ~1.22 GeV. Config:
  `piN.conf`.

  **Wave keys for baryons:** the config wave key is the numeric `2J+1` (S/P J=1/2 ->
  `2`, the J=3/2 waves P13/P33 -> `4`), because the S/P/D letters encode `2J+1 =
  2l+1`, valid only for mesons (J=l). For baryons J = l +- 1/2, so distinct waves
  can share `2J+1` (S31/P31 and S11/P11 all have `2J+1=2`); the orbital `l` is
  carried in the synthetic-id excitation (`n n`) slot to keep their codes distinct.
  All Roy-Steiner coefficients and a numerical validation are in
  `piN_RoySteiner_reference.md`.
- `KN_S01`, `KN_P01`, `KN_P03`, `KN_S11`, `KN_P11`, `KN_P13` — **kaon-nucleon**
  (exotic S=+1) S- and P-waves (Gibbs, Arceo, Phys. Rev. C75 (2007) 054005). A
  meson-baryon channel (B=+1, S=+1, **fermionic**); **non-resonant** (no S=+1
  resonances), so synthetic clusters, each its own single-wave channel; elastic
  below the K-pi-N threshold ~1.57 GeV. The I=1 S-wave is strongly repulsive; the
  I=0 P-waves show a spin-orbit splitting. Same baryon wave-key convention as pi-N
  (numeric 2J+1, orbital l in the excitation slot). Coefficients + validation in
  `KN_GibbsArceo_reference.md`. Config: `KN.conf`.

### Subsumption by PDG coincidence
A channel can REUSE a real resonance's PDG code (`memberPdg`, the `-` list/decay
columns in the config) instead of a synthetic id. If that code is already in the
list (sigma, f0(980)) the phase-shift density **overrides** its thermal
contribution while it stays in the list (still a decay product); if absent (the
kappa) it is **created** with the channel's decays. Either way the real resonance
is not separately counted. Because these are real (non-synthetic)
codes, the cheap enable/disable toggle is not exact for them, so the GUI rebuilds
the list when toggling (see `CountOverriddenResonances`).

pi-K (S=+1), pi-N (B=+1) and K-N (B=+1, S=+1) carry a conserved charge, so unlike
pi-pi they are **not** self-conjugate multiplets: every Iz is a distinct member and
the antiparticles form the conjugate (S=-1 / antibaryon) sector. The builder handles
this automatically; synthetic ids encode the multiplet index I+Iz (0..2I) rather
than 2|Iz|. To run several channels together, load their configs together (the GUI
loads `pipi.conf` + `piK.conf` + `piN.conf` + `KN.conf` by default).
