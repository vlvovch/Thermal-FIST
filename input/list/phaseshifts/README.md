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
  contribution. Config: `pipi.conf`.
- `piK_I32` — repulsive pi-K I=3/2 (S-wave; Pelaez, Rodas, Phys. Rev. D93 (2016)
  074025). Non-resonant, elastic up to ~1.74 GeV. Synthetic clusters. Config:
  `piK.conf`.
- `piK_I12` — attractive pi-K I=1/2 (S-wave; same reference), i.e. the
  kappa/K0*(700). Elastic only below the K-eta threshold. It **reuses the real
  kappa codes (9000321, 9000311)**; since the kappa is usually excluded from the
  list, these are created with the isospin-CG decays (if present, overridden).
  Config: `piK.conf`.

### Subsumption by PDG coincidence
A channel can REUSE a real resonance's PDG code (`memberPdg`, the `-` list/decay
columns in the config) instead of a synthetic id. If that code is already in the
list (sigma, f0(980)) the phase-shift density **overrides** its thermal
contribution while it stays in the list (still a decay product); if absent (the
kappa) it is **created** with the channel's decays. Either way the real resonance
is not separately counted. Because these are real (non-synthetic)
codes, the cheap enable/disable toggle is not exact for them, so the GUI rebuilds
the list when toggling (see `CountOverriddenResonances`).

pi-K carries strangeness (S=+1), so unlike pi-pi it is **not** a self-conjugate
multiplet: every Iz is a distinct member (Q = Iz + 1/2 for the S=+1 sector) and
the antiparticles form the S=-1 sector. The builder handles this automatically;
synthetic ids encode the multiplet index I+Iz (0..2I) rather than 2|Iz|. To run
pi-pi + pi-K together, load both `pipi.conf` and `piK.conf`.
