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

The cluster members are added as ordinary hadrons (a pi-pi channel cluster has
B = 0, so it is a meson). They therefore receive excluded-volume / vdW
parameters from the *same* quantum-number-based assignment as every other
species (uniform radius, mass-proportional `(m/m_N)^{1/3}`, or
baryon-content / point-like-meson). No special-casing is applied or needed:
because the clusters are suppressed consistently with the pions they are made
of, total (free + feeddown) particle densities stay consistent — whereas
leaving the clusters EV-free while the pions are suppressed could drive the
total pion density negative.

Notes:
- The cluster's nominal mass (the threshold m1 + m2) is what a mass-proportional
  radius scheme sees; its thermodynamics still come from the phase-shift
  spectral integral, not this nominal mass.
- In the diagonal-EV model clusters never source excluded volume (their density
  is an interaction correction, possibly negative) but are suppressed by the
  global free-volume factor; in the matrix models (Crossterms/VDW/RealGas) they
  participate through their (meson) parameters like any other meson.
- In the GUI, enabling phase shifts adds the clusters once; toggling the
  checkbox afterwards just switches their densities on/off without rebuilding the
  list (PhaseShifts::SetPhaseShiftsEnabled). The clusters keep their meson
  EV/vdW parameters either way; when off their zero density makes them inert, so
  interactions never need manual re-apply.

## Channels

- `pipi_I2` — repulsive pi-pi I=2 (S- and D-waves; Garcia-Martin et al.,
  Phys. Rev. D83 (2011) 074004). Non-resonant: subsumes no resonances.
