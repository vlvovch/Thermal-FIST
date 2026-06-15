# S-matrix / phase-shift list modules

Each scattering channel is a drop-in Thermal-FIST list module — a particle
table (`list-<name>.dat`) plus a decays file (`decays-<name>.dat`) — that adds
the channel's isospin multiplet as effective "cluster" degrees of freedom. The
clusters carry the channel's conserved charges and decay (with isospin
Clebsch-Gordan branchings) into their constituent hadrons.

The synthetic PDG ids follow `99 F (2I) (2|Iz|) nn (2J+1)`: the last digit is
the wave's spin degeneracy `2J+1`, the `99` prefix keeps them clear of real PDG
codes, and the sign selects particle/antiparticle (generated automatically at
load). Charges/baryon/strangeness live in the table columns, not in the id.

These files only define the *structure*. The actual phase shift `delta(M)` — the
dynamics — is attached in code by `PhaseShifts::AttachDensities(TPS, channel,
model)`, where `model` is a selectable `PhaseShiftModel` (a named analytic
parametrization or a tabulated phase shift). This separation lets the same
channel be evaluated with different models.

## Usage

Simplest — one call adds a fully-wired channel (members + antiparticles + decays
+ delta(M) + resonance subsumption) to an already-loaded list:

```cpp
#include "HRGPhaseShifts/PhaseShiftChannel.h"
#include "HRGPhaseShifts/PhaseShiftModel.h"
using namespace thermalfist;

ThermalParticleSystem TPS("input/list/PDG2020/list.dat");   // your usual base list
PhaseShifts::AddPhaseShiftChannel(TPS, PhaseShifts::PiPi_I2_Channel(),
    PhaseShifts::AnalyticModel("pipi_I2:GarciaMartin2011"));
ThermalModelIdeal model(&TPS);
```

Single config file listing the channels you want (one `<channel> <model>` per
line) — wire them all in one call. A ready-made `pipi.conf` ships here:

```cpp
ThermalParticleSystem TPS("input/list/PDG2020/list.dat");
PhaseShifts::AddPhaseShiftChannelsFromFile(TPS, "input/list/phaseshifts/pipi.conf");
```
```
# pipi.conf
pipi_I2   pipi_I2:GarciaMartin2011
# tabulated alternative for the same channel:
# pipi_I2 tab:1=delta_S.dat,5=delta_D.dat
```

Low-level (data-file modules) — load the generated `.dat` files through the
standard list build, then attach the model:

```cpp
ThermalParticleSystem TPS(
  { "input/list/PDG2020/list.dat",  "input/list/phaseshifts/list-pipi_I2.dat" },
  { "input/list/PDG2020/decays.dat","input/list/phaseshifts/decays-pipi_I2.dat" });
PhaseShifts::PhaseShiftChannel ch = PhaseShifts::PiPi_I2_Channel();
PhaseShifts::PhaseShiftModel model = PhaseShifts::AnalyticModel("pipi_I2:GarciaMartin2011");
PhaseShifts::AttachDensities(TPS, ch, model);
PhaseShifts::SubsumeResonances(TPS, ch);
```

The `.dat` files here are generated from the channel definitions via
`PhaseShifts::WritePhaseShiftFiles(channel, model, dir)` and committed for
inspection; regenerate them if a channel's structure changes.

## Channels

- `pipi_I2` — repulsive pi-pi I=2 (S- and D-waves; Garcia-Martin et al.,
  Phys. Rev. D83 (2011) 074004). Non-resonant: subsumes no resonances.
