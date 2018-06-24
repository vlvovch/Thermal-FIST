# Thermal-FIST
[![Build Status](https://travis-ci.org/vlvovch/Thermal-FIST.svg?branch=master)](https://travis-ci.org/vlvovch/Thermal-FIST)

**Thermal-FIST** (or simply **The FIST**) is a C++ package designed for the convenient general-purpose analysis within the family of the hadron resonance gas (HRG) models.
This mainly includes the statistical analysis of particle production in heavy-ion collisions and the phenomenology of the hadronic equation of state. 

Particular emphasis is put on fluctuations and correlations of conserved charges, effects of probabilistic decay, chemical non-equilibrium, and inclusion of hadronic interactions.

Calculations are possible within the grand canonical ensemble, the canonical ensemble, as well as in mixed-canonical ensembles combining canonical treatment of strangeness/charm with the grand-canonical treatment of other conserved numbers.

The package contains a fast thermal event generator, with Blast Wave model based momentum distributions, and possibility of simultaneous inclusion of effects of exact charge conservation and hadronic interactions.

The package also includes **QtThermalFIST** -- a Qt-based graphical user interface designed for a fast and convenient general-purpose thermal model analysis. 

### Notes on the present version
The package is not fully documented yet, but should be perfectly usable. See the [Quick Start Guide](docs/quickstart.md)

Note that calculations may break down for excessive/overstressed parameters.
This should be carefully monitored.

*Copyright (C) 2018  Volodymyr Vovchenko*
