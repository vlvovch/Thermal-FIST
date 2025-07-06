# Changelog

All notable changes to this project will be documented in this file. 

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)

## [Version 1.5]

Date: 2024-03-22

Version 1.5 contains a number of new features and improvements.

## Real gas HRG

This version contains a first complete version of real gas HRG.
Real gas is a generalization of the van der Waals HRG.
It implements a generalization of the excluded volume prescription by considering non-linear EV prescriptions, such as [Carnahan-Starling EoS](http://www.sklogwiki.org/SklogWiki/index.php/Carnahan-Starling_equation_of_state).
It also allows one to incorporate arbitrary mean fields which can be a function of all densities in the system.
This allows one to model the effects of interactions in dense regime, such as in neutron stars.

A real gas model for a single component case has been formulated in [Vovchenko, Phys. Rev. C 96, 015206 (2017)](https://inspirehep.net/literature/1510282) while a paper detailing the formulation of the model for multiple components is in preparation.

Real gas HRG is implemented as a new [**ThermalModelRealGas**](include/HRGRealGas/ThermalModelRealGas.h) class.
It takes as input the [excluded volume model](include/HRGRealGas/ExcludedVolumeModelsMulti.h) and the [mean field model](include/HRGRealGas/MeanFieldModelsMulti.h), both of which can be arbitrary.
The source code contains a number of examples of both models, such as van der Waals/Carnahan-Starling excluded volume, and van der Waals or baryon density dependent mean field.

Real gas HRG fully supports the calculation of all the standard thermodynamic quantities and susceptibilities and has the complete functionality of van der Waals HRG model from before.

See examples [ThermodynamicsBQS.cpp](https://vovchenko.net/project/thermal-fist/doc/_thermodynamics_b_q_s_8cpp-example.html) and [SusceptibilitiesBQS.cpp](https://vovchenko.net/project/thermal-fist/doc/_susceptibilities_b_q_s_8cpp-example.html) which include Real Gas HRG option and [this link](https://vovchenko.net/project/thermal-fist/doc/classthermalfist_1_1_thermal_model_real_gas.html) for documentation.

## Cosmic trajectories

This version contains an implementation of cosmic trajectories in the Early Universe, given baryon, charge, lepton flavor asymmetries.

- Reference: [Vovchenko et al., Phys. Rev. Lett. 126, 012701 (2021)](https://inspirehep.net/literature/1815329).
- Class: [**CosmicEoS**](include/CosmicEos/CosmicEoS.h)
- Documentation: [this link](https://vovchenko.net/project/thermal-fist/doc/classthermalfist_1_1_cosmic_eo_s.html)
- Example: [this link](https://vovchenko.net/project/thermal-fist/doc/_cosmic_trajectory_8cpp-example.html)
- GUI: `Equation of state` tab
- Implements the EoS of HRG + leptons + photons
- Computes trajectories in the space of B, Q, and lepton flavor chemical potentials

## Other new features

- Partial chemical equilibrium with baryon annihilation

  This version contains an implementation of partial chemical equilibrium with baryon annihilation, namely $B\bar{B} \leftrightarrow n\pi$ reaction, where $n$ is a parameter~(typically $n = 5$), which can be utilized in thermal fits.
  
  - Class: [**ThermalModelPCEAnnihilation**](include/HRGPCE/ThermalModelPCEAnnihilation.h).
  - Reference: [Vovchenko, Koch, Phys. Lett. B 835, 137577 (2022)](https://inspirehep.net/literature/2172067).
  - Documentation: [this link](https://vovchenko.net/project/thermal-fist/doc/classthermalfist_1_1_thermal_model_p_c_e_annihilation.html)

- Neutron star matter calculations
  - Supported by adding leptons to the list from [/input/list/electroweak](/input/list/electroweak) and setting charge neutrality and $\mu_S = 0$
  - For a meaningful calculation, one should use an interacting version of the HRG, e.g. real gas
  - Example [here](https://vovchenko.net/project/thermal-fist/doc/_neutron_stars-_c_s_h_r_g_8cpp-example.html)

- Effective mass model for pions
  - Repulsive interactions and Bose-Einstein condensation
  - Described in https://arxiv.org/abs/2004.09004
 
- Magnetic field
  - HRG model in the magnetic field as formulated in https://arxiv.org/abs/1301.1307
  - Additional features: Hadron number densities and fluctuations
  - GUI: `PCE/Saha/Other...` configuration
  - Reference: [Vovchenko, Phys. Rev. C 110, 034914 (2024)](https://inspirehep.net/literature/2790376)

- New thermodynamic quantities
  - Temperature derivatives of the energy, entropy and hadron number densities
  - Temperature derivatives of susceptibilities
  - Hessian matrix of the pressure wrt $T$ and $\mu$'s
  - Adiabatic and isothermal sound velocity (with multiple conserved charges and specific entropies)
  - Specific heat capacity $c_V/T^3$

## Event generator
- Shear and bulk viscous corrections
  - Part of FIST sampler (see also https://github.com/vlvovch/fist-sampler)
  - Currently no excluded volume and bulk viscosity simultaneously
  - Thanks to: [J. Hammelmann](https://github.com/JanHammelmann) and [N. Goetz](https://github.com/ngoetz)
- Support for output to SMASH
- Propagation of particles after (weak) decays and DCA cuts

## New examples
- [Cosmic trajectories](https://vovchenko.net/project/thermal-fist/doc/_cosmic_trajectory_8cpp-example.html)
- [Neutron star matter EoS at T = 0](https://vovchenko.net/project/thermal-fist/doc/_neutron_stars-_c_s_h_r_g_8cpp-example.html)
- [Thermodynamic functions in $(T, \mu_B, \mu_Q, \mu_S)$ phase diagram](https://vovchenko.net/project/thermal-fist/doc/_thermodynamics_b_q_s_8cpp-example.html)
- [Susceptibilities, speed of sound, and heat capacity in $(T, \mu_B, \mu_Q, \mu_S)$ phase diagram](https://vovchenko.net/project/thermal-fist/doc/_susceptibilities_b_q_s_8cpp-example.html)

## Technical improvements

- New logo
- GUI supports and uses Qt6 by default now
  - Qt5 still works but not will not be guaranteed to work in the future
- Linear algebra: Migrated to Eigen-3.4.0
- GitHub actions for automated tests/compilation on various platforms
- Expanded unit tests
- Better warning and error handling
- Various fixes of minor errors

## GUI enhancements
- Calculations within real gas HRG (in the model selection)
- Baryon annihilation in PCE (`PCE/Saha/Other...`)
- Effective mass model for pions (`PCE/Saha/Other...`)
- Tabular view of calculation in the EoS tab and export to file (`Equation of state` tab)
- Calculation of cosmic trajectories (`Equation of state` tab)
- Calculations at fixed conserved charge densities (instead of fixed $\mu$'s) (`Use densities` checkbox in `Thermal model` tab)
- Various enhancements


## [Version 1.4.2] 

Date: 2022-11-07

- Fix bug (introduced in v1.4.1) with the normalization of two-particle susceptibilities in ThermalModelIdeal

## [Version 1.4.1] 

Date: 2022-10-18

- GUI: Added more clarity concerning the loading of the decays list
- GUI: Number of known decays channels shown in the table for unstable hadrons
- Added explicit calculation of Fermi-Dirac integrals at T = 0

## [Version 1.4] 

Date: 2022-08-25

Main new feature: Excluded volume effect in the Monte Carlo event generator

## Enhancements
- Hard-core repulsion in event generator (see [FIST sampler](https://github.com/vlvovch/fist-sampler))
- Crossterms EV-HRG model is now implemented as partial case of QvdW-HRG model (for a = 0)
- GUI: Blast-wave model: Option to set the value of parameter r_max (the maximum transverse radius)
- GUI: Minor improvements for PCE-HRG model thermal fits

## Bugfixes
- GUI: Fix bug with the sign of strongly intensive quantity Delta
- GUI: Options for different treatments of hard-core repulsion in event generator

## [Version 1.3.4] 

Date: 2022-05-05

GUI: Fixed a bug with the calculation of scaled moments and Pearson correlation coefficients that was introduced in the previous version 1.3.3

## [Version 1.3.3] 

Date: 2022-05-03

GUI: Added scaled moments and Pearson correlation coefficients into the correlations dialog

## [Version 1.3.2] 

Date: 2022-02-18

Minor enhancements and bugfixes

## Enhancements
- GUI: Binning options in event generator
- GUI: Ability to add resonances and other particles to online analysis
- GUI: Export Monte Carlo momentum distributions to file
- Event generator: Monte Carlo sampler of ideal HRG from an arbitrary (numerical) hypersurface
- Event generator: Event generator output through a new event writer class
- Event generator: Experimental support for HepMC format
- Lists: Checks the list for particles with non-zero strangeness/charm but zero strangeness/charm quark content and outputs warnings

## Bugfixes
- Fix the ordering of p0 and px,py,pz output in ascii event generator output. Was px,py,pz,p0, now it is p0,px,py,pz which is consistent with the header line
- Strangeness/charm neutrality equations will now be solved properly for lists with incorrect zero |s| and |c| columns

## [Version 1.3.1] 

Date: 2021-05-26

Minor enhancements and bugfixes

## Enhancements
- GUI: Output moments in addition to susceptiblities in the correlations dialog, output nudyn and D measure in equation of state tab
- GUI: Show lattice data for ratios of quantities for the equation of state tab in the GUI, where available
- GUI: Display mean pT in the event generator tab
- Event generator: Sampling of the space-time coordinates in addition to the momenta (*not yet supported for the Siemens-Rasmussen model*)
- Event generator: Write the output in the format tailored for UrQMD afterburner from https://github.com/jbernhard/urqmd-afterburner
- Event generator: Output particle energy (p0) by default when writing to file


## Bugfixes
- Fix the flag setting in ThermalModelPCE for the Saha equation mode and the freeze-out of long-lived resonances
- GUI: Fix the output of strongly intensive quantities for net-particle correlations
- Event generator: Correct resonance mass-energy rescaling when masses of decay products exceed the mass of a resonance
- Add contributions of the K0S decays to the inclusive (weak decay included) yields of pions. Note that decay contributions from K0L are not included

## [Version 1.3] 

Date: 2020-11-07

Version 1.3 contains a number of new features, improvements, and bugfixes

## Partial chemical equilibrium (PCE)

This release implements HRG in partial chemical equilibrium (PCE) [Bebie et al, [Nucl. Phys. B **378**, 95 (1992)](https://inspirehep.net/literature/316591)], as appropriate for modeling the hadronic phase of heavy-ion collisions. PCE is implemented as a new **ThermalModelPCE** class, and has been used in two recent papers, [Phys. Lett. B **800**, 135131 (2020)](https://inspirehep.net/literature/1726464) and [Phys. Rev. C **102**, 024909 (2020)](https://inspirehep.net/literature/1751960)

See [PCE-Saha-LHC.cpp](./src/examples/PCE/PCE-Saha-LHC.cpp) for a usage example and [this link](https://vovchenko.net/project/thermal-fist/doc/classthermalfist_1_1_thermal_model_p_c_e.html) for documentation.

The PCE-HRG model is available in the thermal fitting routine, with Tkin being an extra fitting parameter

The PCE-HRG is also accessible for use in the graphical user interface program through a button `PCE/Saha/Other...`

## Particle lists

- The default particle list is now based on the 2020 edition of the PDG listing. As before, the list contains only hadrons and resonances with an established status. Changes compared to the old 2014 PDG list are mild.

- Version 1.3 now fully supports charm. The charmed hadrons are available in an extended input list file [list-withcharm.dat](./input/list/PDG2020/list-withcharm.dat)

- The excited nuclei with mass number up to A = 5 have been added with this release, they are available in an extended input list file [list-withexcitednuclei.dat](./input/list/PDG2020/list-withexcitednuclei.dat). The reference for this list of excited nuclei is [Phys. Lett. B **809**, 135746 (2020)](https://inspirehep.net/literature/1790683)

- The hadron list converted from [SMASH](https://smash-transport.github.io/) transport model (version 1.8) is available [here](./input/list/SMASH-1.8)

## Other changes

- Mixed-canonical ensembles in Monte Carlo event generator
- Reworked cylindrical blast-wave event generator with much faster initialization
- Cracow freeze-out generator
- Fermi-Dirac and Bose-Einstein momentum samplers for the event generator with quantum statistics
- Correct decay kinematics involving leptons and photons, the decay leptons/photons, if any, are included in event generator output
- Option for modular loading of particle lists
- Functions for computing net-particle/net-charge correlators (see documentation)
- Leptons and photons will appear in decays of some hadrons with correct kinematics, even if they are not included in the particle list

## Changes to the graphical user interface
- Calculation of matrices of various 2nd-order correlators (accessible on the `Thermal model` tab via `Correlations...` button
- For quantum van der Waals and X-terms excluded volume models: the EV/vdW interactions by default are now only for baryon-baryon and antibaryon-antibaryon pairs, in order to make it consistent with the literature
- Enable PCE mode via `PCE/Saha/Other...`, then possible to fit Tkin in the thermal fits mode
- Option for quantum statistical momentum distribution in the event generator mode (note that multiplicity distributions are still Poisson!)
- Mixed-canonical ensembles in the event generator mode
- Cracow model in the event generator mode

## Bugfixes
- Fix for Bose-Einstein condensation issue in event generator with finite resonance widths
- Fix for Fermi-Dirac integrals in case of massless particles
- Fix for decay properties of some charmed hadrons
- Fix for momentum sampling in event generators in case of zero collective/radial flow

## [Version 1.2] 

Date: 2019-07-22

**The version of the code published in [Computer Physics Communications](https://doi.org/10.1016/j.cpc.2019.06.024)**

### Changes
- Improved performace of the event generator. Exact baryon number conservation now taken into account using [Devroye's method](https://doi.org/10.1016/S0167-7152(02)00055-X) of sampling the Bessel distribution
- Possibility to set volume on an event-by-event basis in the event generator
- The beta parameter in cylindrical blast-wave model now corresponds to the surface velocity parameter \beta_s
- Mixed-canonical ensembles are now available through GUI (use the canonical ensemble and conservation laws dialogs)

### Bugfixes
- Mixed-canonical ensembles now work at non-zero densities of conserved charges
- Correct analytic calculation of fluctuations in mixed-canonical ensembles
- Small violation of electric charge conservation in event generator within full canonical ensemble
- Proxy net charge susceptibility now corresponds again to net-charge instead of net-pion number. This fixes the output of the cpc4 example, a bug introduced in the pre-release of version 1.2

## [Version 1.1] 

Date: 2019-01-25

- A pdg code is now held in a 64-bit integer (was 32-bit)
- A pdg code with a trailing zero in thermal fit/yields applications can now be interpreted as a sum of yields of particles and antiparticles, which correspond to the pdg code without this trailing zero
- The fit procedure is improved to catch possible issues with a Bose-Einstein condensation
- cpc2 and cpc3 examples slightly reworked and optimized, calculation of fit number of dof in cpc3 corrected
- Cosmetic changes to the code, mainly in response to compiler warnings

## [Version 1.0] 

Date: 2019-01-16

Matches the code released with the [arXiv documentation paper](https://arxiv.org/abs/1901.05249v1)

- Charm particles
- Possibility to set seed in the random generator

## [Version 0.9] 

Date: 2019-01-09

**Library**

- Possibility to constrain baryochemical potential by entropy-per-baryon ratio
- Support for charm-canonical ensemble in event generator
- Code structure optimization
- Doxygen documentation

**GUI frontend**

- More convenient HRG model and thermal fit configuration in GUI
- Thermal fit plots
- New equation of state tab for calculating temperature dependencies of various observables,
- Possibility to increase/decrease font size

**Other**

- More files with various experimental data
- New examples

## [Version 0.8] 

Date: 2018-12-26

#### Library
- New, more human-readable format of decays input file (old one still supported)
- Broyden method related routines refactored. Non-linear equations are now solved when needed using the generic routines, resulting in a cleaner code
- Fix for calculations of the ideal gas functions in the case of zero particle mass
- Decay feeddown contributions can now be separated into strong/electromagnetic/weak decays
- Fixed bug regarding light quark content calculations for hypernuclei (only relevant if both chemical non-eq. of light quarks  and light hypernuclei are considered)
- All third-party libraries/code moved into the `thirdparty` folder
- Resonance widths: eBW scheme with constant branching ratios added (very similar results to the standard eBW scheme but much faster)

#### GUI
- New tab for editing the particle list on the fly
- Tool-tips
- Possibility to do calculations in the Thermal Model tab using the parameters extracted from a thermal fit in the other tab
- *About Thermal-FIST* dialog window
- Menu, with an option to increase/decrease font size
- Splash screen, app icon (Windows only)

## [Version 0.7] 

Date: 2018-12-10

- eBW scheme now also implemented in event generator
- Quantum statistics by default now computed using quadratures (was cluster expansion)
- Zero-width scheme is now the default one (was energy-independent Breit-Wigner)
- Spectral functions and (energy-dependent) branching ratios visualization in gui
- Some improvements in numerics at low temperatures

## [Version 0.6] 

Date: 2018-08-02

**The first public version of Thermal-FIST**

[Version 1.5]: https://github.com/vlvovch/Thermal-FIST/compare/v1.4.2...v1.5

[Version 1.4.2]: https://github.com/vlvovch/Thermal-FIST/compare/v1.4.1...v1.4.2

[Version 1.4.1]: https://github.com/vlvovch/Thermal-FIST/compare/v1.4...v1.4.1

[Version 1.4]: https://github.com/vlvovch/Thermal-FIST/compare/v1.3.4...v1.4

[Version 1.3.4]: https://github.com/vlvovch/Thermal-FIST/compare/v1.3.3...v1.3.4

[Version 1.3.3]: https://github.com/vlvovch/Thermal-FIST/compare/v1.3.2...v1.3.3

[Version 1.3.2]: https://github.com/vlvovch/Thermal-FIST/compare/v1.3.1...v1.3.2

[Version 1.3.1]: https://github.com/vlvovch/Thermal-FIST/compare/v1.3...v1.3.1

[Version 1.3]: https://github.com/vlvovch/Thermal-FIST/compare/v1.2(cpc)...v1.3

[Version 1.2]: https://github.com/vlvovch/Thermal-FIST/compare/v1.1...v1.2(cpc)

[Version 1.1]: https://github.com/vlvovch/Thermal-FIST/compare/v1.0...v1.1

[Version 1.0]: https://github.com/vlvovch/Thermal-FIST/compare/v0.9...v1.0

[Version 0.9]: https://github.com/vlvovch/Thermal-FIST/compare/v0.8...v0.9

[Version 0.8]: https://github.com/vlvovch/Thermal-FIST/compare/v0.7...v0.8

[Version 0.7]: https://github.com/vlvovch/Thermal-FIST/compare/v0.6...v0.7

[Version 0.6]: https://github.com/vlvovch/Thermal-FIST/releases/tag/v0.6

