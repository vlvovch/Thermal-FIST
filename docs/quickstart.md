# Thermal-FIST: Quick start guide

### External dependencies

The external dependencies are minimal.
The core library has no external dependencies.
The thermal fitting routines use the MINUIT2 library from CERN ROOT which is included as a standalone package. If a ROOT installation with MINUIT2 
is found in the system, MINUIT2 from the installation is used instead.

The QtThermalFIST GUI frontend requires the open source [**Qt5 framework**](http://qt-project.org) to be installed. Otherwise the GUI frontend is not built.

### Building

The preferred way is to use **cmake**.

For example, to download and build the package on a Linux-like system 
run the following commands in a command shell:
~~~.bash
git clone https://github.com/vlvovch/Thermal-FIST.git
cd Thermal-FIST
mkdir build
cd build
cmake ../
make
~~~

This will build the libraries in `build/lib`, the QtThermalFIST GUI in `build/bin`,
and a couple of sample macros in `build/bin/examples`.

Please see also the commands listed in the [**.travis.yml**](../.travis.yml) file for another example of building and running the package under a Linux system.

### QtThermalFIST

The standard analysis can be performed within the GUI.

From the build folder, run with
```bash
./bin/QtThermalFIST
```

### Running in the cloud

You can use [GitHub Codespaces](https://github.com/features/codespaces) to run Thermal-FIST in the cloud from a browser. To run the graphical user interface program QtThermalFIST from a browser, follow the guide [here](https://github.com/devcontainers/features/tree/main/src/desktop-lite#connecting-to-the-desktop).

### Using QtThermalFIST

Using the GUI should be self-explanatory for most part. 

The particle list is read from a file.
The default list is located in the `$(ThermalFIST)/input/list/PDG2014` folder, containing particles consisting of light and strange quarks in accordance with the PDG2014 compilation.
There the `list-withnuclei.dat` file contains hadrons plus light nuclei, while the `list.dat` file contains hadrons only. There are also two files with a `-withcharm` prefix which contain charm hadrons
in addition to light and strange hadrons.
In the same folder there is also the `decays.dat` file where all the decay channels are listed.
Additionally, the `$(ThermalFIST)/input/list/` folder also contains particle lists extracted from the open source THERMUS-2.3/3.0 ([link](http://www.phy.uct.ac.za/phy/people/academic/wheaton/research)) package, which can be used to compare the results between different lists and also to compare Thermal-FIST and THERMUS.

The **QtThermalFIST** GUI program consists of five tabs:

#### *Thermal model* tab

Performs a single thermal model calculation for a given set of thermal parameters.
The model configuration is specified in the right-hand side panel.

**HRG model.**

1. **Ideal** - the standard ideal gas of hadrons and resonances
2. **Excluded-volume (Diagonal)** - the most widely used excluded-volume HRG model. References: [nucl-th/9711062](https://arxiv.org/abs/nucl-th/9711062) and [nucl-th/9808012](https://arxiv.org/abs/nucl-th/9808012) 
3. **Excluded-volume (X-terms)** - generalized ''Crossterms'' EV model. EV parameters specified for each pair of species, instead of for each specie only. In this way one can model differently, e.g., baryon-baryon and baryon-antibaryon repulsion. See [nucl-th/9906068](https://arxiv.org/abs/nucl-th/9906068) and [1606.06218](https://arxiv.org/abs/1606.06218) for details.
4. **QvdW-HRG** - multi-component quantum van der Waals HRG. Both attractive and repulsive interaction parameters can be specified for any pair of species. See [1707.09215](https://arxiv.org/abs/1707.09215)

**Ensemble.** The statistical ensemble (canonical vs grand-canonical)

1. **Grand-canonical** - the most commonly used one. Conserved charges are conserved "on average".
2. **Canonical** - conserved charges are conserved exactly. Only Ideal HRG supported.
3. **Strangeness-canonical** - Canonical treatment of strangeness combined with the grand canonical treatment of baryon number and electric charge. Charm not supported. Only Boltzmann statistics for strange particles.
4. **Charm-canonical** - Canonical treatment of charm combined with the grand canonical treatment of baryon number, electric charge, and strangeness. Only Boltzmann statistics for charm particles.

---

**NOTE**

Analytic calculations within the strangeness-canonical ensemble in the presence of the EV/vdW interactions are implemented, but use approximations. 
More specifically, it is assumed that strange particles form a very small part of the whole system. This is appropriate for low energies (e.g. HADES), but not for high energies (like LHC). Use with caution!

---

**Resonance widths.** The prescription to treat finite resonance widths in spectral functions of resonances.

The possibilities are:

1. Zero-width approximation
2. Energy-independent Breit-Wigner, restricted to a ±2Γ interval around the pole mass
3. Energy-dependent Breit-Wigner (eBW)
4. Energy-dependent Breit-Wigner for the spectral function, but constant branching ratios for evaluating feeddown

**Statistics.** Maxwell-Boltzmann or quantum statistics. 

For quantum statistics the default method (*Use quadratures* checkbox checked) is to use Gauss-Legendre quadratures to compute the relevant integrals numerically. A faster method is the series (cluster) expansion in Bessel functions. This method however does not work well if chemical potential is close to or larger than particle's mass. 

It is possible to include quantum statistics for pions or for mesons only. This can speed up significantly calculations in the canonical ensemble.

**Conservation laws.** 
Chemical potentials in the GCE can be constrained by various conservation laws. The following options are possible.

1. Baryochemical potential can be constrained to reproduce a particular entropy per baryon ratio, S/B.
2. Electric charge chemical potential can be constrained to reproduce a particular electric-to-baryon charge ratio, Q/B. For heavy-ion collisions it is typically Q/B = 0.4.
3. Strange chemical potential can be constrained from the strangeness neutrality condition.
4. Charm chemical potential can be constrained from the charm neutrality condition.

**Excluded volume/van der Waals interactions.** 
Set the interaction parameters for an excluded volume or a van der Waals HRG model. The are two possibilities.

1. Set parameters manually, by providing the values of the a and b parameters. Interpreting these input parameters as those of protons, it is possible to assign the same parameters to all other particles, or impose a mass-proportional (bag model) or baryon conent proportional scaling of parameters. For Crossterms EV or QvdW models it is also possible to switch off interactions for various hadron pairs such as meson-meson, meson-baryon baryon-baryon, and baryon-antibaryon.
2. The interaction parameters can be read from an external file. 
   * For *Diagonal EV* model the external file should contain a list of rows (one for each species) with two columns: 1 - PDGID, 2 - v_i. 
   * For *Crossterms EV*, one row for each pair of species, three columns: 1 - PDGID1, 2 - PDGID2, 3 - b_ij. 
   * For *QvdW-HRG*, one row for each pair of species, four columns: 1 - PDGID1, 2 - PDGID2, 3 - b_ij, 4 - a_ij. 
   * If some species/pair of species is not found in the input file, then the parameters are assumed to be zero for those. A number of sample input files is provided in the `$(ThermalFIST)/input/list/interaction` folder.

The button *Calculate* invokes a calculation of all the hadron yields, equation of state properties, and, optionally, fluctuations, for the specified values of thermal parameters.

The button *Equation of state...* shows the equation of state properties corresponding to a given calculation, as well as some extra quantities, such as the second-order fluctuations incuding contributions from probabilistic decays.

Double-clicking on a particle in the list will open a window with additional useful information about the particle, including, e.g., all the feeddown contributions from different resonance decays.


#### *Thermal fits* tab

Thermal fitting. The HRG model specification is done in the same manner as in the previous tab.

It is possible to specify which parameters should be fitted, their initial/fixed values, and the value range where the fit is performed.
For fits in the (partially) canonical picture it also possible to fix the value of the correlation volume $V_c$ relative to the total volume $V$.

The experimental data can be loaded from an external file. Some samples are provided in the [$(ThermalFIST)/input/data](../input/data) folder. Possibility to input the data directly in the GUI is provided as well (double-click on the yield in the table to edit it).
It is possible to view the fit results directly within the GUI in a form of some of the common thermal fit plots.
Another useful feature is the analysis of the $χ^2$ profiles (the *Chi2 profile...* button), which can be done after a global fit was performed.

#### *Equation of state* tab

Offers a possibility to study the temperature dependence at fixed μ_B (or at fixed S/B) 
of some common equation of state observables, conserved charges susceptibilities, and number densities (primordial or with feeddown). It is also possible to consider ratios of any pair of these observables.
For a finite μ_B or S/B value, 
the chemical potentials μ_Q, μ_S, and/or μ_C can be fixed by conservation laws, or set to zero otherwise.
The thermal model specification is done similarly to the previous two tabs, except that here only the grand-canonical ensemble is considered.
For a number of observables at $μ_B = 0$ the published lattice QCD data of the Wuppertal-Budapest and/or HotQCD collaborations is plotted along with the calculation results for convenience.

#### *Event generator* tab

A Monte Carlo generator which generates events with particle numbers distributed according to the corresponding multiplicity distribution in a HRG model, while momenta are distributed in accordance with the Blast-Wave model (spherically or cylindrically symmetric). Currently restricted to the case of the Maxwell-Boltzmann statistics.
In theory provides the same yield, fluctuations etc. as analytic calculations (in case of EV/vdW interactions a sufficiently large system volume is necessary), also applicable to study combined effect of EV/vdW interactions and exact charge conservation.
Can be useful, e.g., in detector response simulations.

#### *Particle list editor* tab

Allows to edit the particle list and decay channels on the fly.
Note that antiparticles are not listed, but are generated automatically from particles if at least one of the quantum numbers (baryon number, electric charge, strangeness, charm) is non-zero.
All changes can be saved to a file.

### Other notes

- Sample C++ macros which use Thermal-FIST as a library can be found in the [$(ThermalFIST)/src/examples](../src/examples) folder
  
- An example of using the Thermal-FIST library as a submodule is provided by [this repository](https://github.com/vlvovch/1807.02079)

- Since version 1.2.1 it is possible to run interactive Thermal-FIST sessions in Jupyter Notebook. See the [FIST-jupyter](https://github.com/vlvovch/FIST-jupyter) repository for an example.
  
- The canonical ensemble calculations use the method from [this paper](https://arxiv.org/abs/nucl-th/0112021) to perform analytically the integration over the baryon number fugacity whenever this is possible. This allows to significantly speed up the calculations. This method cannot be used if the quantum statisics is applied for baryons, or if there are multi-baryons (light nuclei) in the particle list. In this case direct numerical integration is performed, which can make calculations very slow.

### Common issues

- CMake cannot find Qt5 and therefore GUI cannot be built. The following message usually appears after running the `cmake` command: 

  > By not providing "FindQt5Widgets.cmake" in CMAKE_MODULE_PATH this project has asked CMake to find a package configuration file provided by "Qt5Widgets", but CMake did not find one.

  Please make sure that Qt5 is installed in the system. If CMake still cannot find it, try specifying the Qt5 directory to CMake explicitly

  ```bash
  cmake -DCMAKE_PREFIX_PATH=<path-to-qt5> ../
  ```
  
  The above solution can also be applied in the case when CMake finds an incompatible version of Qt5, e.g. an x86 Qt build for a x64 Thermal-FIST build configuration.

- Program outputs warnings about the Bose-Einstein condensation. This occurs if chemical potentials of bosons at some point exceed their mass. This sometimes happens in the process of constaining electric/strange chemical potentials through conservation laws. If these warnings eventually disappear and final values of chemical potentials are reasonable, things should be fine.

- Accuracy of the calculations in the canonical ensemble will eventually break down as the system volume is increased to a large value, such as that corresponding to central heavy-ion collisions. Accuracy of a canonical ensemble calculation can be verified by comparing the values of the exactly conserved charges computed in the model with the given ones.

- Canonical ensemble calculations seem to take forever. Canonical ensemble calculations are slow if the particle list contains multi-baryon states (light nuclei) or if quantum statistical effects for baryons are included. Unless required, it is recommended to use particle list without light nuclei and use Boltzmann approximation for baryons in the canonical ensemble.
