## Thermal-FIST: Quick start guide

### External dependencies

- Qt5 [for GUI (QtThermalFIST) only]: in theory, any of the Qt5 versions should be compatible
- ROOT [optional]: Fits are done using Minuit2. If ROOT (containing Minuit2) is found, then the Minuit2 library from ROOT is used. Otherwise Minuit2 is built inside the Thermal-FIST package.

### Building

The preferred way is to use **cmake**.

For example, run the following commands from root project folder to build the package on a Linux-like system:
```bash
mkdir build
cd build
cmake ../
make
make install
```

This will build the libraries in `build/lib`, the QtThermalFIST GUI in `build/bin`,
and a couple of sample macros in `build/bin/routines`.

Please see also the commands listed in the [.travis.yml](../.travis.yml) file for another example of building and running the package under the Linux system.

### QtThermalFIST

The standard analysis can be performed within the GUI.

From the build folder, run with
```bash
./bin/QtThermalFIST
```

Using the GUI should be self-explanatory for most part. 

The particle list is read from file.
The standard list is in the `$(ThermalFIST)/input/list/PDG2014` folder, containing particles consisting of light and strange quarks in accordance with the PDG2014 compilation.
There the `list-withnuclei.dat` file contains hadrons plus multibaryons, while the `list.dat` file contains hadrons only. In the same folder there is also the `decays.dat` file specifying all the decay channels.
Additionally, the `$(ThermalFIST)/input/list/` folder also contains particle lists extracted from the open source THERMUS-2.3 and THERMUS-3.0 ([link](http://www.phy.uct.ac.za/phy/people/academic/wheaton/research)) codes, which can be used to compare the results between different lists and also compare Thermal-FIST to THERMUS.

The **QtThermalFIST** GUI program consists of four tabs:

#### Thermal model tab

Performs a single thermal model calculation for a given set of thermal parameters.

**HRG model:**

1. **Ideal** - the standard non-interacting gas of hadrons and resonances
2. **Diagonal EV** - the most widely used excluded-volume HRG. References: [nucl-th/9711062](https://arxiv.org/abs/nucl-th/9711062) and [nucl-th/9808012](https://arxiv.org/abs/nucl-th/9808012) 
3. **Crossterms EV** - generalized EV model. EV parameters specified for each pair of species, instead of for each specie only. In this way one can model differently, e.g., baryon-baryon and baryon-antibaryon repulsion. See [nucl-th/9906068](https://arxiv.org/abs/nucl-th/9906068) and [1606.06218](https://arxiv.org/abs/1606.06218) for details.
4. **QvdW-HRG** - multi-component vdW HRG. Both attractive and repulsive interaction parameters can be specified for any pair of species. See [1707.09215](https://arxiv.org/abs/1707.09215)

**Statistics:** Boltzmann or quantum statistics. For quantum statistics the default method (*Use quadratures* checkbox checked) is to use Gauss-Legendre quadratures to compute the relevant integrals numerically. A faster method is the series expansion in Bessel functions. This method however does not work well if chemical potential is close to or larger than particle's mass. 

**Excluded volume/van der Waals:** 
Set the interaction parameters for EV/QvdW HRG models:

1. **Same for all**: constant radius parameter for all species.
2. **Bag-like**: mass-proportional eigenvolumes. The *Radius (fm)* value determines the nucleon radius in such a scheme.
3. **Point-like mesons**: eigenvolume is proportional to the absolute baryon number, i.e. mesons are point-like
4. **Custom**: The interaction parameters are read from external file. 
   * For *Diagonal EV* the external file should contain a list of rows (one for each species) with two columns: 1 - PDGID, 2 - v_i. 
   * For *Crossterms EV*, one row for each pair of species, three columns: 1 - PDGID1, 2 - PDGID2, 3 - b_ij. 
   * For *QvdW-HRG*, one row for each pair of species, four columns: 1 - PDGID1, 2 - PDGID2, 3 - b_ij, 4 - a_ij. 
   * If some species/pair of species is not found in the input file, then the parameters are zero for that. Some sample input files are provided in `$(ThermalFIST)/input/list/PDG2014/interaction` folder.

The button *Equation of state...* shows the equation of state properties corresponding to a given calculation, as well as some extra quantities, such as the second-order fluctuations incuding the contributions from probabilistic decays.

Double-clicking on a particle in the list will open a window with additional useful information about the particle, including e.g. all the contributions from different decays.

**Important notice:** Analytic calculations within the strangeness-canonical ensemble in the presence of the EV/vdW interactions are done approximately, assuming that strange particles are a very small part of the whole system. This is appropriate for low energies (e.g. HADES), but not for high energies (like LHC). Use with caution!

#### Thermal fits tab

Thermal fitting. HRG model specification same as in previous tab.
If $mu_Q$ and $mu_S$ are not marked to be constrained, they are fitted as free parameters.

The experimental data can be loaded from external file. Some samples are provided in the [$(ThermalFIST)/input/data](../input/data) folder. Possibility to input the data directly in GUI is provided as well (double-click on the yield in the table).
Another useful feature is the analysis of the $chi^2$ profiles (the *Chi2 profile...* button), done in a separate window. Should be opened after a global fit was performed.


#### Event generator tab

A Monte Carlo generator which generates events with particle numbers distributed according to the corresponding multiplicity distribution in a HRG model, while momenta are distributed in accordance with the Blast-Wave model (spherically or cylindrically symmetric). Currently restricted to the case of the Maxwell-Boltzmann statistics.
In theory provides the same yield, fluctuations etc. as analytic calculations (in case of EV/vdW interactions a sufficiently large volume is necessary), also applicable to study combined effect of EV/vdW interactions and exact charge conservation.

#### Particle list editor tab

Allows to edit the particle list and decay channels on the fly.
Note that antiparticles are not listed, but are generated autmatically from particles is at least one of the quantum numbers (baryon number, electric charge, strangeness, charm) is non-zero.
All changes can be saved to file.

### Other notes

- Sample C++ macros which use Thermal-FIST as a library can be found in the [$(ThermalFIST)/src/routines](../src/routines) folder
- An example of using the Thermal-FIST library as a submodule can be found [here](https://github.com/vlvovch/1807.02079)

### Common issues

- CMake cannot find Qt5 and therefore GUI cannot be built. Usually this type of message appears after running the `cmake` command: 

> By not providing "FindQt5Widgets.cmake" in CMAKE_MODULE_PATH this project has asked CMake to find a package configuration file provided by "Qt5Widgets", but CMake did not find one.

Make sure that Qt5 is installed in the system. If CMake still cannot find it, try specifying the Qt5 directory explicitly

```bash
cmake -DCMAKE_PREFIX_PATH=<path-to-qt5> ../
```
