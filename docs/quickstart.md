## CuteHRG 0.1

### External dependencies

- Qt5 (for gui only): in theory, any of the Qt5 versions should be compatible
- ROOT (optional): Fits are done using Minuit2. If ROOT is found then the Minuit2 library from ROOT is used. Otherwise Minuit2 is compiled along with the package.

### Building

Project is generated with cmake.

Example, run from root project folder:
```bash
mkdir build
cd build
cmake ../
make
```

This will build the libraries in `build/lib`, the QtHRG gui in `build/bin`,
and a couple of sample macro in `build/bin/routines`.

### QtHRG

The standard analysis can be performed within the gui.

From the build folder, run with
```bash
./bin/QtHRG
```

For most part using the gui should be self-explanatory. 

The particle list is read from file.
Standard list is in the `$(CuteHRG)/input/list/PDG2014` folder.
There `list-mulibaryons.dat` contains hadrons plus multibaryons, while `list.dat` contains hadrons only. In the same there should also be a `decays.dat` file specifying all the decay channels.

The gui consists of three tabs:

#### Thermal model tab

Performs a single thermal model calculation for a given set of thermal parameters.

**HRG model:**

1. **Ideal** - the standard non-interacting gas of hadrons and resonances
2. **Diagonal EV** - the most widely used excluded-volume HRG. References: [nucl-th/9711062](https://arxiv.org/abs/nucl-th/9711062) and [nucl-th/9808012](https://arxiv.org/abs/nucl-th/9808012) 
3. **Crossterms EV** - generalized EV model. EV parameters specified for each pair of species, instead of for each specie only. In this way one can model differently, e.g., baryon-baryon and baryon-antibaryon repulsion. See [nucl-th/9906068](https://arxiv.org/abs/nucl-th/9906068) and [1606.06218](https://arxiv.org/abs/1606.06218) for details.
4. **QvdW-HRG** - multi-component vdW HRG. Both attractive and repulsive interaction parameters can be specified for any pair of species. See [1707.09215](https://arxiv.org/abs/1707.09215)

**Statistics:** Boltzmann or quantum statistics. For quantum statistics the default (and faster) method is to use the series expansion in Bessel functions. This method does not work well if chemical potential is close to or larger than particle's mass. In this case *Use quadratures* option should be used, to calculate the integrals numerically.

**Excluded volume/van der Waals:** 
Set the interaction parameters for EV/QvdW HRG models:

1. **Same for all**: constant radius parameter for all species.
2. **Bag-like**: mass-proportional eigenvolumes. The *Radius (fm)* value determines the nucleon radius in such a scheme.
3. **Point-like mesons**: eigenvolume is proportional to the absolute baryon number, i.e. mesons are point-like
4. **Custom**: The interaction parameters are read from external file. For *Diagonal EV* the external file should contain a list of rows (one for each species) with two columns: 1 - PDGID, 2 - v_i. For *Crossterms EV*, one row for each pair of species, three columns: 1 - PDGID1, 2 - PDGID2, 3 - b_ij. For *QvdW-HRG*, one row for each pair of species, four columns: 1 - PDGID1, 2 - PDGID2, 3 - b_ij, 4 - a_ij. If some species/pair of species is not found in the input file, then the parameters are zero for that. Some sample input files are provided in `$(CuteHRG)/input/list/PDG2014/interaction` folder.


#### Thermal fits tab

Thermal fitting. HRG model specification same as in previous tab.
If $mu_Q$ and $mu_S$ are not constrained, they are fitted as well.

The experimental data can be loaded from external file. Some samples are provided in `$(CuteHRG)/input/data`. Limited possibility for manual editing is provided as well.
Another useful feature is the analysis of the $chi^2$ profiles, done in separate window. Should be opened after a global fit was performed.

#### Event generator tab

A Monte Carlo generator which generates events with particle numbers distributed according to the corresponding multiplicity distribution in a HRG model, while momenta are distributed in accordance with the Blast-Wave model. An experimental feature, in theory should provide the same yield, fluctuations etc. as analytic calculations, still under development.