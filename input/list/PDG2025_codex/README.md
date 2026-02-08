# PDG2025 hadron list for Thermal-FIST

This directory contains a `PDG2025` hadron list generated from the PDG API
SQLite database (`pdglisting/pdg-2025-v0.2.2.sqlite`) using the `PDG2020`
Thermal-FIST charm-inclusive list (`../PDG2020/list-withcharm.dat`) as
baseline.

## Files

- `aux/generate_pdg2025.py`: generator script
- `aux/check_nndc_nuclei.py`: optional NuDat nuclei checker
- `aux/check_nndc_walletcards.py`: checker for NNDC Wallet Cards nuclei data
- `list.dat`: PDG2025 established-only hadron list (default)
- `decays.dat`: single canonical decays file containing all generated decays
- `list-withcharm.dat`: PDG2025 established-only hadron list including charm
- `list-extra.dat`: PDG2025 extra hadron list (baseline + expanded PDG2025 states)
- `list-withnuclei.dat`: established hadrons plus stable light nuclei
- `list-withexcitednuclei.dat`: established hadrons plus stable and excited light nuclei
- `pdglisting/`: auxiliary source files:
  - `mass_width_2025.mcd`
  - `pdg-2025-v0.2.2.sqlite`
  - `pdglive-established-nodes-2025.txt`
- `local-diagnostics/` (git-ignored): local audit artifacts and detailed
  cross-check tables used during list construction

## Generation policy

The generator uses a conservative baseline plus controlled expansion:

1. **Baseline state membership**: start from `PDG2020/list-withcharm.dat`, but remove
   baseline states that are absent in both PDG API 2025 `PART` listings
   (MCID- or name-matched) and the official 2025
   `pdglisting/mass_width_2025.mcd`
   Monte Carlo table.
2. **Masses and widths**:
   - prefer values from `pdglisting/mass_width_2025.mcd` when available,
   - fallback to PDG API 2025 summary values for states not covered by MCD
     (including MCID-less states via unique normalized-name match to `PART`
     parent nodes).
3. **State expansion**: add additional light/strange/charmed `PART` states from
   PDG2025 that have MCIDs and map to meson/baryon content built from
   `u,d,s,c` quarks
   (excluding `K_L^0` and `K_S^0` to avoid duplicating neutral kaon weak
   eigenstates in this list convention).
4. **Decays**:
   - update BRs for channels with explicit, mappable `BFX` decay modes in the
     PDG API data,
   - keep baseline channels for generic/ambiguous modes,
   - for added states, use nearest-state decay templates with matching
     conserved charges (`B,Q,S,C`) when explicit mappings are incomplete,
   - apply an explicit `K(1)(1650)` override to the Viewer-style channel set
     (`K pi pi` and `K phi`) with charge-conserving representatives,
   - renormalize to 100% for each decaying parent.
5. **Thresholds**:
   - keep baseline threshold conventions for existing states,
   - for newly added unstable hadrons, recompute `threshold` from the lightest
     hadronic channel that satisfies exact strong-conservation (`B,Q,S,C`);
     channels with leptons/photons are excluded from threshold evaluation.
6. **Output variants**:
   - `established-with-charm`: `list-withcharm.dat`,
     subset built from `pdglisting/mass_width_2025.mcd` plus states marked established in
     `pdgLive` group tables
     (`pdglisting/pdglive-established-nodes-2025.txt`) using the table markers:
     bullet (`•`) or baryon star rating `***`/`****`; written to
     `list-withcharm.dat`,
  - `established-only` (default): `list.dat`, fully charm-free projection of
    the established-with-charm set (excludes open charm and hidden charm),
   - `extra`: `list-extra.dat`, includes states found in
     PDG API `PART` listings (MCID or name match) and in
     `pdglisting/mass_width_2025.mcd`,
   - `with stable nuclei`: `list-withnuclei.dat`, formed by adding the
     nucleus-only delta from `PDG2020/list-withnuclei-withcharm.dat` to the
     established PDG2025 hadron list,
   - `with excited nuclei`: `list-withexcitednuclei.dat`, formed by adding the
     excited-nucleus delta from `PDG2020/list-withexcitednuclei.dat`.
7. **Single decays file**: `decays.dat` contains decays for all generated
   hadrons (`list-extra.dat`) plus nuclei/excited nuclei states sourced from
   `PDG2020/decays.dat`; each block is normalized to 100%.
8. **Runtime compatibility**: Thermal-FIST decays loader applies decay blocks
   only for parents present in the loaded particle list, so one superset
   `decays.dat` can be used with all list variants.
9. **Validation workflow**: detailed pdgViewer/PDG-code cross-check tables are
   produced during development and kept in `local-diagnostics/` (not included in
   public release payload).

This keeps established Thermal-FIST list semantics while applying PDG2025
updates where API mapping is unambiguous and fills missing decay structure using
charge-conserving nearest-neighbor templates.

## Recent changes

- Introduced a charm-aware established list split:
  - `list.dat`: established charm-free states (no open or hidden charm).
  - `list-withcharm.dat`: established states including charm.
- Added open-charm expansion from PDG2025 and updated charm-sector masses,
  widths, and decay BRs from PDG API 2025 where mappable.
- Added fallback to `mass_width_2025.mcd` for established states whose API
  parent rows have missing summary mass/width entries (e.g. `f2(2150)0`).
- Added MCID-less mass/width updates via unique name-to-parent matching
  (e.g. `phi(2170)`, `pi2(1880)`, `Omega(2012)`).
- Fixed charm-content bookkeeping for generated states: `|C|` now uses MCID
  valence charm content (so hidden-charm states such as charmonia carry
  `|C|=2` while keeping net `C=0`).
- Switched mass/width MCD parser to accept Fortran-style values like `6.E-02`
  and prioritize MCD values whenever present.
- Added PDG upper-limit enforcement for widths (`limit_type='U'`), with
  reference-list rounding where applicable (e.g. `Lambda(c)(2625)+`).
- Switched threshold evaluation to a hadronic/strong policy for unstable
  newly added hadrons (exclude leptonic/radiative channels when computing
  `threshold`).
- Replaced `K(1)(1650)0/+` resonance-resolved decay blocks with an explicit
  Viewer-style two-channel model (`K pi pi`, `K phi`), normalized to 100%.
- Updated output ordering to follow a PDG2020-like layout: isospin partners are
  grouped together (e.g. `pi0`/`pi+`), non-charm sectors are listed first, and
  charm sectors are moved to the end.
- `list-withcharm.dat` excludes the double-charm state `Xi(cc)++` to keep the
  canonical charm list focused on hidden/open single-charm hadrons;
  `list-extra.dat` still retains it.
- Switched to a single canonical `decays.dat` (all decays) and removed
  per-variant decay files.
- Added local pdgViewer and charm-PDG-code audit workflow; detailed outputs are
  stored in `local-diagnostics/` and excluded from release tracking.
- Improved existing-state BR updates to avoid low-coverage partial-BR inflation
  when no new mapped channels are available.
- Added established Viewer states without DB MCIDs via controlled overrides in
  three sectors: charmonium, `D_s` excitations, and charmed baryons
  (including `D3*(2750)+`/`D3*(2750)0`).
- Aligned synthetic charm IDs and charge multiplets for manual overrides to the
  `PDG2025/list-withcharm.dat` conventions (e.g. `psi(4230)`, `chi(c1)(3872)`,
  `Sigma(c)(2800)0/+ /++`, `Xi(c)(2970)0/+`, `Xi(c)(3080)0/+`).
- For charm manual additions, if a matching ID exists in
  `PDG2025/list-withcharm.dat`, keep that reference threshold and use its width
  when API/MCD width is missing.
- For manual charm additions with reference IDs, align mass to
  `PDG2025/list-withcharm.dat` to avoid name-inference drift
  (e.g. `Sigma(c)(2800)`, `Xi(c)(2970)+`).
- Added an explicit `D(s1)(2536)` decay override that includes
  `D*(2007)0 K+` and keeps a conservative, normalized channel set.
- For charm overrides with missing official `J^P`, assigned likely quark-model
  `J` values and synchronized synthetic MCID last digits to `2J+1`:
  `Xi(c)(3055)0`/`Xi(c)(3080)0` -> `5/2`,
  `Sigma(c)(2800)0` -> `3/2`,
  `Omega(c)(3000)`/`Omega(c)(3050)`/`Omega(c)(3185)` -> `1/2`,
  `Omega(c)(3065)`/`Omega(c)(3090)`/`Omega(c)(3327)` -> `3/2`,
  `Omega(c)(3120)` -> `5/2`,
  `Xi(cc)++` -> `1/2`.
- For established Viewer baryons lacking extractable API summary masses,
  masses are initialized from the state label (e.g. `Xi(c)(2970)0 -> 2.970 GeV`)
  and widths default to zero.

## Nuclei consistency check

- The local PDG2025 API SQLite (`pdglisting/pdg-2025-v0.2.2.sqlite`) contains
  no entries with MCIDs matching the 38 nuclei/excited-nuclei IDs used here.
- The local PDG2025 Monte Carlo table (`pdglisting/mass_width_2025.mcd`) also
  contains no matching nuclei/excited-nuclei MCIDs.
- Therefore nuclei properties are taken from existing local reference databases:
  `PDG2020/list-withnuclei-withcharm.dat`,
  `PDG2020/list-withexcitednuclei.dat`, and `PDG2020/decays.dat`.
- Cross-check after generation: nuclei field values (`name`, `stable`, `mass`,
  `degeneracy`, `statistics`, `B`, `Q`, `S`, `C`, `|S|`, `|C|`, `width`,
  `threshold`) are identical to PDG2020 reference entries, and excited-nuclei
  decay-channel signatures match exactly.
- NNDC validation utilities are in `aux/`:
  `check_nndc_nuclei.py` and `check_nndc_walletcards.py`.
- Intermediate NNDC cache/input files are not kept in this folder after
  cleanup; the checkers regenerate transient artifacts under `pdglisting/`
  when rerun.
- Wallet Cards are the primary reference for ground/isomeric properties.

## Summary vs PDG2020

Reference baselines for comparison: `../PDG2020/list.dat` and
`../PDG2020/list-withcharm.dat`.

- `list.dat` (established-only, default): 242 -> 248 states (`+11`, `-5`)
- `list-withcharm.dat` (established): 291 -> 330 states (`+44`, `-5`)
- `list-extra.dat` (extra): 291 -> 354 states (`+63`, `-0`)
- `list-withnuclei.dat`: 248 -> 256 states (`+8` stable nuclei)
- `list-withexcitednuclei.dat`: 256 -> 286 states (`+30` excited nuclei)
- Charm coverage in `list-withcharm.dat`: 60 net-charm states (39 charm
  baryons) plus 22 hidden-charm states (`C=0`);
  `list.dat` excludes both net and hidden charm.
- Detailed per-state additions/removals and validation tables are kept in
  `local-diagnostics/` (git-ignored).

## Source reference

- PDG API docs: <https://pdg.lbl.gov/2025/api/index.html>
- API schema/docs: <https://pdgapi.lbl.gov/doc/>
- Citation used in the dataset metadata:
  S. Navas et al. (Particle Data Group), Phys. Rev. D 110, 030001 (2024) and
  2025 update.
