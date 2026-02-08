# PDG2025 Particle List for Thermal-FIST

This directory contains the PDG2025 hadron resonance list and decay channels for use with Thermal-FIST.

## Contents

| File | Description |
|------|-------------|
| `list.dat` | PDG2025 particle list — 248 light+strange hadrons (444 with antiparticles) |
| `list-withcharm.dat` | Extended list with charm — 330 particles (586 with antiparticles). C ∈ {-1, 0, 1}, compatible with charm-canonical ensemble |
| `list-withnuclei.dat` | Extended list with stable light nuclei — 256 particles |
| `list-withexcitednuclei.dat` | Extended list with stable and excited light nuclei — 286 particles |
| `decays.dat` | Decay channels for all unstable hadrons |
| `generate_pdg2025.py` | Generation script: takes PDG2020 as base and applies PDG2025 updates |
| `pdglisting/mass_width_2025.txt` | PDG2025 mass-width summary data used by `generate_pdg2025.py` |
| `crosscheck/generate_pdg2025_from_pdfs.py` | Independent cross-check script: encodes all particle properties directly from PDG2025 Summary Tables |
| `crosscheck/list_from_pdfs.dat` | Output of cross-check script (validated to match `list.dat`) |
| `crosscheck/decays_from_pdfs.dat` | Decay output of cross-check script (validated to match `decays.dat`) |

## How the PDG2025 List Was Generated

### Primary script: `generate_pdg2025.py`

This script takes the PDG2020 list (`../PDG2020/list.dat`, `../PDG2020/decays.dat`) as a base and applies the following updates:

1. **Removed 5 particles** not established in PDG2025:
   - pi(1)(1400)0 (9000113), pi(1)(1400)+ (9000213) -- removed from PDG summary tables
   - Sigma(2250)-, Sigma(2250)0, Sigma(2250)+ (33114, 33214, 33224) -- demoted from summary tables

2. **Added 11 new particles** (7 new states, some with multiple charge states):
   - f(2)(1565) -- new I=0 tensor meson
   - f(0)(2020) -- new I=0 scalar meson
   - f(2)(2150) -- new I=0 tensor meson
   - K(0)*(700) / kappa -- light scalar kaon (K0, K+)
   - K(1)(1650) -- excited axial kaon (K0, K+)
   - K(0)*(1950) -- excited scalar kaon (K0, K+)
   - K(2)*(1980) -- excited tensor kaon (K0, K+)

3. **Updated masses and widths** from PDG2025 data for all existing particles, using `pdglisting/mass_width_2025.txt`.

4. **Updated decay branching fractions** for select particles based on PDG2025 values, including refined isospin decomposition of baryon decay channels with Clebsch-Gordan coefficients.

### Cross-check script: `generate_pdg2025_from_pdfs.py`

This script independently encodes all 248 particle properties (mass, width, quantum numbers, PDG IDs) directly from the PDG 2025 Summary Table PDFs, without referencing the PDG2020 base list. It serves as a validation that the primary generation script produced correct results.

Running it confirms:
- All 248 particles match in mass, width, degeneracy, and quantum numbers (B, Q, S)
- Thermodynamic quantities (density, pressure, energy) computed with Thermal-FIST match exactly between the two lists

## Particle Count by Sector (list.dat)

| Sector | Particles | Notes |
|--------|-----------|-------|
| Light unflavored mesons | 70 | pi through f(2)(2340), including hidden strangeness (eta, phi, ...) |
| Strange mesons | 32 | K through K(4)*(2045), neutral and charged states |
| N baryons | 40 | p, n, N(1440) through N(2600), with 0/+ charge states |
| Delta baryons | 48 | Delta(1232) through Delta(2420), with -/0/+/++ charge states |
| Lambda baryons | 14 | Lambda through Lambda(2350) |
| Sigma baryons | 27 | Sigma+/0/- through Sigma(2030), with 3 charge states |
| Xi baryons | 12 | Xi0/- through Xi(2030), with 2 charge states |
| Omega baryons | 3 | Omega-, Omega(2012), Omega(2250) |
| **Total** | **248** | Antiparticles generated automatically (444 total) |

---

## list-withcharm.dat — Charm-Extended Particle List

### Overview

`list-withcharm.dat` extends `list.dat` with 82 charmed hadrons from the PDG2025 Review of Particle Physics Summary Tables. All entries satisfy C ∈ {-1, 0, 1}, making the list compatible with charm-canonical ensemble calculations.

The doubly charmed baryon Xi_cc++ (C=2, pdgid 4422) is intentionally excluded to maintain this compatibility.

| | Count |
|---|-------|
| Light + strange hadrons | 248 |
| Charmed hadrons | 82 |
| **Total particles** | **330** |
| **Total with antiparticles** | **586** |

### Charm Sector Breakdown

| Sector | Entries | Notes |
|--------|---------|-------|
| Charmonium (cc̄) | 22 | eta_c through psi(4660); C=0, \|C\|=2 |
| D mesons | 13 | D through D*_3(2750); C=±1, S=0 |
| D_s mesons | 8 | D_s through D*_s3(2860); C=±1, S=±1 |
| Lambda_c | 6 | Lambda_c+ through Lambda_c(2940)+ |
| Sigma_c | 9 | Sigma_c(2455) through Sigma_c(2800), 3 charges each |
| Xi_c | 15 | Xi_c through Xi_c(3080), 1–2 charges |
| Omega_c | 9 | Omega_c through Omega_c(3327) |
| **Total charm** | **82** | |

### States with Assigned J^P

Several charmed baryons appear in the PDG2025 RPP with J^P = ??. Degeneracy (2J+1) was assigned using quark model predictions, following the same approach used for non-charm states with unknown J^P (e.g., Xi(1950), Xi(2030), Omega(2250)). Reasoning is documented in `#` comment blocks in the data file.

| State | Assigned J^P | deg | Reasoning |
|-------|-------------|-----|-----------|
| Sigma_c(2800) | 3/2- | 4 | 1P excitation; broad width (~70 MeV) favors lower J over 5/2-. Analog: Sigma(1670) |
| Xi_c(3055)+ | 3/2+ | 4 | D-wave (L=2) excitation; narrow width (7.8 MeV) |
| Xi_c(3080) | 5/2+ | 6 | D-wave partner of Xi_c(3055); narrower width (3.6–5.6 MeV) |
| Omega_c(3000) | 1/2- | 2 | 1P quintuplet (LHCb 2017); QCD sum rules (Wang et al., EPJC 77, 2017) |
| Omega_c(3050) | 1/2- | 2 | 1P quintuplet |
| Omega_c(3065) | 3/2- | 4 | 1P quintuplet |
| Omega_c(3090) | 3/2- | 4 | 1P quintuplet |
| Omega_c(3120) | 5/2- | 6 | 1P quintuplet |
| Omega_c(3185) | 1/2+ | 2 | 2S radial excitation; broad width (50 MeV) favors lower J |
| Omega_c(3327) | 3/2+ | 4 | 1D excitation; moderate width (20 MeV) |

**Note on the Omega_c P-wave quintuplet:** Some quark models swap J^P assignments between individual states within the quintuplet; the total degeneracy sum (2+2+4+4+6=18) across all five states is robust regardless of the specific ordering.

### MC ID Convention for Charm Entries

- **Standard PDG IDs** (e.g., 441, 443, 4122, 4232): official MC numbering
- **PDG IEXTRA IDs** (9000000–9099999, e.g., 9000443 for psi(4040)): official PDG assignments for overflow states
- **Custom IDs** (90XXXXX): for states without official MC IDs. Format: `9` + counter + quark-content digits + `2J+1` suffix. Examples: `9004222` = Sigma_c(2800)++, `9104332` = Omega_c(3000), `9204324` = Xi_c(3080)+

---

## PDG Monte Carlo ID Numbering Scheme

### Standard PDG IDs (218 particles in list.dat)

Most particles use the official PDG Monte Carlo numbering scheme, which encodes quark content and quantum numbers in digit positions:

**Mesons:** `n_r n_L n_q1 n_q2 2S+1` (up to 6 digits)
- Example: `111` = pi0, `113` = rho(770)0, `10113` = b(1)(1235)0, `100113` = rho(1450)0

**Baryons:** `n_r n_q1 n_q2 n_q3 2J+1` (up to 5 digits, 6 for excited states)
- Example: `2212` = proton, `2112` = neutron, `12212` = N(1440)+, `1114` = Delta(1232)-

**PDG "IEXTRA" range (9000000-9099999):** Used for particles that cannot be accommodated in the standard digit-encoding scheme (e.g., very broad states, additional states with same quantum numbers). These IDs are assigned by PDG.
- Examples: `9000221` = f(0)(500), `9010221` = f(0)(980), `9030221` = f(0)(1500)

### Charge state patterns for standard baryons

For baryons with multiple charge states, the standard PDG scheme increments specific digit positions:

| Particle | I | Q=-1 | Q=0 | Q=+1 | Q=+2 |
|----------|---|------|-----|------|------|
| N(939) | 1/2 | -- | 2112 | 2212 | -- |
| Delta(1232) | 3/2 | 1114 | 2114 | 2214 | 2224 |
| Sigma | 1 | 3112 | 3212 | 3222 | -- |
| Xi | 1/2 | 3312 | 3322 | -- | -- |

For excited states with the same pattern:
| N(1440) | 1/2 | -- | 12112 | 12212 | -- |
| N(1520) | 1/2 | -- | 1214 | 2124 | -- |
| Delta(1600) | 3/2 | 31114 | 32114 | 32214 | 32224 |

The pattern: changing the leading quark-content digits by `1000` or `100` shifts between u and d quarks, encoding the charge state.

### Custom Thermal-FIST IDs for light/strange baryons (30 particles)

Thirty high-mass baryon resonances cannot be accommodated in the standard 5-digit PDG scheme because all available digit combinations for their quantum numbers are already used by lower-mass states. These use 7-8 digit IDs in the `99xxxxx` range:

**Nucleon resonances (20 IDs, 10 states x 2 charges):**

| State | Q=0 ID | Q=+1 ID | Pattern |
|-------|--------|---------|---------|
| N(1875) | 9902114 | 9902214 | `99` + `0` + `2112/2212` + `4` |
| N(1880) | 9902112 | 9902212 | `99` + `0` + `2112/2212` + `2` |
| N(1895) | 9912112 | 9912212 | `99` + `1` + `2112/2212` + `2` |
| N(1900) | 9912114 | 9912214 | `99` + `1` + `2112/2212` + `4` |
| N(2060) | 9902116 | 9902216 | `99` + `0` + `2112/2212` + `6` |
| N(2100) | 9922112 | 9922212 | `99` + `2` + `2112/2212` + `2` |
| N(2120) | 9922114 | 9922214 | `99` + `2` + `2112/2212` + `4` |
| N(2220) | 99021110 | 99022110 | `99` + `0` + `21110/22110` |
| N(2250) | 99121110 | 99122110 | `99` + `1` + `21110/22110` |
| N(2600) | 99021112 | 99022112 | `99` + `0` + `21112/22112` |

**Delta resonances (8 IDs, 2 states x 4 charges):**

| State | Q=-1 | Q=0 | Q=+1 | Q=+2 |
|-------|------|-----|------|------|
| Delta(2200) | 9901118 | 9901218 | 9902128 | 9902228 |
| Delta(2420) | 99011112 | 99012112 | 99021212 | 99022212 |

**Strange baryons (2 IDs):**

| State | ID |
|-------|-----|
| Lambda(2350) | 99031210 |
| Omega(2012) | 9903334 |

The construction follows the pattern: `99` prefix + radial counter + standard baryon quark-content digits + `2J+1` suffix. The charge state is encoded using the same quark-digit shifting as in the standard scheme (`2112` <-> `2212` for N, `1112` <-> `1212` <-> `2122` <-> `2222` for Delta). All 30 custom IDs were inherited unchanged from the PDG2020 list used in Thermal-FIST.

## File Formats

### list.dat format (14 columns)

```
pdgid  name  stable  mass[GeV]  degeneracy  statistics  B  Q  S  C  |S|  |C|  width[GeV]  threshold[GeV]
```

- `stable`: 1 if treated as stable (very narrow width), 0 if decays are resolved
- `degeneracy`: 2J+1 (exception: f(0)(500) has degeneracy 0 to suppress its thermal contribution)
- `statistics`: -1 for mesons (Bose-Einstein), +1 for baryons (Fermi-Dirac)
- `|S|`, `|C|`: absolute strangeness/charm content (not net strangeness). E.g., eta has |S|=1.33, phi has |S|=2
- `threshold`: lightest decay channel mass sum (GeV)

### decays.dat format

```
pdgid                    # particle name
N                        # number of decay channels
BR  daughter1 daughter2  # branching ratio and daughter PDG IDs
...
```

Antiparticle decays are not listed; they are generated automatically by charge conjugation.

## Post-Generation Corrections

The following corrections were applied manually after generation, based on cross-checks against the PDG2025 API database and MC data files.

### Property fixes

| File | Particle | Property | Old | New | PDG2025 value |
|------|----------|----------|-----|-----|---------------|
| `list.dat`, `list-withcharm.dat` | pi2(1880) | mass | 1.895 GeV | 1.874 GeV | 1874 +26/-5 MeV |
| `list-withcharm.dat` | D(1)(2420)+ | mass | 2.4261 GeV | 2.4221 GeV | 2422.1 ± 0.5 MeV |
| `list-withcharm.dat` | h_c(1P) | width | 0.7 MeV | 0.78 MeV | 0.78 ± 0.28 MeV |

### Branching ratio fixes

| Particle | Channel | Old BR | New BR | PDG2025 value |
|----------|---------|--------|--------|---------------|
| D\*(2007)0 | D0 pi0 | 0.619 | 0.647 | 64.7 ± 0.9% |
| D\*(2007)0 | D0 gamma | 0.381 | 0.353 | 35.3 ± 0.9% |

### Charm decay catch-all channels

Branching ratios from the PDG are often incomplete: only experimentally measured exclusive modes are listed, so the sum of known BRs is frequently less than 100%. For thermal model feeddown calculations, BR sums must equal 1.0 to ensure proper normalization.

To handle this, an explicit **catch-all decay channel** is added to each charm particle, carrying the unmeasured fraction BR = 1.0 - (sum of known BRs). For particles with no measured exclusive modes, the catch-all carries 100%.

The catch-all channel uses the lightest kinematically allowed final state that conserves the appropriate quantum numbers:

| Particle type | Catch-all channel | Daughters | Rationale |
|---------------|-------------------|-----------|-----------|
| Charmonia (cc&#x0304;, C=0) | pi+ pi- pi0 | 211 -211 111 | Lightest hadronic state with correct quantum numbers |
| D0 (weak decay) | K- pi+ pi0 | -321 211 111 | Cabibbo-favored (c&#x2192;s) hadronic final state |
| D+ (weak decay) | K&#x0304;0 pi+ pi0 | -311 211 111 | Cabibbo-favored (c&#x2192;s) hadronic final state |
| D resonances (C=1, S=0) | D + pi | 421 211 or 421 111 | Lightest open-charm + pion (charge-appropriate) |
| D_s resonances (C=1, S=1) | D_s + pi0 | 431 111 | Lightest open-charm-strange + pion |
| Lambda_c, Sigma_c (C=1, S=0) | Lambda_c + pi | 4122 + pi | Lightest charmed baryon + pion |
| Xi_c (C=1, S=-1) | Xi_c + pi | 4232/4132 + pi | Lightest charmed-strange baryon + pion |
| Omega_c (C=1, S=-2) | Omega_c + pi0 | 4332 111 | Lightest charmed-double-strange baryon + pion |

**Summary of catch-all additions in `decays.dat`:**

- 13 charmonia with existing PDG decay channels: catch-all added for the unmeasured fraction
- 9 charmonia with no PDG exclusive modes: 100% catch-all
- 2 ground-state D mesons (D0, D+): catch-all for unmeasured weak decay fraction (2.1% and 12.0%)
- 7 D/D_s meson resonances: 100% catch-all
- 17 charmed baryon resonances (Lambda_c, Sigma_c, Xi_c, Omega_c): 100% catch-all

After these additions, all 82 charm particles (55 unstable + 27 stable) have BR sums = 1.0.

---

## Validation

### list.dat
1. Cross-checked all particle properties against independently encoded PDG2025 Summary Table data
2. Computed HRG thermodynamics at T=155 MeV, muB=0 with Thermal-FIST and confirmed exact agreement between original and independently generated lists

### list-withcharm.dat
1. Verified quantum number consistency (B, Q, S, C, |C|, statistics, degeneracy) for all 82 charm entries
2. Verified conservation of B, Q, S, C in all threshold decay channels
3. Confirmed threshold mass values match lightest allowed decay products within 1 MeV
4. Loaded all 586 particles in Thermal-FIST C++ and computed thermal densities at T=155 MeV
5. Confirmed C ∈ {-1, 0, 1} for all entries (charm-canonical compatibility)

### decays.dat (charm sector)
1. Verified all 55 unstable charm particles have BR sum = 1.000 (within 0.001 tolerance)
2. Verified all daughter particles exist in the particle list
3. Verified all decay thresholds ≤ parent mass
4. Charge (Q) conservation verified in all catch-all channels
