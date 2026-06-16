# π-N phase shifts — Roy–Steiner reference (for the future πN channel)

Reference data for the pion–nucleon S-matrix / Beth–Uhlenbeck clusters.
**Implemented** in `MesonBaryonPhaseShifts.{h,cpp}` and the `piN_*` channels:
the Δ(1232) (`piN_Delta`, P33, reuses the real Δ codes) plus all five non-resonant
background waves (`piN_S31`, `piN_S11`, `piN_P31`, `piN_P11`, `piN_P13`, synthetic
clusters). Config: `piN.conf`. This file is the vetted source for the coefficients
below and a numerical validation.

**Source:** M. Hoferichter, J. Ruiz de Elvira, B. Kubis, U.-G. Meißner,
*Roy–Steiner-equation analysis of pion–nucleon scattering*, Phys. Rept. 625 (2016) 1
[arXiv:1510.06039]. Parameters are Appendix D, Tables D.1–D.2 (transcribed from the
arXiv LaTeX source `piN_RS_v3.tex`, not OCR — exact).

## Scope and validity
- Six lowest s-channel partial waves solved: **S11, S31, P13, P11, P31, P33**
  (notation P_{2I,2J}; the `l±` of the paper: `+` → J=l+1/2, `−` → J=l−1/2).
- Valid from the πN threshold up to the matching point **√s_m = 1.38 GeV**
  (covers the Δ(1232) at √s = 1.232 GeV).
- The P33 (Δ) uses a **conformal** parametrization; all others a **Schenk** form.

## Conventions
- `q` = CM momentum [GeV]; `s` in GeV²; threshold `s₊ = (m_N + M_π)²`.
- Masses that reproduce the published curves: `m_N = 0.938272`, `M_π = 0.13957` GeV
  (so √s₊ = 1.0778 GeV). q(s) = √([s−(m_N+M_π)²][s−(m_N−M_π)²]) / (2√s).
- Threshold meaning of the leading coefficients: `A₀₊ = a₀₊` (S-wave scattering
  length), `Ã₁₊ = a₁₊` (P33 scattering volume).
- All coefficients are "in appropriate powers of GeV" (so `q`, `s` in GeV units).

## Formulas

S-waves (Eq. D.1):

    tan δ₀₊ = |q| · (A₀₊ + B₀₊ q² + C₀₊ q⁴ + D₀₊ q⁶ + E₀₊ q⁸) · (s₊ − s₀₊)/(s − s₀₊)

P-waves P13, P11, P31 — Schenk (Eq. D.2):

    tan δ₁± = |q|³ · (A₁± + B₁± q² + C₁± q⁴ + D₁± q⁶) · (s₊ − s₁±)/(s − s₁±)

P33 = Δ(1232) — conformal (Eq. D.3). NOTE the `1/Ã` is the **constant term inside
the bracket**, not a prefactor:

    cot δ = (1/|q|³) · (s − s₁₊)/(s₊ − s₁₊) · ( 1/Ã₁₊ + B̃₁₊ [w(s) − w₊] + C̃₁₊ [w(s) − w₊]² )
    w(s) = (√s − √(s̄₁₊ − s)) / (√s + √(s̄₁₊ − s)),   w₊ = w(s₊)

Branch-track the P33 phase through 90° with `atan2(1, cot δ)` (resonant), like K*(892).
The Schenk waves are non-resonant: `atan(tan δ)` is fine (δ stays in (−90°, 90°)).

## Parameters — S-waves (Table D.1)

| coeff | S11 (I=1/2) | S31 (I=3/2) |
|-------|-------------|-------------|
| A₀₊   | 1.217       | −0.6183     |
| B₀₊   | −1.879×10   | −1.831×10   |
| C₀₊   | 1.958×10²   | 3.090×10²   |
| D₀₊   | −1.235×10³  | −2.846×10³  |
| E₀₊   | 3.350×10³   | 9.529×10³   |
| s₀₊   | 2.494       | −1.809×10³  |

## Parameters — P-waves Schenk (Table D.2)

| coeff | P13 (1+, I=1/2) | P11 (1−, I=1/2) | P31 (1−, I=3/2) |
|-------|-----------------|-----------------|-----------------|
| A₁±   | −1.085×10       | −2.569×10       | −1.477×10       |
| B₁±   | −1.145×10       | 8.062×10²       | 1.467×10²       |
| C₁±   | 3.651×10²       | −4.214×10³      | −1.633×10³      |
| D₁±   | −1.052×10³      | 3.986×10⁴       | 6.508×10³       |
| s₁±   | 0.9639          | 0.9340          | 0.4081          |

## Parameters — P33 conformal (Table D.2, middle)

| Ã₁₊       | B̃₁₊         | C̃₁₊     | s₁₊    | √s̄₁₊ (s̄₁₊)     |
|-----------|-------------|---------|--------|------------------|
| 7.781×10  | −3.986×10⁻² | −0.3098 | 0.4509 | 1.540 (2.3716)   |

## Error-band parameters (optional; Eq. D.4, Gaussian propagation of A,B with ρ_AB)

- S11/S31: ΔA₀₊ = 1.433×10⁻² / 1.289×10⁻², ΔB₀₊ = 8.592×10⁻² / 0.1744, ρ₀₊ = 1.000 / −0.2584
- P13:     ΔA₁₊ = 0.5649, ΔB₁₊ = 1.896, ρ₁₊ = −1.000
- P33:     ΔÃ₁₊ = 1.257, ΔB̃₁₊ = 1.113×10⁻⁴, ρ₁₊ = 0.2094
- P11/P31: ΔA₁₋ = 1.800 / 0.7257, ΔB₁₋ = 2.650×10 / 5.507, ρ₁₋ = −0.2510 / −0.9882

## Elasticity (Beth–Uhlenbeck cutoff Mmax)
Inelasticity model (Eq. 5.1), s_inel = (m_N + 2M_π)² → √s_inel ≈ 1.217 GeV:
- **η ≈ 1 (elastic to 1.38 GeV): S31, P13, P33** — the Δ is elastic across its peak.
- **S11, P31, P11**: elastic only below √s_inel ≈ 1.217 GeV (ππN opens; S11→ηN;
  Roper P11(1440)). Cap Mmax at ~1.217 GeV, or treat the rest as inelastic (crude).

## Validation (these formulas + parameters, δ in degrees)
Reproduces the known πN phase shifts — S31 repulsive, S11 attractive, P-waves small,
P11 rising toward the Roper, and **P33 through 90° at the Δ (√s=1.232) → ~151° at 1.38**:

| √s [GeV] | 1.10 | 1.15 | 1.20 | 1.232 | 1.28 | 1.32 | 1.38 |
|----------|------|------|------|-------|------|------|------|
| S11      |  5.0 |  8.4 | 10.1 | 10.9  | 12.0 | 12.8 | 14.9 |
| S31      | −3.1 | −7.3 |−11.2 |−13.6  |−17.1 |−20.3 |−24.7 |
| P13      | −0.2 | −1.0 | −2.0 | −2.5  | −3.3 | −3.9 | −4.5 |
| P11      | −0.4 | −1.0 |  0.3 |  2.7  |  9.4 | 19.0 | 40.8 |
| P31      | −0.3 | −1.8 | −3.6 | −4.8  | −7.0 | −9.0 |−11.5 |
| P33 (Δ)  |  1.9 | 15.7 | 55.5 | 93.7  |127.2 |140.0 |151.1 |

### Reproduce the table
```python
import math
mN, Mpi = 0.938272, 0.13957
sp = (mN+Mpi)**2; sm_ = (mN-Mpi)**2
def qcm(s): return math.sqrt((s-sp)*(s-sm_))/(2*math.sqrt(s))
def schenkS(s,A,B,C,D,E,s0):
    q=qcm(s); q2=q*q
    return math.atan(q*(A+B*q2+C*q2**2+D*q2**3+E*q2**4)*(sp-s0)/(s-s0))
def schenkP(s,A,B,C,D,s1):
    q=qcm(s); q2=q*q
    return math.atan(q**3*(A+B*q2+C*q2**2+D*q2**3)*(sp-s1)/(s-s1))
def conf(s,At,Bt,Ct,s1,sbar):              # 1/At is the CONSTANT term
    q=qcm(s)
    w=lambda x:(math.sqrt(x)-math.sqrt(sbar-x))/(math.sqrt(x)+math.sqrt(sbar-x))
    x=w(s)-w(sp)
    return math.atan2(1.0,(1/q**3)*(s-s1)/(sp-s1)*(1/At + Bt*x + Ct*x*x))
waves={
 'S11':('S',dict(A=1.217,B=-1.879e1,C=1.958e2,D=-1.235e3,E=3.350e3,s0=2.494)),
 'S31':('S',dict(A=-0.6183,B=-1.831e1,C=3.090e2,D=-2.846e3,E=9.529e3,s0=-1.809e3)),
 'P13':('P',dict(A=-1.085e1,B=-1.145e1,C=3.651e2,D=-1.052e3,s1=0.9639)),
 'P11':('P',dict(A=-2.569e1,B=8.062e2,C=-4.214e3,D=3.986e4,s1=0.9340)),
 'P31':('P',dict(A=-1.477e1,B=1.467e2,C=-1.633e3,D=6.508e3,s1=0.4081)),
 'P33':('C',dict(At=7.781e1,Bt=-3.986e-2,Ct=-0.3098,s1=0.4509,sbar=1.540**2)),
}
r2d=180/math.pi
for nm,(t,p) in waves.items():
    f={'S':schenkS,'P':schenkP,'C':conf}[t]
    print(nm,[round(f(W*W,**p)*r2d,1) for W in (1.10,1.15,1.20,1.232,1.28,1.32,1.38)])
```

## Implementation notes (maps onto the existing machinery)
- P33 → resonant conformal wave like `PiK_delta_I12_P` (K*(892)): same `atan2`
  branch tracking; only the bracket form `1/Ã + B̃ x + C̃ x²` and `w(s)` differ.
- S/P backgrounds → one small new Schenk helper `tan δ = |q|^{2l+1} poly(q²) · pole`.
- Cluster: B=1, **Fermi statistics** (`statistics = +1`) → the `mn` alternating sign
  in `DensityCluster` is automatic; constituents are π + N (set `m1=M_π`, `m2=m_N`,
  and the constituent fugacity metadata gives gammaq^? gammaS^0 — N has 3 light
  quarks, π has 2 → gammaq^5 for a πN cluster).
- Decays: πN (isospin Clebsch-Gordan), feeddown to N + π.
- Reuse real PDG codes for the Δ multiplet (Δ++ 2224, Δ+ 2214, Δ0 2114, Δ− 1114)
  via `memberPdg`, like the K*(892); the background S/P waves are synthetic clusters.
