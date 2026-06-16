# K-N phase shifts — Gibbs–Arceo reference

Reference data for the kaon–nucleon (exotic S=+1) S-matrix / Beth–Uhlenbeck
clusters. **Implemented** in `MesonBaryonPhaseShifts.{h,cpp}` and the `KN_*`
channels (S- and P-waves, both isospins). Config: `KN.conf`.

**Source:** W. R. Gibbs, R. Arceo, *Partial-wave analysis of K⁺ nucleon
scattering*, Phys. Rev. C75 (2007) 054005 [arXiv:nucl-th/0611095]. Parameters are
Tables of that paper, transcribed from the arXiv LaTeX source (exact).

## Scope and physics
- K⁺N is exotic (S=+1, no q-qbar content): **no resonances**, so all waves are
  synthetic clusters. The I=1 S-wave is strongly repulsive; the I=0 S-wave is
  weakly repulsive (small negative scattering length); the I=0 P-waves show a
  spin-orbit splitting (P01 attractive, P03 repulsive) — the paper notes the T=0
  ℓ>0 phases are consistent with a pure spin-orbit potential.
- Implemented: S- and P-waves for both isospins (I=0: S01, P01, P03; I=1: S11,
  P11, P13). The small D/F waves (also in the paper) are omitted.
- A B=+1, S=+1 fermionic cluster; constituent fugacity γq⁴γS (K = 1 light + 1
  strange, N = 3 light).

## Conventions
- `q` = CM momentum [GeV/c]; threshold √s₊ = m_K + m_N ≈ 1.434 GeV.
- Coefficients are "in (GeV/c) powers" as tabulated below.
- **Elastic range:** below the K-π-N threshold m_K + m_π + m_N ≈ 1.574 GeV
  (q ≈ 0.31 GeV/c); the integration of every K-N wave is cut there (`KN_Mmax`).

## Formula
S- and P-waves (paper Eq. for the S/P phase shifts):

    delta_{l±} = atan[ a q^(2l+1) / (1 + b1 q² + b2 q⁴ + b3 q⁶) ]

(S-wave uses b1,b2; P1/2 uses b1; P3/2 uses b1,b2,b3.) Non-resonant: the phases
stay within (−90°, 90°), so a plain `atan` (no branch tracking) is used.

## Parameters — a [(GeV/c)^−(2l+1)], b_i [(GeV/c)^−2i]

| wave | I | l | a | b1 | b2 | b3 |
|------|---|---|---|----|----|----|
| S01  | 0 | 0 | −0.531  | −1.206 | 1.362    |        |
| P01  | 0 | 1 | 23.765  | 3.690  |          |        |
| P03  | 0 | 1 | −3.808  | 2.919  | −10.042  | 212.021 |
| S11  | 1 | 0 | −1.562  | −1.108 | 0.217    |        |
| P11  | 1 | 1 | −12.002 | 31.139 |          |        |
| P13  | 1 | 1 | 13.357  | 126.676| −666.951 | 1276.123 |

(The omitted D-waves use δ = q⁵(c + d q² + e q⁴), F-waves δ = q⁷ f; coefficients
in the paper's Tables if ever needed.)

## Validation (these formulas + parameters, δ in degrees)
I=1 S-wave strongly repulsive; I=0 P-waves spin-orbit split (P01 +, P03 −):

| q [GeV/c] | √s [GeV] | S01 | S11 | P01 | P03 | P11 | P13 |
|-----------|----------|-----|-----|-----|-----|-----|-----|
| 0.10 | 1.450 | −3.1 | −9.0  | 1.3  | −0.2 | −0.5 | 0.4 |
| 0.20 | 1.494 | −6.3 | −18.1 | 9.4  | −1.6 | −2.5 | 1.2 |
| 0.30 | 1.565 | −10.0| −27.5 | 25.7 | −4.4 | −4.9 | 2.6 |

### Reproduce
```python
import math
mK, mN = 0.496, 0.938272
sthr=(mK+mN)**2; sdif=(mK-mN)**2
def qcm(s): return math.sqrt((s-sthr)*(s-sdif))/(2*math.sqrt(s))
def sqrts(q): return math.sqrt(mK*mK+q*q)+math.sqrt(mN*mN+q*q)
def dSP(q,l,a,b1=0,b2=0,b3=0):
    q2=q*q; ql=q if l==0 else q*q*q
    return math.atan(a*ql/(1.+q2*(b1+q2*(b2+q2*b3))))
waves={'S01':(0,-0.531,-1.206,1.362,0),'S11':(0,-1.562,-1.108,0.217,0),
 'P01':(1,23.765,3.690,0,0),'P11':(1,-12.002,31.139,0,0),
 'P03':(1,-3.808,2.919,-10.042,212.021),'P13':(1,13.357,126.676,-666.951,1276.123)}
r2d=180/math.pi
for q in (0.1,0.2,0.3):
    print(round(sqrts(q),3),{n:round(dSP(q,*p)*r2d,1) for n,p in waves.items()})
```

## Implementation notes
- Six single-wave channels `KN_S01/P01/P03/S11/P11/P13` (synthetic), `FamilyKN`.
- Constituents K{+1:321,−1:311} ⊗ N{+1:2212,−1:2112}; I=0 and I=1 via isospin CG.
- Orbital l in the `excitation` slot disambiguates S/P1/2 (both 2J+1=2); P3/2 has
  2J+1=4. Same baryon-encoding pattern as the π-N background waves.
