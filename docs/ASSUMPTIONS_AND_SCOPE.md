# Assumptions and scope

What the certificates establish, stated as narrowly as the calculation warrants.

## Established (symbolic identities, arbitrary smooth q)

- For **constant `M`** and **arbitrary smooth `q(v)`**, the displayed fields
  solve all eight reduced Euler–Lagrange equations on `r > 0` patches, with **no
  added matter**. Proved as a symbolic identity in independent regulator jets,
  not from samples. (`proofs/verify_exact_solution.py`)
- About a **static background**, the perturbation `q = q0 + ε p(v)` (with
  `δRf = 0`) changes the scalar `(∇r)²` by a nonzero first variation, so — since
  `δRf = 0` forces `ξ^r = 0` — it is **not** pure gauge. Requires `M > 0`.
  (`proofs/invariant_distinction.py`)
- The vacuum statement is **specific to constant `M`**: with `M -> M(v)`,
  exactly the `g_rr` equation gains a source `= -2 M'(v)` (reduced normalization)
  and the other seven stay vacuum. (`proofs/verify_mass_scope.py`)

## Explicit assumptions carried

- `q(v) > 0`, `ell != 0`, `r > 0`.
- Statements are **local**, on patches with `a0 b0 != 0`. They need not contain
  the center, and a horizon `f = 0` is not by itself excluded.
- A time-dependent `M(v)` sources matter; in the source model this is
  `T_vv = M'(v) / (4 pi r^2)`. This repository certifies the reduced-density
  identity `E[g11] = -2 M'(v)` (source in the `g_rr` equation only), not the
  physical `1/(4 pi)` normalization, which we do not independently fix here.
- Domain restrictions here are those actually needed for the exact-family proof
  (`r > 0`, `ell != 0`, real EF branch); we do not import the tighter loci from a
  separate hyperbolicity calculation. `a0 b0 != 0` is a conservative patch
  choice for the local framing, not a denominator required by these identities.

## NOT claimed by this repository

- No claim of a global or non-spherical well-posed evolution.
- No stability, instability, or mass-inflation claim.
- No statement about singularity resolution, quantum horizon structure, or any
  observational signal.
- No refutation of any result in the literature. See
  `RELATIONSHIP_TO_PRIOR_WORK.md` for the (deliberately neutral) relationship to
  the source paper's perturbative analysis.
- No priority claim. See `PROVENANCE.md`.
