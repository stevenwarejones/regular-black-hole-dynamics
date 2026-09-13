# Assumptions and scope

What the certificates establish, stated as narrowly as the calculation warrants.

## Established (certified in exact arithmetic)

- For **constant `M`** and a smooth function `q(v)`, the displayed fields solve
  all eight reduced Euler–Lagrange equations on `r > 0` patches, with **no added
  matter**. (`proofs/verify_exact_solution.py`)
- The `q(v)` variation changes the scalar invariant `(nabla r)^2`, so it is
  **not** a pure coordinate transformation. (`proofs/invariant_distinction.py`)
- The vacuum statement is **specific to constant `M`**: with `M -> M(v)` the
  residuals are proportional to `M'(v)` and no longer vanish.
  (`proofs/verify_mass_scope.py`)

## Explicit assumptions carried

- `q(v) > 0`, `ell != 0`, `r > 0`.
- Statements are **local**, on patches with `a0 b0 != 0`. They need not contain
  the center, and a horizon `f = 0` is not by itself excluded.
- A time-dependent `M(v)` sources matter; in the source model this is
  `T_vv = M'(v) / (4 pi r^2)`. This repository certifies only the qualitative
  constant-`M` boundary of that statement, not the `1/(4 pi)` normalization,
  which we do not independently fix here.

## NOT claimed by this repository

- No claim of a global or non-spherical well-posed evolution.
- No stability, instability, or mass-inflation claim.
- No statement about singularity resolution, quantum horizon structure, or any
  observational signal.
- No refutation of any result in the literature. See
  `RELATIONSHIP_TO_PRIOR_WORK.md` for the (deliberately neutral) relationship to
  the source paper's perturbative analysis.
- No priority claim. See `PROVENANCE.md`.
