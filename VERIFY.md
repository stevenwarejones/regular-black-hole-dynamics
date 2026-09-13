# Verification

```bash
pip install -r requirements.txt
python run_checks.py     # symbolic-identity certificates (~2 min)
pytest -q                # all certificates + negative regression tests
```

The pipeline runs three certificates as independent subprocesses; any nonzero
exit fails the whole run. The primary checks are **symbolic identities**, not
sample evaluations.

Three distinct kinds of check appear, and are labeled as such in the output:

- **Symbolic identity (the proof).** Residuals are reduced to a single fraction
  and the numerator is required to be the zero polynomial in *independent* jet
  symbols. No tolerance, no sampling.
- **Sample regression (tripwire, not a proof).** Evaluation at fixed rational
  points with one fixed `q` profile — only catches gross breakage.
- **Numerical `theta`-independence (sanity).** A generous-tolerance float check
  that the reduced density carries no residual `theta` dependence.

## 1. `proofs/verify_exact_solution.py` — primary

For arbitrary smooth `q(v)` (independent jets `Q0, Q1, …`) and constant `M`, all
eight reduced Euler–Lagrange residuals (`g00, g01, g11, Rf, A0, A1, B0, B1`) are
**identically zero**. **Expected:** each `IDENTICALLY ZERO`; final `PASS`.

## 2. `proofs/verify_mass_scope.py` — scope

With `M → M(v)` and independent mass jets `P0, P1, …`: `E[g11] − (−2·M′(v))` is
**identically zero** (the full source, derived not assumed; no higher-derivative
terms), and every **other** equation is identically zero even for `M′(v) ≠ 0`.
So the source lands in the **`g11` (= `g_rr`) equation**, proportional to
`M′(v)`. **Expected:** `PASS`. The coefficient `−2` is in the reduced-density
normalization; the physical `1/(4π)` of `T_vv` is *not* asserted here.

## 3. `proofs/invariant_distinction.py` — static-background gauge argument

Perturbing a static background `q = q0 + ε p(v)` with `δRf = 0`, the first
variation of `(∇r)²` equals `4 M ℓ² r² / (r³ + 2ℓ²q0)² · p`, checked
symbolically. Because `δRf = 0` forces `ξ^r = 0` and the static scalar is
`v`-independent, a pure-gauge mode would give zero; the nonzero variation (for
`M > 0`) shows the perturbation is not pure gauge. **Expected:** `PASS`. A
non-gauge perturbation is not by itself a propagating mode or an instability.

## Negative regression tests (`tests/`)

`pytest` runs all three certificates **and** negative tests confirming the
identity checker *rejects* residuals that vanish only on the old fixed sample
profile — `q'' − 1/5`, `q'''`, a nonlinear `M'^2`, and a higher derivative `M''`.
This is what distinguishes a genuine identity proof from finite sampling.

## Reproducibility notes

- `requirements.txt` gives **lower bounds**, not a lock. For a citable run,
  record `pip freeze` (and the Python version) alongside the log; a pinned
  constraints file can be added if exact reproduction is required.
- The exact steps use `sympy.Rational` / exact polynomial arithmetic; there is
  no floating-point tolerance to tune except in the `theta`-independence sanity
  check, which is explicitly not part of the proof.
- `VERIFICATION_LOG.txt` in the repo root is a captured run of `run_checks.py`
  and `pytest` (regenerate it with `python run_checks.py > VERIFICATION_LOG.txt`
  and appending `pytest -q`). It is a convenience transcript, not a proof.
- `hashes.txt` is a transfer check only, not evidence for the equations;
  regenerate with `python make_hashes.py` after content changes.
