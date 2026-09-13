# Verification

```bash
pip install -r requirements.txt
python run_checks.py
```

The pipeline runs three certificates as independent subprocesses; any nonzero
exit fails the whole run. It checks the **science**, not file hashes.

## 1. `proofs/verify_exact_solution.py` — primary

Rebuilds the reduced Euler–Lagrange operator from the action and evaluates every
one of the eight residuals (`g00, g01, g11, Rf, A0, A1, B0, B1`) on the solution
at seven rational sample points in exact arithmetic. **Expected:** each residual
`ALL ZERO`; final `PASS`. A `theta`-independence sanity check runs first.

## 2. `proofs/verify_mass_scope.py` — scope

Promotes `M -> M(v)` and confirms, in exact arithmetic, that the `g00` residual
is `0` when `M'(v)=0` and proportional to `M'(v)` otherwise. **Expected:**
`residual with M'(v)=0 : 0`, `d(residual)/dM'(v) : nonzero`, `PASS`. This is why
the vacuum claim is scoped to constant `M`.

## 3. `proofs/invariant_distinction.py` — physicality

Computes `(nabla r)^2 = g^{rr}` on the solution and shows `d/dq (nabla r)^2 !=
0`, so the `q(v)` variation changes a scalar invariant and is not a coordinate
transformation. **Expected:** `PASS`.

## Reproducibility notes

- Pin versions with `requirements.txt`; record `pip freeze` alongside any run
  you want to cite.
- Certificates use exact `sympy.Rational` arithmetic; there is no floating-point
  tolerance to tune. The one numerical step is the `theta`-independence sanity
  check, which has a generous relative tolerance and is not part of the exact
  proof.
- To harden further: increase the number of sample points in
  `src/rbh/certificate.py::DEFAULT_POINTS`, or replace point-evaluation with a
  full symbolic `simplify` to zero (slower).
