# regular-black-hole-dynamics

Exact **time-dependent** spherical solutions of a two-vector model of regular
black holes, with machine-checkable certificates in exact rational arithmetic.

> **Status.** AI-produced research artifact under adversarial cross-review. Not
> reviewed by a human domain expert. No historical priority claimed. The
> certificates are the evidence — run them. See
> [`docs/PROVENANCE.md`](docs/PROVENANCE.md).

## The result

For the action (Eichhorn & Fernandes, arXiv:2508.00686v2, Eqs. 3–4; `G=c=1`,
signature `-+++`)

$$
S=\frac1{16\pi}\int\!\sqrt{-g}\,\big[R+\ell^2(\mathcal L[A]-\mathcal L[B])\big],
\quad
\mathcal L[W]=4G^{\mu\nu}W_\mu W_\nu+8W^2\nabla_\mu W^\mu+6(W^2)^2,
$$

take the `q_A=q`, `q_B=-q` static branch and **promote the constant regulator to
a function `q(v)`**. In ingoing Eddington–Finkelstein form,

$$
D=r^3+2\ell^2 q(v),\quad f=1-\frac{2Mr^2}{D},\quad
ds^2=-f\,dv^2+2\,dv\,dr+r^2 d\Omega^2,
$$
$$
A=\Big(\tfrac{1-f}{2r}-\tfrac{q(v)}{4r^2}\Big)dv,\qquad
B=\Big(\tfrac{1-f}{2r}+\tfrac{q(v)}{4r^2}\Big)dv .
$$

**Claim (certified).** For **constant `M`** and arbitrary smooth `q(v)`, these
fields solve *all eight* reduced Euler–Lagrange equations on `r>0` patches, with
**no added matter** — including the angular equation, because the areal radius
is kept as an independent field and is not gauge-fixed before variation.

The equations of motion are **rebuilt from the action** in
[`src/rbh/model.py`](src/rbh/model.py); nothing is imported from an external
derivation. The certificate is a **symbolic identity**, not a sample check:
each residual is derived, the solution substituted, and `q(v)` together with all
its `v`-derivatives replaced by **independent** symbols `Q0, Q1, …`; the reduced
numerator is then required to be the zero polynomial in those jets and the
symbolic parameters `r, M, ℓ`. That is a proof for *all* smooth `q(v)` on the
patch. (Fixed-profile rational sample evaluations are retained only as a labeled
secondary regression tripwire — they do not establish an identity.)

## What this is and is not

Scope, assumptions, and the (deliberately neutral) relationship to the source
paper's perturbative analysis are in
[`docs/ASSUMPTIONS_AND_SCOPE.md`](docs/ASSUMPTIONS_AND_SCOPE.md) and
[`docs/RELATIONSHIP_TO_PRIOR_WORK.md`](docs/RELATIONSHIP_TO_PRIOR_WORK.md). In
short: this repository reports **one exact solution and its certificates**. It
makes no stability, instability, singularity-resolution, quantum, or
observational claim, and does not assert that any published result is incorrect.
The "vacuum" statement is specific to constant `M` (a time-dependent `M(v)`
sources matter).

## Run it

```bash
pip install -r requirements.txt
python run_checks.py        # full verification pipeline (~1 minute)
pytest -q                   # same certificates, as tests
```

See [`VERIFY.md`](VERIFY.md) for exactly what each certificate checks and the
expected output.

## Layout

```
src/rbh/        model (action -> EL operator), solution, exact-arithmetic certificate
proofs/         runnable certificates (exact solution, mass scope, invariant distinction)
tests/          pytest wrappers for CI
docs/           model & conventions, scope, provenance, prior-work relationship, sources
run_checks.py   single-command pipeline (failures propagate, nonzero exit)
```

## Citing / license

MIT (see [`LICENSE`](LICENSE)). Please read [`CITATION.cff`](CITATION.cff) and
`docs/PROVENANCE.md` before citing — particularly the no-priority and
not-peer-reviewed statements.
