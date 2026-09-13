#!/usr/bin/env python3
"""GAUGE CERTIFICATE: the q perturbation about a STATIC background is not pure gauge.

The general slogan "a coordinate change cannot alter a scalar invariant" is not a
sufficient perturbative gauge argument: a scalar's coordinate representation
changes under a diffeomorphism by its Lie derivative, and for an already
time-dependent solution a constant time translation gives a variation
proportional to q'(v) while preserving the EF ansatz. We therefore make the
precise, valid statement about a STATIC background.

Setup: background with constant q0, perturbation q = q0 + eps p(v), delta Rf = 0.
For a spherical diffeomorphism xi = xi^v d_v + xi^r d_r, preserving delta Rf = 0
about Rf = r forces xi^r = 0. The static scalar I0 = (nabla r)^2 = f0(r) is
independent of advanced time, so a pure-gauge perturbation obeying that condition
has delta_xi I0 = xi^r f0'(r) = 0.

But the proposed perturbation gives (this script checks the first variation
symbolically)

    delta I = 4 M ell^2 r^2 / (r^3 + 2 ell^2 q0)^2 * p(v),

nonzero for M > 0 wherever p(v) != 0 on the nondegenerate positive-mass patch.
Hence this perturbation about the static background is NOT pure gauge.

Caveats kept explicit: M > 0 is required for the nonvanishing coefficient (the
metric-invariant argument alone does not settle M = 0); a non-gauge perturbation
is not automatically a propagating healthy mode or an instability.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import sympy as sp
from rbh import solution
from rbh.model import ell, M, r


def main() -> int:
    print("[invariant_distinction] first variation of (nabla r)^2 about a static background\n")
    dI, p, q0 = solution.delta_grad_r_squared_static()
    expected = 4 * M * ell ** 2 * r ** 2 / (r ** 3 + 2 * ell ** 2 * q0) ** 2 * p
    matches = sp.simplify(dI - expected) == 0
    print(f"  computed delta I : {dI}")
    print(f"  expected         : {expected}")
    print(f"  match            : {matches}")

    # nonvanishing on the positive-mass nondegenerate patch (M > 0, p != 0)
    sample = {M: sp.Rational(1), r: sp.Rational(3), ell: sp.Rational(1), q0: sp.Rational(3), p: sp.Rational(1)}
    val = sp.simplify(dI.subs(sample))
    nonzero = val != 0
    print(f"  value at (M,r,ell,q0,p)=(1,3,1,3,1): {val}  (nonzero: {nonzero})")

    print()
    if matches and nonzero:
        print("PASS: the first variation matches the displayed expression and is "
              "nonzero for M > 0; the perturbation about the static background is "
              "not pure gauge (delta Rf = 0 forces xi^r = 0, so a gauge mode would "
              "give zero).")
        return 0
    print("FAIL: gauge certificate not established.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
