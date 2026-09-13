#!/usr/bin/env python3
"""INVARIANT CERTIFICATE: the q(v) variation is physical, not pure gauge.

Computes the scalar invariant (nabla r)^2 = g^{mu nu} d_mu r d_nu r on the
solution and shows its dependence on q is nonzero. A pure coordinate
transformation cannot change a scalar invariant, so a nonzero d/dq of (nabla r)^2
certifies that varying q is not a residual coordinate freedom.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import sympy as sp
from rbh import solution
from rbh.model import ell, M, r, v


def main() -> int:
    print("[invariant_distinction] (nabla r)^2 on the solution\n")
    qsym = sp.Symbol("q0", positive=True)      # treat q as a parameter here
    inv = solution.grad_r_squared(qfunc=qsym)  # = f = 1 - 2 M r^2 / (r^3 + 2 ell^2 q0)
    inv = sp.simplify(inv)
    dq = sp.simplify(sp.diff(inv, qsym))
    print(f"  (nabla r)^2 = {inv}")
    print(f"  d/dq (nabla r)^2 = {dq}")

    # nonzero for generic M, r, ell, q0
    sample = {M: sp.Rational(1), r: sp.Rational(3), ell: sp.Rational(1), qsym: sp.Rational(3)}
    val = sp.simplify(dq.subs(sample))
    print(f"  value at (M,r,ell,q0)=(1,3,1,3): {val}")

    ok = (dq != 0) and (val != 0)
    print()
    if ok:
        print("PASS: (nabla r)^2 depends on q, so the q(v) variation changes a "
              "scalar invariant and is not a pure coordinate transformation.")
        return 0
    print("FAIL: invariant is q-independent.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
