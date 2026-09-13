#!/usr/bin/env python3
"""SCOPE CERTIFICATE: the vacuum statement is constant-M only.

Promotes M -> M(v) and shows exactly one field equation stops being satisfied:
the variation with respect to g_rr (field 'g11'), whose residual becomes
proportional to M'(v). All seven other equations remain satisfied. This is the
reduced-action image of the Vaidya-type source T_vv proportional to M'(v):
varying g_vv gives an upper-index equation (G^{vv}, and g^{vv}=0 in EF
coordinates), so the source instead appears in the g_rr variation, which picks
up G_{vv}.

We certify the qualitative, correctly-located statement in exact arithmetic:
  * E[g11] = 0 when M'(v)=0 and is linear-and-nonzero in M'(v);
  * every other E[field] = 0 even for M'(v) != 0.
We do not assert the 1/(4 pi) normalization; the reduced action suppresses the
common angular factor (see docs/MODEL_AND_CONVENTIONS.md).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import sympy as sp
from rbh import model, solution
from rbh.model import ell, v, r, FIELDS


def _q_rat(order, vv):
    if order == 0:
        return sp.Rational(3) + sp.Rational(1, 10) * vv ** 2
    if order == 1:
        return sp.Rational(1, 5) * vv
    if order == 2:
        return sp.Rational(1, 5)
    return sp.Integer(0)


def _residual(field_name, v0, rval, ellval):
    M0, Mp = sp.symbols("M0 Mp", real=True)
    mass = M0 + Mp * (v - sp.Rational(v0))       # M(v0)=M0, M'(v0)=Mp
    sol_jet = solution.solution_jet(mass=mass)
    E = model.euler_lagrange(field_name).xreplace(sol_jet)
    qf = solution.q
    subs = {(sp.Derivative(qf, (v, k)) if k > 0 else qf): _q_rat(k, sp.Rational(v0))
            for k in range(7)}
    E = E.xreplace(subs).subs({v: sp.Rational(v0), r: sp.Rational(rval), ell: sp.Rational(ellval)})
    return sp.expand(sp.cancel(E)), Mp


def main() -> int:
    print("[verify_mass_scope] promoting M -> M(v) (local Taylor M0 + Mp*(v-v0))\n")
    v0, rval, ellval = 2, 3, 1
    ok = True

    # 1. the sourced equation
    e, Mp = _residual("g11", v0, rval, ellval)
    static = sp.simplify(e.subs(Mp, 0))
    slope = sp.simplify(sp.diff(e, Mp))
    print(f"  sourced equation E[g11] at (v0,r,ell)=({v0},{rval},{ellval}):")
    print(f"    residual with M'(v0)=0 : {static}")
    print(f"    d(residual)/dM'(v0)    : {slope}  ({'nonzero' if slope != 0 else 'zero'})")
    ok = ok and (static == 0) and (slope != 0)

    # 2. all other equations remain satisfied even for M'(v) != 0
    print("\n  other equations (must stay 0 for all M'):")
    for fn in FIELDS:
        if fn == "g11":
            continue
        e2, Mp2 = _residual(fn, v0, rval, ellval)
        vanishes = sp.simplify(e2) == 0
        print(f"    E[{fn:4s}] identically 0 in M' : {vanishes}")
        ok = ok and vanishes

    print()
    if ok:
        print("PASS: M -> M(v) sources exactly the g_rr equation (proportional to "
              "M'(v)); all others stay vacuum. The no-matter claim is correctly "
              "scoped to constant M.")
        return 0
    print("FAIL: unexpected mass-scope behavior.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
