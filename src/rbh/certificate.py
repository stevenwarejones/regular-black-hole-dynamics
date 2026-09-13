"""Exact-arithmetic certificates.

The Euler-Lagrange residuals of the reduced action are rational functions of the
fields and their derivatives. Evaluating them on the solution at generic RATIONAL
sample points in exact arithmetic gives exactly 0 iff the residual is identically
0 on the solution (a nonzero rational function cannot vanish at generic rational
points). This is the same "exact certificate, no numerical tolerance" approach
used elsewhere in the author's repositories.
"""
from __future__ import annotations
import sympy as sp
from . import model, solution
from .model import ell, M, v, r, FIELDS

# A default rational q(v) used for point evaluation: q = 3 + v^2/10 (> 0 on the
# sampled domain). Only its value and derivatives at the sample v enter.
def _q_rational(order: int, vval: sp.Rational) -> sp.Rational:
    if order == 0:
        return sp.Rational(3) + sp.Rational(1, 10) * vval ** 2
    if order == 1:
        return sp.Rational(1, 5) * vval
    if order == 2:
        return sp.Rational(1, 5)
    return sp.Integer(0)


# Rational sample points (v, r, M, ell); r > 0, ell != 0, chosen to avoid the
# excluded loci (a0 b0 = 0, horizons) at generic values.
DEFAULT_POINTS = [
    (1, 2, 1, 1), (2, 3, 2, 1), (-1, 5, 1, 2),
    (3, 7, 5, 3), (0, 4, 1, 1), (5, 2, 3, 2), (4, 9, 2, 5),
]


def _eval_residual_at_point(E, sol_jet, vv, rv, Mv, ev, max_order=6):
    e = E.xreplace(sol_jet)  # inject solution jets (symbolic in v,r,M,ell,q-derivs)
    subs = {}
    qf = solution.q
    for k in range(max_order + 1):
        der = sp.Derivative(qf, (v, k)) if k > 0 else qf
        subs[der] = _q_rational(k, sp.Rational(vv))
    e = e.xreplace(subs)
    e = e.subs({v: sp.Rational(vv), r: sp.Rational(rv),
                M: sp.Rational(Mv), ell: sp.Rational(ev)})
    return sp.cancel(sp.together(e))


def verify_exact_solution(points=None, verbose=True):
    """Check every Euler-Lagrange residual vanishes on the solution (constant M).

    Returns (ok: bool, results: dict[field] -> list of exact residual values).
    """
    points = points or DEFAULT_POINTS
    sol_jet = solution.solution_jet()  # constant M
    results, ok = {}, True
    for fn in FIELDS:
        E = model.euler_lagrange(fn)
        vals = [_eval_residual_at_point(E, sol_jet, *pt) for pt in points]
        vals = [sp.nsimplify(x) if not x.free_symbols else sp.simplify(x) for x in vals]
        field_ok = all(x == 0 for x in vals)
        ok = ok and field_ok
        results[fn] = vals
        if verbose:
            status = "ALL ZERO" if field_ok else f"NONZERO: {vals}"
            print(f"  E[{fn:4s}] at {len(points)} rational points -> {status}", flush=True)
    return ok, results


def check_theta_independence(trials=6, seed=0, verbose=True):
    """Sanity check that the reduced density (with sin(theta) removed) carries no
    residual theta dependence, as spherical symmetry requires. Evaluates at
    random PHYSICALLY-VALID metrics (Lorentzian gamma: g00 < 0 so det gamma < 0)
    for several theta values and confirms the density does not vary with theta.
    """
    import numpy as np
    from .model import theta as theta_sym
    L = model.reduced_lagrangian_full()
    used = sorted([k for k in model.jet if model.jet[k] in L.free_symbols], key=str)
    fn = sp.lambdify([model.jet[k] for k in used] + [ell, M, theta_sym], L, "numpy")
    rng = np.random.default_rng(seed)
    ok = True
    for _ in range(trials):
        vals = []
        for k in used:
            if k == ("g00", 0, 0):
                vals.append(-rng.uniform(0.4, 1.3))
            elif k == ("g01", 0, 0):
                vals.append(rng.uniform(-0.2, 0.2))
            else:
                vals.append(rng.uniform(0.4, 1.3))
        vals += [1.0, 1.0]
        outs = [float(fn(*(vals + [t]))) for t in (0.7, 1.3, 2.2)]
        spread = max(abs(o - outs[0]) for o in outs)
        if spread > 1e-7 * max(1.0, abs(outs[0])):
            ok = False
    if verbose:
        print(f"  theta-independence over {trials} random Lorentzian states -> "
              f"{'OK' if ok else 'FAILED'}", flush=True)
    return ok
