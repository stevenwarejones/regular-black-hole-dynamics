"""Certificates for the exact time-dependent family.

The primary certificates are SYMBOLIC IDENTITIES, not sample checks. For each
field, the Euler-Lagrange expression is derived from the action, the proposed
solution is substituted, and q(v) (and, for the mass-scope certificate, M(v))
together with all their v-derivatives are replaced by INDEPENDENT symbols
Q0, Q1, ... (and P0, P1, ...). The reduced numerator must be identically zero
as a polynomial in the independent jets and the symbolic parameters r, M, ell.
That is a proof for arbitrary smooth q(v) on the stated patch --- not an
inference from finitely many samples.

Sample-point evaluations are retained only as clearly-labeled secondary
regression checks; they are NOT the proof and do not establish an identity.
"""
from __future__ import annotations
import sympy as sp
from . import model, solution
from .model import ell, M, v, r, FIELDS

# Independent jet symbols for q(v) and M(v). The reduced EL operator carries
# field jets to order 4, so q/M enter at most to that order; we allocate a few
# extra and assert below that none beyond what appears is needed.
QMAX = 6
Q = sp.symbols(f"Q0:{QMAX + 1}")
P = sp.symbols(f"P0:{QMAX + 1}")

ALLOWED_PARAMS = {r, M, ell, *Q, *P}


def _replace_by_independent_jets(expr, func, symbols):
    reps = {func: symbols[0]}
    for k in range(1, QMAX + 1):
        reps[sp.Derivative(func, (v, k))] = symbols[k]
    return expr.xreplace(reps)


def _reduced_numerator(expr):
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    return sp.expand(num), den


def _fully_resolved(expr) -> bool:
    """True if expr contains no unresolved model jets, solution functions, or
    derivatives --- only the allowed independent symbols and parameters."""
    if expr.atoms(sp.Derivative) or expr.atoms(sp.core.function.AppliedUndef):
        return False
    return expr.free_symbols <= ALLOWED_PARAMS


def is_identically_zero(expr) -> bool:
    """Reduce expr to a single fraction and test that its numerator is the zero
    polynomial in the independent symbols. Used by both the proof scripts and
    the negative regression tests. No tolerances, no sampling, no nsimplify."""
    num, _ = _reduced_numerator(expr)
    return num == 0 and _fully_resolved(num)


# ---------------------------------------------------------------------------
# Vacuum identity (constant M): all eight EL residuals vanish for arbitrary q.
# ---------------------------------------------------------------------------
def vacuum_residual(field):
    """Euler-Lagrange residual E[field] on the solution (constant M), with q and
    its v-derivatives replaced by independent symbols Q0, Q1, .... Returns the
    residual as a sympy expression in {r, M, ell, Q0, Q1, ...}."""
    sol_jet = solution.solution_jet()
    E = model.euler_lagrange(field).xreplace(sol_jet)
    return _replace_by_independent_jets(E, solution.q, Q)


def verify_vacuum_symbolic(verbose=True):
    ok, results = True, {}
    for fn in FIELDS:
        E = vacuum_residual(fn)
        num, _ = _reduced_numerator(E)
        resolved = _fully_resolved(num)
        jets = sorted(k for k in range(QMAX + 1) if E.has(Q[k]))
        field_ok = (num == 0) and resolved
        ok = ok and field_ok
        results[fn] = num
        if verbose:
            status = "IDENTICALLY ZERO" if field_ok else f"NONZERO: {num}"
            print(f"  E[{fn:4s}] (indep. q-jets {jets or 'none'}) -> {status}", flush=True)
    return ok, results


# ---------------------------------------------------------------------------
# Mass-scope identity: with M -> M(v), exactly E[g11] = -2 M'(v); others vanish.
# ---------------------------------------------------------------------------
def mass_residual(field):
    """E[field] on the solution with M -> M(v), q -> q(v), all v-derivatives of
    both replaced by independent symbols P0,P1,... (mass) and Q0,Q1,...
    (regulator). Returns the residual in {r, ell, P*, Q*}."""
    Mv = sp.Function("Mfun")(v)
    sol_jet = solution.solution_jet(mass=Mv)
    E = model.euler_lagrange(field).xreplace(sol_jet)
    E = _replace_by_independent_jets(E, Mv, P)
    E = _replace_by_independent_jets(E, solution.q, Q)
    return E


def verify_mass_source_symbolic(verbose=True):
    """Certify, with independent mass/regulator jets:
      * E[g11] - (-2 * P1) is identically zero  (the full source, not just slope);
      * E[g11] contains no mass jet P_k for k >= 2 (no higher-derivative source);
      * every other E[field] is identically zero even with M'(v) != 0.
    The coefficient -2 is in the reduced-density normalization (the common
    angular/action factor is suppressed; see docs/MODEL_AND_CONVENTIONS.md), so
    we do NOT assert the physical 1/(4 pi) of T_vv here."""
    ok = True
    Eg11 = mass_residual("g11")
    source = -2 * P[1]
    full_ok = is_identically_zero(Eg11 - source)
    higher = sorted(k for k in range(2, QMAX + 1) if Eg11.has(P[k]))
    q_dep = sorted(k for k in range(QMAX + 1) if Eg11.has(Q[k]))
    ok = ok and full_ok and not higher
    if verbose:
        print(f"  E[g11] - (-2*M'(v)) identically zero : {full_ok}")
        print(f"  higher mass jets in E[g11] (want none): {higher or 'none'}")
        print(f"  regulator jets in E[g11]              : {q_dep or 'none'}")
    for fn in FIELDS:
        if fn == "g11":
            continue
        z = is_identically_zero(mass_residual(fn))
        ok = ok and z
        if verbose:
            print(f"  E[{fn:4s}] identically zero (indep. M',q jets): {z}")
    return ok


# ---------------------------------------------------------------------------
# Secondary regression checks (NOT proofs): the solution at fixed rational
# points with a fixed rational q profile. Kept only to catch gross breakage.
# ---------------------------------------------------------------------------
def _q_rational(order, vv):
    if order == 0:
        return sp.Rational(3) + sp.Rational(1, 10) * vv ** 2
    if order == 1:
        return sp.Rational(1, 5) * vv
    if order == 2:
        return sp.Rational(1, 5)
    return sp.Integer(0)


DEFAULT_POINTS = [(1, 2, 1, 1), (2, 3, 2, 1), (-1, 5, 1, 2),
                  (3, 7, 5, 3), (0, 4, 1, 1), (5, 2, 3, 2), (4, 9, 2, 5)]


def sample_regression(verbose=True):
    """Secondary: evaluate each vacuum residual at rational points with the fixed
    profile q = 3 + v^2/10. A regression tripwire, explicitly NOT an identity
    proof (a residual can vanish on this profile without vanishing in general)."""
    sol_jet = solution.solution_jet()
    ok = True
    for fn in FIELDS:
        E = model.euler_lagrange(fn).xreplace(sol_jet)
        for (vv, rv, Mv, ev) in DEFAULT_POINTS:
            reps = {(sp.Derivative(solution.q, (v, k)) if k else solution.q):
                    _q_rational(k, sp.Rational(vv)) for k in range(QMAX + 1)}
            val = E.xreplace(reps).subs({v: sp.Rational(vv), r: sp.Rational(rv),
                                         M: sp.Rational(Mv), ell: sp.Rational(ev)})
            if sp.cancel(sp.together(val)) != 0:
                ok = False
    if verbose:
        print(f"  sample regression over {len(DEFAULT_POINTS)} points -> "
              f"{'ok' if ok else 'FAILED'}")
    return ok


def check_theta_independence(trials=6, seed=0, verbose=True):
    """Numerical sanity check (NOT part of the exact proof): the reduced density,
    with sin(theta) removed, does not vary with theta, as spherical symmetry
    requires. Uses random physically-valid (Lorentzian) metrics."""
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
        if max(abs(o - outs[0]) for o in outs) > 1e-7 * max(1.0, abs(outs[0])):
            ok = False
    if verbose:
        print(f"  theta-independence (numerical sanity) over {trials} states -> "
              f"{'OK' if ok else 'FAILED'}")
    return ok
