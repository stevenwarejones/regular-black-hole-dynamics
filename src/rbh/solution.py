"""The proposed exact time-dependent spherical family.

In ingoing Eddington-Finkelstein form with lapse N = 1 and c = d = 0:

    D(v, r) = r^3 + 2 ell^2 q(v),
    f(v, r) = 1 - 2 M r^2 / D,
    A0 = (1 - f) / (2 r) - q(v) / (4 r^2),
    B0 = (1 - f) / (2 r) + q(v) / (4 r^2),

    ds^2 = -f dv^2 + 2 dv dr + r^2 dOmega^2,   A = A0 dv,   B = B0 dv.

M is the mass parameter; q(v) is the (time-dependent) regularization function.
This is the q_A = q, q_B = -q branch of the static background of
arXiv:2508.00686v2, with the constant regulator promoted to a function q(v).

SCOPE (see docs/ASSUMPTIONS_AND_SCOPE.md):
  * The vacuum (no-added-matter) statement is for CONSTANT M. A time-dependent
    M(v) sources T_vv = M'(v) / (4 pi r^2); see proofs/verify_mass_scope.py.
  * Statements are local on r > 0 patches with a0 b0 != 0; q(v) > 0, ell != 0.
"""
from __future__ import annotations
import sympy as sp
from . import model
from .model import ell, M, v, r, FIELDS

# q as an abstract smooth function of v (constant M assumed for the vacuum claim)
q = sp.Function("q", positive=True)(v)


def field_expressions(mass=M, qfunc=q):
    """Explicit (v, r) expressions for the eight ansatz fields."""
    D = r ** 3 + 2 * ell ** 2 * qfunc
    f = 1 - 2 * mass * r ** 2 / D
    a0 = (1 - f) / (2 * r) - qfunc / (4 * r ** 2)
    b0 = (1 - f) / (2 * r) + qfunc / (4 * r ** 2)
    return {
        "g00": -f, "g01": sp.Integer(1), "g11": sp.Integer(0), "Rf": r,
        "A0": a0, "A1": sp.Integer(0), "B0": b0, "B1": sp.Integer(0),
    }


def solution_jet(mass=M, qfunc=q):
    """Map every jet symbol to the corresponding derivative of the solution."""
    exprs = field_expressions(mass=mass, qfunc=qfunc)
    sj = {}
    for fn, expr in exprs.items():
        for (name, i, j), sym in model.jet.items():
            if name == fn:
                sj[sym] = sp.diff(expr, v, i, r, j)
    return sj


def grad_r_squared(mass=M, qfunc=q):
    """The invariant (nabla r)^2 = g^{rr} on the solution.

    For the EF metric with N = 1 this equals f = 1 - 2 M r^2 / D. See
    proofs/invariant_distinction.py for the static-background gauge argument that
    uses it.
    """
    exprs = field_expressions(mass=mass, qfunc=qfunc)
    g = sp.Matrix([
        [exprs["g00"], exprs["g01"]],
        [exprs["g01"], exprs["g11"]],
    ])
    return sp.simplify(g.inv()[1, 1])


def delta_grad_r_squared_static():
    """First variation of the invariant (nabla r)^2 about a STATIC background.

    Perturb q = q0 + eps * p (delta Rf = 0) about constant q0. Because
    (nabla r)^2 depends on v only through q, its first variation is

        delta I = d/deps (nabla r)^2 |_{eps=0}
                = 4 M ell^2 r^2 / (r^3 + 2 ell^2 q0)^2 * p .

    Returns (delta_I, p, q0). Compare against the displayed expression in
    proofs/invariant_distinction.py. A spherical pure-gauge perturbation
    preserving delta Rf = 0 has xi^r = 0 and cannot change this static scalar,
    so a nonzero delta_I (for M > 0) certifies the perturbation is not pure gauge.
    """
    eps, p, q0 = sp.symbols("eps p q0", real=True)
    I = grad_r_squared(qfunc=q0 + eps * p)
    dI = sp.simplify(sp.diff(I, eps).subs(eps, 0))
    return dI, p, q0
