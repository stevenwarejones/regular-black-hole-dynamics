"""Reduced action and Euler-Lagrange operator for the two-vector model.

Action (G = c = 1, signature -,+,+,+):

    S = (1/16 pi) \\int d^4x sqrt(-g) [ R + ell^2 ( L[A] - L[B] ) ],
    L[W] = 4 G^{mu nu} W_mu W_nu + 8 W^2 nabla_mu W^mu + 6 (W^2)^2.

The action is that of Eichhorn & Fernandes, arXiv:2508.00686v2, Eqs. (3)-(4);
see docs/RELATIONSHIP_TO_PRIOR_WORK.md and docs/SOURCES.md.

Spherically symmetric ansatz with EIGHT independent fields --- the areal
function Rf(v, r) is kept as an independent field and is NOT gauge-fixed before
variation, so the angular equation is obtained honestly rather than assumed:

    ds^2 = g00 dv^2 + 2 g01 dv dr + g11 dr^2 + Rf^2 dOmega^2,
    A = A0 dv + A1 dr,   B = B0 dv + B1 dr.

This module builds the reduced Lagrangian density on a jet space and provides
the per-field Euler-Lagrange operator. It performs NO simplification of the
huge intermediate expressions; certificates evaluate the residuals at rational
points (see rbh.certificate), which is exact and fast.
"""
from __future__ import annotations
import sympy as sp

FIELDS = ("g00", "g01", "g11", "Rf", "A0", "A1", "B0", "B1")
_ORD = 4  # jet order carried (enough for a 2nd-order Lagrangian's EL operator)

theta = sp.Symbol("theta")
ell, M = sp.symbols("ell M", positive=True)
v, r = sp.symbols("v r")

# jet coordinates: (field, d/dv order, d/dr order)
jet: dict[tuple[str, int, int], sp.Symbol] = {}
for _f in FIELDS:
    for _i in range(_ORD + 1):
        for _j in range(_ORD + 1 - _i):
            jet[(_f, _i, _j)] = sp.Symbol(f"{_f}_{_i}_{_j}")


def field(name: str) -> sp.Symbol:
    """The undifferentiated jet symbol for a field."""
    return jet[(name, 0, 0)]


def total_derivative(expr, coord: int):
    """Total derivative of a jet expression w.r.t. v (coord=0) or r (coord=1)."""
    expr = sp.sympify(expr)
    out = sp.Integer(0)
    fs = expr.free_symbols
    for (f, i, j), s in jet.items():
        if s in fs:
            ni, nj = (i + 1, j) if coord == 0 else (i, j + 1)
            out += jet[(f, ni, nj)] * sp.diff(expr, s)
    return out


def _coord_derivative(expr, mu: int):
    if mu == 0:
        return total_derivative(expr, 0)
    if mu == 1:
        return total_derivative(expr, 1)
    if mu == 2:  # theta
        return sp.diff(sp.sympify(expr), theta)
    return sp.Integer(0)  # phi: nothing depends on it


def _build_reduced_lagrangian():
    S = field
    g = sp.Matrix([
        [S("g00"), S("g01"), 0, 0],
        [S("g01"), S("g11"), 0, 0],
        [0, 0, S("Rf") ** 2, 0],
        [0, 0, 0, S("Rf") ** 2 * sp.sin(theta) ** 2],
    ])
    ginv = g.inv()
    detgam = S("g00") * S("g11") - S("g01") ** 2
    sqrtmg = sp.sqrt(-detgam) * S("Rf") ** 2 * sp.sin(theta)

    n = 4
    Gam = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                s = 0
                for e in range(n):
                    if ginv[a, e] == 0:
                        continue
                    s += ginv[a, e] * (_coord_derivative(g[e, b], c)
                                       + _coord_derivative(g[e, c], b)
                                       - _coord_derivative(g[b, c], e))
                Gam[a][b][c] = Gam[a][c][b] = s / 2

    Ric = sp.zeros(n, n)
    for b in range(n):
        for c in range(b, n):
            s = 0
            for a in range(n):
                s += _coord_derivative(Gam[a][b][c], a) - _coord_derivative(Gam[a][b][a], c)
                for e in range(n):
                    s += Gam[a][a][e] * Gam[e][b][c] - Gam[a][c][e] * Gam[e][b][a]
            Ric[b, c] = Ric[c, b] = s

    Rscalar = sum(ginv[a, b] * Ric[a, b] for a in range(n) for b in range(n))

    Gup = sp.zeros(n, n)
    for a in range(n):
        for b in range(n):
            s = 0
            for c in range(n):
                for e in range(n):
                    if ginv[a, c] == 0 or ginv[b, e] == 0:
                        continue
                    s += ginv[a, c] * ginv[b, e] * Ric[c, e]
            Gup[a, b] = s - ginv[a, b] * Rscalar / 2

    def Lvec(W0, W1):
        W = [W0, W1, sp.Integer(0), sp.Integer(0)]
        W2 = sum(ginv[a, b] * W[a] * W[b] for a in range(2) for b in range(2))
        Wup = [sum(ginv[a, b] * W[b] for b in range(2)) for a in range(2)]
        div = sum(_coord_derivative(sqrtmg * Wup[a], a) for a in range(2)) / sqrtmg
        kinetic = 4 * sum(Gup[a, b] * W[a] * W[b] for a in range(2) for b in range(2))
        return kinetic + 8 * W2 * div + 6 * W2 ** 2

    Lfull = sqrtmg * (Rscalar + ell ** 2 * (Lvec(S("A0"), S("A1")) - Lvec(S("B0"), S("B1")))) / sp.sin(theta)
    # After removing the common sin(theta) the density is theta-independent
    # (spherical symmetry); certificate.check_theta_independence() verifies this.
    return Lfull


_LFULL = None
_LRED = None


def reduced_lagrangian_full():
    """Cached reduced Lagrangian density with the theta dependence retained
    (used only by the theta-independence sanity check)."""
    global _LFULL
    if _LFULL is None:
        _LFULL = _build_reduced_lagrangian()
    return _LFULL


def reduced_lagrangian():
    """Cached reduced Lagrangian density evaluated at theta = pi/2."""
    global _LRED
    if _LRED is None:
        _LRED = reduced_lagrangian_full().subs(theta, sp.pi / 2)
    return _LRED


def euler_lagrange(name: str):
    """Euler-Lagrange expression E[field] for a second-order Lagrangian:

        E = L_phi - D_v L_{phi_v} - D_r L_{phi_r}
              + D_v^2 L_{phi_vv} + D_r^2 L_{phi_rr} + D_v D_r L_{phi_vr}.
    """
    L = reduced_lagrangian()
    D = total_derivative
    E = sp.diff(L, jet[(name, 0, 0)])
    E -= D(sp.diff(L, jet[(name, 1, 0)]), 0)
    E -= D(sp.diff(L, jet[(name, 0, 1)]), 1)
    E += D(D(sp.diff(L, jet[(name, 2, 0)]), 0), 0)
    E += D(D(sp.diff(L, jet[(name, 0, 2)]), 1), 1)
    E += D(D(sp.diff(L, jet[(name, 1, 1)]), 0), 1)
    return E
