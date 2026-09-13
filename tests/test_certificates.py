"""pytest wrappers around the certificates, for CI."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rbh import certificate


def test_theta_independence():
    assert certificate.check_theta_independence(verbose=False)


def test_exact_solution_all_fields_zero():
    ok, results = certificate.verify_exact_solution(verbose=False)
    assert ok, {k: v for k, v in results.items() if any(x != 0 for x in v)}


def test_grad_r_squared_depends_on_q():
    import sympy as sp
    from rbh import solution
    q0 = sp.Symbol("q0", positive=True)
    inv = solution.grad_r_squared(qfunc=q0)
    assert sp.simplify(sp.diff(inv, q0)) != 0
