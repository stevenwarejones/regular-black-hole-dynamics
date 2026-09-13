"""pytest wrappers and NEGATIVE regression tests for the certificates.

The negative tests are the point of the review response: they confirm the
identity checker actually rejects residuals that vanish only on the old fixed
sample profile (e.g. q'' - 1/5, q''', M'^2, M''), rather than silently
specializing them to zero.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import sympy as sp
from rbh import certificate
from rbh.certificate import Q, P


def test_theta_independence():
    assert certificate.check_theta_independence(verbose=False)


def test_vacuum_identities_all_fields():
    ok, results = certificate.verify_vacuum_symbolic(verbose=False)
    assert ok, {k: v for k, v in results.items() if v != 0}


def test_mass_source_identity():
    assert certificate.verify_mass_source_symbolic(verbose=False)


def test_invariant_first_variation():
    from rbh import solution
    from rbh.model import ell, M, r
    dI, p, q0 = solution.delta_grad_r_squared_static()
    expected = 4 * M * ell ** 2 * r ** 2 / (r ** 3 + 2 * ell ** 2 * q0) ** 2 * p
    assert sp.simplify(dI - expected) == 0


# --- negative tests: the checker must REJECT these, not zero them out ---

def test_checker_rejects_qpp_minus_const():
    # q'' - 1/5 vanishes for the old sample profile but is not identically zero.
    assert not certificate.is_identically_zero(Q[2] - sp.Rational(1, 5))


def test_checker_rejects_qppp():
    # q''' vanishes for the old sample profile but is not identically zero.
    assert not certificate.is_identically_zero(Q[3])


def test_checker_rejects_nonlinear_mass_term():
    # A quadratic mass term is not a linear source; must be rejected.
    assert not certificate.is_identically_zero(P[1] ** 2)


def test_checker_rejects_higher_mass_derivative():
    assert not certificate.is_identically_zero(P[2])


def test_checker_accepts_true_zero():
    assert certificate.is_identically_zero(sp.Integer(0))


def _load_verify_exact_solution():
    import importlib.util
    path = os.path.join(os.path.dirname(__file__), "..", "proofs", "verify_exact_solution.py")
    spec = importlib.util.spec_from_file_location("verify_exact_solution", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_pipeline_fails_if_sample_regression_fails(monkeypatch):
    # The primary script must incorporate the secondary regression result: if the
    # sample check fails, main() must exit nonzero (reproduces the review finding).
    mod = _load_verify_exact_solution()
    monkeypatch.setattr(mod.certificate, "check_theta_independence", lambda *a, **k: True)
    monkeypatch.setattr(mod.certificate, "verify_vacuum_symbolic", lambda *a, **k: (True, {}))
    monkeypatch.setattr(mod.certificate, "sample_regression", lambda *a, **k: False)
    assert mod.main() != 0


def test_pipeline_passes_when_all_subchecks_pass(monkeypatch):
    mod = _load_verify_exact_solution()
    monkeypatch.setattr(mod.certificate, "check_theta_independence", lambda *a, **k: True)
    monkeypatch.setattr(mod.certificate, "verify_vacuum_symbolic", lambda *a, **k: (True, {}))
    monkeypatch.setattr(mod.certificate, "sample_regression", lambda *a, **k: True)
    assert mod.main() == 0
