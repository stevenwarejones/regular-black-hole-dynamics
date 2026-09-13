#!/usr/bin/env python3
"""PRIMARY CERTIFICATE.

Verifies that the proposed time-dependent family (constant M, free q(v)) solves
every Euler-Lagrange equation of the reduced two-vector action, in exact
rational arithmetic, at several sample points. Exits nonzero on any failure.

The Euler-Lagrange operator is derived from the action in src/rbh/model.py; the
solution is defined in src/rbh/solution.py. Nothing here is imported from any
external result --- the equations are rebuilt from the action.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rbh import certificate


def main() -> int:
    print("[verify_exact_solution] reduced-action Euler-Lagrange residuals")
    print("Eight independent fields including the ungauged areal radius Rf.\n")

    print("theta-independence sanity check:")
    if not certificate.check_theta_independence():
        print("FAIL: reduced density retains theta dependence.")
        return 1

    print("\nexact-rational residuals on the solution (constant M):")
    ok, results = certificate.verify_exact_solution()

    print()
    if ok:
        print("PASS: all eight Euler-Lagrange residuals are EXACTLY ZERO on the "
              "solution at every rational sample point.")
        return 0
    print("FAIL: a residual did not vanish. Offending fields:")
    for fn, vals in results.items():
        if any(x != 0 for x in vals):
            print(f"  {fn}: {vals}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
