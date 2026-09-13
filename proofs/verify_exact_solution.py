#!/usr/bin/env python3
"""PRIMARY CERTIFICATE (symbolic identity).

Proves that the time-dependent family (constant M, arbitrary smooth q(v)) solves
every reduced Euler-Lagrange equation, by showing each residual is IDENTICALLY
ZERO as a polynomial in independent regulator jets Q0, Q1, ... and the symbolic
parameters r, M, ell. This is a proof for all smooth q on the stated patch, not
an inference from samples.

The equations of motion are rebuilt from the action in src/rbh/model.py; nothing
is imported from an external derivation.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rbh import certificate


def main() -> int:
    print("[verify_exact_solution] symbolic Euler-Lagrange identities")
    print("Eight independent fields including the ungauged areal radius Rf.")
    print("q(v) and its derivatives are INDEPENDENT symbols (no fixed profile).\n")

    print("numerical theta-independence sanity check:")
    if not certificate.check_theta_independence():
        print("FAIL: reduced density retains theta dependence.")
        return 1

    print("\nsymbolic residuals on the solution (constant M, arbitrary q(v)):")
    ok, results = certificate.verify_vacuum_symbolic()

    print("\nsecondary regression check (not a proof):")
    samples_ok = certificate.sample_regression()
    ok = ok and samples_ok

    print()
    if ok:
        print("PASS: all eight reduced Euler-Lagrange residuals are IDENTICALLY "
              "ZERO for arbitrary smooth q(v) on the stated patch.")
        return 0
    print("FAIL:")
    if not samples_ok:
        print("  secondary sample regression did not vanish (see above).")
    for fn, num in results.items():
        if num != 0:
            print(f"  residual E[{fn}] did not reduce to zero: {num}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
