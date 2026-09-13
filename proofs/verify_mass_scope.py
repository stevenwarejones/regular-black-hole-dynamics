#!/usr/bin/env python3
"""SCOPE CERTIFICATE (symbolic identity): the vacuum statement is constant-M only.

With M -> M(v) and independent mass jets P0, P1, ... (and independent regulator
jets Q0, Q1, ...), we prove:

  * E[g11] - (-2 * M'(v)) is IDENTICALLY ZERO   (the full source term, derived
    not assumed --- the whole residual equals -2 P1, with no higher-derivative
    or regulator dependence);
  * every OTHER field equation is IDENTICALLY ZERO even with M'(v) != 0.

So promoting M to M(v) sources exactly the g_rr equation, proportional to M'(v)
-- the reduced-density image of the Vaidya-type stress T_vv proportional to
M'(v). Varying g_vv gives an upper-index equation (G^{vv}, and g^{vv}=0 in EF
coordinates), which is why the source lands in the g_rr variation.

The coefficient -2 is in the reduced-density normalization (the common angular /
action factor is suppressed; see docs/MODEL_AND_CONVENTIONS.md). We do NOT assert
the physical 1/(4 pi) normalization of T_vv here.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rbh import certificate


def main() -> int:
    print("[verify_mass_scope] M -> M(v) with independent mass and regulator jets\n")
    ok = certificate.verify_mass_source_symbolic()
    print()
    if ok:
        print("PASS: M -> M(v) sources exactly E[g11] = -2 M'(v) (full expression, "
              "no higher-derivative terms); all other equations stay vacuum. The "
              "no-matter claim is correctly scoped to constant M.")
        return 0
    print("FAIL: mass-source identity not established as stated.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
