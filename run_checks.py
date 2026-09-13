#!/usr/bin/env python3
"""Single-command verification pipeline.

Runs every certificate as a subprocess and reports pass/fail. A nonzero exit
from any certificate fails the whole run (assertions are real; failures
propagate). This checks the SCIENCE, not just file hashes.

    python run_checks.py
"""
import subprocess, sys, os, time

HERE = os.path.dirname(os.path.abspath(__file__))
PROOFS = [
    ("exact solution (primary)", "proofs/verify_exact_solution.py"),
    ("mass scope (constant-M)",  "proofs/verify_mass_scope.py"),
    ("invariant distinction",    "proofs/invariant_distinction.py"),
]


def main() -> int:
    print("=" * 70)
    print("regular-black-hole-dynamics : verification pipeline")
    print("=" * 70)
    failures = []
    for label, path in PROOFS:
        print(f"\n>>> {label}  ({path})")
        print("-" * 70)
        t0 = time.time()
        rc = subprocess.call([sys.executable, os.path.join(HERE, path)])
        dt = time.time() - t0
        print(f"--- {label}: {'PASS' if rc == 0 else 'FAIL'}  ({dt:.1f}s)")
        if rc != 0:
            failures.append(label)

    print("\n" + "=" * 70)
    if failures:
        print(f"RESULT: FAIL ({len(failures)} of {len(PROOFS)}): {', '.join(failures)}")
        return 1
    print(f"RESULT: PASS (all {len(PROOFS)} certificates)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
