"""regular-black-hole-dynamics: exact time-dependent spherical solutions of a
two-vector model, with machine-checkable certificates.

See docs/ for model, conventions, scope, and provenance. This package is a
self-contained verifier; it derives the reduced Euler-Lagrange operator from
the action and checks the proposed solution against it.
"""
__all__ = ["model", "solution", "certificate"]
