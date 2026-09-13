# Provenance and status

- This is an **AI-produced research artifact under adversarial cross-review**.
  It has **not been reviewed by a human domain expert**.
- **No historical priority is claimed** for any result. "Novel" statements mean
  only that a bounded literature search did not surface a duplicate; that is not
  a priority certificate.
- The mathematics is checked by machine: the equations of motion are rebuilt
  from the action, and the solution is verified against them in **exact rational
  arithmetic** (no numerical tolerance). The certificates are the evidence; run
  them (`python run_checks.py`) rather than trusting this document.
- The action and static background originate in prior published work
  (`RELATIONSHIP_TO_PRIOR_WORK.md`, `SOURCES.md`). This repository does not
  independently re-derive that prior work.
- Scope limits are collected in `ASSUMPTIONS_AND_SCOPE.md`. The repository makes
  no stability, instability, singularity-resolution, quantum, or observational
  claim, and does not assert that any prior result is incorrect.

## How to attack this repository (invited)

The point of releasing certificates rather than prose is that disagreement can
be made concrete. Reviewers are invited to:

1. Re-derive the Euler–Lagrange operator independently (e.g. a 4D covariant
   computation, or a different symbolic pipeline) and compare on the solution.
2. Add sample points, or replace the rational point-evaluation with a full
   symbolic simplification to zero.
3. Vary the assumptions in `ASSUMPTIONS_AND_SCOPE.md` (signed perturbations,
   `F <= 0`, the excluded loci) and report where the statements stop holding.
4. Challenge the model transcription in `src/rbh/model.py` against the source
   action.
