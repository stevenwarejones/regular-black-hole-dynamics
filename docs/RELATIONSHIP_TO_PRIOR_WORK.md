# Relationship to prior work

## Origin of the model

The action and the static background specialized here are from:

> A. Eichhorn and P. G. S. Fernandes, *Regular black holes without
> mass-inflation instability and gravastars from modified gravity*,
> arXiv:2508.00686v2; Phys. Rev. D 113, L081501 (2026).

That paper introduces the two-vector action `S = (1/16 pi) \int sqrt(-g)[R +
ell^2(L[A] - L[B])]` and the static regular black hole used as the starting
point. This repository does not reproduce or re-derive their static analysis; it
takes their action and static branch as given and studies one time-dependent
extension.

## What is new here

Promoting the constant regulator `q` to a function `q(v)` and certifying that,
for constant `M`, the resulting time-dependent fields solve every equation of
motion with no added matter. A bounded literature search did not find this exact
time-dependent family in the inspected sources; this is a search finding, not a
priority certificate.

## Neutral note on the perturbative sector

The time-dependent family lies in the sector addressed by the source paper's
supplemental perturbative analysis. A detailed comparison between this exact
family and that analysis is **left to future work** and is not attempted here.
This repository makes **no claim** that any result of the source paper is
incorrect; it reports an exact solution and its certificates, and nothing more.

## Related lineage (context, not competition)

The vector–tensor construction and related primary hair appear in Charmousis,
Fernandes & Hassaine (arXiv:2504.13084). Dynamical regular black holes,
dynamical regularization scales, and vector contributions to gravitational
charges have close precedents in the literature; those broader ideas are **not**
claimed as novel here. This repository is scoped to the single exact solution
and its certificates.
