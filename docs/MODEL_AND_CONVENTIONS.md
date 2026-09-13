# Model and conventions

Units: `G = c = 1`. Signature `(-, +, +, +)`.

## Action

$$
S=\frac1{16\pi}\int d^4x\sqrt{-g}\,\big[R+\ell^2(\mathcal L[A]-\mathcal L[B])\big],
\qquad
\mathcal L[W]=4G^{\mu\nu}W_\mu W_\nu+8W^2\nabla_\mu W^\mu+6(W^2)^2 .
$$

Here `A` and `B` are two covector fields, `G^{mu nu}` is the Einstein tensor,
and `W^2 = g^{mu nu} W_mu W_nu`. There is **no Maxwell term**. The covariant
vector components are varied independently of the metric. This action is that of
Eichhorn & Fernandes, arXiv:2508.00686v2, Eqs. (3)–(4); see
`RELATIONSHIP_TO_PRIOR_WORK.md` and `SOURCES.md`.

## Spherical ansatz (eight independent fields)

$$
ds^2 = g_{00}\,dv^2 + 2g_{01}\,dv\,dr + g_{11}\,dr^2 + R_f^2\,d\Omega^2,
\qquad A = A_0\,dv + A_1\,dr,\quad B = B_0\,dv + B_1\,dr .
$$

The areal function `Rf(v, r)` is retained as an **independent field** and is not
gauge-fixed before variation. This matters: fixing the areal radius to `r` a
priori removes a scalar and forces the angular equation to be recovered
indirectly (by a Noether identity). Keeping `Rf` free means the angular
(`Rf`) equation is checked on the same footing as the others. The verifier
therefore certifies **eight** Euler–Lagrange equations, `g00, g01, g11, Rf, A0,
A1, B0, B1`.

## Reduced Lagrangian

`src/rbh/model.py` builds the Christoffel symbols, Ricci tensor, scalar
curvature, and Einstein tensor for the ansatz directly, assembles the density
`sqrt(-g)[R + ell^2(L[A] - L[B])]`, and removes the common `sin(theta)` factor.
The resulting density is `theta`-independent (spherical symmetry); this is
checked numerically at random Lorentzian field values in
`certificate.check_theta_independence()`. The per-field Euler–Lagrange operator
for this second-order Lagrangian is

$$
E[\phi]=L_\phi-D_v L_{\phi_v}-D_r L_{\phi_r}
        +D_v^2 L_{\phi_{vv}}+D_r^2 L_{\phi_{rr}}+D_vD_r L_{\phi_{vr}} .
$$

## Static background specialized here

With `N = 1`, `A_1 = B_1 = 0`,

$$
D=r^3+2\ell^2 q,\quad f=1-\frac{2Mr^2}{D},\quad
A_0=\frac{1-f}{2r}-\frac{q}{4r^2},\quad
B_0=\frac{1-f}{2r}+\frac{q}{4r^2}.
$$

This is the `q_A = q`, `q_B = -q` branch. The contribution of this repository is
to promote the constant `q` to a function `q(v)` and certify that, for constant
`M`, the resulting time-dependent fields still solve every equation of motion.
