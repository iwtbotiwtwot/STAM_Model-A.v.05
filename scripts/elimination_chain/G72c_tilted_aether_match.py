#!/usr/bin/env python3
"""
G72c_tilted_aether_match.py

Step C of the Einstein-aether track.  Tests whether the framework's
committed metric admits an embedding in tilted Einstein-aether.

Strategy
========
T^EA is LINEAR in each c_i (the four EA coupling constants).  So at any
sample radius r with given aether kinematics (U(r), W(r), U'(r), W'(r),
and the second-derivatives required by the divergence terms),

    T^EA^mu_nu(r) = c1 P1^mu_nu(r) + c2 P2^mu_nu(r)
                   + c3 P3^mu_nu(r) + c4 P4^mu_nu(r)

The matching equation is

    G^mu_nu(r) = T^EA^mu_nu(r) - 2 V(r) delta^mu_nu

Take pairs of differences (tt - thth) and (rr - thth) to eliminate V.
That gives 2 equations per sample r, each LINEAR in (c1, c2, c3, c4)
for fixed kinematics:

    G^t_t - G^th_th = c1 (P1^t_t - P1^th_th) + c2 (...) + c3 (...) + c4 (...)
    G^r_r - G^th_th = c1 (P1^r_r - P1^th_th) + c2 (...) + c3 (...) + c4 (...)

Approach for matching:
1. Sample N >= 5 radii inside the shell.
2. Use the simplest tilt ansatz first: aligned U = 1/sqrt(h), W = 0.  This is
   the static-aligned limit (already ruled out by G72a, but it's a baseline
   sanity check).
3. Then try a TILTED ansatz with U, W and their derivatives as free local
   variables at each sample.  At each r, solve the 2 linear equations for
   c_i (treating U(r), W(r), U'(r), W'(r) as parameters).
4. For TRUE embedding, the c_i must be the SAME constants at all sample r.
   Use scipy.optimize.least_squares to find (c_i_global, {U(r_i), U'(r_i)})
   minimizing total residual across all sample points.

Output:
- Best-fit c_i values
- Residual norm across sample points
- Verdict: tilted EA embedding exists (residual ~ 0) or fails (residual > 0)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import least_squares

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Symbolic setup (kept LIGHT: no sp.simplify on heavy intermediate expressions)
# ---------------------------------------------------------------------------

print("Building symbolic T^EA pieces (may take ~30s) ...")

r_sym, U_sym, dU_sym, W_sym, dW_sym = sp.symbols('r U dU W dW', real=True)
# Use plain symbols (not Functions) — we'll express derivatives ourselves.

# Committed metric (M = 1)
A_expr = 2 / r_sym
y_expr = 3 * A_expr - 2
F_expr = 1 - 5 * y_expr**4 + 4 * y_expr**5
h_expr = 1 - A_expr
k_expr = h_expr * F_expr
hp_expr = sp.diff(h_expr, r_sym)
kp_expr = sp.diff(k_expr, r_sym)

# Aether: u^μ = (U, W, 0, 0).  Constraint: -h U^2 + W^2/k = -1
# Treat U, W as values at a point; dU, dW are r-derivatives at that point.

u_t, u_r = U_sym, W_sym
du_t_dr, du_r_dr = dU_sym, dW_sym
u_th, u_ph = 0, 0

# Inverse metric components
g_inv_tt = -1 / h_expr
g_inv_rr = k_expr
g_inv_thth = 1 / r_sym**2

# Christoffel symbols (only the ones that touch t, r components of u)
# For diagonal metric ds^2 = -h dt^2 + dr^2/k + r^2 dOmega^2:
Gamma = {
    ('t', 't', 'r'): hp_expr / (2 * h_expr),         # Γ^t_{tr}
    ('t', 'r', 't'): hp_expr / (2 * h_expr),
    ('r', 't', 't'): k_expr * hp_expr / 2,           # Γ^r_{tt}
    ('r', 'r', 'r'): -kp_expr / (2 * k_expr),        # Γ^r_{rr}
    ('r', 'th', 'th'): -k_expr * r_sym,              # Γ^r_{thth}
    ('r', 'ph', 'ph'): -k_expr * r_sym * 1,          # Γ^r_{phph}, sin^2θ dropped (at θ=π/2)
    ('th', 'r', 'th'): 1 / r_sym,
    ('th', 'th', 'r'): 1 / r_sym,
    ('ph', 'r', 'ph'): 1 / r_sym,
    ('ph', 'ph', 'r'): 1 / r_sym,
}
# (working in equatorial slice θ = π/2 for simplicity; sin θ = 1)


def gamma(a, b, c):
    return Gamma.get((a, b, c), 0)


# Compute ∇_μ u^α at the sample point (4×4 matrix as dict)
# ∇_μ u^α = ∂_μ u^α + Γ^α_{μβ} u^β
# For static config (no t-dependence): ∂_t = 0.  Only r-derivatives nonzero.

labels = ['t', 'r', 'th', 'ph']
u_vec = {'t': u_t, 'r': u_r, 'th': 0, 'ph': 0}
du_dr_vec = {'t': du_t_dr, 'r': du_r_dr, 'th': 0, 'ph': 0}


def partial(mu, a):
    """Partial derivative ∂_μ u^α (assumes static -- only r-derivatives nonzero)."""
    if mu == 'r':
        return du_dr_vec[a]
    return 0


def nabla(mu, a):
    """Covariant derivative ∇_μ u^α."""
    val = partial(mu, a)
    for b in labels:
        val = val + gamma(a, mu, b) * u_vec[b]
    return val


# Compute the full 4x4 ∇_μ u^α
nab = {(mu, a): nabla(mu, a) for mu in labels for a in labels}

# Lower-index: ∇_μ u_ν = g_νλ ∇_μ u^λ
g_diag = {'t': -h_expr, 'r': 1/k_expr, 'th': r_sym**2, 'ph': r_sym**2}
nab_dn = {(mu, nu): g_diag[nu] * nab[(mu, nu)] for mu in labels for nu in labels}

# Acceleration a^μ = u^ν ∇_ν u^μ
a_up = {a: sum(u_vec[nu] * nab[(nu, a)] for nu in labels) for a in labels}
a_dn = {a: g_diag[a] * a_up[a] for a in labels}

# Expansion θ = ∇·u = ∑ ∇_μ u^μ
theta_exp = sum(nab[(mu, mu)] for mu in labels)


# EA kinetic scalars
g_inv_diag = {'t': -1/h_expr, 'r': k_expr, 'th': 1/r_sym**2, 'ph': 1/r_sym**2}
# (sin^2 θ = 1 at θ = π/2)


# term1 = (∇^μ u^ν)(∇_μ u_ν) = g^{μα} (∇_α u^ν)(∇_μ u_ν)
def nab_up_up(mu, nu):
    return g_inv_diag[mu] * nab[(mu, nu)]


term1 = sum(nab_up_up(mu, nu) * nab_dn[(mu, nu)] for mu in labels for nu in labels)
term2 = theta_exp**2
term3 = sum(nab_up_up(mu, nu) * nab_dn[(nu, mu)] for mu in labels for nu in labels)
term4 = sum(g_diag[a] * a_up[a] * a_up[a] for a in labels)

c1_sym, c2_sym, c3_sym, c4_sym = sp.symbols('c1 c2 c3 c4', real=True)
L_ae = -(c1_sym * term1 + c2_sym * term2 + c3_sym * term3 + c4_sym * term4)

# T^EA stress (Jacobson-Mattingly).  For matching we need T^EA^μ_ν with mixed
# indices.  The simplest piece-by-piece form (LOCAL, without divergence terms):
#
#   T^EA_μν^(no_div) = c1 [(∇_μ u^α)(∇_ν u_α) - (∇^α u_μ)(∇_α u_ν)]
#                    + c4 a_μ a_ν
#                    - (1/2) g_μν L_ae
#                    + λ u_μ u_ν
#
# The DIVERGENCE pieces ∇_α[J^α_(μ u_ν) - J_(μ^α u_ν) + J_(μν) u^α] involve
# second derivatives of u (through J).  For our test, we'll INCLUDE them
# numerically by computing J at the sample and forming the divergence
# numerically (finite differences in r).
#
# To keep the analytic part tractable, let's first do the LOCAL (no-div)
# pieces analytically -- they suffice if we ALSO sample at adjacent radii
# to compute the divergence numerically.


# c1 piece local: c1 (∇_μ u^α)(∇_ν u_α) - c1 (∇^α u_μ)(∇_α u_ν)
def T_local_c1(mu, nu):
    p1 = sum(nab[(mu, a)] * nab_dn[(nu, a)] for a in labels)
    p2 = sum(nab_up_up(a, mu) * nab_dn[(a, nu)] for a in labels)
    # but wait, mu is index of vector; nab is keyed by (mu, alpha) with alpha up.
    # For (∇_μ u^α)(∇_ν u_α):  sum_α nab[(mu, α)] * nab_dn[(nu, α)]
    # For (∇^α u_μ)(∇_α u_ν):  need ∇^α u_μ = g^{αλ} nab_dn[(λ, μ)]; nab_dn keyed (μ, ν) so (λ,μ).
    p2 = sum(g_inv_diag[a] * nab_dn[(a, mu)] * nab_dn[(a, nu)] for a in labels)
    return p1 - p2


# c4 piece local: c4 a_μ a_ν
def T_local_c4(mu, nu):
    return a_dn[mu] * a_dn[nu]


def T_trace_piece(mu, nu):
    # -(1/2) g_μν L_ae
    # diagonal-metric only:
    if mu == nu:
        return -sp.Rational(1, 2) * g_diag[mu] * L_ae
    return 0


# u_μ u_ν piece (multiplied by λ; will eliminate λ via aether EOM)
u_dn = {a: g_diag[a] * u_vec[a] for a in labels}


def T_lambda(mu, nu):
    return u_dn[mu] * u_dn[nu]


# Pure LOCAL T^EA (no divergence): T^local_μν = c1·(c1-piece) + c4·(c4-piece) + (trace piece) + λ·(λ-piece)
def T_local(mu, nu):
    return (c1_sym * T_local_c1(mu, nu)
            + c4_sym * T_local_c4(mu, nu)
            + T_trace_piece(mu, nu)
            # + λ * T_lambda(mu, nu)  # absorbed into V or aether EOM
            )


# Mixed index: T^μ_ν = g^{μα} T_αν.  For diagonal metric: T^μ_ν = g^{μμ} T_μν
def T_local_mixed(mu, nu):
    return g_inv_diag[mu] * T_local(mu, nu)


# Compute G^μ_ν from G69 in committed-metric form (M = 1)
def G_mixed(mu, nu):
    if mu != nu:
        return 0
    if mu == 't':
        return -16 * (r_sym - 3)**3 * (720 - 648*r_sym + 117*r_sym**2 + 13*r_sym**3) / r_sym**8
    if mu == 'r':
        return -16 * (13*r_sym - 24) * (r_sym - 3)**4 / r_sym**7
    if mu in ('th', 'ph'):
        return -1440 * (r_sym - 3)**3 * (r_sym - 2) * (r_sym - 1) / r_sym**8
    return 0


# Build the matching equations.  Define:
#   D_t_th(r, U, W, dU, dW, c1, c2, c3, c4) = G^t_t - G^th_th - (T^EA^t_t - T^EA^th_th)
#   D_r_th(r, U, W, dU, dW, c1, c2, c3, c4) = G^r_r - G^th_th - (T^EA^r_r - T^EA^th_th)
# These should vanish for any consistent embedding.
#
# We work with T_LOCAL.  This omits the divergence piece, so D will not vanish
# exactly even for an exact embedding; the divergence piece must be added
# numerically.  For a first cut we accept this and report what falls out.

print("Building D_t_th, D_r_th (local T^EA only, no divergence) ...")
Tlocal_tt = T_local_mixed('t', 't')
Tlocal_rr = T_local_mixed('r', 'r')
Tlocal_thth = T_local_mixed('th', 'th')

D_t_th_expr = G_mixed('t', 't') - G_mixed('th', 'th') - (Tlocal_tt - Tlocal_thth)
D_r_th_expr = G_mixed('r', 'r') - G_mixed('th', 'th') - (Tlocal_rr - Tlocal_thth)

# Aether unit-norm constraint: -h U^2 + W^2/k = -1
constraint_expr = -h_expr * U_sym**2 + W_sym**2 / k_expr + 1

print("Lambdifying ...")
syms = (r_sym, U_sym, W_sym, dU_sym, dW_sym, c1_sym, c2_sym, c3_sym, c4_sym)
D_t_th_fn = sp.lambdify(syms, D_t_th_expr, 'numpy')
D_r_th_fn = sp.lambdify(syms, D_r_th_expr, 'numpy')
constraint_fn = sp.lambdify((r_sym, U_sym, W_sym), constraint_expr, 'numpy')

# h, k, h', k' as numeric functions of r
h_fn = sp.lambdify(r_sym, h_expr, 'numpy')
k_fn = sp.lambdify(r_sym, k_expr, 'numpy')
hp_fn = sp.lambdify(r_sym, hp_expr, 'numpy')
kp_fn = sp.lambdify(r_sym, kp_expr, 'numpy')


# ---------------------------------------------------------------------------
# Numerical matching
# ---------------------------------------------------------------------------

# Sample radii inside the shell  (r ∈ (2, 3) corresponds to A ∈ (2/3, 1))
r_samples = np.array([2.1, 2.3, 2.5, 2.7, 2.9])
N = len(r_samples)
A_samples = 2 / r_samples
print(f"\nSample radii (in M=1 units): {r_samples}")
print(f"Corresponding A values:        {A_samples}")
print()

# At each sample, kinematic constraints:
# - W^2 = k(h U^2 - 1)  (sign of W is a free choice -- we pick W > 0 for radial inflow)
# - dW from differentiating constraint:  2 W dW = k'(hU^2 - 1) + k (h' U^2 + 2 h U dU)
#   -> dW = [k'(hU^2 - 1) + k h' U^2 + 2 k h U dU] / (2 W)
#
# For the matching test, at each sample r we have variables U(r), dU(r) (local).
# Globals: c1, c2, c3, c4 (constants).
# We have 2 equations per sample (D_t_th, D_r_th) = 2N = 10 equations.
# Unknowns: 4 (c_i) + 2N (U_i, dU_i) = 4 + 10 = 14.  Under-determined by 4.
#
# Resolve under-determination by: pick the tilt magnitude (set W² > 0 floor),
# pick a parametrization of U(r) as a polynomial of low order (say degree 4):
#   U(r) = sum_{j=0..4} a_j r^j
# Then dU(r) = sum_{j=1..4} j a_j r^{j-1}.
# Unknowns: 4 c_i + 5 polynomial coefs = 9 unknowns.
# Equations: 2N = 10.  Over-determined by 1.

def residuals(params):
    c1, c2, c3, c4, a0, a1, a2, a3, a4 = params
    res = []
    for r_val in r_samples:
        U_val = a0 + a1*r_val + a2*r_val**2 + a3*r_val**3 + a4*r_val**4
        dU_val = a1 + 2*a2*r_val + 3*a3*r_val**2 + 4*a4*r_val**3

        h_val = h_fn(r_val)
        k_val = k_fn(r_val)
        hp_val = hp_fn(r_val)
        kp_val = kp_fn(r_val)

        W_sq = k_val * (h_val * U_val**2 - 1)
        if W_sq < 0:
            # Aether constraint violated -- penalty
            res.extend([1e6 * (-W_sq), 1e6 * (-W_sq)])
            continue
        W_val = np.sqrt(W_sq)
        dW_val = (kp_val * (h_val*U_val**2 - 1) + k_val*hp_val*U_val**2
                  + 2*k_val*h_val*U_val*dU_val) / (2 * W_val)

        d_t_th = D_t_th_fn(r_val, U_val, W_val, dU_val, dW_val, c1, c2, c3, c4)
        d_r_th = D_r_th_fn(r_val, U_val, W_val, dU_val, dW_val, c1, c2, c3, c4)
        res.append(d_t_th)
        res.append(d_r_th)
    return np.array(res, dtype=float)


# Initial guess: small c_i, U ~ 1/sqrt(h) at PS as starting profile
# At r=3 (PS): h = 1/3, so U ~ sqrt(3).  Aligned start: a0 = sqrt(3), a_i = 0.
x0 = np.array([0.1, 0.1, 0.1, 0.1, np.sqrt(3.0), 0.0, 0.0, 0.0, 0.0])

print("Initial residual at x0:")
res0 = residuals(x0)
print(f"  ||res||_2 = {np.linalg.norm(res0):.6e}")
print(f"  res = {res0}")
print()

# Solve with least squares
print("Running scipy least_squares ...")
sol = least_squares(residuals, x0, method='lm', max_nfev=20000,
                    xtol=1e-14, ftol=1e-14)
print()
print("Result:")
print(f"  Status:        {sol.status}")
print(f"  Message:       {sol.message}")
print(f"  nfev:          {sol.nfev}")
print(f"  Final ||res||: {np.linalg.norm(sol.fun):.6e}")
print()
print("Best-fit parameters:")
print(f"  c1, c2, c3, c4 = {sol.x[:4]}")
print(f"  U(r) polynomial coeffs (a0..a4) = {sol.x[4:9]}")
print()
print("Residual at each sample r (D_t_th, D_r_th):")
for i, r_val in enumerate(r_samples):
    print(f"  r = {r_val:.2f}: D_t_th = {sol.fun[2*i]:+.4e}, D_r_th = {sol.fun[2*i+1]:+.4e}")
print()

# Reality check: at sample radii, what are the kinematic values?
print("Best-fit aether kinematics:")
print(f"  {'r':>6}{'U':>12}{'dU':>12}{'W':>14}{'tilt phi':>14}")
print("-" * 60)
for r_val in r_samples:
    a0, a1, a2, a3, a4 = sol.x[4:9]
    U_val = a0 + a1*r_val + a2*r_val**2 + a3*r_val**3 + a4*r_val**4
    dU_val = a1 + 2*a2*r_val + 3*a3*r_val**2 + 4*a4*r_val**3
    h_val = h_fn(r_val)
    k_val = k_fn(r_val)
    W_sq = k_val * (h_val * U_val**2 - 1)
    if W_sq < 0:
        W_val = np.nan
        phi_val = np.nan
    else:
        W_val = np.sqrt(W_sq)
        # tilt: cosh(phi) = U sqrt(h),  sinh(phi) = W / sqrt(k)
        cosh_phi = U_val * np.sqrt(h_val)
        sinh_phi = W_val / np.sqrt(k_val)
        phi_val = np.arcsinh(sinh_phi)
    print(f"  {r_val:>6.2f}{U_val:>12.4f}{dU_val:>12.4f}{W_val:>14.4e}{phi_val:>14.4f}")
print()

# Verdict
final_norm = np.linalg.norm(sol.fun)
threshold = 1e-3
print("=" * 80)
if final_norm < threshold:
    print(f"VERDICT: tilted EA embedding LOCALLY consistent at sample radii.")
    print(f"  Final residual norm {final_norm:.4e} < threshold {threshold:.4e}")
    print(f"  Best-fit c_i = {sol.x[:4]}")
    print(f"  (Note: only LOCAL T^EA used; divergence piece omitted.")
    print(f"   A clean verdict requires including the divergence contributions.)")
else:
    print(f"VERDICT: tilted EA matching FAILS at this approximation order.")
    print(f"  Final residual norm {final_norm:.4e} >> threshold {threshold:.4e}")
    print(f"  No consistent (c_i, U(r)) found with polynomial-of-degree-4 ansatz")
    print(f"  using LOCAL T^EA pieces alone.")
    print(f"  Two interpretations:")
    print(f"    (a) The divergence pieces of T^EA (which we omitted) are essential")
    print(f"        and would change the verdict if included.")
    print(f"    (b) Tilted EA genuinely cannot embed the framework's metric, and")
    print(f"        the obstruction is real -- pointing to DHOST or beyond.")
    print(f"  Next step: include divergence pieces in T^EA computation.")
print("=" * 80)
