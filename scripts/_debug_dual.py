"""Test dual extraction: numerical vs analytical."""
from __future__ import annotations

import numpy as np
import cvxpy as cp


def build_solve(A_perturb_idx: int = -1, perturb: float = 0.0, V_max: float = 1.0e8):
    rng = np.random.default_rng(0)
    T = 24 * 30
    t = np.arange(T)
    D = 50_000 + 5_000 * np.sin(2 * np.pi * t / 24) + 500 * rng.standard_normal(T)
    g_ren = np.maximum(0, 12_000 * np.sin(2 * np.pi * (t / 24) - np.pi / 3))
    A = 25_000 + 5_000 * np.sin(2 * np.pi * t / (24 * 7))
    if A_perturb_idx >= 0:
        A[A_perturb_idx] += perturb

    c1, c2 = 50.0, 1e-3
    Gh_max, Gt_max = 80_000.0, 60_000.0
    inv_V = 1.0 / V_max

    gh = cp.Variable(T, nonneg=True)
    c = cp.Variable(T, nonneg=True)
    d = cp.Variable(T, nonneg=True)
    spill = cp.Variable(T, nonneg=True)
    Vn = cp.Variable(T + 1)
    gt = D - gh - g_ren + c - d
    obj = c1 * cp.sum(gt) + 0.5 * c2 * cp.sum_squares(gt) + 1e4 * cp.sum(d) + 10 * cp.sum(c) + 1 * cp.sum(spill)
    dyn = Vn[1:] == Vn[:-1] + (A - gh - spill) * inv_V
    cons = [
        gh <= Gh_max,
        gt >= 0, gt <= Gt_max,
        Vn >= 0, Vn <= 1.0,
        Vn[0] == 0.50, Vn[-1] >= 0.20, Vn[-1] <= 0.95,
        dyn,
    ]
    prob = cp.Problem(cp.Minimize(obj), cons)
    prob.solve(solver=cp.CLARABEL, verbose=False)
    return prob.value, dyn.dual_value, gh.value, gt.value


print("Test 1: V_max = 1e7 (apertado)")
v0, lam0, gh0, gt0 = build_solve(V_max=1.0e7)
print(f"  cost={v0:.2f}")
print(f"  lam stats: min={np.min(np.abs(lam0)):.4f}  max={np.max(np.abs(lam0)):.4f}  mean={np.mean(np.abs(lam0)):.4f}")
print(f"  lam * inv_V = water_value: mean={np.mean(np.abs(lam0)) / 1.0e7:.6f}")
print(f"  gt mean: {np.mean(gt0):.2f}, expected water_value c1+c2*gt: {50 + 1e-3*np.mean(gt0):.4f}")

# Sanity: numerical perturbation at index 100
v_pert, _, _, _ = build_solve(A_perturb_idx=100, perturb=1.0, V_max=1.0e7)
print(f"  numeric water_value at idx 100: {v0 - v_pert:.4f}")
print(f"  dual at idx 100 * inv_V: {lam0[100] / 1.0e7:.4f}  (raw lam: {lam0[100]:.4f})")

print()
print("Test 2: V_max = 1e8 (folgado)")
v1, lam1, gh1, gt1 = build_solve(V_max=1.0e8)
print(f"  cost={v1:.2f}")
print(f"  lam stats: min={np.min(np.abs(lam1)):.4f}  max={np.max(np.abs(lam1)):.4f}  mean={np.mean(np.abs(lam1)):.4f}")
print(f"  water_value = lam / V_max: mean={np.mean(np.abs(lam1)) / 1.0e8:.6f}")
v_pert, _, _, _ = build_solve(A_perturb_idx=100, perturb=1.0, V_max=1.0e8)
print(f"  numeric water_value at idx 100: {v1 - v_pert:.4f}")
print(f"  dual at idx 100: {lam1[100]:.4f}, lam/V_max: {lam1[100] / 1.0e8:.6f}")
