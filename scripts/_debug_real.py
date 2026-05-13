"""Inspect water_value on real SIN data, with diagnostics."""
from __future__ import annotations

import numpy as np
import pandas as pd

from validate_model.pipeline import SINPaths, build_panel, load_hidrologia
from validate_model.dispatch_models import AnnualParams, solve_annual_dispatch


paths = SINPaths(root="validate_model", year=2025)
panel = build_panel(paths)

sin = (
    panel.groupby("din_instante", as_index=False)
    .agg({"D": "sum", "gs": "sum", "gr": "sum", "gh": "sum", "gn_obs": "sum"})
    .sort_values("din_instante")
    .reset_index(drop=True)
)
sin["g_ren"] = sin["gs"] + sin["gr"]
sin = sin.dropna(subset=["D", "gh", "g_ren"]).reset_index(drop=True)

hydro_daily = load_hidrologia(paths.hidro_di_path)
hydro_daily["date"] = pd.to_datetime(hydro_daily["din_instante"]).dt.normalize()
sin["date"] = sin["din_instante"].dt.normalize()
sin = sin.merge(hydro_daily[["date", "A_MW"]], on="date", how="left")
sin["A_MW"] = sin["A_MW"].interpolate(method="linear").bfill().ffill()

Gh_max = float(np.quantile(sin["gh"].dropna(), 0.995))
Gt_max = float(np.quantile(sin["gn_obs"].dropna(), 0.999))
alpha = float(hydro_daily["vol_util_pct"].iloc[0] / 100.0)

params = AnnualParams(
    Gh_max=Gh_max, Gt_max=Gt_max, V_max=1.0e8,
    alpha=alpha, terminal_mode="band", beta_min=0.20, beta_max=0.95,
    c1=50.0, c2=1e-3,
    pi_d=1e4, pi_c=10.0, pi_s=1.0,
)
import cvxpy as cp

D = sin["D"].to_numpy(dtype=float)
g_ren = sin["g_ren"].to_numpy(dtype=float)
A = sin["A_MW"].to_numpy(dtype=float)
T = len(D)
gh = cp.Variable(T); c_var = cp.Variable(T); d_var = cp.Variable(T); spill = cp.Variable(T)
s = cp.Variable(T+1)
gt = D - gh - g_ren + c_var - d_var
cost = params.c1*cp.sum(gt) + 0.5*params.c2*cp.sum_squares(gt) + params.pi_d*cp.sum(d_var) + params.pi_c*cp.sum(c_var) + params.pi_s*cp.sum(spill)
V = cp.Variable(T+1)
V_max = params.V_max
V_0 = params.alpha * params.V_max
dyn = V[1:] == V[:-1] + (A - gh - spill)
cons = [
    gh >= 0, gh <= params.Gh_max,
    gt >= 0, gt <= params.Gt_max,
    c_var >= 0, d_var >= 0, spill >= 0,
    V >= 0, V <= V_max,
    V[0] == V_0,
    dyn,
    V[-1] >= params.beta_min*V_max, V[-1] <= params.beta_max*V_max,
]
s = V  # rename for compatibility
prob = cp.Problem(cp.Minimize(cost), cons)
prob.solve(solver=cp.CLARABEL, verbose=False)
print(f"status={prob.status} cost={prob.value:.2f}")

dval = dyn.dual_value
print("dyn dual_value type:", type(dval))
print("dyn dual_value shape:", None if dval is None else np.asarray(dval).shape)
print("dyn dual_value abs stats: min=", np.min(np.abs(dval)), "max=", np.max(np.abs(dval)), "mean=", np.mean(np.abs(dval)))
print("dyn dual_value head:", np.asarray(dval).ravel()[:5])
print("dyn dual_value tail:", np.asarray(dval).ravel()[-5:])

# For comparison, with sample interior gh point, expected lambda = scale * (c1 + c2*gt)
gh_v = gh.value; gt_v = (D - gh_v - g_ren + c_var.value - d_var.value)
scale = 1.0
print("expected (c1 + c2*gt) at t=100:", (params.c1 + params.c2*gt_v[100]))
print("actual dyn dual at t=100:", np.asarray(dval).ravel()[100] if dval is not None else None)

res = {"gh": gh_v, "gt": gt_v, "c": c_var.value, "d": d_var.value, "spill": spill.value, "V": s.value, "water_value": np.asarray(dval).ravel() if dval is not None else None, "total_cost": prob.value, "status": prob.status}
print(f"status={res['status']}, cost={res['total_cost']:.2f}, T={len(sin)}")

gh = res["gh"]; gt = res["gt"]; c = res["c"]; d = res["d"]; spill = res["spill"]
wv = res["water_value"]

print(f"gh:  mean={gh.mean():.0f}  Gh_max={Gh_max:.0f}  hours_at_max={(gh > Gh_max - 1).sum()}  hours_at_zero={(gh < 1).sum()}")
print(f"gt:  mean={gt.mean():.2f}  hours_zero={(gt < 1).sum()}")
print(f"c:   mean={c.mean():.0f}  hours_pos={(c > 1).sum()}")
print(f"d:   mean={d.mean():.0f}  hours_pos={(d > 1).sum()}")
print(f"V:   min={res['V'].min()/params.V_max*100:.1f}%  max={res['V'].max()/params.V_max*100:.1f}%  hours_at_floor={(res['V'] < 1).sum()}  hours_at_ceil={(res['V'] > params.V_max - 1).sum()}")
print(f"\nwater_value: mean={np.nanmean(wv):.4f}  min={np.nanmin(wv):.6f}  max={np.nanmax(wv):.4f}")
print(f"  hours nonzero: {(wv > 1e-3).sum()} / {len(wv)}")
print(f"  expected (c1+c2*gt_mean): {params.c1 + params.c2*gt.mean():.4f}")

# Decomposition: are most hours at Gh_max binding or gh=0 binding?
hist_gh = pd.Series(gh).describe()
print("\n=== gh series describe ===")
print(hist_gh)
