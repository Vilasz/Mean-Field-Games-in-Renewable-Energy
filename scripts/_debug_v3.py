"""Reproduz o erro do clearing_locacional do v3."""
import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import cvxpy as cp

ROOT = Path(r"c:/Users/joaof/Mean-Field-Games-in-Renewable-Energy")
sys.path.insert(0, str(ROOT))
from validate_model.pipeline import SINPaths, canonical_subsys

# Mínimo necessário
CACHE = ROOT / "outputs" / "cache"
panel = pd.read_parquet(CACHE / "panel_hourly_v3.parquet")
DATE_MIN = pd.Timestamp("2025-01-01 00:00:00")
DATE_MAX = pd.Timestamp("2025-10-31 23:00:00")
panel = panel[(panel["din_instante"] >= DATE_MIN) & (panel["din_instante"] <= DATE_MAX)].copy()

# Reproduz construção de variáveis chave
panel = panel.sort_values(["id_subsistema", "din_instante"]).reset_index(drop=True)
panel["G_s_total"] = panel["solar"].fillna(0.0)
panel["G_w_total"] = panel["wind"].fillna(0.0)
panel["G_h_total"] = panel["hydro"].fillna(0.0)
panel["G_th_total"] = panel["thermal"].fillna(0.0)
panel["G_nuc_inflex"] = panel["nuclear"].fillna(0.0)
panel["L_global"] = panel["D"]
panel["hora"] = pd.to_datetime(panel["din_instante"]).dt.hour

SUBSYS_ORDER = ["N", "NE", "SE", "S"]

typ_loc = (panel.groupby(["id_subsistema", "hora"])
                .agg(L=("L_global", "mean"),
                     Gs=("G_s_total", "mean"),
                     Gw=("G_w_total", "mean"),
                     Gh=("G_h_total", "mean"),
                     Gth=("G_th_total", "mean"),
                     Gnuc=("G_nuc_inflex", "mean"))
                .reset_index())

K_eff = panel.groupby("id_subsistema").agg(
    K_s=("G_s_total", lambda x: float(np.nanquantile(x, 0.99))),
    K_w=("G_w_total", lambda x: float(np.nanquantile(x, 0.99))),
    K_h=("G_h_total", lambda x: float(np.nanquantile(x, 0.99))),
    K_th=("G_th_total", lambda x: float(np.nanquantile(x, 0.99))),
).reindex(SUBSYS_ORDER)

HOURS = 24
def _arr(sub, col):
    return (typ_loc.query("id_subsistema == @sub").sort_values("hora")[col]
                   .reindex(range(HOURS)).fillna(method="ffill").fillna(method="bfill").values)

L_lt    = {s: _arr(s, "L")  for s in SUBSYS_ORDER}
Gw_lt   = {s: _arr(s, "Gw") for s in SUBSYS_ORDER}
Gnuc_lt = {s: _arr(s, "Gnuc") for s in SUBSYS_ORDER}
Havail_lt = {s: _arr(s, "Gh") for s in SUBSYS_ORDER}
a_s_loc = {s: (_arr(s, "Gs") / max(K_eff.loc[s, "K_s"], 1e-6)).clip(0, 1.5)
           for s in SUBSYS_ORDER}
K_TH_loc = {s: float(K_eff.loc[s, "K_th"]) for s in SUBSYS_ORDER}

# Lambda obs por subsistema-hora
import sys
print("L_lt N shape:", L_lt["N"].shape, "type:", type(L_lt["N"]))
print("a_s_loc N shape:", a_s_loc["N"].shape)

# Construir clearing simplificado
ALL_LINKS = [("N", "NE"), ("N", "SE"), ("NE", "SE"), ("SE", "S")]
F_MAX = {k: 5000.0 for k in ALL_LINKS}

# c local
c1_loc = {s: 100.0 for s in SUBSYS_ORDER}
c2_loc = {s: 0.005 for s in SUBSYS_ORDER}

# m0 inicial uniforme
K_S_OBS_TOTAL = float(K_eff["K_s"].sum())
G_s_lt = {ell: float(K_eff.loc[ell, "K_s"]) * a_s_loc[ell] for ell in SUBSYS_ORDER}

# Tenta clearing
locs = SUBSYS_ORDER
n = {l: cp.Variable(HOURS, nonneg=True) for l in locs}
h = {l: cp.Variable(HOURS, nonneg=True) for l in locs}
curt = {l: cp.Variable(HOURS, nonneg=True) for l in locs}
defc = {l: cp.Variable(HOURS, nonneg=True) for l in locs}
F = {(o, d): cp.Variable(HOURS) for (o, d) in ALL_LINKS}

cost = 0
for l in locs:
    cost += (c1_loc[l] * cp.sum(n[l]) + 0.5 * c2_loc[l] * cp.sum_squares(n[l])
             + 80.0 * cp.sum(h[l])
             + 3500.0 * cp.sum(defc[l]) + 30.0 * cp.sum(curt[l]))

bal = {}
cons = []
for l in locs:
    inflow  = sum(F[(o, d)] for (o, d) in ALL_LINKS if d == l)
    outflow = sum(F[(o, d)] for (o, d) in ALL_LINKS if o == l)
    bal[l] = (G_s_lt[l] + Gw_lt[l] + Gnuc_lt[l] + h[l] + n[l] + inflow - outflow) \
             == (L_lt[l] + curt[l] - defc[l])
    cons.append(bal[l])
    cons.append(n[l] <= K_TH_loc[l])
    cons.append(h[l] <= Havail_lt[l])
    cons.append(cp.sum(h[l]) <= Havail_lt[l].sum() * 0.95)
for (o, d) in ALL_LINKS:
    cons.append(F[(o, d)] <=  F_MAX[(o, d)])
    cons.append(F[(o, d)] >= -F_MAX[(o, d)])

prob = cp.Problem(cp.Minimize(cost), cons)
try:
    prob.solve(solver=cp.CLARABEL)
    print("status:", prob.status, "cost:", prob.value)
except Exception as e:
    print("ERROR:", type(e).__name__, e)
    import traceback
    traceback.print_exc()
