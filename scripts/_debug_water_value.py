"""Inspect water_value (dual of reservoir dynamics) on a small instance."""
from __future__ import annotations

import numpy as np

from validate_model.dispatch_models import AnnualParams, solve_annual_dispatch


def run(label: str, params: AnnualParams, *, T: int = 24 * 14) -> None:
    rng = np.random.default_rng(0)
    t = np.arange(T)
    D = 50_000 + 5_000 * np.sin(2 * np.pi * t / 24) + 500 * rng.standard_normal(T)
    g_ren = np.maximum(0, 12_000 * np.sin(2 * np.pi * (t / 24) - np.pi / 3))
    A = 25_000 + 5_000 * np.sin(2 * np.pi * t / (24 * 7))

    res = solve_annual_dispatch(D=D, g_ren=g_ren, A=A, params=params)
    wv = res["water_value"]
    gt = res["gt"]
    print(f"=== {label} ===")
    print(f"  status={res['status']}  cost={res['total_cost']:.2f}")
    print(f"  gh: mean={np.nanmean(res['gh']):.0f}  max={np.nanmax(res['gh']):.0f}  Gh_max={params.Gh_max}")
    print(f"  gt: mean={np.nanmean(gt):.2f}  max={np.nanmax(gt):.2f}")
    print(f"  V_min={np.nanmin(res['V'])/params.V_max*100:.1f}%  V_max={np.nanmax(res['V'])/params.V_max*100:.1f}%")
    if wv is not None:
        print(f"  water_value: mean={np.nanmean(wv):.4f}  min={np.nanmin(wv):.4f}  max={np.nanmax(wv):.4f}")
        print(f"  expected ~ c1+c2*gt_mean = {params.c1 + params.c2*np.nanmean(gt):.4f}")
    else:
        print("  water_value: None")
    print()


if __name__ == "__main__":
    base = AnnualParams(
        Gh_max=80_000.0, Gt_max=60_000.0, V_max=2.0e7,
        alpha=0.50, terminal_mode="band", beta_min=0.20, beta_max=0.95,
        c1=50.0, c2=1e-3,
    )
    run("Folgado (V_max=20 TWh)", base)
    run("Apertado (Gh_max=20 GW)", base.with_overrides(Gh_max=20_000.0))
    run("Apertado (V_max=2 TWh)", base.with_overrides(V_max=2.0e6))
