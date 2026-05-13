"""Smoke test do módulo dispatch_models — usa séries sintéticas curtas."""
from __future__ import annotations

import numpy as np
import pandas as pd

from validate_model.dispatch_models import (
    AnnualParams,
    DailyParams,
    audit_annual_solution,
    fit_polynomial_regression,
    run_permutation_scenario,
    solve_annual_dispatch,
    solve_daily_dispatch,
    solar_profit,
)


def main() -> None:
    rng = np.random.default_rng(0)
    T = 24 * 30  # 30 dias horários
    t = np.arange(T)
    D = 50_000 + 5_000 * np.sin(2 * np.pi * t / 24) + 1_000 * rng.standard_normal(T)
    g_ren = np.maximum(0, 8_000 + 12_000 * np.sin(2 * np.pi * (t / 24) - np.pi / 3))
    g_ren = np.where(np.cos(2 * np.pi * t / 24) > 0, g_ren, 0)
    A = 25_000 + 5_000 * np.sin(2 * np.pi * t / (24 * 7))

    params = AnnualParams(
        Gh_max=80_000.0,
        Gt_max=60_000.0,
        V_max=2.9e8,
        alpha=0.40,
        terminal_mode="band",
        beta_min=0.20,
        beta_max=0.95,
    )
    res = solve_annual_dispatch(D=D, g_ren=g_ren, A=A, params=params)
    print("annual status:", res["status"], "cost:", res["total_cost"])
    audit = audit_annual_solution(res, A=A, D=D, g_ren=g_ren, tol=1e-4)
    print(audit)

    sin_df = pd.DataFrame({
        "din_instante": pd.date_range("2025-01-01", periods=T, freq="h"),
        "D": D, "gs": g_ren, "gr": np.zeros_like(g_ren),
        "A_MW": A,
    })
    perm = run_permutation_scenario(sin_df, params=params, variables_to_swap={"A_MW": 7})
    print("permutation status:", perm["status"], "cost:", perm["total_cost"])

    fit = fit_polynomial_regression(np.linspace(0, 100, 50), np.linspace(50, 200, 50) + rng.standard_normal(50), 2)
    print("fit:", fit["status"], "r2:", fit["r2"])

    # Daily dispatch
    day = pd.DataFrame({
        "din_instante": pd.date_range("2025-01-01", periods=24, freq="h"),
        "D": 50_000 + 5_000 * np.sin(2 * np.pi * np.arange(24) / 24),
        "gs": np.maximum(0, 12_000 * np.sin(2 * np.pi * np.arange(24) / 24 - np.pi / 3)),
        "gh": np.full(24, 8_000.0),
        "gr": np.zeros(24),
    })
    dparams = DailyParams(K_S1=8_000.0, K_S2=12_000.0, line_limit=4_000.0, line_loss=0.04)
    sol = solve_daily_dispatch(day, day, dparams)
    print("daily status:", sol["status"], "obj:", sol["objective"])
    print(sol["dispatch"][["lambda", "cmo_theoretical", "gt"]].head())
    profit = solar_profit(sol["dispatch"])
    print("profit s1:", profit["pi_s1_total"], "s2:", profit["pi_s2_total"])


if __name__ == "__main__":
    main()
