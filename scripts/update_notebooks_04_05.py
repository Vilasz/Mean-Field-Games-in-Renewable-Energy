from __future__ import annotations

import subprocess
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook


ROOT = Path(__file__).resolve().parents[1]
NB4_PATH = ROOT / "04_despacho_hidrotermico.ipynb"
NB5_PATH = ROOT / "05_prob_diario_cmo_termico.ipynb"


CELL_12 = """@dataclass
class HydroThermalParams:
    Gh_max: float = 0.0
    Gt_max: float = 0.0
    V_max: float = 0.0
    alpha: float = 0.30
    c1: float = 50.0
    c2: float = 1e-3
    pi_d: float = 1e4
    pi_c: float = 10.0
    pi_s: float = 1.0
    beta_final: float | None = None
    final_volume_min: float = 0.05
    final_volume_max: float = 0.95
    terminal_tol: float = 0.01


def calibrate_params(sin_df: pd.DataFrame, hydro_daily: pd.DataFrame) -> HydroThermalParams:
    \"\"\"Calibra parÃ¢metros do modelo a partir dos dados observados.\"\"\"
    Gh_max = float(np.quantile(sin_df[\"gh\"].dropna(), 0.995))
    Gt_max = float(np.quantile(sin_df[\"gn_obs\"].dropna(), 0.999))

    vol_series = sin_df[\"vol_util_pct\"].dropna() if \"vol_util_pct\" in sin_df.columns else pd.Series(dtype=float)
    if len(vol_series) == 0:
        vol_series = hydro_daily[\"vol_util_pct\"].dropna()

    alpha = float(vol_series.iloc[0] / 100.0)
    beta_final = float(vol_series.iloc[-1] / 100.0)

    hd = hydro_daily.dropna(subset=[\"vol_util_pct\", \"A_MW\", \"T_MW\", \"Spill_MW\"]).copy()
    net_daily_MWh = (hd[\"A_MW\"] - hd[\"T_MW\"] - hd[\"Spill_MW\"]) * 24
    E_cum = net_daily_MWh.cumsum().values
    dV_frac = (hd[\"vol_util_pct\"].values - hd[\"vol_util_pct\"].values[0]) / 100.0

    mask = np.abs(E_cum) > 1e3
    if mask.sum() > 10:
        V_max = float(np.sum(E_cum[mask] ** 2) / np.sum(E_cum[mask] * dV_frac[mask]))
        V_max = abs(V_max)
    else:
        cum_surplus = np.cumsum(sin_df[\"A_MW\"].values - sin_df[\"A_MW\"].mean())
        V_max = 2.0 * (cum_surplus.max() - cum_surplus.min())

    V_max = max(V_max, Gh_max * 24 * 7)

    return HydroThermalParams(
        Gh_max=Gh_max,
        Gt_max=Gt_max,
        V_max=V_max,
        alpha=alpha,
        beta_final=beta_final,
    )


params = calibrate_params(sin, hydro_daily)

print(\"ParÃ¢metros calibrados:\")
print(f\"  Gh_max      = {params.Gh_max:>12,.0f} MW\")
print(f\"  Gt_max      = {params.Gt_max:>12,.0f} MW\")
print(f\"  V_max       = {params.V_max:>12,.0f} MWh  ({params.V_max/1e6:.2f} TWh)\")
print(f\"  V_0         = {params.alpha * params.V_max:>12,.0f} MWh  ({params.alpha*100:.1f}% de V_max)\")
print(f\"  V_T alvo    = {params.beta_final * params.V_max:>12,.0f} MWh  ({params.beta_final*100:.1f}% de V_max)\")
print(f\"  Î±           = {params.alpha:.3f}\")
print(f\"  Î²_final     = {params.beta_final:.3f}  (meta terminal observada no perÃ­odo)\")
print(f\"  c1          = {params.c1} R$/MWh\")
print(f\"  c2          = {params.c2} R$/MWÂ²h\")
print(f\"  Ï€_d         = {params.pi_d:,.0f} R$/MWh\")
print(f\"  Ï€_c         = {params.pi_c} R$/MWh\")
print(f\"  Ï€_s         = {params.pi_s} R$/MWh (vertimento)\")
print(\"  Estado s_t  = V_t / 1e6  (armazenamento em TWh para estabilidade numerica)\")
"""


CELL_14 = """def solve_hydrothermal(
    D: np.ndarray,
    g_ren: np.ndarray,
    A: np.ndarray,
    params: HydroThermalParams,
    verbose: bool = False,
) -> dict:
    \"\"\"Resolve o despacho hidrotÃ©rmico via QP convexo (CVXPY).\"\"\"
    T = len(D)
    D = np.asarray(D, dtype=float)
    g_ren = np.asarray(g_ren, dtype=float)
    A = np.asarray(A, dtype=float)

    gh = cp.Variable(T, name=\"gh\")
    c_var = cp.Variable(T, name=\"c\")
    d_var = cp.Variable(T, name=\"d\")
    spill = cp.Variable(T, name=\"spill\")
    s = cp.Variable(T + 1, name=\"storage_twh\")

    gt = D - gh - g_ren + c_var - d_var

    cost = params.c1 * cp.sum(gt) + 0.5 * params.c2 * cp.sum_squares(gt)
    cost += params.pi_d * cp.sum(d_var)
    cost += params.pi_c * cp.sum(c_var)
    cost += params.pi_s * cp.sum(spill)

    twh = 1e6
    S_max = params.V_max / twh
    S_0 = params.alpha * params.V_max / twh
    reservoir_dyn = (s[1:] == s[:-1] + (A - gh - spill) / twh)

    terminal_floor = params.final_volume_min
    terminal_cap = params.final_volume_max
    if params.beta_final is not None:
        terminal_floor = max(terminal_floor, params.beta_final - params.terminal_tol)
        terminal_cap = min(terminal_cap, params.beta_final + params.terminal_tol)
    S_floor = terminal_floor * params.V_max / twh
    S_cap = terminal_cap * params.V_max / twh

    cons = [
        gh >= 0, gh <= params.Gh_max,
        gt >= 0, gt <= params.Gt_max,
        c_var >= 0, d_var >= 0, spill >= 0,
        s >= 0, s <= S_max,
        s[0] == S_0,
        reservoir_dyn,
        s[-1] >= S_floor,
        s[-1] <= S_cap,
    ]

    prob = cp.Problem(cp.Minimize(cost), cons)

    solvers = [cp.OSQP, cp.SCS]
    if hasattr(cp, \"CLARABEL\"):
        solvers.append(cp.CLARABEL)
    if hasattr(cp, \"ECOS\"):
        solvers.append(cp.ECOS)

    for solver in solvers:
        try:
            kw = {\"verbose\": verbose}
            if solver == cp.SCS:
                kw.update({\"max_iters\": 100_000, \"eps\": 1e-5})
            prob.solve(solver=solver, **kw)
            if prob.status in (\"optimal\", \"optimal_inaccurate\"):
                break
        except Exception:
            continue

    def _val(x, n):
        return np.asarray(x.value).ravel() if x.value is not None else np.full(n, np.nan)

    S = _val(s, T + 1)
    V = S * twh
    v_frac = V / params.V_max

    wv = None
    if reservoir_dyn.dual_value is not None:
        wv = -np.asarray(reservoir_dyn.dual_value).ravel() / twh

    return {
        \"gh\": _val(gh, T),
        \"gt\": np.asarray(gt.value).ravel() if gt.value is not None else np.full(T, np.nan),
        \"c\": _val(c_var, T),
        \"d\": _val(d_var, T),
        \"spill\": _val(spill, T),
        \"V\": V,
        \"v_frac\": v_frac,
        \"water_value\": wv,
        \"total_cost\": float(prob.value) if prob.value is not None else np.nan,
        \"status\": prob.status,
    }


print(\"Solver pronto.\")
"""


CELL_29 = """alphas = [0.10, 0.20, 0.30, 0.50, 0.70, 0.90]
sens_a = []


def safe_stat(x, func):
    arr = np.asarray(x, dtype=float)
    arr = arr[np.isfinite(arr)]
    return func(arr) if arr.size else np.nan


for a in alphas:
    p = HydroThermalParams(
        Gh_max=params.Gh_max,
        Gt_max=params.Gt_max,
        V_max=params.V_max,
        alpha=a,
        c1=params.c1,
        c2=params.c2,
        pi_d=params.pi_d,
        pi_c=params.pi_c,
        pi_s=params.pi_s,
        beta_final=None,
        final_volume_min=params.final_volume_min,
        final_volume_max=params.final_volume_max,
        terminal_tol=params.terminal_tol,
    )
    r = solve_hydrothermal(sin[\"D\"].values, sin[\"g_ren\"].values, sin[\"A_MW\"].values, p)
    gt = np.asarray(r[\"gt\"], dtype=float)
    gh = np.asarray(r[\"gh\"], dtype=float)
    d = np.maximum(np.asarray(r[\"d\"], dtype=float), 0)
    sp = np.maximum(np.asarray(r[\"spill\"], dtype=float), 0)
    v_frac = np.asarray(r[\"v_frac\"], dtype=float)

    sens_a.append({
        \"alpha\": a,
        \"status\": r[\"status\"],
        \"cost_bi\": r[\"total_cost\"] / 1e9,
        \"gt_mean\": safe_stat(gt, np.mean),
        \"gt_p95\": safe_stat(gt, lambda z: np.quantile(z, 0.95)),
        \"gh_mean\": safe_stat(gh, np.mean),
        \"deficit_MWh\": safe_stat(d, np.sum),
        \"spill_MWh\": safe_stat(sp, np.sum),
        \"V_min_pct\": safe_stat(v_frac, np.min) * 100,
        \"V_final_pct\": safe_stat(v_frac[-1:], np.mean) * 100,
        \"V_traj\": r[\"V\"],
    })

    gt_mean_val = safe_stat(gt, np.mean)
    gt_mean_txt = f\"{gt_mean_val:,.0f} MW\" if np.isfinite(gt_mean_val) else \"n/a\"
    print(
        f\"  alpha={a:.2f} -> {r['total_cost']/1e9:.3f} bi, \"
        f\"gt={gt_mean_txt}, V_T={safe_stat(v_frac[-1:], np.mean)*100:.1f}% [{r['status']}]\" + "\\n"
    )

df_sa = pd.DataFrame([{k: v for k, v in s.items() if k != \"V_traj\"} for s in sens_a])
display(df_sa.round(3))

print(
    \"Leitura rapida: para a sensibilidade de alpha, a meta terminal foi relaxada de proposito, \"
    \"isolando o efeito do estoque inicial e eliminando inviabilidades artificiais.\"
)
"""


CELL_30 = """fig, axes = plt.subplots(1, 3, figsize=(17, 5))

ax = axes[0]
ax.plot(df_sa[\"alpha\"] * 100, df_sa[\"cost_bi\"], \"o-\", color=\"indianred\", lw=2, ms=8)
ax.set_xlabel(\"Î± (%)\"); ax.set_ylabel(\"Custo (R$ bi)\")
ax.set_title(\"Custo Total vs Î±\"); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(df_sa[\"alpha\"] * 100, df_sa[\"gt_mean\"], \"s-\", color=\"indianred\", lw=2, ms=8, label=\"TÃ©rmica mÃ©dia\")
ax.plot(df_sa[\"alpha\"] * 100, df_sa[\"gt_p95\"], \"d--\", color=\"firebrick\", lw=1.6, ms=6, label=\"TÃ©rmica p95\")
ax.plot(df_sa[\"alpha\"] * 100, df_sa[\"gh_mean\"], \"^-\", color=\"cyan\", lw=2, ms=8, label=\"Hidro mÃ©dia\")
ax.set_xlabel(\"Î± (%)\"); ax.set_ylabel(\"GeraÃ§Ã£o (MW)\")
ax.set_title(\"GeraÃ§Ã£o MÃ©dia e p95 vs Î±\")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[2]
for s in sens_a:
    vp = np.asarray(s[\"V_traj\"])[:len(sin)] / params.V_max * 100
    ax.plot(sin[\"din_instante\"], vp, lw=1.2, label=f\"Î±={s['alpha']*100:.0f}%\")
ax.set_ylabel(\"Volume (%)\"); ax.set_title(\"TrajetÃ³ria do ReservatÃ³rio por Î±\")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

fig.tight_layout(); plt.show()
"""


PERM_MD = """### 8.2 PermutaÃ§Ãµes de Dias â€” Hidro e Solar

Nesta etapa fazemos quatro permutaÃ§Ãµes simples em blocos inteiros de dias para os perfis de afluÃªncia hidrelÃ©trica e geraÃ§Ã£o solar, preservando a estrutura intradiÃ¡ria de cada fonte. A ideia Ã© deslocar a oferta renovÃ¡vel ao redor do calendÃ¡rio observado e medir como a geraÃ§Ã£o tÃ©rmica Ã³tima responde.

ConvenÃ§Ã£o usada abaixo:

- `hidro +10d`: a afluÃªncia do dia *t* passa a usar o perfil observado 10 dias Ã  frente.
- `solar -10d`: a geraÃ§Ã£o solar do dia *t* passa a usar o perfil observado 10 dias para trÃ¡s.
- As permutaÃ§Ãµes sÃ£o circulares e feitas em passos de 24 horas.
"""


PERM_CODE_1 = """def shift_profile_days(x, days):
    steps = int(days * 24)
    return np.roll(np.asarray(x, dtype=float), steps)


perm_specs = [
    {\"scenario\": \"Base\", \"shift_hydro_d\": 0, \"shift_solar_d\": 0},
    {\"scenario\": \"H+10 | S-10\", \"shift_hydro_d\": 10, \"shift_solar_d\": -10},
    {\"scenario\": \"H+15 | S-15\", \"shift_hydro_d\": 15, \"shift_solar_d\": -15},
    {\"scenario\": \"H-10 | S+10\", \"shift_hydro_d\": -10, \"shift_solar_d\": 10},
    {\"scenario\": \"H-15 | S+15\", \"shift_hydro_d\": -15, \"shift_solar_d\": 15},
]

perm_results = []
for spec in perm_specs:
    A_perm = shift_profile_days(sin[\"A_MW\"].values, spec[\"shift_hydro_d\"])
    gs_perm = shift_profile_days(sin[\"gs\"].values, spec[\"shift_solar_d\"])
    g_ren_perm = gs_perm + sin[\"gr\"].values

    r = solve_hydrothermal(
        D=sin[\"D\"].values,
        g_ren=g_ren_perm,
        A=A_perm,
        params=params,
    )

    gt = np.asarray(r[\"gt\"], dtype=float)
    perm_results.append({
        \"scenario\": spec[\"scenario\"],
        \"shift_hydro_d\": spec[\"shift_hydro_d\"],
        \"shift_solar_d\": spec[\"shift_solar_d\"],
        \"status\": r[\"status\"],
        \"cost_bi\": r[\"total_cost\"] / 1e9,
        \"gt_mean\": safe_stat(gt, np.mean),
        \"gt_p95\": safe_stat(gt, lambda z: np.quantile(z, 0.95)),
        \"gt_peak\": safe_stat(gt, np.max),
        \"gt_energy_GWh\": safe_stat(gt, np.sum) / 1e3,
        \"V_final_pct\": safe_stat(r[\"v_frac\"][-1:], np.mean) * 100,
        \"gt_series\": gt,
    })

df_perm = pd.DataFrame([{k: v for k, v in x.items() if k != \"gt_series\"} for x in perm_results])
display(df_perm.round(3))
"""


PERM_CODE_2 = """fig, axes = plt.subplots(1, 3, figsize=(18, 5))

ax = axes[0]
ax.bar(df_perm[\"scenario\"], df_perm[\"gt_mean\"], color=\"indianred\", alpha=0.8)
ax.set_ylabel(\"MW\"); ax.set_title(\"GeraÃ§Ã£o TÃ©rmica MÃ©dia\")
ax.tick_params(axis=\"x\", rotation=20); ax.grid(True, alpha=0.3, axis=\"y\")

ax = axes[1]
ax.bar(df_perm[\"scenario\"], df_perm[\"gt_p95\"], color=\"firebrick\", alpha=0.8, label=\"p95\")
ax.bar(df_perm[\"scenario\"], df_perm[\"gt_peak\"], color=\"salmon\", alpha=0.55, label=\"pico\")
ax.set_ylabel(\"MW\"); ax.set_title(\"TÃ©rmica p95 e Pico\")
ax.tick_params(axis=\"x\", rotation=20); ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, axis=\"y\")

ax = axes[2]
for row in perm_results:
    q = np.sort(np.asarray(row[\"gt_series\"], dtype=float))[::-1]
    ax.plot(q, lw=1.6, label=row[\"scenario\"])
ax.set_xlabel(\"Horas ordenadas\"); ax.set_ylabel(\"MW\")
ax.set_title(\"Curvas de DuraÃ§Ã£o da TÃ©rmica\")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

fig.tight_layout(); plt.show()

base_gt = next(row[\"gt_series\"] for row in perm_results if row[\"scenario\"] == \"Base\")
fig, ax = plt.subplots(figsize=(16, 4))
ax.plot(sin[\"din_instante\"].iloc[:24*21], base_gt[:24*21], lw=2.0, color=\"black\", label=\"Base\")
for row in perm_results[1:]:
    ax.plot(sin[\"din_instante\"].iloc[:24*21], row[\"gt_series\"][:24*21], lw=1.3, label=row[\"scenario\"])
ax.set_ylabel(\"MW\"); ax.set_title(\"Primeiras 3 Semanas â€” Resposta da TÃ©rmica Ã s PermutaÃ§Ãµes\")
ax.legend(ncol=3, fontsize=8); ax.grid(True, alpha=0.3)
plt.show()

best_perm = df_perm.loc[df_perm[\"gt_mean\"].idxmin()]
worst_perm = df_perm.loc[df_perm[\"gt_mean\"].idxmax()]
print(f\"Melhor cenÃ¡rio para a tÃ©rmica mÃ©dia: {best_perm['scenario']} ({best_perm['gt_mean']:.1f} MW)\")
print(f\"Pior cenÃ¡rio para a tÃ©rmica mÃ©dia:  {worst_perm['scenario']} ({worst_perm['gt_mean']:.1f} MW)\")
"""


CELL_38 = """print(\"=\" * 72)
print(\"  DESPACHO HIDROTERMICO OTIMO - RESUMO\")
print(\"=\" * 72)
print(f\"\\nPeriodo: {sin['din_instante'].min()} -> {sin['din_instante'].max()}\")
print(f\"T = {len(sin):,} horas | Resolucao: horaria\")
print(f\"Dados de afluencia: hidrologicos_di_2025.csv (175 reservatorios, diario)\")
print(\"\\nObservacao metodologica: o estado do reservatorio foi reescalonado para v_t = V_t / V_max,\")
print(\"o que melhora a estabilidade numerica e elimina inviabilidades artificiais na sensibilidade de alpha.\")

print(f\"\\n{'PARAMETRO':<28} {'VALOR':>15}\")
print(\"-\" * 45)
print(f\"{'Gh_max (MW)':<28} {params.Gh_max:>15,.0f}\")
print(f\"{'Gt_max (MW)':<28} {params.Gt_max:>15,.0f}\")
print(f\"{'V_max (TWh)':<28} {params.V_max/1e6:>15.2f}\")
print(f\"{'alpha inicial':<28} {params.alpha:>15.1%}\")
print(f\"{'beta final alvo':<28} {params.beta_final:>15.1%}\")
print(f\"{'c1 (R$/MWh)':<28} {params.c1:>15.1f}\")
print(f\"{'c2 (R$/MW2h)':<28} {params.c2:>15.1e}\")
print(f\"{'pi_d (R$/MWh)':<28} {params.pi_d:>15,.0f}\")
print(f\"{'pi_c (R$/MWh)':<28} {params.pi_c:>15.1f}\")
print(f\"{'pi_s (R$/MWh)':<28} {params.pi_s:>15.1f}\")

print(f\"\\n{'METRICA':<38} {'OTIMO':>12} {'OBSERVADO':>12}\")
print(\"-\" * 64)
print(f\"{'Geracao termica media (MW)':<38} {sin['gt_opt'].mean():>12,.0f} {sin['gn_obs'].mean():>12,.0f}\")
print(f\"{'Geracao hidro media (MW)':<38} {sin['gh_opt'].mean():>12,.0f} {sin['gh'].mean():>12,.0f}\")
print(f\"{'Custo termico (R$ bi)':<38} {sin['cost_th_opt'].sum()/1e9:>12.3f} {sin['cost_th_obs'].sum()/1e9:>12.3f}\")
print(f\"{'Horas c/ deficit':<38} {(sin['d_opt']>1e-2).sum():>12,} {'-':>12}\")
print(f\"{'Horas c/ curtailment':<38} {(sin['c_opt']>1e-2).sum():>12,} {'-':>12}\")
print(f\"{'Horas c/ vertimento':<38} {(sin['spill_opt']>1e-2).sum():>12,} {'-':>12}\")
print(f\"{'Vertimento total (% afluencia)':<38} {sin['spill_opt'].sum()/sin['A_MW'].sum()*100:>12.1f}% {'-':>12}\")
print(f\"{'Volume min (% V_max)':<38} {sin['V_opt'].min()/params.V_max*100:>12.1f} {'-':>12}\")
print(f\"{'Volume final (% V_max)':<38} {sin['V_opt'].iloc[-1]/params.V_max*100:>12.1f} {'-':>12}\")

econ = (sin['cost_th_obs'].sum() - sin['cost_th_opt'].sum()) / 1e9
econ_pct = econ / (sin['cost_th_obs'].sum()/1e9) * 100
print(f\"\\nEconomia potencial: R$ {econ:.3f} bilhoes ({econ_pct:.1f}%)\")
print(\"=\" * 72)
"""


def update_notebook_4() -> None:
    nb = nbformat.read(NB4_PATH, as_version=4)
    base_raw = subprocess.check_output(["git", "show", "HEAD:04_despacho_hidrotermico.ipynb"], cwd=ROOT)
    base_nb = nbformat.reads(base_raw.decode("utf-8"), as_version=4)
    nb.cells[12].source = CELL_12
    nb.cells[14].source = CELL_14
    nb.cells[29].source = CELL_29
    nb.cells[30].source = CELL_30

    insertion = [
        new_markdown_cell(PERM_MD),
        new_code_cell(PERM_CODE_1),
        new_code_cell(PERM_CODE_2),
    ]
    tail = [
        base_nb.cells[31],
        base_nb.cells[32],
        base_nb.cells[33],
        base_nb.cells[34],
        base_nb.cells[35],
        base_nb.cells[36],
        base_nb.cells[37],
        new_code_cell(CELL_38),
        base_nb.cells[39],
    ]
    nb.cells = nb.cells[:31] + insertion + tail
    nbformat.write(nb, NB4_PATH)


def create_notebook_5() -> None:
    meta = nbformat.read(NB4_PATH, as_version=4).metadata
    cells = []

    cells.append(new_markdown_cell("""# Problema Diario, CMO Termico e Permutacoes Sazonais

Este notebook constroi um estudo diario simples e direto, mas ancorado nos dados do SIN 2025. A ideia e observar a termica como variavel central de decisao e interpretar o CMO a partir do custo marginal da geracao termica.

Modelos usados:

- **Visao do operador/planner:** minimizar `sum(gt_t^2)` sujeito ao balanco `gt_t + gh_t + gs1_t + gs2_t = D_t`.
- **Visao do agente termico:** maximizar `pi_t gt_t - gt_t^2`, o que gera a oferta otima `gt_t = pi_t / 2`.
- **Ligacao com o CMO:** quando a termica e marginal, temos `CMO_t ~= pi_t = 2 gt_t`.

O sistema estudado tera:

- dois paineis solares com capacidades distintas;
- uma usina hidreletrica agregada com orcamento diario de agua;
- uma usina termica agregada;
- um bloco simples de transmissao para o segundo painel solar.
"""))

    cells.append(new_code_cell("""import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cvxpy as cp

from validate_model.pipeline import SINPaths, build_panel, load_hidrologia

pd.set_option("display.max_columns", 200)
plt.rcParams.update({
    "figure.dpi": 110,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "font.size": 10,
})

paths = SINPaths(root="validate_model", year=2025)
"""))

    cells.append(new_markdown_cell("## 1. Base SIN 2025 e selecao de dias representativos"))

    cells.append(new_code_cell("""panel = build_panel(paths)

sin = (
    panel
    .groupby("din_instante", as_index=False)
    .agg({
        "D": "sum",
        "gs": "sum",
        "gr": "sum",
        "gh": "sum",
        "gn_obs": "sum",
    })
    .sort_values("din_instante")
    .reset_index(drop=True)
)

hydro_daily = load_hidrologia(paths.hidro_di_path)
hydro_daily["date"] = pd.to_datetime(hydro_daily["din_instante"]).dt.normalize()
sin["date"] = sin["din_instante"].dt.normalize()
sin = sin.merge(hydro_daily[["date", "A_MW"]], on="date", how="left")
sin["A_MW"] = sin["A_MW"].interpolate(method="linear").bfill().ffill()

season_map = {
    12: "Verao", 1: "Verao", 2: "Verao",
    3: "Outono", 4: "Outono", 5: "Outono",
    6: "Inverno", 7: "Inverno", 8: "Inverno",
    9: "Primavera", 10: "Primavera", 11: "Primavera",
}
sin["season"] = sin["din_instante"].dt.month.map(season_map)
sin["day_type"] = np.where(sin["din_instante"].dt.weekday < 5, "Dia util", "Fim de semana")

daily_meta = (
    sin.groupby("date")
    .agg(
        demand_GWh=("D", "sum"),
        solar_GWh=("gs", "sum"),
        hydro_GWh=("gh", "sum"),
        season=("season", "first"),
        day_type=("day_type", "first"),
        n_hours=("D", "size"),
    )
    .reset_index()
)
daily_meta = daily_meta[daily_meta["n_hours"] == 24].copy()


def representative_day(meta: pd.DataFrame, season: str, day_type: str) -> pd.Timestamp:
    group = meta[(meta["season"] == season) & (meta["day_type"] == day_type)].copy()
    if group.empty:
        raise ValueError(f"Sem dados para {season} | {day_type}")
    med = group["demand_GWh"].median()
    idx = (group["demand_GWh"] - med).abs().idxmin()
    return pd.Timestamp(group.loc[idx, "date"])


rep_rows = []
for season in ["Verao", "Outono", "Inverno", "Primavera"]:
    for day_type in ["Dia util", "Fim de semana"]:
        day = representative_day(daily_meta, season, day_type)
        row = daily_meta.loc[daily_meta["date"] == day].iloc[0]
        rep_rows.append({
            "label": f"{season} | {day_type}",
            "date": day,
            "demand_GWh": row["demand_GWh"] / 1e3,
            "solar_GWh": row["solar_GWh"] / 1e3,
            "hydro_GWh": row["hydro_GWh"] / 1e3,
        })

rep_days = pd.DataFrame(rep_rows)
display(rep_days.round(3))
"""))

    cells.append(new_code_cell("""fig, axes = plt.subplots(1, 2, figsize=(14, 4))

ax = axes[0]
sns.barplot(data=rep_days, x="label", y="demand_GWh", hue="label", dodge=False, ax=ax, palette="crest")
ax.set_title("Demanda diaria dos dias representativos")
ax.set_ylabel("GWh")
ax.set_xlabel("")
ax.tick_params(axis="x", rotation=35)
if ax.legend_ is not None:
    ax.legend_.remove()

ax = axes[1]
tmp = rep_days.melt(
    id_vars=["label"],
    value_vars=["solar_GWh", "hydro_GWh"],
    var_name="fonte",
    value_name="energia_GWh",
)
sns.barplot(data=tmp, x="label", y="energia_GWh", hue="fonte", ax=ax, palette=["gold", "deepskyblue"])
ax.set_title("Oferta solar e hidrica por dia representativo")
ax.set_ylabel("GWh")
ax.set_xlabel("")
ax.tick_params(axis="x", rotation=35)

fig.tight_layout(); plt.show()
"""))

    cells.append(new_markdown_cell("""## 2. Formulacao diaria do problema

Vamos usar uma versao diaria muito simples:

- `D_t`: demanda horaria;
- `gs1_t` e `gs2_t`: producao dos dois paineis solares;
- `gh_t`: geracao hidreletrica;
- `gt_t`: geracao termica;
- `shed_t`: deficit, penalizado fortemente;
- `curt_t`: corte solar, com penalidade baixa.

Problema diario:

`min sum(gt_t^2) + M * sum(shed_t) + eps * sum(curt_t)`

sujeito a:

- `gt_t + gh_t + gs1_t + gs2_t + shed_t - curt_t = D_t`
- `0 <= gh_t <= Gh_cap`
- `sum(gh_t) <= H_budget`
- `0 <= gt_t <= Gt_max`

Depois de resolver o despacho, inferimos `pi_t = 2 gt_t`, logo `CMO_t = pi_t`.
"""))

    cells.append(new_code_cell("""K_S1 = 10_000.0
K_S2 = 16_000.0
GT_MAX = 50_000.0
HYDRO_SCALE = 1.00
BIG_M = 2e5
EPS_CURT = 5.0


def make_solar_profiles(day_supply: pd.DataFrame, k1: float = K_S1, k2: float = K_S2) -> tuple[np.ndarray, np.ndarray]:
    raw = day_supply["gs"].to_numpy(dtype=float)
    base_cf = raw / max(raw.max(), 1.0)
    cf1 = np.clip(base_cf, 0, 1)
    cf2 = np.clip(0.90 * np.roll(base_cf, 1) + 0.10 * base_cf, 0, 1)
    return k1 * cf1, k2 * cf2


def solve_daily_dispatch(
    day_demand: pd.DataFrame,
    day_supply: pd.DataFrame,
    line_limit: float = np.inf,
    line_loss: float = 0.0,
    hydro_scale: float = HYDRO_SCALE,
    thermal_cap: float = GT_MAX,
) -> dict:
    D = day_demand["D"].to_numpy(dtype=float)
    gs1_raw, gs2_raw = make_solar_profiles(day_supply)
    gs2_export = np.minimum(gs2_raw, line_limit)
    gs2_delivered = (1 - line_loss) * gs2_export
    solar_total = gs1_raw + gs2_delivered

    hydro_energy_budget = float(day_supply["gh"].sum() * hydro_scale)
    gh_cap = np.full(24, max(day_supply["gh"].max() * 1.05, hydro_energy_budget / 24 * 1.10))

    gt = cp.Variable(24, nonneg=True, name="gt")
    gh = cp.Variable(24, nonneg=True, name="gh")
    shed = cp.Variable(24, nonneg=True, name="shed")
    curt = cp.Variable(24, nonneg=True, name="curt")

    cons = [
        gt <= thermal_cap,
        gh <= gh_cap,
        cp.sum(gh) <= hydro_energy_budget,
        curt <= solar_total,
        gt + gh + solar_total + shed - curt == D,
    ]

    obj = cp.Minimize(cp.sum_squares(gt) + BIG_M * cp.sum(shed) + EPS_CURT * cp.sum(curt))
    prob = cp.Problem(obj, cons)

    for solver in [cp.OSQP, cp.SCS]:
        try:
            kw = {"verbose": False}
            if solver == cp.SCS:
                kw.update({"max_iters": 50_000, "eps": 1e-5})
            prob.solve(solver=solver, **kw)
            if prob.status in ("optimal", "optimal_inaccurate"):
                break
        except Exception:
            continue

    gt_val = np.asarray(gt.value).ravel()
    gh_val = np.asarray(gh.value).ravel()
    shed_val = np.asarray(shed.value).ravel()
    curt_val = np.asarray(curt.value).ravel()
    cmo = 2 * gt_val

    out = pd.DataFrame({
        "din_instante": day_demand["din_instante"].to_numpy(),
        "D": D,
        "gs1": gs1_raw,
        "gs2": gs2_delivered,
        "gh": gh_val,
        "gt": gt_val,
        "CMO": cmo,
        "shed": shed_val,
        "curt": curt_val,
        "gs2_raw": gs2_raw,
        "gs2_export": gs2_export,
    })

    return {
        "status": prob.status,
        "objective": float(prob.value),
        "dispatch": out,
        "hydro_budget_GWh": hydro_energy_budget / 1e3,
        "solar2_spilled_GWh": float(np.sum(gs2_raw - gs2_export) / 1e3),
        "line_limit": line_limit,
        "line_loss": line_loss,
    }
"""))

    cells.append(new_markdown_cell("## 3. Cenarios base por estacao e por tipo de dia"))

    cells.append(new_code_cell("""def day_frame(day: pd.Timestamp) -> pd.DataFrame:
    return sin[sin["date"] == pd.Timestamp(day)].copy().sort_values("din_instante").reset_index(drop=True)


base_runs = []
base_dispatch = {}
for _, row in rep_days.iterrows():
    label = row["label"]
    date = pd.Timestamp(row["date"])
    day = day_frame(date)
    solved = solve_daily_dispatch(day, day)
    disp = solved["dispatch"]
    base_dispatch[label] = disp
    base_runs.append({
        "label": label,
        "date": date,
        "status": solved["status"],
        "gt_mean": disp["gt"].mean(),
        "gt_peak": disp["gt"].max(),
        "gt_energy_GWh": disp["gt"].sum() / 1e3,
        "gh_energy_GWh": disp["gh"].sum() / 1e3,
        "solar_energy_GWh": (disp["gs1"].sum() + disp["gs2"].sum()) / 1e3,
        "CMO_mean": disp["CMO"].mean(),
        "CMO_peak": disp["CMO"].max(),
        "shed_MWh": disp["shed"].sum(),
        "curt_GWh": disp["curt"].sum() / 1e3,
    })

df_base_daily = pd.DataFrame(base_runs)
display(df_base_daily.round(3))
"""))

    cells.append(new_code_cell("""plot_df = df_base_daily.copy()
plot_df["season"] = plot_df["label"].str.split(" | ", regex=False).str[0]
plot_df["day_type"] = plot_df["label"].str.split(" | ", regex=False).str[1]

pivot_gt = plot_df.pivot(index="season", columns="day_type", values="gt_energy_GWh")
pivot_cmo = plot_df.pivot(index="season", columns="day_type", values="CMO_mean")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
sns.heatmap(pivot_gt, annot=True, fmt=".1f", cmap="Reds", ax=axes[0])
axes[0].set_title("Energia termica diaria (GWh)")

sns.heatmap(pivot_cmo, annot=True, fmt=".0f", cmap="Oranges", ax=axes[1])
axes[1].set_title("CMO medio diario (R$/MWh)")

fig.tight_layout(); plt.show()
"""))

    cells.append(new_markdown_cell("""## 4. Permutacoes de dias

Agora separamos **demanda** e **oferta**. A demanda vem de um dia representativo e a oferta hidrica/solar vem de outro. Isso cria um experimento simples e util:

- se o perfil de oferta for deslocado para uma condicao mais favoravel, a termica cai;
- se a oferta for deslocada para uma condicao menos favoravel, a termica sobe e o CMO acompanha.

As quatro permutacoes abaixo cobrem trocas sazonais e tambem trocas entre dia util e fim de semana.
"""))

    cells.append(new_code_cell("""perm_day_specs = [
    {"scenario": "Demanda Verao util | Oferta Inverno util", "demand_label": "Verao | Dia util", "supply_label": "Inverno | Dia util"},
    {"scenario": "Demanda Inverno util | Oferta Verao util", "demand_label": "Inverno | Dia util", "supply_label": "Verao | Dia util"},
    {"scenario": "Demanda Primavera util | Oferta Primavera FDS", "demand_label": "Primavera | Dia util", "supply_label": "Primavera | Fim de semana"},
    {"scenario": "Demanda Outono FDS | Oferta Outono util", "demand_label": "Outono | Fim de semana", "supply_label": "Outono | Dia util"},
]

perm_day_results = []
perm_dispatch = {}

for spec in perm_day_specs:
    demand_date = pd.Timestamp(rep_days.loc[rep_days["label"] == spec["demand_label"], "date"].iloc[0])
    supply_date = pd.Timestamp(rep_days.loc[rep_days["label"] == spec["supply_label"], "date"].iloc[0])

    demand_day = day_frame(demand_date)
    supply_day = day_frame(supply_date)
    solved = solve_daily_dispatch(demand_day, supply_day)
    disp = solved["dispatch"]
    perm_dispatch[spec["scenario"]] = disp

    perm_day_results.append({
        "scenario": spec["scenario"],
        "demand_date": demand_date.date(),
        "supply_date": supply_date.date(),
        "status": solved["status"],
        "gt_mean": disp["gt"].mean(),
        "gt_peak": disp["gt"].max(),
        "gt_energy_GWh": disp["gt"].sum() / 1e3,
        "CMO_mean": disp["CMO"].mean(),
        "CMO_peak": disp["CMO"].max(),
        "shed_MWh": disp["shed"].sum(),
        "curt_GWh": disp["curt"].sum() / 1e3,
    })

df_perm_daily = pd.DataFrame(perm_day_results)
display(df_perm_daily.round(3))
"""))

    cells.append(new_code_cell("""fig, axes = plt.subplots(2, 2, figsize=(15, 8), sharex=True, sharey=True)

for ax, (scenario, disp) in zip(axes.ravel(), perm_dispatch.items()):
    ax.plot(disp["din_instante"], disp["D"], color="black", lw=1.7, label="Demanda")
    ax.fill_between(disp["din_instante"], 0, disp["gs1"], color="gold", alpha=0.65, label="Solar 1")
    ax.fill_between(disp["din_instante"], disp["gs1"], disp["gs1"] + disp["gs2"], color="orange", alpha=0.55, label="Solar 2")
    ax.fill_between(disp["din_instante"], disp["gs1"] + disp["gs2"], disp["gs1"] + disp["gs2"] + disp["gh"], color="deepskyblue", alpha=0.55, label="Hidro")
    ax.fill_between(disp["din_instante"], disp["gs1"] + disp["gs2"] + disp["gh"], disp["gs1"] + disp["gs2"] + disp["gh"] + disp["gt"], color="indianred", alpha=0.60, label="Termica")
    ax.set_title(scenario)
    ax.grid(True, alpha=0.25)
    ax.tick_params(axis="x", rotation=25)

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles[:5], labels[:5], loc="upper center", ncol=5, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.95)); plt.show()

fig, ax = plt.subplots(figsize=(12, 4))
for scenario, disp in perm_dispatch.items():
    ax.plot(disp["din_instante"], disp["CMO"], lw=1.8, label=scenario)
ax.set_title("CMO horario nas permutacoes diarias")
ax.set_ylabel("R$/MWh")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, ncol=2)
plt.show()
"""))

    cells.append(new_markdown_cell("""## 5. Transmissao: efeito simples no sistema

Para incorporar transmissao da maneira mais simples possivel, vamos supor que o **Painel Solar 2** e remoto e precisa atravessar uma linha ate o centro de carga.

Representacao adotada:

- `line_limit`: limite de exportacao da linha;
- `line_loss`: perda percentual na energia exportada;
- quanto menor a capacidade da linha ou maior a perda, menor a solar entregue e maior a necessidade de termica.

Isso nao substitui um modelo de fluxo em rede, mas ja mostra o mecanismo economico principal: restricao de transmissao desloca o equilibrio para mais termica e CMO maior.
"""))

    cells.append(new_code_cell("""stress_case = df_perm_daily.sort_values("gt_peak", ascending=False).iloc[0]
stress_name = stress_case["scenario"]
stress_spec = next(spec for spec in perm_day_specs if spec["scenario"] == stress_name)

demand_date = pd.Timestamp(rep_days.loc[rep_days["label"] == stress_spec["demand_label"], "date"].iloc[0])
supply_date = pd.Timestamp(rep_days.loc[rep_days["label"] == stress_spec["supply_label"], "date"].iloc[0])

demand_day = day_frame(demand_date)
supply_day = day_frame(supply_date)

tx_specs = [
    {"name": "Sem restricao", "line_limit": np.inf, "line_loss": 0.00},
    {"name": "Linha folgada", "line_limit": 5_000.0, "line_loss": 0.02},
    {"name": "Linha intermediaria", "line_limit": 3_500.0, "line_loss": 0.04},
    {"name": "Linha restrita", "line_limit": 2_000.0, "line_loss": 0.06},
]

tx_results = []
tx_dispatch = {}
for spec in tx_specs:
    solved = solve_daily_dispatch(
        demand_day,
        supply_day,
        line_limit=spec["line_limit"],
        line_loss=spec["line_loss"],
    )
    disp = solved["dispatch"]
    tx_dispatch[spec["name"]] = disp
    tx_results.append({
        "cenario": spec["name"],
        "gt_mean": disp["gt"].mean(),
        "gt_peak": disp["gt"].max(),
        "gt_energy_GWh": disp["gt"].sum() / 1e3,
        "CMO_mean": disp["CMO"].mean(),
        "CMO_peak": disp["CMO"].max(),
        "shed_MWh": disp["shed"].sum(),
        "curt_GWh": disp["curt"].sum() / 1e3,
        "solar2_spilled_GWh": solved["solar2_spilled_GWh"],
    })

df_tx = pd.DataFrame(tx_results)
display(df_tx.round(3))
"""))

    cells.append(new_code_cell("""fig, axes = plt.subplots(1, 2, figsize=(13, 4))

ax = axes[0]
for name, disp in tx_dispatch.items():
    ax.plot(disp["din_instante"], disp["gt"], lw=1.8, label=name)
ax.set_title(f"Termica horaria no caso critico\\n{stress_name}")
ax.set_ylabel("MW")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

ax = axes[1]
for name, disp in tx_dispatch.items():
    ax.plot(disp["din_instante"], disp["CMO"], lw=1.8, label=name)
ax.set_title("CMO horario com diferentes restricoes de transmissao")
ax.set_ylabel("R$/MWh")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

fig.tight_layout(); plt.show()
"""))

    cells.append(new_markdown_cell("## 6. Leitura final do estudo"))

    cells.append(new_code_cell("""base_best = df_base_daily.loc[df_base_daily["gt_energy_GWh"].idxmin()]
base_worst = df_base_daily.loc[df_base_daily["gt_energy_GWh"].idxmax()]
perm_best = df_perm_daily.loc[df_perm_daily["gt_energy_GWh"].idxmin()]
perm_worst = df_perm_daily.loc[df_perm_daily["gt_energy_GWh"].idxmax()]
tx_worst = df_tx.loc[df_tx["gt_energy_GWh"].idxmax()]

print("Principais leituras:")
print(f"  1. Entre os dias representativos, o menor uso termico foi em {base_best['label']} ({base_best['gt_energy_GWh']:.2f} GWh).")
print(f"  2. O maior uso termico foi em {base_worst['label']} ({base_worst['gt_energy_GWh']:.2f} GWh).")
print(f"  3. A melhor permutacao foi '{perm_best['scenario']}' ({perm_best['gt_energy_GWh']:.2f} GWh termicos).")
print(f"  4. A pior permutacao foi '{perm_worst['scenario']}' ({perm_worst['gt_energy_GWh']:.2f} GWh termicos).")
print(f"  5. Na transmissao, o cenario mais restritivo levou a {tx_worst['gt_energy_GWh']:.2f} GWh termicos e CMO medio de {tx_worst['CMO_mean']:.1f} R$/MWh.")
print()
print("Conclusao curta:")
print("  Quando a oferta solar/hidrica e deslocada para perfis menos aderentes a demanda, a termica cresce e o CMO sobe quase linearmente com gt.")
print("  Restricoes de transmissao agravam esse efeito porque reduzem a energia renovavel efetivamente entregue ao centro de carga.")
"""))

    nb = new_notebook(cells=cells, metadata=meta)
    nbformat.write(nb, NB5_PATH)


if __name__ == "__main__":
    update_notebook_4()
    create_notebook_5()

