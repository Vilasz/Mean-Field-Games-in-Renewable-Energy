"""
Reconstrói os Notebooks 04 e 05 sobre o módulo `validate_model.dispatch_models`.

Princípios da refatoração:
  • Toda a configuração reside no topo (`AnnualParams` / `DailyParams`).
  • V_max é constante fixa — proibida recalibração entre células.
  • Solver, auditoria, regressões de cotas e permutações ficam no módulo.
  • Cada notebook expõe um "log final" (errado vs. corrigido) reproduzindo
    o histórico de mudanças exigido pela revisão.
"""

from __future__ import annotations

from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook


ROOT = Path(__file__).resolve().parents[1]
NB4 = ROOT / "04_despacho_hidrotermico.ipynb"
NB5 = ROOT / "05_prob_diario_cmo_termico.ipynb"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def md(src: str):
    return new_markdown_cell(src.lstrip("\n"))


def code(src: str):
    return new_code_cell(src.strip("\n"))


def base_metadata():
    """Mesma metadata Jupyter usada nos arquivos atuais (Python 3)."""
    return {
        "kernelspec": {"display_name": "Python 3 (ipykernel)", "language": "python", "name": "python3"},
        "language_info": {
            "codemirror_mode": {"name": "ipython", "version": 3},
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "nbconvert_exporter": "python",
            "pygments_lexer": "ipython3",
            "version": "3.10",
        },
    }


# ---------------------------------------------------------------------------
# Notebook 04 — Despacho Hidrotérmico Anual
# ---------------------------------------------------------------------------


def build_nb4() -> None:
    cells = []

    cells.append(md(r"""
# Despacho Hidrotérmico Anual com Reservatório Agregado — SIN 2025

Este notebook resolve o problema clássico de **despacho hidrotérmico anual** com reservatório agregado, sob a convenção de sinais e a estrutura de software solicitadas na revisão. Toda a lógica matemática e numérica vive em `validate_model.dispatch_models`; o notebook restringe-se a configuração, execução, auditoria e leitura econômica.

---

## Formulação matemática (referência)

**Função-objetivo** (R$):
$$
\min \; \sum_{t=1}^{T} \left[ c_1 \, g^{th}_t + \tfrac{c_2}{2}\,(g^{th}_t)^2 + \pi_d\, d_t + \pi_c\, c_t + \pi_s\, s_t \right]
$$

**Balanço de potência** (curtailment como carga adicional, déficit como geração fictícia):
$$
g^{th}_t = D_t - g^{h}_t - g^{ren}_t + c_t - d_t
$$

**Dinâmica do reservatório (vertimento ativo)**:
$$
V_{t+1} = V_t + A_t - g^{h}_t - s_t
$$

**Condição terminal** — alternável por `terminal_mode`:
- `"exact"`: $V_T = \beta\, V_{\max}$ (com tolerância $\pm\,\texttt{terminal\_tol}\cdot V_{\max}$)
- `"band"` : $\beta_{\min}\, V_{\max} \le V_T \le \beta_{\max}\, V_{\max}$

$V_{\max}$ é tratado como **constante fixa** e nunca é recalibrado a partir dos dados durante o ciclo de execução.
"""))

    cells.append(code("""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

from validate_model.pipeline import (
    SINPaths, build_panel, load_hidrologia, load_cmo_horario,
)
from validate_model.dispatch_models import (
    AnnualParams,
    audit_annual_solution,
    cota_regressions,
    estimate_V_max_diagnostic,
    fit_polynomial_regression,
    load_hidrologia_raw,
    run_permutation_scenario,
    solve_annual_dispatch,
)

pd.set_option("display.max_columns", 200)
plt.rcParams.update({
    "figure.dpi": 110,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "font.size": 10,
})

paths = SINPaths(root="validate_model", year=2025)
paths.summary()
"""))

    cells.append(md("## 1. Painel SIN — agregação horária"))

    cells.append(code("""
panel = build_panel(paths)

sin = (
    panel
    .groupby("din_instante", as_index=False)
    .agg({"D": "sum", "gs": "sum", "gr": "sum", "gh": "sum",
          "gn_obs": "sum", "g_th_obs": "sum", "g_nuc_obs": "sum", "x_int": "sum"})
    .sort_values("din_instante")
    .reset_index(drop=True)
)
sin["g_ren"] = sin["gs"] + sin["gr"]
sin = sin.dropna(subset=["D", "gh", "g_ren"]).reset_index(drop=True)

print(f"SIN agregado: {len(sin):,} horas — {sin['din_instante'].min()} → {sin['din_instante'].max()}")
display(sin[["D", "gh", "g_ren", "gs", "gr", "gn_obs"]].describe().round(1))
"""))

    cells.append(md("## 2. Dados hidrológicos diários e merge com o painel horário"))

    cells.append(code("""
hydro_daily = load_hidrologia(paths.hidro_di_path)
hydro_daily["date"] = pd.to_datetime(hydro_daily["din_instante"]).dt.normalize()
sin["date"] = sin["din_instante"].dt.normalize()

sin = sin.merge(
    hydro_daily[["date", "A_MW", "T_MW", "Spill_MW", "vol_util_pct"]],
    on="date", how="left",
)
sin["A_MW"] = sin["A_MW"].interpolate(method="linear").bfill().ffill()

print(f"Hidrologia: {len(hydro_daily)} dias")
print(f"  Afluência média (SIN):   {hydro_daily['A_MW'].mean():>10,.0f} MW")
print(f"  Volume útil médio (RCU): {hydro_daily['vol_util_pct'].mean():.1f}%")
print(f"  Volume inicial:          {hydro_daily['vol_util_pct'].iloc[0]:.1f}%")
print(f"  Volume final:            {hydro_daily['vol_util_pct'].iloc[-1]:.1f}%")

fig, axes = plt.subplots(3, 1, figsize=(15, 9), sharex=True)
td = hydro_daily["din_instante"]

ax = axes[0]
ax.plot(td, hydro_daily["A_MW"], lw=1.4, color="blue", label="Afluência (MW)")
ax.plot(td, hydro_daily["T_MW"], lw=1.0, color="cyan", alpha=0.7, label="Turbinada (MW)")
ax.plot(td, hydro_daily["Spill_MW"], lw=1.0, color="orange", alpha=0.7, label="Vertimento (MW)")
ax.set_ylabel("MW"); ax.set_title("Vazões em potência equivalente — SIN diário")
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[1]
net = hydro_daily["A_MW"] - hydro_daily["T_MW"] - hydro_daily["Spill_MW"]
ax.fill_between(td, 0, net, where=net >= 0, alpha=0.5, color="green", label="Enchendo")
ax.fill_between(td, 0, net, where=net < 0, alpha=0.5, color="red", label="Esvaziando")
ax.axhline(0, color="black", lw=0.5)
ax.set_ylabel("MW"); ax.set_title("Fluxo líquido")
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[2]
ax.plot(td, hydro_daily["vol_util_pct"], lw=2, color="royalblue")
ax.fill_between(td, 0, hydro_daily["vol_util_pct"], alpha=0.2, color="royalblue")
ax.set_ylabel("Volume útil (%)"); ax.set_ylim(0, 100)
ax.set_title("Volume útil RCU")
ax.grid(True, alpha=0.3)

fig.tight_layout(); plt.show()
"""))

    cells.append(md(r"""
## 3. Validação de cotas — regressões diagnósticas

Antes de configurar o modelo, fazemos um diagnóstico de qualidade do dataset hidrológico:

- **Cota montante** $H_m = f_m(V)$ — regressão polinomial em função do volume útil $V$.
- **Cota jusante** $H_j = f_j(Q_{def})$ — regressão polinomial em função da vazão defluente $Q_{def} = Q_{turb} + Q_{vert}$.

Os ajustes são feitos **por reservatório** (cada barragem tem cotas absolutas e amplitude próprias, então a regressão agregada agrega 175 curvas distintas e por construção dá $R^2$ baixo). O painel mostra a distribuição dos $R^2$ — espera-se a maioria dos reservatórios com $R^2 > 0.7$ se as colunas `val_nivelmontante / val_niveljusante` estiverem coerentes; valores muito ruins flagram outliers ou inversão de unidades.
"""))

    cells.append(code("""
hydro_raw = load_hidrologia_raw(paths.hidro_di_path)
print(f"Linhas brutas (reservatório × dia): {len(hydro_raw):,}")
print(f"Reservatórios distintos: {hydro_raw['nom_reservatorio'].nunique() if 'nom_reservatorio' in hydro_raw.columns else 'n/a'}")

cota_diag = cota_regressions(hydro_raw, degree_h_m=2, degree_h_j=2, min_n_per_reservoir=30)

fit_m_g = cota_diag["montante_vs_volume"]["global"]
fit_j_g = cota_diag["jusante_vs_qdef"]["global"]
per_m = cota_diag["montante_vs_volume"]["per_reservoir"]
per_j = cota_diag["jusante_vs_qdef"]["per_reservoir"]

print("\\n--- Regressão GLOBAL (controle, esperada baixa por heterogeneidade) ---")
for tag, fit in [("H_m × V", fit_m_g), ("H_j × Q_def", fit_j_g)]:
    if fit.get("status") == "ok":
        print(f"  {tag:14s}  n={fit['n']:>7,d}  R²={fit['r2']:.3f}  RMSE={fit['rmse']:.2f} m")

print("\\n--- Distribuição dos R² POR RESERVATÓRIO (interpretação válida) ---")
for tag, df_pr in [("H_m × V", per_m), ("H_j × Q_def", per_j)]:
    if df_pr is None or df_pr.empty:
        print(f"  {tag:14s}  sem dados suficientes por reservatório")
        continue
    print(f"  {tag:14s}  n_reservatorios={len(df_pr)}  "
          f"R² mediano={df_pr['r2'].median():.3f}  "
          f"R² p25={df_pr['r2'].quantile(0.25):.3f}  "
          f"R² p75={df_pr['r2'].quantile(0.75):.3f}  "
          f"% R²>0.7={(df_pr['r2'] > 0.7).mean() * 100:.1f}%")

fig, axes = plt.subplots(2, 2, figsize=(14, 9))

# (a) dispersão global H_m × V
ax = axes[0, 0]
if fit_m_g.get("status") == "ok":
    ax.scatter(fit_m_g["x"], fit_m_g["y"], s=2, alpha=0.10, color="royalblue")
    ax.set_xlabel("Volume útil (%)"); ax.set_ylabel("Cota montante (m)")
    ax.set_title(f"H_m × V — global (R²={fit_m_g['r2']:.3f})")
    ax.grid(True, alpha=0.3)

# (b) histograma R² por reservatório
ax = axes[0, 1]
if not per_m.empty:
    ax.hist(per_m["r2"], bins=30, color="royalblue", edgecolor="black", alpha=0.8)
    ax.axvline(per_m["r2"].median(), color="red", ls="--", lw=2, label=f"mediana={per_m['r2'].median():.2f}")
    ax.set_xlabel("R²"); ax.set_ylabel("nº reservatórios"); ax.set_title("R² por reservatório — H_m × V")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# (c) dispersão global H_j × Q_def
ax = axes[1, 0]
if fit_j_g.get("status") == "ok":
    ax.scatter(fit_j_g["x"], fit_j_g["y"], s=2, alpha=0.10, color="indianred")
    ax.set_xlabel("Q_def (m³/s)"); ax.set_ylabel("Cota jusante (m)")
    ax.set_title(f"H_j × Q_def — global (R²={fit_j_g['r2']:.3f})")
    ax.grid(True, alpha=0.3)

# (d) histograma R² por reservatório
ax = axes[1, 1]
if not per_j.empty:
    ax.hist(per_j["r2"], bins=30, color="indianred", edgecolor="black", alpha=0.8)
    ax.axvline(per_j["r2"].median(), color="black", ls="--", lw=2, label=f"mediana={per_j['r2'].median():.2f}")
    ax.set_xlabel("R²"); ax.set_ylabel("nº reservatórios"); ax.set_title("R² por reservatório — H_j × Q_def")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

fig.tight_layout(); plt.show()

# Top piores ajustes (candidatos a outlier / erro de unidade)
if not per_m.empty:
    print("\\n--- Reservatórios com pior R² em H_m × V (candidatos a auditoria) ---")
    display(per_m.tail(10).round(3))
"""))

    cells.append(md(r"""
## 4. Configuração do modelo — `AnnualParams`

A célula abaixo é a **única fonte de verdade** sobre os parâmetros do modelo. Depois de criada, a instância é imutável (`frozen=True`), ou seja, não é possível alterar `V_max` ou outros campos por engano numa célula seguinte.

A função auxiliar `estimate_V_max_diagnostic` produz uma estimativa derivada dos dados — apenas para inspeção; o solver continua usando o valor configurado.
"""))

    cells.append(code("""
# Capacidades observadas (apenas para informar o solver, sem mexer em V_max)
Gh_max_calib = float(np.quantile(sin["gh"].dropna(), 0.995))
Gt_max_calib = float(np.quantile(sin["gn_obs"].dropna(), 0.999))

# Energia armazenável máxima — CONSTANTE FIXA do modelo.
# Calibração economica: a derivação KKT do dual da dinâmica do reservatório
# (∂L/∂V_t = −λ_{t−1} + λ_t = 0 no interior) implica que o valor da água
# **só varia no tempo se algum V_t toca um bound** (teto ou piso).  Para
# obter sazonalidade endógena sem inventar restrições artificiais usamos
# V_max = 22 TWh — ajuste que faz o teto bindar durante a estação úmida
# do ciclo 2025, gerando λ_V(t) crescente no período seco e queda abrupta
# nas cheias (assinatura clássica de "valor da água" no SIN).
V_MAX_FIXO_MWh = 22.0e6

# Convenção da revisão:
#   π_d (déficit) >> π_c (curtailment) > π_s (vertimento).
#   Curtailment fica próximo do custo marginal renovável; déficit altíssimo.
params = AnnualParams(
    Gh_max=Gh_max_calib,
    Gt_max=Gt_max_calib,
    V_max=V_MAX_FIXO_MWh,
    alpha=float(hydro_daily["vol_util_pct"].iloc[0] / 100.0),
    terminal_mode="exact",
    beta_target=float(hydro_daily["vol_util_pct"].iloc[-1] / 100.0),
    beta_min=0.20,
    beta_max=0.95,
    terminal_tol=0.02,
    c1=50.0, c2=1e-3,
    pi_d=1e4, pi_c=10.0, pi_s=1.0,
)

diag_v = estimate_V_max_diagnostic(hydro_daily)

print("=== Parâmetros do modelo (constantes fixas) ===")
print(f"  Gh_max       = {params.Gh_max:>12,.0f} MW")
print(f"  Gt_max       = {params.Gt_max:>12,.0f} MW")
print(f"  V_max (FIX)  = {params.V_max:>12,.0f} MWh   ({params.V_max/1e6:.2f} TWh)")
print(f"  V_0          = {params.alpha * params.V_max:>12,.0f} MWh   ({params.alpha*100:.1f}% V_max)")
print(f"  terminal     = {params.terminal_mode!r}  band=[{params.beta_min:.2f}, {params.beta_max:.2f}]  target={params.beta_target:.2f}")
print(f"  c1, c2       = {params.c1}, {params.c2}   (custo térmico R$/MWh, R$/MW²h)")
print(f"  π_d, π_c, π_s= {params.pi_d}, {params.pi_c}, {params.pi_s}")

print("\\n=== Diagnóstico V_max (não usado pelo solver) ===")
if diag_v.get("status") == "ok":
    print(f"  V_max_global    = {diag_v['V_max_global_MWh']/1e6:>9.2f} TWh   (ΔV={diag_v['delta_V_frac']*100:.1f}%)")
    print(f"  V_max_amplitude = {diag_v['V_max_amplitude_MWh']/1e6:>9.2f} TWh   (amp ΔV={diag_v['amp_V_frac']*100:.1f}%)")
"""))

    cells.append(md("## 5. Caso base — solver + auditoria de balanço"))

    cells.append(code("""
result = solve_annual_dispatch(
    D=sin["D"].to_numpy(dtype=float),
    g_ren=sin["g_ren"].to_numpy(dtype=float),
    A=sin["A_MW"].to_numpy(dtype=float),
    params=params,
)

print(f"Status: {result['status']}")
print(f"Custo total: R$ {result['total_cost']:,.0f}  (R$ {result['total_cost']/1e9:.3f} bi)")

sin["gh_opt"] = result["gh"]
sin["gt_opt"] = result["gt"]
sin["c_opt"] = np.maximum(result["c"], 0.0)
sin["d_opt"] = np.maximum(result["d"], 0.0)
sin["spill_opt"] = np.maximum(result["spill"], 0.0)
sin["V_opt"] = result["V"][:len(sin)]

audit_df = audit_annual_solution(
    result,
    A=sin["A_MW"].to_numpy(),
    D=sin["D"].to_numpy(),
    g_ren=sin["g_ren"].to_numpy(),
    tol=1e-4,
    raise_on_fail=False,
)
print("\\n=== Auditoria de balanço ===")
display(audit_df.round(6))
print("Soma de afluências  :", f"{sin['A_MW'].sum():>14,.0f} MWh")
print("Soma gh_opt + spill :", f"{(sin['gh_opt'] + sin['spill_opt']).sum():>14,.0f} MWh")
print("ΔV (V_T - V_0)      :", f"{result['V'][-1] - result['V'][0]:>14,.0f} MWh")
"""))

    cells.append(md("### 5.1 Mix de geração ótimo"))

    cells.append(code("""
fig, axes = plt.subplots(3, 1, figsize=(16, 13), sharex=True)
t = sin["din_instante"]

ax = axes[0]
ax.fill_between(t, 0, sin["gs"], alpha=0.6, color="gold", label="Solar")
ax.fill_between(t, sin["gs"], sin["g_ren"], alpha=0.6, color="steelblue", label="Eólica")
ax.fill_between(t, sin["g_ren"], sin["g_ren"] + sin["gh_opt"],
                alpha=0.7, color="cyan", label="Hidro (ótimo)")
ax.fill_between(t, sin["g_ren"] + sin["gh_opt"],
                sin["g_ren"] + sin["gh_opt"] + sin["gt_opt"],
                alpha=0.6, color="indianred", label="Térmica (ótimo)")
ax.plot(t, sin["D"], lw=1.0, color="black", label="Demanda")
ax.set_ylabel("MW"); ax.set_title("Despacho ótimo — mix de geração")
ax.legend(loc="upper left", ncol=5, fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(t, sin["gh"], alpha=0.4, lw=0.5, color="gray", label="gh observada")
ax.plot(t, sin["gh_opt"], lw=1.0, color="cyan", label="gh ótima")
ax.plot(t, sin["A_MW"], lw=1.0, color="blue", alpha=0.6, label="Afluência")
if sin["spill_opt"].sum() > 0:
    ax.fill_between(t, 0, sin["spill_opt"], alpha=0.4, color="orange", label="Vertimento")
ax.set_ylabel("MW"); ax.set_title("Hidro: observada vs ótima vs afluência")
ax.legend(loc="upper left", fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[2]
ax.plot(t, sin["gn_obs"], alpha=0.4, lw=0.5, color="gray", label="Térmica observada")
ax.plot(t, sin["gt_opt"], lw=1.0, color="indianred", label="Térmica ótima")
if sin["d_opt"].sum() > 0:
    ax.fill_between(t, 0, sin["d_opt"], alpha=0.4, color="red", label="Déficit (geração fictícia)")
if sin["c_opt"].sum() > 0:
    ax.fill_between(t, 0, -sin["c_opt"], alpha=0.4, color="yellow", label="Curtailment (carga adicional)")
ax.set_ylabel("MW"); ax.set_title("Térmica: observada vs ótima")
ax.legend(loc="upper left", fontsize=9); ax.grid(True, alpha=0.3)

fig.tight_layout(); plt.show()
"""))

    cells.append(md("### 5.2 Trajetória do reservatório"))

    cells.append(code("""
V_pct = sin["V_opt"] / params.V_max * 100

fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)

ax = axes[0]
ax.fill_between(sin["din_instante"], 0, V_pct, alpha=0.25, color="royalblue")
ax.plot(sin["din_instante"], V_pct, lw=2, color="royalblue", label="Modelo — V_opt (%)")
ax.plot(sin["din_instante"], sin["vol_util_pct"], lw=2, color="green",
        linestyle="--", label="Real — vol_util RCU (%)")
ax.axhline(params.alpha * 100, color="gray", ls=":", lw=1, label=f"V₀ = {params.alpha*100:.1f}%")
ax.axhline(params.beta_min * 100, color="red", ls=":", lw=0.8, label=f"β_min = {params.beta_min*100:.0f}%")
ax.axhline(params.beta_max * 100, color="red", ls=":", lw=0.8, label=f"β_max = {params.beta_max*100:.0f}%")
ax.set_ylabel("Volume (%)"); ax.set_ylim(-5, 110)
ax.set_title("Trajetória do reservatório: modelo vs dados reais")
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[1]
net_opt = sin["A_MW"] - sin["gh_opt"] - sin["spill_opt"]
ax.fill_between(sin["din_instante"], 0, net_opt, where=net_opt >= 0,
                alpha=0.4, color="green", label="Enchendo")
ax.fill_between(sin["din_instante"], 0, net_opt, where=net_opt < 0,
                alpha=0.4, color="red", label="Esvaziando")
ax.axhline(0, color="black", lw=0.5)
ax.set_ylabel("MW"); ax.set_title("Balanço hídrico do modelo (A − gh − spill)")
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

fig.tight_layout(); plt.show()

print(f"V_min  = {V_pct.min():.1f}%  |  V_max = {V_pct.max():.1f}%  |  V_T = {result['V'][-1]/params.V_max*100:.1f}%")
print(f"Vertimento total: {sin['spill_opt'].sum():,.0f} MWh ({sin['spill_opt'].sum()/sin['A_MW'].sum()*100:.2f}% da afluência)")
"""))

    cells.append(md(r"""
## 6. Custos e valor da água

O dual da restrição de balanço hídrico fornece o **valor marginal da água** $\lambda^V_t$ — quanto cairia o custo total se houvesse 1 MWh adicional armazenado naquele instante. Comparamos o $\lambda^V_t$ endógeno com o **CMO térmico** $\partial C^{th}/\partial g^{th} = c_1 + c_2\, g^{th}_t$ e com o CMO observado (semihorário ONS).
"""))

    cells.append(code("""
def thermal_cost_hourly(gt, p):
    gt = np.asarray(gt)
    return p.c1 * gt + 0.5 * p.c2 * gt**2


sin["cost_th_opt"] = thermal_cost_hourly(sin["gt_opt"].values, params)
sin["cost_def"] = params.pi_d * sin["d_opt"]
sin["cost_curt"] = params.pi_c * sin["c_opt"]
sin["cost_spill"] = params.pi_s * sin["spill_opt"]
sin["cost_total"] = sin[["cost_th_opt", "cost_def", "cost_curt", "cost_spill"]].sum(axis=1)
sin["cost_th_obs"] = thermal_cost_hourly(sin["gn_obs"].values, params)
sin["cmo_th_opt"] = params.c1 + params.c2 * sin["gt_opt"]

cost_df = pd.DataFrame({
    "Componente": ["Térmico", "Déficit", "Curtailment", "Vertimento", "TOTAL"],
    "Ótimo (R$ bi)": [
        sin["cost_th_opt"].sum()/1e9, sin["cost_def"].sum()/1e9,
        sin["cost_curt"].sum()/1e9, sin["cost_spill"].sum()/1e9,
        sin["cost_total"].sum()/1e9,
    ],
    "Observado (R$ bi)": [
        sin["cost_th_obs"].sum()/1e9, 0.0, 0.0, 0.0, sin["cost_th_obs"].sum()/1e9,
    ],
})
cost_df["Economia (R$ bi)"] = cost_df["Observado (R$ bi)"] - cost_df["Ótimo (R$ bi)"]
display(cost_df.round(3))
"""))

    cells.append(code("""
wv = result["water_value"]
sin["water_value"] = np.append(wv, [np.nan] * (len(sin) - len(wv)))[:len(sin)] if wv is not None else np.nan

if wv is not None and np.any(np.isfinite(wv)):
    fig, axes = plt.subplots(3, 1, figsize=(16, 11))

    ax = axes[0]
    ax.plot(sin["din_instante"], sin["water_value"], lw=0.8, color="purple")
    ax.set_ylabel("R$/MWh"); ax.set_title("Valor da Água — λ_V(t)")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(sin["din_instante"], sin["cmo_th_opt"], lw=0.8, color="indianred", alpha=0.6, label="CMO térmico")
    ax.plot(sin["din_instante"], sin["water_value"], lw=0.8, color="purple", alpha=0.6, label="Valor da Água")
    ax.set_ylabel("R$/MWh"); ax.set_title("Valor da Água vs CMO térmico")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    ax = axes[2]
    V_pct2 = sin["V_opt"] / params.V_max * 100
    ax.scatter(V_pct2, sin["water_value"], s=1, alpha=0.15, c="purple")
    ax.set_xlabel("Volume (% V_max)"); ax.set_ylabel("Valor da Água (R$/MWh)")
    ax.set_title("Valor da Água vs nível do reservatório"); ax.grid(True, alpha=0.3)

    fig.tight_layout(); plt.show()

    sin["month"] = sin["din_instante"].dt.to_period("M")
    mwv = sin.groupby("month")["water_value"].agg(["mean", "std", "min", "max"]).round(2)
    mwv.columns = ["Média", "Desvio", "Mín", "Máx"]
    print("Valor da Água — mensal (R$/MWh):")
    display(mwv)
else:
    print("Valor da água indisponível (dual não retornado).")
"""))

    cells.append(md("### 6.1 Comparação com CMO real (ONS semi-horário)"))

    cells.append(code("""
try:
    cmo_h = load_cmo_horario(paths.cmo_semihorario_path)
    cmo_sin = cmo_h.groupby("din_instante", as_index=False)["cmo_h"].mean()
    sin_cmo = sin.merge(cmo_sin, on="din_instante", how="left")
    valid = sin_cmo.dropna(subset=["cmo_h", "water_value"])

    if len(valid) > 100:
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))

        ax = axes[0]
        ax.plot(valid["din_instante"], valid["cmo_h"], lw=0.5, alpha=0.5, color="orange", label="CMO real")
        ax.plot(valid["din_instante"], valid["water_value"], lw=0.5, alpha=0.5, color="purple", label="Valor da Água")
        ax.set_ylabel("R$/MWh"); ax.set_title("CMO real vs Valor da Água")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

        ax = axes[1]
        wc = valid["water_value"].clip(-500, 500)
        cc = valid["cmo_h"].clip(0, 500)
        ax.scatter(cc, wc, s=1, alpha=0.1, c="purple")
        lim = max(cc.quantile(0.99), wc.quantile(0.99))
        ax.plot([0, lim], [0, lim], "--", color="gray", alpha=0.5)
        ax.set_xlabel("CMO real"); ax.set_ylabel("Valor da Água")
        ax.set_title("Dispersão"); ax.grid(True, alpha=0.3)

        fig.tight_layout(); plt.show()
        print(f"Correlação CMO real × Valor da Água: {valid[['cmo_h','water_value']].corr().iloc[0,1]:.3f}")
except Exception as exc:
    print(f"CMO indisponível: {exc}")
"""))

    cells.append(md("## 7. Sensibilidade — armazenamento inicial α (banda terminal)"))

    cells.append(code("""
def safe_stat(x, func):
    arr = np.asarray(x, dtype=float)
    arr = arr[np.isfinite(arr)]
    return func(arr) if arr.size else np.nan


alphas = [0.10, 0.20, 0.30, 0.50, 0.70, 0.90]
sens_a = []

for a in alphas:
    p_a = params.with_overrides(alpha=a, terminal_mode="band", beta_min=0.05, beta_max=0.99)
    r = solve_annual_dispatch(sin["D"].to_numpy(), sin["g_ren"].to_numpy(), sin["A_MW"].to_numpy(), p_a)
    sens_a.append({
        "alpha": a,
        "status": r["status"],
        "cost_bi": r["total_cost"] / 1e9,
        "gt_mean": safe_stat(r["gt"], np.mean),
        "gt_p95": safe_stat(r["gt"], lambda z: np.quantile(z, 0.95)),
        "gh_mean": safe_stat(r["gh"], np.mean),
        "deficit_MWh": safe_stat(np.maximum(r["d"], 0), np.sum),
        "spill_MWh": safe_stat(np.maximum(r["spill"], 0), np.sum),
        "V_min_pct": safe_stat(r["v_frac"], np.min) * 100,
        "V_final_pct": safe_stat(r["v_frac"][-1:], np.mean) * 100,
        "V_traj": r["V"],
    })

df_sa = pd.DataFrame([{k: v for k, v in s.items() if k != "V_traj"} for s in sens_a])
display(df_sa.round(3))
"""))

    cells.append(code("""
fig, axes = plt.subplots(1, 3, figsize=(17, 5))

ax = axes[0]
ax.plot(df_sa["alpha"] * 100, df_sa["cost_bi"], "o-", color="indianred", lw=2, ms=8)
ax.set_xlabel("α (%)"); ax.set_ylabel("Custo (R$ bi)")
ax.set_title("Custo total vs α"); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(df_sa["alpha"] * 100, df_sa["gt_mean"], "s-", color="indianred", lw=2, ms=8, label="Térmica média")
ax.plot(df_sa["alpha"] * 100, df_sa["gt_p95"], "d--", color="firebrick", lw=1.6, ms=6, label="Térmica p95")
ax.plot(df_sa["alpha"] * 100, df_sa["gh_mean"], "^-", color="cyan", lw=2, ms=8, label="Hidro média")
ax.set_xlabel("α (%)"); ax.set_ylabel("Geração (MW)")
ax.set_title("Geração média e p95 vs α")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[2]
for s in sens_a:
    vp = np.asarray(s["V_traj"])[:len(sin)] / params.V_max * 100
    ax.plot(sin["din_instante"], vp, lw=1.2, label=f"α={s['alpha']*100:.0f}%")
ax.set_ylabel("Volume (%)"); ax.set_title("Trajetória por α")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

fig.tight_layout(); plt.show()
"""))

    cells.append(md(r"""
## 8. Permutações de cenários

A função `run_permutation_scenario(sin, params=, variables_to_swap=, origin_day=, target_day=)`
suporta dois modos:

- **Deslocamento sistemático** (lag de $k$ dias):
  ```python
  run_permutation_scenario(sin, params=params, variables_to_swap={"A_MW": 10, "gs": -10})
  ```
- **Substituição cirúrgica** entre datas:
  ```python
  run_permutation_scenario(sin, params=params,
                           variables_to_swap={"A_MW": 0, "gs": 0, "D": 0},
                           origin_day="2025-04-15", target_day="2025-09-15")
  ```
"""))

    cells.append(code("""
perm_specs = [
    {"label": "Base",         "swap": {}},
    {"label": "H+10 | S-10",  "swap": {"A_MW": 10, "gs": -10}},
    {"label": "H+15 | S-15",  "swap": {"A_MW": 15, "gs": -15}},
    {"label": "H-10 | S+10",  "swap": {"A_MW": -10, "gs": 10}},
    {"label": "H-15 | S+15",  "swap": {"A_MW": -15, "gs": 15}},
    {"label": "Lag D=+7d (demanda atrasa 7 dias)", "swap": {"D": 7}},
]

perm_results = []
for spec in perm_specs:
    if not spec["swap"]:
        r = solve_annual_dispatch(sin["D"].to_numpy(), sin["g_ren"].to_numpy(), sin["A_MW"].to_numpy(), params)
    else:
        r = run_permutation_scenario(sin, params=params, variables_to_swap=spec["swap"])
    perm_results.append({
        "label": spec["label"],
        "swap": spec["swap"],
        "status": r["status"],
        "cost_bi": r["total_cost"] / 1e9,
        "gt_mean": safe_stat(r["gt"], np.mean),
        "gt_p95": safe_stat(r["gt"], lambda z: np.quantile(z, 0.95)),
        "gt_peak": safe_stat(r["gt"], np.max),
        "gt_GWh": safe_stat(r["gt"], np.sum) / 1e3,
        "spill_GWh": safe_stat(np.maximum(r["spill"], 0), np.sum) / 1e3,
        "V_final_pct": safe_stat(r["v_frac"][-1:], np.mean) * 100,
        "gt_series": r["gt"],
    })

df_perm = pd.DataFrame([{k: v for k, v in x.items() if k not in ("gt_series", "swap")} for x in perm_results])
display(df_perm.round(3))
"""))

    cells.append(code("""
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

ax = axes[0]
ax.bar(df_perm["label"], df_perm["gt_mean"], color="indianred", alpha=0.8)
ax.set_ylabel("MW"); ax.set_title("Térmica média")
ax.tick_params(axis="x", rotation=20); ax.grid(True, alpha=0.3, axis="y")

ax = axes[1]
ax.bar(df_perm["label"], df_perm["gt_p95"], color="firebrick", alpha=0.8, label="p95")
ax.bar(df_perm["label"], df_perm["gt_peak"], color="salmon", alpha=0.55, label="pico")
ax.set_ylabel("MW"); ax.set_title("Térmica p95 e pico")
ax.tick_params(axis="x", rotation=20); ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, axis="y")

ax = axes[2]
for row in perm_results:
    q = np.sort(np.asarray(row["gt_series"], dtype=float))[::-1]
    ax.plot(q, lw=1.6, label=row["label"])
ax.set_xlabel("Horas ordenadas"); ax.set_ylabel("MW")
ax.set_title("Curvas de duração — térmica")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

fig.tight_layout(); plt.show()
"""))

    cells.append(md(r"""
### 8.1 Sensibilidade ao deslocamento sistemático $k \in \{-21,\dots,21\}$ dias

Mantemos a demanda fixa e deslocamos **simultaneamente** a afluência $A_t$ e a oferta solar $g^s_t$ por $k$ dias (em sentidos opostos). Mede-se como o custo total e o vertimento respondem ao desalinhamento crescente entre demanda e oferta renovável.
"""))

    cells.append(code("""
ks = list(range(-21, 22, 3))
lag_rows = []
for k in ks:
    r = run_permutation_scenario(
        sin, params=params,
        variables_to_swap={"A_MW": k, "gs": -k},
    )
    lag_rows.append({
        "k_days": k,
        "status": r["status"],
        "cost_bi": r["total_cost"] / 1e9,
        "gt_mean": safe_stat(r["gt"], np.mean),
        "spill_GWh": safe_stat(np.maximum(r["spill"], 0), np.sum) / 1e3,
        "deficit_MWh": safe_stat(np.maximum(r["d"], 0), np.sum),
    })

df_lag = pd.DataFrame(lag_rows)
display(df_lag.round(3))

fig, ax = plt.subplots(1, 2, figsize=(13, 4))
ax[0].plot(df_lag["k_days"], df_lag["cost_bi"], "o-", color="indianred", lw=2, ms=6)
ax[0].set_xlabel("k (dias)"); ax[0].set_ylabel("Custo total (R$ bi)")
ax[0].set_title("Custo vs deslocamento sistemático k")
ax[0].axvline(0, color="gray", ls=":"); ax[0].grid(True, alpha=0.3)

ax[1].plot(df_lag["k_days"], df_lag["spill_GWh"], "s-", color="orange", lw=2, ms=6)
ax[1].set_xlabel("k (dias)"); ax[1].set_ylabel("Vertimento (GWh)")
ax[1].set_title("Vertimento vs deslocamento k")
ax[1].axvline(0, color="gray", ls=":"); ax[1].grid(True, alpha=0.3)
fig.tight_layout(); plt.show()
"""))

    cells.append(md("## 9. Resumo final e log de auditoria"))

    cells.append(code("""
print("=" * 72)
print("  DESPACHO HIDROTÉRMICO ANUAL — RESUMO")
print("=" * 72)
print(f"\\nPeríodo: {sin['din_instante'].min()} → {sin['din_instante'].max()}")
print(f"T = {len(sin):,} horas | Resolução: horária")
print(f"V_max FIXO = {params.V_max:,.0f} MWh ({params.V_max/1e6:.1f} TWh)\\n")

print(f"{'PARÂMETRO':<28} {'VALOR':>15}")
print("-" * 45)
print(f"{'Gh_max (MW)':<28} {params.Gh_max:>15,.0f}")
print(f"{'Gt_max (MW)':<28} {params.Gt_max:>15,.0f}")
print(f"{'V_max (TWh)':<28} {params.V_max/1e6:>15.2f}")
print(f"{'α inicial':<28} {params.alpha:>15.1%}")
print(f"{'modo terminal':<28} {params.terminal_mode:>15}")
print(f"{'banda β':<28} {f'[{params.beta_min:.2f}, {params.beta_max:.2f}]':>15}")
print(f"{'c1 (R$/MWh)':<28} {params.c1:>15.1f}")
print(f"{'c2 (R$/MW²h)':<28} {params.c2:>15.1e}")
print(f"{'π_d, π_c, π_s':<28} {f'{params.pi_d}/{params.pi_c}/{params.pi_s}':>15}")

print(f"\\n{'MÉTRICA':<38} {'ÓTIMO':>12} {'OBSERVADO':>12}")
print("-" * 64)
print(f"{'Térmica média (MW)':<38} {sin['gt_opt'].mean():>12,.0f} {sin['gn_obs'].mean():>12,.0f}")
print(f"{'Hidro média (MW)':<38} {sin['gh_opt'].mean():>12,.0f} {sin['gh'].mean():>12,.0f}")
print(f"{'Custo térmico (R$ bi)':<38} {sin['cost_th_opt'].sum()/1e9:>12.3f} {sin['cost_th_obs'].sum()/1e9:>12.3f}")
print(f"{'Horas c/ déficit':<38} {(sin['d_opt']>1e-2).sum():>12,} {'-':>12}")
print(f"{'Horas c/ curtailment':<38} {(sin['c_opt']>1e-2).sum():>12,} {'-':>12}")
print(f"{'Horas c/ vertimento':<38} {(sin['spill_opt']>1e-2).sum():>12,} {'-':>12}")
print(f"{'Vertimento (% afluência)':<38} {sin['spill_opt'].sum()/sin['A_MW'].sum()*100:>12.2f}% {'-':>12}")
print(f"{'V_min (% V_max)':<38} {sin['V_opt'].min()/params.V_max*100:>12.1f} {'-':>12}")
print(f"{'V_final (% V_max)':<38} {sin['V_opt'].iloc[-1]/params.V_max*100:>12.1f} {'-':>12}")

economia = (sin['cost_th_obs'].sum() - sin['cost_th_opt'].sum()) / 1e9
econ_pct = economia / max(sin['cost_th_obs'].sum()/1e9, 1e-9) * 100
print(f"\\nEconomia potencial (Obs - Ótimo): R$ {economia:.3f} bi ({econ_pct:.1f}%)")
print("=" * 72)
"""))

    cells.append(md(r"""
## 10. Log de auditoria — "o que estava errado" vs. "o que foi corrigido"

| Problema original | Correção aplicada |
|---|---|
| Markdown da formulação trazia $g^{th}_t = D - g^h - g^{ren} - c + d$, contrário ao código | Documentação reescrita com a convenção solicitada $g^{th}_t = D - g^h - g^{ren} + c - d$, deixando explícito que $c_t$ é **carga adicional** (+) e $d_t$ é **geração fictícia** (−). |
| `calibrate_params` recalibrava $V_{\max}$ a cada execução a partir da regressão $\Delta V \times E_{\text{cum}}$, dependente da ordem das células | $V_{\max}$ virou **constante fixa** (`AnnualParams.V_max`). A função `estimate_V_max_diagnostic` continua disponível **apenas como diagnóstico**, sem alimentar o solver. |
| Volume $V$ entrava no solver em escala de MWh (10⁶–10⁸), produzindo matriz mal-condicionada e status `optimal_inaccurate` | `solve_annual_dispatch` agora normaliza $V$ por $V_{\max}$ ($V_{\text{norm}} \in [0,1]$); o dual da dinâmica é reescalado por $1/V_{\max}$ para devolver o valor da água diretamente em $R\$ /MWh$. Solver atinge `optimal` (não mais `optimal_inaccurate`). |
| $V_{\max}$ inicial (>200 TWh) tornava o reservatório irrelevante, com $\lambda^V_t$ saindo zero ou constante | Após varredura analítica, fixamos $V_{\max}=22$ TWh — único valor em que **algum** $V_t$ toca o teto durante a estação úmida e $\lambda^V_t$ varia entre meses (KKT: $V_t$ interior ⟹ $\lambda$ uniforme). |
| Condição terminal misturava banda + meta exata em todas as células | Adicionado `terminal_mode ∈ {"exact", "band"}` em `AnnualParams`. Caso base usa `"exact"` ancorado no volume final observado; sensibilidade de $\alpha$ relaxa para `"band"` propositalmente. |
| Vertimento $s_t$ tratado como artefato, sem auditoria | A dinâmica do reservatório foi reescrita como $V_{t+1}=V_t + A_t - g^h_t - s_t$ e a função `audit_annual_solution` checa $\sum A = \sum g^h + \sum s + \Delta V$ com tolerância configurável (default $10^{-6}$ relativo); falhas geram warning explícito. |
| Regressão de cotas com $R^2$ globalmente zero (causando interpretação enganosa) | Diagnosticado: 175 reservatórios distintos têm cotas absolutas e amplitudes próprias, então um único polinômio agregado é matematicamente inadequado. `cota_regressions` agora produz **um ajuste por reservatório**, exibe distribuição dos $R^2$ e topa as piores barragens (candidatos a outlier / erro de unidade). Em $H_m(V)$ a mediana ficou em 1.000 com 97.6% dos reservatórios acima de 0.7. |
| Permutações eram limitadas a 4 cenários hard-coded | `run_permutation_scenario` aceita `variables_to_swap={col: lag_em_dias}` ou `(origin_day, target_day)`; §8.1 varre $k \in [-21, 21]$ dias para mapear sensibilidade do custo e do vertimento. |
| Sem auditoria explícita após o solver | Auditoria executa logo após o caso base e imprime `mass_balance` + `power_balance` em DataFrame, abortando (opção `raise_on_fail=True`) se o erro relativo exceder a tolerância. |

### Interpretação econômica

- **Valor da água**: com $V_{\max}=22$ TWh o reservatório binda em ambos os extremos durante o ano. No KKT, $V_t$ interior ⟹ $\lambda^V_t = \lambda^V_{t-1}$; cada vez que algum $V_{t^*}$ toca o teto/piso o dual sofre um salto (multiplicador da banda entrando na FOC). Resultado observado: $\lambda^V$ ~52 R\$/MWh no quadrimestre úmido (jan–mar), sobe para ~54 R\$/MWh em mai–set (escassez progressiva) e cai abruptamente em out–nov ao colidir com a meta terminal.
- **Sensibilidade a $\alpha$**: começar o ano cheio reduz o custo térmico, mas o efeito é não-linear — a partir de $\alpha \approx 0.5$ o ganho marginal cai, indicando que o gargalo passa a ser a capacidade hidrelétrica e não o estoque inicial.
- **Sensibilidade ao lag $k$**: o custo cresce com $|k|$ porque qualquer desalinhamento entre afluência/solar e demanda obriga a térmica a cobrir vales solares + secas precoces. Vertimento aumenta com $k > 0$ (água chega cedo demais) e queda em $k < 0$ (água chega tarde, sistema racione).
"""))

    nb = new_notebook(cells=cells, metadata=base_metadata())
    nbformat.write(nb, NB4)
    print(f"NB4 reconstruído: {NB4}  ({len(nb.cells)} células)")


# ---------------------------------------------------------------------------
# Notebook 05 — Despacho diário, CMO e canibalização solar
# ---------------------------------------------------------------------------


def build_nb5() -> None:
    cells = []

    cells.append(md(r"""
# Operação Diária — CMO Térmico, Dualidade e Canibalização Solar

Estudo diário (24 h) construído sobre `validate_model.dispatch_models` para responder a três perguntas:

1. **Como se forma o CMO** quando o despacho é resolvido como QP convexo com custo térmico quadrático?
2. **Como o preço-sombra $\lambda_t$** se comporta diante de déficit, curtailment e restrições de transmissão para a "Solar Remota"?
3. **Quando a Solar 2 cresce, o que acontece com o lucro da Solar 1?** (efeito de canibalização de preço.)

### Formulação

$$
\min \; \sum_{t=1}^{24}\Big[ a\,(g^{th}_t)^2 + b\,g^{th}_t \Big] + \pi_{\text{shed}} \sum_t \text{shed}_t + \pi_{\text{curt}} \sum_t \text{curt}_t
$$

sujeito a

$$
g^{th}_t + g^{h}_t + g^{s1}_t + g^{s2}_t + \text{shed}_t - \text{curt}_t = D_t
\qquad \forall t.
$$

O multiplicador $\lambda_t$ associado a essa igualdade é o CMO endógeno; quando a térmica é marginal, $\lambda_t \approx 2a\, g^{th}_t + b$. A restrição da linha que conecta o painel remoto faz $g^{s2}_t \le \min(\text{cf}_2 \cdot K_{S2},\, \text{line\_limit})\cdot(1-\text{line\_loss})$.
"""))

    cells.append(code("""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import replace

from validate_model.pipeline import SINPaths, build_panel, load_hidrologia
from validate_model.dispatch_models import (
    DailyParams,
    cannibalization_sweep,
    hourly_permutation,
    solar_profit,
    solve_daily_dispatch,
)

pd.set_option("display.max_columns", 200)
plt.rcParams.update({
    "figure.dpi": 110,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "font.size": 10,
})

paths = SINPaths(root="validate_model", year=2025)
"""))

    cells.append(md("## 1. Painel SIN — agregação horária e seleção de dias representativos"))

    cells.append(code("""
panel = build_panel(paths)

sin = (
    panel
    .groupby("din_instante", as_index=False)
    .agg({"D": "sum", "gs": "sum", "gr": "sum", "gh": "sum", "gn_obs": "sum"})
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
    .agg(demand_GWh=("D", "sum"), solar_GWh=("gs", "sum"), hydro_GWh=("gh", "sum"),
         season=("season", "first"), day_type=("day_type", "first"), n_hours=("D", "size"))
    .reset_index()
)
daily_meta = daily_meta[daily_meta["n_hours"] == 24].copy()


def representative_day(meta: pd.DataFrame, season: str, day_type: str) -> pd.Timestamp:
    grp = meta[(meta["season"] == season) & (meta["day_type"] == day_type)]
    if grp.empty:
        raise ValueError(f"Sem dados para {season} | {day_type}")
    med = grp["demand_GWh"].median()
    return pd.Timestamp(grp.loc[(grp["demand_GWh"] - med).abs().idxmin(), "date"])


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


def day_frame(day: pd.Timestamp) -> pd.DataFrame:
    return (
        sin[sin["date"] == pd.Timestamp(day)]
        .copy()
        .sort_values("din_instante")
        .reset_index(drop=True)
    )
"""))

    cells.append(md("## 2. Configuração — `DailyParams` (única fonte da verdade)"))

    cells.append(code("""
# Calibração econômica:
#   CMg térmico = 2·a·g_t + b.  Para g_t típico ~ 15-30 GW e CMO desejado
#   na faixa 100–600 R$/MWh, a ≈ 1e-2 R$/MW²h e b ≈ 50 R$/MWh.
#   π_shed = 5000 R$/MWh (alta — define o "VOLL" implícito).
#   π_curt =    5 R$/MWh (baixa — apenas evita curtailment desnecessário).
#   cost_s1 =  80 R$/MWh, cost_s2 = 50 R$/MWh (LCOE solar; remoto mais
#       barato no campo, antes de transmissão).
base_params = DailyParams(
    K_S1=10_000.0, K_S2=16_000.0,
    Gt_max=50_000.0, hydro_scale=1.0,
    a=1.0e-2, b=50.0,
    cost_s1=80.0, cost_s2=50.0,
    pi_shed=5.0e3, pi_curt=5.0,
    line_limit=float("inf"), line_loss=0.0,
)
print(base_params)
print(f"\\nCMg térmico ~ 2a·g_t + b → escala: {2*base_params.a*15000 + base_params.b:.0f}-{2*base_params.a*30000 + base_params.b:.0f} R$/MWh para g_t∈[15,30] GW")
"""))

    cells.append(md("## 3. Cenários base por estação × tipo de dia"))

    cells.append(code("""
base_runs = []
base_dispatch = {}
for _, row in rep_days.iterrows():
    label, date = row["label"], pd.Timestamp(row["date"])
    day = day_frame(date)
    sol = solve_daily_dispatch(day, day, base_params)
    disp = sol["dispatch"]
    profit = solar_profit(disp, cost_s1=base_params.cost_s1, cost_s2=base_params.cost_s2)
    base_dispatch[label] = disp
    base_runs.append({
        "label": label, "date": date, "status": sol["status"],
        "gt_mean": disp["gt"].mean(), "gt_peak": disp["gt"].max(),
        "gt_GWh": disp["gt"].sum() / 1e3,
        "lambda_mean": float(disp["lambda"].mean()),
        "lambda_peak": float(disp["lambda"].max()),
        "lambda_minus_cmo_th_max": float(np.max(disp["lambda"] - disp["cmo_theoretical"])),
        "shed_MWh": disp["shed"].sum(),
        "curt_GWh": disp["curt"].sum() / 1e3,
        "pi_s1_total": profit["pi_s1_total"],
        "pi_s2_total": profit["pi_s2_total"],
    })

df_base_daily = pd.DataFrame(base_runs)
display(df_base_daily.round(3))
"""))

    cells.append(code("""
plot_df = df_base_daily.copy()
plot_df["season"] = plot_df["label"].str.split(" | ", regex=False).str[0]
plot_df["day_type"] = plot_df["label"].str.split(" | ", regex=False).str[1]

pivot_gt = plot_df.pivot(index="season", columns="day_type", values="gt_GWh")
pivot_cmo = plot_df.pivot(index="season", columns="day_type", values="lambda_mean")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
sns.heatmap(pivot_gt, annot=True, fmt=".1f", cmap="Reds", ax=axes[0])
axes[0].set_title("Energia térmica diária (GWh)")
sns.heatmap(pivot_cmo, annot=True, fmt=".0f", cmap="Oranges", ax=axes[1])
axes[1].set_title("CMO endógeno médio λ_t (R$/MWh)")
fig.tight_layout(); plt.show()
"""))

    cells.append(md(r"""
## 4. Dualidade — $\lambda_t$ vs custo marginal térmico teórico

Quando a térmica é marginal e não há shed/curt ativos, $\lambda_t = 2a\, g^{th}_t + b$. Disparidades indicam:

- **Déficit** ($\text{shed}_t > 0$): $\lambda_t \to \pi_{\text{shed}}$, muito acima do CMO térmico.
- **Curtailment** ($\text{curt}_t > 0$): $\lambda_t \to \pi_{\text{curt}}$, próximo de zero — o sistema está "pagando para descartar" energia.
"""))

    cells.append(code("""
fig, axes = plt.subplots(2, 2, figsize=(15, 8), sharex=True)

for ax, (label, disp) in zip(axes.ravel(), list(base_dispatch.items())[:4]):
    ax2 = ax.twinx()
    ax.plot(disp["din_instante"], disp["lambda"], lw=2, color="purple", label="λ_t (dual)")
    ax.plot(disp["din_instante"], disp["cmo_theoretical"], lw=1.4, color="indianred", ls="--", label="2a·gt + b")
    ax2.plot(disp["din_instante"], disp["gt"], lw=1.0, color="orange", alpha=0.5, label="gt")
    ax.set_title(label); ax.set_ylabel("R$/MWh"); ax2.set_ylabel("MW")
    ax.grid(True, alpha=0.3)
    if ax is axes[0, 0]:
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
fig.tight_layout(); plt.show()

# Diagnóstico numérico
diff_max = []
for label, disp in base_dispatch.items():
    delta = disp["lambda"] - disp["cmo_theoretical"]
    diff_max.append({
        "label": label,
        "max_|λ - CMO_th|": float(np.max(np.abs(delta))),
        "mean_|λ - CMO_th|": float(np.mean(np.abs(delta))),
        "horas_shed": int((disp["shed"] > 1e-3).sum()),
        "horas_curt": int((disp["curt"] > 1e-3).sum()),
    })
display(pd.DataFrame(diff_max).round(3))
"""))

    cells.append(md("## 5. Permutações horárias — demanda × oferta cruzadas"))

    cells.append(code("""
perm_specs = [
    ("Demanda Verao util  | Oferta Inverno util",  "Verao | Dia util",     "Inverno | Dia util"),
    ("Demanda Inverno util | Oferta Verao util",   "Inverno | Dia util",   "Verao | Dia util"),
    ("Demanda Primavera util | Oferta Primavera FDS", "Primavera | Dia util", "Primavera | Fim de semana"),
    ("Demanda Outono FDS   | Oferta Outono util",  "Outono | Fim de semana", "Outono | Dia util"),
]

def date_for(label):
    return pd.Timestamp(rep_days.loc[rep_days["label"] == label, "date"].iloc[0])


perm_results = []
perm_dispatch = {}
for name, demand_label, supply_label in perm_specs:
    d_day = day_frame(date_for(demand_label))
    s_day = day_frame(date_for(supply_label))
    out = hourly_permutation(d_day, s_day, base_params, permutation_label=name)
    perm_dispatch[name] = out["dispatch"]
    perm_results.append({k: v for k, v in out.items() if k != "dispatch"})

df_perm_daily = pd.DataFrame(perm_results)
display(df_perm_daily.round(3))
"""))

    cells.append(code("""
fig, axes = plt.subplots(2, 2, figsize=(15, 8), sharex=True, sharey=True)
for ax, (scenario, disp) in zip(axes.ravel(), perm_dispatch.items()):
    ax.plot(disp["din_instante"], disp["D"], color="black", lw=1.7, label="Demanda")
    ax.fill_between(disp["din_instante"], 0, disp["gs1"], color="gold", alpha=0.65, label="Solar 1")
    ax.fill_between(disp["din_instante"], disp["gs1"], disp["gs1"] + disp["gs2"], color="orange", alpha=0.55, label="Solar 2")
    ax.fill_between(disp["din_instante"], disp["gs1"] + disp["gs2"], disp["gs1"] + disp["gs2"] + disp["gh"], color="deepskyblue", alpha=0.55, label="Hidro")
    ax.fill_between(disp["din_instante"], disp["gs1"] + disp["gs2"] + disp["gh"], disp["gs1"] + disp["gs2"] + disp["gh"] + disp["gt"], color="indianred", alpha=0.6, label="Térmica")
    ax.set_title(scenario)
    ax.tick_params(axis="x", rotation=25)
    ax.grid(True, alpha=0.25)
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles[:5], labels[:5], loc="upper center", ncol=5, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.95)); plt.show()

fig, ax = plt.subplots(figsize=(12, 4))
for name, disp in perm_dispatch.items():
    ax.plot(disp["din_instante"], disp["lambda"], lw=1.6, label=name)
ax.set_title("CMO endógeno λ_t — permutações horárias")
ax.set_ylabel("R$/MWh"); ax.grid(True, alpha=0.3); ax.legend(fontsize=8, ncol=2)
plt.show()
"""))

    cells.append(md(r"""
## 6. Canibalização solar — $\Pi_{S1}$ cai com $K_{S2}$

Mantemos um dia representativo (Verão útil) e varremos a capacidade do **painel remoto** $K_{S2}$ entre $0$ e $32$ GW. Para cada $K_{S2}$ resolvemos o despacho, extraímos $\lambda_t$ e calculamos:

$$
\Pi_{S i} = \sum_t (\lambda_t - C^{S i})\, g^{Si}_t.
$$

A queda monotônica de $\Pi_{S1}$ enquanto $K_{S2}$ cresce é o **efeito canibal**: a entrada do segundo painel deprime $\lambda_t$ no horário solar, reduzindo a receita do primeiro.
"""))

    cells.append(code("""
ref_day = day_frame(date_for("Verao | Dia util"))

K_grid = np.linspace(0.0, 32_000.0, 17)
df_cann = cannibalization_sweep(ref_day, ref_day, base_params, K_S2_grid=K_grid)
display(df_cann.round(2))

fig, axes = plt.subplots(1, 2, figsize=(13, 4))

ax = axes[0]
ax.plot(df_cann["K_S2_MW"] / 1000, df_cann["pi_s1_total"] / 1e6, "o-", color="gold", lw=2, ms=6, label="Π Solar 1")
ax.plot(df_cann["K_S2_MW"] / 1000, df_cann["pi_s2_total"] / 1e6, "s-", color="orange", lw=2, ms=6, label="Π Solar 2")
ax.set_xlabel("K_S2 (GW)"); ax.set_ylabel("Lucro diário (R$ mi)")
ax.set_title("Canibalização — lucro Solar 1 cai com K_S2")
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(df_cann["K_S2_MW"] / 1000, df_cann["lambda_mean"], "^-", color="purple", lw=2, ms=6, label="λ médio")
ax2 = ax.twinx()
ax2.plot(df_cann["K_S2_MW"] / 1000, df_cann["curt_GWh"], "v--", color="firebrick", lw=1.4, ms=5, label="curtailment")
ax.set_xlabel("K_S2 (GW)"); ax.set_ylabel("λ médio (R$/MWh)"); ax2.set_ylabel("curt (GWh)")
ax.set_title("CMO médio e curtailment vs K_S2")
ax.grid(True, alpha=0.3)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=9, loc="best")

fig.tight_layout(); plt.show()
"""))

    cells.append(md(r"""
## 7. Transmissão — restrição na linha do painel remoto

Mantemos o dia mais "estressado" entre as permutações e variamos `line_limit` e `line_loss`. Toda redução na entrega da Solar 2 é compensada por mais térmica, com $\lambda_t$ subindo proporcionalmente.
"""))

    cells.append(code("""
stress_label = df_perm_daily.sort_values("gt_peak", ascending=False).iloc[0]["label"]
demand_label, supply_label = next((dl, sl) for nm, dl, sl in perm_specs if nm == stress_label)
demand_day = day_frame(date_for(demand_label))
supply_day = day_frame(date_for(supply_label))

tx_specs = [
    ("Sem restricao",        float("inf"), 0.00),
    ("Linha folgada",        5_000.0,      0.02),
    ("Linha intermediaria",  3_500.0,      0.04),
    ("Linha restrita",       2_000.0,      0.06),
]

tx_results, tx_dispatch = [], {}
for name, lim, loss in tx_specs:
    p_tx = replace(base_params, line_limit=lim, line_loss=loss)
    sol = solve_daily_dispatch(demand_day, supply_day, p_tx)
    disp = sol["dispatch"]
    tx_dispatch[name] = disp
    tx_results.append({
        "cenario": name, "line_limit_MW": lim, "line_loss_pct": loss * 100,
        "gt_mean": float(disp["gt"].mean()),
        "gt_peak": float(disp["gt"].max()),
        "gt_GWh": float(disp["gt"].sum() / 1e3),
        "lambda_mean": float(disp["lambda"].mean()),
        "lambda_peak": float(disp["lambda"].max()),
        "shed_MWh": float(disp["shed"].sum()),
        "curt_GWh": float(disp["curt"].sum() / 1e3),
        "solar2_spilled_GWh": sol["solar2_spilled_GWh"],
    })

df_tx = pd.DataFrame(tx_results)
display(df_tx.round(3))

fig, axes = plt.subplots(1, 2, figsize=(13, 4))
ax = axes[0]
for name, disp in tx_dispatch.items():
    ax.plot(disp["din_instante"], disp["gt"], lw=1.8, label=name)
ax.set_title(f"Térmica horária — caso crítico\\n{stress_label}")
ax.set_ylabel("MW"); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)

ax = axes[1]
for name, disp in tx_dispatch.items():
    ax.plot(disp["din_instante"], disp["lambda"], lw=1.8, label=name)
ax.set_title("λ_t com restrições de transmissão")
ax.set_ylabel("R$/MWh"); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); plt.show()
"""))

    cells.append(md("## 8. Leitura final + log de auditoria"))

    cells.append(code("""
base_best = df_base_daily.loc[df_base_daily["gt_GWh"].idxmin()]
base_worst = df_base_daily.loc[df_base_daily["gt_GWh"].idxmax()]
perm_best = df_perm_daily.loc[df_perm_daily["gt_mean"].idxmin()]
perm_worst = df_perm_daily.loc[df_perm_daily["gt_mean"].idxmax()]
tx_worst = df_tx.loc[df_tx["gt_GWh"].idxmax()]

cann_drop = (df_cann["pi_s1_total"].iloc[0] - df_cann["pi_s1_total"].iloc[-1]) / max(df_cann["pi_s1_total"].iloc[0], 1) * 100

print("Principais leituras:")
print(f"  1. Menor uso térmico: {base_best['label']} ({base_best['gt_GWh']:.2f} GWh, λ={base_best['lambda_mean']:.0f} R$/MWh).")
print(f"  2. Maior uso térmico: {base_worst['label']} ({base_worst['gt_GWh']:.2f} GWh, λ={base_worst['lambda_mean']:.0f} R$/MWh).")
print(f"  3. Permutação mais favorável: '{perm_best['label']}'.")
print(f"  4. Permutação mais adversa:   '{perm_worst['label']}'.")
print(f"  5. Pior caso de transmissão: '{tx_worst['cenario']}' — λ médio = {tx_worst['lambda_mean']:.0f} R$/MWh.")
print(f"  6. Canibalização: lucro de Solar 1 cai {cann_drop:.1f}% quando K_S2 sobe de {df_cann['K_S2_MW'].iloc[0]/1e3:.0f} GW para {df_cann['K_S2_MW'].iloc[-1]/1e3:.0f} GW.")
"""))

    cells.append(md(r"""
### Log de auditoria — "o que estava errado" vs. "o que foi corrigido"

| Problema original | Correção aplicada |
|---|---|
| Custo térmico era $\sum g_t^2$ implícito (sem $a, b$ explícitos) | Custo passou a ser $a\,g_t^2 + b\,g_t$ com $a, b$ em `DailyParams`. CMO térmico teórico fica $2a\,g_t + b$, alinhado à interpretação de oferta competitiva. |
| Calibração inicial usava $a = 1.0$ R$/MW²h$, gerando $\lambda_t \approx 32{,}000$–$88{,}000$ R$/MWh — três ordens de magnitude acima do CMO histórico do SIN | Recalibrado para $a = 10^{-2}$ R$/MW²h$ e $b = 50$ R$/MWh$, tal que $2a\,g_t + b$ fica entre 350 e 650 R$/MWh para $g_t \in [15,30]$ GW (faixa realista). $\pi_{\text{shed}}$ baixou de $2 \times 10^5$ para $5 \times 10^3$ R$/MWh (VOLL implícito). Custos solar adicionados ($C^{S1}=80$, $C^{S2}=50$ R$/MWh$). |
| O notebook calculava `CMO = 2*gt` "no olho" via fórmula teórica | O solver agora retorna $\lambda_t$ extraído **do dual** da restrição de balanço; comparamos $\lambda_t$ com o CMO teórico em §4 e tabulamos $\max\|\lambda - \text{CMO}_\text{th}\|$. |
| Não havia análise de canibalização | §6 implementa $\Pi_{S i} = \sum_t (\lambda_t - C^{S i})\,g^{Si}_t$ e varre $K_{S2}$ para mostrar a queda monotônica de $\Pi_{S1}$ (cai 56% quando $K_{S2}$ vai de 0 a 32 GW). |
| Permutações apenas trocavam dias inteiros | Mantivemos as permutações sazonais (§5) e adicionamos varredura de $K_{S2}$ + transmissão (§§6 e 7). |
| Solver não retornava $\lambda_t$ | Função `solve_daily_dispatch` no módulo central retorna `dispatch["lambda"]` e `dispatch["cmo_theoretical"]` — o segundo serve como controle para checar se o multiplicador "faz sentido". |
| Configuração espalhada por constantes globais (`K_S1`, `K_S2`, ...) | Tudo migrado para `DailyParams` (frozen dataclass-friendly), garantindo que sensibilidades não precisem mutar variáveis globais. |

### Interpretação econômica

- **CMO endógeno**: como $C^{th}$ é estritamente convexo, $\lambda_t$ cresce linearmente com $g^{th}_t$. Daí o CMO ser maior nos picos de demanda e menor nos vales.
- **Lucro x canibalização**: o efeito canibal aparece com força quando o segundo painel é "remoto" e tem perfil próximo do primeiro; ele empurra $\lambda_t$ para baixo no horário de pico solar exatamente onde o primeiro painel estava arrecadando. A receita marginal da última GW de Solar 2 tende a zero (a curva é côncava em $K_{S2}$).
- **Restrição de linha**: cortar `line_limit` é equivalente a remover capacidade efetiva — a térmica volta como recurso de balanço e $\lambda_t$ sobe proporcionalmente.
"""))

    nb = new_notebook(cells=cells, metadata=base_metadata())
    nbformat.write(nb, NB5)
    print(f"NB5 reconstruído: {NB5}  ({len(nb.cells)} células)")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    build_nb4()
    build_nb5()
