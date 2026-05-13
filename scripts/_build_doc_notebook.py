"""Construi o notebook 06_implementacao_documento_mfg.ipynb.

Estrategia: usar r'''...''' para wrappers e manter qualquer triple-double-quote
dentro das celulas (docstrings de funcoes nao causam conflito).

Estrutura segue Apendice C do documento final.
"""
from __future__ import annotations
import nbformat as nbf
from pathlib import Path

NB_PATH = Path(r"c:\Users\joaof\Mean-Field-Games-in-Renewable-Energy\06_implementacao_documento_mfg.ipynb")

nb = nbf.v4.new_notebook()
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.13"},
}
cells: list = []


def md(src: str) -> None:
    cells.append(nbf.v4.new_markdown_cell(src.strip("\n")))


def code(src: str) -> None:
    cells.append(nbf.v4.new_code_cell(src.strip("\n")))


# =========================================================================
md(r'''
# Intermitência renovável, despacho hidrotérmico e equilíbrio locacional

## Implementação completa do documento MFG — Final

**Autor da implementação:** Claude (com base no manuscrito de J. F. Vilas Boas)
**Período de dados:** 2025‑01‑01 → **2025‑10‑31** (Sistema Interligado Nacional)
**Subsistemas:** N · NE · SE · S

Este notebook reproduz, em uma única execução end‑to‑end, **todos os blocos
empíricos e modelos** sugeridos no documento (Apêndice C):

1. Importação e padronização dos dados horários do ONS
2. Construção de calendário (estação × tipo de dia × hora)
3. Limpeza de carga, geração por fonte, preço (CMO/PLD) e hidrologia
4. Cálculo de fatores de capacidade e perfis típicos
5. Estimativa de demanda líquida `D^net`, residual `D^res` e rampas
6. Curva do pato e contrafactual de expansão solar (`α = 1, 1.5, 2, 3`)
7. Fator de captura (H1 — canibalização)
8. Modelo hidrelétrico reduzido (estado de armazenamento, valor da água)
9. Despacho centralizado (planejador social) — benchmark `Q(K)`
10. Modelo situacional finito (2 solares + térmica + hidrelétrica)
11. Análise de sensibilidade de `K_2`, `F_i`, `q_i`, hidrologia (H2–H5)
12. Iteração MFG locacional — preço‑sombra λ_{ℓ,t} por subsistema
13. Comparação ótimo social × jogo finito × MFG
14. Apêndice B — tabela consolidada de métricas

> **Saídas:** todas as figuras (`figs/`) e tabelas (`tables/`) são exportadas
> para `outputs_doc/`. As figuras seguem a numeração do manuscrito
> (Figura 1 → Figura 19).
''')


md(r'''
## 0. Setup, imports e diretórios de saída
''')

code(r'''
from __future__ import annotations
import os, sys, math, json, warnings, time
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import cvxpy as cp

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 140,
    "axes.titleweight": "bold",
    "axes.spines.right": False,
    "axes.spines.top": False,
    "font.size": 10,
})

ROOT      = Path(r"c:/Users/joaof/Mean-Field-Games-in-Renewable-Energy")
DATA      = ROOT / "validate_model" / "data"
OUT       = ROOT / "outputs_doc"
FIGS      = OUT / "figs"
TABLES    = OUT / "tables"
CACHE     = OUT / "cache"
for d in (OUT, FIGS, TABLES, CACHE):
    d.mkdir(parents=True, exist_ok=True)

DATE_MIN = pd.Timestamp("2025-01-01 00:00:00")
DATE_MAX = pd.Timestamp("2025-10-31 23:00:00")

SUBSYS_ORDER = ["N", "NE", "SE", "S"]
SUB_LABEL = {"N": "Norte", "NE": "Nordeste", "SE": "Sudeste/CO", "S": "Sul"}
SUB_COLOR = {"N": "#1f77b4", "NE": "#d62728", "SE": "#2ca02c", "S": "#9467bd"}

sys.path.insert(0, str(ROOT))
from validate_model.pipeline import (
    SINPaths, build_panel, load_cmo_semihorario, load_cmo_horario,
    load_hidrologia, canonical_subsys, load_demanda_efetiva,
    load_generation, load_intercambio_interno,
)

print("Pacotes OK · cvxpy", cp.__version__, "· pandas", pd.__version__)
print("Diretórios prontos em:", OUT)
''')


md(r'''
## 1. Importação e padronização (Apêndice C, passo i)

Usamos `validate_model.pipeline` para padronizar separadores, encodings e
nomes de colunas. Toda observação horária recebe (i) `id_subsistema` em
código canônico (`N|NE|SE|S`); (ii) `din_instante` truncado à hora; (iii)
valores numéricos no padrão brasileiro convertidos para `float`.

Para evitar releitura dos CSVs (>500 MB no total), o painel unificado é
persistido em **parquet** em `outputs_doc/cache/`.
''')

code(r'''
PANEL_CACHE = CACHE / "panel_hourly.parquet"
CMO_CACHE   = CACHE / "cmo_hourly.parquet"
HYDRO_CACHE = CACHE / "hidrologia.parquet"

paths = SINPaths(root=str(ROOT / "validate_model"), year=2025)

if PANEL_CACHE.exists():
    panel = pd.read_parquet(PANEL_CACHE)
    print(f"[CACHE] painel: {len(panel):,} linhas")
else:
    panel = build_panel(paths)
    panel.to_parquet(PANEL_CACHE, index=False)

if CMO_CACHE.exists():
    cmo = pd.read_parquet(CMO_CACHE)
else:
    cmo = load_cmo_horario(paths.cmo_semihorario_path)
    cmo.to_parquet(CMO_CACHE, index=False)

if HYDRO_CACHE.exists():
    hydro = pd.read_parquet(HYDRO_CACHE)
else:
    hydro = load_hidrologia(paths.hidro_di_path)
    hydro.to_parquet(HYDRO_CACHE, index=False)

panel = panel[(panel["din_instante"] >= DATE_MIN) & (panel["din_instante"] <= DATE_MAX)].copy()
cmo   = cmo[(cmo["din_instante"] >= DATE_MIN) & (cmo["din_instante"] <= DATE_MAX)].copy()
hydro = hydro[(hydro["din_instante"] >= DATE_MIN) & (hydro["din_instante"] <= DATE_MAX)].copy()

print(f"Painel ........: {len(panel):>10,} linhas | {panel['din_instante'].min()} -> {panel['din_instante'].max()}")
print(f"CMO horário ..: {len(cmo):>10,} linhas | {cmo['din_instante'].min()} -> {cmo['din_instante'].max()}")
print(f"Hidrologia ...: {len(hydro):>10,} linhas | {hydro['din_instante'].min()} -> {hydro['din_instante'].max()}")
panel.head()
''')


md(r'''
## 2. Calendário horário, estação e tipo de dia (Seção 5.1)

Para cada timestamp:

- `mes` ∈ {1..12}
- `estacao` ∈ {`verao`, `outono`, `inverno`, `primavera`} — convenção do
  hemisfério sul: DJF = verão, MAM = outono, JJA = inverno, SON = primavera.
- `dia_semana` ∈ {0=Seg..6=Dom}
- `tipo_dia` ∈ {`util`, `sabado`, `domingo`, `feriado`}
- `hora` ∈ {0..23}
''')

code(r'''
SOUTH_HEMI_SEASON = {
    12: "verao", 1: "verao", 2: "verao",
    3: "outono", 4: "outono", 5: "outono",
    6: "inverno", 7: "inverno", 8: "inverno",
    9: "primavera", 10: "primavera", 11: "primavera",
}
SEASON_ORDER = ["verao", "outono", "inverno", "primavera"]

FERIADOS_2025 = {
    "2025-01-01", "2025-03-03", "2025-03-04", "2025-04-18",
    "2025-04-21", "2025-05-01", "2025-06-19", "2025-09-07",
    "2025-10-12", "2025-11-02", "2025-11-15", "2025-11-20",
    "2025-12-25",
}

def add_calendar(df, dt_col="din_instante"):
    df = df.copy()
    t = df[dt_col]
    df["data"]       = t.dt.date
    df["ano"]        = t.dt.year
    df["mes"]        = t.dt.month
    df["dia_semana"] = t.dt.dayofweek
    df["hora"]       = t.dt.hour
    df["estacao"]    = df["mes"].map(SOUTH_HEMI_SEASON)
    feriados = pd.to_datetime(list(FERIADOS_2025)).date
    is_feriado = df["data"].isin(feriados)
    df["tipo_dia"] = np.select(
        [is_feriado, df["dia_semana"] == 5, df["dia_semana"] == 6],
        ["feriado", "sabado", "domingo"],
        default="util",
    )
    return df

panel = add_calendar(panel)
cmo   = add_calendar(cmo)
panel.groupby("estacao", observed=False).size().reindex(SEASON_ORDER)
''')


md(r'''
## 3. Limpeza e definição das variáveis empíricas (Seções 3.1, 4.1, 5)

$$
D^{\text{net}}_{\ell,t} = D_{\ell,t} - G^{s}_{\ell,t}
\qquad
D^{\text{res}}_{\ell,t} = D_{\ell,t} - G^{s}_{\ell,t} - G^{w}_{\ell,t}
\qquad
\Delta D^{\text{res}}_{\ell,t} = D^{\text{res}}_{\ell,t+1} - D^{\text{res}}_{\ell,t}
$$

`D_t` é a carga observada da Rede Básica (já líquida de MMGD).
Para `D^{global}` adicionamos a geração distribuída (PQU MMGD).
''')

code(r'''
print("Re-lendo fotovoltaicas para isolar MMGD vs centralizada...")
sol_raw = pd.read_csv(paths.solar_path, sep=";", parse_dates=["din_instante"])
sol_raw["din_instante"] = sol_raw["din_instante"].dt.floor("h")
sol_raw["id_subsistema"] = sol_raw["nom_subsistema"].astype(str).map(canonical_subsys)
sol_raw = sol_raw[(sol_raw["din_instante"] >= DATE_MIN) & (sol_raw["din_instante"] <= DATE_MAX)]
is_mmgd = sol_raw["cod_modalidadeoperacao"].astype(str).str.contains("MMGD", case=False, regex=False)

sol_mmgd = (sol_raw[is_mmgd]
    .groupby(["din_instante", "id_subsistema"], as_index=False)["val_geracao"].sum()
    .rename(columns={"val_geracao": "g_mmgd"}))
sol_cen  = (sol_raw[~is_mmgd]
    .groupby(["din_instante", "id_subsistema"], as_index=False)["val_geracao"].sum()
    .rename(columns={"val_geracao": "g_solar_cen"}))

panel = (panel
    .merge(sol_mmgd, on=["din_instante", "id_subsistema"], how="left")
    .merge(sol_cen,  on=["din_instante", "id_subsistema"], how="left"))
panel["g_mmgd"]      = panel["g_mmgd"].fillna(0.0)
panel["g_solar_cen"] = panel["g_solar_cen"].fillna(0.0)

panel = panel.rename(columns={"D": "D_rb"})
panel["D_global"] = panel["D_rb"] + panel["g_mmgd"]
panel["D_net"]    = panel["D_global"] - panel["gs"]
panel["D_res"]    = panel["D_global"] - panel["gs"] - panel["gr"]

panel = panel.sort_values(["id_subsistema", "din_instante"])
panel["dD_res"]   = panel.groupby("id_subsistema")["D_res"].diff()
panel["dD_net"]   = panel.groupby("id_subsistema")["D_net"].diff()
panel["dD_glob"]  = panel.groupby("id_subsistema")["D_global"].diff()

print(f"\nLinhas finais .........: {len(panel):,}")
print(f"Subsistemas ...........: {sorted(panel['id_subsistema'].dropna().unique())}")
print(f"Horizonte temporal ....: {panel['din_instante'].min()} -> {panel['din_instante'].max()}")

panel.groupby("id_subsistema").agg(
    D_rb_mean=("D_rb", "mean"),
    D_global_mean=("D_global", "mean"),
    Gs_mean=("gs", "mean"),
    Gw_mean=("gr", "mean"),
    Gh_mean=("gh", "mean"),
    Gth_mean=("g_th_obs", "mean"),
).round(1)
''')


md(r'''
### 3.1 Capacidades instaladas estimadas (proxy `K_ℓ^r`)

Estimamos `K_ℓ^r` pelo percentil 99 da geração horária observada.
''')

code(r'''
K_inst = (panel.groupby("id_subsistema")
          .agg(K_s=("gs", lambda x: float(np.nanquantile(x, 0.99))),
               K_w=("gr", lambda x: float(np.nanquantile(x, 0.99))),
               K_h=("gh", lambda x: float(np.nanquantile(x, 0.99))),
               K_t=("g_th_obs", lambda x: float(np.nanquantile(x, 0.99))))
          .reindex(SUBSYS_ORDER))
K_inst.to_csv(TABLES / "capacidades_estimadas_p99.csv", float_format="%.1f")
print("Capacidades estimadas (MW, p99 horário):")
display(K_inst.round(0))

for fonte, kcol in [("gs", "K_s"), ("gr", "K_w")]:
    a = panel[fonte] / panel["id_subsistema"].map(K_inst[kcol])
    panel[f"a_{fonte}"] = a.clip(lower=0.0, upper=1.5)

print("\nFator de capacidade médio por subsistema:")
panel.groupby("id_subsistema").agg(a_s_mean=("a_gs", "mean"),
                                   a_w_mean=("a_gr", "mean")).round(3)
''')


md(r'''
## 4. Perfis horários por estação e tipo de dia (Seções 6.1, 3.1)

### Figura 1 — `D^global`, `D^RB` e `D^res` por estação × tipo de dia
''')

code(r'''
def avg_profile(df, cols, by=("id_subsistema", "estacao", "tipo_dia_g", "hora")):
    d = df.copy()
    d["tipo_dia_g"] = np.where(d["tipo_dia"].isin(["util"]), "util", "fim_semana_feriado")
    return d.groupby(list(by), observed=False)[list(cols)].mean().reset_index()

prof_load = avg_profile(panel, ["D_global", "D_rb", "D_res", "gs", "gr", "gh", "g_th_obs"])
prof_load.to_csv(TABLES / "perfis_horarios_carga.csv", index=False, float_format="%.1f")

fig, axes = plt.subplots(4, 2, figsize=(14, 14), sharex=True)
for j, td in enumerate(["util", "fim_semana_feriado"]):
    for i, sub in enumerate(SUBSYS_ORDER):
        ax = axes[i, j]
        for est, ls in zip(SEASON_ORDER, ["-", "--", "-.", ":"]):
            g = prof_load.query("id_subsistema == @sub and estacao == @est and tipo_dia_g == @td")
            if len(g) == 0:
                continue
            ax.plot(g["hora"], g["D_global"]/1000, color=SUB_COLOR[sub], ls=ls, lw=2, label=f"D_global · {est}")
            ax.plot(g["hora"], g["D_rb"]/1000,     color="gray",           ls=ls, lw=1)
            ax.plot(g["hora"], g["D_res"]/1000,    color="black",          ls=ls, lw=1)
        ax.set_title(f"[{sub}] {SUB_LABEL[sub]} - {td}")
        ax.set_ylabel("GW")
        ax.set_xticks(range(0, 24, 3))
        ax.grid(alpha=.3)

axes[0, 0].legend(fontsize=8, ncols=2, loc="upper right")
axes[-1, 0].set_xlabel("hora")
axes[-1, 1].set_xlabel("hora")
fig.suptitle("Figura 1 - Carga global (D_global), Rede Basica (D_RB) e residual (D_res) por subsistema x estacao x tipo_dia (2025)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig01_carga_global_rb_res.png", bbox_inches="tight")
plt.show()
''')


md(r'''
### Figura 3 — Perfis horários de carga, solar e eólica
''')

code(r'''
def quantile_band(df, val, sub):
    g = df[df["id_subsistema"] == sub]
    q = g.groupby("hora")[val].quantile([0.1, 0.5, 0.9]).unstack().reindex(range(24))
    return q.rename(columns={0.1: "q10", 0.5: "q50", 0.9: "q90"})

fig, axes = plt.subplots(4, 3, figsize=(16, 14), sharex=True)
metrics = [("D_global", "Carga global (MW)"),
           ("gs",       "Geracao solar (MW)"),
           ("gr",       "Geracao eolica (MW)")]

for i, sub in enumerate(SUBSYS_ORDER):
    for j, (col, title) in enumerate(metrics):
        ax = axes[i, j]
        q = quantile_band(panel, col, sub)
        ax.fill_between(q.index, q["q10"], q["q90"], color=SUB_COLOR[sub], alpha=.25, label="P10-P90")
        ax.plot(q.index, q["q50"], color=SUB_COLOR[sub], lw=2, label="Mediana")
        ax.set_title(f"[{sub}] {title}")
        ax.grid(alpha=.3)
        if i == 0 and j == 0:
            ax.legend()
        if j == 0:
            ax.set_ylabel(SUB_LABEL[sub])
axes[-1, 0].set_xlabel("hora")
axes[-1, 1].set_xlabel("hora")
axes[-1, 2].set_xlabel("hora")
fig.suptitle("Figura 3 - Perfis horarios de carga, solar e eolica (mediana + banda P10-P90, 2025)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig03_perfis_horarios_carga_solar_eolica.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 5. Diagnóstico operacional do SIN (Figura 2)

Heatmaps por (subsistema × estação × tipo de dia):
horário de carga residual mínima, rampa máxima e fração de horas com excedente.
''')

code(r'''
d = panel.copy()
d["tipo_dia_g"] = np.where(d["tipo_dia"] == "util", "util", "fim_semana_feriado")

ag = d.groupby(["id_subsistema", "estacao", "tipo_dia_g", "data"])
hora_min = ag.apply(lambda x: x.loc[x["D_res"].idxmin(), "hora"] if x["D_res"].notna().any() else np.nan,
                    include_groups=False).reset_index(name="hora_min_dres")
hora_min_avg = hora_min.groupby(["id_subsistema", "estacao", "tipo_dia_g"])["hora_min_dres"].mean().reset_index()

ramp = ag["dD_res"].apply(lambda x: x.abs().max()).reset_index(name="rampa_max_dres")
ramp_avg = ramp.groupby(["id_subsistema", "estacao", "tipo_dia_g"])["rampa_max_dres"].mean().reset_index()

d["pot_exc"] = ((d["gs"] + d["gr"]) > 0.7 * d["D_global"]).astype(int)
exc = d.groupby(["id_subsistema", "estacao", "tipo_dia_g"])["pot_exc"].mean().reset_index()

fig, axes = plt.subplots(1, 3, figsize=(16, 4))
def heatmap(df, val, title, ax, fmt=".1f", cmap="viridis"):
    piv = df.pivot_table(index="id_subsistema", columns=["tipo_dia_g", "estacao"], values=val)
    piv = piv.reindex(index=SUBSYS_ORDER)
    sns.heatmap(piv, annot=True, fmt=fmt, cmap=cmap, ax=ax, cbar_kws={"shrink": .7})
    ax.set_title(title)
    ax.set_xlabel("")

heatmap(hora_min_avg, "hora_min_dres", "Hora de carga residual minima", axes[0], fmt=".1f", cmap="YlOrRd")
heatmap(ramp_avg,     "rampa_max_dres", "Rampa maxima |dD_res| (MW/h)", axes[1], fmt=".0f", cmap="Reds")
heatmap(exc,          "pot_exc",        "Fracao de horas com excedente (>70% D)", axes[2], fmt=".2f", cmap="Blues")
fig.suptitle("Figura 2 - Diagnostico operacional do SIN por estacao x tipo_dia (2025)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig02_diagnostico_operacional.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 6. Curva do pato e contrafactual de expansão solar (Figuras 4 e 14)

Para `α ∈ {1, 1.5, 2, 3}`:
$$
D^{\text{net}(\alpha)}_t = D_t - \alpha\, G^{s}_t
$$
''')

code(r'''
ALPHAS = [1.0, 1.5, 2.0, 3.0]

d = panel.copy()
d["tipo_dia_g"] = np.where(d["tipo_dia"] == "util", "util", "fim_semana_feriado")
typ = d.groupby(["id_subsistema", "estacao", "tipo_dia_g", "hora"], observed=False).agg(
    D_global=("D_global", "mean"),
    gs=("gs", "mean"),
    gr=("gr", "mean"),
    gh=("gh", "mean"),
    g_th_obs=("g_th_obs", "mean"),
).reset_index()

ramp_summary = []
for alpha in ALPHAS:
    typ[f"D_net_a{alpha}"] = typ["D_global"] - alpha * typ["gs"]
    typ[f"dD_net_a{alpha}"] = (typ.sort_values("hora")
                                  .groupby(["id_subsistema", "estacao", "tipo_dia_g"])[f"D_net_a{alpha}"]
                                  .diff())
    r = typ.groupby(["id_subsistema", "estacao", "tipo_dia_g"])[f"dD_net_a{alpha}"].agg(
        rampa_max_pos="max", rampa_max_neg="min").reset_index()
    r["alpha"] = alpha
    ramp_summary.append(r)
ramp_summary = pd.concat(ramp_summary, ignore_index=True)
ramp_summary.to_csv(TABLES / "rampa_max_por_alpha.csv", index=False, float_format="%.0f")

fig, axes = plt.subplots(4, 2, figsize=(14, 14), sharex=True)
for j, td in enumerate(["util", "fim_semana_feriado"]):
    for i, sub in enumerate(SUBSYS_ORDER):
        ax = axes[i, j]
        g = typ.query("id_subsistema == @sub and tipo_dia_g == @td")
        gp = g.groupby("hora").mean(numeric_only=True).reset_index()
        ax.plot(gp["hora"], gp["D_global"]/1000, "k-", lw=2, label="D_global (obs.)")
        for alpha, col in zip(ALPHAS, ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]):
            ax.plot(gp["hora"], gp[f"D_net_a{alpha}"]/1000, color=col, lw=1.5,
                    label=f"D_net alpha={alpha}")
        ax.set_title(f"[{sub}] {SUB_LABEL[sub]} - {td}")
        ax.set_ylabel("GW")
        ax.grid(alpha=.3)
axes[0, 0].legend(fontsize=8, ncols=2)
axes[-1, 0].set_xlabel("hora"); axes[-1, 1].set_xlabel("hora")
fig.suptitle("Figura 4 / 14 - Curva do pato observada e contrafactual de expansao solar (alpha em {1, 1.5, 2, 3})",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig04_14_curva_do_pato_contrafactual.png", bbox_inches="tight")
plt.show()

display(ramp_summary.pivot_table(index=["id_subsistema", "tipo_dia_g"],
                                 columns="alpha", values="rampa_max_pos").round(0))
''')


md(r'''
## 7. Preços (CMO) e fator de captura solar (Seções 4.2, 5.3, 6.3)

### 7.1 Dispersão CMO × demanda líquida (Figura 5)
''')

code(r'''
panel = panel.merge(cmo[["din_instante", "id_subsistema", "cmo_h"]],
                    on=["din_instante", "id_subsistema"], how="left")

fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=False)
for ax, sub in zip(axes.flat, SUBSYS_ORDER):
    g = panel.dropna(subset=["cmo_h", "D_net"])[panel["id_subsistema"] == sub]
    sc = ax.scatter(g["D_net"]/1000, g["cmo_h"], c=g["hora"], cmap="twilight",
                    s=4, alpha=.4, rasterized=True)
    ax.set_title(f"[{sub}] {SUB_LABEL[sub]}: CMO horario x D_net")
    ax.set_xlabel("D_net (GW)")
    ax.set_ylabel("CMO (R$/MWh)")
    ax.grid(alpha=.3)
fig.colorbar(sc, ax=axes.ravel().tolist(), label="hora do dia", shrink=.7)
fig.suptitle("Figura 5 - CMO horario x demanda liquida por subsistema (2025)",
             fontsize=12, fontweight="bold")
fig.savefig(FIGS / "fig05_cmo_vs_dnet.png", bbox_inches="tight")
plt.show()
''')


md(r'''
### 7.2 Fator de captura solar/eólica (H1 — Figuras 13 & 15)

$$
F^{\text{capture}}_{\ell,r} = \frac{P^{\text{cap}}_{\ell,r}}{P^{\text{avg}}_{\ell}}
$$

Também simulamos um contrafactual usando `P̂(D_net)` polinômio de 2º grau.
''')

code(r'''
from numpy.polynomial.polynomial import polyfit, polyval

def fcap_obs(df, sub, gen_col):
    g = df[(df["id_subsistema"] == sub) & df["cmo_h"].notna() & df[gen_col].notna()]
    P_avg  = g["cmo_h"].mean()
    Es     = g[gen_col].sum()
    if Es <= 0:
        return np.nan, P_avg, np.nan
    P_cap  = (g["cmo_h"] * g[gen_col]).sum() / Es
    return P_cap / P_avg if P_avg > 0 else np.nan, P_avg, P_cap

records = []
for sub in SUBSYS_ORDER:
    for col, src in [("gs", "solar"), ("gr", "wind")]:
        Fcap, P_avg, P_cap = fcap_obs(panel, sub, col)
        records.append({"id_subsistema": sub, "fonte": src,
                        "F_capture": Fcap, "P_avg": P_avg, "P_cap": P_cap})

fcap_df = pd.DataFrame(records)
fcap_df.to_csv(TABLES / "fator_captura_observado.csv", index=False, float_format="%.4f")
display(fcap_df.round(3))

def fit_price_curve(df, sub):
    g = df[(df["id_subsistema"] == sub)].dropna(subset=["cmo_h", "D_net"])
    return polyfit(g["D_net"].values, g["cmo_h"].values, 2)

fc_rows = []
for sub in SUBSYS_ORDER:
    coefs = fit_price_curve(panel, sub)
    g = panel[panel["id_subsistema"] == sub].dropna(subset=["D_global", "gs"]).copy()
    for alpha in [1.0, 1.5, 2.0, 3.0]:
        Dn = (g["D_global"] - alpha * g["gs"]).values
        Phat = np.clip(polyval(Dn, coefs), 0, None)
        Gs   = (alpha * g["gs"]).values
        if Gs.sum() <= 0:
            continue
        Fc = (Phat * Gs).sum() / Gs.sum() / max(Phat.mean(), 1e-9)
        fc_rows.append({"id_subsistema": sub, "alpha": alpha, "F_capture_hat": Fc,
                        "P_avg_hat": float(Phat.mean()),
                        "P_cap_hat": float((Phat * Gs).sum() / Gs.sum())})
fc_alpha = pd.DataFrame(fc_rows)
fc_alpha.to_csv(TABLES / "fator_captura_contrafactual.csv", index=False, float_format="%.4f")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
pivot = fcap_df.pivot(index="id_subsistema", columns="fonte", values="F_capture").reindex(SUBSYS_ORDER)
pivot.plot.bar(ax=ax, color=["#fdae61", "#74add1"])
ax.axhline(1.0, color="black", ls="--", lw=1)
ax.set_ylabel("Fator de captura F")
ax.set_title("Figura 13 - F_capture observado por fonte/subsistema (2025)")
ax.set_ylim(0, max(1.2, pivot.max().max()*1.1))
ax.legend(title="fonte")

ax = axes[1]
for sub in SUBSYS_ORDER:
    g = fc_alpha[fc_alpha["id_subsistema"] == sub]
    ax.plot(g["alpha"], g["F_capture_hat"], "o-", color=SUB_COLOR[sub], label=f"{sub}")
ax.axhline(1.0, color="black", ls="--", lw=1)
ax.set_xlabel("alpha - escala de expansao solar")
ax.set_ylabel("F_capture solar (preco estimado)")
ax.set_title("Figura 15 - Canibalizacao contrafactual da receita solar")
ax.legend()
ax.grid(alpha=.3)
fig.tight_layout()
fig.savefig(FIGS / "fig13_15_fator_captura.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 8. Hidrologia e valor de oportunidade (Seções 6.4, 7) — Figura 6
''')

code(r'''
hydro_sin = hydro.copy()
hydro_sin["Storage_MW"] = hydro_sin["A_MW"].cumsum() - hydro_sin["Outflow_MW"].cumsum()

if hydro_sin["vol_util_pct"].notna().any():
    s = hydro_sin["vol_util_pct"]
else:
    s = hydro_sin["Storage_MW"]
hydro_sin["Z_hyd"] = (s - s.min()) / max(s.max() - s.min(), 1e-9)

fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)
ax = axes[0, 0]
ax.plot(hydro_sin["din_instante"], hydro_sin["A_MW"]/1000, color="#1f77b4", label="Afluencia")
ax.plot(hydro_sin["din_instante"], hydro_sin["T_MW"]/1000, color="#d62728", label="Turbinada")
ax.set_ylabel("GWmed")
ax.set_title("Afluencia vs turbinamento (SIN, diario)")
ax.legend(); ax.grid(alpha=.3)

ax = axes[0, 1]
ax.plot(hydro_sin["din_instante"], hydro_sin["Spill_MW"]/1000, color="#9467bd", label="Vertimento")
ax.plot(hydro_sin["din_instante"], hydro_sin["Outflow_MW"]/1000, color="black", label="Defluente")
ax.set_ylabel("GWmed")
ax.set_title("Vertimento e defluencia")
ax.legend(); ax.grid(alpha=.3)

ax = axes[1, 0]
ax.plot(hydro_sin["din_instante"], hydro_sin["vol_util_pct"], color="#2ca02c")
ax.set_ylabel("% vol. util")
ax.set_title("Estado de armazenamento medio (RCU)")
ax.grid(alpha=.3)

ax = axes[1, 1]
ax.plot(hydro_sin["din_instante"], hydro_sin["Z_hyd"], color="#ff7f0e")
ax.set_ylabel("Z_hyd em [0, 1]")
ax.set_title("Estado hidrologico normalizado Z_hyd")
ax.grid(alpha=.3)

for ax in axes.flat:
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))

fig.suptitle("Figura 6 - Series hidrologicas SIN (diarias) - 2025",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig06_series_hidrologicas.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 9. Fatos estilizados de 2025 (Figura 13)
''')

code(r'''
fig, axes = plt.subplots(2, 2, figsize=(14, 9))

ax = axes[0, 0]
g = panel.groupby(["id_subsistema", "hora"])["D_global"].mean().reset_index()
for sub in SUBSYS_ORDER:
    s = g[g["id_subsistema"] == sub]
    ax.plot(s["hora"], s["D_global"]/1000, color=SUB_COLOR[sub], lw=2, label=sub)
ax.set_title("(a) Carga global media horaria"); ax.set_ylabel("GW"); ax.legend(ncols=4)
ax.grid(alpha=.3)

ax = axes[0, 1]
g = panel.groupby(["id_subsistema", "hora"])["a_gs"].mean().reset_index()
for sub in SUBSYS_ORDER:
    s = g[g["id_subsistema"] == sub]
    ax.plot(s["hora"], s["a_gs"], color=SUB_COLOR[sub], lw=2, label=sub)
ax.set_title("(b) Fator de capacidade solar a_s,h"); ax.set_ylabel("fracao")
ax.grid(alpha=.3); ax.legend(ncols=4)

ax = axes[1, 0]
g = panel.groupby(["id_subsistema", "hora"])["a_gr"].mean().reset_index()
for sub in SUBSYS_ORDER:
    s = g[g["id_subsistema"] == sub]
    ax.plot(s["hora"], s["a_gr"], color=SUB_COLOR[sub], lw=2, label=sub)
ax.set_title("(c) Fator de capacidade eolica a_w,h"); ax.set_ylabel("fracao")
ax.set_xlabel("hora"); ax.grid(alpha=.3); ax.legend(ncols=4)

ax = axes[1, 1]
g = panel.groupby(["id_subsistema", "hora"])["cmo_h"].mean().reset_index()
for sub in SUBSYS_ORDER:
    s = g[g["id_subsistema"] == sub]
    ax.plot(s["hora"], s["cmo_h"], color=SUB_COLOR[sub], lw=2, label=sub)
ax.set_title("(d) CMO medio horario (R$/MWh)"); ax.set_ylabel("R$/MWh")
ax.set_xlabel("hora"); ax.grid(alpha=.3); ax.legend(ncols=4)

fig.suptitle("Figura 13 - Fatos estilizados dos dados de 2025 (4 paineis)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig13_fatos_estilizados.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 10. Calibração de parâmetros do despacho (Figura 7)
''')

code(r'''
g = panel.dropna(subset=["cmo_h"]).copy()
g["g_disp"] = g["g_th_obs"] + g["g_nuc_obs"]
y = g["cmo_h"].values
x = g["g_disp"].values

A = np.column_stack([np.ones_like(x), x])
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
c1_lin, c2_lin = float(coef[0]), float(coef[1])

c1 = max(c1_lin, 50.0)
c2 = max(c2_lin, 1e-4)

agg_panel = panel.groupby("din_instante", as_index=False).agg(
    D_global=("D_global", "sum"),
    gs=("gs", "sum"), gr=("gr", "sum"),
    gh=("gh", "sum"), g_th=("g_th_obs", "sum"),
    g_nuc=("g_nuc_obs", "sum"),
)
K_T = float(np.nanquantile(agg_panel["g_th"], 0.995))
K_H = float(np.nanquantile(agg_panel["gh"],  0.995))

dG_th = agg_panel["g_th"].diff().abs().dropna()
rho_up = rho_dn = float(np.nanquantile(dG_th, 0.99))

PARAMS = {
    "c1_RS_per_MWh":      round(c1, 2),
    "c2_RS_per_MWh2":     round(c2, 6),
    "K_T_MW":             round(K_T, 0),
    "K_H_MW":             round(K_H, 0),
    "rho_up_MW_per_h":    round(rho_up, 0),
    "rho_dn_MW_per_h":    round(rho_dn, 0),
    "pi_c_RS_per_MWh":    0.0,
    "pi_u_RS_per_MWh":    3500.0,
    "VOLL":               3500.0,
    "kappa_hydro":        0.9 * 9.81 / 1000,
    "alpha_S_inicial":    0.55,
    "beta_S_target":      0.55,
}
pd.DataFrame.from_dict(PARAMS, orient="index", columns=["valor"]).to_csv(
    TABLES / "params_dispatch_calibrados.csv")
print("Parametros calibrados (Figura 7):")
display(pd.DataFrame.from_dict(PARAMS, orient="index", columns=["valor"]).round(3))
''')


md(r'''
## 11. Despacho centralizado (Seção 8)

Para cada (subsistema × estação × tipo_dia), resolvemos o problema
uninodal com 24 horas. `λ_t` é o multiplicador do balanço.
''')

code(r'''
def typical_day(panel, sub, season, td):
    g = panel.query("id_subsistema == @sub").copy()
    g["tipo_dia_g"] = np.where(g["tipo_dia"] == "util", "util", "fim_semana_feriado")
    g = g.query("estacao == @season and tipo_dia_g == @td")
    prof = g.groupby("hora").agg(
        D_global=("D_global", "mean"),
        gs=("gs", "mean"),
        gr=("gr", "mean"),
        gh=("gh", "mean"),
        g_th=("g_th_obs", "mean"),
        g_nuc=("g_nuc_obs", "mean"),
        cmo=("cmo_h", "mean"),
    ).reindex(range(24))
    return prof

def solve_dispatch_day(D, g_ren, K_T, K_H, c1, c2, rho_up, rho_dn,
                      pi_c=0.0, pi_u=3500.0, H_avail=None, energy_budget_H=None):
    T = len(D)
    n = cp.Variable(T, nonneg=True)
    h = cp.Variable(T, nonneg=True)
    c = cp.Variable(T, nonneg=True)
    u = cp.Variable(T, nonneg=True)

    cost = (c1 * cp.sum(n) + 0.5 * c2 * cp.sum_squares(n)
            + pi_c * cp.sum(c) + pi_u * cp.sum(u))

    bal = (g_ren + h + n) - (D + c - u) == 0
    cons = [
        bal, n <= K_T, h <= K_H,
        n[1:] - n[:-1] <= rho_up,
        n[:-1] - n[1:] <= rho_dn,
    ]
    if H_avail is not None:
        cons.append(h <= H_avail)
    if energy_budget_H is not None:
        cons.append(cp.sum(h) <= energy_budget_H)

    prob = cp.Problem(cp.Minimize(cost), cons)
    prob.solve(solver=cp.CLARABEL)
    lam = -bal.dual_value if bal.dual_value is not None else np.full(T, np.nan)
    return {
        "status": prob.status, "n": np.array(n.value), "h": np.array(h.value),
        "c": np.array(c.value), "u": np.array(u.value),
        "lam": np.array(lam), "cost": float(prob.value),
    }

rows = []
prof_records = {}
for sub in SUBSYS_ORDER:
    for est in SEASON_ORDER:
        for td in ["util", "fim_semana_feriado"]:
            prof = typical_day(panel, sub, est, td)
            if prof["D_global"].isna().all():
                continue
            prof = prof.ffill().bfill().fillna(0.0)
            D     = prof["D_global"].values
            g_ren = (prof["gs"] + prof["gr"]).values
            K_T_l = float(K_inst.loc[sub, "K_t"])
            K_H_l = float(K_inst.loc[sub, "K_h"])
            d_th = panel.loc[panel["id_subsistema"] == sub, "g_th_obs"].diff().abs()
            rho_l = float(np.nanquantile(d_th.dropna(), 0.99))
            res = solve_dispatch_day(
                D, g_ren, K_T_l, K_H_l,
                c1=PARAMS["c1_RS_per_MWh"], c2=PARAMS["c2_RS_per_MWh2"],
                rho_up=rho_l, rho_dn=rho_l,
            )
            prof_records[(sub, est, td)] = (prof, res)
            rows.append({
                "id_subsistema": sub, "estacao": est, "tipo_dia": td,
                "status": res["status"], "cost": res["cost"],
                "G_th": float(res["n"].sum()),
                "G_h":  float(res["h"].sum()),
                "curt": float(res["c"].sum()),
                "def":  float(res["u"].sum()),
                "lam_mean": float(np.nanmean(res["lam"])),
                "lam_p95":  float(np.nanpercentile(res["lam"], 95)),
            })

despacho_central = pd.DataFrame(rows)
despacho_central.to_csv(TABLES / "despacho_centralizado_dias_tipicos.csv",
                        index=False, float_format="%.2f")
print("Despacho centralizado por (subsistema x estacao x tipo_dia):")
display(despacho_central.head(12))
''')


md(r'''
## 12. Modelo situacional finito (Seção 9) — Figuras 9 e 16

Sistema 1‑nó SIN com:
- Solar 1 (incumbente, capacidade `K_1`)
- Solar 2 (entrante, capacidade `K_2 = κ · K_1`)
- Eólica como exógena `G^w_t`
- Térmica e Hidrelétrica com reservatório

Horizonte 192 h (8 blocos de 24h: 4 estações × 2 tipos de dia), com
balanço hídrico contínuo entre blocos.
''')

code(r'''
agg = panel.groupby("din_instante", as_index=False).agg(
    D_global=("D_global", "sum"),
    gs=("gs", "sum"),
    gr=("gr", "sum"),
    gh=("gh", "sum"),
    g_th=("g_th_obs", "sum"),
    g_nuc=("g_nuc_obs", "sum"),
)
agg = add_calendar(agg)
agg["tipo_dia_g"] = np.where(agg["tipo_dia"] == "util", "util", "fim_semana_feriado")
typ_sin = (agg.groupby(["estacao", "tipo_dia_g", "hora"], observed=False)
              .agg(D=("D_global", "mean"),
                   gs=("gs", "mean"),
                   gr=("gr", "mean"),
                   gh=("gh", "mean"),
                   g_th=("g_th", "mean"),
                   g_nuc=("g_nuc", "mean"))
              .reindex(pd.MultiIndex.from_product(
                  [SEASON_ORDER, ["util", "fim_semana_feriado"], range(24)],
                  names=["estacao", "tipo_dia_g", "hora"])))
typ_sin = typ_sin.ffill().bfill().reset_index()

afl_season = (hydro.assign(estacao=lambda d: d["din_instante"].dt.month.map(SOUTH_HEMI_SEASON))
                   .groupby("estacao")["A_MW"].mean())

K1 = float(K_inst.loc[:, "K_s"].sum())
KAPPAS = [0.0, 0.25, 0.5, 1.0, 2.0]

y_sin = panel.dropna(subset=["cmo_h"])["cmo_h"].values
x_sin = panel.dropna(subset=["cmo_h"])[["g_th_obs", "g_nuc_obs"]].sum(axis=1).values
A_mat = np.column_stack([np.ones_like(x_sin), x_sin])
coef_sin, *_ = np.linalg.lstsq(A_mat, y_sin, rcond=None)
c1_sin, c2_sin = max(float(coef_sin[0]), 50.0), max(float(coef_sin[1]), 1e-6)

K_T_SIN = float(np.nanquantile(agg["g_th"], 0.995))
K_H_SIN = float(np.nanquantile(agg["gh"],   0.995))

# Geração hidro média observada por (estacao, tipo_dia, hora) — usaremos como
# limite hidráulico realista (capacidade instalada inclui restrições não
# modeladas: defluência ecológica, mínimo operativo, irrigação, etc.)
H_obs_block = typ_sin["gh"].values  # já alinhado em SEASON × tipo × hora
# Energia diária média de hidro por bloco (24h cada)
H_energy_daily = float(typ_sin.groupby(["estacao", "tipo_dia_g"])["gh"].sum().mean())

HOURS = 8 * 24

def build_typical_year_8days(typ_sin):
    rows = []
    for est in SEASON_ORDER:
        for td in ["util", "fim_semana_feriado"]:
            d = typ_sin.query("estacao == @est and tipo_dia_g == @td").sort_values("hora")
            d = d.copy()
            d["bloco"] = f"{est}-{td}"
            rows.append(d)
    return pd.concat(rows, ignore_index=True)

template = build_typical_year_8days(typ_sin)
K_s_obs_total = float(K_inst.loc[:, "K_s"].sum())
K_w_obs_total = float(K_inst.loc[:, "K_w"].sum())
template["a_s"] = template["gs"] / max(K_s_obs_total, 1e-6)
template["a_w"] = template["gr"] / max(K_w_obs_total, 1e-6)

D_t   = template["D"].values
a_s_t = template["a_s"].clip(0, 1.5).values
a_w_t = template["a_w"].clip(0, 1.5).values
G_w_t = a_w_t * K_w_obs_total

A_h_list = []
for est in SEASON_ORDER:
    for td in ["util", "fim_semana_feriado"]:
        Av = float(afl_season.get(est, np.nan))
        if not np.isfinite(Av):
            Av = float(afl_season.mean())
        A_h_list.append(np.full(24, Av))
A_h_base = np.concatenate(A_h_list)

# Volume armazenável proxy: 30 dias de afluência média (realista para o SIN
# em escala plurimensal, mas reduzido para que o modelo de 8 dias tenha trade-off
# intertemporal não-degenerado).
V_max_energy = float(afl_season.mean() * 24 * 30)
V_init  = 0.55 * V_max_energy
V_target = 0.55 * V_max_energy

# Valor terminal da água: calibrado para que seja vinculante na ordem
# de grandeza do CMO médio observado.
theta_water = float(np.nanmean(panel["cmo_h"].dropna()))   # ~ 165 R$/MWh

# Aumentamos c2 (componente quadrático) para evitar lambda colapsado: a
# regressão CMO ~ G_th dá c2 muito pequeno porque a correlação intra-dia é
# baixa. Aqui forçamos curvatura mínima compatível com elasticidade térmica.
c2_sin_eff = max(c2_sin, 5e-3)

# Penalidade simbólica em curtailment para evitar despacho gratuito
PI_CURT = 30.0

def solve_situational(K1, kappa, A_h=None, theta=None, H_block_cap=None):
    if A_h is None:
        A_h = A_h_base
    if theta is None:
        theta = theta_water
    if H_block_cap is None:
        # limite hidráulico horário = mín(capacidade instalada, 1.2× geração observada)
        Hb = np.nan_to_num(H_obs_block, nan=0.0)
        H_block_cap = np.minimum(K_H_SIN, 1.2 * Hb)
    # Sanitização: se K1 for NaN/Inf/negativo, retornamos falha controlada
    if not np.isfinite(K1) or K1 < 0:
        return {"K1": K1, "K2": kappa * K1, "kappa": kappa, "status": "infeasible_input",
                "n": np.full(HOURS, np.nan), "h": np.full(HOURS, np.nan),
                "c": np.full(HOURS, np.nan), "u": np.full(HOURS, np.nan),
                "s": np.full(HOURS, np.nan), "V": np.full(HOURS+1, np.nan),
                "lam": np.full(HOURS, np.nan), "mu_water": np.full(HOURS, np.nan),
                "g1": np.full(HOURS, np.nan), "g2": np.full(HOURS, np.nan),
                "cost": np.nan}
    K2 = kappa * K1
    T  = HOURS
    n  = cp.Variable(T, nonneg=True)
    h  = cp.Variable(T, nonneg=True)
    c  = cp.Variable(T, nonneg=True)
    u  = cp.Variable(T, nonneg=True)
    s  = cp.Variable(T, nonneg=True)
    V  = cp.Variable(T + 1, nonneg=True)

    g1 = K1 * a_s_t
    g2 = K2 * a_s_t

    cost = (c1_sin * cp.sum(n) + 0.5 * c2_sin_eff * cp.sum_squares(n)
            + PARAMS["pi_u_RS_per_MWh"] * cp.sum(u)
            + PI_CURT * cp.sum(c)
            - theta * V[-1])

    bal = g1 + g2 + G_w_t + h + n - (D_t + c - u) == 0
    res = V[1:] == V[:-1] + A_h - h - s

    cons = [
        bal, res,
        V[0] == V_init,
        V <= V_max_energy,
        V[-1] >= V_target - 0.05 * V_max_energy,
        n <= K_T_SIN,
        h <= H_block_cap,           # limite hidráulico horário realista
        n[1:] - n[:-1] <= rho_up,
        n[:-1] - n[1:] <= rho_dn,
    ]
    prob = cp.Problem(cp.Minimize(cost), cons)
    prob.solve(solver=cp.CLARABEL)

    lam = -bal.dual_value if bal.dual_value is not None else np.full(T, np.nan)
    mu  = res.dual_value if res.dual_value is not None else np.full(T, np.nan)
    return {
        "K1": K1, "K2": K2, "kappa": kappa, "status": prob.status,
        "n": np.array(n.value), "h": np.array(h.value),
        "c": np.array(c.value), "u": np.array(u.value),
        "s": np.array(s.value), "V": np.array(V.value),
        "lam": np.array(lam), "mu_water": np.array(mu),
        "g1": K1 * a_s_t, "g2": K2 * a_s_t,
        "cost": float(prob.value),
    }

results_sit = {kappa: solve_situational(K1, kappa) for kappa in KAPPAS}
print("Status por kappa:", {k: v["status"] for k, v in results_sit.items()})

rb = results_sit[0.5]
fig, axes = plt.subplots(4, 2, figsize=(14, 14), sharex=False)
labels_block = template[["estacao", "tipo_dia_g"]].drop_duplicates().reset_index(drop=True)
for i in range(8):
    ax = axes[i//2, i%2]
    sl = slice(i*24, (i+1)*24)
    hours = np.arange(24)
    ax.fill_between(hours, 0, rb["g1"][sl]/1000, color="#fdae61", alpha=.7, label="Solar 1")
    ax.fill_between(hours, rb["g1"][sl]/1000, (rb["g1"][sl]+rb["g2"][sl])/1000,
                    color="#fee08b", alpha=.7, label="Solar 2")
    base = (rb["g1"][sl]+rb["g2"][sl])
    ax.fill_between(hours, base/1000, (base + G_w_t[sl])/1000, color="#74add1", alpha=.6, label="Eolica")
    base2 = base + G_w_t[sl]
    ax.fill_between(hours, base2/1000, (base2 + rb["h"][sl])/1000, color="#abd9e9", alpha=.7, label="Hidro")
    base3 = base2 + rb["h"][sl]
    ax.fill_between(hours, base3/1000, (base3 + rb["n"][sl])/1000, color="#f46d43", alpha=.7, label="Termica")
    ax.plot(hours, D_t[sl]/1000, "k-", lw=2, label="Demanda")
    ax2 = ax.twinx()
    ax2.plot(hours, rb["lam"][sl], "g--", lw=1.5, label="lambda_t")
    ax2.set_ylabel("lambda R$/MWh", color="g")
    ax.set_xticks(range(0, 24, 3))
    bl = labels_block.iloc[i]
    ax.set_title(f"Bloco {bl['estacao']} - {bl['tipo_dia_g']}")
    ax.set_ylabel("GW")
    ax.grid(alpha=.3)
    if i == 0:
        ax.legend(loc="upper left", fontsize=8, ncols=2)
fig.suptitle("Figura 9 / 16 - Despacho situacional (kappa=0.5): solar1, solar2, eolica, hidro, termica e lambda_t",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig09_16_despacho_situacional.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 13. Canibalização e sensibilidade (Figuras 10 e 17)

Para `κ = K_2/K_1 ∈ {0, 0.25, 0.5, 1, 2}`:
- `R_cap_i` (preço médio capturado)
- `Π_i = Σ λ_t · g_i,t − I_i(K_i)` (lucro)
- `C_curt`, `G_therm`, `R_max`

Custo de investimento: `I_i(K) = F + q · K`, com `q = 150 R$/kW·ano` anualizado
(CAPEX ~ 1 BRL/W, vida 25y, WACC 8% → ~150 R$/MWh × 8760 h).
''')

code(r'''
HOURS_YEAR = 8760
F_i = 0.0
# CAPEX solar realista: ~4 BRL/Wp = 4 milhões R$/MW
# Anualizado em 25y a WACC 8% → fator 0.094 → 376 mil R$/MW/ano
# Para 8 dias típicos (192h): q_i = 376000 * 192/8760 ≈ 8240 R$/MW
q_i_annual = 376000.0
q_i = q_i_annual * (HOURS / HOURS_YEAR)
eta_i = 0.0

def Inv(K):
    return F_i + q_i * K + 0.5 * eta_i * K**2

rows = []
for kappa in KAPPAS:
    r = results_sit[kappa]
    g1 = r["g1"]; g2 = r["g2"]; lam = r["lam"]
    E1, E2 = float(g1.sum()), float(g2.sum())
    R1 = float((lam * g1).sum())
    R2 = float((lam * g2).sum())
    Rcap1 = R1 / E1 if E1 > 0 else np.nan
    Rcap2 = R2 / E2 if E2 > 0 else np.nan
    Pi1 = R1 - Inv(r["K1"])
    Pi2 = R2 - Inv(r["K2"])
    rows.append({
        "kappa": kappa, "K1": r["K1"], "K2": r["K2"],
        "R1": R1, "R2": R2,
        "R_cap_1": Rcap1, "R_cap_2": Rcap2,
        "Pi_1": Pi1, "Pi_2": Pi2,
        "G_therm": float(r["n"].sum()),
        "C_curt": float(r["c"].sum()),
        "Deficit": float(r["u"].sum()),
    })
canib = pd.DataFrame(rows)
canib.to_csv(TABLES / "canibalizacao_kappa.csv", index=False, float_format="%.2f")
display(canib.round(2))

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

ax = axes[0]
ax.plot(canib["kappa"], canib["R_cap_1"], "o-", color="#fdae61", lw=2, label="Solar 1 (incumbente)")
ax.plot(canib["kappa"], canib["R_cap_2"], "s-", color="#d73027", lw=2, label="Solar 2 (entrante)")
ax.set_xlabel("kappa = K_2/K_1")
ax.set_ylabel("R_cap (R$/MWh)")
ax.set_title("(a) Receita capturada por agente"); ax.legend(); ax.grid(alpha=.3)

ax = axes[1]
ax.plot(canib["kappa"], canib["Pi_1"]/1e6, "o-", color="#fdae61", lw=2, label="Pi_1 Solar 1")
ax.plot(canib["kappa"], canib["Pi_2"]/1e6, "s-", color="#d73027", lw=2, label="Pi_2 Solar 2")
ax.axhline(0, color="black", lw=1)
ax.set_xlabel("kappa"); ax.set_ylabel("Lucro Pi (R$ x10^6)")
ax.set_title("(b) Lucro privado dos produtores solares"); ax.legend(); ax.grid(alpha=.3)

ax = axes[2]
ax.bar(np.arange(len(KAPPAS)) - 0.2, canib["G_therm"]/1e6, width=0.4, color="#f46d43", label="G_therm (TWh)")
ax.bar(np.arange(len(KAPPAS)) + 0.2, canib["C_curt"]/1e6,  width=0.4, color="#74add1", label="Curtailment (TWh)")
ax.set_xticks(np.arange(len(KAPPAS)))
ax.set_xticklabels([f"{k}" for k in KAPPAS])
ax.set_xlabel("kappa"); ax.set_ylabel("Energia (TWh equivalente)")
ax.set_title("(c) Termica e curtailment vs kappa"); ax.legend(); ax.grid(alpha=.3)

fig.suptitle("Figuras 10 / 17 - Canibalizacao: lucro, R_cap e operacao vs tamanho do entrante",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig10_17_canibalizacao.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 14. Valor marginal da água (H3 — Figura 18)
''')

code(r'''
SCEN_H = {"seca (-30%)": 0.7, "base": 1.0, "umida (+30%)": 1.3}
hidro_results = {name: solve_situational(K1, 0.5, A_h=A_h_base * sc)
                 for name, sc in SCEN_H.items()}

rows = []
for name, r in hidro_results.items():
    rows.append({
        "cenario": name,
        "G_therm (TWh)": r["n"].sum()/1e6,
        "G_h (TWh)":     r["h"].sum()/1e6,
        "Spill (TWh)":   r["s"].sum()/1e6,
        "Vend (MWh)":    r["V"][-1]/1e3,
        "mu_w medio":    float(np.nanmean(r["mu_water"])),
        "lambda medio":  float(np.nanmean(r["lam"])),
    })
hydro_summary = pd.DataFrame(rows)
hydro_summary.to_csv(TABLES / "valor_agua_cenarios.csv", index=False, float_format="%.3f")
display(hydro_summary.round(3))

fig, axes = plt.subplots(2, 2, figsize=(14, 9))

ax = axes[0, 0]
for (name, r), ls in zip(hidro_results.items(), ["-", "--", "-."]):
    ax.plot(r["V"]/1e3, ls=ls, label=name, lw=2)
ax.set_title("(a) Estado de armazenamento V_t")
ax.set_xlabel("hora"); ax.set_ylabel("GWh-equiv.")
ax.legend(); ax.grid(alpha=.3)

ax = axes[0, 1]
for (name, r), ls in zip(hidro_results.items(), ["-", "--", "-."]):
    ax.plot(np.abs(r["mu_water"]), ls=ls, label=name, lw=1.5)
ax.set_title("(b) Valor marginal da agua |mu_t^w|")
ax.set_xlabel("hora"); ax.set_ylabel("R$/MWh"); ax.legend(); ax.grid(alpha=.3)

ax = axes[1, 0]
ax.bar(hydro_summary["cenario"], hydro_summary["G_therm (TWh)"], color="#f46d43")
ax.set_title("(c) Geracao termica total")
ax.set_ylabel("TWh"); ax.grid(alpha=.3, axis="y")

ax = axes[1, 1]
ax.bar(hydro_summary["cenario"], hydro_summary["Spill (TWh)"], color="#9467bd")
ax.set_title("(d) Vertimento total")
ax.set_ylabel("TWh"); ax.grid(alpha=.3, axis="y")

fig.suptitle("Figura 18 - Valor marginal da agua por estado hidrologico (H3)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig18_valor_da_agua.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 15. Mean Field Games locacional (Seções 10–12)

Versão simplificada com 4 localizações (N, NE, SE, S), 4 tipos
tecnológicos e horizonte T = 24 h. Iteração de ponto fixo:
- λ⁰ ← CMO médio
- Resolver clearing locacional → λ̂ atualizado
- λ ← (1−ω) λ + ω λ̂   (relaxação)
- Parar quando ‖λ − λ_prev‖∞ < ε
''')

code(r'''
typ_loc = (panel.groupby(["id_subsistema", "hora"])
                .agg(D=("D_global", "mean"),
                     gs=("gs", "mean"),
                     gr=("gr", "mean"),
                     gh=("gh", "mean"),
                     g_th=("g_th_obs", "mean"))
                .reset_index())

intc_path = paths.interc_interno_path
intc_raw = pd.read_csv(intc_path, sep=";", parse_dates=["din_instante"])
intc_raw = intc_raw[intc_raw["din_instante"].between(DATE_MIN, DATE_MAX)]
intc_raw["o"] = intc_raw["nom_subsistema_origem"].astype(str).map(canonical_subsys)
intc_raw["d"] = intc_raw["nom_subsistema_destino"].astype(str).map(canonical_subsys)
flow_lim = (intc_raw.groupby(["o", "d"])["val_intercambiomwmed"]
            .apply(lambda x: float(np.nanquantile(x.abs(), 0.99))).reset_index()
            .rename(columns={"val_intercambiomwmed": "F_max"}))
display(flow_lim)
''')

code(r'''
ALL_LINKS = [("N", "NE"), ("N", "SE"), ("NE", "SE"), ("SE", "S")]
F_MAX = {(o, d): float(flow_lim.query("o == @o and d == @d")["F_max"].iloc[0])
         if len(flow_lim.query("o == @o and d == @d")) else 0.0
         for (o, d) in ALL_LINKS}
print("Limites de fluxo F_max (MW):", F_MAX)


def clearing_locacional(D_lt, G_ren_lt, K_T_l, K_H_l, c1_l, c2_l, H_avail_lt,
                       pi_u=3500.0):
    T = 24
    locs = SUBSYS_ORDER
    n = {l: cp.Variable(T, nonneg=True) for l in locs}
    h = {l: cp.Variable(T, nonneg=True) for l in locs}
    c = {l: cp.Variable(T, nonneg=True) for l in locs}
    u = {l: cp.Variable(T, nonneg=True) for l in locs}
    F = {(o, d): cp.Variable(T) for (o, d) in ALL_LINKS}

    cost = 0
    # theta_hydro = custo de oportunidade da água (valor sombra plurimensal)
    # calibrado para que o preço-sombra fique na faixa observada do CMO
    theta_hydro = 80.0
    for l in locs:
        # Forçamos curvatura mínima para evitar lambda degenerado
        c2_eff = max(c2_l[l], 5e-3)
        cost += c1_l[l] * cp.sum(n[l]) + 0.5 * c2_eff * cp.sum_squares(n[l]) \
                + theta_hydro * cp.sum(h[l]) \
                + pi_u * cp.sum(u[l]) + 30.0 * cp.sum(c[l])

    cons = []
    bal = {}
    for l in locs:
        inflow  = sum( F[(o, d)] for (o, d) in ALL_LINKS if d == l )
        outflow = sum( F[(o, d)] for (o, d) in ALL_LINKS if o == l )
        bal[l] = (G_ren_lt[l] + h[l] + n[l] + inflow - outflow) == (D_lt[l] + c[l] - u[l])
        cons.append(bal[l])
        cons.append(n[l] <= K_T_l[l])
        cons.append(h[l] <= H_avail_lt[l])
        # orçamento diário de energia hidráulica = ~95% do dia típico observado
        # (dá flexibilidade intra‑dia mas mantém balanço energético plausível)
        cons.append(cp.sum(h[l]) <= H_avail_lt[l].sum() * 0.95)
    for (o, d) in ALL_LINKS:
        cons.append(F[(o, d)] <=  F_MAX[(o, d)])
        cons.append(F[(o, d)] >= -F_MAX[(o, d)])
    prob = cp.Problem(cp.Minimize(cost), cons)
    prob.solve(solver=cp.CLARABEL)

    lam_l = {l: -bal[l].dual_value if bal[l].dual_value is not None else np.full(T, np.nan)
             for l in locs}
    return {
        "n": {l: n[l].value for l in locs},
        "h": {l: h[l].value for l in locs},
        "c": {l: c[l].value for l in locs},
        "u": {l: u[l].value for l in locs},
        "F": {k: F[k].value for k in F},
        "lam": lam_l,
        "cost": float(prob.value),
        "status": prob.status,
    }

D_lt    = {l: typ_loc.query("id_subsistema == @l").sort_values("hora")["D"].values for l in SUBSYS_ORDER}
Gren_lt = {l: (typ_loc.query("id_subsistema == @l").sort_values("hora")[["gs", "gr"]].sum(axis=1).values)
           for l in SUBSYS_ORDER}
# Disponibilidade hidráulica horária = 1.5× geração média observada
# (limita o quanto a hidro pode "subir" além do despacho típico)
Havail_lt = {l: 1.5 * typ_loc.query("id_subsistema == @l").sort_values("hora")["gh"].values
             for l in SUBSYS_ORDER}
K_T_loc = {l: float(K_inst.loc[l, "K_t"]) for l in SUBSYS_ORDER}
K_H_loc = {l: float(K_inst.loc[l, "K_h"]) for l in SUBSYS_ORDER}

c1_loc, c2_loc = {}, {}
for l in SUBSYS_ORDER:
    g = panel.dropna(subset=["cmo_h"])[panel["id_subsistema"] == l].copy()
    g["g_disp"] = g["g_th_obs"] + g["g_nuc_obs"]
    if len(g) > 100:
        A_loc = np.column_stack([np.ones_like(g["g_disp"]), g["g_disp"]])
        coef_loc, *_ = np.linalg.lstsq(A_loc, g["cmo_h"].values, rcond=None)
        c1_loc[l] = max(float(coef_loc[0]), 30.0)
        c2_loc[l] = max(float(coef_loc[1]), 1e-6)
    else:
        c1_loc[l] = 80.0; c2_loc[l] = 1e-4

print("Custos termicos locais:")
display(pd.DataFrame({"c1": c1_loc, "c2": c2_loc}))
''')


md(r'''
### Iteração MFG
''')

code(r'''
N_ITER = 30
OMEGA  = 0.6
EPS    = 1e-2

lam_obs = (panel.groupby(["id_subsistema", "hora"])["cmo_h"].mean().unstack().reindex(SUBSYS_ORDER))
lam = {l: lam_obs.loc[l].values for l in SUBSYS_ORDER}
hist = {"err": [], "cost": []}

for k in range(N_ITER):
    res = clearing_locacional(D_lt, Gren_lt, K_T_loc, K_H_loc, c1_loc, c2_loc, Havail_lt)
    lam_new = res["lam"]
    err = max(np.nanmax(np.abs(lam_new[l] - lam[l])) for l in SUBSYS_ORDER)
    hist["err"].append(err)
    hist["cost"].append(res["cost"])
    for l in SUBSYS_ORDER:
        lam[l] = (1 - OMEGA) * lam[l] + OMEGA * lam_new[l]
    if err < EPS:
        print(f"convergiu em {k+1} iteracoes | erro = {err:.4f}")
        break

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(hist["err"], "o-")
axes[0].set_yscale("log"); axes[0].set_xlabel("iteracao"); axes[0].set_ylabel("||lambda_{k+1} - lambda_k||_inf")
axes[0].set_title("(a) Erro maximo de preco")
axes[0].grid(alpha=.3)

axes[1].plot(hist["cost"], "s-", color="#d73027")
axes[1].set_xlabel("iteracao"); axes[1].set_ylabel("Custo total (R$)")
axes[1].set_title("(b) Custo total por iteracao")
axes[1].grid(alpha=.3)

fig.suptitle("Figura 12 - Convergencia do algoritmo MFG locacional",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig12_convergencia_mfg.png", bbox_inches="tight")
plt.show()
''')


md(r'''
### Figura 19 — Equilíbrio MFG locacional
''')

code(r'''
fig = plt.figure(figsize=(16, 12))
gs  = fig.add_gridspec(3, 2)

ax = fig.add_subplot(gs[0, :])
for l in SUBSYS_ORDER:
    ax.plot(np.arange(24), res["lam"][l], "o-", color=SUB_COLOR[l], label=l)
ax.set_title("(a) Preco locacional lambda_{l,t} no equilibrio MFG")
ax.set_xlabel("hora"); ax.set_ylabel("R$/MWh")
ax.legend(ncols=4); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[1, 0])
for l in SUBSYS_ORDER:
    ax.plot(np.arange(24), res["n"][l]/1000, "-", color=SUB_COLOR[l], label=f"{l} termica")
ax.set_title("(b) Geracao termica n_{l,t}")
ax.set_xlabel("hora"); ax.set_ylabel("GW"); ax.legend(); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[1, 1])
for l in SUBSYS_ORDER:
    ax.plot(np.arange(24), res["h"][l]/1000, "-", color=SUB_COLOR[l], label=f"{l} hidro")
ax.set_title("(c) Geracao hidro h_{l,t}")
ax.set_xlabel("hora"); ax.set_ylabel("GW"); ax.legend(); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[2, 0])
for l in SUBSYS_ORDER:
    ax.plot(np.arange(24), res["c"][l], "-", color=SUB_COLOR[l], label=f"{l}")
ax.set_title("(d) Curtailment c_{l,t}"); ax.set_xlabel("hora"); ax.set_ylabel("MW"); ax.legend(); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[2, 1])
for (o, d) in ALL_LINKS:
    ax.plot(np.arange(24), res["F"][(o, d)]/1000, "-", lw=2, label=f"{o}->{d}")
    ax.axhline(F_MAX[(o, d)]/1000, color="gray", ls="--", lw=0.5)
ax.set_title("(e) Fluxos inter-subsistema F_{lk,t}")
ax.set_xlabel("hora"); ax.set_ylabel("GW"); ax.legend(); ax.grid(alpha=.3)

fig.suptitle("Figura 19 - Equilibrio MFG locacional (24 h, dia tipico anual)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig19_mfg_locacional.png", bbox_inches="tight")
plt.show()
''')


md(r'''
### Figura 11 — Diagrama conceitual do ponto fixo MFG
''')

code(r'''
fig, ax = plt.subplots(figsize=(11, 5.5))
ax.axis("off")
box = dict(boxstyle="round,pad=0.5", facecolor="#fff5b8", edgecolor="black")
boxd = dict(boxstyle="round,pad=0.5", facecolor="#d3e8ff", edgecolor="black")
boxc = dict(boxstyle="round,pad=0.5", facecolor="#d8f0c8", edgecolor="black")

ax.text(0.05, 0.5, "lambda^k_{l,t}\n(preco locacional)", ha="center", va="center", bbox=boxd, fontsize=11)
ax.annotate("", xy=(0.25, 0.5), xytext=(0.15, 0.5),
            arrowprops=dict(arrowstyle="->", lw=2))
ax.text(0.35, 0.5, "pi^tau_{l,t}\n(controle otimo\npor tipo tau)", ha="center", va="center", bbox=box, fontsize=11)
ax.annotate("", xy=(0.55, 0.5), xytext=(0.45, 0.5),
            arrowprops=dict(arrowstyle="->", lw=2))
ax.text(0.65, 0.5, "m^tau_{l,t}\n(distribuicao\npopulacional)", ha="center", va="center", bbox=box, fontsize=11)
ax.annotate("", xy=(0.85, 0.5), xytext=(0.75, 0.5),
            arrowprops=dict(arrowstyle="->", lw=2))
ax.text(0.95, 0.5, "G_{l,t}\n(oferta agregada)", ha="center", va="center", bbox=boxc, fontsize=11)
ax.annotate("clearing locacional ->", xy=(0.05, 0.55), xytext=(0.95, 0.55),
            arrowprops=dict(arrowstyle="->", lw=2, color="red"),
            ha="center", color="red", fontsize=10)
ax.set_title("Figura 11 - Diagrama de ponto fixo MFG locacional",
             fontsize=12, fontweight="bold")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
fig.tight_layout()
fig.savefig(FIGS / "fig11_diagrama_mfg.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 16. Comparação ótimo social × equilíbrio privado (H5)
''')

code(r'''
def Q_of_K(Ks_total, theta=theta_water):
    K1_local = Ks_total
    r = solve_situational(K1_local, kappa=0.0, theta=theta)
    return r["cost"], r

K_grid = np.array([0.3, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0]) * K1
soc_rows = []
for K in K_grid:
    cost_op, r = Q_of_K(K)
    soc_rows.append({"K": K, "cost_op": cost_op,
                     "Inv": Inv(K), "total": cost_op + Inv(K)})
soc = pd.DataFrame(soc_rows)
soc.to_csv(TABLES / "expansao_social.csv", index=False, float_format="%.2f")
display(soc.round(0))

K_soc = float(soc.loc[soc["total"].idxmin(), "K"])
print(f"\nCapacidade socialmente otima estimada: K_soc = {K_soc/1000:.1f} GW")
print(f"Capacidade observada (proxy p99):       K_obs = {K1/1000:.1f} GW")

def best_response_priv(K_init, n_iter=10, q_i=q_i, eta_i=eta_i):
    K = float(K_init)
    K_min = 0.1 * K_init   # piso
    K_max = 5.0 * K_init   # teto
    for k in range(n_iter):
        K = float(np.clip(K, K_min, K_max))
        r = solve_situational(K, kappa=0.0)
        lam = r["lam"]
        if np.all(np.isnan(lam)) or np.isnan(lam).any():
            print(f"  iter {k}: lambda com NaN, parando.")
            break
        rev_per_K = float((lam * a_s_t).sum())
        if not np.isfinite(rev_per_K) or rev_per_K <= 0:
            break
        K_new = K * rev_per_K / max(q_i, 1e-6)
        K_new = float(np.clip(K_new, K_min, K_max))
        if abs(K_new - K) / max(K, 1.0) < 0.01:
            break
        # relaxação para estabilidade
        K = 0.5 * K + 0.5 * K_new
    return K

K_priv = best_response_priv(K1)
print(f"Capacidade privada de equilibrio:        K_priv = {K_priv/1000:.1f} GW")

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(soc["K"]/1000, soc["cost_op"]/1e6, "o-", label="Custo operacional Q(K)", color="#1f77b4")
ax.plot(soc["K"]/1000, soc["Inv"]/1e6,     "s-", label="Custo investimento I(K)", color="#d62728")
ax.plot(soc["K"]/1000, soc["total"]/1e6,   "^-", label="Total Q(K)+I(K)", color="black", lw=2)
ax.axvline(K_soc/1000,  color="green", ls="--", label=f"K_soc = {K_soc/1000:.1f} GW")
ax.axvline(K_priv/1000, color="orange", ls=":",  label=f"K_priv = {K_priv/1000:.1f} GW")
ax.set_xlabel("K_solar (GW)"); ax.set_ylabel("R$ x10^6")
ax.set_title("H5 - Comparacao otimo social x equilibrio privado")
ax.legend(); ax.grid(alpha=.3)
fig.tight_layout()
fig.savefig(FIGS / "fig_H5_social_vs_privado.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 17. Apêndice B — Métricas consolidadas
''')

code(r'''
metrics = panel.groupby("id_subsistema").agg(
    PLD_avg=("cmo_h", "mean"),
    PLD_p95=("cmo_h", lambda x: float(np.nanpercentile(x.dropna(), 95))),
    D_mean=("D_global", "mean"),
    D_max=("D_global", "max"),
    D_net_mean=("D_net", "mean"),
    G_solar_total=("gs", "sum"),
    G_wind_total=("gr", "sum"),
    G_hydro_total=("gh", "sum"),
    G_thermal_total=("g_th_obs", "sum"),
    G_nuclear_total=("g_nuc_obs", "sum"),
    Ramp_max_pos=("dD_res", "max"),
    Ramp_max_neg=("dD_res", "min"),
).reindex(SUBSYS_ORDER)

fc_solar = fcap_df.query("fonte == 'solar'").set_index("id_subsistema")["F_capture"].reindex(SUBSYS_ORDER)
fc_wind  = fcap_df.query("fonte == 'wind' ").set_index("id_subsistema")["F_capture"].reindex(SUBSYS_ORDER)
metrics["F_capture_solar"] = fc_solar
metrics["F_capture_wind"]  = fc_wind

metrics["frac_horas_excedente"] = (
    panel.assign(exc=lambda d: (d["gs"]+d["gr"] > 0.7*d["D_global"]).astype(int))
         .groupby("id_subsistema")["exc"].mean()
         .reindex(SUBSYS_ORDER)
)
metrics["taxa_curtailment_estimada"] = (
    panel.assign(disp=lambda d: d["gs"]+d["gr"]+d["g_th_obs"]+d["g_nuc_obs"]+d["gh"])
         .assign(excesso=lambda d: np.maximum(d["disp"]-d["D_global"], 0))
         .groupby("id_subsistema").apply(lambda g: g["excesso"].sum() / g["disp"].sum(), include_groups=False)
         .reindex(SUBSYS_ORDER)
)

metrics_out = metrics.round(2)
metrics_out.to_csv(TABLES / "apendice_B_metricas_consolidadas.csv")
print("APENDICE B - Metricas consolidadas (2025, ate 2025-10-31):")
display(metrics_out)
''')


# =========================================================================
# SEÇÃO 18 — Modelagem física da geração solar (latitude)
# =========================================================================
md(r'''
## 18. Modelagem física da geração solar (latitude · irradiância clear‑sky)

O documento usa o fator de capacidade observado `a^s_{ℓ,t}`. Para entender
**por que** `a^s` difere entre subsistemas, modelamos a irradiância
teórica usando o modelo *clear‑sky* dependente da latitude e do dia do ano:

$$
\begin{aligned}
\delta(d) &= 23.45^\circ \sin\!\left(\tfrac{2\pi(d-81)}{365}\right) \quad \text{(declinação solar)} \\
H(t) &= 15^\circ\,(t - 12) \quad \text{(ângulo horário em graus)} \\
\cos\theta_z &= \sin\phi\sin\delta + \cos\phi\cos\delta\cos H \quad \text{(ângulo zenital)} \\
I_{\text{horiz}}(t) &= I_{\text{sc}}\,(1-\rho_{\text{cloud}})\,\max(0, \cos\theta_z) \quad (\text{kW/m}^2)
\end{aligned}
$$

A integral anual `H_a(φ, ρ)` reproduz a queda esperada de irradiância
com `|lat|` crescente, modulada por nebulosidade regional. Os **centroides
operacionais** dos 4 subsistemas são:

| Subsistema | Centroide (lat, lon) | Cidades de referência |
| --- | --- | --- |
| N (Norte) | (−5°, −55°) | Manaus / Belém |
| NE (Nordeste) | (−9°, −40°) | Petrolina / Salvador |
| SE (Sudeste/CO) | (−18.5°, −47°) | Brasília / São Paulo |
| S (Sul) | (−28°, −52°) | Porto Alegre / Curitiba |
''')

code(r'''
SUBSYS_CENTROIDS = {
    "N":  (-5.0,  -55.0),
    "NE": (-9.0,  -40.0),
    "SE": (-18.5, -47.0),
    "S":  (-28.0, -52.0),
}
# Climatologia simplificada de cobertura média de nuvem (Atlas Brasileiro
# de Energia Solar 2017 — INPE/LABREN, ordens de grandeza)
SUBSYS_CLOUDINESS = {
    "N":  0.45,   # Amazônia: alta nebulosidade
    "NE": 0.15,   # semi-árido: céu claro
    "SE": 0.30,   # Cerrado/litoral: moderado
    "S":  0.35,   # frentes frias: moderado-alto
}

def declinacao(d):
    """Declinação solar (graus) para dia do ano d (1..365)."""
    return 23.45 * np.sin(np.deg2rad(360.0 * (d - 81) / 365.0))

def irradiancia_horizontal(lat, d, hora_solar, I_sc=1.36, cloud=0.0):
    """Irradiância (kW/m²) em superfície horizontal, modelo clear-sky reduzido.
    lat, declin em graus; hora_solar entre 0..24 (12 = meio-dia)."""
    delta = declinacao(d)
    H = 15.0 * (hora_solar - 12.0)
    cos_z = (np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(delta))
             + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(delta)) * np.cos(np.deg2rad(H)))
    return I_sc * (1 - cloud) * np.maximum(0.0, cos_z)

# Cálculo anual H_a por subsistema (resolução horária)
days = np.arange(1, 366)
hours = np.arange(0, 24)
DAYS, HRS = np.meshgrid(days, hours, indexing="ij")

lat_table = []
for sub, (lat, lon) in SUBSYS_CENTROIDS.items():
    cloud = SUBSYS_CLOUDINESS[sub]
    I = irradiancia_horizontal(lat, DAYS, HRS, cloud=cloud)
    H_diario = I.sum(axis=1)                # kWh/m²/dia
    H_anual  = H_diario.sum()                # kWh/m²/ano
    H_medio_pico = I.max(axis=1).mean()      # kW/m² (irradiância de pico médio)
    # CF teórico ≈ H_anual / (8760 * 1 kW/m²) — fração de horas-equivalentes
    CF_teorico = H_anual / 8760.0

    # CF observado (a_gs) — média anual
    g = panel[panel["id_subsistema"] == sub].dropna(subset=["a_gs"])
    CF_obs = float(g["a_gs"].mean()) if len(g) else np.nan

    lat_table.append({
        "subsistema": sub, "lat": lat, "lon": lon,
        "nebulosidade": cloud,
        "H_anual_kWh_m2": H_anual,
        "H_pico_kW_m2":   H_medio_pico,
        "CF_teorico":     CF_teorico,
        "CF_observado":   CF_obs,
        "razao_obs_teor": CF_obs / CF_teorico if CF_teorico > 0 else np.nan,
    })
lat_df = pd.DataFrame(lat_table)
lat_df.to_csv(TABLES / "latitude_irradiancia_solar.csv", index=False, float_format="%.4f")
print("Modelo de latitude (clear-sky):")
display(lat_df.round(3))

# Figura: perfil diário de irradiância para o equinócio (d=80, ~21/mar)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
for sub, (lat, _lon) in SUBSYS_CENTROIDS.items():
    cloud = SUBSYS_CLOUDINESS[sub]
    I = irradiancia_horizontal(lat, 80, hours, cloud=cloud)
    ax.plot(hours, I, lw=2, color=SUB_COLOR[sub], label=f"{sub} (lat={lat}°, nuvem={cloud:.0%})")
ax.set_xlabel("hora solar local"); ax.set_ylabel("Irradiância (kW/m²)")
ax.set_title("(a) Perfil clear-sky — equinócio (21 mar)")
ax.legend(); ax.grid(alpha=.3)

ax = axes[1]
sub_list = list(SUBSYS_CENTROIDS.keys())
x = np.arange(len(sub_list))
ax.bar(x - 0.2, lat_df["CF_teorico"],   width=0.4, label="CF teórico (clear-sky)", color="#1f77b4")
ax.bar(x + 0.2, lat_df["CF_observado"], width=0.4, label="CF observado 2025",       color="#fdae61")
ax.set_xticks(x); ax.set_xticklabels(sub_list)
ax.set_ylabel("Fator de capacidade")
ax.set_title("(b) CF teórico × observado por subsistema")
ax.legend(); ax.grid(alpha=.3, axis="y")

fig.suptitle("Figura 20 — Modelo físico de irradiância solar (latitude + nebulosidade)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig20_latitude_irradiancia.png", bbox_inches="tight")
plt.show()
''')


# =========================================================================
# SEÇÃO 19 — Mapas geográficos do SIN
# =========================================================================
md(r'''
## 19. Mapas geográficos do SIN

Visualizações esquemáticas do Brasil com os 4 subsistemas. Cada mapa
mostra uma variável-chave por região e a topologia dos corredores de
transmissão (fluxos máximos `F^{max}_{ℓk}`).

> **Nota.** Polígonos são esquemáticos (não usam shapefiles oficiais).
> Centroides correspondem às coordenadas reais.
''')

code(r'''
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm

# Polígonos esquemáticos dos 4 subsistemas (lon, lat)
SUBSYS_POLYGONS = {
    "N":  [(-74, -10), (-74,   5), (-50,  5), (-44, -2), (-44, -10)],
    "NE": [(-46, -18), (-46,  -2), (-34, -3), (-34, -17), (-40, -18)],
    "SE": [(-58, -25), (-58, -15), (-48, -10), (-44, -10), (-39, -20),
           (-40, -24), (-50, -25)],
    "S":  [(-58, -34), (-58, -25), (-50, -25), (-48, -29), (-50, -34)],
}

def draw_brazil_map(ax, value_per_subsys, cmap_name="YlOrRd", title="",
                     label_fmt=".1f", vmin=None, vmax=None, label_unit=""):
    """Desenha mapa esquemático do Brasil com cor por subsistema."""
    cmap = cm.get_cmap(cmap_name)
    vals = np.array([value_per_subsys.get(s, np.nan) for s in SUBSYS_ORDER])
    vmin = vmin if vmin is not None else np.nanmin(vals)
    vmax = vmax if vmax is not None else np.nanmax(vals)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    patches = []
    colors = []
    for sub in SUBSYS_ORDER:
        poly = SUBSYS_POLYGONS[sub]
        patches.append(MplPolygon(poly, closed=True))
        v = value_per_subsys.get(sub, np.nan)
        colors.append(cmap(norm(v)) if np.isfinite(v) else (0.9, 0.9, 0.9, 1.0))

    pc = PatchCollection(patches, facecolor=colors, edgecolor="black", linewidth=1.0, alpha=0.85)
    ax.add_collection(pc)

    # Centroides + valores
    for sub in SUBSYS_ORDER:
        lat, lon = SUBSYS_CENTROIDS[sub]
        v = value_per_subsys.get(sub, np.nan)
        ax.plot(lon, lat, "ko", markersize=5)
        ax.annotate(f"{sub}\n{v:{label_fmt}}{label_unit}",
                    (lon, lat), xytext=(0, -15), textcoords="offset points",
                    ha="center", fontsize=9, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))

    ax.set_xlim(-76, -32); ax.set_ylim(-36, 8)
    ax.set_aspect("equal")
    ax.set_xlabel("longitude (°)"); ax.set_ylabel("latitude (°)")
    ax.set_title(title)
    ax.grid(alpha=.2)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, shrink=.6, label=label_unit)


# Construir variáveis por subsistema
K_inst_dict = {s: float(K_inst.loc[s, "K_s"]) for s in SUBSYS_ORDER}
CF_solar    = {s: float(panel[panel["id_subsistema"] == s]["a_gs"].mean()) for s in SUBSYS_ORDER}
PLD_avg     = {s: float(panel[panel["id_subsistema"] == s]["cmo_h"].mean()) for s in SUBSYS_ORDER}
Fcap_solar  = {s: float(fcap_df.query("id_subsistema == @s and fonte == 'solar'")["F_capture"].iloc[0])
               for s in SUBSYS_ORDER}

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
draw_brazil_map(axes[0, 0], {s: v/1000 for s, v in K_inst_dict.items()},
                cmap_name="YlOrRd", title="(a) Capacidade solar instalada (proxy p99)",
                label_fmt=".1f", label_unit=" GW")
draw_brazil_map(axes[0, 1], CF_solar,
                cmap_name="YlOrRd", title="(b) Fator de capacidade solar médio",
                label_fmt=".2f", label_unit="")
draw_brazil_map(axes[1, 0], PLD_avg,
                cmap_name="Reds", title="(c) PLD/CMO médio (R$/MWh)",
                label_fmt=".0f", label_unit=" R$/MWh")
draw_brazil_map(axes[1, 1], Fcap_solar,
                cmap_name="RdYlGn_r", title="(d) Fator de captura solar (H1)",
                label_fmt=".2f", label_unit="")

# Sobrepor corredores de transmissão na sub-figura (d)
for ax in axes.flat:
    for (o, d) in ALL_LINKS:
        la_o, lo_o = SUBSYS_CENTROIDS[o]
        la_d, lo_d = SUBSYS_CENTROIDS[d]
        fmw = F_MAX[(o, d)]
        lw = max(1, fmw / 3000)
        ax.annotate("", xy=(lo_d, la_d), xytext=(lo_o, la_o),
                    arrowprops=dict(arrowstyle="->", color="#0a3d62",
                                    lw=lw, alpha=0.6))

fig.suptitle("Figura 21 — Mapas do SIN: capacidade solar, fator de capacidade, PLD e canibalização",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig21_mapas_sin.png", bbox_inches="tight")
plt.show()
''')


# =========================================================================
# SEÇÃO 20 — Análise detalhada da transmissão
# =========================================================================
md(r'''
## 20. Análise detalhada da transmissão entre subsistemas

A transmissão entre subsistemas é o elo entre o "ótimo local" e o
"equilíbrio sistêmico" do documento (Seção 11.4). Quando um corredor
satura (`|F_{ℓk,t}| = F^{max}_{ℓk}`), o preço-sombra local diverge —
é o **diferencial locacional** que justifica a formulação MFG.

A análise abaixo mede:

1. **Limites observados** `F^{max}_{ℓk}` (p99 do fluxo absoluto)
2. **Utilização horária** `|F_{ℓk,t}| / F^{max}_{ℓk}`
3. **Frequência de congestionamento** (% de horas com utilização > 80 %)
''')

code(r'''
# Reler intercâmbio interno (já existe intc_raw)
intc_use = intc_raw[["din_instante", "o", "d", "val_intercambiomwmed"]].copy()
intc_use = intc_use.rename(columns={"val_intercambiomwmed": "F"})
intc_use["F_abs"] = intc_use["F"].abs()
intc_use["par"]   = intc_use["o"] + "->" + intc_use["d"]
intc_use = add_calendar(intc_use)

# Utilização horária por corredor
F_MAX_S = pd.Series({f"{o}->{d}": F_MAX[(o, d)] for (o, d) in ALL_LINKS})
intc_use["F_max"] = intc_use["par"].map(F_MAX_S)
intc_use["util"]  = intc_use["F_abs"] / intc_use["F_max"]
intc_use["congest"] = (intc_use["util"] > 0.8).astype(int)

trans_summary = (intc_use.groupby("par")
                 .agg(F_max=("F_max", "first"),
                      F_abs_mean=("F_abs", "mean"),
                      util_mean=("util", "mean"),
                      util_p95=("util", lambda x: float(np.nanpercentile(x, 95))),
                      frac_congest=("congest", "mean"))
                 .reset_index())
trans_summary.to_csv(TABLES / "transmissao_utilizacao.csv", index=False, float_format="%.4f")
display(trans_summary.round(3))

# Figura: 3 painéis — F_max, utilização por hora do dia, frequência de congestão
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

ax = axes[0]
y_pos = np.arange(len(trans_summary))
ax.barh(y_pos, trans_summary["F_max"]/1000, color="#1f77b4")
ax.set_yticks(y_pos)
ax.set_yticklabels(trans_summary["par"])
ax.set_xlabel("F_max (GW)")
ax.set_title("(a) Capacidade máxima de fluxo por corredor"); ax.grid(alpha=.3, axis="x")
for i, v in enumerate(trans_summary["F_max"]/1000):
    ax.text(v + 0.1, i, f"{v:.1f}", va="center")

ax = axes[1]
util_hora = (intc_use.groupby(["par", "hora"])["util"].mean().unstack().T)
for par in util_hora.columns:
    ax.plot(util_hora.index, util_hora[par]*100, "o-", lw=1.5, label=par)
ax.axhline(80, color="red", ls="--", lw=1, label="80% (congestão)")
ax.set_xlabel("hora"); ax.set_ylabel("Utilização média (%)")
ax.set_title("(b) Utilização média por hora")
ax.legend(loc="best", fontsize=8); ax.grid(alpha=.3)

ax = axes[2]
ax.bar(trans_summary["par"], trans_summary["frac_congest"]*100, color="#d62728")
ax.set_ylabel("% horas com utilização > 80%")
ax.set_title("(c) Frequência de congestionamento")
ax.grid(alpha=.3, axis="y")
for i, v in enumerate(trans_summary["frac_congest"]*100):
    ax.text(i, v + 1, f"{v:.1f}%", ha="center", fontweight="bold")

fig.suptitle("Figura 22 — Transmissão entre subsistemas: capacidade, utilização e congestão (2025)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig22_transmissao_detalhada.png", bbox_inches="tight")
plt.show()
''')


# =========================================================================
# SEÇÃO 21 — Custo de entrada detalhado por subsistema
# =========================================================================
md(r'''
## 21. Custo de entrada para o agente solar — detalhamento regional

O documento (Seção 9.3) usa `I_i(K) = F_i + q_i K + η_i/2 K²`. Para
tornar o modelo defensável **com dados brasileiros**, separamos `q_ℓ`
em **três componentes**, todos por MW instalado e anualizados:

| Componente | Valor base | Fonte/justificativa |
| --- | --- | --- |
| `CAPEX_solar` | 3.5 a 4.5 milhões R$/MW | EPE 2024, ABSOLAR Boletim 2025 |
| `Custo de conexão` | 0.15 a 0.40 milhões R$/MW | TUST + EUST por subsistema |
| `O&M anual` | 0.06 a 0.08 milhões R$/MW/ano | ~2% do CAPEX |
| `WACC × vida` | (8%, 25 anos) → fator 0.094 | TLP-CDI + prêmio renovável |
| **`q_ℓ` anualizado total** | **350 a 500 mil R$/MW/ano** | soma anualizada |

Diferenças regionais:
- **NE**: CAPEX baixo (terreno barato), conexão alta (rede saturada), mas
  irradiância ótima → q_efetivo competitivo.
- **SE**: CAPEX médio, conexão barata (proximidade aos centros de carga),
  irradiância menor.
- **N**: CAPEX alto (logística amazônica), conexão difícil (rede esparsa),
  nebulosidade alta — solar menos competitiva.
- **S**: CAPEX médio, conexão moderada, irradiância baixa.

Para `q_ℓ` em **base 8 dias típicos** (período do modelo), aplicamos
`q_ℓ * (192/8760)`.
''')

code(r'''
# Premissas regionais (milhões R$/MW)
capex_reg = {"N": 4.2, "NE": 3.6, "SE": 3.8, "S": 4.0}     # CAPEX bruto
conn_reg  = {"N": 0.40, "NE": 0.30, "SE": 0.15, "S": 0.22} # conexão
om_reg    = {"N": 0.084, "NE": 0.072, "SE": 0.076, "S": 0.080}  # O&M anual

WACC = 0.08
VIDA = 25
FATOR_ANUAL = WACC * (1 + WACC)**VIDA / ((1 + WACC)**VIDA - 1)  # ≈ 0.0937

q_reg = {}
for s in SUBSYS_ORDER:
    capex_total = (capex_reg[s] + conn_reg[s]) * 1e6   # R$/MW
    q_anual = capex_total * FATOR_ANUAL + om_reg[s] * 1e6
    q_reg[s] = q_anual

# Versão por horizonte (192 h)
q_reg_192h = {s: q_reg[s] * (HOURS / HOURS_YEAR) for s in SUBSYS_ORDER}

inv_table = pd.DataFrame({
    "CAPEX_MR$/MW": capex_reg,
    "Conexao_MR$/MW": conn_reg,
    "O&M_MR$/MW/ano": om_reg,
    "Total_capex_MR$/MW": {s: capex_reg[s] + conn_reg[s] for s in SUBSYS_ORDER},
    "q_anualizado_R$/MW/ano": q_reg,
    "q_horizonte_R$/MW (192h)": q_reg_192h,
}).T
inv_table.columns = SUBSYS_ORDER
inv_table.to_csv(TABLES / "custo_entrada_solar_regional.csv", float_format="%.0f")
print("Custo de entrada para agente solar — premissas regionais:")
display(inv_table.round(2))

# Figura: composição do custo + comparação receita média
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
df_stack = pd.DataFrame({
    "CAPEX": [capex_reg[s] for s in SUBSYS_ORDER],
    "Conexão": [conn_reg[s] for s in SUBSYS_ORDER],
    "O&M (25y)": [om_reg[s] * VIDA for s in SUBSYS_ORDER],  # somado em 25y
}, index=SUBSYS_ORDER)
df_stack.plot.bar(stacked=True, ax=ax, color=["#1f77b4", "#ff7f0e", "#2ca02c"])
ax.set_ylabel("milhões R$/MW (vida útil)")
ax.set_title("(a) Composição do custo de entrada por subsistema")
ax.legend(); ax.grid(alpha=.3, axis="y")

ax = axes[1]
# Receita média horária = PLD_avg * CF_solar
rcap_avg = [PLD_avg[s] * CF_solar[s] * 8760 / 1e6 for s in SUBSYS_ORDER]
qreg_y = [q_reg[s] / 1e6 for s in SUBSYS_ORDER]
x = np.arange(len(SUBSYS_ORDER))
ax.bar(x - 0.2, rcap_avg, width=0.4, color="#2ca02c", label="Receita esperada (R$/MW/ano)")
ax.bar(x + 0.2, qreg_y,   width=0.4, color="#d62728", label="Custo q_ℓ (R$/MW/ano)")
ax.set_xticks(x); ax.set_xticklabels(SUBSYS_ORDER)
ax.set_ylabel("milhões R$/MW/ano")
ax.set_title("(b) Receita esperada × custo q_ℓ (sem canibalização)")
ax.legend(); ax.grid(alpha=.3, axis="y")
for i, (r, q) in enumerate(zip(rcap_avg, qreg_y)):
    margin = r - q
    ax.text(i, max(r, q) + 0.02, f"Δ={margin:+.2f}", ha="center", fontsize=9,
            color="green" if margin > 0 else "red")

fig.suptitle("Figura 23 — Custo de entrada regional e margem ex‑ante",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig23_custo_entrada_regional.png", bbox_inches="tight")
plt.show()
''')


# =========================================================================
# SEÇÃO 22 — PLD como receita do agente
# =========================================================================
md(r'''
## 22. PLD horário — proxy a partir do CMO com tetos regulatórios

A receita do **agente solar** é determinada pelo **PLD horário** (Preço
de Liquidação das Diferenças), não pelo CMO interno do modelo. Em 2025
o PLD horário CCEE opera como `PLD = clip(CMO, PLD^{min}, PLD^{max})`
com:

- `PLD^{max}_{2025} ≈ 939.96 R$/MWh` (Resolução Homologatória ANEEL nº 3.342/2024)
- `PLD^{min}_{2025} ≈ 65.69 R$/MWh`

Construímos `PLD̂_t = clip(CMO_t, PLD^{min}, PLD^{max})` e calculamos a
receita capturada `R^{cap,PLD}_{ℓ,r}` substituindo `λ_t` por `PLD̂_t`
nas fórmulas da §4.2 do documento.
''')

code(r'''
PLD_MAX_2025 = 939.96
PLD_MIN_2025 = 65.69

panel["pld_hat"] = panel["cmo_h"].clip(lower=PLD_MIN_2025, upper=PLD_MAX_2025)

# Receita capturada PLD-ponderada para solar e eólica
pld_rev = []
for sub in SUBSYS_ORDER:
    g = panel[panel["id_subsistema"] == sub].dropna(subset=["pld_hat"])
    for col, src in [("gs", "solar"), ("gr", "wind")]:
        Es = g[col].sum()
        if Es <= 0:
            continue
        R_PLD = (g["pld_hat"] * g[col]).sum() / Es      # R$/MWh capturados
        R_CMO = (g["cmo_h"] * g[col]).sum() / Es
        pld_rev.append({"id_subsistema": sub, "fonte": src,
                        "PLD_avg": float(g["pld_hat"].mean()),
                        "P_cap_PLD": R_PLD,
                        "P_cap_CMO": R_CMO,
                        "F_cap_PLD": R_PLD / float(g["pld_hat"].mean())})
pld_rev = pd.DataFrame(pld_rev)
pld_rev.to_csv(TABLES / "receita_capturada_PLD.csv", index=False, float_format="%.2f")
display(pld_rev.round(2))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
# CMO vs PLD ao longo do tempo
agg_price = panel.groupby("din_instante").agg(cmo=("cmo_h", "mean"), pld=("pld_hat", "mean"))
ax.plot(agg_price.index, agg_price["cmo"], "-", color="#1f77b4", lw=1, label="CMO (média SIN)", alpha=.6)
ax.plot(agg_price.index, agg_price["pld"], "-", color="#d62728", lw=1, label="PLD̂ = clip(CMO)", alpha=.8)
ax.axhline(PLD_MIN_2025, color="green", ls="--", lw=0.8, label=f"PLD_min={PLD_MIN_2025:.0f}")
ax.axhline(PLD_MAX_2025, color="red",   ls="--", lw=0.8, label=f"PLD_max={PLD_MAX_2025:.0f}")
ax.set_ylabel("R$/MWh"); ax.set_title("(a) CMO vs PLD̂ (com tetos regulatórios)")
ax.legend(loc="upper left", fontsize=8); ax.grid(alpha=.3)

ax = axes[1]
piv = pld_rev.pivot(index="id_subsistema", columns="fonte", values="F_cap_PLD").reindex(SUBSYS_ORDER)
piv.plot.bar(ax=ax, color=["#fdae61", "#74add1"])
ax.axhline(1.0, color="black", ls="--", lw=1)
ax.set_ylabel("F_capture PLD")
ax.set_title("(b) Fator de captura calculado com PLD̂ (versus PLD médio)")
ax.legend(title="fonte"); ax.grid(alpha=.3, axis="y")

fig.suptitle("Figura 24 — PLD̂ horário (proxy via CMO + tetos) e receita capturada",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig24_pld_receita.png", bbox_inches="tight")
plt.show()
''')


# =========================================================================
# SEÇÃO 23 — MFG com decisão de localização do agente solar
# =========================================================================
md(r'''
## 23. MFG com **decisão de localização** do agente solar

Aqui o agente solar **escolhe onde investir**. A variável de decisão
é a alocação `K_ℓ` (capacidade solar no subsistema ℓ), com restrição
agregada `Σ_ℓ K_ℓ ≤ K^{total}_{disp}` (orçamento de investimento).

O equilíbrio de **best‑response** é caracterizado por:

$$
\frac{\partial \Pi_{\ell}}{\partial K_\ell}
= \underbrace{\sum_t \mathrm{PLD}_{\ell,t}(K)\, a^{s}_{\ell,t}}_{\text{receita marginal}}
- q_\ell = \mu \quad \forall \ell \text{ com } K_\ell > 0
$$

isto é, a **margem por MW** é igualada entre subsistemas no ótimo
privado (com `μ ≥ 0` multiplicador do orçamento). O preço local
`PLD_{ℓ,t}(K)` é **endógeno**: depende de quanto solar entra em cada
local via o clearing locacional com transmissão.

### Algoritmo: iteração best‑response + clearing locacional

1. Inicializar `K_ℓ⁰` ← distribuição observada.
2. Clearing locacional com transmissão → preços `PLD_{ℓ,t}(K^k)`.
3. Para cada subsistema, calcular margem `m_ℓ = Σ_t PLD_{ℓ,t} a^s_{ℓ,t} − q_ℓ`.
4. Realocar capacidade para subsistemas com maior margem (proporcional).
5. Repetir até convergência.
''')

code(r'''
# Total de capacidade solar a alocar
K_TOT = float(K_inst.loc[:, "K_s"].sum())   # ~37 GW (proxy observado)

# O clearing locacional é em base 24h. Convertemos q_anual → q por dia típico
# anualizando a receita: rev_dia * 365 = rev_anual a ser comparada com q_anual.
# Equivalente: usar q_dia = q_anual / 365.
q_loc_dia = {s: q_reg[s] / 365.0 for s in SUBSYS_ORDER}
print(f"Orçamento total K_TOT = {K_TOT/1000:.1f} GW")
print("q_ℓ (R$/MW · dia):", {s: f"{v:.0f}" for s, v in q_loc_dia.items()})

# Perfis horários por subsistema (dia típico anual)
a_solar_l = {s: typ_loc.query("id_subsistema == @s").sort_values("hora")["gs"].values
             / max(K_inst.loc[s, "K_s"], 1e-6)
             for s in SUBSYS_ORDER}

# Demanda residual sem solar (para calcular preço endógeno)
D_lt_nosol = {s: (typ_loc.query("id_subsistema == @s").sort_values("hora")["D"].values
                  - typ_loc.query("id_subsistema == @s").sort_values("hora")["gr"].values)
              for s in SUBSYS_ORDER}

def clearing_with_solar_alloc(K_alloc):
    """Clearing locacional dado alocação de capacidade solar K_alloc[ℓ]."""
    Gren_with = {s: a_solar_l[s] * K_alloc[s]
                    + typ_loc.query("id_subsistema == @s").sort_values("hora")["gr"].values
                 for s in SUBSYS_ORDER}
    return clearing_locacional(D_lt, Gren_with, K_T_loc, K_H_loc,
                                c1_loc, c2_loc, Havail_lt)

def best_response_alloc(K_total, n_iter=15, eta=0.3, eps=1e-3):
    """Iteração best-response para alocação solar entre subsistemas.

    Receita por MW (R$/MW/ano) = 365 × Σ_t PLD_t · a^s_t (clearing 24h)
    Margem por MW = receita_anual − q_anual
    A realocação favorece subsistemas com maior margem **relativa** —
    se todas margens são negativas, alocamos proporcionalmente à receita
    (best-response continua válida mesmo com lucro negativo, porque o
    agente compara entre alternativas).
    """
    K_obs = {s: float(K_inst.loc[s, "K_s"]) for s in SUBSYS_ORDER}
    scale = K_total / sum(K_obs.values())
    K = {s: K_obs[s] * scale for s in SUBSYS_ORDER}

    hist = []
    for it in range(n_iter):
        res = clearing_with_solar_alloc(K)
        margens, receitas = {}, {}
        for s in SUBSYS_ORDER:
            pld = np.clip(res["lam"][s], PLD_MIN_2025, PLD_MAX_2025)
            rev_dia    = float((pld * a_solar_l[s]).sum())     # R$/MW · dia
            rev_anual  = rev_dia * 365.0                         # R$/MW · ano
            margens[s]  = rev_anual - q_reg[s]
            receitas[s] = rev_anual
        # Realocação proporcional à receita (best-response sob orçamento total)
        # com componente de margem (incentivo extra à região mais lucrativa)
        r_arr = np.array([receitas[s] for s in SUBSYS_ORDER])
        if r_arr.sum() <= 0:
            break
        m_arr = np.array([margens[s] for s in SUBSYS_ORDER])
        # softmax-like weight: combina nível de receita + bônus para margem positiva
        m_norm = (m_arr - m_arr.min()) / max(m_arr.max() - m_arr.min(), 1.0)
        weights = r_arr * (1 + m_norm)
        w_new = weights / weights.sum()
        K_new = {s: (1 - eta) * K[s] + eta * w_new[i] * K_total
                 for i, s in enumerate(SUBSYS_ORDER)}
        err = max(abs(K_new[s] - K[s]) for s in SUBSYS_ORDER)
        hist.append({"it": it, "err": err, "margens": margens, "receitas": receitas,
                     "K": dict(K_new)})
        K = K_new
        if err < eps * K_total:
            break
    return K, hist

K_eq, hist_br = best_response_alloc(K_TOT, n_iter=20)
print("\nDistribuição de capacidade solar no equilíbrio MFG:")
last = hist_br[-1] if hist_br else {"margens": {s: np.nan for s in SUBSYS_ORDER},
                                    "receitas": {s: np.nan for s in SUBSYS_ORDER}}
df_eq = pd.DataFrame({
    "K_obs (GW)":          {s: K_inst.loc[s, "K_s"]/1000 for s in SUBSYS_ORDER},
    "K_eq (GW)":           {s: K_eq[s]/1000 for s in SUBSYS_ORDER},
    "Δ (GW)":              {s: (K_eq[s] - K_inst.loc[s, "K_s"])/1000 for s in SUBSYS_ORDER},
    "margem R$/MW/ano":    last["margens"],
    "receita R$/MW/ano":   last["receitas"],
    "q_ℓ R$/MW/ano":       q_reg,
}).T
df_eq.columns = SUBSYS_ORDER
df_eq.to_csv(TABLES / "mfg_alocacao_solar.csv", float_format="%.2f")
display(df_eq.round(2))

# Figura: convergência + comparação observada vs equilíbrio
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

ax = axes[0]
errs = [h["err"]/K_TOT for h in hist_br]
ax.plot(range(len(errs)), errs, "o-", color="#1f77b4")
ax.set_yscale("log"); ax.set_xlabel("iteração"); ax.set_ylabel("‖ΔK‖∞ / K_TOT")
ax.set_title("(a) Convergência do best-response")
ax.grid(alpha=.3)

ax = axes[1]
sub_x = np.arange(len(SUBSYS_ORDER))
K_obs_arr = np.array([K_inst.loc[s, "K_s"]/1000 for s in SUBSYS_ORDER])
K_eq_arr  = np.array([K_eq[s]/1000 for s in SUBSYS_ORDER])
ax.bar(sub_x - 0.2, K_obs_arr, width=0.4, label="K_observada", color="#1f77b4")
ax.bar(sub_x + 0.2, K_eq_arr,  width=0.4, label="K_equilíbrio MFG", color="#d62728")
ax.set_xticks(sub_x); ax.set_xticklabels(SUBSYS_ORDER)
ax.set_ylabel("GW"); ax.set_title("(b) Capacidade solar por subsistema")
ax.legend(); ax.grid(alpha=.3, axis="y")
for i, (o, e) in enumerate(zip(K_obs_arr, K_eq_arr)):
    delta = (e - o)
    ax.text(i, max(o, e) + 0.5, f"Δ={delta:+.1f}", ha="center", fontsize=9,
            color="green" if delta > 0 else "red")

ax = axes[2]
m_arr = np.array([last["margens"][s]/1e3 for s in SUBSYS_ORDER])
ax.bar(SUBSYS_ORDER, m_arr, color=["green" if m > 0 else "red" for m in m_arr])
ax.axhline(0, color="black", lw=1)
ax.set_ylabel("Margem (R$ ×10³/MW/ano)")
ax.set_title("(c) Margem (receita − q_ℓ) por subsistema no equilíbrio")
ax.grid(alpha=.3, axis="y")

fig.suptitle("Figura 25 — MFG com decisão de localização: best-response e alocação ótima privada",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig25_mfg_alocacao.png", bbox_inches="tight")
plt.show()
''')


md(r'''
### Mapa do equilíbrio de localização

Mostra **graficamente** onde o agente privado prefere investir versus
onde a capacidade já está instalada hoje.
''')

code(r'''
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

draw_brazil_map(axes[0], {s: K_inst.loc[s, "K_s"]/1000 for s in SUBSYS_ORDER},
                cmap_name="YlOrRd", title="(a) Distribuição observada (2025)",
                label_fmt=".1f", label_unit=" GW")
draw_brazil_map(axes[1], {s: K_eq[s]/1000 for s in SUBSYS_ORDER},
                cmap_name="YlOrRd", title="(b) Distribuição de equilíbrio MFG (best-response)",
                label_fmt=".1f", label_unit=" GW",
                vmin=0, vmax=max(K_inst["K_s"])/1000)

# Sobrepor corredores de transmissão
for ax in axes:
    for (o, d) in ALL_LINKS:
        la_o, lo_o = SUBSYS_CENTROIDS[o]
        la_d, lo_d = SUBSYS_CENTROIDS[d]
        ax.annotate("", xy=(lo_d, la_d), xytext=(lo_o, la_o),
                    arrowprops=dict(arrowstyle="->", color="#0a3d62",
                                    lw=max(1, F_MAX[(o, d)]/3000), alpha=0.5))

fig.suptitle("Figura 26 — Mapa: distribuição observada vs equilíbrio MFG da capacidade solar",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGS / "fig26_mapa_equilibrio_mfg.png", bbox_inches="tight")
plt.show()
''')


md(r'''
## 24. Erros / inconsistências encontradas no documento

### (E1) Resumo — typo "avaliarssa"
> *"… define hipóteses empiricamente testáveis para **avaliarssa** a ser
> limitada por flexibilidade …"*

**Sugestão:** *"… para **avaliar essa nova fronteira**, que tende a ser
limitada por …"*

### (E2) Inconsistência D_t / L_t / D_global
A Seção 4.1 declara que `L_t` é a carga observada e `D_t` é demanda
econômica. Porém:
- A Seção 1 usa `D_t` para carga: `g + ... = D_t + c_t - u_t`.
- A Seção 3.1 introduz `D^global, D^RB, D^res` (carga observada).
- A Seção 4.1 distingue `L_t ≠ D_t`, mas todas as fórmulas posteriores
  seguem usando `D_t` para carga.

**Sugestão:** Eliminar `L_t`. Manter `D_t` para carga observada e usar
`D̃(p)` apenas se introduzir função de demanda econômica.

### (E3) Balanço situacional sem eólica (Seção 9.1)
> `g^s_{1,t} + g^s_{2,t} + h_t + n_t = D_t + c_t - u_t`

A equação omite eólica. Para comparação empírica com o Brasil, tratar
eólica como exógena `G^w_t`:
> `g^s_{1,t} + g^s_{2,t} + G^w_t + h_t + n_t = D_t + c_t - u_t`

Implementado dessa forma neste notebook.

### (E4) Banda de preço na Seção 11.4
> `-π_c ≤ λ_{ℓ,t} ≤ π_u`

O sinal negativo é correto na convenção `g = D + c - u` mas precisa ser
explicitado em §11.4 / Apêndice A.

### (E5) "Outras saídas" no balanço hídrico
A Seção 4.4 introduz `q^out_t` sem defini-lo. Os arquivos do ONS expõem
`val_vazaoevaporacaoliquida` e `val_vazaousoconsuntivo`, que podem
preencher: `q^out_t = q^evap_t + q^cons_t`.

### (E6) Unidades de q_i
A unidade de `q_i` em `I(K) = F + qK + ηK²/2` não é definida. Deve
estar em R$/(MW·ano) para que `I` esteja em R$/ano.

### (E7) Convenção do dual no clearing
Em CVXPY o sinal do `dual_value` depende do solver. Usamos
`lam = -bal.dual_value` para que `λ ≥ 0` em horas de escassez,
consistente com a banda do documento.

### (E8) PLD horário em 2025
O PLD horário começou em 2024. Para 2025 já há série oficial. Usamos
CMO semi-horário agregado (mais granular).

### (E9) Inconsistência entre nomenclatura `Δ` (Seção 3.1) e `R` (Seção 4.3)
A Seção 3.1 usa `ΔD^res_t` para rampa. A Seção 4.3 usa `R^net_t`. São o
mesmo objeto.

### (E10) Disponibilidade de dados
Os dados de carga do ONS de 2025 cessam em **2025-10-24**, não
2025-10-31. As demais variáveis vão até 2025-11-05 (geração) ou
2025-12-31 (CMO/hidrologia). Filtramos tudo no janela conjunta.
''')


md(r'''
## 25. Resumo executivo e inventário de saídas
''')

code(r'''
print("== Figuras geradas ==")
for p in sorted(FIGS.glob("*.png")):
    print(" -", p.relative_to(ROOT))
print()
print("== Tabelas geradas ==")
for p in sorted(TABLES.glob("*.csv")):
    print(" -", p.relative_to(ROOT))
''')


# ==========================================================================
nb["cells"] = cells
NB_PATH.write_text(nbf.writes(nb), encoding="utf-8")
print(f"Notebook escrito em {NB_PATH} ({len(cells)} celulas)")
