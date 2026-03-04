"""
Pipeline de dados do SIN (Sistema Interligado Nacional).

Módulo compartilhado que centraliza:
  - Leitura robusta de CSVs (detecção de separador, encoding, formato numérico BR)
  - Normalização de nomes de colunas e subsistemas
  - Carregamento de demanda, geração por fonte e intercâmbio
  - Construção do painel unificado (hora × subsistema)
  - Funções auxiliares de visualização (heatmaps, curvas de duração)
"""

from __future__ import annotations

import os
import re
import glob
import unicodedata
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

SUBSYS_MAP: dict[str, str] = {
    "norte": "N", "nordeste": "NE", "sudeste": "SE", "sul": "S",
    "n": "N", "ne": "NE", "se": "SE", "s": "S",
}

SUBSYS_ORDER: list[str] = ["N", "NE", "SE", "S"]

MESES_PT: dict[int, str] = {
    1: "Jan", 2: "Fev", 3: "Mar", 4: "Abr", 5: "Mai", 6: "Jun",
    7: "Jul", 8: "Ago", 9: "Set", 10: "Out", 11: "Nov", 12: "Dez",
}

# ---------------------------------------------------------------------------
# Normalização de texto
# ---------------------------------------------------------------------------

def norm(s: str) -> str:
    """Remove acentos, converte para minúsculas e substitui caracteres
    não-alfanuméricos por underscore."""
    s = str(s).strip()
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    s = s.lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def canonical_subsys(x: str | None) -> str | None:
    """Converte nome livre de subsistema para código canônico (N/NE/SE/S)."""
    if x is None:
        return None
    k = norm(x)
    if k in SUBSYS_MAP:
        return SUBSYS_MAP[k]
    for kk in sorted(SUBSYS_MAP, key=len, reverse=True):
        if k == kk or k.endswith("_" + kk):
            return SUBSYS_MAP[kk]
    return None


def subsys_from_filename(path: str) -> str | None:
    """Extrai subsistema a partir do nome do arquivo."""
    fn = norm(os.path.basename(path))
    for token in ["nordeste", "sudeste", "norte", "sul", "ne", "se", "n", "s"]:
        if token in fn:
            return canonical_subsys(token)
    return None


# ---------------------------------------------------------------------------
# Leitura robusta de CSV
# ---------------------------------------------------------------------------

def _sniff_sep(sample: str) -> str | None:
    candidates = [";", ",", "\t", "|"]
    counts = {c: sample.count(c) for c in candidates}
    best = max(counts, key=counts.get)
    return best if counts[best] > 0 else None


def read_csv_robust(path: str) -> pd.DataFrame:
    """Lê CSV com detecção automática de separador, fallback de encoding,
    normalização de colunas e correção de artefatos de encoding em strings."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")

    with open(path, "rb") as f:
        raw = f.read(4096)
    try:
        sample = raw.decode("utf-8", errors="ignore")
    except Exception:
        sample = raw.decode("latin1", errors="ignore")
    sep = _sniff_sep(sample)

    last_err = None
    df = None
    for enc in ["utf-8", "latin1"]:
        try:
            kw = dict(engine="python", encoding=enc)
            df = pd.read_csv(path, sep=sep, **kw) if sep else pd.read_csv(path, sep=None, **kw)
            break
        except Exception as e:
            last_err = e
            df = None

    if df is None:
        raise RuntimeError(f"Falha ao ler {path}. Último erro: {last_err}")

    for c in df.select_dtypes(include=["object"]).columns:
        s = df[c].astype(str)
        try:
            fixed = s.str.encode("latin1", errors="ignore").str.decode("utf-8", errors="ignore")
            if fixed.str.contains("Ã").mean() < s.str.contains("Ã").mean():
                df[c] = fixed
        except Exception:
            pass

    df.columns = [norm(c) for c in df.columns]
    return df


def to_float(x: pd.Series) -> pd.Series:
    """Converte série de strings numéricas para float, tratando formato BR
    (vírgula decimal, ponto de milhar) de forma robusta."""
    s = x.astype(str).str.strip().str.replace(r"\s+", "", regex=True)

    has_comma = s.str.contains(",", regex=False)
    s.loc[has_comma] = (
        s.loc[has_comma]
        .str.replace(".", "", regex=False)
        .str.replace(",", ".", regex=False)
    )

    no_comma = ~has_comma
    multi_dot = no_comma & (s.str.count(r"\.") >= 2)

    def _keep_last_dot(val: str) -> str:
        parts = val.split(".")
        return "".join(parts[:-1]) + "." + parts[-1]

    s.loc[multi_dot] = s.loc[multi_dot].apply(_keep_last_dot)
    return pd.to_numeric(s, errors="coerce")


# ---------------------------------------------------------------------------
# Detecção automática de colunas
# ---------------------------------------------------------------------------

def detect_datetime_col(df: pd.DataFrame) -> str | None:
    cands = [c for c in df.columns if any(k in c for k in
             ["din_instante", "data", "hora", "datetime", "timestamp"])]
    return cands[0] if cands else None


def detect_load_col(df: pd.DataFrame) -> str | None:
    if "val_cargaenergiahomwmed" in df.columns:
        return "val_cargaenergiahomwmed"
    cands = [c for c in df.columns if any(k in c for k in
             ["carga", "demanda", "mwmed", "load", "mw"])]
    return cands[0] if cands else None


def detect_gen_col(df: pd.DataFrame) -> str | None:
    if "val_geracao" in df.columns:
        return "val_geracao"
    cands = [c for c in df.columns if any(k in c for k in ["geracao", "generation"])]
    return cands[0] if cands else None


def detect_subsys_col(df: pd.DataFrame) -> str | None:
    for c in ["id_subsistema", "cod_subsistema", "nom_subsistema", "subsistema"]:
        if c in df.columns:
            return c
    return None


# ---------------------------------------------------------------------------
# Loaders de dados
# ---------------------------------------------------------------------------

def load_demanda_efetiva(paths: list[str]) -> pd.DataFrame:
    """Carrega curvas de carga reais (demanda efetiva) e agrega por hora."""
    out = []
    for p in paths:
        df = read_csv_robust(p)
        dtc = detect_datetime_col(df)
        lc = detect_load_col(df)
        if dtc is None or lc is None:
            raise ValueError(f"[{p}] colunas datetime/carga não encontradas. cols={df.columns.tolist()}")
        ss = subsys_from_filename(p)
        d = pd.DataFrame({
            "din_instante": pd.to_datetime(df[dtc], errors="coerce").dt.floor("h"),
            "D": to_float(df[lc]),
            "id_subsistema": ss,
        })
        d = d.dropna(subset=["din_instante", "id_subsistema", "D"])
        d = d.groupby(["din_instante", "id_subsistema"], as_index=False)["D"].mean()
        out.append(d)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["din_instante", "id_subsistema", "D"])


def load_demanda_prevista(paths: list[str]) -> pd.DataFrame:
    """Carrega previsões dia-seguinte de demanda por patamar."""
    out = []
    for p in paths:
        df = read_csv_robust(p)
        day_col = next((c for c in df.columns if any(k in c for k in
                        ["din_programacaodia", "programacao", "dia", "data"])), None)
        pat_col = next((c for c in df.columns if "patamar" in c), None)
        val_col = ("val_demanda" if "val_demanda" in df.columns else
                   next((c for c in df.columns if "demanda" in c or "carga" in c), None))
        if day_col is None or pat_col is None or val_col is None:
            raise ValueError(f"[{p}] colunas dia/patamar/demanda não encontradas. cols={df.columns.tolist()}")
        ss = subsys_from_filename(p)
        dd = pd.DataFrame({
            "dia": pd.to_datetime(df[day_col], errors="coerce").dt.date,
            "patamar": pd.to_numeric(df[pat_col], errors="coerce"),
            "D_prev_pat": to_float(df[val_col]),
            "id_subsistema": ss,
        })
        dd = dd.dropna(subset=["dia", "patamar", "D_prev_pat", "id_subsistema"])
        dd = dd.groupby(["id_subsistema", "dia", "patamar"], as_index=False)["D_prev_pat"].mean()
        out.append(dd)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["id_subsistema", "dia", "patamar", "D_prev_pat"])


def load_generation(path: str, fonte: str) -> pd.DataFrame:
    """Carrega geração de uma fonte (solar, wind, hydro, nuclear, thermal)
    e agrega por hora e subsistema."""
    df = read_csv_robust(path)
    dtc = detect_datetime_col(df)
    vc = detect_gen_col(df)
    sc = detect_subsys_col(df)
    if dtc is None or vc is None or sc is None:
        raise ValueError(f"[{path}] colunas (datetime, geracao, subsistema) não encontradas. cols={df.columns.tolist()}")
    g = pd.DataFrame({
        "din_instante": pd.to_datetime(df[dtc], errors="coerce").dt.floor("h"),
        "val_geracao": to_float(df[vc]),
        "id_subsistema": df[sc].astype(str).str.strip().apply(canonical_subsys),
        "fonte": fonte,
    })
    g = g.dropna(subset=["din_instante", "id_subsistema", "val_geracao"])
    g = g.groupby(["din_instante", "id_subsistema", "fonte"], as_index=False)["val_geracao"].sum()
    return g


def load_intercambio_interno(path: str) -> pd.DataFrame:
    """Carrega fluxos de potência entre subsistemas e retorna importação
    líquida x_int(t, s) em MWmed por hora."""
    if not os.path.exists(path):
        return pd.DataFrame(columns=["din_instante", "id_subsistema", "x_int"])

    df = read_csv_robust(path)
    dtc = detect_datetime_col(df)
    valc = "val_intercambiomwmed" if "val_intercambiomwmed" in df.columns else None
    oc = "id_subsistema_origem" if "id_subsistema_origem" in df.columns else None
    dc = "id_subsistema_destino" if "id_subsistema_destino" in df.columns else None
    if dtc is None or valc is None or oc is None or dc is None:
        raise ValueError(f"[{path}] colunas esperadas não encontradas. cols={df.columns.tolist()}")

    t = pd.to_datetime(df[dtc], errors="coerce").dt.floor("h")
    f = to_float(df[valc])
    o = df[oc].astype(str).str.strip().apply(canonical_subsys)
    d = df[dc].astype(str).str.strip().apply(canonical_subsys)

    base = pd.DataFrame({"din_instante": t, "o": o, "d": d, "f": f}).dropna()

    out_o = base[["din_instante", "o", "f"]].rename(columns={"o": "id_subsistema"}).copy()
    out_o["x_int"] = -out_o["f"]
    out_o = out_o.drop(columns=["f"])

    out_d = base[["din_instante", "d", "f"]].rename(columns={"d": "id_subsistema"}).copy()
    out_d["x_int"] = +out_d["f"]
    out_d = out_d.drop(columns=["f"])

    result = pd.concat([out_o, out_d], ignore_index=True)
    result = result.groupby(["din_instante", "id_subsistema"], as_index=False)["x_int"].sum()
    return result


def load_intercambio_sin(path: str) -> pd.DataFrame:
    """Carrega intercâmbio externo agregado do SIN (por país)."""
    if not os.path.exists(path):
        return pd.DataFrame(columns=["din_instante", "nom_paisdestino", "val_intercambiomwmed"])

    df = read_csv_robust(path)
    dtc = detect_datetime_col(df)
    valc = "val_intercambiomwmed" if "val_intercambiomwmed" in df.columns else None
    pc = "nom_paisdestino" if "nom_paisdestino" in df.columns else None
    if dtc is None or valc is None or pc is None:
        raise ValueError(f"[{path}] colunas esperadas não encontradas. cols={df.columns.tolist()}")

    result = pd.DataFrame({
        "din_instante": pd.to_datetime(df[dtc], errors="coerce").dt.floor("h"),
        "nom_paisdestino": df[pc].astype(str).str.strip(),
        "val_intercambiomwmed": to_float(df[valc]),
    })
    return result.dropna()


# ---------------------------------------------------------------------------
# Loaders de preços (CMO)
# ---------------------------------------------------------------------------

def load_cmo_semanal(path: str) -> pd.DataFrame:
    """Carrega CMO semanal (resultado DECOMP) com patamares de carga.

    Retorna DataFrame com colunas:
        din_instante, id_subsistema, cmo_medio, cmo_leve, cmo_media, cmo_pesada
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")
    df = read_csv_robust(path)
    dtc = detect_datetime_col(df)
    sc = detect_subsys_col(df)
    if dtc is None or sc is None:
        raise ValueError(f"[{path}] colunas datetime/subsistema não encontradas. cols={df.columns.tolist()}")

    result = pd.DataFrame({
        "din_instante": pd.to_datetime(df[dtc], errors="coerce"),
        "id_subsistema": df[sc].astype(str).str.strip().apply(canonical_subsys),
    })
    col_map = {
        "val_cmomediasemanal": "cmo_medio",
        "val_cmoleve": "cmo_leve",
        "val_cmomedia": "cmo_media",
        "val_cmopesada": "cmo_pesada",
    }
    for orig, dest in col_map.items():
        if orig in df.columns:
            result[dest] = to_float(df[orig])
        else:
            result[dest] = np.nan

    return result.dropna(subset=["din_instante", "id_subsistema"])


def load_cmo_semihorario(path: str) -> pd.DataFrame:
    """Carrega CMO semi-horário (30 em 30 min) com valor único de CMO.

    Retorna DataFrame com colunas:
        din_instante, id_subsistema, cmo
    Opcionalmente agrega para resolução horária com cmo_h (média de 2 meias-horas).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")
    df = read_csv_robust(path)
    dtc = detect_datetime_col(df)
    sc = detect_subsys_col(df)
    val_col = "val_cmo" if "val_cmo" in df.columns else None
    if dtc is None or sc is None or val_col is None:
        raise ValueError(f"[{path}] colunas esperadas não encontradas. cols={df.columns.tolist()}")

    result = pd.DataFrame({
        "din_instante": pd.to_datetime(df[dtc], errors="coerce"),
        "id_subsistema": df[sc].astype(str).str.strip().apply(canonical_subsys),
        "cmo": to_float(df[val_col]),
    })
    return result.dropna(subset=["din_instante", "id_subsistema", "cmo"])


def load_cmo_horario(path: str) -> pd.DataFrame:
    """Carrega CMO semi-horário e agrega para resolução horária (média)."""
    df = load_cmo_semihorario(path)
    df["din_instante_h"] = df["din_instante"].dt.floor("h")
    hourly = (df.groupby(["din_instante_h", "id_subsistema"], as_index=False)["cmo"]
              .mean()
              .rename(columns={"din_instante_h": "din_instante", "cmo": "cmo_h"}))
    return hourly


# ---------------------------------------------------------------------------
# Construção do painel unificado
# ---------------------------------------------------------------------------

class SINPaths:
    """Encapsula todos os caminhos de dados do SIN para um dado ano."""

    def __init__(self, root: str = "validate_model", year: int = 2025):
        self.root = root
        self.data = os.path.join(root, "data")
        self.outputs = os.path.join(root, "outputs")
        os.makedirs(self.outputs, exist_ok=True)

        self.curva_paths = sorted(glob.glob(
            os.path.join(self.data, "demanda_efetiva", f"CURVA_CARGA_*_{year}.csv")))
        self.prev_paths = sorted(glob.glob(
            os.path.join(self.data, "demanda_esperada", f"DEMANDA_*_{year}.csv")))
        self.solar_path = os.path.join(self.data, "producao_solar", f"fotovoltaicas_{year}.csv")
        self.wind_path = os.path.join(self.data, "producao_eolica", f"eolicas_{year}.csv")
        self.nuclear_path = os.path.join(self.data, "producao_non_renewable", f"nuclear_{year}.csv")
        self.term_path = os.path.join(self.data, "producao_non_renewable", f"TERMICAS_{year}.csv")
        self.hydro_paths = sorted(glob.glob(
            os.path.join(self.data, "producao_*", f"*hidro*_{year}.csv")))
        self.interc_sin_path = os.path.join(self.data, "intercambio", f"Intercambio_do_SIN_{year}.csv")
        self.interc_interno_path = os.path.join(self.data, "intercambio", f"intercambio_interno_{year}.csv")

        self.cmo_semanal_path = os.path.join(self.data, "precos", f"cmo_semanal{year}.csv")
        self.cmo_semihorario_path = os.path.join(self.data, "precos", f"cmo_semihorario{year}.csv")

    def summary(self) -> None:
        """Imprime resumo dos caminhos e existência dos arquivos."""
        print(f"ROOT: {self.root}")
        print(f"Curva carga ({len(self.curva_paths)} arquivos): {self.curva_paths}")
        print(f"Previsões   ({len(self.prev_paths)} arquivos): {self.prev_paths}")
        for name, p in [("Solar", self.solar_path), ("Eólica", self.wind_path),
                        ("Nuclear", self.nuclear_path), ("Térmica", self.term_path),
                        ("Intercâmbio SIN", self.interc_sin_path),
                        ("Intercâmbio interno", self.interc_interno_path),
                        ("CMO semanal", self.cmo_semanal_path),
                        ("CMO semi-horário", self.cmo_semihorario_path)]:
            print(f"  {name:20s}: {'✓' if os.path.exists(p) else '✗'} {p}")
        if self.hydro_paths:
            print(f"  Hidrelétrica       : ✓ {self.hydro_paths}")
        else:
            print(f"  Hidrelétrica       : ✗ (nenhum arquivo encontrado)")


def _load_gen_safe(path: str, fonte: str) -> pd.DataFrame:
    """Carrega geração se o arquivo existir, caso contrário retorna DataFrame vazio."""
    empty = pd.DataFrame(columns=["din_instante", "id_subsistema", "fonte", "val_geracao"])
    if not os.path.exists(path):
        print(f"  ⚠ {fonte}: arquivo não encontrado ({path})")
        return empty
    g = load_generation(path, fonte)
    print(f"  ✓ {fonte}: {len(g):,} linhas")
    return g


def build_panel(paths: SINPaths) -> pd.DataFrame:
    """Constrói o painel unificado (hora × subsistema) a partir dos dados brutos.

    Colunas do painel resultante:
        din_instante, id_subsistema, D, x_int, D_net,
        gs (solar), gr (eólica), gh (hidro),
        g_nuc_obs (nuclear), g_th_obs (térmica), gn_obs (controlável total)
    """
    print("Carregando demanda efetiva...")
    D_eff = load_demanda_efetiva(paths.curva_paths)
    print(f"  ✓ {len(D_eff):,} linhas")

    print("Carregando geração por fonte...")
    G_sol = _load_gen_safe(paths.solar_path, "solar")
    G_win = _load_gen_safe(paths.wind_path, "wind")
    G_nuc = _load_gen_safe(paths.nuclear_path, "nuclear")
    G_ter = _load_gen_safe(paths.term_path, "thermal")

    if paths.hydro_paths:
        G_hyd = _load_gen_safe(paths.hydro_paths[0], "hydro")
    else:
        print("  ⚠ hydro: nenhum arquivo encontrado — usando zeros")
        G_hyd = pd.DataFrame(columns=["din_instante", "id_subsistema", "fonte", "val_geracao"])

    G_all = pd.concat([G_sol, G_win, G_hyd, G_nuc, G_ter], ignore_index=True)
    print(f"  Total de gerações: {len(G_all):,} linhas")

    G_piv = (G_all.pivot_table(
        index=["din_instante", "id_subsistema"], columns="fonte",
        values="val_geracao", aggfunc="sum",
    ).fillna(0.0).reset_index())

    print("Carregando intercâmbio interno...")
    X_int = load_intercambio_interno(paths.interc_interno_path)
    print(f"  ✓ {len(X_int):,} linhas")

    panel = (G_piv
             .merge(D_eff, on=["din_instante", "id_subsistema"], how="outer")
             .merge(X_int, on=["din_instante", "id_subsistema"], how="left")
             .sort_values(["id_subsistema", "din_instante"])
             .reset_index(drop=True))

    for col in ["solar", "wind", "hydro", "nuclear", "thermal"]:
        if col not in panel.columns:
            panel[col] = 0.0

    panel["D"] = pd.to_numeric(panel["D"], errors="coerce")
    panel["x_int"] = panel["x_int"].fillna(0.0)

    panel["gs"] = panel["solar"].fillna(0.0)
    panel["gr"] = panel["wind"].fillna(0.0)
    panel["gh"] = panel["hydro"].fillna(0.0)
    panel["g_nuc_obs"] = panel["nuclear"].fillna(0.0)
    panel["g_th_obs"] = panel["thermal"].fillna(0.0)
    panel["gn_obs"] = panel["g_nuc_obs"] + panel["g_th_obs"]
    panel["D_net"] = panel["D"] - panel["x_int"]

    print(f"\nPainel montado: {len(panel):,} linhas")
    print(f"Subsistemas: {sorted(panel['id_subsistema'].dropna().unique())}")
    print(f"Período: {panel['din_instante'].min()} → {panel['din_instante'].max()}")
    return panel


# ---------------------------------------------------------------------------
# Visualização: Heatmaps e Curvas de Duração
# ---------------------------------------------------------------------------

def heatmap_day_hour(
    df: pd.DataFrame,
    value_col: str,
    ss: str,
    title: str,
    agg: str = "mean",
    cmap: str = "viridis",
    figsize: tuple[int, int] = (14, 4),
) -> None:
    """Plota heatmap dia × hora para uma variável e subsistema."""
    g = df[df["id_subsistema"] == ss].copy()
    if len(g) == 0:
        print(f"[{ss}] sem dados para {value_col}")
        return
    g = g.dropna(subset=[value_col])
    if len(g) == 0:
        print(f"[{ss}] {value_col} todo NaN")
        return

    g["dia"] = g["din_instante"].dt.date
    g["hora"] = g["din_instante"].dt.hour
    piv = g.pivot_table(index="hora", columns="dia", values=value_col, aggfunc=agg).sort_index()

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(piv.values, aspect="auto", interpolation="nearest", cmap=cmap)
    plt.colorbar(im, ax=ax, label=value_col)
    ax.set_yticks(range(0, 24, 2))
    ax.set_yticklabels(range(0, 24, 2))
    ax.set_title(f"[{ss}] {title}")
    ax.set_xlabel("Dia")
    ax.set_ylabel("Hora")
    fig.tight_layout()
    plt.show()


def plot_duration_curve(
    df: pd.DataFrame,
    ss: str,
    use_net: bool = True,
    figsize: tuple[int, int] = (14, 5),
) -> None:
    """Plota curva de duração (horas ordenadas por demanda decrescente)
    com geração empilhada."""
    g = df[df["id_subsistema"] == ss].copy()
    col_dem = "D_net" if use_net else "D"
    label_dem = "D líquida (D − x_int)" if use_net else "Demanda (D)"

    g = g.dropna(subset=[col_dem])
    if len(g) == 0:
        print(f"[{ss}] Sem {col_dem}")
        return

    g = g.sort_values(col_dem, ascending=False).reset_index(drop=True)
    x = np.arange(len(g))

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(x, g[col_dem].values, linewidth=2.5, label=label_dem, color="black")
    ax.fill_between(x, 0, g["gs"].values, alpha=0.6, label="Solar", color="gold")
    ax.fill_between(x, g["gs"].values,
                    g["gs"].values + g["gr"].values,
                    alpha=0.6, label="Eólica", color="steelblue")
    ax.fill_between(x, g["gs"].values + g["gr"].values,
                    g["gs"].values + g["gr"].values + g["gh"].values,
                    alpha=0.7, label="Hidrelétrica", color="cyan")
    ax.plot(x, g["gn_obs"].values, linewidth=1.8, label="gn_obs (não-renov.)", color="red")
    ax.set_title(f"[{ss}] Curva de duração (ordenada por {label_dem})")
    ax.set_xlabel("Horas ordenadas")
    ax.set_ylabel("MWmed")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    plt.show()


def plot_temporal_window(
    df: pd.DataFrame,
    ss: str,
    start: Optional[str] = None,
    end: Optional[str] = None,
    figsize: tuple[int, int] = (15, 6),
) -> None:
    """Plota janela temporal com demanda, geração e despacho otimizado."""
    g = df[df["id_subsistema"] == ss].copy().sort_values("din_instante")
    if start:
        g = g[g["din_instante"] >= pd.to_datetime(start)]
    if end:
        g = g[g["din_instante"] <= pd.to_datetime(end)]
    if len(g) == 0:
        print(f"[{ss}] Janela temporal vazia")
        return

    if "D_net" not in g.columns:
        g["D_net"] = g["D"] - g.get("x_int", 0.0) if "D" in g.columns else 0.0
    if "x_int" not in g.columns:
        g["x_int"] = 0.0

    t = g["din_instante"]
    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(t, g["D_net"], label="D líquida (D − x_int)", linewidth=2.2, color="black")
    ax.plot(t, g["x_int"], label="Intercâmbio (+ = import)", linewidth=1.2, alpha=0.7, color="gray")
    ax.plot(t, g["gs"] + g["gr"], label="Renováveis (gs + gr)", linewidth=1.8, alpha=0.8, color="orange")
    ax.plot(t, g["gh"], label="Hidrelétrica (gh)", linewidth=1.6, alpha=0.85, color="cyan")

    if "gn_hat" in g.columns:
        ax.plot(t, g["gn_hat"], label="gn_hat (modelo)", linewidth=2.0, color="blue")
    if "gn_obs" in g.columns:
        ax.plot(t, g["gn_obs"], label="gn_obs (observado)", linewidth=1.5, alpha=0.6, color="red", linestyle="--")
    if "u_hat" in g.columns:
        ax.fill_between(t, 0, g["u_hat"], alpha=0.3, color="red", label="Déficit (u_hat)")
    if "c_hat" in g.columns:
        ax.fill_between(t, 0, g["c_hat"], alpha=0.3, color="yellow", label="Curtailment (c_hat)")

    ax.set_title(f"[{ss}] Despacho social — {t.min().date()} a {t.max().date()}")
    ax.set_xlabel("Data/Hora")
    ax.set_ylabel("Potência (MWmed)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", ncol=2)
    fig.tight_layout()
    plt.show()
