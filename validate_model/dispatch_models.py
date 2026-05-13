"""
Módulo central de modelagem hidrotérmica — SIN 2025.

Centraliza:
  • Dataclasses de configuração (anual, diário) — sem mágica entre células.
  • Solvers QP (CVXPY) com extração explícita de duais.
  • Auditoria de balanço de massa/energia (asserts e warnings).
  • Regressões de cotas (montante/jusante) com diagnóstico.
  • Permutações genéricas (origem/destino, variáveis a permutar, lag).
  • Lucro solar e canibalização (NB5).

Convenções de sinal (NB4 — anual)
---------------------------------
Balanço de potência:
    g^th_t = D_t - g^h_t - g^ren_t + c_t - d_t

  onde c_t (curtailment) entra como **carga fictícia adicional** (>=0)
  e d_t (déficit) entra como **geração fictícia** (>=0).

Dinâmica do reservatório (com vertimento ativo):
    V_{t+1} = V_t + A_t - g^h_t - spill_t

V_max é tratado como **constante fixa** lida do dataclass de configuração;
nunca é recalibrado entre células.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field, replace
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

try:
    import cvxpy as cp
except ImportError:  # pragma: no cover - cvxpy é dependência declarada
    cp = None  # type: ignore[assignment]

from .pipeline import read_csv_robust, to_float, detect_datetime_col


# ---------------------------------------------------------------------------
# Loader auxiliar: hidrologia bruta (linha por reservatório-dia)
# ---------------------------------------------------------------------------


def load_hidrologia_raw(path: str) -> pd.DataFrame:
    """Carrega hidrologia diária **sem agregação**, para regressões de cota.

    Mantém uma linha por (reservatório, dia) com colunas numéricas limpas:
        val_nivelmontante, val_niveljusante,
        val_vazaoafluente, val_vazaoturbinada, val_vazaovertida, val_vazaodefluente,
        val_volumeutilcon, head, din_instante.
    """
    df = read_csv_robust(path)
    dtc = detect_datetime_col(df)
    if dtc is None:
        raise ValueError(f"[{path}] coluna datetime não encontrada.")
    df["din_instante"] = pd.to_datetime(df[dtc], errors="coerce")

    num_cols = [
        "val_nivelmontante", "val_niveljusante",
        "val_vazaoafluente", "val_vazaoturbinada",
        "val_vazaovertida", "val_vazaodefluente",
        "val_volumeutilcon",
    ]
    for c in num_cols:
        if c in df.columns:
            if df[c].dtype == object:
                df[c] = to_float(df[c])
            df[c] = pd.to_numeric(df[c], errors="coerce")

    if {"val_nivelmontante", "val_niveljusante"}.issubset(df.columns):
        df["head"] = df["val_nivelmontante"] - df["val_niveljusante"]
    return df


# ---------------------------------------------------------------------------
# Config — Despacho anual (NB4)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AnnualParams:
    """Parâmetros do despacho hidrotérmico anual (reservatório agregado).

    Todos os parâmetros são imutáveis (`frozen=True`) para impedir
    recalibrações acidentais entre células do notebook.
    """

    # Capacidades (MW)
    Gh_max: float
    Gt_max: float

    # Reservatório agregado (MWh) — CONSTANTE FIXA
    V_max: float

    # Estado inicial (fração de V_max)
    alpha: float = 0.30

    # Custos térmicos (R$/MWh, R$/MW²h)
    c1: float = 50.0
    c2: float = 1e-3

    # Penalidades
    pi_d: float = 1e4   # déficit (R$/MWh) — alta
    pi_c: float = 10.0  # curtailment (R$/MWh) — baixa
    pi_s: float = 1.0   # vertimento (R$/MWh) — baixa

    # Condição terminal
    # terminal_mode ∈ {"exact", "band"}
    #   "exact": V_T = beta_target * V_max  (com tolerância terminal_tol)
    #   "band":  beta_min*V_max ≤ V_T ≤ beta_max*V_max
    terminal_mode: str = "band"
    beta_target: float = 0.30
    beta_min: float = 0.20
    beta_max: float = 0.95
    terminal_tol: float = 0.005

    def __post_init__(self) -> None:
        if self.terminal_mode not in {"exact", "band"}:
            raise ValueError("terminal_mode deve ser 'exact' ou 'band'.")
        if not 0.0 < self.alpha < 1.0:
            raise ValueError("alpha deve estar em (0, 1).")
        if self.V_max <= 0:
            raise ValueError("V_max deve ser positivo.")
        if self.terminal_mode == "exact" and not 0.0 < self.beta_target < 1.0:
            raise ValueError("beta_target deve estar em (0, 1).")
        if self.terminal_mode == "band" and not 0.0 <= self.beta_min < self.beta_max <= 1.0:
            raise ValueError("Banda terminal inválida (beta_min < beta_max).")

    def with_overrides(self, **kwargs) -> "AnnualParams":
        """Retorna nova instância com substituições explícitas — útil em
        análises de sensibilidade sem mutar o objeto-base."""
        return replace(self, **kwargs)


# ---------------------------------------------------------------------------
# Config — Despacho diário (NB5)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DailyParams:
    """Parâmetros do despacho diário (24 h)."""

    # Capacidades nominais
    K_S1: float = 10_000.0   # MW — pico do painel solar 1
    K_S2: float = 16_000.0   # MW — pico do painel solar 2 (remoto)
    Gt_max: float = 50_000.0  # MW — capacidade térmica
    hydro_scale: float = 1.0  # multiplicador no orçamento diário hidro

    # Custo térmico quadrático: a*gt^2 + b*gt   [R$/MW²h, R$/MWh]
    # Defaults calibrados para que CMg = 2a·gt + b fique na faixa 100-600
    # R$/MWh para gt entre 15-30 GW (compatível com CMO histórico do SIN).
    a: float = 1.0e-2
    b: float = 50.0

    # Custos de oportunidade do solar (R$/MWh) — LCOE de referência
    cost_s1: float = 80.0
    cost_s2: float = 50.0

    # Penalidades (R$/MWh)
    pi_shed: float = 5.0e3   # alta — limite superior implícito do CMO
    pi_curt: float = 5.0     # baixa — evita curtailment desnecessário

    # Transmissão para Solar 2 (remoto)
    line_limit: float = float("inf")
    line_loss: float = 0.0


# ---------------------------------------------------------------------------
# Estimador OFFLINE de V_max (apenas diagnóstico)
# ---------------------------------------------------------------------------


def estimate_V_max_diagnostic(hydro_daily: pd.DataFrame, min_dV_pct: float = 2.0) -> dict:
    """Estima V_max a partir do balanço observado.

    NÃO use o resultado para alimentar o solver — V_max do modelo deve ser
    tratado como **constante fixa**. Esta função existe apenas para diagnosticar
    se a magnitude assumida está em ordem de grandeza compatível com os dados.

    Retorna dict com: V_max_global, V_max_amplitude, status, residuals.
    """
    hd = hydro_daily.dropna(subset=["vol_util_pct", "A_MW", "T_MW", "Spill_MW"]).copy()
    if hd.empty:
        return {"status": "no_data"}

    net_daily_MWh = (hd["A_MW"] - hd["T_MW"] - hd["Spill_MW"]).to_numpy() * 24.0
    E_cum = np.cumsum(net_daily_MWh)
    dV_frac = (hd["vol_util_pct"].to_numpy() - hd["vol_util_pct"].iloc[0]) / 100.0

    delta_E = float(E_cum[-1])
    delta_V = float(dV_frac[-1])
    amp_E = float(E_cum.max() - E_cum.min())
    amp_V = float(dV_frac.max() - dV_frac.min())

    if abs(delta_V) > min_dV_pct / 100.0:
        v_global = abs(delta_E / delta_V)
    else:
        v_global = float("nan")
    v_amp = abs(amp_E / amp_V) if amp_V > 0 else float("nan")

    return {
        "V_max_global_MWh": v_global,
        "V_max_amplitude_MWh": v_amp,
        "delta_V_frac": delta_V,
        "amp_V_frac": amp_V,
        "status": "ok",
    }


# ---------------------------------------------------------------------------
# Solver anual (NB4) — QP convexo via CVXPY
# ---------------------------------------------------------------------------


def _solve_with_fallbacks(prob, verbose: bool = False, prefer_duals: bool = True) -> str:
    """Tenta resolver `prob` com solvers que entregam **duais confiáveis**.

    Ordem padrão (`prefer_duals=True`): CLARABEL → ECOS → SCS → OSQP.
    Solvers IPM (CLARABEL/ECOS) entregam multiplicadores Lagrangianos
    consistentes na otimalidade; OSQP é mantido apenas como fallback.
    """
    if cp is None:
        raise ImportError("cvxpy não está disponível.")

    solvers: list = []
    if prefer_duals:
        if hasattr(cp, "CLARABEL"):
            solvers.append(cp.CLARABEL)
        if hasattr(cp, "ECOS"):
            solvers.append(cp.ECOS)
        solvers.extend([cp.SCS, cp.OSQP])
    else:
        solvers = [cp.OSQP, cp.SCS]
        if hasattr(cp, "CLARABEL"):
            solvers.append(cp.CLARABEL)
        if hasattr(cp, "ECOS"):
            solvers.append(cp.ECOS)

    last_status = "unknown"
    for solver in solvers:
        try:
            kw = {"verbose": verbose}
            if solver == cp.SCS:
                kw.update({"max_iters": 200_000, "eps": 1e-7})
            prob.solve(solver=solver, **kw)
            last_status = prob.status
            if prob.status in ("optimal", "optimal_inaccurate"):
                return prob.status
        except Exception as exc:  # pragma: no cover — log defensivo
            last_status = f"error:{type(exc).__name__}"
            continue
    return last_status


def solve_annual_dispatch(
    D: np.ndarray,
    g_ren: np.ndarray,
    A: np.ndarray,
    params: AnnualParams,
    verbose: bool = False,
) -> dict:
    """Resolve o despacho hidrotérmico anual.

    Convenções:
        g^th_t = D_t - g^h_t - g^ren_t + c_t - d_t
        V_{t+1} = V_t + A_t - g^h_t - spill_t

    Retorna dict com séries (gh, gt, c, d, spill, V), valor da água,
    status, custo total e duals utilizados pelo módulo de auditoria.
    """
    if cp is None:
        raise ImportError("cvxpy não está disponível.")

    D = np.asarray(D, dtype=float)
    g_ren = np.asarray(g_ren, dtype=float)
    A = np.asarray(A, dtype=float)
    T = len(D)
    if not (len(g_ren) == len(A) == T):
        raise ValueError("D, g_ren e A precisam ter o mesmo tamanho.")

    # Normalizamos V em [0,1] (V_norm = V_phys / V_max) e (A-gh-spill) em
    # frações de V_max — equivale ao problema original e gera matriz de
    # restrições bem-condicionada, mantendo a recuperação do dual em R$/MWh.
    inv_Vmax = 1.0 / params.V_max

    gh = cp.Variable(T, name="gh")
    c_var = cp.Variable(T, name="c")
    d_var = cp.Variable(T, name="d")
    spill = cp.Variable(T, name="spill")
    V_norm = cp.Variable(T + 1, name="V_norm")  # adimensional ∈ [0, 1]

    # Balanço (variável afim de g^th_t — nunca declarado como variável solta)
    gt = D - gh - g_ren + c_var - d_var

    cost = (
        params.c1 * cp.sum(gt)
        + 0.5 * params.c2 * cp.sum_squares(gt)
        + params.pi_d * cp.sum(d_var)
        + params.pi_c * cp.sum(c_var)
        + params.pi_s * cp.sum(spill)
    )

    # Para fechar o ciclo do dual: ∂L/∂gh_t = -(c1+c2·gt) + λ_t · inv_Vmax = 0
    # ⇒ λ_t = (c1 + c2·gt) / inv_Vmax = (c1 + c2·gt) · V_max.  O valor da água
    # em R$/MWh é então water_value = λ_t · inv_Vmax.
    reservoir_dyn = V_norm[1:] == V_norm[:-1] + (A - gh - spill) * inv_Vmax

    cons = [
        gh >= 0, gh <= params.Gh_max,
        gt >= 0, gt <= params.Gt_max,
        c_var >= 0, d_var >= 0, spill >= 0,
        V_norm >= 0.0, V_norm <= 1.0,
        V_norm[0] == params.alpha,
        reservoir_dyn,
    ]

    if params.terminal_mode == "exact":
        cons += [
            V_norm[-1] >= params.beta_target - params.terminal_tol,
            V_norm[-1] <= params.beta_target + params.terminal_tol,
        ]
    else:
        cons += [
            V_norm[-1] >= params.beta_min,
            V_norm[-1] <= params.beta_max,
        ]

    prob = cp.Problem(cp.Minimize(cost), cons)
    status = _solve_with_fallbacks(prob, verbose=verbose)

    def _val(x, n):
        if x.value is None:
            return np.full(n, np.nan)
        return np.asarray(x.value).ravel()

    v_frac = _val(V_norm, T + 1)
    V_vec = v_frac * params.V_max

    water_value = None
    if reservoir_dyn.dual_value is not None:
        # CVXPY devolve λ_t referente à equação adimensional.  Para R$/MWh
        # multiplicamos por inv_Vmax (ver derivação acima).  abs() padroniza
        # o sinal entre solvers (CLARABEL/ECOS/SCS adotam convenções distintas).
        lam_norm = np.asarray(reservoir_dyn.dual_value, dtype=float).ravel()
        water_value = np.abs(lam_norm * inv_Vmax)

    gt_val = (D - _val(gh, T) - g_ren + _val(c_var, T) - _val(d_var, T)) if status.startswith("optimal") else np.full(T, np.nan)

    return {
        "status": status,
        "total_cost": float(prob.value) if prob.value is not None else np.nan,
        "gh": _val(gh, T),
        "gt": gt_val,
        "c": _val(c_var, T),
        "d": _val(d_var, T),
        "spill": _val(spill, T),
        "V": V_vec,
        "v_frac": v_frac,
        "water_value": water_value,
        "params": params,
    }


# ---------------------------------------------------------------------------
# Auditoria — balanço de massa/energia
# ---------------------------------------------------------------------------


def audit_annual_solution(
    result: Mapping,
    A: np.ndarray,
    D: np.ndarray | None = None,
    g_ren: np.ndarray | None = None,
    tol: float = 1e-6,
    raise_on_fail: bool = False,
) -> pd.DataFrame:
    """Verifica balanço de massa (água) e potência da solução anual.

    Critérios:
      • mass:    sum(A) - sum(gh) - sum(spill) - (V_T - V_0) ≈ 0
      • power:   max|D - gh - g_ren - gt + c - d| ≈ 0  (resíduo do balanço hora-a-hora)

    Retorna DataFrame de diagnóstico; se algum erro relativo exceder `tol`,
    emite warning (ou levanta `AssertionError` se raise_on_fail=True).
    """
    A = np.asarray(A, dtype=float)
    gh = np.asarray(result["gh"], dtype=float)
    spill = np.asarray(result["spill"], dtype=float)
    V = np.asarray(result["V"], dtype=float)

    sum_A = float(np.nansum(A))
    sum_gh = float(np.nansum(gh))
    sum_spill = float(np.nansum(spill))
    delta_V = float(V[-1] - V[0])
    mass_residual = sum_A - sum_gh - sum_spill - delta_V
    mass_rel = abs(mass_residual) / max(sum_A, 1.0)

    rows = [{
        "check": "mass_balance",
        "value": mass_residual,
        "scale": sum_A,
        "rel_error": mass_rel,
        "ok": mass_rel < tol,
    }]

    if D is not None and g_ren is not None:
        D_arr = np.asarray(D, dtype=float)
        g_ren_arr = np.asarray(g_ren, dtype=float)
        gt = np.asarray(result["gt"], dtype=float)
        c = np.asarray(result["c"], dtype=float)
        d = np.asarray(result["d"], dtype=float)
        residual_h = D_arr - gh - g_ren_arr - gt + c - d
        max_abs = float(np.nanmax(np.abs(residual_h)))
        rel = max_abs / max(np.nanmean(np.abs(D_arr)), 1.0)
        rows.append({
            "check": "power_balance_max",
            "value": max_abs,
            "scale": float(np.nanmean(np.abs(D_arr))),
            "rel_error": rel,
            "ok": rel < tol,
        })

    audit = pd.DataFrame(rows)
    bad = audit[~audit["ok"]]
    if len(bad) > 0:
        msg = (
            "Auditoria de balanço falhou:\n"
            + bad.to_string(index=False)
        )
        if raise_on_fail:
            raise AssertionError(msg)
        warnings.warn(msg, stacklevel=2)
    return audit


# ---------------------------------------------------------------------------
# Regressões de cotas — diagnóstico de qualidade do dataset
# ---------------------------------------------------------------------------


def _rsquared(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ss_res = np.nansum((y_true - y_pred) ** 2)
    ss_tot = np.nansum((y_true - np.nanmean(y_true)) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def _rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.sqrt(np.nanmean((np.asarray(y_true) - np.asarray(y_pred)) ** 2)))


def fit_polynomial_regression(
    x: np.ndarray,
    y: np.ndarray,
    degree: int = 2,
) -> dict:
    """Regressão polinomial 1D com diagnóstico (R², RMSE, resíduos)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x_f, y_f = x[mask], y[mask]
    if len(x_f) <= degree + 1:
        return {"status": "insufficient_data", "n": int(mask.sum())}
    coeffs = np.polyfit(x_f, y_f, degree)
    y_hat = np.polyval(coeffs, x_f)
    return {
        "status": "ok",
        "n": int(mask.sum()),
        "coeffs": coeffs,
        "degree": degree,
        "r2": _rsquared(y_f, y_hat),
        "rmse": _rmse(y_f, y_hat),
        "x": x_f,
        "y": y_f,
        "y_hat": y_hat,
        "residuals": y_f - y_hat,
    }


def cota_regressions(
    df_raw: pd.DataFrame,
    *,
    montante_col: str = "val_nivelmontante",
    jusante_col: str = "val_niveljusante",
    volume_col: str = "val_volumeutilcon",
    q_turb_col: str = "val_vazaoturbinada",
    q_vert_col: str = "val_vazaovertida",
    q_def_col: str = "val_vazaodefluente",
    reservoir_col: str | None = "nom_reservatorio",
    degree_h_m: int = 2,
    degree_h_j: int = 2,
    min_n_per_reservoir: int = 30,
) -> dict:
    """Ajusta, **por reservatório** (e também agregado):
        H_m = f_m(V)              — cota montante × volume útil
        H_j = f_j(Q_def)          — cota jusante × vazão defluente

    Cada reservatório tem sua própria curva característica (cotas absolutas
    em metros e volumes nominais distintos). O ajuste agregado existe apenas
    como controle; a interpretação válida está na **distribuição dos R²**
    por reservatório.

    Retorna dict com:
        {
          "montante_vs_volume": {"global": fit_global, "per_reservoir": DataFrame},
          "jusante_vs_qdef":    {"global": fit_global, "per_reservoir": DataFrame},
        }
    """
    cols_required = [montante_col, jusante_col]
    missing = [c for c in cols_required if c not in df_raw.columns]
    if missing:
        raise KeyError(f"Colunas ausentes para regressão de cotas: {missing}")

    df = df_raw.copy()

    if q_def_col in df.columns:
        df["__q_def"] = pd.to_numeric(df[q_def_col], errors="coerce")
    else:
        df["__q_def"] = (
            pd.to_numeric(df.get(q_turb_col, np.nan), errors="coerce").fillna(0)
            + pd.to_numeric(df.get(q_vert_col, np.nan), errors="coerce").fillna(0)
        )

    has_volume = volume_col in df.columns

    fit_m_global = fit_polynomial_regression(
        df[volume_col].to_numpy(dtype=float, na_value=np.nan) if has_volume
        else np.full(len(df), np.nan),
        df[montante_col].to_numpy(dtype=float, na_value=np.nan),
        degree=degree_h_m,
    )
    fit_j_global = fit_polynomial_regression(
        df["__q_def"].to_numpy(dtype=float, na_value=np.nan),
        df[jusante_col].to_numpy(dtype=float, na_value=np.nan),
        degree=degree_h_j,
    )

    per_m_rows: list[dict] = []
    per_j_rows: list[dict] = []
    if reservoir_col is not None and reservoir_col in df.columns:
        for name, sub in df.groupby(reservoir_col, dropna=True):
            if has_volume:
                fit = fit_polynomial_regression(
                    sub[volume_col].to_numpy(dtype=float, na_value=np.nan),
                    sub[montante_col].to_numpy(dtype=float, na_value=np.nan),
                    degree=degree_h_m,
                )
                if fit.get("status") == "ok" and fit["n"] >= min_n_per_reservoir:
                    per_m_rows.append({
                        "reservatorio": name, "n": fit["n"],
                        "r2": fit["r2"], "rmse": fit["rmse"],
                    })
            fit = fit_polynomial_regression(
                sub["__q_def"].to_numpy(dtype=float, na_value=np.nan),
                sub[jusante_col].to_numpy(dtype=float, na_value=np.nan),
                degree=degree_h_j,
            )
            if fit.get("status") == "ok" and fit["n"] >= min_n_per_reservoir:
                per_j_rows.append({
                    "reservatorio": name, "n": fit["n"],
                    "r2": fit["r2"], "rmse": fit["rmse"],
                })
    per_m = pd.DataFrame(per_m_rows).sort_values("r2", ascending=False) if per_m_rows else pd.DataFrame()
    per_j = pd.DataFrame(per_j_rows).sort_values("r2", ascending=False) if per_j_rows else pd.DataFrame()

    return {
        "montante_vs_volume": {"global": fit_m_global, "per_reservoir": per_m},
        "jusante_vs_qdef":   {"global": fit_j_global, "per_reservoir": per_j},
    }


# ---------------------------------------------------------------------------
# Permutações — análise de cenários (NB4 e NB5)
# ---------------------------------------------------------------------------


def shift_circular(x: np.ndarray, hours: int) -> np.ndarray:
    return np.roll(np.asarray(x, dtype=float), int(hours))


def run_permutation_scenario(
    sin: pd.DataFrame,
    *,
    params: AnnualParams,
    variables_to_swap: Mapping[str, int] | None = None,
    origin_day: int | None = None,
    target_day: int | None = None,
) -> dict:
    """Executa uma permutação genérica antes de resolver o despacho anual.

    Dois modos de uso:
      1. `variables_to_swap = {"A_MW": 10, "gs": -10, ...}`
         Aplica deslocamento circular em **dias** para cada coluna.
         Equivalente ao "deslocamento sistemático" de k dias.

      2. `origin_day, target_day`: copia o perfil de 24 h de `origin_day`
         para o lugar de `target_day` em todas as variáveis listadas em
         `variables_to_swap` (chaves apenas; valores ignorados).
    """
    sin = sin.sort_values("din_instante").reset_index(drop=True)
    df = sin.copy()
    variables_to_swap = dict(variables_to_swap or {})

    if origin_day is not None and target_day is not None:
        if not variables_to_swap:
            raise ValueError("Informe variables_to_swap quando origin_day/target_day forem dados.")
        for col in variables_to_swap.keys():
            if col not in df.columns:
                continue
            origin_idx = (df["din_instante"].dt.normalize() == pd.Timestamp(origin_day))
            target_idx = (df["din_instante"].dt.normalize() == pd.Timestamp(target_day))
            if origin_idx.sum() != target_idx.sum() or origin_idx.sum() == 0:
                raise ValueError(f"Days {origin_day}/{target_day} sem 24 h em '{col}'.")
            df.loc[target_idx.values, col] = df.loc[origin_idx.values, col].to_numpy()
    else:
        for col, days in variables_to_swap.items():
            if col not in df.columns:
                continue
            df[col] = shift_circular(df[col].to_numpy(dtype=float), days * 24)

    g_ren_perm = df.get("gs", pd.Series(np.zeros(len(df)))).to_numpy(dtype=float) + \
                 df.get("gr", pd.Series(np.zeros(len(df)))).to_numpy(dtype=float)
    res = solve_annual_dispatch(
        D=df["D"].to_numpy(dtype=float),
        g_ren=g_ren_perm,
        A=df["A_MW"].to_numpy(dtype=float),
        params=params,
    )
    res["sin_perm"] = df
    res["spec"] = {
        "variables_to_swap": dict(variables_to_swap),
        "origin_day": origin_day,
        "target_day": target_day,
    }
    return res


# ---------------------------------------------------------------------------
# Solver diário (NB5) com extração explícita de duais
# ---------------------------------------------------------------------------


def make_solar_profiles(
    day_supply: pd.DataFrame,
    *,
    K_S1: float,
    K_S2: float,
    solar_col: str = "gs",
) -> tuple[np.ndarray, np.ndarray]:
    """Constrói perfis horários (24 h) dos dois painéis a partir do dia base."""
    raw = day_supply[solar_col].to_numpy(dtype=float)
    if len(raw) != 24:
        raise ValueError("Dia base precisa ter exatamente 24 horas.")
    base_cf = raw / max(raw.max(), 1.0)
    cf1 = np.clip(base_cf, 0.0, 1.0)
    cf2 = np.clip(0.90 * np.roll(base_cf, 1) + 0.10 * base_cf, 0.0, 1.0)
    return K_S1 * cf1, K_S2 * cf2


def solve_daily_dispatch(
    day_demand: pd.DataFrame,
    day_supply: pd.DataFrame,
    params: DailyParams,
    *,
    verbose: bool = False,
) -> dict:
    """Resolve o despacho diário com:

        min  Σ_t [ a*gt² + b*gt + π_shed*shed + π_curt*curt ]
        s.a. gt + gh + gs1 + gs2 + shed - curt = D
             0 ≤ gh, gt; capacidades; orçamento hidro;
             gs2 limitada por linha (line_limit) e perda (line_loss).

    Retorna o dispatch hora-a-hora **e** o preço-sombra (λ_t) da restrição
    de balanço de potência — o CMO endógeno do problema.
    """
    if cp is None:
        raise ImportError("cvxpy não está disponível.")

    if len(day_demand) != 24:
        raise ValueError("day_demand precisa ter exatamente 24 linhas (horárias).")

    D = day_demand["D"].to_numpy(dtype=float)
    gs1_raw, gs2_raw = make_solar_profiles(
        day_supply, K_S1=params.K_S1, K_S2=params.K_S2,
    )
    # Restrições de transmissão para o painel remoto (Solar 2)
    gs2_export = np.minimum(gs2_raw, params.line_limit)
    gs2_delivered = (1.0 - params.line_loss) * gs2_export

    hydro_budget = float(day_supply["gh"].sum() * params.hydro_scale)
    gh_cap = float(max(day_supply["gh"].max() * 1.05, hydro_budget / 24.0 * 1.10))

    gt = cp.Variable(24, nonneg=True, name="gt")
    gh = cp.Variable(24, nonneg=True, name="gh")
    shed = cp.Variable(24, nonneg=True, name="shed")
    curt = cp.Variable(24, nonneg=True, name="curt")

    solar_total = gs1_raw + gs2_delivered

    balance = gt + gh + solar_total + shed - curt == D
    cons = [
        gt <= params.Gt_max,
        gh <= gh_cap,
        cp.sum(gh) <= hydro_budget,
        curt <= solar_total,
        balance,
    ]

    obj = cp.Minimize(
        params.a * cp.sum_squares(gt) + params.b * cp.sum(gt)
        + params.pi_shed * cp.sum(shed) + params.pi_curt * cp.sum(curt)
    )
    prob = cp.Problem(obj, cons)
    status = _solve_with_fallbacks(prob, verbose=verbose)

    def _val(x):
        return np.asarray(x.value).ravel() if x.value is not None else np.full(24, np.nan)

    gt_val = _val(gt)
    gh_val = _val(gh)
    shed_val = _val(shed)
    curt_val = _val(curt)

    # λ_t — preço-sombra do balanço
    if balance.dual_value is not None:
        lam = np.asarray(balance.dual_value).ravel()
        # CVXPY retorna o sinal coerente com a derivada do Lagrangiano:
        # para uma restrição de igualdade, λ é o ganho marginal de
        # relaxar o lado direito (carga). Mantemos como vem.
        lam = np.abs(lam)
    else:
        lam = 2.0 * params.a * gt_val + params.b  # fallback teórico

    cmo_theoretical = 2.0 * params.a * gt_val + params.b

    out = pd.DataFrame({
        "din_instante": day_demand["din_instante"].to_numpy(),
        "D": D,
        "gs1": gs1_raw,
        "gs2": gs2_delivered,
        "gs2_raw": gs2_raw,
        "gs2_export": gs2_export,
        "gh": gh_val,
        "gt": gt_val,
        "shed": shed_val,
        "curt": curt_val,
        "lambda": lam,
        "cmo_theoretical": cmo_theoretical,
    })

    return {
        "status": status,
        "objective": float(prob.value) if prob.value is not None else np.nan,
        "dispatch": out,
        "hydro_budget_GWh": hydro_budget / 1e3,
        "solar2_spilled_GWh": float(np.sum(gs2_raw - gs2_export) / 1e3),
        "params": params,
    }


# ---------------------------------------------------------------------------
# Lucratividade solar e canibalização (NB5)
# ---------------------------------------------------------------------------


def solar_profit(
    dispatch: pd.DataFrame,
    *,
    cost_s1: float = 0.0,
    cost_s2: float = 0.0,
    lambda_col: str = "lambda",
) -> dict:
    """Lucro horário e total dos dois painéis solares."""
    lam = dispatch[lambda_col].to_numpy(dtype=float)
    pi_s1 = (lam - cost_s1) * dispatch["gs1"].to_numpy(dtype=float)
    pi_s2 = (lam - cost_s2) * dispatch["gs2"].to_numpy(dtype=float)
    return {
        "pi_s1_total": float(np.nansum(pi_s1)),
        "pi_s2_total": float(np.nansum(pi_s2)),
        "pi_s1_hourly": pi_s1,
        "pi_s2_hourly": pi_s2,
        "lambda_mean": float(np.nanmean(lam)),
        "lambda_p95": float(np.nanquantile(lam, 0.95)),
    }


def cannibalization_sweep(
    day_demand: pd.DataFrame,
    day_supply: pd.DataFrame,
    base_params: DailyParams,
    *,
    K_S2_grid: Sequence[float],
) -> pd.DataFrame:
    """Varre `K_S2` (capacidade do painel 2) e mede como o lucro do painel 1 cai."""
    rows = []
    for k2 in K_S2_grid:
        p = replace(base_params, K_S2=float(k2))
        sol = solve_daily_dispatch(day_demand, day_supply, p)
        disp = sol["dispatch"]
        prof = solar_profit(disp, cost_s1=p.cost_s1, cost_s2=p.cost_s2)
        rows.append({
            "K_S2_MW": float(k2),
            "lambda_mean": prof["lambda_mean"],
            "pi_s1_total": prof["pi_s1_total"],
            "pi_s2_total": prof["pi_s2_total"],
            "gt_mean": float(disp["gt"].mean()),
            "curt_GWh": float(disp["curt"].sum() / 1e3),
            "status": sol["status"],
        })
    return pd.DataFrame(rows)


def hourly_permutation(
    day_demand: pd.DataFrame,
    day_supply: pd.DataFrame,
    base_params: DailyParams,
    *,
    permutation_label: str,
) -> dict:
    """Combina demanda de um dia com oferta de outro (24 h)."""
    sol = solve_daily_dispatch(day_demand, day_supply, base_params)
    disp = sol["dispatch"]
    prof = solar_profit(disp, cost_s1=base_params.cost_s1, cost_s2=base_params.cost_s2)
    return {
        "label": permutation_label,
        "dispatch": disp,
        "lambda_mean": prof["lambda_mean"],
        "lambda_peak": float(np.nanmax(disp["lambda"].to_numpy())),
        "gt_mean": float(disp["gt"].mean()),
        "gt_peak": float(disp["gt"].max()),
        "shed_MWh": float(disp["shed"].sum()),
        "curt_GWh": float(disp["curt"].sum() / 1e3),
        "pi_s1_total": prof["pi_s1_total"],
        "pi_s2_total": prof["pi_s2_total"],
        "status": sol["status"],
    }


__all__ = [
    "AnnualParams",
    "DailyParams",
    "load_hidrologia_raw",
    "estimate_V_max_diagnostic",
    "solve_annual_dispatch",
    "audit_annual_solution",
    "fit_polynomial_regression",
    "cota_regressions",
    "shift_circular",
    "run_permutation_scenario",
    "make_solar_profiles",
    "solve_daily_dispatch",
    "solar_profit",
    "cannibalization_sweep",
    "hourly_permutation",
]
