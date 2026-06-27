"""
Diagnostic experiments for T-zone enrichment mechanisms.

This script is designed to be run next to one of the user's simulation files, for example:
    mortcellulaire_Tzone_sweep_lambda_d_beta_d_only(4).py
    mortcellulaire_feedback_ERK_1000_stats_unified_burnin_Tzone_binomial(6).py

It performs three diagnostic experiments:

Experiment 1A: fixed death-count lambda_d sweep
    - fixed activation regime: lambda_a_T=0.5, lambda_a_c=0.005, beta_a_R=2.5, beta_a_T=3.0
    - ERK regimes:
        strong_ERK_protection: beta_d_R=1.0, beta_d_T=0.4
        reference_ERK:         beta_d_R=2.0, beta_d_T=0.8
        weak_ERK_protection:   beta_d_R=4.0, beta_d_T=1.6
    - sweep lambda_d = 0.25, 0.5, 1.0, 2.0, 4.0
    - stopping rule: burn-in 500 deaths, analyze next 1000 deaths

Experiment 1B: fixed physical-time lambda_d control
    - same parameter grid as Experiment 1A
    - first estimates the reference physical window from Experiment 1A at
      lambda_d=1.0, beta_d_R=2.0, beta_d_T=0.8
    - then after burn-in, observes all deaths during the same physical time window
      for every lambda_d and ERK regime

Experiment 2: death/ERK robustness conditional on activation localization
    - activation regimes:
        non_local: beta_a_R=1.0, beta_a_T=1.2
        medium:    beta_a_R=2.5, beta_a_T=3.0
        local:     beta_a_R=5.0, beta_a_T=3.0
      with lambda_a_T=0.5 and lambda_a_c=0.005 fixed in all regimes.
    - in each activation regime, sweep lambda_d and the three ERK regimes.
    - performs pooled binomial tests and Bonferroni correction over all settings.

Outputs are written to:
    diagnostic_outputs_Tzone_mechanisms/

Recommended workflow:
    1. First set N_REPLICATES = 3 or 5 to test runtime.
    2. For final results, use N_REPLICATES = 20 if runtime allows.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Tuple
import importlib.util
import math
import sys

import numpy as np
import pandas as pd

try:
    from scipy.stats import binomtest, spearmanr
except Exception as exc:  # pragma: no cover
    raise RuntimeError("This diagnostic script requires scipy: binomtest and spearmanr.") from exc


# ============================================================
# 0. User-adjustable settings
# ============================================================

# For a quick smoke test use 3 or 5. For final results use 20 if runtime allows.
N_REPLICATES = 10
SEED0 = 7001

BURN_IN_DEATHS = 500
TARGET_DEATHS = 1000
MAX_PROPOSALS = 120_000_000

LAMBDA_D_VALUES = [0.25, 0.5, 1.0, 2.0, 4.0]

ERK_REGIMES = [
    {
        "erk_regime": "strong_ERK_protection",
        "beta_d_R": 1.0,
        "beta_d_T": 0.4,
        "meaning": "large ERK radius, long ERK lifetime",
    },
    {
        "erk_regime": "reference_ERK",
        "beta_d_R": 2.0,
        "beta_d_T": 0.8,
        "meaning": "reference ERK protection",
    },
    {
        "erk_regime": "weak_ERK_protection",
        "beta_d_R": 4.0,
        "beta_d_T": 1.6,
        "meaning": "small ERK radius, short ERK lifetime",
    },
]

ACTIVATION_REGIMES = [
    {
        "activation_regime": "non_local_activation",
        "beta_a_R": 1.0,
        "beta_a_T": 1.2,
        "meaning": "large and long-lived active disks",
    },
    {
        "activation_regime": "medium_activation",
        "beta_a_R": 2.5,
        "beta_a_T": 3.0,
        "meaning": "current visible reference regime",
    },
    {
        "activation_regime": "local_activation",
        "beta_a_R": 5.0,
        "beta_a_T": 3.0,
        "meaning": "small and short-lived active disks",
    },
]

# Fixed T-zone activation contrast for these diagnostic experiments.
FIXED_LAMBDA_A_T = 0.5
FIXED_LAMBDA_A_C = 0.005

OUTPUT_DIR = Path("diagnostic_outputs_Tzone_mechanisms")

# The script will try to import the first existing file from this list.
SIMULATION_FILE_CANDIDATES = [
    "mortcellulaire_Tzone_sweep_lambda_d_beta_d_only.py",
    "mortcellulaire_Tzone_sweep_lambda_d_beta_d_only(4).py",
    "mortcellulaire_Tzone_sweep_lambda_d_beta_d_only (4).py",
    "mortcellulaire_Tzone_sweep_lambda_d_beta_d_only(3).py",
    "mortcellulaire_Tzone_sweep_lambda_d_beta_d_only (3).py",
    "mortcellulaire_Tzone_sweep_lambda_d_beta_d_only(2).py",
    "mortcellulaire_Tzone_sweep_lambda_d_beta_d_only (1).py",
    "mortcellulaire_feedback_ERK_1000_stats_unified_burnin_Tzone_binomial.py",
    "mortcellulaire_feedback_ERK_1000_stats_unified_burnin_Tzone_binomial(6).py",
    "mortcellulaire_feedback_ERK_1000_stats_unified_burnin_Tzone_binomial (6).py",
]


# ============================================================
# 1. Load the user's simulation module
# ============================================================


def load_simulation_module() -> Any:
    here = Path(__file__).resolve().parent
    candidates = []
    for name in SIMULATION_FILE_CANDIDATES:
        candidates.append(here / name)
    # Also search the current working directory, in case the script is run elsewhere.
    cwd = Path.cwd().resolve()
    if cwd != here:
        for name in SIMULATION_FILE_CANDIDATES:
            candidates.append(cwd / name)

    for path in candidates:
        if path.exists():
            spec = importlib.util.spec_from_file_location("user_simulation_module", path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            sys.modules["user_simulation_module"] = module
            spec.loader.exec_module(module)
            print(f"Loaded simulation module: {path}")
            return module

    tried = "\n".join(str(p) for p in candidates)
    raise FileNotFoundError(
        "Could not find a simulation file. Put this diagnostic script next to one of:\n"
        + "\n".join(SIMULATION_FILE_CANDIDATES)
        + "\n\nTried:\n"
        + tried
    )


# ============================================================
# 2. Shared utilities
# ============================================================


def pstr(x: float) -> str:
    return str(x).replace(".", "p").replace("-", "m")


def seed_for(rep_idx: int, offset: int = 0) -> int:
    return SEED0 + offset + rep_idx


def make_base_parameters(sim: Any) -> Any:
    return sim.Parameters(
        burn_in_deaths=BURN_IN_DEATHS,
        target_deaths=TARGET_DEATHS,
        max_proposals=MAX_PROPOSALS,
    )


def make_parameters(
    sim: Any,
    base: Any,
    *,
    seed: int,
    lambda_d: float,
    beta_d_R: float,
    beta_d_T: float,
    beta_a_R: float = 2.5,
    beta_a_T: float = 3.0,
    lambda_a_T: float = FIXED_LAMBDA_A_T,
    lambda_a_c: float = FIXED_LAMBDA_A_C,
) -> Any:
    return replace(
        base,
        seed=seed,
        lambda_a_T=lambda_a_T,
        lambda_a_c=lambda_a_c,
        beta_a_R=beta_a_R,
        beta_a_T=beta_a_T,
        lambda_d=lambda_d,
        beta_d_R=beta_d_R,
        beta_d_T=beta_d_T,
        max_proposals=MAX_PROPOSALS,
    )


def t_zone_stats_from_arrays(sim: Any, x: np.ndarray, y: np.ndarray, p: Any) -> Dict[str, float | int]:
    n = int(len(x))
    if n == 0:
        area_T = float(sim.t_zone_area(p))
        area_W = float(p.area_W)
        p0 = area_T / area_W if area_W > 0 else np.nan
        return {
            "n_observed_deaths_analyzed": 0,
            "t_zone_binomial_test_k": 0,
            "t_zone_binomial_test_n": 0,
            "t_zone_area_fraction": p0,
            "t_zone_observed_fraction": np.nan,
            "t_zone_density_ratio": np.nan,
            "t_zone_binomial_p_value_greater": np.nan,
        }

    in_T = np.array([sim.is_inside_T_zone(float(xx), float(yy), p) for xx, yy in zip(x, y)], dtype=bool)
    k = int(np.sum(in_T))
    n_out = n - k

    area_T = float(sim.t_zone_area(p))
    area_W = float(p.area_W)
    area_out = area_W - area_T
    p0 = area_T / area_W if area_W > 0 else np.nan

    frac = k / n
    density_T = k / area_T if area_T > 0 else np.nan
    density_out = n_out / area_out if area_out > 0 else np.nan
    ratio = density_T / density_out if density_out > 0 else np.nan
    pval = binomtest(k=k, n=n, p=p0, alternative="greater").pvalue if n > 0 else np.nan

    return {
        "n_observed_deaths_analyzed": n,
        "t_zone_binomial_test_k": k,
        "t_zone_binomial_test_n": n,
        "t_zone_area_fraction": p0,
        "t_zone_observed_fraction": frac,
        "t_zone_density_ratio": ratio,
        "t_zone_binomial_p_value_greater": pval,
    }


def activations_during_window(result: Any, t0: float, t1: float) -> int:
    act_t = np.asarray(result.activation_t, dtype=float)
    if act_t.size == 0:
        return 0
    return int(np.sum((act_t > t0) & (act_t <= t1)))


def row_from_fixed_death_result(
    sim: Any,
    result: Any,
    *,
    experiment: str,
    base_scenario: str,
    replicate: int,
    activation_regime: str,
    erk_regime: str,
    meaning: str = "",
) -> Dict[str, Any]:
    p = result.params
    stats = t_zone_stats_from_arrays(sim, np.asarray(result.death_x), np.asarray(result.death_y), p)
    analysis_start = float(result.burn_in_time)
    analysis_end = float(result.raw_final_time if hasattr(result, "raw_final_time") else result.final_time)
    # result has final_time; raw_final_time is only a CSV column.
    analysis_end = float(result.final_time)
    return {
        "experiment": experiment,
        "base_scenario": base_scenario,
        "replicate": replicate,
        "seed": p.seed,
        "activation_regime": activation_regime,
        "erk_regime": erk_regime,
        "regime_meaning": meaning,
        "lambda_a_T": p.lambda_a_T,
        "lambda_a_c": p.lambda_a_c,
        "beta_a_R": p.beta_a_R,
        "beta_a_T": p.beta_a_T,
        "lambda_d": p.lambda_d,
        "beta_d_R": p.beta_d_R,
        "beta_d_T": p.beta_d_T,
        "burn_in_deaths": p.burn_in_deaths,
        "target_deaths": p.target_deaths,
        "burn_in_time": float(result.burn_in_time),
        "analysis_time_span": float(result.analysis_time_span),
        "observed_death_rate_after_burnin": (
            float(result.n_observed_deaths_analyzed) / float(result.analysis_time_span)
            if float(result.analysis_time_span) > 0
            else np.nan
        ),
        "accepted_activations_total": int(result.n_accepted_activations),
        "accepted_activations_during_analysis": activations_during_window(result, analysis_start, analysis_end),
        "death_candidates": int(result.n_death_candidates),
        "death_rejected": int(result.n_death_rejected),
        "activation_candidates": int(result.n_activation_candidates),
        "activation_rejected": int(result.n_activation_rejected),
        "active_center_expirations": int(result.n_active_expirations),
        "erk_zone_expirations": int(result.n_erk_expirations),
        "all_proposals": int(result.n_proposals),
        "stopped_by_target": bool(result.stopped_by_target),
        **stats,
    }


# ============================================================
# 3. Fixed physical-time simulation
# ============================================================


def simulate_fixed_time_after_burnin(sim: Any, scenario: str, p: Any, observation_time: float) -> Dict[str, Any]:
    """Simulate until burn-in deaths, then observe all deaths for a fixed physical time."""
    rng = np.random.default_rng(p.seed)
    t = 0.0

    active_x: List[float] = []
    active_y: List[float] = []
    active_r: List[float] = []

    erk_x: List[float] = []
    erk_y: List[float] = []
    erk_r: List[float] = []

    activation_t: List[float] = []
    death_x_all: List[float] = []
    death_y_all: List[float] = []
    death_t_all: List[float] = []
    death_x_an: List[float] = []
    death_y_an: List[float] = []
    death_t_an: List[float] = []

    n_activation_candidates = 0
    n_death_candidates = 0
    n_activation_rejected = 0
    n_death_rejected = 0
    n_active_expirations = 0
    n_erk_expirations = 0
    n_proposals = 0

    burn_in_done = p.burn_in_deaths == 0
    burn_in_time = 0.0 if burn_in_done else np.nan
    t_stop = observation_time if burn_in_done else np.inf

    def is_inside_active_zone(x: float, y: float) -> bool:
        return sim.is_inside_union_of_disks(x, y, active_x, active_y, active_r)

    def is_inside_erk_zone(x: float, y: float) -> bool:
        return sim.is_inside_union_of_disks(x, y, erk_x, erk_y, erk_r)

    def activation_intensity(x: float, y: float) -> float:
        if is_inside_active_zone(x, y):
            return p.lambda_a_1
        if sim.is_inside_T_zone(x, y, p):
            return p.lambda_a_T
        return p.lambda_a_c

    while n_proposals < p.max_proposals:
        if burn_in_done and t >= t_stop:
            break

        n_active = len(active_x)
        n_erk = len(erk_x)

        activation_proposal_rate = p.lambda_a_1 * p.area_W
        death_proposal_rate = p.lambda_d * p.area_W
        active_expiration_rate = n_active * p.beta_a_T
        erk_expiration_rate = n_erk * p.beta_d_T
        a0 = activation_proposal_rate + death_proposal_rate + active_expiration_rate + erk_expiration_rate
        if a0 <= 0.0:
            break

        dt = float(rng.exponential(scale=1.0 / a0))
        if burn_in_done and t + dt > t_stop:
            t = float(t_stop)
            break

        t += dt
        n_proposals += 1

        u = rng.uniform(0.0, 1.0)
        p_activation_event = activation_proposal_rate / a0
        p_death_event = death_proposal_rate / a0
        p_active_expiration_event = active_expiration_rate / a0

        if u <= p_activation_event:
            n_activation_candidates += 1
            x = float(rng.uniform(0.0, p.Lx))
            y = float(rng.uniform(0.0, p.Ly))
            accept_prob = activation_intensity(x, y) / p.lambda_a_1
            if rng.uniform(0.0, 1.0) <= accept_prob:
                activation_t.append(t)
                r = float(rng.exponential(scale=1.0 / p.beta_a_R))
                active_x.append(x)
                active_y.append(y)
                active_r.append(r)
            else:
                n_activation_rejected += 1

        elif u <= p_activation_event + p_death_event:
            n_death_candidates += 1
            x = float(rng.uniform(0.0, p.Lx))
            y = float(rng.uniform(0.0, p.Ly))
            inside_active = is_inside_active_zone(x, y)
            inside_erk = is_inside_erk_zone(x, y)

            if inside_active and not inside_erk:
                death_x_all.append(x)
                death_y_all.append(y)
                death_t_all.append(t)

                if burn_in_done:
                    death_x_an.append(x)
                    death_y_an.append(y)
                    death_t_an.append(t)
                elif len(death_t_all) >= p.burn_in_deaths:
                    burn_in_done = True
                    burn_in_time = t
                    t_stop = burn_in_time + observation_time

                # ERK feedback after accepted observed death.
                r_E = float(rng.exponential(scale=1.0 / p.beta_d_R))
                erk_x.append(x)
                erk_y.append(y)
                erk_r.append(r_E)
            else:
                n_death_rejected += 1

        elif u <= p_activation_event + p_death_event + p_active_expiration_event:
            if active_x:
                idx = int(rng.integers(0, len(active_x)))
                active_x.pop(idx)
                active_y.pop(idx)
                active_r.pop(idx)
                n_active_expirations += 1
        else:
            if erk_x:
                idx = int(rng.integers(0, len(erk_x)))
                erk_x.pop(idx)
                erk_y.pop(idx)
                erk_r.pop(idx)
                n_erk_expirations += 1

    if not burn_in_done:
        # Simulation failed to reach burn-in before max_proposals.
        burn_in_time = float(death_t_all[-1]) if death_t_all else 0.0
        t_stop = burn_in_time

    death_x_arr = np.asarray(death_x_an, dtype=float)
    death_y_arr = np.asarray(death_y_an, dtype=float)
    death_t_arr = np.asarray(death_t_an, dtype=float)
    stats = t_zone_stats_from_arrays(sim, death_x_arr, death_y_arr, p)

    activation_t_arr = np.asarray(activation_t, dtype=float)
    activations_during_analysis = int(np.sum((activation_t_arr > burn_in_time) & (activation_t_arr <= t_stop)))

    return {
        "scenario": scenario,
        "seed": p.seed,
        "burn_in_time": float(burn_in_time),
        "observation_time_after_burnin": float(observation_time),
        "final_time": float(t),
        "analysis_time_span": float(max(0.0, t_stop - burn_in_time)),
        "accepted_activations_total": int(len(activation_t_arr)),
        "accepted_activations_during_analysis": activations_during_analysis,
        "death_candidates": int(n_death_candidates),
        "death_rejected": int(n_death_rejected),
        "activation_candidates": int(n_activation_candidates),
        "activation_rejected": int(n_activation_rejected),
        "active_center_expirations": int(n_active_expirations),
        "erk_zone_expirations": int(n_erk_expirations),
        "all_proposals": int(n_proposals),
        "stopped_by_time": bool(burn_in_done),
        **stats,
    }


# ============================================================
# 4. Aggregation utilities
# ============================================================


def grouped_summary(df: pd.DataFrame, group_cols: List[str]) -> pd.DataFrame:
    metrics = [
        "n_observed_deaths_analyzed",
        "analysis_time_span",
        "observed_death_rate_after_burnin",
        "accepted_activations_during_analysis",
        "t_zone_observed_fraction",
        "t_zone_density_ratio",
        "t_zone_binomial_p_value_greater",
        "death_candidates",
        "death_rejected",
        "activation_candidates",
        "activation_rejected",
    ]
    use_metrics = [m for m in metrics if m in df.columns]
    agg = df.groupby(group_cols, dropna=False)[use_metrics].agg(["mean", "std", "min", "max"])
    agg.columns = [f"{a}_{b}" for a, b in agg.columns]
    agg = agg.reset_index()
    return agg


def pooled_bonferroni(df: pd.DataFrame, group_cols: List[str], alpha: float = 0.05) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for keys, g in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {col: val for col, val in zip(group_cols, keys)}
        pooled_k = int(g["t_zone_binomial_test_k"].sum())
        pooled_n = int(g["t_zone_binomial_test_n"].sum())
        p0 = float(g["t_zone_area_fraction"].iloc[0])
        pooled_fraction = pooled_k / pooled_n if pooled_n > 0 else np.nan
        pval = binomtest(k=pooled_k, n=pooled_n, p=p0, alternative="greater").pvalue if pooled_n > 0 else np.nan
        n_out = pooled_n - pooled_k
        density_ratio = (pooled_k / p0) / (n_out / (1.0 - p0)) if pooled_n > 0 and n_out > 0 else np.nan
        row.update(
            {
                "n_replicates": int(len(g)),
                "pooled_t_zone_deaths_k": pooled_k,
                "pooled_total_deaths_n": pooled_n,
                "t_zone_area_fraction_p0": p0,
                "pooled_t_zone_fraction": pooled_fraction,
                "pooled_t_zone_density_ratio": density_ratio,
                "pooled_binomial_p_value_greater": pval,
            }
        )
        rows.append(row)

    out = pd.DataFrame(rows)
    m = len(out)
    bonf_alpha = alpha / m if m > 0 else np.nan
    out["bonferroni_alpha"] = bonf_alpha
    out["bonferroni_adjusted_p_value"] = np.minimum(out["pooled_binomial_p_value_greater"] * m, 1.0)
    out["significant_after_bonferroni_0p05"] = out["pooled_binomial_p_value_greater"] < bonf_alpha
    out = out.sort_values(["significant_after_bonferroni_0p05", "pooled_binomial_p_value_greater"], ascending=[False, True])
    return out


def spearman_trend_tests(df: pd.DataFrame, group_cols: List[str], x_col: str, y_cols: List[str], label: str) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for keys, g in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = {col: val for col, val in zip(group_cols, keys)}
        for y_col in y_cols:
            sub = g[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(sub) < 3 or sub[x_col].nunique() < 2 or sub[y_col].nunique() < 2:
                rho, pval = np.nan, np.nan
            else:
                rho, pval = spearmanr(sub[x_col], sub[y_col])
            rows.append(
                {
                    "experiment": label,
                    **base,
                    "x": x_col,
                    "y": y_col,
                    "n": int(len(sub)),
                    "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                    "spearman_p_value": float(pval) if np.isfinite(pval) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def write_markdown_report(out_dir: Path, sections: Dict[str, str]) -> None:
    text = "# Diagnostic experiments for T-zone mechanisms\n\n"
    for title, body in sections.items():
        text += f"## {title}\n\n{body.strip()}\n\n"
    (out_dir / "diagnostic_experiments_report.md").write_text(text, encoding="utf-8")


# ============================================================
# 5. Experiments
# ============================================================


def run_experiment_1(sim: Any, out_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame, float]:
    """Experiment 1A and 1B: lambda_d mechanism with ERK blocking regimes."""
    base = make_base_parameters(sim)
    rows_A: List[Dict[str, Any]] = []

    fixed_activation = {
        "activation_regime": "medium_activation",
        "beta_a_R": 2.5,
        "beta_a_T": 3.0,
    }

    print("\n=== Experiment 1A: fixed 1000 deaths lambda_d sweep across ERK regimes ===")
    for er in ERK_REGIMES:
        for lambda_d in LAMBDA_D_VALUES:
            base_scenario = f"E1A_{er['erk_regime']}_ld{pstr(lambda_d)}"
            for rep in range(1, N_REPLICATES + 1):
                p = make_parameters(
                    sim,
                    base,
                    seed=seed_for(rep, offset=10_000 + 1000 * ERK_REGIMES.index(er) + 100 * LAMBDA_D_VALUES.index(lambda_d)),
                    lambda_d=lambda_d,
                    beta_d_R=er["beta_d_R"],
                    beta_d_T=er["beta_d_T"],
                    beta_a_R=fixed_activation["beta_a_R"],
                    beta_a_T=fixed_activation["beta_a_T"],
                )
                scenario = f"{base_scenario}__rep{rep:02d}"
                print(scenario)
                result = sim.simulate_one(scenario, p)
                rows_A.append(
                    row_from_fixed_death_result(
                        sim,
                        result,
                        experiment="E1A_fixed_1000_deaths",
                        base_scenario=base_scenario,
                        replicate=rep,
                        activation_regime=fixed_activation["activation_regime"],
                        erk_regime=er["erk_regime"],
                        meaning=er["meaning"],
                    )
                )

    df_A = pd.DataFrame(rows_A)
    df_A.to_csv(out_dir / "experiment1A_fixed_1000_deaths_replicates.csv", index=False)
    summary_A = grouped_summary(df_A, ["erk_regime", "beta_d_R", "beta_d_T", "lambda_d"])
    summary_A.to_csv(out_dir / "experiment1A_fixed_1000_deaths_group_summary.csv", index=False)
    pooled_A = pooled_bonferroni(df_A, ["erk_regime", "beta_d_R", "beta_d_T", "lambda_d"])
    pooled_A.to_csv(out_dir / "experiment1A_fixed_1000_deaths_pooled_bonferroni.csv", index=False)

    ref_mask = (
        (df_A["erk_regime"] == "reference_ERK")
        & np.isclose(df_A["lambda_d"].astype(float), 1.0)
    )
    reference_observation_time = float(df_A.loc[ref_mask, "analysis_time_span"].mean())
    if not np.isfinite(reference_observation_time) or reference_observation_time <= 0:
        reference_observation_time = float(df_A["analysis_time_span"].median())
    print(f"Reference fixed physical observation window: {reference_observation_time:.6f}")

    print("\n=== Experiment 1B: fixed physical-time lambda_d control ===")
    rows_B: List[Dict[str, Any]] = []
    for er in ERK_REGIMES:
        for lambda_d in LAMBDA_D_VALUES:
            base_scenario = f"E1B_{er['erk_regime']}_ld{pstr(lambda_d)}"
            for rep in range(1, N_REPLICATES + 1):
                p = make_parameters(
                    sim,
                    base,
                    seed=seed_for(rep, offset=20_000 + 1000 * ERK_REGIMES.index(er) + 100 * LAMBDA_D_VALUES.index(lambda_d)),
                    lambda_d=lambda_d,
                    beta_d_R=er["beta_d_R"],
                    beta_d_T=er["beta_d_T"],
                    beta_a_R=fixed_activation["beta_a_R"],
                    beta_a_T=fixed_activation["beta_a_T"],
                )
                scenario = f"{base_scenario}__rep{rep:02d}"
                print(scenario)
                row = simulate_fixed_time_after_burnin(sim, scenario, p, reference_observation_time)
                row.update(
                    {
                        "experiment": "E1B_fixed_physical_time",
                        "base_scenario": base_scenario,
                        "replicate": rep,
                        "activation_regime": fixed_activation["activation_regime"],
                        "erk_regime": er["erk_regime"],
                        "regime_meaning": er["meaning"],
                        "lambda_a_T": p.lambda_a_T,
                        "lambda_a_c": p.lambda_a_c,
                        "beta_a_R": p.beta_a_R,
                        "beta_a_T": p.beta_a_T,
                        "lambda_d": p.lambda_d,
                        "beta_d_R": p.beta_d_R,
                        "beta_d_T": p.beta_d_T,
                    }
                )
                row["observed_death_rate_after_burnin"] = (
                    row["n_observed_deaths_analyzed"] / row["analysis_time_span"]
                    if row["analysis_time_span"] > 0
                    else np.nan
                )
                rows_B.append(row)

    df_B = pd.DataFrame(rows_B)
    df_B.to_csv(out_dir / "experiment1B_fixed_physical_time_replicates.csv", index=False)
    summary_B = grouped_summary(df_B, ["erk_regime", "beta_d_R", "beta_d_T", "lambda_d"])
    summary_B.to_csv(out_dir / "experiment1B_fixed_physical_time_group_summary.csv", index=False)
    pooled_B = pooled_bonferroni(df_B, ["erk_regime", "beta_d_R", "beta_d_T", "lambda_d"])
    pooled_B.to_csv(out_dir / "experiment1B_fixed_physical_time_pooled_bonferroni.csv", index=False)

    trend_A = spearman_trend_tests(
        df_A,
        group_cols=["erk_regime", "beta_d_R", "beta_d_T"],
        x_col="lambda_d",
        y_cols=[
            "analysis_time_span",
            "accepted_activations_during_analysis",
            "t_zone_observed_fraction",
            "t_zone_density_ratio",
        ],
        label="E1A_fixed_1000_deaths",
    )
    trend_B = spearman_trend_tests(
        df_B,
        group_cols=["erk_regime", "beta_d_R", "beta_d_T"],
        x_col="lambda_d",
        y_cols=[
            "n_observed_deaths_analyzed",
            "accepted_activations_during_analysis",
            "t_zone_observed_fraction",
            "t_zone_density_ratio",
        ],
        label="E1B_fixed_physical_time",
    )
    pd.concat([trend_A, trend_B], ignore_index=True).to_csv(out_dir / "experiment1_spearman_trend_tests.csv", index=False)

    return df_A, df_B, reference_observation_time


def run_experiment_2(sim: Any, out_dir: Path) -> pd.DataFrame:
    """Experiment 2: death/ERK robustness conditional on activation localization."""
    base = make_base_parameters(sim)
    rows: List[Dict[str, Any]] = []

    print("\n=== Experiment 2: death/ERK sweep within each activation localization regime ===")
    for ar in ACTIVATION_REGIMES:
        for er in ERK_REGIMES:
            for lambda_d in LAMBDA_D_VALUES:
                base_scenario = (
                    f"E2_{ar['activation_regime']}"
                    f"_{er['erk_regime']}"
                    f"_ld{pstr(lambda_d)}"
                )
                for rep in range(1, N_REPLICATES + 1):
                    p = make_parameters(
                        sim,
                        base,
                        seed=seed_for(
                            rep,
                            offset=(
                                30_000
                                + 5000 * ACTIVATION_REGIMES.index(ar)
                                + 1000 * ERK_REGIMES.index(er)
                                + 100 * LAMBDA_D_VALUES.index(lambda_d)
                            ),
                        ),
                        lambda_d=lambda_d,
                        beta_d_R=er["beta_d_R"],
                        beta_d_T=er["beta_d_T"],
                        beta_a_R=ar["beta_a_R"],
                        beta_a_T=ar["beta_a_T"],
                    )
                    scenario = f"{base_scenario}__rep{rep:02d}"
                    print(scenario)
                    result = sim.simulate_one(scenario, p)
                    rows.append(
                        row_from_fixed_death_result(
                            sim,
                            result,
                            experiment="E2_activation_localization_x_deathERK",
                            base_scenario=base_scenario,
                            replicate=rep,
                            activation_regime=ar["activation_regime"],
                            erk_regime=er["erk_regime"],
                            meaning=f"activation: {ar['meaning']}; ERK: {er['meaning']}",
                        )
                    )

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "experiment2_activation_localization_deathERK_replicates.csv", index=False)

    group_cols = [
        "activation_regime",
        "beta_a_R",
        "beta_a_T",
        "erk_regime",
        "beta_d_R",
        "beta_d_T",
        "lambda_d",
    ]
    summary = grouped_summary(df, group_cols)
    summary.to_csv(out_dir / "experiment2_activation_localization_deathERK_group_summary.csv", index=False)

    pooled = pooled_bonferroni(df, group_cols)
    pooled.to_csv(out_dir / "experiment2_activation_localization_deathERK_pooled_bonferroni.csv", index=False)

    # Compact regime-level robustness table.
    robust = pooled.groupby(["activation_regime", "beta_a_R", "beta_a_T"]).agg(
        n_settings=("pooled_t_zone_fraction", "size"),
        n_significant_after_bonferroni=("significant_after_bonferroni_0p05", "sum"),
        min_pooled_fraction=("pooled_t_zone_fraction", "min"),
        mean_pooled_fraction=("pooled_t_zone_fraction", "mean"),
        max_pooled_fraction=("pooled_t_zone_fraction", "max"),
        min_density_ratio=("pooled_t_zone_density_ratio", "min"),
        mean_density_ratio=("pooled_t_zone_density_ratio", "mean"),
        max_density_ratio=("pooled_t_zone_density_ratio", "max"),
    ).reset_index()
    robust.to_csv(out_dir / "experiment2_activation_regime_robustness_summary.csv", index=False)

    trend = spearman_trend_tests(
        df,
        group_cols=["activation_regime", "erk_regime", "beta_d_R", "beta_d_T"],
        x_col="lambda_d",
        y_cols=[
            "analysis_time_span",
            "accepted_activations_during_analysis",
            "t_zone_observed_fraction",
            "t_zone_density_ratio",
        ],
        label="E2_by_activation_and_ERK_regime",
    )
    trend.to_csv(out_dir / "experiment2_lambda_d_trend_tests_by_activation_and_ERK.csv", index=False)

    return df


# ============================================================
# 6. Main
# ============================================================


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sim = load_simulation_module()

    df_A, df_B, reference_time = run_experiment_1(sim, OUTPUT_DIR)
    df_E2 = run_experiment_2(sim, OUTPUT_DIR)

    sections = {
        "Experiment 1A: fixed 1000 deaths": f"""
This experiment sweeps lambda_d across three ERK regimes while keeping the activation regime fixed at
lambda_a_T={FIXED_LAMBDA_A_T}, lambda_a_c={FIXED_LAMBDA_A_C}, beta_a_R=2.5, beta_a_T=3.0.
Each run uses the original stopping rule: burn-in {BURN_IN_DEATHS} deaths, then analyze the next {TARGET_DEATHS} deaths.
Use experiment1A_fixed_1000_deaths_group_summary.csv and experiment1_spearman_trend_tests.csv to check whether lambda_d has a stable main effect across ERK regimes.
""",
        "Experiment 1B: fixed physical time": f"""
The fixed physical observation window was estimated from Experiment 1A at the reference ERK regime and lambda_d=1.0:
T_obs = {reference_time:.6f} time units after burn-in.
All lambda_d values are then compared over the same physical time window.
If lambda_d mainly changes the speed of reaching 1000 deaths, its effect on T-zone fraction should be weaker here than in Experiment 1A,
while the number of observed deaths should increase with lambda_d.
""",
        "Experiment 2: activation localization x death/ERK robustness": f"""
This experiment keeps lambda_a_T={FIXED_LAMBDA_A_T} and lambda_a_c={FIXED_LAMBDA_A_C} fixed, then compares three activation-localization regimes:
non-local (beta_a_R=1.0, beta_a_T=1.2), medium (beta_a_R=2.5, beta_a_T=3.0), and local (beta_a_R=5.0, beta_a_T=3.0).
Inside each activation regime, it sweeps the same lambda_d values and ERK regimes.
Use experiment2_activation_regime_robustness_summary.csv and experiment2_activation_localization_deathERK_pooled_bonferroni.csv to decide whether
activation localization controls the existence of T-zone enrichment, while death/ERK parameters modulate strength and speed.
""",
        "Key output files": """
- experiment1A_fixed_1000_deaths_replicates.csv
- experiment1A_fixed_1000_deaths_group_summary.csv
- experiment1A_fixed_1000_deaths_pooled_bonferroni.csv
- experiment1B_fixed_physical_time_replicates.csv
- experiment1B_fixed_physical_time_group_summary.csv
- experiment1B_fixed_physical_time_pooled_bonferroni.csv
- experiment1_spearman_trend_tests.csv
- experiment2_activation_localization_deathERK_replicates.csv
- experiment2_activation_localization_deathERK_group_summary.csv
- experiment2_activation_localization_deathERK_pooled_bonferroni.csv
- experiment2_activation_regime_robustness_summary.csv
- experiment2_lambda_d_trend_tests_by_activation_and_ERK.csv
""",
    }
    write_markdown_report(OUTPUT_DIR, sections)

    print("\nDone.")
    print(f"Outputs are in: {OUTPUT_DIR.resolve()}")
    print("Recommended first checks:")
    print("  1. experiment1_spearman_trend_tests.csv")
    print("  2. experiment1B_fixed_physical_time_group_summary.csv")
    print("  3. experiment2_activation_regime_robustness_summary.csv")
    print("  4. experiment2_activation_localization_deathERK_pooled_bonferroni.csv")


if __name__ == "__main__":
    main()
