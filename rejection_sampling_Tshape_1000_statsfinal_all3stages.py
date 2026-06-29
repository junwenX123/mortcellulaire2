"""
Rejection-sampling / thinning simulation for a T-shaped activation model
with ERK negative feedback, using the same statistical analysis workflow
as the Gillespie script.

Main points:
  1. Keep the rejection-sampling / thinning construction:
       - generate activation candidates from a dominating Poisson process;
       - generate death candidates from a homogeneous Poisson process;
       - accept/reject candidates using the active-zone and ERK-zone rules.
  2. Burn-in first: simulate 500 observed deaths but do not use them in statistics.
  3. Then collect and analyze the next 1000 observed deaths.
  4. For each parameter setting, run several independent seeds.
  5. Save the same kind of statistics as the Gillespie file:
       - number of analyzed observed deaths,
       - time needed for the analyzed 1000 deaths after burn-in,
       - mean / median inter-death time after burn-in,
       - 2D spatial histogram for visualization,
       - one-sided exact binomial test for T-zone enrichment,
       - T-zone density ratio,
       - replicate aggregation,
       - pooled Bonferroni correction across parameter settings,
       - optional diagnostic experiments matching the Gillespie diagnostic script.

"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Dict, List, Tuple
import csv
import math

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

try:
    from scipy.stats import binomtest as scipy_binomtest
except Exception:  # pragma: no cover - fallback if scipy is unavailable
    scipy_binomtest = None

try:
    from scipy.stats import spearmanr as scipy_spearmanr
except Exception:  # pragma: no cover - diagnostic trend tests become unavailable
    scipy_spearmanr = None


# ============================================================
# 1. Parameters
# ============================================================


@dataclass(frozen=True)
class Parameters:
    # Reproducibility
    seed: int = 123

    # Spatial domain W = [0, Lx] x [0, Ly]
    Lx: float = 10.0
    Ly: float = 10.0

    # Burn-in and stopping rule.
    # The simulation runs until burn_in_deaths + target_deaths observed deaths.
    # Statistics are computed only on the target_deaths after burn-in.
    burn_in_deaths: int = 500
    target_deaths: int = 1000

    # Activation intensities with fixed T-shaped extension.
    # lambda_a_1 is the dominating intensity and active-zone intensity.
    lambda_a_1: float = 5.00
    lambda_a_T: float = 0.50
    lambda_a_c: float = 0.05

    # Candidate death intensity on W.
    lambda_d: float = 1.0

    # Activation marks: R^a ~ Exp(beta_a_R), T^a ~ Exp(beta_a_T).
    beta_a_R: float = 2.5
    beta_a_T: float = 1.2

    # ERK marks after accepted deaths: R^d ~ Exp(beta_d_R), T^d ~ Exp(beta_d_T).
    # In this rejection-sampling version, T^d is pre-sampled as an exponential lifetime.
    beta_d_R: float = 2.0
    beta_d_T: float = 0.8

    # Plotting/statistics choices.
    bins_time: int = 30
    bins_space: int = 20       # used for the 2D plot and descriptive spatial CV
    bins_chi2: int = 10        # kept for backward compatibility; not used as central evidence

    # Safety stop, only to avoid infinite runs if bad parameters are chosen.
    max_proposals: int = 10_000_000

    # Rejection-sampling-specific block length.
    # Candidates are generated block by block from Poisson processes on W x [t0, t1].
    block_duration: float = 1.0
    max_blocks: int = 200_000

    @property
    def area_W(self) -> float:
        return self.Lx * self.Ly

    @property
    def total_deaths_to_simulate(self) -> int:
        return self.burn_in_deaths + self.target_deaths


@dataclass
class SimulationResult:
    scenario: str
    params: Parameters

    # All accepted observed deaths, including burn-in deaths.
    death_x_all: np.ndarray
    death_y_all: np.ndarray
    death_t_all: np.ndarray

    # Analyzed observed deaths, after burn-in.
    death_x: np.ndarray
    death_y: np.ndarray
    death_t: np.ndarray
    death_t_since_burnin: np.ndarray

    # Accepted activation centers, historical.
    activation_x: np.ndarray
    activation_y: np.ndarray
    activation_t: np.ndarray

    final_time: float
    burn_in_time: float
    analysis_time_span: float

    n_observed_deaths_total: int
    n_observed_deaths_analyzed: int
    n_burn_in_deaths: int
    n_accepted_activations: int
    n_activation_candidates: int
    n_death_candidates: int
    n_activation_rejected: int
    n_death_rejected: int
    n_death_rejected_outside_active: int
    n_death_rejected_by_erk: int
    n_proposals: int
    n_blocks: int
    stopped_by_target: bool


# ============================================================
# 2. Geometry: fixed T-shaped zone and current disk unions
# ============================================================


def validate_params(p: Parameters) -> None:
    if not (0.0 <= p.lambda_a_c <= p.lambda_a_T <= p.lambda_a_1):
        raise ValueError("Need 0 <= lambda_a_c <= lambda_a_T <= lambda_a_1 for thinning.")
    if min(p.lambda_d, p.beta_a_R, p.beta_a_T, p.beta_d_R, p.beta_d_T) <= 0:
        raise ValueError("Rates lambda_d, beta_a_R, beta_a_T, beta_d_R, beta_d_T must be positive.")
    if p.Lx <= 0 or p.Ly <= 0:
        raise ValueError("Need positive domain size.")
    if p.burn_in_deaths < 0 or p.target_deaths <= 0:
        raise ValueError("Need burn_in_deaths >= 0 and target_deaths > 0.")
    if p.bins_time <= 0 or p.bins_space <= 0 or p.bins_chi2 <= 0:
        raise ValueError("Need positive histogram bin numbers.")
    if p.block_duration <= 0 or p.max_blocks <= 0 or p.max_proposals <= 0:
        raise ValueError("Need positive block_duration, max_blocks, and max_proposals.")


def t_zone_bounds(p: Parameters) -> Tuple[float, float, float, float]:
    """Return the 3x3-grid cut points defining the fixed T-shaped zone."""
    x1 = p.Lx / 3.0
    x2 = 2.0 * p.Lx / 3.0
    y1 = p.Ly / 3.0
    y2 = 2.0 * p.Ly / 3.0
    return x1, x2, y1, y2


def t_zone_area(p: Parameters) -> float:
    """Area of the fixed T-zone: middle vertical column plus left middle arm."""
    x1, x2, y1, y2 = t_zone_bounds(p)
    middle_column_area = (x2 - x1) * p.Ly
    left_arm_area = x1 * (y2 - y1)
    return middle_column_area + left_arm_area


def is_inside_T_zone(x: float, y: float, p: Parameters) -> bool:
    """Fixed T-zone: middle vertical column plus left middle arm."""
    x1, x2, y1, y2 = t_zone_bounds(p)
    inside_middle_column = (x1 <= x <= x2) and (0.0 <= y <= p.Ly)
    inside_left_arm = (0.0 <= x <= x1) and (y1 <= y <= y2)
    return inside_middle_column or inside_left_arm


def is_inside_current_disks(x: float, y: float, centers_x: List[float], centers_y: List[float], radii: List[float]) -> bool:
    """Return True if (x,y) is inside at least one currently alive disk."""
    for cx, cy, r in zip(centers_x, centers_y, radii):
        if (x - cx) ** 2 + (y - cy) ** 2 <= r ** 2:
            return True
    return False


def prune_expired_disks(t: float, centers_x: List[float], centers_y: List[float], radii: List[float], end_times: List[float]) -> None:
    """
    Remove disks whose lifetime ended strictly before t.

    We keep disks with end_time >= t, corresponding to the convention
        S < t <= S + T.
    Since accepted disks are inserted only after their birth time, all current disks
    automatically have S < t at later candidate events with probability one.
    """
    i = len(end_times) - 1
    while i >= 0:
        if end_times[i] < t:
            centers_x.pop(i)
            centers_y.pop(i)
            radii.pop(i)
            end_times.pop(i)
        i -= 1


def add_t_zone_overlay(ax: plt.Axes, p: Parameters) -> None:
    """Add dashed rectangles showing the fixed T-zone."""
    x1, x2, y1, y2 = t_zone_bounds(p)
    ax.add_patch(Rectangle((x1, 0.0), x2 - x1, p.Ly, fill=False, linewidth=1.2, linestyle="--"))
    ax.add_patch(Rectangle((0.0, y1), x1, y2 - y1, fill=False, linewidth=1.2, linestyle="--"))


# ============================================================
# 3. Rejection sampling / thinning simulation
# ============================================================


def simulate_one(scenario: str, p: Parameters) -> SimulationResult:
    """
    Simulate until p.burn_in_deaths + p.target_deaths observed deaths are accepted.

    In each block [t0, t1], we generate:
        - activation candidates from PPP(lambda_a_1) on W x [t0, t1],
        - death candidates from PPP(lambda_d) on W x [t0, t1].

    Candidate events are processed in chronological order inside the block.

    Activation thinning probability:
        1                                      if x is inside current active zones,
        lambda_a_T / lambda_a_1               if x is in the fixed T-zone but outside active zones,
        lambda_a_c / lambda_a_1               otherwise.

    Death candidates are accepted iff they are inside active zones and outside ERK zones.
    Accepted deaths create ERK protection disks with pre-sampled exponential lifetimes.
    """
    validate_params(p)
    rng = np.random.default_rng(p.seed)

    # Current active disks V^a_t.
    active_x: List[float] = []
    active_y: List[float] = []
    active_r: List[float] = []
    active_end: List[float] = []

    # Current ERK protection disks V^p_t.
    erk_x: List[float] = []
    erk_y: List[float] = []
    erk_r: List[float] = []
    erk_end: List[float] = []

    # Historical accepted activations.
    activation_x: List[float] = []
    activation_y: List[float] = []
    activation_t: List[float] = []

    # Accepted observed deaths, including burn-in.
    death_x: List[float] = []
    death_y: List[float] = []
    death_t: List[float] = []

    n_activation_candidates = 0
    n_death_candidates = 0
    n_activation_rejected = 0
    n_death_rejected_outside_active = 0
    n_death_rejected_by_erk = 0
    n_proposals = 0

    t0 = 0.0
    n_blocks_used = 0

    for block_index in range(p.max_blocks):
        n_blocks_used = block_index + 1
        t1 = t0 + p.block_duration

        # Candidate activations in W x [t0,t1].
        n_act = int(rng.poisson(p.lambda_a_1 * p.area_W * p.block_duration))
        act_t = rng.uniform(t0, t1, size=n_act)
        act_x = rng.uniform(0.0, p.Lx, size=n_act)
        act_y = rng.uniform(0.0, p.Ly, size=n_act)
        act_u = rng.uniform(0.0, 1.0, size=n_act)
        act_r = rng.exponential(scale=1.0 / p.beta_a_R, size=n_act)
        act_tau = rng.exponential(scale=1.0 / p.beta_a_T, size=n_act)

        # Candidate deaths in W x [t0,t1].
        n_death = int(rng.poisson(p.lambda_d * p.area_W * p.block_duration))
        death_cand_t = rng.uniform(t0, t1, size=n_death)
        death_cand_x = rng.uniform(0.0, p.Lx, size=n_death)
        death_cand_y = rng.uniform(0.0, p.Ly, size=n_death)
        death_cand_r = rng.exponential(scale=1.0 / p.beta_d_R, size=n_death)
        death_cand_tau = rng.exponential(scale=1.0 / p.beta_d_T, size=n_death)

        n_activation_candidates += n_act
        n_death_candidates += n_death
        n_proposals += n_act + n_death

        if n_proposals > p.max_proposals:
            raise RuntimeError(
                f"Safety stop: processed {n_proposals} candidates before reaching "
                f"{p.total_deaths_to_simulate} observed deaths. Increase max_proposals."
            )

        # Combine event times and process the two independent candidate PPPs chronologically.
        # type 0 = activation candidate, type 1 = death candidate.
        event_t = np.concatenate([act_t, death_cand_t])
        event_type = np.concatenate([np.zeros(n_act, dtype=np.int8), np.ones(n_death, dtype=np.int8)])
        event_index = np.concatenate([np.arange(n_act, dtype=np.int64), np.arange(n_death, dtype=np.int64)])

        for pos in np.argsort(event_t):
            if len(death_t) >= p.total_deaths_to_simulate:
                break

            t = float(event_t[pos])
            kind = int(event_type[pos])
            idx = int(event_index[pos])
            
            prune_expired_disks(t, active_x, active_y, active_r, active_end)
            prune_expired_disks(t, erk_x, erk_y, erk_r, erk_end)

            if kind == 0:
                x = float(act_x[idx])
                y = float(act_y[idx])
                u = float(act_u[idx])
                r = float(act_r[idx])
                tau = float(act_tau[idx])

                inside_active = is_inside_current_disks(x, y, active_x, active_y, active_r)
                if inside_active:
                    accept_prob = 1.0
                elif is_inside_T_zone(x, y, p):
                    accept_prob = p.lambda_a_T / p.lambda_a_1
                else:
                    accept_prob = p.lambda_a_c / p.lambda_a_1

                if u <= accept_prob:
                    activation_x.append(x)
                    activation_y.append(y)
                    activation_t.append(t)
                    active_x.append(x)
                    active_y.append(y)
                    active_r.append(r)
                    active_end.append(t + tau)
                else:
                    n_activation_rejected += 1

            else:
                x = float(death_cand_x[idx])
                y = float(death_cand_y[idx])
                r = float(death_cand_r[idx])
                tau = float(death_cand_tau[idx])

                inside_active = is_inside_current_disks(x, y, active_x, active_y, active_r)
                inside_erk = is_inside_current_disks(x, y, erk_x, erk_y, erk_r)

                if inside_active and not inside_erk:
                    death_x.append(x)
                    death_y.append(y)
                    death_t.append(t)
                    erk_x.append(x)
                    erk_y.append(y)
                    erk_r.append(r)
                    erk_end.append(t + tau)
                else:
                    if inside_erk:
                        n_death_rejected_by_erk += 1
                    else:
                        n_death_rejected_outside_active += 1

        if len(death_t) >= p.total_deaths_to_simulate:
            break

        t0 = t1
    else:
        raise RuntimeError(
            f"Only {len(death_t)} observed deaths after {p.max_blocks} blocks. "
            "Increase max_blocks, max_proposals, lambda_d, lambda_a_T/lambda_a_c, "
            "or active-zone coverage."
        )

    death_x_all = np.asarray(death_x[:p.total_deaths_to_simulate])
    death_y_all = np.asarray(death_y[:p.total_deaths_to_simulate])
    death_t_all = np.asarray(death_t[:p.total_deaths_to_simulate])

    if len(death_t_all) > p.burn_in_deaths:
        burn_in_time = 0.0 if p.burn_in_deaths == 0 else float(death_t_all[p.burn_in_deaths - 1])
        death_x_an = death_x_all[p.burn_in_deaths:]
        death_y_an = death_y_all[p.burn_in_deaths:]
        death_t_an = death_t_all[p.burn_in_deaths:]
        death_t_since_burnin = death_t_an - burn_in_time
        analysis_time_span = float(death_t_all[-1] - burn_in_time) if len(death_t_an) else 0.0
    else:
        burn_in_time = float(death_t_all[-1]) if len(death_t_all) else 0.0
        death_x_an = np.array([])
        death_y_an = np.array([])
        death_t_an = np.array([])
        death_t_since_burnin = np.array([])
        analysis_time_span = 0.0

    n_death_rejected = n_death_rejected_outside_active + n_death_rejected_by_erk

    return SimulationResult(
        scenario=scenario,
        params=p,
        death_x_all=death_x_all,
        death_y_all=death_y_all,
        death_t_all=death_t_all,
        death_x=death_x_an,
        death_y=death_y_an,
        death_t=death_t_an,
        death_t_since_burnin=death_t_since_burnin,
        activation_x=np.asarray(activation_x),
        activation_y=np.asarray(activation_y),
        activation_t=np.asarray(activation_t),
        final_time=float(death_t_all[-1]) if len(death_t_all) else 0.0,
        burn_in_time=burn_in_time,
        analysis_time_span=analysis_time_span,
        n_observed_deaths_total=len(death_t_all),
        n_observed_deaths_analyzed=len(death_t_an),
        n_burn_in_deaths=min(p.burn_in_deaths, len(death_t_all)),
        n_accepted_activations=len(activation_t),
        n_activation_candidates=n_activation_candidates,
        n_death_candidates=n_death_candidates,
        n_activation_rejected=n_activation_rejected,
        n_death_rejected=n_death_rejected,
        n_death_rejected_outside_active=n_death_rejected_outside_active,
        n_death_rejected_by_erk=n_death_rejected_by_erk,
        n_proposals=n_proposals,
        n_blocks=n_blocks_used,
        stopped_by_target=(len(death_t_all) >= p.total_deaths_to_simulate),
    )


# ============================================================
# 4. Statistics and plots on the observed death process
# ============================================================


def spatial_histogram(result: SimulationResult, bins: int) -> np.ndarray:
    p = result.params
    H, _, _ = np.histogram2d(
        result.death_x,
        result.death_y,
        bins=bins,
        range=[[0.0, p.Lx], [0.0, p.Ly]],
    )
    return H



# ============================================================
# 4. Statistics and plots on the observed death process
#    same central evidence as the Gillespie scripts
# ============================================================


def binomial_upper_tail_p_value(k: int, n: int, p0: float) -> float:
    """
    One-sided exact binomial p-value for enrichment.

    H0: pi = p0, HA: pi > p0, with X ~ Binomial(n, p0).
    The p-value is P(X >= k).
    """
    if n <= 0 or k < 0 or k > n or not np.isfinite(p0) or p0 < 0.0 or p0 > 1.0:
        return np.nan

    if scipy_binomtest is not None:
        return float(scipy_binomtest(k=k, n=n, p=p0, alternative="greater").pvalue)

    # Stable fallback without scipy, useful for n around 1000 or 5000.
    if p0 == 0.0:
        return 1.0 if k == 0 else 0.0
    if p0 == 1.0:
        return 1.0
    if k <= 0:
        return 1.0

    log_terms: List[float] = []
    log_p = math.log(p0)
    log_1mp = math.log1p(-p0)
    for i in range(k, n + 1):
        log_choose = math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
        log_terms.append(log_choose + i * log_p + (n - i) * log_1mp)
    m = max(log_terms)
    return float(math.exp(m) * sum(math.exp(v - m) for v in log_terms))


def t_zone_density_stats(result: SimulationResult) -> Dict[str, float | int | bool]:
    """T-zone density ratio and one-sided exact binomial enrichment test."""
    p = result.params
    in_T = np.array([
        is_inside_T_zone(float(x), float(y), p)
        for x, y in zip(result.death_x, result.death_y)
    ], dtype=bool)
    n = int(len(in_T))
    n_T = int(np.sum(in_T))
    n_out = int(n - n_T)

    area_T = float(t_zone_area(p))
    area_out = float(p.area_W - area_T)
    p0 = area_T / float(p.area_W) if p.area_W > 0 else np.nan

    observed_fraction = n_T / n if n > 0 else np.nan
    density_T = n_T / area_T if area_T > 0 else np.nan
    density_out = n_out / area_out if area_out > 0 else np.nan
    ratio = density_T / density_out if density_out > 0 else np.nan
    p_value_greater = binomial_upper_tail_p_value(k=n_T, n=n, p0=p0)

    return {
        "t_zone_area": area_T,
        "outside_t_zone_area": area_out,
        "t_zone_observed_deaths": n_T,
        "outside_t_zone_observed_deaths": n_out,
        "t_zone_observed_fraction": observed_fraction,
        "t_zone_area_fraction": p0,
        "t_zone_death_density": density_T,
        "outside_t_zone_death_density": density_out,
        "t_zone_density_ratio": ratio,
        "t_zone_binomial_test_k": n_T,
        "t_zone_binomial_test_n": n,
        "t_zone_binomial_null_probability_area_fraction": p0,
        "t_zone_binomial_observed_fraction": observed_fraction,
        "t_zone_binomial_p_value_greater": p_value_greater,
        "t_zone_binomial_significant_0p05": bool(np.isfinite(p_value_greater) and p_value_greater < 0.05),
    }


def split_scenario_replicate(scenario: str) -> Tuple[str, int]:
    """Split names such as 'baseline__rep03' into ('baseline', 3)."""
    if "__rep" not in scenario:
        return scenario, 1
    base_name, rep_part = scenario.rsplit("__rep", 1)
    try:
        return base_name, int(rep_part)
    except ValueError:
        return base_name, 1


def death_statistics(result: SimulationResult) -> Dict[str, float | int | str | bool]:
    """Replicate-level summary, with T-zone binomial evidence as the central statistic."""
    p = result.params
    base_scenario, replicate = split_scenario_replicate(result.scenario)

    if len(result.death_t) > 0:
        inter_death_after_burnin = np.diff(np.concatenate([[result.burn_in_time], result.death_t]))
    else:
        inter_death_after_burnin = np.array([])

    activation_acceptance_rate = (
        result.n_accepted_activations / result.n_activation_candidates
        if result.n_activation_candidates > 0 else np.nan
    )
    death_rejection_rate = (
        result.n_death_rejected / result.n_death_candidates
        if result.n_death_candidates > 0 else np.nan
    )
    death_rejection_rate_outside_active = (
        result.n_death_rejected_outside_active / result.n_death_candidates
        if result.n_death_candidates > 0 else np.nan
    )
    death_rejection_rate_by_erk = (
        result.n_death_rejected_by_erk / result.n_death_candidates
        if result.n_death_candidates > 0 else np.nan
    )
    death_acceptance_rate_total = (
        result.n_observed_deaths_total / result.n_death_candidates
        if result.n_death_candidates > 0 else np.nan
    )

    t_stats = t_zone_density_stats(result)

    return {
        "scenario": result.scenario,
        "base_scenario": base_scenario,
        "replicate": replicate,
        "seed": p.seed,
        "burn_in_deaths": p.burn_in_deaths,
        "analyzed_target_deaths": p.target_deaths,
        "total_deaths_simulated": result.n_observed_deaths_total,
        "n_observed_deaths_analyzed": result.n_observed_deaths_analyzed,
        "stopped_by_target": result.stopped_by_target,
        "burn_in_time": result.burn_in_time,
        "raw_final_time": result.final_time,
        "total_time_for_analyzed_1000_deaths_after_burnin": result.analysis_time_span,
        "observed_death_rate_after_burnin": result.n_observed_deaths_analyzed / result.analysis_time_span if result.analysis_time_span > 0 else np.nan,
        "mean_inter_death_time_after_burnin": float(np.mean(inter_death_after_burnin)) if len(inter_death_after_burnin) else np.nan,
        "median_inter_death_time_after_burnin": float(np.median(inter_death_after_burnin)) if len(inter_death_after_burnin) else np.nan,
        "accepted_activations_total": result.n_accepted_activations,
        "activation_candidates": result.n_activation_candidates,
        "activation_rejected": result.n_activation_rejected,
        "activation_acceptance_rate": activation_acceptance_rate,
        "death_candidates": result.n_death_candidates,
        "death_rejected": result.n_death_rejected,
        "death_rejected_outside_active": result.n_death_rejected_outside_active,
        "death_rejected_by_erk": result.n_death_rejected_by_erk,
        "death_rejection_rate": death_rejection_rate,
        "death_rejection_rate_outside_active": death_rejection_rate_outside_active,
        "death_rejection_rate_by_erk": death_rejection_rate_by_erk,
        "death_acceptance_rate_total": death_acceptance_rate_total,
        "all_proposals": result.n_proposals,
        "n_blocks": result.n_blocks,
        "spatial_bins_per_axis_for_plot": p.bins_space,
        **t_stats,
        "lambda_a_1": p.lambda_a_1,
        "lambda_a_T": p.lambda_a_T,
        "lambda_a_c": p.lambda_a_c,
        "lambda_a_T_over_lambda_a_c": p.lambda_a_T / p.lambda_a_c if p.lambda_a_c > 0 else np.inf,
        "lambda_d": p.lambda_d,
        "beta_a_R": p.beta_a_R,
        "beta_a_T": p.beta_a_T,
        "beta_d_R": p.beta_d_R,
        "beta_d_T": p.beta_d_T,
    }


# ============================================================
# 5. Output helpers
# ============================================================


def save_death_events_csv(result: SimulationResult, out_dir: Path) -> Path:
    path = out_dir / f"{result.scenario}_observed_deaths_after_burnin.csv"
    p = result.params
    in_T = [is_inside_T_zone(float(x), float(y), p) for x, y in zip(result.death_x, result.death_y)]
    raw_indices = np.arange(p.burn_in_deaths + 1, p.burn_in_deaths + result.n_observed_deaths_analyzed + 1)
    analysis_indices = np.arange(1, result.n_observed_deaths_analyzed + 1)
    data = np.column_stack([
        analysis_indices,
        raw_indices,
        result.death_t,
        result.death_t_since_burnin,
        result.death_x,
        result.death_y,
        np.asarray(in_T, dtype=int),
    ])
    header = "analysis_death_index,raw_death_index,raw_death_time,time_since_burnin,x,y,in_T_zone"
    np.savetxt(
        path,
        data,
        delimiter=",",
        header=header,
        comments="",
        fmt=["%d", "%d", "%.10f", "%.10f", "%.10f", "%.10f", "%d"],
    )
    return path


def wrap_underscore_name(name: str, width: int = 55) -> str:
    """Wrap a long scenario name at underscores so that figure titles are not clipped."""
    base_name, replicate = split_scenario_replicate(name)
    parts = base_name.split("_")
    lines: List[str] = []
    current = ""
    for part in parts:
        candidate = part if current == "" else current + "_" + part
        if len(candidate) > width and current != "":
            lines.append(current)
            current = part
        else:
            current = candidate
    if current:
        lines.append(current)
    if "__rep" in name:
        lines.append(f"rep{replicate:02d}")
    return "\n".join(lines)


def plot_time_histogram(result: SimulationResult, out_dir: Path) -> Path:
    p = result.params
    path = out_dir / f"{result.scenario}_time_histogram_after_burnin.png"
    scenario_title = wrap_underscore_name(result.scenario, width=55)

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    ax.hist(result.death_t_since_burnin, bins=p.bins_time, range=(0.0, result.analysis_time_span))
    ax.set_xlabel("time after burn-in")
    ax.set_ylabel("number of observed deaths")
    ax.set_title(
        f"Observed cell death times after burn-in\n"
        f"{scenario_title}\n"
        f"burn-in = {p.burn_in_deaths}, N = {result.n_observed_deaths_analyzed}, "
        f"time span = {result.analysis_time_span:.3f}",
        fontsize=11,
        pad=12,
    )
    fig.savefig(path, dpi=180, bbox_inches="tight", pad_inches=0.20)
    plt.close(fig)
    return path


def plot_spatial_histogram_2d(result: SimulationResult, out_dir: Path, colorbar_vmax: int) -> Path:
    p = result.params
    path = out_dir / f"{result.scenario}_space_histogram_2d_after_burnin.png"
    H = spatial_histogram(result, p.bins_space)
    scenario_title = wrap_underscore_name(result.scenario, width=55)

    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    image = ax.imshow(
        H.T,
        origin="lower",
        extent=[0.0, p.Lx, 0.0, p.Ly],
        aspect="equal",
        interpolation="nearest",
        vmin=0,
        vmax=colorbar_vmax,
    )
    add_t_zone_overlay(ax, p)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        f"2D histogram of observed death locations\n"
        f"{scenario_title}\n"
        f"after burn-in {p.burn_in_deaths}, {p.bins_space} x {p.bins_space} bins, "
        f"N = {result.n_observed_deaths_analyzed}",
        fontsize=11,
        pad=12,
    )
    cbar = fig.colorbar(image, ax=ax, label="number of observed deaths", shrink=0.88)
    cbar.set_ticks(np.arange(0, colorbar_vmax + 1, max(1, colorbar_vmax // 7)))
    fig.savefig(path, dpi=180, bbox_inches="tight", pad_inches=0.20)
    plt.close(fig)
    return path


def plot_cumulative_deaths(result: SimulationResult, out_dir: Path) -> Path:
    path = out_dir / f"{result.scenario}_cumulative_deaths_after_burnin.png"
    scenario_title = wrap_underscore_name(result.scenario, width=65)

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    ax.step(result.death_t_since_burnin, np.arange(1, result.n_observed_deaths_analyzed + 1), where="post")
    ax.set_xlabel("time after burn-in")
    ax.set_ylabel("cumulative observed deaths after burn-in")
    ax.set_title(f"Cumulative observed deaths after burn-in\n{scenario_title}", fontsize=11, pad=12)
    fig.savefig(path, dpi=180, bbox_inches="tight", pad_inches=0.20)
    plt.close(fig)
    return path


def save_summary_csv(rows: List[Dict[str, Any]], out_dir: Path) -> Path:
    path = out_dir / "parameter_sweep_summary_after_burnin.csv"
    keys = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return np.nan


def _mean_std_min_max(values: List[float]) -> Tuple[float, float, float, float]:
    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return np.nan, np.nan, np.nan, np.nan
    mean = float(np.mean(arr))
    std = float(np.std(arr, ddof=1)) if arr.size >= 2 else 0.0
    return mean, std, float(np.min(arr)), float(np.max(arr))


def aggregate_replicate_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Aggregate replicate-level rows by base_scenario."""
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        base_name = str(row.get("base_scenario", row.get("scenario", "unknown")))
        groups.setdefault(base_name, []).append(row)

    metrics = [
        "total_time_for_analyzed_1000_deaths_after_burnin",
        "observed_death_rate_after_burnin",
        "t_zone_observed_deaths",
        "t_zone_observed_fraction",
        "t_zone_density_ratio",
        "t_zone_binomial_p_value_greater",
        "death_rejected",
        "death_rejection_rate",
        "death_rejection_rate_outside_active",
        "death_rejection_rate_by_erk",
        "accepted_activations_total",
        "activation_acceptance_rate",
    ]

    aggregated_rows: List[Dict[str, Any]] = []
    for base_name, group in groups.items():
        first = group[0]
        out: Dict[str, Any] = {
            "base_scenario": base_name,
            "n_replicates": len(group),
            "n_finished_target": sum(bool(r.get("stopped_by_target", False)) for r in group),
            "lambda_a_1": first.get("lambda_a_1", np.nan),
            "lambda_a_T": first.get("lambda_a_T", np.nan),
            "lambda_a_c": first.get("lambda_a_c", np.nan),
            "lambda_a_T_over_lambda_a_c": first.get("lambda_a_T_over_lambda_a_c", np.nan),
            "lambda_d": first.get("lambda_d", np.nan),
            "beta_a_R": first.get("beta_a_R", np.nan),
            "beta_a_T": first.get("beta_a_T", np.nan),
            "beta_d_R": first.get("beta_d_R", np.nan),
            "beta_d_T": first.get("beta_d_T", np.nan),
        }

        p_values = [_safe_float(r.get("t_zone_binomial_p_value_greater", np.nan)) for r in group]
        finite_p_values = [pv for pv in p_values if np.isfinite(pv)]
        out["fraction_t_zone_binomial_p_value_below_0p05"] = (
            float(np.mean([pv < 0.05 for pv in finite_p_values])) if finite_p_values else np.nan
        )

        for metric in metrics:
            values = [_safe_float(r.get(metric, np.nan)) for r in group]
            mean, std, min_v, max_v = _mean_std_min_max(values)
            out[f"{metric}_mean"] = mean
            out[f"{metric}_std"] = std
            out[f"{metric}_min"] = min_v
            out[f"{metric}_max"] = max_v

        aggregated_rows.append(out)

    return aggregated_rows


def save_aggregate_summary_csv(rows: List[Dict[str, Any]], out_dir: Path) -> Path:
    path = out_dir / "parameter_sweep_summary_by_parameter_setting_after_burnin.csv"
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    keys = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    return path


def pooled_bonferroni_from_rows(rows: List[Dict[str, Any]], alpha: float = 0.05) -> List[Dict[str, Any]]:
    """Pool independent replicates by parameter setting and apply Bonferroni correction."""
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(str(row["base_scenario"]), []).append(row)

    pooled_rows: List[Dict[str, Any]] = []
    for base_scenario, group in groups.items():
        first = group[0]
        pooled_k = int(sum(int(r["t_zone_binomial_test_k"]) for r in group))
        pooled_n = int(sum(int(r["t_zone_binomial_test_n"]) for r in group))
        p0 = float(first["t_zone_binomial_null_probability_area_fraction"])
        pooled_fraction = pooled_k / pooled_n if pooled_n > 0 else np.nan
        pooled_p = binomial_upper_tail_p_value(pooled_k, pooled_n, p0)
        n_out = pooled_n - pooled_k
        pooled_ratio = (pooled_k / p0) / (n_out / (1.0 - p0)) if pooled_n > 0 and n_out > 0 and 0 < p0 < 1 else np.nan

        pooled_rows.append({
            "base_scenario": base_scenario,
            "n_replicates": int(len(group)),
            "pooled_t_zone_deaths_k": pooled_k,
            "pooled_total_deaths_n": pooled_n,
            "t_zone_area_fraction_p0": p0,
            "pooled_t_zone_fraction": pooled_fraction,
            "pooled_t_zone_density_ratio": pooled_ratio,
            "pooled_binomial_p_value_greater": pooled_p,
            "lambda_a_1": first.get("lambda_a_1", np.nan),
            "lambda_a_T": first.get("lambda_a_T", np.nan),
            "lambda_a_c": first.get("lambda_a_c", np.nan),
            "lambda_d": first.get("lambda_d", np.nan),
            "beta_a_R": first.get("beta_a_R", np.nan),
            "beta_a_T": first.get("beta_a_T", np.nan),
            "beta_d_R": first.get("beta_d_R", np.nan),
            "beta_d_T": first.get("beta_d_T", np.nan),
        })

    m = len(pooled_rows)
    bonf_alpha = alpha / m if m > 0 else np.nan
    for row in pooled_rows:
        pval = float(row["pooled_binomial_p_value_greater"])
        row["bonferroni_alpha"] = bonf_alpha
        row["bonferroni_adjusted_p_value"] = min(pval * m, 1.0) if np.isfinite(pval) else np.nan
        row["significant_after_bonferroni_0p05"] = bool(np.isfinite(pval) and pval < bonf_alpha)

    pooled_rows.sort(key=lambda r: (not bool(r["significant_after_bonferroni_0p05"]), _safe_float(r["pooled_binomial_p_value_greater"])))
    return pooled_rows


def save_pooled_bonferroni_outputs(rows: List[Dict[str, Any]], out_dir: Path, alpha: float = 0.05) -> Tuple[Path, Path, List[Dict[str, Any]]]:
    pooled_rows = pooled_bonferroni_from_rows(rows, alpha=alpha)
    csv_path = out_dir / "pooled_bonferroni_summary.csv"
    md_path = out_dir / "pooled_bonferroni_report.md"

    if pooled_rows:
        keys = list(pooled_rows[0].keys())
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(pooled_rows)
    else:
        csv_path.write_text("", encoding="utf-8")

    m = len(pooled_rows)
    bonf_alpha = alpha / m if m > 0 else np.nan
    lines: List[str] = []
    lines.append("# Pooled T-zone binomial test with Bonferroni correction: rejection-sampling version\n\n")
    lines.append("## Method\n\n")
    lines.append(
        "For each parameter setting, independent replicates are pooled. If a parameter setting has five replicates "
        "of 1000 analyzed observed deaths, then the pooled sample size is N = 5000. The pooled T-zone count is tested "
        "against p0 = |T|/|W| by a one-sided exact binomial test.\n\n"
    )
    lines.append("The hypotheses are:\n\n")
    lines.append("```latex\nH_0:\\pi=p_0,\\qquad H_A:\\pi>p_0.\n```\n\n")
    lines.append(
        f"Since {m} parameter settings are tested, the Bonferroni threshold at global alpha={alpha} is "
        f"alpha_Bonferroni = {bonf_alpha:.6g}.\n\n"
    )
    lines.append("## Results\n\n")
    lines.append(
        "| parameter setting | replicates | pooled k/n | pooled fraction | pooled density ratio | pooled p-value | Bonferroni adjusted p-value | significant |\n"
        "|---|---:|---:|---:|---:|---:|---:|---:|\n"
    )
    for r in pooled_rows:
        lines.append(
            f"| {r['base_scenario']} | "
            f"{int(r['n_replicates'])} | "
            f"{int(r['pooled_t_zone_deaths_k'])}/{int(r['pooled_total_deaths_n'])} | "
            f"{float(r['pooled_t_zone_fraction']):.4f} | "
            f"{float(r['pooled_t_zone_density_ratio']):.4f} | "
            f"{float(r['pooled_binomial_p_value_greater']):.3g} | "
            f"{float(r['bonferroni_adjusted_p_value']):.3g} | "
            f"{bool(r['significant_after_bonferroni_0p05'])} |\n"
        )
    md_path.write_text("".join(lines), encoding="utf-8")
    return csv_path, md_path, pooled_rows


def write_analysis_report(
    rows: List[Dict[str, Any]],
    out_dir: Path,
    colorbar_vmax: int,
    aggregated_rows: List[Dict[str, Any]] | None = None,
    pooled_rows: List[Dict[str, Any]] | None = None,
) -> Path:
    path = out_dir / "t_zone_binomial_analysis_report.md"
    if not rows:
        path.write_text("No results.\n", encoding="utf-8")
        return path

    first = rows[0]
    lines: List[str] = []
    lines.append("# T-zone binomial enrichment analysis after burn-in: rejection-sampling version\n\n")
    lines.append(f"All 2D histograms use the same colorbar scale: 0 to {colorbar_vmax} observed deaths per bin.\n\n")
    lines.append(
        f"Burn-in removes the initial transient period: the first {int(first['burn_in_deaths'])} observed deaths "
        f"are simulated but not used for statistics. The reported statistics use the next "
        f"{int(first['n_observed_deaths_analyzed'])} observed deaths.\n\n"
    )
    lines.append("## Interpretation rules\n\n")
    lines.append("- `t_zone_density_ratio`: values > 1 mean deaths are denser in the T-zone than outside.\n")
    lines.append("- `t_zone_binomial_p_value_greater`: one-sided exact binomial p-value for `H0: pi = |T|/|W|` against `HA: pi > |T|/|W|`.\n")
    lines.append("- `pooled_binomial_p_value_greater`: same test after pooling independent replicates of the same parameter setting.\n")
    lines.append("- Old `spatial CV` and `chi-square uniformity` are intentionally not used as the central evidence, because they test global non-uniformity rather than enrichment in the fixed T-zone.\n\n")

    lines.append("## Replicate-level results\n\n")
    lines.append(
        "| scenario | time for analyzed 1000 deaths | T-zone deaths k/n | area fraction p0 | observed fraction | binomial p-value greater | T-zone density ratio | death rejection | by ERK |\n"
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|\n"
    )
    for r in rows:
        lines.append(
            f"| {r['scenario']} | "
            f"{float(r['total_time_for_analyzed_1000_deaths_after_burnin']):.3f} | "
            f"{int(r['t_zone_binomial_test_k'])}/{int(r['t_zone_binomial_test_n'])} | "
            f"{float(r['t_zone_binomial_null_probability_area_fraction']):.4f} | "
            f"{float(r['t_zone_binomial_observed_fraction']):.4f} | "
            f"{float(r['t_zone_binomial_p_value_greater']):.3g} | "
            f"{float(r['t_zone_density_ratio']):.4f} | "
            f"{float(r['death_rejection_rate']):.4f} | "
            f"{float(r['death_rejection_rate_by_erk']):.4f} |\n"
        )

    if aggregated_rows:
        lines.append("\n## Aggregated replicate summary\n\n")
        lines.append(
            "| parameter setting | replicates | mean T-zone fraction | mean T-zone ratio | mean binomial p-value | fraction p<0.05 | mean time | mean death rejection |\n"
            "|---|---:|---:|---:|---:|---:|---:|---:|\n"
        )
        for r in aggregated_rows:
            lines.append(
                f"| {r['base_scenario']} | "
                f"{int(r['n_replicates'])} | "
                f"{float(r['t_zone_observed_fraction_mean']):.4f} | "
                f"{float(r['t_zone_density_ratio_mean']):.4f} | "
                f"{float(r['t_zone_binomial_p_value_greater_mean']):.3g} | "
                f"{float(r['fraction_t_zone_binomial_p_value_below_0p05']):.3f} | "
                f"{float(r['total_time_for_analyzed_1000_deaths_after_burnin_mean']):.3f} | "
                f"{float(r['death_rejection_rate_mean']):.4f} |\n"
            )

    if pooled_rows:
        lines.append("\n## Pooled Bonferroni summary\n\n")
        lines.append(
            "| parameter setting | pooled k/n | pooled fraction | pooled ratio | pooled p-value | Bonferroni alpha | significant |\n"
            "|---|---:|---:|---:|---:|---:|---:|\n"
        )
        for r in pooled_rows:
            lines.append(
                f"| {r['base_scenario']} | "
                f"{int(r['pooled_t_zone_deaths_k'])}/{int(r['pooled_total_deaths_n'])} | "
                f"{float(r['pooled_t_zone_fraction']):.4f} | "
                f"{float(r['pooled_t_zone_density_ratio']):.4f} | "
                f"{float(r['pooled_binomial_p_value_greater']):.3g} | "
                f"{float(r['bonferroni_alpha']):.6g} | "
                f"{bool(r['significant_after_bonferroni_0p05'])} |\n"
            )

    baseline = rows[0]
    strongest_by_ratio = max(rows, key=lambda r: _safe_float(r["t_zone_density_ratio"]))
    strongest_by_pvalue = min(rows, key=lambda r: _safe_float(r["t_zone_binomial_p_value_greater"]))
    lines.append("\n## Short conclusion\n\n")
    lines.append(
        f"The first setting `{baseline['scenario']}` has observed T-zone fraction "
        f"{float(baseline['t_zone_binomial_observed_fraction']):.4f}, area fraction "
        f"{float(baseline['t_zone_binomial_null_probability_area_fraction']):.4f}, "
        f"one-sided binomial p-value {float(baseline['t_zone_binomial_p_value_greater']):.3g}, "
        f"and density ratio {float(baseline['t_zone_density_ratio']):.4f}.\n"
    )
    lines.append(
        f"The largest replicate-level T-zone density ratio occurs in `{strongest_by_ratio['scenario']}` "
        f"with ratio {float(strongest_by_ratio['t_zone_density_ratio']):.4f}.\n"
    )
    lines.append(
        f"The smallest replicate-level one-sided binomial p-value occurs in `{strongest_by_pvalue['scenario']}` "
        f"with p-value {float(strongest_by_pvalue['t_zone_binomial_p_value_greater']):.3g}.\n"
    )

    path.write_text("".join(lines), encoding="utf-8")
    return path


# ============================================================
# 6. Parameter sweep: same parameter families as the Gillespie analysis
# ============================================================


def make_scenarios(base: Parameters) -> Dict[str, Parameters]:
    """
    Build a replicated parameter sweep.

    The statistics reported later are aligned with
    the Gillespie files.
    """
    n_replicates = 5
    replicate_seeds = [1001 + i for i in range(n_replicates)]

    parameter_sets: Dict[str, Parameters] = {
        "baseline": base,
        "higher_death_rate": replace(base, lambda_d=1.5),
        "larger_ERK_radius": replace(base, beta_d_R=1.0),
        "shorter_ERK_duration": replace(base, beta_d_T=2.0),
        "stronger_T_zone_activation": replace(base, lambda_a_T=1.0),
        "strong_visible_T_zone": replace(
            base,
            lambda_a_T=2.0,
            lambda_a_c=0.005,
            beta_a_R=5.0,
            beta_a_T=3.0,
            max_proposals=50_000_000,
        ),
        "T0p5_c0p005_ratio_1e2": replace(base, lambda_a_T=5e-1, lambda_a_c=5e-3),
        "T0p5_c0p0005_ratio_1e3": replace(base, lambda_a_T=5e-1, lambda_a_c=5e-4),
        "T0p5_c0p00005_ratio_1e4": replace(base, lambda_a_T=5e-1, lambda_a_c=5e-5),
        "T0p5_c0p025_ratio_20": replace(base, lambda_a_T=5e-1, lambda_a_c=2.5e-2),
        "T0p5_c0p01_ratio_50": replace(base, lambda_a_T=5e-1, lambda_a_c=1e-2),
        "T1e-1_c1e-2_ratio_1e1": replace(base, lambda_a_T=1e-1, lambda_a_c=1e-2),
        "T1e-1_c1e-3_ratio_1e2": replace(base, lambda_a_T=1e-1, lambda_a_c=1e-3),
        "T1e-1_c1e-4_ratio_1e3": replace(base, lambda_a_T=1e-1, lambda_a_c=1e-4),
        "T1e-2_c1e-3_ratio_1e1": replace(base, lambda_a_T=1e-2, lambda_a_c=1e-3, max_proposals=30_000_000),
        "T1e-2_c1e-4_ratio_1e2": replace(base, lambda_a_T=1e-2, lambda_a_c=1e-4, max_proposals=30_000_000),
        "T1e-3_c1e-4_ratio_1e1": replace(base, lambda_a_T=1e-3, lambda_a_c=1e-4, max_proposals=50_000_000),
        "T0p5_c0p005_ratio_1e2_local_activation": replace(
            base, lambda_a_T=5e-1, lambda_a_c=5e-3, beta_a_R=5.0, beta_a_T=3.0, max_proposals=50_000_000
        ),
        "T0p5_c0p0005_ratio_1e3_local_activation": replace(
            base, lambda_a_T=5e-1, lambda_a_c=5e-4, beta_a_R=5.0, beta_a_T=3.0, max_proposals=50_000_000
        ),
        "T1e-1_c1e-4_ratio_1e3_local_activation": replace(
            base, lambda_a_T=1e-1, lambda_a_c=1e-4, beta_a_R=5.0, beta_a_T=3.0, max_proposals=80_000_000
        ),
    }

    def pstr(x: float) -> str:
        return str(x).replace(".", "p").replace("-", "m")

    beta_a_R_values = [2.5, 3.5, 5.0, 7.5]
    beta_a_T_values = [1.2, 2.0, 3.0, 5.0]

    for beta_R in beta_a_R_values:
        for beta_T in beta_a_T_values:
            name = f"T0p5_c0p005_ratio_1e2_betaR{pstr(beta_R)}_betaT{pstr(beta_T)}"
            parameter_sets[name] = replace(
                base,
                lambda_a_T=5e-1,
                lambda_a_c=5e-3,
                beta_a_R=beta_R,
                beta_a_T=beta_T,
                max_proposals=80_000_000,
            )

    for beta_R in beta_a_R_values:
        for beta_T in beta_a_T_values:
            name = f"T1e-1_c1e-4_ratio_1e3_betaR{pstr(beta_R)}_betaT{pstr(beta_T)}"
            parameter_sets[name] = replace(
                base,
                lambda_a_T=1e-1,
                lambda_a_c=1e-4,
                beta_a_R=beta_R,
                beta_a_T=beta_T,
                max_proposals=120_000_000,
            )

    scenarios: Dict[str, Parameters] = {}
    for parameter_name, p0 in parameter_sets.items():
        for rep_idx, seed in enumerate(replicate_seeds, start=1):
            scenarios[f"{parameter_name}__rep{rep_idx:02d}"] = replace(p0, seed=seed)
    return scenarios



# ============================================================
# 7. Stage 2: Death/ERK-only parameter sweep
# ============================================================


def make_stage2_death_erk_scenarios(base: Parameters) -> Dict[str, Parameters]:
    """
    Build the Stage 2 death/ERK-only sweep from the report.

    Fixed activation regime:
        lambda_a_T = 0.5,
        lambda_a_c = 0.005,
        beta_a_R   = 2.5,
        beta_a_T   = 3.0.

    Swept death / ERK parameters:
        lambda_d  in {0.5, 1.0, 2.0},
        beta_d_R  in {1.0, 2.0, 4.0},
        beta_d_T  in {0.4, 0.8, 1.6}.

    This gives 3 x 3 x 3 = 27 parameter settings.
    Each setting is repeated over 5 independent seeds, so the pooled
    binomial test for each setting uses N = 5 x 1000 = 5000 analyzed deaths.
    """
    n_replicates = 5
    replicate_seeds = [1001 + i for i in range(n_replicates)]

    activation_reference = replace(
        base,
        lambda_a_T=0.5,
        lambda_a_c=0.005,
        beta_a_R=2.5,
        beta_a_T=3.0,
        max_proposals=80_000_000,
    )

    lambda_d_values = [0.5, 1.0, 2.0]
    beta_d_R_values = [1.0, 2.0, 4.0]
    beta_d_T_values = [0.4, 0.8, 1.6]

    parameter_sets: Dict[str, Parameters] = {}
    for lambda_d in lambda_d_values:
        for beta_d_R in beta_d_R_values:
            for beta_d_T in beta_d_T_values:
                name = (
                    f"stage2_deathERK_ld{pstr(lambda_d)}"
                    f"_bdR{pstr(beta_d_R)}"
                    f"_bdT{pstr(beta_d_T)}"
                )
                parameter_sets[name] = replace(
                    activation_reference,
                    lambda_d=lambda_d,
                    beta_d_R=beta_d_R,
                    beta_d_T=beta_d_T,
                )

    scenarios: Dict[str, Parameters] = {}
    for parameter_name, p0 in parameter_sets.items():
        for rep_idx, seed in enumerate(replicate_seeds, start=1):
            scenarios[f"{parameter_name}__rep{rep_idx:02d}"] = replace(p0, seed=seed)

    return scenarios


def run_parameter_sweep(
    *,
    scenarios: Dict[str, Parameters],
    out_dir: Path,
    run_label: str,
) -> Tuple[Path, Path, Path, Path, Path]:
    """
    Run a parameter sweep and save the same statistical outputs for every sweep.

    This helper is shared by:
        - Stage 1 exploratory sweep;
        - Stage 2 death/ERK-only sweep.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {run_label} ===")
    print("Burn-in 500 deaths, then analyze 1000 deaths per replicate.\n")

    results: List[SimulationResult] = []
    for name, params in scenarios.items():
        print(f"Scenario: {name}")
        result = simulate_one(name, params)
        results.append(result)
        print(f"  total observed deaths simulated: {result.n_observed_deaths_total}")
        print(f"  analyzed observed deaths:        {result.n_observed_deaths_analyzed}")
        print(f"  burn-in time:                    {result.burn_in_time:.6f}")
        print(f"  time for analyzed 1000 deaths:   {result.analysis_time_span:.6f}")
        print("")

    global_colorbar_max = int(max(np.max(spatial_histogram(result, result.params.bins_space)) for result in results))
    global_colorbar_max = max(global_colorbar_max, 1)

    summary_rows: List[Dict[str, Any]] = []
    for result in results:
        stats = death_statistics(result)
        summary_rows.append(stats)

        save_death_events_csv(result, out_dir)
        plot_time_histogram(result, out_dir)
        plot_spatial_histogram_2d(result, out_dir, global_colorbar_max)
        plot_cumulative_deaths(result, out_dir)

        print(f"Scenario: {result.scenario}")
        print(f"  activation acceptance:              {stats['activation_acceptance_rate']:.4f}")
        print(f"  death rejection rate:               {stats['death_rejection_rate']:.4f}")
        print(f"    outside active:                   {stats['death_rejection_rate_outside_active']:.4f}")
        print(f"    by ERK protection:                {stats['death_rejection_rate_by_erk']:.4f}")
        print(f"  T-zone deaths k/n:                  {int(stats['t_zone_binomial_test_k'])}/{int(stats['t_zone_binomial_test_n'])}")
        print(f"  T-zone observed fraction:           {stats['t_zone_binomial_observed_fraction']:.4f}")
        print(f"  T-zone area fraction p0:            {stats['t_zone_binomial_null_probability_area_fraction']:.4f}")
        print(f"  binomial p-value greater:           {stats['t_zone_binomial_p_value_greater']:.4g}")
        print(f"  T-zone density ratio:               {stats['t_zone_density_ratio']:.4f}")
        print("")

    summary_path = save_summary_csv(summary_rows, out_dir)
    aggregated_rows = aggregate_replicate_rows(summary_rows)
    aggregate_summary_path = save_aggregate_summary_csv(aggregated_rows, out_dir)
    pooled_csv_path, pooled_md_path, pooled_rows = save_pooled_bonferroni_outputs(summary_rows, out_dir)
    report_path = write_analysis_report(summary_rows, out_dir, global_colorbar_max, aggregated_rows, pooled_rows)

    print(f"{run_label} done.")
    print(f"Common colorbar vmax: {global_colorbar_max}")
    print(f"Replicate-level summary CSV: {summary_path}")
    print(f"Aggregated-by-parameter summary CSV: {aggregate_summary_path}")
    print(f"Pooled Bonferroni CSV: {pooled_csv_path}")
    print(f"Pooled Bonferroni report: {pooled_md_path}")
    print(f"T-zone binomial analysis report: {report_path}")
    print(f"Plots and death-event CSV files are in: {out_dir.resolve()}")

    return summary_path, aggregate_summary_path, pooled_csv_path, pooled_md_path, report_path


# ============================================================
# 7. Optional diagnostic experiments matching the Gillespie diagnostic file
# ============================================================


RUN_DIAGNOSTIC_EXPERIMENTS = False
N_DIAGNOSTIC_REPLICATES = 10
DIAGNOSTIC_SEED0 = 7001
DIAGNOSTIC_OUTPUT_DIR = Path("diagnostic_outputs_rejection_Tzone_mechanisms")
DIAGNOSTIC_MAX_PROPOSALS = 120_000_000

DIAGNOSTIC_LAMBDA_D_VALUES = [0.25, 0.5, 1.0, 2.0, 4.0]
DIAGNOSTIC_ERK_REGIMES = [
    {"erk_regime": "strong_ERK_protection", "beta_d_R": 1.0, "beta_d_T": 0.4, "meaning": "large ERK radius, long ERK lifetime"},
    {"erk_regime": "reference_ERK", "beta_d_R": 2.0, "beta_d_T": 0.8, "meaning": "reference ERK protection"},
    {"erk_regime": "weak_ERK_protection", "beta_d_R": 4.0, "beta_d_T": 1.6, "meaning": "small ERK radius, short ERK lifetime"},
]
DIAGNOSTIC_ACTIVATION_REGIMES = [
    {"activation_regime": "non_local_activation", "beta_a_R": 1.0, "beta_a_T": 1.2, "meaning": "large and long-lived active disks"},
    {"activation_regime": "medium_activation", "beta_a_R": 2.5, "beta_a_T": 3.0, "meaning": "current visible reference regime"},
    {"activation_regime": "local_activation", "beta_a_R": 5.0, "beta_a_T": 3.0, "meaning": "small and short-lived active disks"},
]
DIAGNOSTIC_FIXED_LAMBDA_A_T = 0.5
DIAGNOSTIC_FIXED_LAMBDA_A_C = 0.005


def pstr(x: float) -> str:
    return str(x).replace(".", "p").replace("-", "m")


def seed_for(rep_idx: int, offset: int = 0) -> int:
    return DIAGNOSTIC_SEED0 + offset + rep_idx


def t_zone_stats_from_arrays(x: np.ndarray, y: np.ndarray, p: Parameters) -> Dict[str, float | int]:
    n = int(len(x))
    area_T = float(t_zone_area(p))
    area_W = float(p.area_W)
    p0 = area_T / area_W if area_W > 0 else np.nan
    if n == 0:
        return {
            "n_observed_deaths_analyzed": 0,
            "t_zone_binomial_test_k": 0,
            "t_zone_binomial_test_n": 0,
            "t_zone_area_fraction": p0,
            "t_zone_observed_fraction": np.nan,
            "t_zone_density_ratio": np.nan,
            "t_zone_binomial_p_value_greater": np.nan,
        }

    in_T = np.array([is_inside_T_zone(float(xx), float(yy), p) for xx, yy in zip(x, y)], dtype=bool)
    k = int(np.sum(in_T))
    n_out = n - k
    area_out = area_W - area_T
    frac = k / n
    density_T = k / area_T if area_T > 0 else np.nan
    density_out = n_out / area_out if area_out > 0 else np.nan
    ratio = density_T / density_out if density_out > 0 else np.nan
    pval = binomial_upper_tail_p_value(k=k, n=n, p0=p0)
    return {
        "n_observed_deaths_analyzed": n,
        "t_zone_binomial_test_k": k,
        "t_zone_binomial_test_n": n,
        "t_zone_area_fraction": p0,
        "t_zone_observed_fraction": frac,
        "t_zone_density_ratio": ratio,
        "t_zone_binomial_p_value_greater": pval,
    }


def activations_during_window(result: SimulationResult, t0: float, t1: float) -> int:
    act_t = np.asarray(result.activation_t, dtype=float)
    if act_t.size == 0:
        return 0
    return int(np.sum((act_t > t0) & (act_t <= t1)))


def diagnostic_row_from_fixed_death_result(
    result: SimulationResult,
    *,
    experiment: str,
    base_scenario: str,
    replicate: int,
    activation_regime: str,
    erk_regime: str,
    meaning: str = "",
) -> Dict[str, Any]:
    p = result.params
    stats = t_zone_stats_from_arrays(np.asarray(result.death_x), np.asarray(result.death_y), p)
    analysis_start = float(result.burn_in_time)
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
            if float(result.analysis_time_span) > 0 else np.nan
        ),
        "accepted_activations_total": int(result.n_accepted_activations),
        "accepted_activations_during_analysis": activations_during_window(result, analysis_start, analysis_end),
        "death_candidates": int(result.n_death_candidates),
        "death_rejected": int(result.n_death_rejected),
        "death_rejected_outside_active": int(result.n_death_rejected_outside_active),
        "death_rejected_by_erk": int(result.n_death_rejected_by_erk),
        "activation_candidates": int(result.n_activation_candidates),
        "activation_rejected": int(result.n_activation_rejected),
        "all_proposals": int(result.n_proposals),
        "n_blocks": int(result.n_blocks),
        "stopped_by_target": bool(result.stopped_by_target),
        **stats,
    }


def simulate_fixed_time_after_burnin(scenario: str, p: Parameters, observation_time: float) -> Dict[str, Any]:
    """Rejection-sampling version of fixed physical-time control after burn-in."""
    validate_params(p)
    rng = np.random.default_rng(p.seed)

    active_x: List[float] = []
    active_y: List[float] = []
    active_r: List[float] = []
    active_end: List[float] = []

    erk_x: List[float] = []
    erk_y: List[float] = []
    erk_r: List[float] = []
    erk_end: List[float] = []

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
    n_death_rejected_outside_active = 0
    n_death_rejected_by_erk = 0
    n_proposals = 0
    n_blocks_used = 0

    burn_in_done = p.burn_in_deaths == 0
    burn_in_time = 0.0 if burn_in_done else np.nan
    t_stop = observation_time if burn_in_done else np.inf

    t0 = 0.0
    for block_index in range(p.max_blocks):
        n_blocks_used = block_index + 1
        t1 = t0 + p.block_duration

        n_act = int(rng.poisson(p.lambda_a_1 * p.area_W * p.block_duration))
        act_t = rng.uniform(t0, t1, size=n_act)
        act_x = rng.uniform(0.0, p.Lx, size=n_act)
        act_y = rng.uniform(0.0, p.Ly, size=n_act)
        act_u = rng.uniform(0.0, 1.0, size=n_act)
        act_r = rng.exponential(scale=1.0 / p.beta_a_R, size=n_act)
        act_tau = rng.exponential(scale=1.0 / p.beta_a_T, size=n_act)

        n_death = int(rng.poisson(p.lambda_d * p.area_W * p.block_duration))
        death_cand_t = rng.uniform(t0, t1, size=n_death)
        death_cand_x = rng.uniform(0.0, p.Lx, size=n_death)
        death_cand_y = rng.uniform(0.0, p.Ly, size=n_death)
        death_cand_r = rng.exponential(scale=1.0 / p.beta_d_R, size=n_death)
        death_cand_tau = rng.exponential(scale=1.0 / p.beta_d_T, size=n_death)

        n_activation_candidates += n_act
        n_death_candidates += n_death
        n_proposals += n_act + n_death
        if n_proposals > p.max_proposals:
            break

        event_t = np.concatenate([act_t, death_cand_t])
        event_type = np.concatenate([np.zeros(n_act, dtype=np.int8), np.ones(n_death, dtype=np.int8)])
        event_index = np.concatenate([np.arange(n_act, dtype=np.int64), np.arange(n_death, dtype=np.int64)])

        for pos in np.argsort(event_t):
            t = float(event_t[pos])
            if burn_in_done and t > t_stop:
                break

            kind = int(event_type[pos])
            idx = int(event_index[pos])
            prune_expired_disks(t, active_x, active_y, active_r, active_end)
            prune_expired_disks(t, erk_x, erk_y, erk_r, erk_end)

            if kind == 0:
                x = float(act_x[idx])
                y = float(act_y[idx])
                u = float(act_u[idx])
                r = float(act_r[idx])
                tau = float(act_tau[idx])
                inside_active = is_inside_current_disks(x, y, active_x, active_y, active_r)
                if inside_active:
                    accept_prob = 1.0
                elif is_inside_T_zone(x, y, p):
                    accept_prob = p.lambda_a_T / p.lambda_a_1
                else:
                    accept_prob = p.lambda_a_c / p.lambda_a_1
                if u <= accept_prob:
                    activation_t.append(t)
                    active_x.append(x)
                    active_y.append(y)
                    active_r.append(r)
                    active_end.append(t + tau)
                else:
                    n_activation_rejected += 1
            else:
                x = float(death_cand_x[idx])
                y = float(death_cand_y[idx])
                r = float(death_cand_r[idx])
                tau = float(death_cand_tau[idx])
                inside_active = is_inside_current_disks(x, y, active_x, active_y, active_r)
                inside_erk = is_inside_current_disks(x, y, erk_x, erk_y, erk_r)
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
                    erk_x.append(x)
                    erk_y.append(y)
                    erk_r.append(r)
                    erk_end.append(t + tau)
                else:
                    if inside_erk:
                        n_death_rejected_by_erk += 1
                    else:
                        n_death_rejected_outside_active += 1

        if burn_in_done and t1 >= t_stop:
            break
        t0 = t1

    if not burn_in_done:
        burn_in_time = float(death_t_all[-1]) if death_t_all else 0.0
        t_stop = burn_in_time

    death_x_arr = np.asarray(death_x_an, dtype=float)
    death_y_arr = np.asarray(death_y_an, dtype=float)
    death_t_arr = np.asarray(death_t_an, dtype=float)
    stats = t_zone_stats_from_arrays(death_x_arr, death_y_arr, p)
    activation_t_arr = np.asarray(activation_t, dtype=float)
    activations_window = int(np.sum((activation_t_arr > burn_in_time) & (activation_t_arr <= t_stop)))
    n_death_rejected = n_death_rejected_outside_active + n_death_rejected_by_erk

    return {
        "scenario": scenario,
        "seed": p.seed,
        "burn_in_time": float(burn_in_time),
        "observation_time_after_burnin": float(observation_time),
        "final_time": float(t_stop),
        "analysis_time_span": float(max(0.0, t_stop - burn_in_time)),
        "accepted_activations_total": int(len(activation_t_arr)),
        "accepted_activations_during_analysis": activations_window,
        "death_candidates": int(n_death_candidates),
        "death_rejected": int(n_death_rejected),
        "death_rejected_outside_active": int(n_death_rejected_outside_active),
        "death_rejected_by_erk": int(n_death_rejected_by_erk),
        "activation_candidates": int(n_activation_candidates),
        "activation_rejected": int(n_activation_rejected),
        "all_proposals": int(n_proposals),
        "n_blocks": int(n_blocks_used),
        "stopped_by_time": bool(burn_in_done),
        **stats,
    }


def grouped_summary(rows: List[Dict[str, Any]], group_cols: List[str]) -> List[Dict[str, Any]]:
    groups: Dict[Tuple[Any, ...], List[Dict[str, Any]]] = {}
    for row in rows:
        key = tuple(row.get(c) for c in group_cols)
        groups.setdefault(key, []).append(row)

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
    out_rows: List[Dict[str, Any]] = []
    for key, group in groups.items():
        out = {col: val for col, val in zip(group_cols, key)}
        for metric in metrics:
            vals = [_safe_float(r.get(metric, np.nan)) for r in group]
            mean, std, min_v, max_v = _mean_std_min_max(vals)
            out[f"{metric}_mean"] = mean
            out[f"{metric}_std"] = std
            out[f"{metric}_min"] = min_v
            out[f"{metric}_max"] = max_v
        out_rows.append(out)
    return out_rows


def save_rows_csv(rows: List[Dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def pooled_bonferroni_generic(rows: List[Dict[str, Any]], group_cols: List[str], alpha: float = 0.05) -> List[Dict[str, Any]]:
    groups: Dict[Tuple[Any, ...], List[Dict[str, Any]]] = {}
    for row in rows:
        key = tuple(row.get(c) for c in group_cols)
        groups.setdefault(key, []).append(row)

    out: List[Dict[str, Any]] = []
    for key, group in groups.items():
        row = {col: val for col, val in zip(group_cols, key)}
        pooled_k = int(sum(int(r["t_zone_binomial_test_k"]) for r in group))
        pooled_n = int(sum(int(r["t_zone_binomial_test_n"]) for r in group))
        p0 = float(group[0]["t_zone_area_fraction"])
        pooled_fraction = pooled_k / pooled_n if pooled_n > 0 else np.nan
        pval = binomial_upper_tail_p_value(k=pooled_k, n=pooled_n, p0=p0)
        n_out = pooled_n - pooled_k
        density_ratio = (pooled_k / p0) / (n_out / (1.0 - p0)) if pooled_n > 0 and n_out > 0 and 0 < p0 < 1 else np.nan
        row.update({
            "n_replicates": int(len(group)),
            "pooled_t_zone_deaths_k": pooled_k,
            "pooled_total_deaths_n": pooled_n,
            "t_zone_area_fraction_p0": p0,
            "pooled_t_zone_fraction": pooled_fraction,
            "pooled_t_zone_density_ratio": density_ratio,
            "pooled_binomial_p_value_greater": pval,
        })
        out.append(row)

    m = len(out)
    bonf_alpha = alpha / m if m > 0 else np.nan
    for row in out:
        pval = float(row["pooled_binomial_p_value_greater"])
        row["bonferroni_alpha"] = bonf_alpha
        row["bonferroni_adjusted_p_value"] = min(pval * m, 1.0) if np.isfinite(pval) else np.nan
        row["significant_after_bonferroni_0p05"] = bool(np.isfinite(pval) and pval < bonf_alpha)
    out.sort(key=lambda r: (not bool(r["significant_after_bonferroni_0p05"]), _safe_float(r["pooled_binomial_p_value_greater"])))
    return out


def spearman_trend_tests(rows: List[Dict[str, Any]], group_cols: List[str], x_col: str, y_cols: List[str], label: str) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[Any, ...], List[Dict[str, Any]]] = {}
    for row in rows:
        key = tuple(row.get(c) for c in group_cols)
        grouped.setdefault(key, []).append(row)

    out: List[Dict[str, Any]] = []
    for key, group in grouped.items():
        base = {col: val for col, val in zip(group_cols, key)}
        for y_col in y_cols:
            pairs = []
            for r in group:
                x = _safe_float(r.get(x_col, np.nan))
                y = _safe_float(r.get(y_col, np.nan))
                if np.isfinite(x) and np.isfinite(y):
                    pairs.append((x, y))
            if len(pairs) < 3 or len(set(x for x, _ in pairs)) < 2 or len(set(y for _, y in pairs)) < 2 or scipy_spearmanr is None:
                rho, pval = np.nan, np.nan
            else:
                arr = np.asarray(pairs, dtype=float)
                rho, pval = scipy_spearmanr(arr[:, 0], arr[:, 1])
            out.append({
                "experiment": label,
                **base,
                "x": x_col,
                "y": y_col,
                "n": int(len(pairs)),
                "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "spearman_p_value": float(pval) if np.isfinite(pval) else np.nan,
            })
    return out


def run_diagnostic_experiments() -> None:
    """Run the same conceptual diagnostics as the Gillespie diagnostic script, using rejection sampling."""
    out_dir = DIAGNOSTIC_OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    base = Parameters(
        burn_in_deaths=500,
        target_deaths=1000,
        lambda_a_T=DIAGNOSTIC_FIXED_LAMBDA_A_T,
        lambda_a_c=DIAGNOSTIC_FIXED_LAMBDA_A_C,
        max_proposals=DIAGNOSTIC_MAX_PROPOSALS,
    )

    rows_A: List[Dict[str, Any]] = []
    fixed_activation = {"activation_regime": "medium_activation", "beta_a_R": 2.5, "beta_a_T": 3.0}
    print("\n=== Rejection diagnostic Experiment 1A: fixed 1000 deaths lambda_d sweep ===")
    for er_idx, er in enumerate(DIAGNOSTIC_ERK_REGIMES):
        for ld_idx, lambda_d in enumerate(DIAGNOSTIC_LAMBDA_D_VALUES):
            base_scenario = f"E1A_{er['erk_regime']}_ld{pstr(lambda_d)}"
            for rep in range(1, N_DIAGNOSTIC_REPLICATES + 1):
                p = replace(
                    base,
                    seed=seed_for(rep, offset=10_000 + 1000 * er_idx + 100 * ld_idx),
                    lambda_d=lambda_d,
                    beta_d_R=er["beta_d_R"],
                    beta_d_T=er["beta_d_T"],
                    beta_a_R=fixed_activation["beta_a_R"],
                    beta_a_T=fixed_activation["beta_a_T"],
                )
                scenario = f"{base_scenario}__rep{rep:02d}"
                print(scenario)
                result = simulate_one(scenario, p)
                rows_A.append(diagnostic_row_from_fixed_death_result(
                    result,
                    experiment="E1A_fixed_1000_deaths",
                    base_scenario=base_scenario,
                    replicate=rep,
                    activation_regime=fixed_activation["activation_regime"],
                    erk_regime=er["erk_regime"],
                    meaning=er["meaning"],
                ))

    save_rows_csv(rows_A, out_dir / "experiment1A_fixed_1000_deaths_replicates.csv")
    save_rows_csv(grouped_summary(rows_A, ["erk_regime", "beta_d_R", "beta_d_T", "lambda_d"]), out_dir / "experiment1A_fixed_1000_deaths_group_summary.csv")
    save_rows_csv(pooled_bonferroni_generic(rows_A, ["erk_regime", "beta_d_R", "beta_d_T", "lambda_d"]), out_dir / "experiment1A_fixed_1000_deaths_pooled_bonferroni.csv")

    ref_rows = [r for r in rows_A if r["erk_regime"] == "reference_ERK" and math.isclose(float(r["lambda_d"]), 1.0)]
    reference_observation_time = float(np.mean([_safe_float(r["analysis_time_span"]) for r in ref_rows])) if ref_rows else np.nan
    if not np.isfinite(reference_observation_time) or reference_observation_time <= 0:
        reference_observation_time = float(np.median([_safe_float(r["analysis_time_span"]) for r in rows_A]))
    print(f"Reference fixed physical observation window: {reference_observation_time:.6f}")

    rows_B: List[Dict[str, Any]] = []
    print("\n=== Rejection diagnostic Experiment 1B: fixed physical-time lambda_d control ===")
    for er_idx, er in enumerate(DIAGNOSTIC_ERK_REGIMES):
        for ld_idx, lambda_d in enumerate(DIAGNOSTIC_LAMBDA_D_VALUES):
            base_scenario = f"E1B_{er['erk_regime']}_ld{pstr(lambda_d)}"
            for rep in range(1, N_DIAGNOSTIC_REPLICATES + 1):
                p = replace(
                    base,
                    seed=seed_for(rep, offset=20_000 + 1000 * er_idx + 100 * ld_idx),
                    lambda_d=lambda_d,
                    beta_d_R=er["beta_d_R"],
                    beta_d_T=er["beta_d_T"],
                    beta_a_R=fixed_activation["beta_a_R"],
                    beta_a_T=fixed_activation["beta_a_T"],
                )
                scenario = f"{base_scenario}__rep{rep:02d}"
                print(scenario)
                row = simulate_fixed_time_after_burnin(scenario, p, reference_observation_time)
                row.update({
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
                })
                row["observed_death_rate_after_burnin"] = row["n_observed_deaths_analyzed"] / row["analysis_time_span"] if row["analysis_time_span"] > 0 else np.nan
                rows_B.append(row)

    save_rows_csv(rows_B, out_dir / "experiment1B_fixed_physical_time_replicates.csv")
    save_rows_csv(grouped_summary(rows_B, ["erk_regime", "beta_d_R", "beta_d_T", "lambda_d"]), out_dir / "experiment1B_fixed_physical_time_group_summary.csv")
    save_rows_csv(pooled_bonferroni_generic(rows_B, ["erk_regime", "beta_d_R", "beta_d_T", "lambda_d"]), out_dir / "experiment1B_fixed_physical_time_pooled_bonferroni.csv")

    trend_rows: List[Dict[str, Any]] = []
    trend_rows.extend(spearman_trend_tests(
        rows_A,
        group_cols=["erk_regime", "beta_d_R", "beta_d_T"],
        x_col="lambda_d",
        y_cols=["analysis_time_span", "accepted_activations_during_analysis", "t_zone_observed_fraction", "t_zone_density_ratio"],
        label="E1A_fixed_1000_deaths",
    ))
    trend_rows.extend(spearman_trend_tests(
        rows_B,
        group_cols=["erk_regime", "beta_d_R", "beta_d_T"],
        x_col="lambda_d",
        y_cols=["n_observed_deaths_analyzed", "accepted_activations_during_analysis", "t_zone_observed_fraction", "t_zone_density_ratio"],
        label="E1B_fixed_physical_time",
    ))
    save_rows_csv(trend_rows, out_dir / "experiment1_spearman_trend_tests.csv")

    rows_E2: List[Dict[str, Any]] = []
    print("\n=== Rejection diagnostic Experiment 2: activation localization x death/ERK ===")
    for ar_idx, ar in enumerate(DIAGNOSTIC_ACTIVATION_REGIMES):
        for er_idx, er in enumerate(DIAGNOSTIC_ERK_REGIMES):
            for ld_idx, lambda_d in enumerate(DIAGNOSTIC_LAMBDA_D_VALUES):
                base_scenario = f"E2_{ar['activation_regime']}_{er['erk_regime']}_ld{pstr(lambda_d)}"
                for rep in range(1, N_DIAGNOSTIC_REPLICATES + 1):
                    p = replace(
                        base,
                        seed=seed_for(rep, offset=30_000 + 5000 * ar_idx + 1000 * er_idx + 100 * ld_idx),
                        lambda_d=lambda_d,
                        beta_d_R=er["beta_d_R"],
                        beta_d_T=er["beta_d_T"],
                        beta_a_R=ar["beta_a_R"],
                        beta_a_T=ar["beta_a_T"],
                    )
                    scenario = f"{base_scenario}__rep{rep:02d}"
                    print(scenario)
                    result = simulate_one(scenario, p)
                    rows_E2.append(diagnostic_row_from_fixed_death_result(
                        result,
                        experiment="E2_activation_localization_x_deathERK",
                        base_scenario=base_scenario,
                        replicate=rep,
                        activation_regime=ar["activation_regime"],
                        erk_regime=er["erk_regime"],
                        meaning=f"activation: {ar['meaning']}; ERK: {er['meaning']}",
                    ))

    save_rows_csv(rows_E2, out_dir / "experiment2_activation_localization_deathERK_replicates.csv")
    save_rows_csv(grouped_summary(rows_E2, ["activation_regime", "beta_a_R", "beta_a_T", "erk_regime", "beta_d_R", "beta_d_T", "lambda_d"]), out_dir / "experiment2_activation_localization_deathERK_group_summary.csv")
    save_rows_csv(pooled_bonferroni_generic(rows_E2, ["activation_regime", "beta_a_R", "beta_a_T", "erk_regime", "beta_d_R", "beta_d_T", "lambda_d"]), out_dir / "experiment2_activation_localization_deathERK_pooled_bonferroni.csv")
    save_rows_csv(spearman_trend_tests(
        rows_E2,
        group_cols=["activation_regime", "erk_regime", "beta_d_R", "beta_d_T"],
        x_col="lambda_d",
        y_cols=["analysis_time_span", "accepted_activations_during_analysis", "t_zone_observed_fraction", "t_zone_density_ratio"],
        label="E2_activation_localization_x_deathERK",
    ), out_dir / "experiment2_lambda_d_trend_tests_by_activation_and_ERK.csv")

    report = (
        "# Rejection-sampling diagnostic experiments for T-zone mechanisms\n\n"
        "This folder reproduces the Gillespie diagnostic logic with the rejection-sampling simulator.\n\n"
        "- Experiment 1A: fixed 1000 analyzed deaths after burn-in, sweeping lambda_d and ERK regimes.\n"
        "- Experiment 1B: fixed physical-time control after burn-in using the reference observation window from Experiment 1A.\n"
        "- Experiment 2: death/ERK robustness conditional on activation localization.\n\n"
        "Central evidence is always the T-zone density ratio, one-sided exact binomial test, pooled binomial test, "
        "Bonferroni correction, and Spearman trend diagnostics.\n"
    )
    (out_dir / "diagnostic_experiments_report.md").write_text(report, encoding="utf-8")
    print(f"Diagnostic outputs saved in: {out_dir.resolve()}")


# ============================================================
# 8. Main
# ============================================================

# Default behavior:
#   - Run all three stages in one file:
#       Stage 1: exploratory rejection-sampling sweep;
#       Stage 2: death/ERK-only rejection-sampling sweep;
#       Stage 3: diagnostic experiments.
#
RUN_STAGE1_EXPLORATORY_SWEEP = True
RUN_STAGE2_DEATH_ERK_SWEEP = True
RUN_DIAGNOSTIC_EXPERIMENTS = True


def main() -> None:
    base = Parameters(burn_in_deaths=500, target_deaths=1000)

    if RUN_STAGE1_EXPLORATORY_SWEEP:
        stage1_scenarios = make_scenarios(base)
        run_parameter_sweep(
            scenarios=stage1_scenarios,
            out_dir=Path("simulation_outputs_1000_rejection_stage1_exploratory"),
            run_label="Stage 1 exploratory rejection-sampling sweep",
        )
    else:
        print("Stage 1 exploratory sweep is included but not run by default.")
        print("Set RUN_STAGE1_EXPLORATORY_SWEEP = True to run it.\n")

    if RUN_STAGE2_DEATH_ERK_SWEEP:
        stage2_scenarios = make_stage2_death_erk_scenarios(base)
        run_parameter_sweep(
            scenarios=stage2_scenarios,
            out_dir=Path("simulation_outputs_1000_rejection_stage2_death_ERK_sweep"),
            run_label="Stage 2 death/ERK-only rejection-sampling sweep",
        )
    else:
        print("Stage 2 death/ERK sweep is included but not run.")
        print("Set RUN_STAGE2_DEATH_ERK_SWEEP = True to run it.\n")

    if RUN_DIAGNOSTIC_EXPERIMENTS:
        run_diagnostic_experiments()
    else:
        print("\nDiagnostic experiments are included but not run by default.")
        print("Set RUN_DIAGNOSTIC_EXPERIMENTS = True to reproduce fixed-time / activation-localization diagnostics.")


if __name__ == "__main__":
    main()
