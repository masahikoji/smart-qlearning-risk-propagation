#!/usr/bin/env python3
"""Paper-profile-aligned exact-integration validation for smooth boundaries and exact ties.

This program reproduces the two-stage smooth-boundary and exact-tie diagnostics using the
same statistical specification as the principal SMART simulations:

* scenarios: 2s_smooth_terminal_local and 2s_exact_tie_terminal;
* sample sizes: 250, 500, and 1000 by default;
* 20,000 independently simulated fitted data sets by default;
* project base seed: 260906;
* the same terminal DGP, narrow/wide candidate bases, AIC rule, and Akaike
  weights as the R simulation engine;
* smooth-boundary local-signal calibration delta=1; the exact-tie contrast is identically zero.

The evaluation-law expectation is computed by exact conditional integration
over X2|X1 followed by Gauss-Legendre integration over X1.  Thus N_eval=100000
from the principal risk simulations is replaced deliberately by deterministic
quadrature for this boundary diagnostic.  The default quadrature order 48 is
the paper-profile value and is numerically indistinguishable from order 128 in
validation checks.

Replicate-specific seeds use the same stable scenario/n/replicate schedule as
the R project.  NumPy and R do not generate identical random streams, so this
is an independent implementation, but its results are invariant to the order
of --n values, --batch, and --jobs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss

DEFAULT_N = (250, 500, 1000)
DEFAULT_REPS = 20_000
DEFAULT_BASE_SEED = 260_906
DEFAULT_QUADRATURE_ORDER = 48
STABILIZE_TOL = 1.0e-10


@dataclass(frozen=True)
class Scenario:
    name: str
    scenario_id: str
    boundary: str
    rho_x: float = 0.30
    xi_x: float = 0.00
    beta0_terminal: float = 0.25
    beta1_terminal: float = 0.50
    beta_past_a1: float = 0.20
    tau0_terminal: float = 0.00
    tau1_terminal: float = 1.00
    sigma_terminal: float = 1.00
    target_delta_terminal: float = 1.00


def scenario_by_name(name: str) -> Scenario:
    if name == "smooth":
        return Scenario(
            name="smooth boundary", scenario_id="2s_smooth_terminal_local",
            boundary="smooth", tau1_terminal=1.0, target_delta_terminal=1.0
        )
    if name == "tie":
        return Scenario(
            name="exact tie", scenario_id="2s_exact_tie_terminal",
            boundary="exact_tie", tau1_terminal=0.0, target_delta_terminal=0.0
        )
    raise ValueError(name)


def scenario_hash(identifier: str) -> int:
    """Exact translation of simulation/main/R/01_utils.R::scenario_hash."""
    h = 104_729
    for value in identifier.encode("utf-8"):
        h = (h * 1009 + value + 97) % 2_147_483_000
    return int(h)


def replicate_seed(base_seed: int, scenario_id: str, n: int, rep_id: int) -> int:
    """Exact integer seed schedule used by the R simulation project."""
    z = base_seed + 1_000_003 * scenario_hash(scenario_id) + 1009 * n + rep_id
    return int(z % 2_147_483_000)


def even_moments_x2(rho_x: float, xi_x: float = 0.0) -> tuple[float, float]:
    """Second/fourth moments from R/00_config.R for the two-stage state."""
    h = 1.0 - rho_x
    m2 = rho_x**2 / 3.0 + h**2 / 3.0 + xi_x**2
    m4 = (
        rho_x**4 / 5.0
        + h**4 / 5.0
        + xi_x**4
        + 6.0 * rho_x**2 * h**2 / 9.0
        + 6.0 * rho_x**2 * xi_x**2 / 3.0
        + 6.0 * h**2 * xi_x**2 / 3.0
    )
    return m2, m4


def calibrated_b(scenario: Scenario) -> float:
    if scenario.boundary == "exact_tie":
        return 0.0
    m2, m4 = even_moments_x2(scenario.rho_x, scenario.xi_x)
    residual_variance = m4 - m2 * m2
    if not residual_variance > 0:
        raise ValueError("Nonpositive residual variance for the added AX2 regressor")
    return (
        scenario.target_delta_terminal
        * scenario.sigma_terminal
        / math.sqrt(residual_variance)
    )


def real_roots(a0: float, a1: float, a2: float, tol: float = 1e-13) -> list[float]:
    if abs(a2) < tol:
        return [] if abs(a1) < tol else [-a0 / a1]
    disc = a1 * a1 - 4.0 * a2 * a0
    if disc < -tol:
        return []
    if abs(disc) <= tol:
        return [-a1 / (2.0 * a2)]
    root_disc = math.sqrt(max(0.0, disc))
    return sorted(((-a1 - root_disc) / (2.0 * a2), (-a1 + root_disc) / (2.0 * a2)))


def quadratic_antiderivative(x: np.ndarray, a0: float, a1: float, a2: float) -> np.ndarray:
    return a0 * x + 0.5 * a1 * x**2 + (a2 / 3.0) * x**3


def mean_switching_uniform(
    centers: np.ndarray,
    halfwidth: float,
    fitted_half_contrast: tuple[float, float, float],
    true_half_contrast: tuple[float, float, float],
) -> np.ndarray:
    """E[|C_hat(X)|-sign(C_star(X)) C_hat(X) | center] exactly."""
    a0, a1, a2 = fitted_half_contrast
    b0, b1, b2 = true_half_contrast
    roots = sorted(
        set(
            root
            for root in (*real_roots(a0, a1, a2), *real_roots(b0, b1, b2))
            if math.isfinite(root)
        )
    )
    cuts = (-math.inf, *roots, math.inf)
    lower_support = centers - halfwidth
    upper_support = centers + halfwidth
    result = np.zeros_like(centers, dtype=float)

    for left, right in zip(cuts[:-1], cuts[1:]):
        lower = np.maximum(lower_support, left)
        upper = np.minimum(upper_support, right)
        active = upper > lower
        if not np.any(active):
            continue
        if math.isfinite(left) and math.isfinite(right):
            midpoint = 0.5 * (left + right)
        elif math.isfinite(right):
            midpoint = right - 1.0
        elif math.isfinite(left):
            midpoint = left + 1.0
        else:
            midpoint = 0.0

        fitted_mid = a0 + a1 * midpoint + a2 * midpoint**2
        true_mid = b0 + b1 * midpoint + b2 * midpoint**2
        fitted_sign = 1 if fitted_mid >= 0 else -1
        true_sign = 1 if true_mid >= 0 else -1
        multiplier = fitted_sign - true_sign
        if multiplier == 0:
            continue
        result[active] += multiplier * (
            quadratic_antiderivative(upper[active], a0, a1, a2)
            - quadratic_antiderivative(lower[active], a0, a1, a2)
        )

    return np.maximum(result / (2.0 * halfwidth), 0.0)


def design_matrices(
    x2: np.ndarray, a1: np.ndarray, a2: np.ndarray, scenario: Scenario
) -> tuple[np.ndarray, np.ndarray]:
    """Exact scenario-specific narrow/wide bases from R/02_features.R."""
    common = [np.ones_like(x2), x2, x2**2]
    if scenario.boundary == "smooth":
        common.append(x2**3)
    common.extend([a1, a2, a2 * x2])
    narrow = np.stack(tuple(common), axis=-1)
    wide = np.concatenate((narrow, (a2 * x2**2)[..., None]), axis=-1)
    return narrow, wide


def batched_ols(X: np.ndarray, y: np.ndarray, stabilize_tol: float) -> tuple[np.ndarray, np.ndarray]:
    """OLS with the same empirical-Gram threshold as the R production code."""
    n = X.shape[1]
    xtx = np.einsum("bni,bnj->bij", X, X, optimize=True)
    gram = xtx / n
    minimum_eigenvalue = np.linalg.eigvalsh(gram)[:, 0]
    invalid = (~np.isfinite(minimum_eigenvalue)) | (minimum_eigenvalue < stabilize_tol)
    if np.any(invalid):
        bad = np.flatnonzero(invalid)
        raise RuntimeError(
            f"Spectral stabilization would be activated in {len(bad)} replicate(s); "
            f"first batch index={int(bad[0])}, min eigenvalue={minimum_eigenvalue[bad[0]]:.3e}."
        )
    xty = np.einsum("bni,bn->bi", X, y, optimize=True)
    beta = np.linalg.solve(xtx, xty[..., None])[..., 0]
    residual = y - np.einsum("bni,bi->bn", X, beta, optimize=True)
    rss = np.einsum("bn,bn->b", residual, residual, optimize=True)
    return beta, np.maximum(rss, np.finfo(float).tiny)


def aggregate_aic(
    beta_narrow: np.ndarray,
    rss_narrow: np.ndarray,
    beta_wide: np.ndarray,
    rss_wide: np.ndarray,
    n: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Hard AIC selection and normalized Akaike weighting from R/03_fit.R."""
    p0 = beta_narrow.shape[1]
    p1 = beta_wide.shape[1]
    embedded_narrow = np.zeros((beta_narrow.shape[0], p1), dtype=float)
    embedded_narrow[:, :p0] = beta_narrow
    aic_narrow = n * np.log(rss_narrow / n) + 2.0 * p0
    aic_wide = n * np.log(rss_wide / n) + 2.0 * p1
    choose_wide = aic_wide < aic_narrow
    selected = np.where(choose_wide[:, None], beta_wide, embedded_narrow)
    wide_weight = 1.0 / (
        1.0 + np.exp(np.clip(0.5 * (aic_wide - aic_narrow), -700.0, 700.0))
    )
    averaged = (1.0 - wide_weight)[:, None] * embedded_narrow + wide_weight[:, None] * beta_wide
    return selected, averaged


def make_replicate_batch(
    scenario: Scenario,
    n: int,
    rep_ids: Iterable[int],
    base_seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rep_ids = list(rep_ids)
    size = len(rep_ids)
    x1 = np.empty((size, n), dtype=float)
    u2 = np.empty((size, n), dtype=float)
    a1 = np.empty((size, n), dtype=float)
    a2 = np.empty((size, n), dtype=float)
    epsilon = np.empty((size, n), dtype=float)
    actions = np.array((-1.0, 1.0))

    for row, rep_id in enumerate(rep_ids):
        seed = replicate_seed(base_seed, scenario.scenario_id, n, rep_id)
        # Explicit MT19937 documents the independent Python RNG engine.
        rng = np.random.Generator(np.random.MT19937(seed))
        x1[row] = rng.uniform(-1.0, 1.0, n)
        u2[row] = rng.uniform(-1.0, 1.0, n)
        a1[row] = rng.choice(actions, n)
        a2[row] = rng.choice(actions, n)
        epsilon[row] = rng.normal(0.0, scenario.sigma_terminal, n)
    return x1, u2, a1, a2, epsilon


def summarize(values: np.ndarray) -> tuple[float, float]:
    return float(np.mean(values)), float(np.std(values, ddof=1) / math.sqrt(values.size))


def simulate_one_n(
    n: int,
    reps: int,
    batch_size: int,
    scenario: Scenario,
    base_seed: int,
    quadrature_order: int,
    raw_directory: Path | None,
) -> list[dict[str, float | int | str]]:
    nodes, weights = leggauss(quadrature_order)
    evaluation_weights = 0.5 * weights  # X1 ~ Uniform(-1,1)
    centers = scenario.rho_x * nodes  # xi_x=0, hence no A1 shift
    halfwidth = 1.0 - scenario.rho_x
    b = calibrated_b(scenario)
    true_half_contrast = (
        scenario.tau0_terminal,
        scenario.tau1_terminal,
        b / math.sqrt(n),
    )

    selection_values = np.empty(reps, dtype=float)
    averaging_values = np.empty(reps, dtype=float)

    for start in range(1, reps + 1, batch_size):
        stop = min(reps, start + batch_size - 1)
        rep_ids = range(start, stop + 1)
        x1, u2, a1, a2, epsilon = make_replicate_batch(
            scenario, n, rep_ids, base_seed
        )
        x2 = scenario.rho_x * x1 + (1.0 - scenario.rho_x) * u2 + scenario.xi_x * a1
        local_treatment = (
            scenario.tau0_terminal
            + scenario.tau1_terminal * x2
            + (b / math.sqrt(n)) * x2**2
        )
        y2 = (
            scenario.beta0_terminal
            + scenario.beta1_terminal * x2
            + scenario.beta_past_a1 * a1
            + a2 * local_treatment
            + epsilon
        )

        X_narrow, X_wide = design_matrices(x2, a1, a2, scenario)
        beta_narrow, rss_narrow = batched_ols(X_narrow, y2, STABILIZE_TOL)
        beta_wide, rss_wide = batched_ols(X_wide, y2, STABILIZE_TOL)
        beta_selection, beta_averaging = aggregate_aic(
            beta_narrow, rss_narrow, beta_wide, rss_wide, n
        )

        offset = start - 1
        for local_index, (coef_selection, coef_averaging) in enumerate(
            zip(beta_selection, beta_averaging)
        ):
            for coef, destination in (
                (coef_selection, selection_values),
                (coef_averaging, averaging_values),
            ):
                if scenario.boundary == "smooth":
                    fitted_half_contrast = (float(coef[5]), float(coef[6]), float(coef[7]))
                else:
                    fitted_half_contrast = (float(coef[4]), float(coef[5]), float(coef[6]))
                transported = mean_switching_uniform(
                    centers,
                    halfwidth,
                    fitted_half_contrast,
                    true_half_contrast,
                )
                scale = n * n if scenario.boundary == "smooth" else n
                destination[offset + local_index] = scale * float(
                    np.dot(evaluation_weights, transported * transported)
                )

    if raw_directory is not None:
        raw_directory.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            raw_directory / f"n{n}.npz",
            rep_id=np.arange(1, reps + 1, dtype=np.int32),
            aic_selection=selection_values,
            akaike_weighting=averaging_values,
        )

    selection_mean, selection_mcse = summarize(selection_values)
    averaging_mean, averaging_mcse = summarize(averaging_values)
    common = {
        "scenario": scenario.name,
        "scenario_id": scenario.scenario_id,
        "n": n,
        "reps": reps,
        "base_seed": base_seed,
        "quadrature_order": quadrature_order,
        "evaluation_method": "exact conditional integration + Gauss-Legendre",
        "stabilize_tol": STABILIZE_TOL,
        "stabilization_count": 0,
        "b": b,
    }
    return [
        {
            **common,
            "procedure": "AIC selection",
            "estimate": selection_mean,
            "mcse": selection_mcse,
        },
        {
            **common,
            "procedure": "Akaike weighting",
            "estimate": averaging_mean,
            "mcse": averaging_mcse,
        },
    ]


def _task(arguments: tuple) -> list[dict[str, float | int | str]]:
    return simulate_one_n(*arguments)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=DEFAULT_REPS)
    parser.add_argument("--batch", type=int, default=50, help="Computational only; results are batch-invariant.")
    parser.add_argument("--n", type=int, nargs="+", default=list(DEFAULT_N))
    parser.add_argument("--scenario", choices=["smooth", "tie", "both"], default="both")
    parser.add_argument(
        "--include-extended-4000",
        action="store_true",
        help="Append n=4000 as an optional large-sample check; not part of the principal paper grid.",
    )
    parser.add_argument("--quadrature-order", type=int, default=DEFAULT_QUADRATURE_ORDER)
    parser.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED)
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 1) - 2))
    parser.add_argument("--out", type=Path, default=Path("boundary_validation_paper.csv"))
    parser.add_argument("--raw-directory", type=Path, default=Path("boundary_validation_raw"))
    args = parser.parse_args()

    if args.reps < 2:
        parser.error("--reps must be at least 2")
    if args.batch < 1:
        parser.error("--batch must be positive")
    if args.quadrature_order < 16:
        parser.error("--quadrature-order must be at least 16")
    if any(n < 20 for n in args.n):
        parser.error("all sample sizes must be at least 20")

    n_values = list(dict.fromkeys(args.n))
    if args.include_extended_4000 and 4000 not in n_values:
        n_values.append(4000)

    scenario_names = ["smooth", "tie"] if args.scenario == "both" else [args.scenario]
    scenarios = [scenario_by_name(name) for name in scenario_names]
    raw_directory = args.raw_directory
    tasks = [
        (n, args.reps, args.batch, scenario, args.base_seed,
         args.quadrature_order, raw_directory / scenario.scenario_id)
        for scenario in scenarios for n in n_values
    ]

    blocks: list[list[dict[str, float | int | str]]]
    if args.jobs <= 1 or len(tasks) == 1:
        blocks = [_task(task) for task in tasks]
    else:
        blocks = []
        with ProcessPoolExecutor(max_workers=min(args.jobs, len(tasks))) as executor:
            futures = {executor.submit(_task, task): task[0] for task in tasks}
            for future in as_completed(futures):
                blocks.append(future.result())

    rows = [row for block in blocks for row in block]
    procedure_order = {"AIC selection": 0, "Akaike weighting": 1}
    scenario_order = {"smooth boundary": 0, "exact tie": 1}
    rows.sort(key=lambda row: (scenario_order[str(row["scenario"])], int(row["n"]), procedure_order[str(row["procedure"])]))

    for row in rows:
        print(
            f"{row['scenario_id']} n={int(row['n']):4d} "
            f"{str(row['procedure']):18s} "
            f"{float(row['estimate']):.4f} ({float(row['mcse']):.4f})",
            flush=True,
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "scenario",
        "scenario_id",
        "n",
        "procedure",
        "reps",
        "estimate",
        "mcse",
        "base_seed",
        "quadrature_order",
        "evaluation_method",
        "stabilize_tol",
        "stabilization_count",
        "b",
    ]
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "program": Path(__file__).name,
        "numpy_version": np.__version__,
        "rng_engine": "numpy.random.MT19937 with R-project replicate seed schedule",
        "scenarios": [asdict(scenario) for scenario in scenarios],
        "n_values": n_values,
        "principal_n_values": list(DEFAULT_N),
        "reps": args.reps,
        "base_seed": args.base_seed,
        "quadrature_order": args.quadrature_order,
        "batch": args.batch,
        "jobs": args.jobs,
        "note": (
            "The evaluation-law expectation is deterministic quadrature rather than "
            "the N_eval=100000 Monte Carlo scoring sample used in principal risk tables."
        ),
    }
    metadata_path = args.out.with_suffix(args.out.suffix + ".json")
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
