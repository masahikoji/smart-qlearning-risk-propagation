#!/usr/bin/env python3
"""Exact-quadrature validation for the two-stage boundary design.

The program simulates the terminal-stage Q-regression used in the manuscript,
constructs AIC hard-selection and normalized Akaike-weighting estimators, and
computes the transported switching remainder P_1 rho_2 by analytic integration
over X_2 | X_1.  It therefore avoids the fixed-descendant-grid error in the
older decomposition diagnostic.
"""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from numpy.polynomial.legendre import leggauss


@dataclass(frozen=True)
class Scenario:
    name: str
    smooth: bool
    rho: float = 0.30
    sigma: float = 1.0


def even_moments_x2(rho: float) -> tuple[float, float]:
    h = 1.0 - rho
    ex2 = (rho * rho + h * h) / 3.0
    ex4 = (rho**4 + h**4) / 5.0 + 6.0 * rho * rho * h * h / 9.0
    return ex2, ex4


def calibrated_b(rho: float) -> float:
    ex2, ex4 = even_moments_x2(rho)
    return 1.0 / math.sqrt(ex4 - ex2 * ex2)


def real_roots(a0: float, a1: float, a2: float, tol: float = 1e-13) -> list[float]:
    if abs(a2) < tol:
        if abs(a1) < tol:
            return []
        return [-a0 / a1]
    disc = a1 * a1 - 4.0 * a2 * a0
    if disc < -tol:
        return []
    if abs(disc) <= tol:
        return [-a1 / (2.0 * a2)]
    s = math.sqrt(max(0.0, disc))
    return sorted([(-a1 - s) / (2.0 * a2), (-a1 + s) / (2.0 * a2)])


def q_antiderivative(x: np.ndarray, a0: float, a1: float, a2: float) -> np.ndarray:
    return a0 * x + 0.5 * a1 * x**2 + (a2 / 3.0) * x**3


def mean_switching_uniform(
    centers: np.ndarray,
    halfwidth: float,
    chat: tuple[float, float, float],
    cstar: tuple[float, float, float] | None,
    tie_sign: int = 1,
) -> np.ndarray:
    """E[ |C_hat(X)| - sign(C_star(X)) C_hat(X) | center ]."""
    a0, a1, a2 = chat
    true_is_tie = cstar is None or max(abs(z) for z in cstar) < 1e-13
    roots = real_roots(a0, a1, a2)
    if not true_is_tie:
        roots += real_roots(*cstar)
    roots = sorted(set(r for r in roots if math.isfinite(r)))
    cuts = [-math.inf, *roots, math.inf]
    lo_all = centers - halfwidth
    hi_all = centers + halfwidth
    ans = np.zeros_like(centers, dtype=float)

    for left, right in zip(cuts[:-1], cuts[1:]):
        lo = np.maximum(lo_all, left)
        hi = np.minimum(hi_all, right)
        ok = hi > lo
        if not np.any(ok):
            continue
        if math.isfinite(left) and math.isfinite(right):
            mid = 0.5 * (left + right)
        elif not math.isfinite(left) and math.isfinite(right):
            mid = right - 1.0
        elif math.isfinite(left) and not math.isfinite(right):
            mid = left + 1.0
        else:
            mid = 0.0
        ch_mid = a0 + a1 * mid + a2 * mid * mid
        sh = 1 if ch_mid >= 0 else -1
        if true_is_tie:
            st = tie_sign
        else:
            b0, b1, b2 = cstar
            cs_mid = b0 + b1 * mid + b2 * mid * mid
            st = tie_sign if abs(cs_mid) < 1e-13 else (1 if cs_mid > 0 else -1)
        k = sh - st
        if k == 0:
            continue
        piece = k * (q_antiderivative(hi[ok], a0, a1, a2) - q_antiderivative(lo[ok], a0, a1, a2))
        ans[ok] += piece
    ans /= 2.0 * halfwidth
    return np.maximum(ans, 0.0)


def design_matrices(x: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    common = [
        np.ones_like(x),
        x,
        x**2,
        x**3,
        a1,
        a2,
        a2 * x,
    ]
    narrow = np.stack(common, axis=-1)
    wide = np.concatenate([narrow, (a2 * x**2)[..., None]], axis=-1)
    return narrow, wide


def batched_ols(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xtx = np.einsum("bni,bnj->bij", X, X, optimize=True)
    xty = np.einsum("bni,bn->bi", X, y, optimize=True)
    beta = np.linalg.solve(xtx, xty[..., None])[..., 0]
    resid = y - np.einsum("bni,bi->bn", X, beta, optimize=True)
    rss = np.einsum("bn,bn->b", resid, resid, optimize=True)
    return beta, rss


def embed_and_aggregate(
    beta0: np.ndarray,
    rss0: np.ndarray,
    beta1: np.ndarray,
    rss1: np.ndarray,
    n: int,
) -> tuple[np.ndarray, np.ndarray]:
    B = beta0.shape[0]
    emb0 = np.zeros((B, 8), dtype=float)
    emb0[:, :7] = beta0
    emb1 = beta1
    aic0 = n * np.log(rss0 / n) + 2.0 * 7
    aic1 = n * np.log(rss1 / n) + 2.0 * 8
    choose_wide = aic1 < aic0
    sel = np.where(choose_wide[:, None], emb1, emb0)
    # Stable two-model softmax.
    w1 = 1.0 / (1.0 + np.exp(np.clip(0.5 * (aic1 - aic0), -700, 700)))
    ma = (1.0 - w1)[:, None] * emb0 + w1[:, None] * emb1
    return sel, ma


def summarize(x: np.ndarray) -> tuple[float, float]:
    return float(np.mean(x)), float(np.std(x, ddof=1) / math.sqrt(len(x)))


def simulate_one_n(
    n: int,
    reps: int,
    batch: int,
    scenario: Scenario,
    seed: int,
    outer_order: int,
) -> list[dict[str, float | int | str]]:
    rng = np.random.default_rng(seed)
    nodes, weights = leggauss(outer_order)
    outer_w = 0.5 * weights  # expectation for X1 ~ Uniform(-1,1)
    centers = scenario.rho * nodes
    halfwidth = 1.0 - scenario.rho
    b = calibrated_b(scenario.rho)
    cstar = (0.0, 1.0, b / math.sqrt(n)) if scenario.smooth else None

    vals_sel: list[np.ndarray] = []
    vals_ma: list[np.ndarray] = []
    completed = 0
    while completed < reps:
        B = min(batch, reps - completed)
        x1 = rng.uniform(-1.0, 1.0, size=(B, n))
        u2 = rng.uniform(-1.0, 1.0, size=(B, n))
        x2 = scenario.rho * x1 + (1.0 - scenario.rho) * u2
        a1 = rng.choice(np.array([-1.0, 1.0]), size=(B, n))
        a2 = rng.choice(np.array([-1.0, 1.0]), size=(B, n))
        eps = rng.normal(0.0, scenario.sigma, size=(B, n))
        if scenario.smooth:
            treatment = x2 + (b / math.sqrt(n)) * x2**2
        else:
            treatment = np.zeros_like(x2)
        y2 = 0.25 + 0.50 * x2 + 0.20 * a1 + a2 * treatment + eps

        X0, X1 = design_matrices(x2, a1, a2)
        beta0, rss0 = batched_ols(X0, y2)
        beta1, rss1 = batched_ols(X1, y2)
        sel, ma = embed_and_aggregate(beta0, rss0, beta1, rss1, n)

        out_sel = np.empty(B, dtype=float)
        out_ma = np.empty(B, dtype=float)
        for j in range(B):
            # union order: int,X,X2,X3,pastA,A,AX,AX2
            for coef, out in ((sel[j], out_sel), (ma[j], out_ma)):
                chat = (float(coef[5]), float(coef[6]), float(coef[7]))
                m = mean_switching_uniform(centers, halfwidth, chat, cstar, tie_sign=1)
                risk = float(np.dot(outer_w, m * m))
                out[j] = (n * n if scenario.smooth else n) * risk
        vals_sel.append(out_sel)
        vals_ma.append(out_ma)
        completed += B

    arr_sel = np.concatenate(vals_sel)
    arr_ma = np.concatenate(vals_ma)
    m_sel, se_sel = summarize(arr_sel)
    m_ma, se_ma = summarize(arr_ma)
    return [
        {
            "scenario": scenario.name,
            "n": n,
            "procedure": "AIC selection",
            "reps": reps,
            "estimate": m_sel,
            "mcse": se_sel,
            "b": b,
        },
        {
            "scenario": scenario.name,
            "n": n,
            "procedure": "Akaike weighting",
            "reps": reps,
            "estimate": m_ma,
            "mcse": se_ma,
            "b": b,
        },
    ]


def _task(args_tuple):
    n, reps, batch, sc, seed, outer_order = args_tuple
    return simulate_one_n(n, reps, batch, sc, seed, outer_order)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20000)
    ap.add_argument("--batch", type=int, default=100)
    ap.add_argument("--n", type=int, nargs="+", default=[250, 500, 1000])
    ap.add_argument("--outer-order", type=int, default=128)
    ap.add_argument("--scenario", choices=["smooth", "tie", "both"], default="both")
    ap.add_argument("--seed", type=int, default=261016)
    ap.add_argument("--jobs", type=int, default=min(6, os.cpu_count() or 1))
    ap.add_argument("--out", type=Path, default=Path("table_s7_exact_quadrature.csv"))
    args = ap.parse_args()

    rows: list[dict[str, float | int | str]] = []
    scenarios_all = {"smooth": Scenario("smooth boundary", True), "tie": Scenario("exact tie", False)}
    scenarios = [scenarios_all["smooth"], scenarios_all["tie"]] if args.scenario == "both" else [scenarios_all[args.scenario]]
    tasks=[]
    for si, sc in enumerate(scenarios):
        for ni, n in enumerate(args.n):
            tasks.append((n,args.reps,args.batch,sc,args.seed+100000*si+1000*ni,args.outer_order))
    if args.jobs <= 1:
        blocks=[_task(z) for z in tasks]
    else:
        blocks=[]
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs={ex.submit(_task,z):z for z in tasks}
            for fut in as_completed(futs):
                blocks.append(fut.result())
    for block in blocks:
        rows.extend(block)
        for row in block:
            print(f"{row['scenario']:15s} n={int(row['n']):4d} {row['procedure']:18s} {row['estimate']:.4f} ({row['mcse']:.4f})", flush=True)
    order={(sc.name,n,p):i for i,(sc,n,p) in enumerate((sc,n,p) for sc in scenarios for n in args.n for p in ["AIC selection","Akaike weighting"])}
    rows.sort(key=lambda r: order[(r["scenario"],r["n"],r["procedure"])])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["scenario", "n", "procedure", "reps", "estimate", "mcse", "b"])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
