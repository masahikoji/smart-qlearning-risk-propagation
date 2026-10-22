#!/usr/bin/env python3
"""Reproduce the risk-ordering and backward-reversal figure."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm


def aic_gap(delta):
    rt2 = np.sqrt(2.0)

    def integrand(x):
        sel = 1.0 if x * x > 2.0 else 0.0
        ma = 1.0 / (1.0 + np.exp(1.0 - x * x / 2.0))
        return ((sel * x - delta) ** 2 - (ma * x - delta) ** 2) * norm.pdf(
            x, loc=delta, scale=1.0
        )

    value, _ = quad(
        integrand,
        -10.0,
        10.0,
        points=[-rt2, rt2],
        epsabs=1e-10,
        epsrel=1e-10,
        limit=200,
    )
    return value


def main():
    outdir = Path(__file__).resolve().parent

    delta_grid = np.linspace(0.0, 4.0, 161)
    gap_grid = np.array([aic_gap(d) for d in delta_grid])
    delta_star = brentq(aic_gap, 2.0, 3.0)

    stage_2 = np.array([1, 2])
    gap_2 = np.array([0.144, -0.195])
    stage_3 = np.array([1, 2, 3])
    gap_3 = np.array([0.246, 0.136, -0.213])

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))

    ax = axes[0]
    ax.plot(delta_grid, gap_grid, linewidth=1.8)
    ax.axhline(0.0, linewidth=0.8)
    ax.axvline(delta_star, linestyle="--", linewidth=0.9)
    ax.text(delta_star + 0.05, 0.40, r"$\delta^*\approx 2.365$", fontsize=8, va="top")
    ax.set_xlabel(r"Local signal $\delta$")
    ax.set_ylabel("Limiting risk gap")
    ax.set_title("(a) Terminal one-stage benchmark", fontsize=9)
    ax.set_xlim(0, 4)
    ax.set_ylim(-0.32, 0.48)
    ax.tick_params(labelsize=8)

    ax = axes[1]
    ax.plot(stage_2, gap_2, marker="o", linewidth=1.5, label="2 stages")
    ax.plot(stage_3, gap_3, marker="s", linewidth=1.5, label="3 stages")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel("Target stage")
    ax.set_ylabel(r"$n\widehat{\Delta}_s$ at $\delta=3$")
    ax.set_title("(b) Backward risk reversal", fontsize=9)
    ax.set_xticks([1, 2, 3])
    ax.set_xlim(0.75, 3.25)
    ax.set_ylim(-0.32, 0.33)
    ax.tick_params(labelsize=8)
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    fig.tight_layout(w_pad=2.2)
    fig.savefig(outdir / "fig_risk_reversal.pdf", bbox_inches="tight")
    fig.savefig(outdir / "fig_risk_reversal.png", dpi=220, bbox_inches="tight")

    print(f"delta_star={delta_star:.10f}")
    print("benchmark values:", [(d, aic_gap(d)) for d in [0.5, 1.0, 2.4, 3.0]])


if __name__ == "__main__":
    main()
