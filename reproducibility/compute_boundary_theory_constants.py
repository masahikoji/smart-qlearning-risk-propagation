#!/usr/bin/env python3
"""Compute the smooth-boundary limiting constants reported in the supplement.

The calculation implements exactly the displayed formula

    f0^2 E[(B_A^r)^4],
    B_A^r = Z + m2*b*{1 - w_r(Xi^2) Xi},

with Z ~ N(0,1), Xi ~ N(1,1), Z independent of Xi,
selection weight w_SEL(x)=1{x>2}, and Akaike model-averaging weight
w_MA(x)={1+exp(1-x/2)}^{-1}.

The outer expectation over Z is evaluated analytically using
E[(Z+c)^4] = c^4 + 6 c^2 + 3; the remaining one-dimensional
expectation over Xi is evaluated by adaptive quadrature.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from scipy.integrate import quad

RHO_X = 0.30
F0 = 1.0 / 1.4


def even_moments_x2(rho_x: float) -> tuple[float, float]:
    h = 1.0 - rho_x
    m2 = rho_x**2 / 3.0 + h**2 / 3.0
    m4 = rho_x**4 / 5.0 + h**4 / 5.0 + 6.0 * rho_x**2 * h**2 / 9.0
    return m2, m4


def normal_pdf(x: float, mean: float = 1.0, sd: float = 1.0) -> float:
    z = (x - mean) / sd
    return math.exp(-0.5 * z * z) / (sd * math.sqrt(2.0 * math.pi))


def w_selection(x2: float) -> float:
    return 1.0 if x2 > 2.0 else 0.0


def w_akaike(x2: float) -> float:
    # Stable logistic evaluation for 1 / {1 + exp(1 - x2/2)}.
    z = 1.0 - 0.5 * x2
    if z >= 0:
        ez = math.exp(-z)
        return ez / (1.0 + ez)
    ez = math.exp(z)
    return 1.0 / (1.0 + ez)


def constant(weight_fn) -> float:
    m2, m4 = even_moments_x2(RHO_X)
    omega = m4 - m2 * m2
    b = 1.0 / math.sqrt(omega)
    m2b = m2 * b

    def integrand(xi: float) -> float:
        c = m2b * (1.0 - weight_fn(xi * xi) * xi)
        fourth = c**4 + 6.0 * c**2 + 3.0
        return fourth * normal_pdf(xi)

    value, error = quad(integrand, -math.inf, math.inf, epsabs=1e-12, epsrel=1e-12, limit=500)
    return F0 * F0 * value, error, m2, m4, omega, b, m2b


def main() -> None:
    sel, sel_err, m2, m4, omega, b, m2b = constant(w_selection)
    ma, ma_err, *_ = constant(w_akaike)
    out = {
        "rho_x": RHO_X,
        "m2": m2,
        "m4": m4,
        "omega": omega,
        "b": b,
        "f0": F0,
        "m2_times_b": m2b,
        "aic_selection_constant": sel,
        "akaike_model_averaging_constant": ma,
        "quadrature_error_selection": sel_err,
        "quadrature_error_averaging": ma_err,
    }
    print(json.dumps(out, indent=2))
    Path("boundary_theory_constants.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
