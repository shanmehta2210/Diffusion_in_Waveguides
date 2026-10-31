"""Steering with a few heaters: least-squares powers and the steering range.

The target is a uniform ramp φ_n = Δφ·x_n/d, which points the beam at
sin θ = Δφ·λ/(2πd). We look for heater powers P ≥ 0 (heaters can only heat)
and a free offset c (a phase shared by every emitter does not move the beam)
that minimise
    ||G·P - (φ_target + c)||².
Writing c = c+ - c- with both parts ≥ 0 turns this into one non-negative
least-squares problem in (P, c+, c-).

Scaling: every constraint has the form "≥ 0", so if (P, c) is the best fit to a
target, (αP, αc) is the best fit to α times that target, for any α > 0. The
residual, the power and the offset all grow in proportion to Δφ. One fit to a
ramp of Δφ = 1 rad per emitter therefore describes every angle in that
direction; the opposite direction is a different target (the ramp flips) and
needs its own fit.
"""
import numpy as np
from scipy.optimize import nnls

from .constants import PITCH_UM, WAVELENGTH_UM

STREHL_MIN = 0.8  # beam-quality threshold for the steering range


def fit_powers(G, target):
    """Best non-negative powers P (W) and free offset c (rad) with G·P ≈ target + c."""
    n, k = G.shape
    ones = np.ones((n, 1))
    solution, _ = nnls(np.hstack([G, -ones, ones]), target)
    return solution[:k], solution[k] - solution[k + 1]


def strehl_vs_step(residual_per_rad, dphi):
    """Strehl ratio toward the target angle, for each ramp step in dphi.

    At step Δφ the phases are Δφ·(target + c + r), where r is the residual of
    the unit fit. Toward the target angle the ramp term cancels exactly (that
    is what aiming means) and c is common to all emitters, so only Δφ·r is left:
        S(Δφ) = |mean_n exp(i Δφ r_n)|².
    This equals farfield.strehl at the target angle (checked in check_pipeline.py).
    """
    return np.abs(np.exp(1j * np.outer(dphi, residual_per_rad)).mean(axis=1)) ** 2


def steer_one_direction(G, x_emit, sign, pitch=PITCH_UM, strehl_min=STREHL_MIN):
    """Fit a unit ramp in one direction (sign = +1 or -1) and find its steering range.

    The range is the largest angle up to which the Strehl ratio stays at or above
    strehl_min, capped at the single-beam field of view arcsin(λ/2d), where Δφ = π.
    """
    target = sign * x_emit / pitch
    powers, offset = fit_powers(G, target)
    residual = G @ powers - target - offset

    dphi = np.linspace(0, np.pi, 1001)
    s = strehl_vs_step(residual, dphi)
    below = np.nonzero(s < strehl_min)[0]
    if below.size == 0:
        dphi_max = np.pi
    else:  # interpolate the first crossing
        j = below[0]
        dphi_max = np.interp(strehl_min, [s[j], s[j - 1]], [dphi[j], dphi[j - 1]])

    u_per_rad = WAVELENGTH_UM / (2 * np.pi * pitch)  # sin θ per rad of Δφ
    return {
        "powers_per_rad": powers,              # W per rad of Δφ, one entry per heater
        "power_per_rad": powers.sum(),         # total W per rad of Δφ
        "offset_per_rad": offset,
        "residual_per_rad": residual,          # rad per rad of Δφ
        "rms_residual_per_rad": np.sqrt(np.mean(residual**2)),
        "dphi_max": dphi_max,
        "range_deg": np.degrees(np.arcsin(dphi_max * u_per_rad)),
        # Near broadside dθ ≈ d(sin θ), so this is the cost of each degree there.
        "mw_per_deg": powers.sum() / u_per_rad * np.pi / 180 * 1e3,
    }


def steer_both_directions(G, x_emit, pitch=PITCH_UM, strehl_min=STREHL_MIN):
    """Fits for both directions. span_deg = range toward +θ plus range toward -θ.

    Broadside needs no power and always has Strehl 1, so the two ranges join
    there and the span is one continuous stretch of angles.
    """
    plus = steer_one_direction(G, x_emit, +1, pitch, strehl_min)
    minus = steer_one_direction(G, x_emit, -1, pitch, strehl_min)
    return {
        "plus": plus,
        "minus": minus,
        "span_deg": plus["range_deg"] + minus["range_deg"],
        "mw_per_deg": 0.5 * (plus["mw_per_deg"] + minus["mw_per_deg"]),
    }
