"""Fabrication errors: what K smooth heaters can correct, against static trim and per-element calibration.

Every waveguide picks up a fixed phase error ε_n from fabrication: width and
thickness vary along its own path, so its effective index does too. Two parts:
  - local: independent from waveguide to waveguide. Along a path of independent
    length L it accumulates like a random walk, σ² ≈ L / L_coh, with
    L_coh = 4.17 mm for strip and 1.61 mm for rib waveguides [Yang 2015];
  - smooth: wafer- and die-level variation, on-chip temperature gradients.
    Modelled here as a Gaussian process across the array with a
    squared-exponential correlation of length ℓ [Xing 2022 for the split].

Ways to deal with ε:
  - uncorrected: steer with the nominal drive; the error stays on the beam.
  - corrected by the K heaters: refit the drive to the target minus ε
    (non-negative least squares with a free offset, as in steering.py). Heat
    only makes smooth profiles, so this removes the part of ε the heaters can
    reach: offsets always, gradients and smooth bumps if the layout can make them.
  - static trim: a one-time per-waveguide correction after fabrication removes
    ε; the heaters then only steer.
  - conventional (one isolated heater per emitter, phases wrapped): corrects any
    ε exactly, and the only cost is power.

Analytic checks (run_fabrication.py tests both):
  - Ruze with the finite-N term: for independent Gaussian errors with nothing
    corrected, E[S] = exp(-σ²) + (1 - exp(-σ²))/N.
  - The exact Gaussian expectation with the design's own residual and the
    least-squares projection (expected_strehl). In words: a fit onto an
    m-dimensional space (K heater shapes plus the offset) removes on average
    only m/N of independent error variance, so few heaters barely help.
"""
import numpy as np

from .steering import fit_powers


def strehl(residual):
    """|mean_n exp(i r_n)|² for each row of residual phases; a shared offset doesn't matter."""
    return np.abs(np.exp(1j * np.atleast_2d(residual)).mean(axis=1)) ** 2


def ruze(sigma, n):
    """Expected Strehl for independent Gaussian phase errors of std sigma on n emitters."""
    s = np.exp(-np.asarray(sigma, float) ** 2)
    return s + (1 - s) / n


def expected_strehl(r0, Q, sigma):
    """Exact mean Strehl for residual r0 + Q·ε with ε independent Gaussian, std sigma.

    For u = r0 + Qε the covariance is Σ = σ²QQᵀ, and for a Gaussian
    E[exp(i(u_n - u_m))] = exp(i(r0_n - r0_m)) · exp(-(Σ_nn + Σ_mm - 2Σ_nm)/2). So
        E[S] = (1/N²) Σ_nm a_n a_m* exp(-(Σ_nn + Σ_mm - 2Σ_nm)/2),  a = exp(i r0).
    Uncorrected: Q = I. Corrected by an unconstrained least-squares fit onto
    the columns of A = [G, 1]: Q = I - A·A⁺, and r0 = -Q·target (the design's own
    residual). The non-negative fit equals it while no heater hits zero power.
    """
    S = sigma**2 * Q @ Q.T
    d = np.diag(S)
    a = np.exp(1j * r0)
    return float(np.real(np.sum(np.outer(a, a.conj()) * np.exp(-(d[:, None] + d[None, :] - 2 * S) / 2)))
                 / len(r0) ** 2)


def phase_errors(x, sigma, ell, n_trials, rng):
    """n_trials draws of ε (rad) on positions x (µm): std sigma, correlation length ell.

    ell = 0 gives independent errors. Otherwise the correlation between two
    waveguides a distance Δx apart is exp(-Δx² / 2ℓ²).
    """
    if ell == 0:
        return rng.normal(0.0, sigma, (n_trials, len(x)))
    C = np.exp(-0.5 * ((x[:, None] - x[None, :]) / ell) ** 2)
    L = np.linalg.cholesky(C + 1e-9 * np.eye(len(x)))   # tiny jitter: C is near singular for large ell
    return sigma * rng.standard_normal((n_trials, len(x))) @ L.T


def uncorrected(G, P_nominal, target, errors):
    """Strehl per trial when the nominal drive is kept and the errors stay."""
    return strehl(G @ P_nominal + errors - target)


def corrected(G, target, errors):
    """Strehl and total drive power (W) per trial after refitting the heaters to target - ε."""
    s, p = np.empty(len(errors)), np.empty(len(errors))
    for i, e in enumerate(errors):
        P, c = fit_powers(G, target - e)
        s[i] = strehl(G @ P + e - target - c)[0]
        p[i] = P.sum()
    return s, p


def wrapped_power(phi, p_pi):
    """Least power (W) for isolated heaters to add the phases phi modulo 2π, over all common offsets.

    The best common offset puts some emitter exactly at 0, so the candidates are
    c = -φ_m mod 2π (as in conventional.py and lens.py).
    """
    offsets = np.append(np.mod(-phi, 2 * np.pi), 0.0)
    return min(p_pi / np.pi * np.mod(phi + c, 2 * np.pi).sum() for c in offsets)
