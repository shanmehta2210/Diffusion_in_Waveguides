"""A thermal cylindrical lens: fitting a centre-hot quadratic phase, and the conventional baselines.

A lens of focal length f across the array needs, in the paraxial limit,
    φ_n = c - k0·x_n²/(2f),
so the centre leads the edges by a = k0·X²/(2f), with X the half-aperture
(32 µm). We call a the lens strength (rad): a·f = k0·X²/2 ≈ 2076 rad·µm, so
a = 1 rad is f ≈ 2.1 mm and a = 2π is f ≈ 0.33 mm. dn/dT > 0, so heating the
centre more than the edges makes a converging lens, like convex glass.

Fitting works exactly as for steering (opa/steering.py): heaters can only heat
and the offset c is free, so the best fit to a·t is a times the best fit to t.
One unit fit gives the power per rad of lens strength and, from its residual,
the Strehl ratio at the focus for every a.

The paraxial form is fine for the weak lenses where this design is
competitive (f ≳ 0.3 mm, NA ≲ 0.1). Below f ≈ 3X ≈ 100 µm it is only approximate.
"""
import numpy as np

from .constants import K0
from .steering import STREHL_MIN, fit_powers, strehl_vs_step


def lens_constant(x_emit):
    """a·f in rad·µm: lens strength a gives focal length f = lens_constant / a."""
    return K0 * np.max(np.abs(x_emit)) ** 2 / 2


def unit_target(x_emit):
    """Centre-hot quadratic with edge-to-centre phase 1 rad (offset is free in the fit)."""
    return -(x_emit / np.max(np.abs(x_emit))) ** 2


def fit_lens(G, x_emit, strehl_min=STREHL_MIN, a_max=60.0):
    """Best non-negative fit of a unit lens; returns power per rad and the largest usable strength.

    a_limit is the strength at which the Strehl at the focus falls to strehl_min.
    """
    target = unit_target(x_emit)
    powers, offset = fit_powers(G, target)
    residual = G @ powers - target - offset
    a = np.linspace(0, a_max, 6001)
    s = strehl_vs_step(residual, a)
    below = np.flatnonzero(s < strehl_min)
    if below.size == 0:
        a_limit = a_max
    else:
        j = below[0]
        a_limit = float(np.interp(strehl_min, [s[j], s[j - 1]], [a[j], a[j - 1]]))
    return {
        "powers_per_rad": powers,            # W per rad of lens strength
        "power_per_rad": powers.sum(),
        "rms_residual_per_rad": float(np.sqrt(np.mean(residual**2))),
        "a_limit": a_limit,                  # rad
        "f_min_um": lens_constant(x_emit) / a_limit if a_limit > 0 else np.inf,
    }


def least_power_lens(G, x_emit, tol=0.02):
    """Least total power per rad of lens strength that keeps every emitter within tol·a of the lens.

    The least-squares fit picks the drive with the smallest residual, which is
    not the cheapest: on a shared membrane many drives make nearly the same
    profile. This linear program asks for the cheapest drive whose phase error
    is at most tol per rad of strength at every emitter (0.02 keeps the Strehl
    above 0.99 up to a = π):
        minimise Σ P_k  subject to  |G·P - t - c| ≤ tol,  P ≥ 0,  c free.
    Returns W per rad, or NaN if the layout can't reach that accuracy.
    """
    from scipy.optimize import linprog

    target = unit_target(x_emit)
    n, k = G.shape
    ones = np.ones((n, 1))
    # Variables: P (k of them, ≥ 0) and c (free). Rows: G·P - c ≤ t + tol and -(G·P - c) ≤ -(t - tol).
    A = np.vstack([np.hstack([G, -ones]), np.hstack([-G, ones])])
    b = np.concatenate([target + tol, -(target - tol)])
    cost = np.concatenate([np.ones(k), [0.0]])
    res = linprog(cost, A_ub=A, b_ub=b, bounds=[(0, None)] * k + [(None, None)], method="highs")
    if res.status != 0:
        return np.nan
    # The tolerance band would let the LP make a slightly weaker lens (up to 2·tol weaker),
    # so divide by the strength it actually makes: the least-squares quadratic of its phase.
    phase = G @ res.x[:k]
    a_eff = np.polyfit(target, phase, 1)[0]
    return res.x[:k].sum() / a_eff


def isolated_power(a, x_emit, p_pi, wrapped):
    """Least total power (W) for one isolated heater per emitter.

    Unwrapped: φ_n = a·(1 - (x_n/X)²) ≥ 0, so the edge emitters need nothing.
    Wrapped (a Fresnel lens): each phase is only needed modulo 2π. As in
    conventional.py, the best common offset puts some emitter exactly at 0, so
    the candidates are c_m = -φ_m (mod 2π).
    """
    phi = a * (1 + unit_target(x_emit))
    if not wrapped:
        return p_pi / np.pi * phi.sum()
    offsets = np.append(np.mod(-phi, 2 * np.pi), 0.0)
    return min(p_pi / np.pi * np.mod(phi + c, 2 * np.pi).sum() for c in offsets)
