"""The conventional design: one heater per emitter, phases wrapped into [0, 2π).

Because exp(iφ) does not change when φ changes by 2π, each emitter only needs
its ramp phase modulo 2π. With one isolated heater per emitter (no thermal
crosstalk), heater n then needs power P_π·φ_n/π, never more than 2·P_π, and the
beam is perfect at every angle.

A common offset c shifts every phase, and so where the 2π resets fall. Raising
c lowers the total power until some emitter's phase wraps from just under 2π
to 0, so the least power always comes with one emitter at exactly 0. The
candidate offsets are therefore c_m = -Δφ·x_m/d (mod 2π), one per emitter.

The real layout is not isolated: at 2 µm pitch a heater warms its neighbours
to 75% of its own rise (Model C). with_crosstalk fits the target with the
real matrix G instead, which shows what the smooth heat does to the sharp 2π
resets.
"""
import numpy as np

from .constants import PITCH_UM
from .farfield import strehl
from .steering import STREHL_MIN, fit_powers


def wrapped_ramp(x_emit, dphi, offset, pitch=PITCH_UM):
    """Ramp phases Δφ·x_n/d + offset, reduced to [0, 2π)."""
    return np.mod(dphi * x_emit / pitch + offset, 2 * np.pi)


def candidate_offsets(x_emit, dphi, pitch=PITCH_UM):
    """The offsets that put one emitter exactly at phase 0, one per emitter."""
    return np.mod(-dphi * x_emit / pitch, 2 * np.pi)


def ideal_power(x_emit, dphi, p_pi, pitch=PITCH_UM):
    """Least total power (W) for isolated heaters, over all common offsets."""
    return min(p_pi / np.pi * wrapped_ramp(x_emit, dphi, c, pitch).sum()
               for c in candidate_offsets(x_emit, dphi, pitch))


def with_crosstalk(G, x_emit, dphi, u_target, pitch=PITCH_UM, strehl_min=STREHL_MIN):
    """Fit the ramp with the real G (heater over every emitter), wrapped or not.

    Candidates: the wrapped ramp at each candidate offset (each puts the 2π
    resets in different places), and the unwrapped ramp (no resets). Keeps the
    least power among those with Strehl ≥ strehl_min, or the best beam if none
    reaches it. Returns (total power W, Strehl, wrapped?).

    When a fit is exact the free offset makes the powers non-unique; at 10° the
    powers found here are within 2% of the true minimum.
    """
    targets = [(wrapped_ramp(x_emit, dphi, c, pitch), True) for c in candidate_offsets(x_emit, dphi, pitch)]
    targets.append((dphi * x_emit / pitch, False))
    best = None
    for target, wrapped in targets:
        powers, _ = fit_powers(G, target)
        s = strehl(G @ powers, x_emit, u_target)
        key = (s >= strehl_min, -powers.sum() if s >= strehl_min else s)
        if best is None or key > best[0]:
            best = (key, powers.sum(), s, wrapped)
    return best[1:]
