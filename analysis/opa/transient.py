"""Step responses from the R4 transient runs.

Each R4 run switches one heater on at t = 0 with 1 W/m and records the rise
T(x_n, t) - T0 at the 33 core centres. The system is linear and time-invariant,
so these step responses are the whole dynamic model: a heater driven with any
power waveform P(t) produces the convolution of dP/dt with its step response,
and several heaters add. Pre-emphasis design and beam-settling estimates are
built on this.

Array layout: S[i_t, k, n] is the rise (K per W/m) at time times[i_t], core n,
for a heater at x_h[k].
"""
import numpy as np

from .comsol_io import column_parameters, read_comsol
from .constants import DATA_DIR


class StepResponses:
    def __init__(self, filename, mirror=False, default_x_h=None, data_dir=DATA_DIR):
        """mirror=True adds negative heater positions using K(x; -x_h) = K(-x; x_h),
        valid for the symmetric (no-slot) structure. default_x_h labels a file
        from a single run, whose headers carry no x_h."""
        names, data = read_comsol(data_dir / filename)
        self.x = data[:, 0]  # core centres, µm
        params = [column_parameters(n) for n in names[2:]]
        t = np.array([p["t"] for p in params])
        x_h = np.array([p.get("x_h", default_x_h) for p in params], dtype=float)
        if np.isnan(x_h).any():
            raise ValueError(f"{filename}: no x_h in the headers; pass default_x_h")
        self.times = np.unique(t)
        self.x_h = np.unique(x_h)
        S = np.empty((len(self.times), len(self.x_h), len(self.x)))
        for col, (ti, xi) in enumerate(zip(t, x_h)):
            S[np.searchsorted(self.times, ti), np.searchsorted(self.x_h, xi)] = data[:, 2 + col]
        if mirror:
            if not np.allclose(self.x, -self.x[::-1]):
                raise ValueError("mirroring needs core positions symmetric about 0")
            neg = self.x_h > 0
            self.x_h = np.concatenate([-self.x_h[neg][::-1], self.x_h])
            S = np.concatenate([S[:, neg][:, ::-1, ::-1], S], axis=1)
        self.S = S

    def column(self, x_h):
        """Step response of one heater: array (n_times, n_cores)."""
        return self.S[:, int(np.flatnonzero(np.isclose(self.x_h, x_h))[0])]

    @property
    def final(self):
        """Rise at the last output time, shape (n_cores, n_heaters): the steady kernel if settled."""
        return self.S[-1].T


def settle_time(times, metric, threshold):
    """Time from which metric(t) stays ≥ threshold, interpolated in log t between outputs."""
    bad = np.flatnonzero(metric < threshold)
    if bad.size == 0:
        return times[0]
    j = bad[-1]
    if j + 1 >= len(times):
        return np.nan
    w = (threshold - metric[j]) / (metric[j + 1] - metric[j])
    return float(np.exp(np.log(times[j]) + w * (np.log(times[j + 1]) - np.log(times[j]))))


def crossing_time(times, y, fraction):
    """First time y(t) reaches fraction × its final value (log-interpolated in t)."""
    target = fraction * y[-1]
    above = np.flatnonzero(y >= target) if y[-1] > 0 else np.flatnonzero(y <= target)
    if above.size == 0:
        return np.nan
    j = above[0]
    if j == 0:
        return times[0]
    w = (target - y[j - 1]) / (y[j] - y[j - 1])
    return float(np.exp(np.log(times[j - 1]) + w * (np.log(times[j]) - np.log(times[j - 1]))))
