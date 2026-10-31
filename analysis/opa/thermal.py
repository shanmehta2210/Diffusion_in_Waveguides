"""Heater-to-phase model: thermal spreading kernels and the matrix G.

The 2D thermal model (Model C) gives, per 1 W/m of heater power, the
temperature rise K(x; x_h) at core height for a heater centred at x_h. The
model is per unit length, so a heater of length L and power P raises the core
by K·P/L, and over that length the phase changes by
    (2π/λ)(dn_eff/dT)(K·P/L)·L = PHASE_COEFF · K · P.
The heater length cancels: phase depends only on total power. Stacking this
for every emitter n and heater k gives φ = G·P with G[n, k] = PHASE_COEFF·K(x_n; x_k).
"""
import numpy as np

from .comsol_io import column_parameters, read_comsol
from .constants import DATA_DIR, N_EMITTERS, PHASE_COEFF, PITCH_UM

# The brief names kernel_checks_1.csv as the final standard kernel, but that
# file is missing. kernel_checks.csv reproduces every value the brief quotes
# (K(0) = 0.1959 K per W/m), so we use it; swap the name if the final file turns up.
STANDARD_KERNEL_FILE = "kernel_checks.csv"
UNDERCUT_KERNEL_FILE = "kernel_undercut.csv"


class StandardKernel:
    """Silicon substrate under the buried oxide.

    The structure looks the same from every heater position, so the rise
    depends only on the distance x - x_h.
    """

    name = "standard"

    def __init__(self, data_dir=DATA_DIR):
        _, data = read_comsol(data_dir / STANDARD_KERNEL_FILE)
        self.x = data[:, 0]
        self.k = data[:, 2]  # first data column: the base model (the others are variants)

    def __call__(self, x, x_h):
        """Rise in K per (W/m) at position(s) x for one heater at x_h."""
        return np.interp(np.subtract(x, x_h), self.x, self.k)


class UndercutKernel:
    """Air cavity (80 µm wide, 15 µm deep) under the array.

    The cavity has edges, so the profile changes shape as the heater moves:
    COMSOL gives one column per heater position x_h = 0, 1, ..., 32 µm, and we
    interpolate linearly between neighbouring columns.
    """

    name = "undercut"

    def __init__(self, data_dir=DATA_DIR):
        names, data = read_comsol(data_dir / UNDERCUT_KERNEL_FILE)
        self.x = data[:, 0]
        self.x_h = np.array([column_parameters(n)["x_h"] for n in names[2:]])
        assert np.array_equal(self.x_h, np.arange(len(self.x_h))), "expected x_h = 0, 1, 2, ... µm"
        self.k = data[:, 2:]  # shape (len(x), 33)
        self.x_h_max = self.x_h[-1]

    def __call__(self, x, x_h):
        """Rise in K per (W/m) at position(s) x for one heater at x_h."""
        if abs(x_h) > self.x_h_max:
            raise ValueError(f"heater at {x_h} µm is outside the modelled range ±{self.x_h_max:g} µm")
        x = np.asarray(x, dtype=float)
        # The structure is mirror-symmetric about x = 0, so K(x; -x_h) = K(-x; x_h).
        if x_h < 0:
            x, x_h = -x, -x_h
        j = min(int(x_h), len(self.x_h) - 2)  # column just below x_h
        w = x_h - self.x_h[j]
        profile = (1 - w) * self.k[:, j] + w * self.k[:, j + 1]
        return np.interp(x, self.x, profile)


def emitter_positions(n=N_EMITTERS, pitch=PITCH_UM):
    """Emitter centres, symmetric about x = 0.

    For an odd count one emitter sits at x = 0: 33 at 2 µm gives 0, ±2, ..., ±32 µm,
    the core positions of the thermal model.
    """
    return (np.arange(n) - (n - 1) / 2) * pitch


def phase_matrix(kernel, x_emit, x_heat):
    """G[n, k]: phase in radians at emitter n per watt in heater k."""
    return PHASE_COEFF * np.column_stack([kernel(x_emit, xh) for xh in x_heat])
