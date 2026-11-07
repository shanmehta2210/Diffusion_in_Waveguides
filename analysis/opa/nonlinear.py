"""Nonlinearity at high temperature: when the linear thermal model stops being exact.

Two effects bend the phase away from "phase ∝ heater power":

1. Heating (measured, Model F). The oxide conducts better when hot, so the
   real rise falls below the linear one. For steady conduction with a
   conductivity k = k0(1 + βΔT), the Kirchhoff transform θ = ∫ k/k0 dT turns the
   heat equation back into the linear one. So the linear model really predicts
   θ, and the true rise follows point by point:
       ΔT = θ - a·θ²,   a = β/2.
   On the undercut membrane this reproduces every Model F profile to within
   0.25% of the peak, and two heaters together to 0.3% (run_nonlinearity.py).
   Because it is pointwise, it applies to any multi-heater drive: superpose the
   linear rises first, then transform.

2. Optics (bracketed, not measured here). Silicon's dn/dT itself grows with
   temperature, at a relative rate b = (1/(dn/dT))·d(dn/dT)/dT. The index
   change is then Δn = (dn/dT)(ΔT + b·ΔT²/2). Our COMSOL models used a constant
   dn/dT (Model A's n_eff is linear in T to 3.5e-6 over 100 K), so b has to come
   from the literature. The source could not be verified, so b is treated as
   a bracket.

To first order the phase goes as φ ≈ φ_lin·(1 + (b/2 - a)·θ): the two effects
pull in opposite directions.
"""
import numpy as np
from scipy.optimize import brentq, least_squares

from .comsol_io import column_parameters, read_comsol
from .constants import DATA_DIR, PHASE_COEFF, WAVELENGTH_UM


def model_f_runs(data_dir=DATA_DIR):
    """Model F profiles: x (µm), the run parameters, and the rises (K) at core height."""
    names, data = read_comsol(data_dir / "nonlinear_line.csv")
    return data[:, 0], [column_parameters(n) for n in names[2:]], data[:, 2:]


def true_rise(theta, a):
    """Real temperature rise (K) from the linear model's rise θ (K): the Kirchhoff law, first order."""
    return theta - a * theta**2


def kirchhoff_theta(dT, a):
    """Inverse of true_rise: the linear-model rise that gives a real rise dT."""
    return (1 - np.sqrt(1 - 4 * a * dT)) / (2 * a)


def fit_heating(x, runs, rises, undercut, x_max=32.0):
    """Best a (1/K) for the Kirchhoff law over the emitter span, from Model F's single-heater runs."""
    lin = rises[:, [i for i, p in enumerate(runs)
                    if int(p["uc"]) == undercut and int(p["nl"]) == 0][0]]
    inside = np.abs(x) <= x_max + 1e-6
    th, dt = [], []
    for i, p in enumerate(runs):
        if int(p["uc"]) == undercut and int(p["nl"]) == 1 and p["Q_l2"] == 0:
            th.append(p["Q_l"] * lin[inside])
            dt.append(rises[inside, i])
    th, dt = np.concatenate(th), np.concatenate(dt)
    return float(np.sum((th - dt) * th**2) / np.sum(th**4))   # least squares for θ - dt = a·θ²


def true_phase(phi_lin, heater_length_m, a, b):
    """Real phase (rad) for linear-model phases phi_lin, with heating law a and dn/dT growth b."""
    theta = phi_lin / (PHASE_COEFF * heater_length_m)          # linear-model rise, K
    dT = true_rise(theta, a)
    return PHASE_COEFF * heater_length_m * (dT + 0.5 * b * dT**2)


def strehl(residual):
    return np.abs(np.exp(1j * np.atleast_2d(residual)).mean(axis=1)) ** 2


def pointing_deg(phase, x_emit):
    """Beam angle from the least-squares slope of the phase across the array."""
    slope = np.polyfit(x_emit, phase, 1)[0]
    return float(np.degrees(np.arcsin(slope * WAVELENGTH_UM / (2 * np.pi))))


def shape_strehl(phase, x_emit):
    """Strehl toward wherever the beam actually points: the error after removing the best tilt."""
    fit = np.polyval(np.polyfit(x_emit, phase, 1), x_emit)
    return float(strehl(phase - fit)[0])


def gain_corrected_drive(G, P_lin, x_emit, theta_target, heater_length_m, a, b):
    """Scale the linear drive by one number s so the real beam points at theta_target.

    This is a one-dimensional calibration: one gain per angle (a lookup table,
    or a smooth curve), on top of the single linear fit.
    """
    def miss(s):
        return pointing_deg(true_phase(G @ (s * P_lin), heater_length_m, a, b), x_emit) - theta_target
    s = brentq(miss, 0.7, 1.4)
    return s * P_lin, s


def refit_drive(G, target, P_start, heater_length_m, a, b):
    """Full recalibration: the non-negative drive whose real phase best matches the target (offset free)."""
    K = G.shape[1]

    def residual(v):
        return true_phase(G @ v[:K], heater_length_m, a, b) - target - v[K]

    c0 = float(np.mean(G @ P_start - target))
    sol = least_squares(residual, np.append(P_start, c0),
                        bounds=(np.append(np.zeros(K), -np.inf), np.full(K + 1, np.inf)))
    return sol.x[:K]
