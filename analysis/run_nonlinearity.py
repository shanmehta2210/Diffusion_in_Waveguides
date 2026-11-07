"""Analysis 8: nonlinearity at high temperature, on the final designs.

Parts:
  1. checks: the Kirchhoff law against every Model F profile, two-heater
     superposition, the brief's "1% at 30 K, 3% at 80 K", and that our COMSOL
     index model had a constant dn/dT;
  2. the final designs with 1 mm heaters: hottest core, and what the
     nonlinearity does to the beam with three levels of calibration:
       linear     the single linear fit, as used everywhere so far;
       gain       the same drive shape, rescaled by one number per angle so the
                  beam points right (a one-dimensional lookup);
       refit      the drive re-optimised against the real phase;
  3. heater length: 0.5, 1 and 2 mm at 10°;
  4. heater temperatures in the slot design, steady and during pre-emphasis.

The dn/dT growth rate b is bracketed (0, 0.7e-3, 1.4e-3 per K): the literature
value could not be verified (opa/nonlinear.py). Writes results/nonlinear_designs.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np

from opa.comsol_io import read_comsol
from opa.constants import DATA_DIR, PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.nonlinear import (fit_heating, gain_corrected_drive, kirchhoff_theta, model_f_runs,
                           pointing_deg, refit_drive, shape_strehl, strehl, true_phase, true_rise)
from opa.steering import steer_one_direction
from opa.thermal import SweptKernel, UndercutKernel, emitter_positions, phase_matrix

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
xe = emitter_positions()
B_BRACKET = (0.0, 0.7e-3, 1.4e-3)   # relative growth of dn/dT per K; 1.4e-3 is the unverified literature value


def dphi(theta_deg):
    return 2 * np.pi * PITCH_UM * np.sin(np.radians(theta_deg)) / WAVELENGTH_UM


def ramp(theta_deg):
    return dphi(theta_deg) * xe / PITCH_UM


# ---------------------------------------------------------------- 1. checks
print("1. Checks against Model F")
x, runs, rises = model_f_runs()
inside = np.abs(x) <= 32 + 1e-6
A = {}
for uc, label in ((1, "undercut"), (0, "standard")):
    a = fit_heating(x, runs, rises, uc)
    A[uc] = a
    lin = rises[:, [i for i, p in enumerate(runs) if int(p["uc"]) == uc and int(p["nl"]) == 0][0]]
    print(f"  {label}: a = {a:.3e} per K (brief: about 3.3e-4)")
    for i, p in enumerate(runs):
        if int(p["uc"]) == uc and int(p["nl"]) == 1 and p["Q_l2"] == 0:
            th = p["Q_l"] * lin[inside]
            real = rises[inside, i]
            print(f"    Q = {p['Q_l']:5g} W/m, peak {real.max():6.1f} K: worst error, % of peak: linear model "
                  f"{100 * np.max(np.abs(th - real)) / real.max():5.2f};  Kirchhoff law "
                  f"{100 * np.max(np.abs(true_rise(th, a) - real)) / real.max():5.2f}")
for uc, (i1, i2, i12), label in ((1, (17, 18, 19), "undercut"), (0, (14, 15, 16), "standard")):
    t1, t2, t12 = rises[inside, i1], rises[inside, i2], rises[inside, i12]
    pred = true_rise(kirchhoff_theta(t1, A[uc]) + kirchhoff_theta(t2, A[uc]), A[uc])
    print(f"  {label}, two heaters together (peak {t12.max():.1f} K): plain sum of the two single runs is off by "
          f"{100 * np.max(np.abs(t1 + t2 - t12)) / t12.max():.2f}% of peak; Kirchhoff law "
          f"{100 * np.max(np.abs(pred - t12)) / t12.max():.2f}%")
for dT in (30, 80):
    print(f"  Heating shortfall at a linear rise of {dT} K: {100 * A[1] * dT:.1f}% (brief: {1 if dT == 30 else 3}%)")
_, neff = read_comsol(DATA_DIR / "neff_vs_T.csv")
T, n = neff[:, 0] - 293.15, neff[:, 1]
c2, c1, _ = np.polyfit(T, n, 2)
print(f"  Model A index: dn_eff/dT = {c1:.4e} per K, curvature 2c2/c1 = {2 * c2 / c1:.1e} per K "
      f"(constant dn/dT in our COMSOL material; real Si grows at about 1e-3 per K, unverified)")

# ---------------------------------------------------------------- 2. designs
slot = SweptKernel("kernel_slot.csv", "slot", False, 67)
vacuum = SweptKernel("kernel_vacuum.csv", "vacuum", False)
designs = [("vacuum, 1 heater", vacuum, [34]), ("slot, 3 heaters", slot, [2, 14, 32]),
           ("undercut, 3 heaters", UndercutKernel(), [-8, 10, 32])]
a = A[1]   # membrane designs: use the undercut fit
L1 = 1e-3  # 1 mm heaters
rows = []
print(f"\n2. Final designs, 1 mm heaters, a = {a:.2e} per K. For each b (×1e-3 per K):")
print("   pointing error with the linear fit;  Strehl toward the target: linear fit / gain-corrected / refit;"
      "  shape Strehl (beam quality wherever it points)")
for name, kernel, heaters in designs:
    G = phase_matrix(kernel, xe, heaters)
    p_unit = steer_one_direction(G, xe, +1)["powers_per_rad"]
    for theta in (2, 5, 10, 15, 20):
        P = dphi(theta) * p_unit
        phi_lin = G @ P
        s_lin0 = strehl(phi_lin - ramp(theta))[0]
        hot_lin = phi_lin.max() / (PHASE_COEFF * L1)
        hot = true_rise(hot_lin, a)
        parts = []
        for b in B_BRACKET:
            phi = true_phase(phi_lin, L1, a, b)
            err = pointing_deg(phi, xe) - theta
            s_lin = strehl(phi - ramp(theta))[0]
            Pg, s = gain_corrected_drive(G, P, xe, theta, L1, a, b)
            s_gain = strehl(true_phase(G @ Pg, L1, a, b) - ramp(theta))[0]
            Pr = refit_drive(G, ramp(theta), Pg, L1, a, b)
            s_ref = strehl(true_phase(G @ Pr, L1, a, b) - ramp(theta))[0]
            rows.append([name, theta, 1.0, b, hot, s_lin0, err, s_lin, s, s_gain, s_ref, shape_strehl(phi, xe)])
            parts.append(f"b={b * 1e3:.1f}: {err:+.2f}°  {s_lin:.3f}/{s_gain:.3f}/{s_ref:.3f}")
        flag = "  > 80 K" if hot > 80 else ("  > 30 K" if hot > 30 else "")
        print(f"  {name:20} {theta:2d}°  hottest core {hot:5.1f} K{flag:8}  linear-model Strehl {s_lin0:.3f} | "
              + " | ".join(parts))

# ---------------------------------------------------------------- 3. heater length
print("\n3. Heater length, slot design at 10° (b = 1.4e-3, the worst of the bracket)")
name, kernel, heaters = designs[1]
G = phase_matrix(kernel, xe, heaters)
P = dphi(10) * steer_one_direction(G, xe, +1)["powers_per_rad"]
for L in (0.5e-3, 1e-3, 2e-3):
    for b in B_BRACKET:
        phi = true_phase(G @ P, L, a, b)
        hot = true_rise((G @ P).max() / (PHASE_COEFF * L), a)
        Pg, s = gain_corrected_drive(G, P, xe, 10, L, a, b)
        s_gain = strehl(true_phase(G @ Pg, L, a, b) - ramp(10))[0]
        rows.append([name, 10, L * 1e3, b, hot, strehl(G @ P - ramp(10))[0], pointing_deg(phi, xe) - 10,
                     strehl(phi - ramp(10))[0], s, s_gain, np.nan, shape_strehl(phi, xe)])
        if b == B_BRACKET[-1] or b == 0:
            print(f"  L = {L * 1e3:.1f} mm, b = {b * 1e3:.1f}e-3: hottest core {hot:6.1f} K, pointing error "
                  f"{pointing_deg(phi, xe) - 10:+.2f}°, Strehl linear fit {strehl(phi - ramp(10))[0]:.3f}, "
                  f"gain-corrected {s_gain:.3f} (gain {s:.3f})")

with open(RESULTS / "nonlinear_designs.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["design", "angle_deg", "heater_mm", "b_per_K", "hottest_core_K", "strehl_linear_model",
                "pointing_error_deg_linear_fit", "strehl_linear_fit", "gain", "strehl_gain_corrected",
                "strehl_refit", "shape_strehl_linear_fit"])
    w.writerows(rows)

# ---------------------------------------------------------------- 4. heater temperatures
print("\n4. Heater temperatures, slot design, 1 mm heaters (heater_Tmax_slot.csv: hottest point in each heater)")
_, tab = read_comsol(DATA_DIR / "heater_Tmax_slot.csv")
tmax = {round(float(r[2]), 3): float(r[3]) for r in tab}     # K per W/m, heater alone


def heater_rises(P_w, L=1e-3):
    """Hottest rise in each heater: its own peak plus the others' rise at its position (linear), then Kirchhoff."""
    q = P_w / L                                                # W/m
    out = []
    for k, xh in enumerate(heaters):
        own = tmax[float(xh)] * q[k]
        others = sum(float(np.interp(xh, kernel.x, kernel.column(xj))) * q[j]
                     for j, xj in enumerate(heaters) if j != k)
        out.append(true_rise(own + others, a))
    return np.array(out)


print(f"  steady 10°: drive {', '.join(f'{p * 1e3:.2f}' for p in P)} mW -> heater rises "
      f"{', '.join(f'{t:.0f}' for t in heater_rises(P))} K (hottest core {true_rise((G @ P).max() / (PHASE_COEFF * L1), a):.0f} K)")
for k, xh in enumerate(heaters):
    burst = np.zeros(len(heaters))
    burst[k] = 3 * P.sum()                                     # the m = 3 pre-emphasis burst, all in one heater
    print(f"  pre-emphasis burst of {burst[k] * 1e3:.1f} mW in the heater at {xh} µm: steady-state bound "
          f"{heater_rises(burst)[k]:.0f} K (the burst lasts about 0.1-0.2 ms, so the real peak is lower)")
