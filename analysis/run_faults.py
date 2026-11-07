"""Fault tolerance: what happens to the steering when something goes wrong after the chip is built.

The four designs used throughout (vacuum 1 heater, slot 3 heaters, undercut 3 heaters,
dual-use 5 heaters, all steering toward +θ), and for comparison the conventional array: a heater
under every emitter on the standard stack, isolated, phases wrapped into [0, 2π).

  1. A dead heater (open circuit). Before recalibration the other heaters keep their drive;
     after it, the remaining heaters are refitted.
  2. Drive noise. Each heater's power wanders by a fraction σ = 1%, either independently or all
     together (supply noise). A heater that makes a ramp turns its error into pointing error; a
     heater per emitter turns it into random phase error, which costs a little Strehl and mostly
     averages out of the pointing.
  3. Drift after calibration, with the drive held fixed: the package pressure rising (vacuum
     design), the air pressure changing (slot in air), the ambient temperature, and radiation,
     which the COMSOL vacuum run left out. These use the 1D membrane model with gas conduction
     (opa/membrane1d.py, calibrated to COMSOL R2/R3b, checked against R3a).

Pointing is read as the least-squares tilt of the phase across the array (the far-field peak
follows it for small errors); Strehl is |mean exp(i r)|² of what is left after the tilt.

Writes results/faults.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np

from opa import membrane1d
from opa.constants import PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.steering import steer_one_direction
from opa.thermal import SweptKernel, UndercutKernel, emitter_positions, phase_matrix

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
rng = np.random.default_rng(7)
xe = emitter_positions()
THETA = 10.0                                        # operating angle for the fault cases, deg
rows = []


def dphi(theta_deg):
    """Phase step per emitter for a beam at theta_deg."""
    return 2 * np.pi * PITCH_UM * np.sin(np.radians(theta_deg)) / WAVELENGTH_UM


def beam(phases):
    """Pointing angle (deg) and Strehl of a phase profile across the array.

    The least-squares line through the phases is the tilt that points the beam; whatever is
    left after removing it is the error that lowers the peak.
    """
    x = xe - xe.mean()
    slope = np.dot(x, phases - phases.mean()) / np.dot(x, x)      # rad per µm
    rest = phases - phases.mean() - slope * x
    theta = np.degrees(np.arcsin(slope * WAVELENGTH_UM / (2 * np.pi)))
    return theta, np.abs(np.exp(1j * rest).mean()) ** 2


def record(section, design, case, **values):
    for quantity, value in values.items():
        rows.append([section, design, case, quantity, value])


class Design:
    def __init__(self, name, G, heaters):
        self.name, self.G, self.heaters = name, G, list(heaters)
        fit = steer_one_direction(G, xe, +1)
        self.range, self.mw_per_deg = fit["range_deg"], fit["mw_per_deg"]
        self.P = dphi(THETA) * fit["powers_per_rad"]                # hold drive at THETA, W

    def phases(self, P):
        return self.G @ P


vacuum = SweptKernel("kernel_vacuum.csv", "vacuum", False)
slot = SweptKernel("kernel_slot.csv", "slot", False, 67)
far = SweptKernel("kernel_undercut_far.csv", "far", True, 41)
undercut = UndercutKernel()
designs = [Design(name, phase_matrix(k, xe, h), h) for name, k, h in (
    ("vacuum, 1 heater", vacuum, [34]),
    ("slot, 3 heaters", slot, [2, 14, 32]),
    ("undercut, 3 heaters", undercut, [-8, 10, 32]),
    ("dual-use, 5 heaters", far, [-33, -16, 0, 16, 33]))]

# The conventional array: emitter n needs the wrapped phase of the ramp. Its error is all that
# matters, so it is added to the unwrapped ramp, which points the beam at THETA exactly.
ramp = dphi(THETA) * xe / PITCH_UM
wrapped = np.mod(ramp, 2 * np.pi)

# ---------------------------------------------------------------- 0. check
print("0. Check: the clean drives point where they should")
for d in designs:
    th, s = beam(d.phases(d.P))
    print(f"  {d.name:22} {th:6.3f}° (target {THETA:g}°), Strehl {s:.4f}; range {d.range:.1f}° at {d.mw_per_deg:.3f} mW/deg")
th, s = beam(ramp)
print(f"  {'conventional':22} {th:6.3f}°, Strehl {s:.4f}")

# ---------------------------------------------------------------- 1. dead heater
print(f"\n1. One heater dies (open circuit), beam at {THETA:g}°")
print("   before recalibration: the others keep their drive;  after: the rest are refitted")
for d in designs:
    for j, xh in enumerate(d.heaters):
        P = d.P.copy()
        P[j] = 0.0
        th, s = beam(d.phases(P))
        rest = [k for k in range(len(d.heaters)) if k != j]
        if rest:
            fit = steer_one_direction(d.G[:, rest], xe, +1)
            after = (fit["range_deg"], fit["mw_per_deg"])
        else:
            after = (0.0, np.nan)
        record("dead heater", d.name, f"heater at {xh:g} µm", pointing_deg=th, strehl=s,
               range_after_deg=after[0], mw_per_deg_after=after[1])
        print(f"  {d.name:22} heater at {xh:4g} µm: beam goes to {th:5.2f}°, Strehl {s:.2f};"
              f"  refitted: range {after[0]:4.1f}° at {after[1]:.3f} mW/deg (was {d.range:.1f}° at {d.mw_per_deg:.3f})")

# A spare for the one-heater vacuum design: a second heater beside the first, wired separately.
fit = steer_one_direction(phase_matrix(vacuum, xe, [32]), xe, +1)
record("dead heater", "vacuum, 1 heater", "spare heater at 32 µm", range_after_deg=fit["range_deg"],
       mw_per_deg_after=fit["mw_per_deg"])
print(f"  vacuum spare at 32 µm alone: range {fit['range_deg']:.1f}° at {fit['mw_per_deg']:.3f} mW/deg")
# The same idea for the slot design: a spare at 34 µm, just outside the edge heater at 32 µm.
fit = steer_one_direction(phase_matrix(slot, xe, [2, 14, 34]), xe, +1)
record("dead heater", "slot, 3 heaters", "spare at 34 µm replacing 32 µm", range_after_deg=fit["range_deg"],
       mw_per_deg_after=fit["mw_per_deg"])
print(f"  slot with a spare at 34 µm in place of 32 µm: range {fit['range_deg']:.1f}° at {fit['mw_per_deg']:.3f} mW/deg")

# Conventional: a dead heater leaves its emitter at the unheated phase, an error of -wrapped[n].
worst, mean_s = 1.0, []
for n in range(len(xe)):
    err = np.zeros_like(xe)
    err[n] = -wrapped[n]
    th, s = beam(ramp + err)
    mean_s.append(s)
    worst = min(worst, s)
record("dead heater", "conventional, 33 heaters", "any one heater", strehl_mean=np.mean(mean_s),
       strehl_worst=worst)
print(f"  conventional, any one of 33: Strehl {np.mean(mean_s):.3f} on average, {worst:.3f} worst; pointing unchanged."
      "\n    Recalibrated: a shared offset makes the dead emitter's target phase 0, so one dead heater costs nothing.")

# ---------------------------------------------------------------- 2. drive noise
SIGMA, TRIALS = 0.01, 4000
print(f"\n2. Drive noise: each heater's power off by {SIGMA:.0%} rms, beam at {THETA:g}° ({TRIALS} trials)")
for kind in ("independent", "common"):
    for d in designs + ["conventional"]:
        th_err, strehls = [], []
        for _ in range(TRIALS):
            if d == "conventional":
                eps = rng.standard_normal(len(xe)) if kind == "independent" else np.full(len(xe), rng.standard_normal())
                th, s = beam(ramp + SIGMA * eps * wrapped)
            else:
                eps = rng.standard_normal(len(d.P)) if kind == "independent" else np.full(len(d.P), rng.standard_normal())
                th, s = beam(d.phases(d.P * (1 + SIGMA * eps)))
            th_err.append(th - THETA)
            strehls.append(s)
        name = d if isinstance(d, str) else d.name
        rms = np.sqrt(np.mean(np.square(np.array(th_err) - np.mean(th_err))))
        record("drive noise", name, kind, pointing_rms_deg=rms, strehl_mean=np.mean(strehls))
        print(f"  {kind:11} {name:24} pointing jitter {rms:.4f}° rms, Strehl {np.mean(strehls):.4f}")

# ---------------------------------------------------------------- 3. drift after calibration
m1 = membrane1d
print("\n3. Drift after calibration, drive held fixed (1D membrane model with gas, est.)")


def drift(heaters, D_um, w_um, p_cal, cases, label):
    """Calibrate the drive for THETA at p_cal, then report the beam for each (case, model kwargs)."""
    G0 = PHASE_COEFF * m1.model(m1.A_H, m1.A_S, D_um, w_um, heaters, p_cal)
    P = dphi(THETA) * steer_one_direction(G0, xe, +1)["powers_per_rad"]
    th0, _ = beam(G0 @ P)
    for case, kw in cases:
        G = PHASE_COEFF * m1.model(m1.A_H, m1.A_S, D_um, w_um, heaters, **kw)
        th, s = beam(G @ P)
        record("drift", label, case, pointing_error_deg=th - th0, strehl=s)
        print(f"  {label:22} {case:30} pointing {th - th0:+7.3f}°, Strehl {s:.3f}")


# 3a. Vacuum design: calibrated at 1 Pa, then the package leaks.
drift([34], 15, 2, 1.0, [(f"pressure {p:g} Pa", {"pressure": p}) for p in (10, 30, 100, 300, 1000, 1e4, 101325)],
      "vacuum, 1 heater")
# 3b. Slot in air: calibrated at 1 atm; weather (±5%) and altitude (2.5 km, about 0.75 atm).
drift([2, 14, 32], 15, 2, 101325, [(f"pressure {f:g} atm", {"pressure": 101325 * f}) for f in (0.75, 0.95, 1.05)],
      "slot in air, 3 heaters")

# 3c. Ambient temperature. Two parts:
#  - the membrane and silicon: conductivity rises as k0(1 + 2a·ΔT) and dn/dT as (1 + b·ΔT), so the
#    phase per watt changes by (b - 2a) per K. a = 3.0e-4/K (our Model F); b is bracketed 0-1.4e-3/K
#    (Cocorullo 1999, unverified), as in notes/nonlinearity.md.
#  - in air, the gas also conducts better when warm: k_air ∝ T^0.8, and the mean free path ∝ T.
A_KIRCHHOFF = 3.0e-4
for b in (0.0, 0.7e-3, 1.4e-3):
    per_10K = THETA * (b - 2 * A_KIRCHHOFF) * 10
    record("drift", "all designs", f"ambient +10 K, b = {b:.1e}", pointing_error_deg=per_10K)
    print(f"  {'all designs':22} ambient +10 K, b = {b * 1e3:.1f}e-3/K   pointing {per_10K:+7.3f}° (membrane k and dn/dT)")
dT = 10.0
drift([2, 14, 32], 15, 2, 101325,
      [("ambient +10 K, gas only", {"pressure": 101325, "k_bulk": m1.K_AIR * (1 + dT / 300) ** 0.8,
                                     "lambda_atm": m1.LAMBDA_ATM * (1 + dT / 300)})],
      "slot in air, 3 heaters")

# 3d. Radiation, left out of the COMSOL vacuum run. Both faces radiate with h = 4εσT³ each, at
# most ε = 1. This is a loss to the surroundings, so it acts like gas under the membrane.
SIGMA_SB, T0 = 5.670e-8, 300.0
G_vac = PHASE_COEFF * m1.profiles(0.0, 0.0, [34])
P = dphi(THETA) * steer_one_direction(G_vac, xe, +1)["powers_per_rad"]
th0, _ = beam(G_vac @ P)
for eps in (0.5, 1.0):
    h_rad = 2 * 4 * eps * SIGMA_SB * T0**3
    G = PHASE_COEFF * m1.profiles(h_rad, 0.0, [34])
    th, s = beam(G @ P)
    fit = steer_one_direction(G, xe, +1)
    record("drift", "vacuum, 1 heater", f"radiation, emissivity {eps:g}", pointing_error_deg=th - th0,
           strehl=s, mw_per_deg_after=fit["mw_per_deg"])
    print(f"  {'vacuum, 1 heater':22} radiation, ε = {eps:g} (h = {h_rad:.1f} W/m²K)   pointing {th - th0:+7.3f}°,"
          f" Strehl {s:.4f}; refitted {fit['mw_per_deg']:.3f} mW/deg against {steer_one_direction(G_vac, xe, +1)['mw_per_deg']:.3f}")

with open(RESULTS / "faults.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["section", "design", "case", "quantity", "value"])
    w.writerows(rows)
print(f"\nWrote {len(rows)} rows to results/faults.csv")
