"""Where the slot design loses heat in air, and how much gas pressure has to drop to recover the vacuum result.

Uses the 1D membrane model with gas conduction in opa/membrane1d.py (its docstring gives the
equations). This script fits the model's two gas-path factors a_h (cavity) and a_s (slot) to
the COMSOL slot runs R2 (2 µm slot, air) and R3b (4 µm slot, air), checks the vacuum result
against R3a, splits the loss between the two gas paths, and sweeps the gas pressure. The fitted
factors are copied into opa/membrane1d.py (A_H, A_S) for other scripts; refit here if the
COMSOL runs change.

Writes results/gas_loss.csv.
"""
import csv
import itertools
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from opa import membrane1d
from opa.constants import PHASE_COEFF
from opa.membrane1d import T_MEM
from opa.steering import steer_one_direction
from opa.thermal import emitter_positions

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
xe = emitter_positions()
GRID = np.arange(-32, 35, 2.0)  # heater positions searched (as the R3 sweeps)


def profiles(h, G_s):
    """Temperature (K per W/m) at the 33 cores for a heater at each GRID position."""
    return membrane1d.profiles(h, G_s, GRID)                # (33 cores, heaters)


def best(Tmat, k):
    """Best one-way layout of k heaters: widest range at Strehl 0.8, then least power."""
    G_all = PHASE_COEFF * Tmat
    top = None
    for combo in itertools.combinations(range(len(GRID)), k):
        if GRID[combo[-1]] < 32:          # every good one-way layout keeps a heater at the slot end
            continue
        r = steer_one_direction(G_all[:, combo], xe, +1)
        key = (round(r["range_deg"], 1), -r["mw_per_deg"])
        if top is None or key > top[0]:
            top = (key, combo, r)
    return [float(GRID[i]) for i in top[1]], top[2]


def one_heater_34(Tmat):
    return steer_one_direction(PHASE_COEFF * Tmat[:, [list(GRID).index(34.0)]], xe, +1)


def model(a_h, a_s, D_um, w_um, pressure=101325.0):
    return membrane1d.model(a_h, a_s, D_um, w_um, GRID, pressure)


# ---------------------------------------------------------------- calibration
COMSOL = {  # (slot width µm, gas) -> one heater at 34: (range°, mW/deg); best 3: mW/deg
    (2.0, "air"): ((3.98, 1.298), 1.383),
    (4.0, "air"): ((3.98, 1.128), 1.224),
}


def misfit(p):
    a_h, a_s = np.exp(p)
    out = []
    for (w, _), ((_, p1), p3) in COMSOL.items():
        T = model(a_h, a_s, 15.0, w)
        out.append(one_heater_34(T)["mw_per_deg"] / p1 - 1)
        out.append(best(T, 3)[1]["mw_per_deg"] / p3 - 1)
    return out


print("Calibrating the gas-loss factors to COMSOL R2 and R3b (this searches heater layouts, about a minute)")
fit = least_squares(misfit, np.log([1.0, 1.0]), diff_step=0.05)
a_h, a_s = np.exp(fit.x)
print(f"  a_h = {a_h:.2f} (cavity), a_s = {a_s:.2f} (slot); residuals {np.round(np.array(fit.fun) * 100, 1)} %")
for (w, gas), ((rng, p1), p3) in COMSOL.items():
    T = model(a_h, a_s, 15.0, w)
    r1, (lay, r3) = one_heater_34(T), best(T, 3)
    print(f"  {w:g} µm slot, air: one heater {r1['range_deg']:.1f}° at {r1['mw_per_deg']:.3f} (COMSOL {rng}° at {p1});"
          f" 3 heaters {lay} {r3['range_deg']:.1f}° at {r3['mw_per_deg']:.3f} (COMSOL {p3})")
T = profiles(0.0, 0.0)
r = one_heater_34(T)
print(f"  vacuum check (no gas at all): one heater {r['range_deg']:.1f}° at {r['mw_per_deg']:.3f} mW/deg "
      f"(COMSOL R3a: 22.8° at 0.532; conduction limit 0.506)")

# ---------------------------------------------------------------- loss split at 1 atm
print("\nWhere the power goes at 1 atm (2 µm slot, 15 µm cavity), best 3-heater layout:")
for label, ah, as_ in (("both gas paths", a_h, a_s), ("slot gas only", 0.0, a_s), ("cavity gas only", a_h, 0.0),
                       ("no gas", 0.0, 0.0)):
    T = model(max(ah, 1e-12), max(as_, 1e-12), 15.0, 2.0)
    lay, r3 = best(T, 3)
    print(f"  {label:16} {r3['range_deg']:5.1f}° at {r3['mw_per_deg']:.3f} mW/deg  (heaters {lay})")

# ---------------------------------------------------------------- pressure sweeps
rows = []
print("\nBest 3-heater layout against gas pressure (air), mW/deg at the full range; vacuum floor as above")
pressures = [101325, 1e4, 3e3, 1e3, 300, 100, 30, 10, 1]
for D, w in ((15, 2), (15, 10), (40, 10), (80, 10)):
    line = []
    for p in pressures:
        T = model(a_h, a_s, D, w, p)
        lay, r3 = best(T, 3)
        rows.append([D, w, p, r3["range_deg"], r3["mw_per_deg"], " ".join(f"{v:g}" for v in lay)])
        line.append(f"{r3['mw_per_deg']:.2f}" + ("" if r3["range_deg"] > 22.7 else f"({r3['range_deg']:.0f}°)"))
    print(f"  cavity {D:3d} µm, slot {w:2d} µm: " + "  ".join(f"{p:>6g} Pa {v}" for p, v in zip(pressures, line)))
with open(RESULTS / "gas_loss.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["cavity_um", "slot_um", "pressure_Pa", "range_deg", "mw_per_deg", "heaters_um"])
    wr.writerows(rows)

# ---------------------------------------------------------------- heavy gases at 1 atm
print("\nA heavy gas sealed in at 1 atm instead of vacuum (conductivity at 300 K, standard values):")
for name, k in (("air", 0.026), ("argon", 0.0177), ("krypton", 0.0094), ("xenon", 0.0055)):
    h = a_h * k / 15e-6
    G_s = a_s * k * T_MEM / 10e-6
    lay, r3 = best(profiles(h, G_s), 3)
    print(f"  {name:8} k = {k:.4f} W/m/K, 15 µm cavity, 10 µm slot: {r3['range_deg']:.1f}° at {r3['mw_per_deg']:.3f} mW/deg")
