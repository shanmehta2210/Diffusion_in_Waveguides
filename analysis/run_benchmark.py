"""Analysis 3: thermal steering against the conventional design.

Designs compared at each target angle:
  conventional        one isolated heater per emitter, wrapped phases, P_π = 19.7 mW
  conventional+xtalk  the same layout with the real standard-stack crosstalk:
                      cheapest of wrapped (each reset placement) and unwrapped
                      that reaches Strehl ≥ 0.8
  thermal designs     the optimised layouts from run_placement.py (results/placement.csv)

For each: total power, Strehl toward the target, and the largest phase on any
emitter. That phase sets the hottest core: ΔT = φ / (PHASE_COEFF · L) for a
heater of length L, about 1.23 K per rad for 1 mm. Writes results/benchmark.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np

from opa.constants import PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.conventional import candidate_offsets, ideal_power, with_crosstalk, wrapped_ramp
from opa.steering import STREHL_MIN, steer_one_direction, strehl_vs_step
from opa.thermal import StandardKernel, UndercutKernel, emitter_positions, phase_matrix

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"

xe = emitter_positions()
std, uc = StandardKernel(), UndercutKernel()
kernels = {"standard": std, "undercut": uc}
p_pi = np.pi / (PHASE_COEFF * std(0, 0))
angles = [1, 2, 5, 8, 10, 12, 15, 16, 20]

with open(RESULTS / "placement.csv") as f:
    layouts = {(r["stack"], int(r["K"])): np.array(r["heater_positions_um"].split(), dtype=float)
               for r in csv.DictReader(f)}
thermal = [("undercut", 3), ("undercut", 4), ("undercut", 8), ("standard", 8)]
g_conv = phase_matrix(std, xe, xe)


def ideal_max_phase(dphi):
    """Largest wrapped phase in the least-power isolated solution."""
    best = min(candidate_offsets(xe, dphi), key=lambda c: wrapped_ramp(xe, dphi, c).sum())
    return wrapped_ramp(xe, dphi, best).max()


def thermal_at(G, fit, dphi):
    """Power, Strehl and largest phase of a thermal design at ramp step dphi (powers scale with dphi)."""
    powers = dphi * fit["powers_per_rad"]
    return powers.sum(), strehl_vs_step(fit["residual_per_rad"], [dphi])[0], (G @ powers).max()


rows = []
names = ["conventional", "conv.+xtalk"] + [f"{s} K={k}" for s, k in thermal]
print("Power in mW (Strehl in brackets; * = below 0.8). Steering toward +θ.\n")
print(f"{'angle':>5} " + " ".join(f"{n:>17}" for n in names))
for deg in angles:
    dphi = 2 * np.pi * PITCH_UM * np.sin(np.radians(deg)) / WAVELENGTH_UM
    results = [(ideal_power(xe, dphi, p_pi), 1.0, ideal_max_phase(dphi))]
    power, s, _ = with_crosstalk(g_conv, xe, dphi, np.sin(np.radians(deg)))
    results.append((power, s, np.nan))
    for stack, k in thermal:
        G = phase_matrix(kernels[stack], xe, layouts[(stack, k)])
        results.append(thermal_at(G, steer_one_direction(G, xe, +1), dphi))
    cells = [f"{p * 1e3:8.0f} ({s:.2f}){'*' if s < STREHL_MIN else ' '}" for p, s, _ in results]
    print(f"{deg:4d}° " + " ".join(f"{c:>17}" for c in cells))
    rows += [[name, deg, p * 1e3, s, phi] for name, (p, s, phi) in zip(names, results)]

print("\nHottest core rise for 1 mm heaters (K), from the largest phase:")
for name in ["conventional"] + names[2:4]:
    rises = [r[4] / (PHASE_COEFF * 1e-3) for r in rows if r[0] == name]
    print(f"  {name:>15}: " + "  ".join(f"{deg}°: {t:4.0f}" for deg, t in zip(angles, rises)))

with open(RESULTS / "benchmark.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["design", "angle_deg", "power_mw", "strehl", "max_phase_rad"])
    writer.writerows(rows)
