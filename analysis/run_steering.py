"""Analysis 1: steering by least squares, with evenly spaced heaters.

For each stack and each heater count K, heaters are spread evenly over the
array (K = 1: one heater at the right-hand end, which can only steer one way).
Prints the steering range each way at Strehl ≥ 0.8, the power per degree, and
full beam metrics at a few angles for two example layouts. Writes the table to
results/steering_even.csv.

Note on the range: with no heating at all, the Strehl toward an angle stays
above 0.8 for the first ~0.35° simply because the beam is that wide. Ranges at
that level mean "does not steer".
"""
import csv
import sys
from pathlib import Path

import numpy as np

from opa.constants import PITCH_UM, WAVELENGTH_UM
from opa.farfield import U_GRID, beam_metrics, element_pattern, far_field
from opa.steering import steer_both_directions
from opa.thermal import StandardKernel, UndercutKernel, emitter_positions, phase_matrix

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)

xe = emitter_positions()
kernels = [StandardKernel(), UndercutKernel()]


def even_layout(n_heaters):
    return np.array([xe[-1]]) if n_heaters == 1 else np.linspace(xe[0], xe[-1], n_heaters)


rows = []
print(f"{'stack':9} {'K':>3} {'+range':>7} {'-range':>7} {'span':>6}  {'mW/deg':>7}  {'rms resid':>9}")
for kernel in kernels:
    for n_heaters in [1, 2, 3, 4, 5, 6, 7, 8, len(xe)]:
        r = steer_both_directions(phase_matrix(kernel, xe, even_layout(n_heaters)), xe)
        rms = r["plus"]["rms_residual_per_rad"]
        print(f"{kernel.name:9} {n_heaters:3d} {r['plus']['range_deg']:6.2f}° {r['minus']['range_deg']:6.2f}° "
              f"{r['span_deg']:5.1f}°  {r['mw_per_deg']:7.2f}  {rms:9.3f}")
        rows.append([kernel.name, n_heaters, r["plus"]["range_deg"], r["minus"]["range_deg"],
                     r["span_deg"], r["mw_per_deg"], rms])
print("(rms resid: rms phase error per rad of ramp step, toward +θ; K = 33 is one heater per emitter)")

with open(RESULTS / "steering_even.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["stack", "K", "range_plus_deg", "range_minus_deg", "span_deg", "mw_per_deg", "rms_residual_per_rad"])
    writer.writerows(rows)

# Full far-field metrics at a few target angles. Powers scale with the ramp
# step, so the unit fit is reused: P(Δφ) = Δφ · P(1 rad).
ep = element_pattern()
for kernel, n_heaters in [(kernels[1], 4), (kernels[0], len(xe))]:
    G = phase_matrix(kernel, xe, even_layout(n_heaters))
    fit = steer_both_directions(G, xe)["plus"]
    print(f"\n{kernel.name}, K = {n_heaters}, evenly spaced:")
    print(f"  {'target':>7} {'beam':>7} {'Strehl':>7} {'SLS':>6} {'grating':>8} {'FWHM':>6} {'power':>9}")
    for target_deg in [2, 5, 10, 15, 20]:
        dphi = 2 * np.pi * PITCH_UM * np.sin(np.radians(target_deg)) / WAVELENGTH_UM
        m = beam_metrics(U_GRID, far_field(dphi * G @ fit["powers_per_rad"], xe, ep), ep)
        power_mw = dphi * fit["power_per_rad"] * 1e3
        print(f"  {target_deg:6.1f}° {m['angle_deg']:6.2f}° {m['strehl']:7.3f} {m['sls_db']:5.1f}dB "
              f"{m['grating_lobe_db']:6.1f}dB {m['fwhm_deg']:5.2f}° {power_mw:7.1f}mW")
