"""Analysis of the R4 transients: checks, rise times, and how fast a steered beam settles.

Cases: R4a undercut without slot (mirrored to negative heater positions),
R4b slot, R4c vacuum (one heater at 34 µm). Heaters step on at t = 0.

Beam settling uses linearity: for heater powers P_k switched on at t = 0, the
phase is φ_n(t) = 814·Σ_k P_k·S_k(x_n, t). The powers are the steady-state
least-squares fit for the target angle, so this is the plain step drive with no
pre-emphasis. Two settling measures:
  - "arrive": Strehl toward the target angle reaches 0.8 (a 0 -> 10° slew);
  - "slew of δ": by linearity a slew of size δ leaves the same error pattern
    scaled by δ/10°, so smaller corrections settle sooner. Reported as the time
    for the Strehl of the remaining error to reach 0.8.

Writes results/R4_settling.csv.
"""
import csv
import itertools
import sys
from pathlib import Path

import numpy as np

from opa.constants import PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.steering import STREHL_MIN, steer_one_direction
from opa.thermal import SweptKernel, UndercutKernel, emitter_positions
from opa.transient import StepResponses, crossing_time, settle_time

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
xe = emitter_positions()


def dphi(theta_deg):
    """Ramp step per emitter (rad) for steering angle θ."""
    return 2 * np.pi * PITCH_UM * np.sin(np.radians(theta_deg)) / WAVELENGTH_UM


def strehl(phase_error):
    """|mean exp(iε)|² over emitters, for each row of phase errors."""
    return np.abs(np.exp(1j * np.atleast_2d(phase_error)).mean(axis=1)) ** 2


cases = [
    ("R4a undercut", StepResponses("transient_R4a.csv", mirror=True), UndercutKernel()),
    ("R4b slot", StepResponses("transient_R4b.csv"), SweptKernel("kernel_slot.csv", "slot", False, 67)),
    ("R4c vacuum", StepResponses("transient_R4vac.csv", default_x_h=34.0),
     SweptKernel("kernel_vacuum.csv", "vacuum", False, 34)),
]

print("Checks: last output time (32 ms) against the steady kernels, overshoot, monotonic rise")
for name, sr, steady in cases:
    assert np.allclose(sr.x, xe), "transient points must be the 33 core centres"
    worst = max(np.max(np.abs(sr.final[:, k] - steady(xe, xh))) / steady(xe, xh).max()
                for k, xh in enumerate(sr.x_h))
    overshoot = np.max(sr.S - sr.S[-1]) / sr.S[-1].max()  # relative to the run's peak
    drops = np.sum(np.diff(sr.S, axis=0) < -1e-6 * sr.S[-1].max())
    print(f"  {name:13} heaters {', '.join(f'{v:g}' for v in sr.x_h)}")
    print(f"  {'':13} final vs steady: worst {worst:.2%} of peak;  overshoot {overshoot:+.2%};  "
          f"falling steps {drops}")
    if name == "R4c vacuum":  # this file has no x_h label: confirm it is the 34 µm run
        match = {xh: np.max(np.abs(sr.final[:, 0] - steady(xe, xh))) for xh in steady.x_h}
        print(f"  {'':13} best-matching steady column: x_h = {min(match, key=match.get):g} µm")

print("\nSingle-heater rise times (63% / 90% of final), heater at the right end")
for name, sr, _ in cases:
    xh = sr.x_h.max()
    s = sr.column(xh)
    parts = []
    for xc in (32, 0, -32):
        n = int(np.argmin(abs(xe - xc)))
        parts.append(f"core {xc:+3d}: {crossing_time(sr.times, s[:, n], 0.63) * 1e3:.2f} / "
                     f"{crossing_time(sr.times, s[:, n], 0.9) * 1e3:.2f} ms")
    print(f"  {name:13} heater {xh:g} µm   " + ";  ".join(parts))


def best_layout(sr, n_heaters):
    """Best one-way (+θ) layout on the run's heater grid, steady kernel = last output time."""
    G_all = PHASE_COEFF * sr.final
    best = None
    for combo in itertools.combinations(range(len(sr.x_h)), n_heaters):
        r = steer_one_direction(G_all[:, combo], xe, +1)
        key = (round(r["range_deg"], 2), -r["mw_per_deg"])
        if best is None or key > best[0]:
            best = (key, combo, r)
    return best[1], best[2]


rows = []
slews = (1, 2, 5, 10)
print("\nBeam settling for a step drive (no pre-emphasis), one-way layouts on each run's heater grid")
print(f"  {'case':13} {'K':>2}  {'heaters (µm)':<14} {'range':>6}  {'target':>6}  {'arrive':>8}  "
      + "  ".join(f"slew {d:>2}°" for d in slews))
for name, sr, _ in cases:
    for n_heaters in ((1,) if len(sr.x_h) == 1 else (1, 2, 3)):
        combo, fit = best_layout(sr, n_heaters)
        target = min(10.0, 0.9 * fit["range_deg"])
        P = dphi(target) * fit["powers_per_rad"]                      # W per heater
        phase = PHASE_COEFF * np.einsum("tkn,k->tn", sr.S[:, list(combo)], P)  # rad, (times, cores)
        ramp = dphi(target) * xe / PITCH_UM
        arrive = settle_time(sr.times, strehl(phase - ramp), STREHL_MIN)
        # Remaining error relative to the final state, scaled to slews of size δ.
        err = phase[-1] - phase
        slew_t = [settle_time(sr.times, strehl(err * dphi(d) / dphi(target)), STREHL_MIN) for d in slews]
        where = ", ".join(f"{sr.x_h[k]:g}" for k in combo)
        print(f"  {name:13} {n_heaters:2d}  {where:<14} {fit['range_deg']:5.1f}°  {target:5.1f}°  "
              f"{arrive * 1e3:6.2f}ms  " + "  ".join(f"{t * 1e3:6.2f}ms" for t in slew_t))
        rows.append([name, n_heaters, where, fit["range_deg"], fit["mw_per_deg"], target, arrive * 1e3,
                     *[t * 1e3 for t in slew_t]])

with open(RESULTS / "R4_settling.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["case", "K", "heaters_um", "range_deg", "mw_per_deg", "target_deg", "arrive_ms",
                *[f"slew_{d}deg_ms" for d in slews]])
    w.writerows(rows)
print("\n(slew δ: time until the error left by a δ-degree step has Strehl ≥ 0.8; "
      "slews larger than the range are not reachable and read as the 'arrive' case would)")
