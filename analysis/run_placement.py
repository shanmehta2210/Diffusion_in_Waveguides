"""Analysis 2: best heater positions for K = 1 to 8, in both stacks.

Ranking: steering span at Strehl ≥ 0.8 first, power per degree second (see
opa/placement.py). Heaters stay within the array, ±32 µm: the undercut kernel
was only computed for heaters there. For K = 1 and 2 an exhaustive 1 µm grid
checks the optimiser. Writes results/placement.csv.
"""
import csv
import sys
import time
from pathlib import Path

import numpy as np

from opa.placement import evaluate_layout, grid_search, optimise_layout, rank_key
from opa.thermal import StandardKernel, UndercutKernel, emitter_positions

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)

xe = emitter_positions()
bounds = (xe[0], xe[-1])
rows = []
for kernel in [StandardKernel(), UndercutKernel()]:
    print(f"\n{kernel.name} stack")
    print(f"{'K':>3} {'+range':>7} {'-range':>7} {'span':>6} {'mW/deg':>7} {'even span':>9}  positions (µm)")
    for n_heaters in range(1, 9):
        t0 = time.time()
        x_best, best = optimise_layout(kernel, xe, n_heaters, bounds, n_random_starts=30)
        even = evaluate_layout(kernel, xe, np.array([xe[-1]]) if n_heaters == 1 else np.linspace(*bounds, n_heaters))
        print(f"{n_heaters:3d} {best['plus']['range_deg']:6.2f}° {best['minus']['range_deg']:6.2f}° "
              f"{best['span_deg']:5.1f}° {best['mw_per_deg']:7.2f} {even['span_deg']:8.1f}°  "
              f"{np.array2string(x_best, precision=2, separator=', ')}  ({time.time() - t0:.0f} s)")
        rows.append([kernel.name, n_heaters, best["plus"]["range_deg"], best["minus"]["range_deg"],
                     best["span_deg"], best["mw_per_deg"], even["span_deg"], " ".join(f"{x:g}" for x in x_best)])
        if n_heaters <= 3:
            # 1 µm grid for K = 1, 2; 2 µm for K = 3 to keep it to seconds.
            x_grid, grid = grid_search(kernel, xe, n_heaters, bounds, step=1.0 if n_heaters <= 2 else 2.0)
            # Pass if the optimiser spans as far (to 0.05°) and, at equal span, needs at most 1% more
            # power. Mirror-image layouts tie, so the grid can legitimately pick the other one.
            as_good = best["span_deg"] >= grid["span_deg"] - 0.05 and (
                best["span_deg"] > grid["span_deg"] + 0.05 or best["mw_per_deg"] <= 1.01 * grid["mw_per_deg"])
            print(f"    grid check: span {grid['span_deg']:.1f}°, {grid['mw_per_deg']:.2f} mW/deg at {x_grid} -> "
                  f"{'optimiser at least as good' if as_good else 'GRID BETTER'}")

with open(RESULTS / "placement.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["stack", "K", "range_plus_deg", "range_minus_deg", "span_deg", "mw_per_deg",
                     "even_span_deg", "heater_positions_um"])
    writer.writerows(rows)
