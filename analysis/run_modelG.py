"""Model G (R1: heaters out to the cavity wall, R2: slot through the membrane,
R3: the slot with its losses removed: vacuum, a 4 µm slot, deeper cavities).

R3 files are optional: each one found in the data folder is checked and steered.

Part 1 checks the COMSOL exports against what they must satisfy (energy
balance, agreement with kernel_undercut.csv) and against what the 1D membrane
model predicted (expect agreement to about 10%, since the membrane model fitted
today's kernel to 7-9%).

Part 2 redoes the steering with the new kernels. The slot design is one-sided
by construction: heat enters at the right-hand end and can only leave through
the left wall, so it steers toward +θ only. A one-sided range R still serves
±R/2 once the static trim builds in a fixed tilt of R/2, so every case here is
ranked by its range toward +θ (Strehl ≥ 0.8), then by power per degree, over
an exhaustive grid of heater positions. Today's kernel goes through the same
search, so the comparison is like for like.

Writes results/modelG_checks.csv and results/modelG_steering.csv.
Usage: python run_modelG.py [data_dir]   (default: the Project folder)
"""
import csv
import itertools
import sys
from pathlib import Path

import numpy as np

from opa.comsol_io import read_comsol
from opa.constants import DATA_DIR, PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.placement import is_feasible
from opa.steering import steer_one_direction
from opa.thermal import SweptKernel, UndercutKernel, emitter_positions

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)
data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else DATA_DIR

xe = emitter_positions()
Y_CORE = 0.11            # µm, where the kernels were exported
HEATER_LENGTH_M = 1e-3   # for absolute temperatures only; phase does not depend on it
DPHI_10DEG = 2 * np.pi * PITCH_UM * np.sin(np.radians(10)) / WAVELENGTH_UM  # ramp step for 10°

checks = []  # rows: [run, quantity, value, expected, verdict]


def check(run, quantity, value, expected, ok):
    verdict = "ok" if ok else "LOOK"
    checks.append([run, quantity, f"{value:.4g}", expected, verdict])
    print(f"  [{verdict:4}] {quantity:<52} {value:9.4g}   (expected {expected})")


def energy_check(run, filename, n_expected):
    """Heat leaving through the substrate bottom must equal the 1 W/m put in."""
    try:
        _, table = read_comsol(data_dir / filename)
    except FileNotFoundError:
        print(f"  {filename} not found")
        return
    flux = np.abs(table[:, -1])
    check(run, f"energy balance: runs in {filename}", len(flux), f"{n_expected}", len(flux) == n_expected)
    check(run, "energy balance: worst |flux - 1 W/m|", np.max(np.abs(flux - 1)), "< 0.01", np.max(np.abs(flux - 1)) < 0.01)


def load(filename, name, mirror, n_columns):
    try:
        return SweptKernel(filename, name, mirror, expected_columns=n_columns, data_dir=data_dir)
    except (FileNotFoundError, ValueError) as err:
        print(f"  cannot use {filename}: {err}")
        return None


def at(kernel, x_h, x):
    """Kernel value at position x for a heater at grid position x_h (K per W/m)."""
    return float(np.interp(x, kernel.x, kernel.column(x_h)))


# ---------------------------------------------------------------- Part 1: checks
today = UndercutKernel(data_dir)

print("R1: heaters out to the wall, slot closed (kernel_undercut_far.csv)")
energy_check("R1", "checks_energy_far.csv", 41)
far = load("kernel_undercut_far.csv", "far (no slot)", mirror=True, n_columns=41)
if far is not None:
    # Columns 0-32 repeat today's kernel; only the slot rectangle (filled with oxide) differs.
    worst = max(np.max(np.abs(far.column(j) - today.k[:, j])) / today.k[:, j].max() for j in range(33))
    check("R1", "x_h = 0-32 vs kernel_undercut.csv, worst diff / peak", worst, "< 0.005", worst < 0.005)
    identical = all(np.array_equal(far.column(j), today.k[:, j]) for j in range(33))
    check("R1", "bit-identical to kernel_undercut.csv (stale export?)", float(identical), "0", not identical)
    check("R1", "centre core, heater at 0", at(far, 0, 0), "2.850 ± 0.5%", abs(at(far, 0, 0) / 2.850 - 1) < 0.005)
    for x_h, expect in [(34, 0.43), (36, 0.29), (38, 0.14)]:
        v = at(far, x_h, 0)
        check("R1", f"centre core, heater at {x_h}", v, f"≈ {expect} (membrane model)", abs(v / expect - 1) < 0.25)
    centre = [at(far, x_h, 0) for x_h in range(32, 41)]
    check("R1", "centre core falls as the heater nears the wall", float(np.all(np.diff(centre) < 0)), "1", np.all(np.diff(centre) < 0))

print("\nR2: slot open at 36-38 µm (kernel_slot.csv)")
energy_check("R2", "checks_energy_slot.csv", 67)
slot = load("kernel_slot.csv", "slot", mirror=False, n_columns=67)
if slot is not None:
    v = at(slot, 0, 0)
    check("R2", "centre core, heater at 0 (2.850 means slot not open)", v, "≈ 4.4", abs(v / 4.4 - 1) < 0.25)
    for x, expect in [(-32, 0.76), (0, 4.05), (32, 8.4)]:
        v = at(slot, 34, x)
        check("R2", f"heater at 34: core at {x:+d} µm", v, f"≈ {expect}", abs(v / expect - 1) < 0.25)
    # The ramp: fit a straight line to the heater-34 profile over the cores and
    # compare the misfit with the rise across the array.
    t = np.interp(xe, slot.x, slot.column(34))
    line = np.polyval(np.polyfit(xe, t, 1), xe)
    rms_frac = np.sqrt(np.mean((t - line) ** 2)) / (line[-1] - line[0])
    check("R2", "heater at 34: rms off a straight line / rise", rms_frac, "≈ 0.022", rms_frac < 0.05)
    stub = at(slot, 34, 39) / at(slot, 34, 32)
    check("R2", "heater at 34: stub (x = 39) / last core (x = 32)", stub, "≪ 1", stub < 0.2)
    # Mirror symmetry must be broken: compare K(x; -20) with K(-x; +20).
    asym = np.max(np.abs(slot.column(-20) - slot.column(20)[::-1])) / slot.column(20).max()
    check("R2", "asymmetry, heater at ±20 (0 means slot not open)", asym, "> 0.05", asym > 0.05)

# R3: copies of the Slot study with the losses removed one at a time. The
# centre-core values come from the membrane model with a leaky slot, scaled by
# R2's COMSOL/model ratio (0.85) for the air cases.
R3_CASES = [  # label, kernel file, energy file, centre core for a heater at 34 µm (K per W/m)
    ("R3a vacuum", "kernel_vacuum.csv", "checks_energy_vacuum.csv", 7.0),
    ("R3b slot 4 µm", "kernel_4um.csv", "checks_energy_4um.csv", 2.8),
    ("R3c slot 4 µm, D 50", "kernel_r3c_slot4_D50.csv", "checks_energy_r3c.csv", 3.7),
    ("R3d slot 4 µm, D 100", "kernel_r3d_slot4_D100.csv", "checks_energy_r3d.csv", 4.0),
]
# With no losses, all of the heater's power q crosses the array to the left wall,
# so the temperature slope is exactly q/(k·t). The ramp step is then
# Δφ = 814·P·d/(k·t) whatever N, the heater length or its position:
# P = k·t·Δφ/(814·d), which is about 0.50 mW/deg.
MEMBRANE_KT = 1.38 * 4.22e-6  # W/K: k_ox × (BOX + core layer + cladding)
U_PER_RAD = WAVELENGTH_UM / (2 * np.pi * PITCH_UM)  # sin θ per rad of ramp step
CONDUCTION_LIMIT = MEMBRANE_KT / (PHASE_COEFF * PITCH_UM * 1e-6) / U_PER_RAD * np.pi / 180 * 1e3  # mW/deg

r3 = {}
for label, kernel_file, energy_file, centre in R3_CASES:
    if not (data_dir / kernel_file).exists():
        continue
    print(f"\n{label} ({kernel_file})")
    kernel = load(kernel_file, label, mirror=False, n_columns=None)
    if kernel is None:
        continue
    n_runs = len(kernel.x_h)
    check(label, "kernel columns (34 for a 2 µm step, 67 for 1 µm)", n_runs, "34 or 67", n_runs in (34, 67))
    energy_check(label, energy_file, n_runs)
    v = at(kernel, 34, 0)
    check(label, "centre core, heater at 34", v, f"≈ {centre}", abs(v / centre - 1) < 0.25)
    t = np.interp(xe, kernel.x, kernel.column(34))
    line = np.polyval(np.polyfit(xe, t, 1), xe)
    rms_frac = np.sqrt(np.mean((t - line) ** 2)) / (line[-1] - line[0])
    check(label, "heater at 34: rms off a straight line / rise", rms_frac,
          "< 0.01" if "vacuum" in label else "< R2's 0.026", rms_frac < (0.01 if "vacuum" in label else 0.026))
    if "vacuum" in label:
        r = steer_one_direction(PHASE_COEFF * t[:, None], xe, +1)
        check(label, "one heater at 34: mW/deg against k·t/(814·d)", r["mw_per_deg"] / CONDUCTION_LIMIT,
              f"≈ 1 (limit {CONDUCTION_LIMIT:.3f} mW/deg)", abs(r["mw_per_deg"] / CONDUCTION_LIMIT - 1) < 0.1)
        check(label, "one heater at 34: range (deg)", r["range_deg"], "22.8 (full)", r["range_deg"] > 22)
    r3[label] = kernel

print("\nR2: hottest point in the heater (heater_Tmax_slot.csv)")
try:
    names, tmax = read_comsol(data_dir / "heater_Tmax_slot.csv")
    x_h_cols = [i for i, n in enumerate(names) if n.strip().startswith("x_h")]
    if not x_h_cols:  # e.g. an X, Y, T profile from Export -> Data
        raise ValueError(f"no x_h column (columns: {names}); expected the Surface Maximum table, one row per x_h")
    col = x_h_cols[0]
    heater_tmax = dict(zip(np.round(tmax[:, col], 6), tmax[:, -1]))
    check("R2", "hottest point in the heater at 34 (K per W/m)", heater_tmax[34.0], "≈ 9-10", 7 < heater_tmax[34.0] < 13)
except (FileNotFoundError, ValueError, KeyError) as err:
    heater_tmax = None
    print(f"  cannot use heater_Tmax_slot.csv: {err}")

with open(RESULTS / "modelG_checks.csv", "w", newline="", encoding="utf-8") as f:
    csv.writer(f).writerows([["run", "quantity", "value", "expected", "verdict"], *checks])


# ---------------------------------------------------------------- Part 2: steering
def best_one_way(kernel, n_heaters, lo, hi, step):
    """Exhaustive search over heater positions on a grid, ranked by +θ range, then mW/deg."""
    grid = np.arange(lo, hi + step / 2, step)
    columns = {xh: PHASE_COEFF * kernel(xe, xh) for xh in grid}  # phase per W, per heater position
    best_key, best = None, None
    for combo in itertools.combinations(grid, n_heaters):
        if not is_feasible(np.array(combo), (lo, hi)):
            continue
        r = steer_one_direction(np.column_stack([columns[xh] for xh in combo]), xe, +1)
        key = (round(r["range_deg"], 2), -r["mw_per_deg"])
        if best_key is None or key > best_key:
            best_key, best = key, (np.array(combo), r)
    return best


def summary(kernel, x_heat, r):
    """Power at 10° (if reachable) and the hottest core then, for 1 mm heaters."""
    if r["range_deg"] < 10:
        return np.nan, np.nan
    powers = DPHI_10DEG * r["powers_per_rad"]  # W
    rise = sum(p * kernel(xe, xh) for p, xh in zip(powers, x_heat)) / HEATER_LENGTH_M
    return powers.sum() * 1e3, rise.max()


rows = []


def report(case, kernel, x_heat, r):
    p10, t10 = summary(kernel, x_heat, r)
    where = ", ".join(f"{x:g}" for x in x_heat)
    print(f"  {case:<26} {len(x_heat):2d}  {where:<18} {r['range_deg']:5.1f}°  {r['mw_per_deg']:6.2f}  "
          f"{p10:7.1f}  {t10:7.1f}")
    rows.append([case, len(x_heat), where, r["range_deg"], r["mw_per_deg"], p10, t10])


header = (f"  {'case':<26} {'K':>2}  {'heaters (µm)':<18} {'range':>6}  {'mW/deg':>6}  "
          f"{'P(10°)':>7}  {'ΔT(10°)':>7}")
print("\nSteering toward +θ, Strehl ≥ 0.8. P(10°) in mW; ΔT(10°) = hottest core in K for 1 mm heaters.")
print(header)

# Reference: today's cavity, heaters within the array.
for n in (1, 2, 3):
    x, r = best_one_way(today, n, -32, 32, 1.0 if n < 3 else 2.0)
    report("today (kernel_undercut)", today, x, r)

if far is not None:
    for x_h in range(32, 41, 2):
        report("R1, single heater", far, [x_h], steer_one_direction(PHASE_COEFF * far(xe, x_h)[:, None], xe, +1))
    for n in (2, 3):
        x, r = best_one_way(far, n, -40, 40, 1.0 if n < 3 else 2.0)
        report("R1, best layout", far, x, r)

if slot is not None:
    report("R2, single heater", slot, [34], steer_one_direction(PHASE_COEFF * slot(xe, 34)[:, None], xe, +1))
    # The prediction: 34 µm plus one interior heater (membrane model said 5 µm).
    report("R2, heaters 34 + 5", slot, [34, 5],
           steer_one_direction(PHASE_COEFF * np.column_stack([slot(xe, 34), slot(xe, 5)]), xe, +1))
    for n in (2, 3):
        x, r = best_one_way(slot, n, -32, 34, 1.0 if n < 3 else 2.0)
        report("R2, best layout", slot, x, r)
    if heater_tmax is not None:
        r = steer_one_direction(PHASE_COEFF * slot(xe, 34)[:, None], xe, +1)
        p_per_deg = r["mw_per_deg"] * 1e-3  # W per degree
        print(f"\n  Heater at 34 µm self-heating: {heater_tmax[34.0]:.2f} K per W/m, so "
              f"{heater_tmax[34.0] * p_per_deg / HEATER_LENGTH_M:.1f} K per degree of steering for a 1 mm heater.")

for label, kernel in r3.items():
    report(f"{label}, one heater", kernel, [34], steer_one_direction(PHASE_COEFF * kernel(xe, 34)[:, None], xe, +1))
    for n in (2, 3):
        x, r = best_one_way(kernel, n, -32, 34, 2.0)  # the R3 sweeps step x_h by 2 µm
        report(f"{label}, best", kernel, x, r)

with open(RESULTS / "modelG_steering.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["case", "K", "heaters_um", "range_plus_deg", "mw_per_deg", "P_10deg_mW", "core_dT_10deg_K_1mm"])
    writer.writerows(rows)
