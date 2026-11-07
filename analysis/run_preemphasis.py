"""Pre-emphasis: how much faster does the beam settle if the heater drives are shaped in time?

Builds on the R4 step responses (opa/transient.py) and the optimiser in
opa/preemphasis.py. Three ways to drive the same heaters, from the same start
to the same target angle:
  - plain step: every heater jumps to its final power (R4);
  - shared waveform: one overdrive waveform scales all heaters together, so the
    profile keeps the steady shape (classic single-channel pre-emphasis);
  - independent waveforms: each heater gets its own waveform (multi-heater
    pre-emphasis). For the one-heater vacuum design the last two coincide.
"Settled" means the Strehl toward the target stays ≥ 0.8 from then on. Every
shaped drive ends on the plain step's hold drive, so all three end in the same
steady state (same final Strehl and pointing) and only the way there differs.

The power limit is on the whole chip: at every moment Σ_k P_k(t) ≤ m·P10,
where P10 is the hold power at 10° for that design. A chip built for m = 2 can
draw twice its 10° hold power for a moment.

Cases: cold start to 10° against m; the same with no core ever more than 5%
above its final temperature; slews from a held angle (where heaters can also be turned
down); more heaters (slot); sensitivity to a wrong thermal time scale; and
robust drives that tolerate it.

Writes results/preemphasis.csv (A-D), preemphasis_robust.csv (E-F) and
preemphasis_traces.csv (A). Sections can be run alone: python run_preemphasis.py E F.
"""
import csv
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from opa.constants import PITCH_UM, WAVELENGTH_UM
from opa.preemphasis import (PulseModel, SIGMA_STREHL, fastest_drive, fine_times, rms_error,
                             strehl, verify)
from opa.steering import STREHL_MIN, steer_one_direction
from opa.thermal import emitter_positions
from opa.transient import StepResponses, settle_time

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
xe = emitter_positions()


def dphi(theta_deg):
    """Ramp step per emitter (rad) for steering angle θ."""
    return 2 * np.pi * PITCH_UM * np.sin(np.radians(theta_deg)) / WAVELENGTH_UM


def ramp(theta_deg):
    return dphi(theta_deg) * xe / PITCH_UM


def pointing_deg(phase):
    """Beam angle from the least-squares slope of the phase across the array, per row."""
    slope = np.polyfit(xe, np.atleast_2d(phase).T, 1)[0]            # rad per µm
    return np.degrees(np.arcsin(slope * WAVELENGTH_UM / (2 * np.pi)))


def shared_waveform_model(sr, heaters, p_unit):
    """One virtual heater = all heaters in the fixed steady ratio; its power is the chip's total."""
    cols = [int(np.flatnonzero(np.isclose(sr.x_h, xh))[0]) for xh in heaters]
    S = np.einsum("tkn,k->tn", sr.S[:, cols], p_unit / p_unit.sum())[:, None, :]
    return PulseModel(SimpleNamespace(times=sr.times, x_h=np.array([0.0]), S=S), [0.0])


def time_stretched(sr, factor):
    """The same step responses on a time axis stretched by factor (>1: slower heat flow)."""
    return SimpleNamespace(times=sr.times * factor, x_h=sr.x_h, S=sr.S)


class Design:
    def __init__(self, name, sr, heaters):
        self.name, self.sr, self.heaters = name, sr, list(heaters)
        self.M = PulseModel(sr, heaters)
        self.p_unit = steer_one_direction(self.M.G, xe, +1)["powers_per_rad"]   # W per rad of Δφ
        self.P10 = dphi(10) * self.p_unit.sum()
        self.shared = shared_waveform_model(sr, heaters, self.p_unit) if len(heaters) > 1 else self.M

    def drive(self, theta):
        """Steady least-squares drive (W per heater) for angle θ."""
        return dphi(theta) * self.p_unit


def step_time(M, P_start, P_end, target):
    t = fine_times(0.03, 6000)
    return settle_time(t, strehl(M.step_drive_phase(t, P_end, P_start) - target), STREHL_MIN)


def optimise(M, P_start, target, cap, T_hi, **kw):
    """Fastest drive, checked on a fine grid. Returns (drive, verification, extra energy over the hold in µJ)."""
    d = fastest_drive(M, P_start, target, cap, T_hi=T_hi, **kw)
    v = verify(M, d, P_start, target)
    dt = np.diff(d["t_switch"])
    P_tot = d["P"].sum(axis=0)
    extra = float(np.sum(dt * (P_tot[:-1] - P_tot[-1]))) * 1e6
    return d, v, extra


rows = []
section = ""   # set by each section below; tags its rows


def record(design, mode, start, end, m, t_step, v=None, extra=0.0, note="", start_peak=0.0):
    """One result row. start_peak: hottest core phase of the starting state (overshoot reference)."""
    peak_ref = max(v["final_peak_phase"], start_peak) if v else np.nan
    t_pre = v["t_strehl"] if v else t_step
    hold = dphi(end) * design.p_unit.sum()
    row = {
        "section": section, "design": design.name,
        "heaters_um": " ".join(f"{h:g}" for h in design.heaters), "mode": mode,
        "start_deg": start, "end_deg": end, "cap_x_P10": m, "P10_mW": design.P10 * 1e3,
        "step_ms": t_step * 1e3, "settle_ms": t_pre * 1e3, "speedup": t_step / t_pre,
        "peak_over_final": v["peak_phase"] / peak_ref if v else 1.0,
        "extra_energy_uJ": extra, "hold_mW": hold * 1e3, "hold_x_settle_mWms": hold * 1e3 * t_pre * 1e3,
        "note": note,
    }
    rows.append(row)
    print(f"  {design.name:8} {mode:12} {start:4.1f}->{end:4.1f}°  m={m:<4g} step {t_step*1e3:6.3f} ms  "
          f"-> {t_pre*1e3:6.3f} ms  (×{row['speedup']:4.1f})  peak/final {row['peak_over_final']:.2f}  "
          f"extra {extra:6.1f} µJ  {note}")
    return row


def run_on(design, drv, factor, P_start, target):
    """Strehl settling time when the drive meets a model whose time scale is factor × nominal."""
    Mf = PulseModel(time_stretched(design.sr, factor), design.heaters)
    t = fine_times(drv["t_switch"][-1])
    phi = Mf.phase(t, drv["t_switch"], drv["P"], P_start)
    return settle_time(t, strehl(phi - target), STREHL_MIN)


drives = {}   # (design, start°, end°, m) -> fastest drive with independent waveforms


def nominal_drive(d, a, b, m):
    key = (d.name, a, b, m)
    if key not in drives:
        t_step = step_time(d.M, d.drive(a), d.drive(b), ramp(b))
        drives[key] = fastest_drive(d.M, d.drive(a), ramp(b), m * d.P10, T_hi=1.1 * t_step, P_final=d.drive(b))
    return drives[key]


def save(path, new_rows, sections):
    """Write rows, keeping any rows already in the file from sections not run this time."""
    old = []
    if path.exists():
        with open(path, encoding="utf-8") as f:
            old = [r for r in csv.DictReader(f) if r.get("section", "") not in sections]
    fields = list(new_rows[0]) if new_rows else list(old[0])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(old + new_rows)


# Sections to run, e.g. "python run_preemphasis.py E F"; all by default. A full run takes over an hour.
SECTIONS = set(sys.argv[1:]) or set("ABCDEF")
t_start = time.time()
designs = [
    Design("vacuum", StepResponses("transient_R4vac.csv", default_x_h=34.0), [34]),
    Design("slot", StepResponses("transient_R4b.csv"), [0, 12, 32]),
    Design("undercut", StepResponses("transient_R4a.csv", mirror=True), [-12, 8, 32]),
]
by_name = {d.name: d for d in designs}
for d in designs:
    print(f"{d.name}: heaters {d.heaters} µm, 10° hold drive "
          f"{', '.join(f'{p*1e3:.2f}' for p in d.drive(10))} mW (total {d.P10*1e3:.2f} mW)")

if "A" in SECTIONS:
    print("\nA. Cold start 0 -> 10°, chip power ≤ m × P10")
    section = "A"
    for d in designs:
        P0, r = np.zeros(d.M.K), ramp(10)
        t_step = step_time(d.M, P0, d.drive(10), r)
        record(d, "plain step", 0, 10, 1, t_step)
        for m in (1, 1.5, 2, 3, 5, 10):
            modes = [("shared", d.shared)] + ([("independent", d.M)] if d.M.K > 1 else [])
            for mode, M in modes:
                if mode == "shared" and m == 1:
                    continue                 # a shared waveform capped at the hold power is the plain step
                Ps = np.zeros(M.K)
                Pf = d.drive(10) if M is d.M else [d.drive(10).sum()]   # end on the step's hold drive
                drv, v, extra = optimise(M, Ps, r, m * d.P10, 1.1 * t_step, P_final=Pf)
                record(d, mode if d.M.K > 1 else "overdrive", 0, 10, m, t_step, v, extra)
                if M is d.M:
                    drives[(d.name, 0, 10, m)] = drv

if "B" in SECTIONS:
    print("\nB. Cold start 0 -> 10°, m = 3, and no core ever more than 5% above the final peak")
    section = "B"
    for d in designs:
        P0, r = np.zeros(d.M.K), ramp(10)
        t_step = step_time(d.M, P0, d.drive(10), r)
        cap_phase = 1.05 * float((d.M.G @ d.drive(10)).max())   # exactly 1.0 makes the LP degenerate
        drv, v, extra = optimise(d.M, P0, r, 3 * d.P10, 1.1 * t_step, P_final=d.drive(10),
                                 phase_cap=cap_phase)
        record(d, "independent" if d.M.K > 1 else "overdrive", 0, 10, 3, t_step, v, extra,
               note="overshoot ≤ 5%")

if "C" in SECTIONS:
    print("\nC. Slews from a held angle (m = 2): heaters can now also be turned down")
    section = "C"
    for d in designs:
        for a, b in ((0, 1), (5, 6), (5, 4), (5, 10), (10, 5), (10, 0)):
            Pa, r = d.drive(a), ramp(b)
            t_step = step_time(d.M, Pa, d.drive(b), r)
            drv, v, extra = optimise(d.M, Pa, r, 2 * d.P10, 1.1 * t_step, P_final=d.drive(b))
            record(d, "independent" if d.M.K > 1 else "overdrive", a, b, 2, t_step, v, extra,
                   start_peak=float((d.M.G @ Pa).max()))
            drives[(d.name, a, b, 2)] = drv

if "D" in SECTIONS:
    print("\nD. Slot: what more heaters buy (cold start 0 -> 10°)")
    section = "D"
    slot3 = by_name["slot"]
    for heaters, caps in (([4, 32], (2, 3)), (list(slot3.sr.x_h), (3,))):   # all 17: slow, so m = 3 only
        d = Design("slot", slot3.sr, heaters)
        for m in caps:
            P0 = np.zeros(d.M.K)
            t_step = step_time(d.M, P0, d.drive(10), ramp(10))
            drv, v, extra = optimise(d.M, P0, ramp(10), m * d.P10, 1.1 * t_step, P_final=d.drive(10))
            record(d, "independent", 0, 10, m, t_step, v, extra, note=f"K = {len(heaters)}")

robust_rows = []
CASES = [("slot", 0, 10, 3), ("vacuum", 0, 10, 3), ("undercut", 0, 10, 3), ("slot", 5, 6, 2), ("vacuum", 5, 6, 2)]

if "E" in SECTIONS:
    print("\nE. Sensitivity: the nominal fastest drive, run on a model whose thermal time scale is off")
    for name, a, b, m in CASES:
        d = by_name[name]
        drv = nominal_drive(d, a, b, m)
        t_step = step_time(d.M, d.drive(a), d.drive(b), ramp(b))
        for factor in (0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.1):
            ts = run_on(d, drv, factor, d.drive(a), ramp(b))
            print(f"  {name:8} {a:g}->{b:g}° m={m}  time scale ×{factor:.2f}: {ts*1e3:6.3f} ms  "
                  f"(step {factor*t_step*1e3:6.3f} ms)")   # a plain step just rescales in time
            robust_rows.append({"section": "E", "design": name, "start_deg": a, "end_deg": b, "cap_x_P10": m,
                                "drive": "nominal", "band": 0, "time_scale": factor,
                                "settle_ms": ts * 1e3, "step_ms": factor * t_step * 1e3})

if "F" in SECTIONS:
    print("\nF. Robust design: one drive that must settle on models with time scale 1-b, 1 and 1+b")
    for name, a, b, m in [c for c in CASES if c[0] != "undercut"]:
        d = by_name[name]
        t_step = step_time(d.M, d.drive(a), d.drive(b), ramp(b))
        nominal = nominal_drive(d, a, b, m)
        for band in (0.02, 0.05, 0.10):
            scen = [PulseModel(time_stretched(d.sr, f), d.heaters) for f in (1 - band, 1.0, 1 + band)]
            drv = fastest_drive(scen, d.drive(a), ramp(b), m * d.P10, T_hi=1.1 * (1 + band) * t_step,
                                P_final=d.drive(b))
            factors = np.linspace(1 - band, 1 + band, 9)      # checked between the design points too
            times = [run_on(d, drv, f, d.drive(a), ramp(b)) for f in factors]
            worst_nominal = max(run_on(d, nominal, f, d.drive(a), ramp(b)) for f in factors)
            print(f"  {name:8} {a:g}->{b:g}° m={m}  ±{band:.0%}: robust drive worst {max(times)*1e3:6.3f} ms, "
                  f"nominal drive worst {worst_nominal*1e3:6.3f} ms, step worst {(1+band)*t_step*1e3:6.3f} ms")
            for f, ts in zip(factors, times):
                robust_rows.append({"section": "F", "design": name, "start_deg": a, "end_deg": b, "cap_x_P10": m,
                                    "drive": f"robust ±{band:.0%}", "band": band, "time_scale": f,
                                    "settle_ms": ts * 1e3, "step_ms": f * t_step * 1e3})

# ---------------------------------------------------------------- outputs
if SECTIONS & set("ABCD"):
    save(RESULTS / "preemphasis.csv", rows, SECTIONS & set("ABCD"))
if SECTIONS & set("EF"):
    save(RESULTS / "preemphasis_robust.csv", robust_rows, SECTIONS & set("EF"))

if "A" in SECTIONS:   # beam path against time, step vs pre-emphasis (m = 3), for plotting
    with open(RESULTS / "preemphasis_traces.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["design", "drive", "t_ms", "total_power_mW", "angle_deg", "strehl", "rms_error_rad",
                    "peak_core_phase_rad"])
        for d in designs:
            drv = drives[(d.name, 0, 10, 3)]
            t = np.linspace(1e-6, 4 * drv["t_switch"][-1], 800)
            edges = np.append(drv["t_switch"], np.inf)
            for label in ("plain step", "pre-emphasis m=3"):
                if label == "plain step":
                    phi = d.M.step_drive_phase(t, d.drive(10), np.zeros(d.M.K))
                    P_tot = np.full_like(t, d.P10)
                else:
                    phi = d.M.phase(t, drv["t_switch"], drv["P"], np.zeros(d.M.K))
                    P_tot = drv["P"].sum(axis=0)[np.searchsorted(edges, t, side="right") - 1]
                for ti, pt, ang, s, e, pk in zip(t, P_tot, pointing_deg(phi), strehl(phi - ramp(10)),
                                                 rms_error(phi - ramp(10)), phi.max(axis=1)):
                    w.writerow([d.name, label, ti * 1e3, pt * 1e3, ang, s, e, pk])

print(f"\nσ* = {SIGMA_STREHL:.4f} rad. Done in {time.time() - t_start:.0f} s.")
