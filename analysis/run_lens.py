"""The thermal lens: how strong a lens each layout makes, at what power, against conventional heaters.

Part 0  checks against analytic limits (membrane model in vacuum, isolated-heater sums)
Part 1  our membrane: best symmetric layouts of K heaters, a uniform wide heater,
        and the best single shaped heater (one wire, graded power density)
Part 2  conventional baselines: one isolated heater per emitter, unwrapped and wrapped
        (Fresnel), on the standard and undercut stacks; one heater per emitter
        on the shared stacks, crosstalk included
Part 3  steering and focusing at the same time
Part 4  how fast the lens settles (R4a step responses)

Lens strength a = edge-to-centre phase; focal length f = 2076 rad·µm / a
(see opa/lens.py). Kernels: today's undercut (heaters within ±32 µm) and R1's
far kernel (heaters out to ±40 µm). The slot kernels are one-sided and can't
make a symmetric dome, so they're left out.

Writes results/lens_layouts.csv and results/lens_comparison.csv.
"""
import csv
import itertools
import sys
from pathlib import Path

import numpy as np

from opa.constants import PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.lens import fit_lens, isolated_power, least_power_lens, lens_constant, unit_target
from opa.steering import STREHL_MIN, fit_powers
from opa.thermal import StandardKernel, SweptKernel, UndercutKernel, emitter_positions, phase_matrix
from opa.transient import StepResponses, settle_time

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
xe = emitter_positions()
AF = lens_constant(xe)            # a·f, rad·µm
A_REF = AF / 1000.0               # lens strength for f = 1 mm (≈ 2.08 rad)


def strehl(phase_error):
    return float(np.abs(np.exp(1j * phase_error).mean()) ** 2)


def ok(flag):
    return "ok  " if flag else "FAIL"


# ---------------------------------------------------------------- Part 0: checks
print("Part 0: checks against analytic limits")
# (a) Isolated heaters: P = (P_π/π)·a·Σ(1 - x²/X²), and Σ(1 - x²/X²) = 33 - 11968/1024 = 21.3125.
s = isolated_power(1.0, xe, np.pi, wrapped=False)
print(f"  [{ok(abs(s - 21.3125) < 1e-9)}] isolated, unwrapped: Σ(1 - x²/X²) = {s:.4f} (analytic 21.3125)")
# (b) Membrane model in vacuum (anchored at ±40 µm, no loss), one heater spread uniformly over
# the array: T'' = -q/(k·t), so the rise is an exact parabola and the lens is perfect, with
# P/a = 2·k·t·W/(814·X²) for a heater of width W = 2X (derivation in the lens write-up).
KT, A_WALL = 5.80e-6, 40.0


def vacuum_membrane(x, x_h):
    """Green's function of k·t·T'' = -δ(x - x_h)·(1 W/m), T(±A) = 0 (K per W/m)."""
    lo, hi = np.minimum(x, x_h), np.maximum(x, x_h)
    return (lo + A_WALL) * (A_WALL - hi) * 1e-6 / (KT * 2 * A_WALL)


wire = np.column_stack([vacuum_membrane(xe, xh) for xh in np.arange(-31, 32, 2.0)]).mean(axis=1)
r = fit_lens(PHASE_COEFF * wire[:, None], xe)
expect = 2 * KT * 64e-6 / (PHASE_COEFF * (32e-6) ** 2)
print(f"  [{ok(abs(r['power_per_rad'] / expect - 1) < 1e-3 and r['rms_residual_per_rad'] < 1e-6)}] "
      f"vacuum membrane, uniform 64 µm heater: {r['power_per_rad'] * 1e3:.4f} mW/rad "
      f"(analytic {expect * 1e3:.4f}), rms residual {r['rms_residual_per_rad']:.1e}")
lp = least_power_lens(PHASE_COEFF * wire[:, None], xe)
print(f"  [{ok(abs(lp / expect - 1) < 1e-3)}] same, least-power drive: {lp * 1e3:.4f} mW/rad "
      "(must equal the exact fit: one wire has nothing to trade)")
# (c) Wrapping can never cost more than not wrapping.
grid = np.linspace(0.1, 40, 400)
worse = [a for a in grid if isolated_power(a, xe, 1.0, True) > isolated_power(a, xe, 1.0, False) + 1e-12]
print(f"  [{ok(not worse)}] wrapped ≤ unwrapped for a = 0.1-40 rad")


# ---------------------------------------------------------------- Part 1: our membrane
def column_cache(kernel, positions):
    return {float(p): PHASE_COEFF * kernel(xe, float(p)) for p in positions}


def symmetric_layouts(n_heaters, x_max, step=1.0):
    """Heater sets symmetric about 0: a centre heater if K is odd, plus mirrored pairs ±x."""
    pair_x = np.arange(step, x_max + step / 2, step)
    n_pairs = n_heaters // 2
    for pairs in itertools.combinations(pair_x, n_pairs):
        x = np.array(sorted([-p for p in pairs] + [p for p in pairs] + ([0.0] if n_heaters % 2 else [])))
        if np.all(np.diff(x) >= 2 - 1e-9):
            yield x


def best_symmetric(cache, n_heaters, x_max):
    best = None
    for x in symmetric_layouts(n_heaters, x_max):
        r = fit_lens(np.column_stack([cache[float(v)] for v in x]), xe)
        key = (round(r["a_limit"], 2), -r["power_per_rad"])
        if best is None or key > best[0]:
            best = (key, x, r)
    return best[1], best[2]


rows = []


def report(stack, design, wires, x, r, G):
    """One table row. Power is the least-power drive within 2% (opa/lens.py); the
    least-squares drive's power is shown too, since the two can differ on a shared membrane."""
    lp = least_power_lens(G, xe) * 1e3
    p_ref = A_REF * lp
    where = "" if x is None else ", ".join(f"{v:g}" for v in x)
    a_txt = "  ≥60" if r["a_limit"] >= 59.99 else f"{r['a_limit']:5.1f}"
    f_txt = "   ≤35" if r["a_limit"] >= 59.99 else f"{r['f_min_um']:6.0f}"
    print(f"  {stack:9} {design:<28} {wires:>5}  {a_txt:>6}  {f_txt:>7}  {lp:7.3f}  "
          f"{r['power_per_rad'] * 1e3:7.3f}  {p_ref:7.2f}   {where}")
    rows.append([stack, design, wires, where, r["a_limit"], r["f_min_um"], lp, r["power_per_rad"] * 1e3, p_ref])


today = UndercutKernel()
far = SweptKernel("kernel_undercut_far.csv", "far", mirror=True, expected_columns=41)
standard = StandardKernel()
print("\nPart 1: lenses on the shared membrane (Strehl ≥ 0.8 at the focus)")
print(f"  {'stack':9} {'design':<28} {'wires':>5}  {'a_max':>6}  {'f_min':>7}  {'mW/rad':>7}  "
      f"{'(lsq)':>7}  {'P(1 mm)':>7}   heaters (µm)")
print(f"  {'':9} {'':<28} {'':>5}  {'(rad)':>6}  {'(µm)':>7}  {'':>7}  {'':>7}  {'(mW)':>7}")
best_far = {}
uniform_far = None
for stack, kernel, x_max in [("undercut", today, 32), ("far", far, 40)]:
    cache = column_cache(kernel, np.arange(-x_max, x_max + 0.5, 1.0))
    for n in (1, 2, 3, 4, 5, 6, 7):
        x, r = best_symmetric(cache, n, x_max)
        G = np.column_stack([cache[float(v)] for v in x])
        report(stack, f"K = {n}, symmetric", n, x, r, G)
        if stack == "far":
            best_far[n] = (x, G)
    # One wire, uniform power density over the array: 32 adjacent 2 µm heaters at equal power.
    G_uniform = np.column_stack([cache[float(v)] for v in np.arange(-31, 32, 2)]).mean(axis=1)[:, None]
    report(stack, "one uniform 64 µm heater", 1, None, fit_lens(G_uniform, xe), G_uniform)
    if stack == "far":
        uniform_far = G_uniform
    # One wire, best fixed shape: the power density over a 1 µm grid that fits the lens best.
    # A graded-width resistor would set this shape; only its overall strength is tunable.
    G_all = np.column_stack([cache[float(v)] for v in np.arange(-x_max, x_max + 0.5, 1.0)])
    shape = fit_lens(G_all, xe)["powers_per_rad"]
    G_shaped = (G_all @ shape / shape.sum())[:, None]
    report(stack, "one shaped heater (best fit)", 1, None, fit_lens(G_shaped, xe), G_shaped)
    # One independently driven heater per emitter.
    G33 = phase_matrix(kernel, xe, xe)
    report(stack, "one heater per emitter", 33, None, fit_lens(G33, xe), G33)

G_std = phase_matrix(standard, xe, xe)
report("standard", "one heater per emitter", 33, None, fit_lens(G_std, xe), G_std)
print("  (mW/rad: cheapest drive with every emitter within 2% of the lens; (lsq): the least-squares drive.\n"
      "   a_max ≥ 60 rad is past the paraxial range, so read it as 'not limited by the heat profile'.)")

with open(RESULTS / "lens_layouts.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["stack", "design", "wires", "heaters_um", "a_max_rad", "f_min_um", "mw_per_rad_least",
                "mw_per_rad_lsq", "P_f1mm_mW"])
    w.writerows(rows)

# ---------------------------------------------------------------- Part 2: conventional baselines
p_pi = {"undercut": np.pi / (PHASE_COEFF * today(0.0, 0.0)), "standard": np.pi / (PHASE_COEFF * standard(0.0, 0.0))}
print(f"\nPart 2: power (mW) against isolated heaters (P_π = {p_pi['undercut'] * 1e3:.3f} mW undercut, "
      f"{p_pi['standard'] * 1e3:.1f} mW standard)")
strengths = [1.0, A_REF, np.pi, 2 * np.pi, 4 * np.pi, 24.0]
print(f"  {'design':<40} " + " ".join(f"{a:6.2f}" for a in strengths) + "   (a, rad)")
print(f"  {'':<40} " + " ".join(f"{AF / a:6.0f}" for a in strengths) + "   (f, µm)")
comp = []


def line(name, powers_mw):
    print(f"  {name:<40} " + " ".join("   n/a" if np.isnan(p) else f"{p:6.1f}" for p in powers_mw))
    comp.append([name, *powers_mw])


for stack in ("undercut", "standard"):
    for wrapped in (False, True):
        line(f"isolated, {stack}, {'wrapped (Fresnel)' if wrapped else 'unwrapped'}",
             [isolated_power(a, xe, p_pi[stack], wrapped) * 1e3 for a in strengths])
# Best literature isolated suspended heater, for a harder baseline: Fang 2011, P_π = 0.49 mW [V].
line("isolated, best published (P_π 0.49 mW), wrapped", [isolated_power(a, xe, 0.49e-3, True) * 1e3 for a in strengths])
DUAL = np.array([-33.0, -16.0, 0.0, 16.0, 33.0])  # steers the full 22.8° and focuses (Part 3)
ours = [(f"ours, K = {n} (lens-best)", best_far[n][1]) for n in (1, 3, 5)]
ours += [("ours, one uniform 64 µm heater", uniform_far),
         ("ours, K = 5 dual-use [±33, ±16, 0]", np.column_stack([PHASE_COEFF * far(xe, v) for v in DUAL]))]
for name, G in ours:
    a_lim = fit_lens(G, xe)["a_limit"]
    per_rad = least_power_lens(G, xe)
    if np.isnan(per_rad):
        per_rad = fit_lens(G, xe)["power_per_rad"]  # K = 1 can't meet 2%: use its least-squares drive
        name += " (lsq)"
    line(name, [a * per_rad * 1e3 if a <= a_lim else np.nan for a in strengths])

with open(RESULTS / "lens_comparison.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["design", *[f"a={a:.3f} (f={AF / a:.0f} um)" for a in strengths]])
    w.writerows(comp)

# ---------------------------------------------------------------- Part 3: steer and focus together
print("\nPart 3: steering and focusing at once (far kernel). Largest lens strength (rad) usable at each angle")


def max_lens_at_angle(G, theta):
    """Largest a with Strehl ≥ 0.8 for the target ramp(θ) + a·lens, by bisection on a."""
    ramp = 2 * np.pi * PITCH_UM * np.sin(np.radians(theta)) / WAVELENGTH_UM * xe / PITCH_UM

    def strehl_at(a):
        target = ramp + a * unit_target(xe)
        P, c = fit_powers(G, target)
        return strehl(G @ P - target - c)

    if strehl_at(0.0) < STREHL_MIN:
        return np.nan
    lo, hi = 0.0, 60.0
    if strehl_at(hi) >= STREHL_MIN:
        return hi
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if strehl_at(mid) >= STREHL_MIN else (lo, mid)
    return lo


angles = [0, 2, 4, 6, 8, 10, 12, 14, 16, 20]
print(f"  {'layout':<44} " + " ".join(f"{t:>5}°" for t in angles))
layouts = [("lens-best K = 5", best_far[5][0]),
           ("steering-best K = 4 (R1)", np.array([-32.75, -6.7, 9.37, 32.59])),
           ("dual-use K = 5", DUAL)]
for name, x in layouts:
    G = np.column_stack([PHASE_COEFF * far(xe, v) for v in x])
    vals = [max_lens_at_angle(G, t) for t in angles]
    print(f"  {name + ' [' + ', '.join(f'{v:g}' for v in x) + ']':<44} "
          + " ".join("    --" if np.isnan(v) else f"{v:6.1f}" for v in vals))
print("  (-- : the angle itself is out of reach; 0.0 : the angle is reachable but no lens on top)")

# ---------------------------------------------------------------- Part 4: speed
print("\nPart 4: lens settling with a step drive (R4a step responses, heaters on its 4 µm grid)")
sr = StepResponses("transient_R4a.csv", mirror=True)
cache4 = {float(v): PHASE_COEFF * sr.final[:, k] for k, v in enumerate(sr.x_h)}
for n in (1, 3, 5):
    best = None
    for x in symmetric_layouts(n, 32, step=4.0):
        r = fit_lens(np.column_stack([cache4[float(v)] for v in x]), xe)
        key = (round(r["a_limit"], 2), -r["power_per_rad"])
        if best is None or key > best[0]:
            best = (key, x, r)
    _, x, r = best
    a = min(np.pi, 0.9 * r["a_limit"])
    idx = [int(np.flatnonzero(np.isclose(sr.x_h, v))[0]) for v in x]
    phase = PHASE_COEFF * np.einsum("tkn,k->tn", sr.S[:, idx], a * r["powers_per_rad"])
    target = a * unit_target(xe)
    # Strehl toward the ideal lens; the offset is free, so remove each time's best common phase.
    s_t = np.array([np.abs(np.exp(1j * (p - target)).mean()) ** 2 for p in phase])
    t_arrive = settle_time(sr.times, s_t, STREHL_MIN)
    err = phase[-1] - phase
    s_err = np.array([np.abs(np.exp(1j * e).mean()) ** 2 for e in err])
    print(f"  K = {n} at [{', '.join(f'{v:g}' for v in x)}], a = {a:.2f} rad (f = {AF / a:.0f} µm): "
          f"focus reaches Strehl 0.8 in {t_arrive * 1e3:.2f} ms; error vs final below the 0.8 level after "
          f"{settle_time(sr.times, s_err, STREHL_MIN) * 1e3:.2f} ms")
