"""Report and slide figures, one function per figure. Writes figures/<name>.pdf and .png.

Run all:            python make_figures.py
Run some:           python make_figures.py fig01 fig03

Every number plotted is read from the COMSOL exports or analysis/results/*.csv, or computed
here with the same opa/ functions the analyses use, so a figure cannot disagree with its
analysis. Literature points (fig05) are typed in, with their source beside them.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from opa import plotstyle as ps
from opa.constants import PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.conventional import ideal_power
from opa.farfield import element_pattern, far_field
from opa.lens import fit_lens, isolated_power, least_power_lens, lens_constant, unit_target
from opa.modes import half_width
from opa.steering import STREHL_MIN, fit_powers, steer_one_direction
from opa.thermal import (StandardKernel, SweptKernel, UndercutKernel, emitter_positions,
                         phase_matrix)

sys.stdout.reconfigure(encoding="utf-8")
ps.apply()
RESULTS = "results/"
xe = emitter_positions()

# Kernels and the four designs used throughout (heater positions in µm)
standard = StandardKernel()
undercut = UndercutKernel()
far = SweptKernel("kernel_undercut_far.csv", "far", True, 41)   # undercut, heaters out to the wall
slot = SweptKernel("kernel_slot.csv", "slot", False, 67)
vacuum = SweptKernel("kernel_vacuum.csv", "vacuum", False)
DESIGNS = {
    "undercut": (undercut, [-8, 10, 32]),
    "slot": (slot, [2, 14, 32]),
    "vacuum": (vacuum, [34]),
}


def dphi(theta_deg):
    """Phase step per emitter for a beam at theta_deg."""
    return 2 * np.pi * PITCH_UM * np.sin(np.radians(theta_deg)) / WAVELENGTH_UM


def drive(name, theta_deg):
    """Heater powers (W) that hold a design at theta_deg, from the usual non-negative fit."""
    kernel, heaters = DESIGNS[name]
    fit = steer_one_direction(phase_matrix(kernel, xe, heaters), xe, +1)
    return dphi(theta_deg) * fit["powers_per_rad"], fit


def temperature(kernel, heaters, powers, x):
    """Temperature rise (K) along x for 1 mm heaters: kernel (K per W/m) times W/m."""
    return sum(kernel(x, xh) * p / 1e-3 for xh, p in zip(heaters, powers))


# ====================================================================== fig01
STACKS = {
    "standard": "Standard stack: a heater on every emitter",
    "undercut": "Undercut: substrate removed under the array",
    "slot": "Slot: the near exit closed",
    "vacuum": "Vacuum: the gas removed",
}


def draw_profile(ax, kind):
    """Temperature across the array at 10°, normalised to its peak, above the cross-section."""
    if kind == "standard":
        # Each emitter needs only the wrapped phase, 0 to 2π: a sawtooth, one tooth per λ/sin θ.
        # (The conventional design's drives already absorb the crosstalk, so this is what the
        # cores reach.)
        T = np.mod(dphi(10) * xe / PITCH_UM, 2 * np.pi) / (2 * np.pi)
        ax.plot(xe, T, color=ps.COLOR[kind], lw=1.2, marker="o", ms=2.5)
        label = "phase each emitter needs at 10° (wrapped, 0–2π)"
    else:
        x = np.linspace(-50, 50, 1001)
        kernel, heaters = DESIGNS[kind]
        T = temperature(kernel, heaters, drive(kind, 10)[0], x)
        T = T / T.max()
        ax.fill_between(x, T, color=ps.COLOR[kind], alpha=0.18, lw=0)
        ax.plot(x, T, color=ps.COLOR[kind], lw=1.4)
        label = "temperature at 10° (one ramp, unwrapped)"
    ax.set_xlim(-50, 50)
    ax.set_ylim(0, 1.4)
    ax.set_yticks([])
    ax.grid(False)
    ax.spines[["left", "bottom"]].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)
    ax.text(-49, 1.38, label, fontsize=7, color=ps.INK2, va="top")


def arrow(ax, x1, y1, x2, y2, lw=1.2, color=ps.HEAT, ls="-", head=7):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls,
                                mutation_scale=head, shrinkA=0, shrinkB=0))


def draw_section(ax, kind, labels=True):
    """A cross-section of the stack (not to scale vertically for the heater and cores)."""
    m = ps.MAT
    ax.add_patch(Rectangle((-50, -20), 100, 18, color=m["si"], lw=0))            # substrate
    if kind != "standard":
        ax.add_patch(Rectangle((-40, -17), 80, 15, fc=m["gas"], ec=ps.AXIS, lw=0.6))  # cavity
    ax.add_patch(Rectangle((-50, -2), 100, 4.22, color=m["oxide"], lw=0))        # BOX + cladding
    for x in xe:
        ax.add_patch(Rectangle((x - 0.3, 0), 0.6, 0.55, color=m["core"], lw=0))
    if kind in ("slot", "vacuum"):
        ax.add_patch(Rectangle((36, -2), 2, 4.22, fc=m["gas"], ec="none", lw=0))  # the slot
    heaters = xe if kind == "standard" else DESIGNS[kind][1]
    for xh in heaters:
        ax.add_patch(Rectangle((xh - 0.8, 2.22), 1.6, 0.8, color=m["heater"], lw=0))

    y = -1.0  # arrows run through the membrane
    if kind == "standard":
        for xh in xe[::4]:
            arrow(ax, xh, 2.0, xh, -6.5, lw=1.0, head=6)
        if labels:
            ax.text(0, -9, "heat goes straight down, 2 µm to the substrate", ha="center",
                    fontsize=7.5, color=ps.INK)
    elif kind == "undercut":
        arrow(ax, 33.5, y, 40, y, lw=2.6)                       # short way out: most of it
        arrow(ax, 30.5, y, -40, y, lw=0.7)                      # the long way: little
        for xh in (-8, 10):                                     # interior heaters: both ways
            arrow(ax, xh - 1.2, 1.4, xh - 7, 1.4, lw=0.7, head=5)
            arrow(ax, xh + 1.2, 1.4, xh + 7, 1.4, lw=0.7, head=5)
        if labels:
            ax.text(27, -4.0, "≈90% of the edge heater's\nheat takes the short way", ha="center",
                    fontsize=7, color=ps.INK, va="top")
            ax.text(0, -12, "cavity (air)", ha="center", fontsize=7.5, color=ps.INK2)
    else:
        xh = DESIGNS[kind][1][-1]
        arrow(ax, xh - 1.5, y, -40, y, lw=2.6)                  # all of it crosses the array
        if kind == "slot":
            arrow(ax, 35.2, 0.6, 38.8, 0.6, lw=0.9, ls=(0, (2, 1.5)), head=5)          # across the slot
            for xl in (-20, 0, 20):
                arrow(ax, xl, -2.1, xl, -6.5, lw=0.9, ls=(0, (2, 1.5)), head=5)         # into the cavity gas
            if labels:
                ax.text(0, -12, "cavity (air): leaks across the slot and downward", ha="center",
                        fontsize=7.5, color=ps.INK2)
        elif labels:
            ax.text(0, -12, "cavity and slot evacuated: no gas path", ha="center", fontsize=7.5,
                    color=ps.INK2)
        if labels:
            ax.text(37, 3.6, "slot", ha="center", fontsize=7, color=ps.INK2)
    if labels and kind != "standard":
        sink = "wall: heat sink" if kind == "undercut" else "far wall: heat sink"
        ax.text(-45, -9.5, sink, rotation=90, ha="center", va="center", fontsize=6.5, color=ps.INK)
        if kind == "undercut":
            ax.text(45, -9.5, sink, rotation=90, ha="center", va="center", fontsize=6.5, color=ps.INK)
    ax.set_xlim(-50, 50)
    ax.set_ylim(-20, 4.5)
    ax.set_aspect("auto")
    ax.set_yticks([])
    ax.grid(False)
    ax.spines[["left"]].set_visible(False)
    ax.set_xticks([-40, -32, 0, 32, 40])
    ax.tick_params(axis="x", labelsize=7)


def fig01():
    """The four stacks: where the heat leaves, with the temperature it makes across the array."""
    fig = plt.figure(figsize=(ps.FULL, 9.0))
    # per stack: profile, cross-section, then an empty spacer row before the next stack
    gs = fig.add_gridspec(12, 1, height_ratios=[0.5, 1.15, 0.32] * 4, hspace=0.06)
    for i, kind in enumerate(STACKS):
        ax_p = fig.add_subplot(gs[3 * i])
        ax_s = fig.add_subplot(gs[3 * i + 1], sharex=ax_p)
        draw_section(ax_s, kind)
        draw_profile(ax_p, kind)
        ax_p.set_title(f"({'abcd'[i]})  {STACKS[kind]}", loc="left", fontsize=9, color=ps.INK, pad=2)
        if i < 3:
            ax_s.tick_params(labelbottom=False)
    ax_s.set_xlabel("x across the array (µm); emitters at −32 … 32, cavity walls at ±40")
    # material key under the last panel
    key = [("si", "silicon"), ("oxide", "oxide (BOX + cladding)"), ("core", "Si cores"),
           ("heater", "heater"), ("gas", "cavity / slot")]
    handles = [Rectangle((0, 0), 1, 1, fc=ps.MAT[k], ec=ps.AXIS, lw=0.5) for k, _ in key]
    handles.append(plt.Line2D([], [], color=ps.HEAT, lw=1.6))
    fig.legend(handles, [t for _, t in key] + ["heat flow (dashed: leaks)"], loc="lower center",
               ncol=6, bbox_to_anchor=(0.5, 0.0), fontsize=7, handlelength=1.2, columnspacing=1.0)
    fig.subplots_adjust(bottom=0.06, top=0.98)
    ps.save(fig, "fig01_stacks")
    # one panel per stack, for the slide sequence
    for kind in STACKS:
        fig = plt.figure(figsize=(ps.FULL, 2.5))
        gs = fig.add_gridspec(2, 1, height_ratios=[0.5, 1.15], hspace=0.08)
        ax_p, ax_s = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        draw_section(ax_s, kind)
        draw_profile(ax_p, kind)
        ax_p.set_title(STACKS[kind], loc="left", fontsize=9, color=ps.INK, pad=2)
        ax_s.set_xlabel("x across the array (µm)")
        ps.save(fig, f"fig01_{kind}")


# ====================================================================== fig02
def fig02():
    """One heater's temperature profile: how far heat spreads, and what the slot does to it."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(ps.FULL, 2.6))
    x = np.linspace(-60, 60, 2401)
    # labels sit in empty space with a thin leader to each curve's half-maximum point
    spots = {"standard": (-58, 0.95, "standard stack"), "undercut": (35, 0.95, "undercut")}
    for kind, kern in (("standard", standard), ("undercut", undercut)):
        T = kern(x, 0.0)
        a.plot(x, T / T.max(), color=ps.COLOR[kind])
        L = half_width(x, T)
        tx, ty, name = spots[kind]
        hx = -L if kind == "standard" else L
        a.annotate(f"{name}\nhalf-width {L:.1f} µm\npeak {T.max():.2f} K/mW", xy=(hx, 0.5),
                   xytext=(tx, ty), fontsize=7, color=ps.INK2, va="top",
                   arrowprops=dict(arrowstyle="-", color=ps.MUTED, lw=0.6, shrinkA=2, shrinkB=1))
        a.plot(hx, 0.5, "o", color=ps.COLOR[kind], ms=3.5)
    a.axvspan(-32, 32, color=ps.GRID, alpha=0.5, lw=0)
    a.set_xlim(-60, 60)
    a.set_ylim(0, 1.05)
    a.set_xlabel("x (µm)")
    a.set_ylabel("temperature rise / peak")
    a.text(0, 1.06, "array", ha="center", fontsize=7, color=ps.INK2)
    ps.panel_label(a, "a")
    a.set_title("Heater at the centre", fontsize=8.5, color=ps.INK, loc="left", pad=10)

    x = np.linspace(-40, 40, 1601)
    profiles = {"undercut": far(x, 34.0),
                "slot": np.interp(x, slot.x, slot.column(34)),
                "vacuum": np.interp(x, vacuum.x, vacuum.column(34))}
    labels = {"undercut": "undercut, no slot", "slot": "slot, in air", "vacuum": "slot, vacuum"}
    for kind, T in profiles.items():
        mask = x <= 36 if kind != "undercut" else np.ones_like(x, bool)
        b.plot(x[mask], T[mask], color=ps.COLOR[kind], label=labels[kind])
    b.axvspan(-32, 32, color=ps.GRID, alpha=0.5, lw=0)
    b.axvspan(36, 38, color=ps.AXIS, alpha=0.6, lw=0)
    b.text(37, b.get_ylim()[1] * 0.02, "slot", rotation=90, ha="center", va="bottom", fontsize=6.5,
           color=ps.INK2)
    b.set_xlim(-40, 40)
    b.set_ylim(0, None)
    b.set_xlabel("x (µm); heater at 34 µm, far wall at −40 µm")
    b.set_ylabel("temperature rise (K per mW)")
    b.legend(loc="upper left")
    ps.panel_label(b, "b")
    b.set_title("Heater at the array's edge", fontsize=8.5, color=ps.INK, loc="left", pad=10)
    fig.tight_layout(w_pad=2.5)
    ps.save(fig, "fig02_kernels")


# ====================================================================== fig03
def fig03():
    """How three heaters build the steering ramp, and the beams it makes."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(ps.FULL, 2.7), gridspec_kw={"width_ratios": [1, 1.15]})
    kernel, heaters = DESIGNS["slot"]
    P, fit = drive("slot", 10)
    x = np.linspace(-32, 32, 641)
    parts = [PHASE_COEFF * kernel(x, xh) * p for xh, p in zip(heaters, P)]
    shades = ps.BLUE_RAMP[:3]
    for xh, part, shade, p in zip(heaters, parts, shades, P):
        a.plot(x, part, color=shade, lw=1.2, label=f"heater at {xh} µm alone: {p * 1e3:.1f} mW")
    total = sum(parts)
    target = dphi(10) * x / PITCH_UM + dphi(10) * fit["offset_per_rad"]
    a.plot(x, total, color=ps.COLOR["slot"], lw=2.4, label="sum of the three")
    a.plot(x, target, color=ps.INK, lw=0.9, ls=(0, (3, 3)), label="ideal ramp for 10°")
    for xh in heaters:
        a.plot(xh, 0, marker="v", color=ps.MAT["heater"], ms=6, clip_on=False, zorder=5)
    a.set_xlim(-34, 34)
    a.set_ylim(0, 66)                    # headroom so the legend sits above the curves
    a.set_xlabel("x across the array (µm); ▼ heaters")
    a.set_ylabel("phase (rad)")
    handles, labels = a.get_legend_handles_labels()
    order = [3, 4, 0, 1, 2]                                   # sum and target first
    a.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper left", fontsize=6.8)
    ps.panel_label(a, "a")
    a.set_title("Slot design held at 10°", fontsize=8.5, color=ps.INK, loc="left", pad=10)

    ep = element_pattern()
    u = np.linspace(-0.9, 0.9, 7201)
    G = phase_matrix(kernel, xe, heaters)
    theta_u = np.degrees(np.arcsin(u))
    b.plot(theta_u, ep(u), color=ps.MUTED, lw=0.9, ls=(0, (4, 3)))
    b.text(-7, 1.04, "one emitter's pattern", fontsize=6.8, color=ps.INK2, ha="right")
    angles = [0, 5, 10, 15, 20]
    for th, shade in zip(angles, ps.BLUE_RAMP[1:]):
        I = far_field(G @ (dphi(th) * fit["powers_per_rad"]), xe, ep, u)
        b.plot(theta_u, I, color=shade, lw=1.2)
        if th % 10 == 0:                                   # label every other beam: no crowding
            k = np.argmax(I)
            b.text(theta_u[k], I[k] + 0.03, f"{th}°", ha="center", fontsize=7, color=ps.INK2)
    edge = np.degrees(np.arcsin(WAVELENGTH_UM / (2 * PITCH_UM)))
    for s in (-1, 1):
        b.axvline(s * edge, color=ps.AXIS, lw=0.8, ls=":")
    b.text(edge + 1.5, 0.42, f"lobe-free\nlimit ±{edge:.1f}°", fontsize=6.5, color=ps.INK2, va="top")
    b.text(-52, 0.86, "grating lobes", fontsize=6.8, color=ps.INK2)
    b.set_xlim(-55, 55)
    b.set_ylim(0, 1.1)
    b.set_xlabel("angle (°)")
    b.set_ylabel("far-field intensity")
    ps.panel_label(b, "b")
    b.set_title("Beams at 0–20° in 5° steps, same three heaters", fontsize=8.5, color=ps.INK,
                loc="left", pad=10)
    fig.tight_layout(w_pad=2.0)
    ps.save(fig, "fig03_ramp_and_beams")


# ====================================================================== fig04
P_PI_STANDARD = np.pi / (PHASE_COEFF * standard(0.0, 0.0))   # our standard-stack heater, 1 mm: 19.7 mW


def fig04():
    """Hold power against angle: our designs against a heater on every emitter."""
    fig, ax = plt.subplots(figsize=(ps.FULL, 3.0))
    th = np.linspace(0.05, 22.8, 300)
    # Heater per emitter, isolated, phases wrapped, at three phase-shifter efficiencies.
    # (The wiggles are real: the best wrap pattern changes with angle on a 33-emitter array.)
    refs = [(P_PI_STANDARD, "-", f"standard stack, P_π {P_PI_STANDARD * 1e3:.1f} mW (our model)"),
            (0.49e-3, (0, (4, 3)), "suspended, P_π 0.49 mW (Fang 2011)"),
            (0.05e-3, (0, (1, 2)), "folded record, P_π 0.05 mW (Lu 2015)")]
    for p_pi, ls, label in refs:
        P = np.array([ideal_power(xe, dphi(t), p_pi) for t in th]) * 1e3
        ax.plot(th, P, color=ps.MUTED, ls=ls, lw=1.0, label=label)
    ax.legend(loc="lower right", title="a heater on every emitter (33 wires):", title_fontsize=7,
              fontsize=6.8, alignment="left")
    names = {"undercut": "undercut, 3 heaters", "slot": "slot in air, 3 heaters", "vacuum": "vacuum, 1 heater"}
    for name, label in names.items():
        _, fit = drive(name, 1.0)
        t = th[th <= fit["range_deg"] + 1e-9]
        P = dphi(t) * fit["power_per_rad"] * 1e3
        ax.plot(t, P, color=ps.COLOR[name], lw=2.0)
        ax.plot(t[-1], P[-1], "o", color=ps.COLOR[name], ms=4)
        ax.text(t[-1] + 0.3, P[-1], label + ("" if t[-1] > 22 else f" (range {t[-1]:.1f}°)"),
                va="center", fontsize=7, color=ps.INK2)
    ax.set_yscale("log")
    ax.set_xlim(0, 25.5)
    ax.set_ylim(0.02, 3000)
    ax.set_xticks([0, 5, 10, 15, 20, 22.8])
    ax.set_xticklabels(["0", "5", "10", "15", "20", "22.8"])
    ax.set_xlabel("beam angle (°)")
    ax.set_ylabel("power to hold the angle (mW)")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ps.save(fig, "fig04_power_vs_angle")


# ====================================================================== fig05
# Published single heaters, plain step, as isolated heaters on 33 emitters at 10° (wrapped);
# time = the slower of rise and fall. Values as in the comparison doc (rev 58).
PUBLISHED = [  # label, settling ms, hold power at 10° mW, label offset (points)
    ("Watts 2013, doped Si", 0.0024, 373, (6, -10)),
    ("Jacques 2019, TiN", 0.0072, 860, (6, 0)),
    ("Miller 2020, multi-pass", 0.0065, 76.4, (6, -2)),
    ("Fang 2011, suspended", 0.144, 14.4, (-6, -10)),
    ("Wang 2020, suspended", 0.38, 14.1, (-4, 6)),
    ("Lu 2015, folded", 0.78, 1.47, (-6, -3)),
    ("Murray 2015, folded", 1.2, 2.79, (6, 2)),
    ("Parsons 2025, foundry undercut", 1.85, 31.5, (6, -8)),
]


def standard_settling_ms():
    """Our standard stack: one heater's own core reaching 92% of its final rise (transient.csv)."""
    from opa.comsol_io import read_comsol
    names, data = read_comsol(ps.FIG_DIR.parent / "transient.csv")
    times = np.array([float(n.split("t=")[1].split(";")[0]) for n in names[2:]])
    uc0 = np.array(["uc=0" in n for n in names[2:]])
    row = data[np.argmin(np.abs(data[:, 0]))]          # the node at x = 0, over the heater
    s = row[2:][uc0] / row[2:][uc0][-1]
    t = times[uc0]
    i = np.argmax(s >= 0.92)
    return 1e3 * np.interp(0.92, s[i - 1:i + 1], t[i - 1:i + 1])


def fig05():
    """The design path on power against settling time, beside published heaters."""
    pe = pd.read_csv(RESULTS + "preemphasis.csv")
    a = pe[(pe.section == "A") & (pe.start_deg == 0) & (pe.end_deg == 10)]
    steps = [("standard", "standard stack, 33 heaters", standard_settling_ms(),
              ideal_power(xe, dphi(10), P_PI_STANDARD) * 1e3, None)]
    for name, label in (("undercut", "undercut, 3 heaters"), ("slot", "slot in air, 3 heaters"),
                        ("vacuum", "vacuum, 1 heater")):
        d = a[a.design == name]
        plain = d[d["mode"] == "plain step"].iloc[0]
        fast = d[(d.cap_x_P10 == 3) & d["mode"].isin(["independent", "overdrive"])].iloc[0]
        steps.append((name, label, plain.settle_ms, plain.hold_mW, fast.settle_ms))

    fig, ax = plt.subplots(figsize=(ps.FULL, 3.6))
    for c, (lx, ly) in ((1, (1.3e-3, 520)), (10, (0.013, 1300)), (100, (0.13, 1300))):
        t = np.array([1e-3, 30])                           # lines of equal power × time
        ax.plot(t, c / t, color=ps.AXIS, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.text(lx, ly, f"{c} mW·ms", fontsize=6.5, color=ps.MUTED, ha="left", va="center")
    for label, t, p, off in PUBLISHED:
        ax.plot(t, p, "o", color=ps.MUTED, ms=4.5, zorder=3)
        ax.annotate(label, (t, p), xytext=off, textcoords="offset points", fontsize=6.5,
                    color=ps.INK2, ha="left" if off[0] > 0 else "right", va="center")
    for (n0, _, t0, p0, _), (n1, _, t1, p1, _) in zip(steps[:-1], steps[1:]):
        ax.annotate("", xy=(t1, p1), xytext=(t0, p0),
                    arrowprops=dict(arrowstyle="-|>", color=ps.INK2, lw=1.0, shrinkA=6, shrinkB=6,
                                    mutation_scale=8), zorder=2)
    for name, label, t, p, t_fast in steps:
        ax.plot(t, p, "o", color=ps.COLOR[name], ms=7, zorder=4)
        if t_fast:
            ax.plot([t_fast, t], [p, p], color=ps.COLOR[name], lw=1.6, alpha=0.5, zorder=2)
            ax.plot(t_fast, p, "o", mfc="white", mec=ps.COLOR[name], mew=1.6, ms=7, zorder=4)
        ax.annotate(f"{label}: {p:.3g} mW", (t, p), xytext=(7, 0), textcoords="offset points",
                    fontsize=7, color=ps.INK, va="center", fontweight="bold")
    ax.text(0.17, 63, "rings: same design with pre-emphasis", fontsize=6.8, color=ps.INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 30)
    ax.set_ylim(0.8, 2000)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("time to settle a 0 → 10° step (ms)")
    ax.set_ylabel("power to hold 10° (mW)")
    ps.save(fig, "fig05_power_vs_speed")


# ====================================================================== fig06
def fig06():
    """Pre-emphasis: the drive and the beam quality against time, plain step against shaped."""
    tr = pd.read_csv(RESULTS + "preemphasis_traces.csv")
    pe = pd.read_csv(RESULTS + "preemphasis.csv")
    a = pe[(pe.section == "A") & (pe.start_deg == 0) & (pe.end_deg == 10)]
    fig, axes = plt.subplots(2, 2, figsize=(ps.FULL, 3.6), sharex="col",
                             gridspec_kw={"height_ratios": [1, 1.15]})
    for col, (name, t_max) in enumerate((("slot", 4.5), ("vacuum", 11))):
        top, bot = axes[0, col], axes[1, col]
        for drive_name, color, lw, label in (("plain step", ps.MUTED, 1.2, "plain step"),
                                             ("pre-emphasis m=3", ps.COLOR[name], 1.8,
                                              "pre-emphasis (up to 3× hold power)")):
            d = tr[(tr.design == name) & (tr.drive == drive_name)]
            top.step(d.t_ms, d.total_power_mW, where="post", color=color, lw=lw, label=label)
            bot.plot(d.t_ms, d.strehl, color=color, lw=lw)
        d = a[a.design == name]
        t_plain = d[d["mode"] == "plain step"].settle_ms.iloc[0]
        t_fast = d[(d.cap_x_P10 == 3) & d["mode"].isin(["independent", "overdrive"])].settle_ms.iloc[0]
        bot.axhline(0.8, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
        for t, color, txt in ((t_fast, ps.COLOR[name], f"{t_fast:.2f} ms"), (t_plain, ps.MUTED, f"{t_plain:.2f} ms")):
            bot.axvline(t, color=color, lw=0.8, ls=":")
            bot.text(t + 0.04 * t_max, 0.35, txt, fontsize=7, color=ps.INK2)
        bot.text(t_max * 0.98, 0.82, "settled: Strehl ≥ 0.8", ha="right", va="bottom", fontsize=6.5,
                 color=ps.INK2)
        top.set_title(f"{'slot in air, 3 heaters' if name == 'slot' else 'vacuum, 1 heater'}: 0 → 10°",
                      loc="left", fontsize=8.5, color=ps.INK)
        top.set_xlim(0, t_max)
        top.set_ylim(0, None)
        bot.set_ylim(0, 1.05)
        bot.set_xlabel("time (ms)")
        top.legend(loc="upper right", fontsize=6.8)
    axes[0, 0].set_ylabel("chip power (mW)")
    axes[1, 0].set_ylabel("Strehl toward 10°")
    for ax, letter in zip(axes.flat, "abcd"):
        ps.panel_label(ax, letter)
    fig.tight_layout(h_pad=0.8, w_pad=2.0)
    ps.save(fig, "fig06_preemphasis")


# ====================================================================== fig07
def fig07():
    """Which phase patterns heat can make: singular values, and the slot design's cheapest modes."""
    modes = pd.read_csv(RESULTS + "modes.csv")
    shapes = pd.read_csv(RESULTS + "modes_shapes.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(ps.FULL, 2.7))
    cases = [("standard, d = 2 µm", "standard", "standard stack"),
             ("undercut", "undercut", "undercut"), ("slot", "slot", "slot in air"),
             ("vacuum", "vacuum", "vacuum")]
    idx = np.arange(1, 13)
    for case, color_key, label in cases:
        r = modes[modes.case == case].iloc[0]
        s = np.array([r[f"s{i}_rad_per_W"] for i in idx])
        a.plot(idx, s / s[0], "o-", color=ps.COLOR[color_key], ms=3.5, lw=1.3)
        a.text(12.4, s[-1] / s[0], f"{label}: {int(r.usable_10x)} usable", fontsize=6.5,
               color=ps.INK2, va="center")
    a.axhline(0.1, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
    a.text(5.3, 0.115, "usable: ≤ 10× the cheapest's power", fontsize=6.5, color=ps.INK2,
           va="bottom")
    a.set_yscale("log")
    a.set_ylim(5e-4, 1.5)
    a.set_xlim(0.5, 17.5)
    a.set_xticks([1, 3, 5, 7, 9, 11])
    a.set_xlabel("pattern (mode) number")
    a.set_ylabel("phase per watt / best pattern's")
    a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ps.panel_label(a, "a")
    a.set_title("Heat passes only the smooth patterns", loc="left", fontsize=8.5, color=ps.INK, pad=10)

    x = np.array([float(c[1:]) for c in shapes.columns if c.startswith("x")])
    slot_modes = shapes[shapes.case == "slot"].sort_values("mode")
    s1 = slot_modes.gain_rad_per_W.iloc[0]
    names = {1: "ramp", 2: "dome", 3: "ripple"}
    for (_, r), shade in zip(slot_modes.iloc[:3].iterrows(), ps.BLUE_RAMP[5:2:-1]):
        u = r[[c for c in shapes.columns if c.startswith("x")]].to_numpy(float)
        u = u * np.sign(u[-1])               # a mode's sign is arbitrary: draw each rising to the right
        b.plot(x, u, color=shade, lw=1.6)
        b.text(33.5, u[-1], f"{r['mode']}: {names[r['mode']]}, {s1 / r.gain_rad_per_W:.1f}×",
               fontsize=6.8, color=ps.INK2, va="center")
    b.axhline(0, color=ps.AXIS, lw=0.7)
    b.set_xlim(-34, 54)
    b.set_xticks([-32, -16, 0, 16, 32])
    b.set_xlabel("x across the array (µm)")
    b.set_ylabel("phase pattern (normalised)")
    b.text(0, b.get_ylim()[0] * 0.97, "×: power for the same\namount of pattern", fontsize=6.5,
           color=ps.INK2, va="bottom")
    ps.panel_label(b, "b")
    b.set_title("Slot design: its three cheapest patterns", loc="left", fontsize=8.5, color=ps.INK, pad=10)
    fig.tight_layout(w_pad=2.0)
    ps.save(fig, "fig07_modes")


# ====================================================================== fig08
def fig08():
    """Gas pressure: what it costs in power, and how far a leak moves the beam."""
    gas = pd.read_csv(RESULTS + "gas_loss.csv")
    faults = pd.read_csv(RESULTS + "faults.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(ps.FULL, 2.8))
    geoms = [(15, 2), (15, 10), (40, 10), (80, 10)]
    for (D, w), shade in zip(geoms, ps.BLUE_RAMP[2:]):
        g = gas[(gas.cavity_um == D) & (gas.slot_um == w)].sort_values("pressure_Pa")
        a.plot(g.pressure_Pa, g.mw_per_deg, "o-", color=shade, ms=3, lw=1.3,
               label=f"cavity {D} µm, slot {w} µm")
    limit = 0.506                                           # k·t·Δφ/(814·d), from notes (Q2)
    a.axhline(limit, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
    a.text(1.2, limit * 0.93, "conduction limit, 0.506 mW/deg", fontsize=6.5, color=ps.INK2, va="top")
    # COMSOL checks: R2 (air, 15 µm cavity, 2 µm slot) and R3a (vacuum), drawn at 1 Pa
    for p, v, txt, off, ha in ((101325, 1.383, "COMSOL R2, air", (-6, 4), "right"),
                               (1, 0.532, "COMSOL R3a, vacuum", (4, 6), "left")):
        a.plot([p], [v], "D", mfc="white", mec=ps.INK, ms=5, zorder=5)
        a.annotate(txt, (p, v), xytext=off, textcoords="offset points", ha=ha, fontsize=6.5,
                   color=ps.INK2)
    a.set_xscale("log")
    a.set_ylim(0.4, 1.6)
    a.set_xlabel("gas pressure (Pa)")
    a.set_ylabel("steering power (mW/deg)")
    a.legend(loc="upper left", fontsize=6.6)
    a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ps.panel_label(a, "a")
    a.set_title("Power: 1D model calibrated to COMSOL", loc="left", fontsize=8.5, color=ps.INK, pad=10)

    d = faults[(faults.section == "drift") & (faults.design == "vacuum, 1 heater")
               & faults.case.str.startswith("pressure") & (faults.quantity == "pointing_error_deg")]
    p = d.case.str.extract(r"pressure ([\d.e+]+) Pa")[0].astype(float)
    order = np.argsort(p.to_numpy())
    b.plot(p.to_numpy()[order], -d.value.to_numpy()[order], "o-", color=ps.COLOR["vacuum"], ms=3.5)
    b.axhline(0.12, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
    b.text(9e4, 0.13, "a tenth of the 1.2° beam", ha="right", fontsize=6.5, color=ps.INK2)
    b.set_xscale("log")
    b.set_yscale("log")
    b.set_xlabel("package pressure (Pa); drive set at 1 Pa")
    b.set_ylabel("beam falls short of 10° by (°)")
    for axis in (b.xaxis, b.yaxis):
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ps.panel_label(b, "b")
    b.set_title("Drift: vacuum design, if the package leaks", loc="left", fontsize=8.5, color=ps.INK,
                pad=10)
    fig.tight_layout(w_pad=2.0)
    ps.save(fig, "fig08_gas_pressure")


# ====================================================================== fig09
def fig09():
    """Fabrication errors: what the heaters can correct, against how smooth the errors are."""
    f = pd.read_csv(RESULTS + "fab_correlation.csv")
    f = f[f.angle_deg == 10]
    fig, ax = plt.subplots(figsize=(ps.FULL * 0.75, 2.9))
    ells = np.sort(f.ell_um.unique())
    xpos = np.where(ells == 0, 1.0, ells)                   # ℓ = 0 (independent errors) drawn at 1
    base = f[f.design == "slot, 3 heaters"].sort_values("ell_um")
    ax.plot(xpos, base.strehl_uncorrected, "o-", color=ps.MUTED, ms=3.5, lw=1.2,
            label="no correction")
    series = [("vacuum, 1 heater", ps.COLOR["vacuum"], "-"), ("slot, 3 heaters", ps.COLOR["slot"], "-"),
              ("dual-use, 5 heaters", ps.COLOR["dual"], "-"),
              ("slot, 17 even", ps.COLOR["slot"], (0, (3, 2)))]
    for name, color, ls in series:
        d = f[f.design == name].sort_values("ell_um")
        ax.plot(xpos, d.strehl_corrected, color=color, ls=ls, marker="o", ms=3.5, lw=1.5,
                label=f"heaters re-fitted: {name.replace('17 even', '17 heaters, evenly spaced')}")
    ax.axhline(0.8, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
    ax.text(1.05, 0.815, "Strehl 0.8", fontsize=6.5, color=ps.INK2)
    ax.set_xscale("log")
    ax.set_xticks(xpos)
    ax.set_xticklabels(["0"] + [f"{e:g}" for e in ells[1:]])
    ax.minorticks_off()
    ax.set_ylim(0.0, 1.02)
    ax.set_xlabel("correlation length of the path errors, ℓ (µm); 0 = independent")
    ax.set_ylabel("Strehl at 10°")
    ax.legend(loc="lower right", fontsize=6.6)
    ax.text(1.0, 1.0, "path errors σ = 1 rad", fontsize=6.8, color=ps.INK2, va="bottom")
    ps.save(fig, "fig09_fabrication")


# ====================================================================== fig10
LENS_AF = lens_constant(xe)        # a·f in rad·µm: a lens of strength a has focal length LENS_AF / a
DUAL = [-33.0, -16.0, 0.0, 16.0, 33.0]
# Lens layouts, all on the undercut membrane's far kernel (heaters out to the walls), as in run_lens.py:
# name: (heater positions in µm, colour, line style, label)
LENS_LAYOUTS = {
    "centre": ([0.0], ps.COLOR["undercut"], (0, (1, 1.4)), "one centre heater"),
    "uniform": (list(np.arange(-31.0, 32.0, 2.0)), ps.COLOR["undercut"], "-", "one 64 µm wide heater"),
    "dual": (DUAL, ps.COLOR["dual"], "-", "dual-use: 0, ±16, ±33 µm"),
    "steer4": ([-32.75, -6.7, 9.37, 32.59], ps.COLOR["undercut"], (0, (4, 2)),
               "steering-best: −32.8, −6.7, 9.4, 32.6 µm"),
}


def lens_matrix(name):
    """Phase matrix (rad per W) of a lens layout. The wide heater is one wire: 32 segments sharing its power."""
    G = phase_matrix(far, xe, LENS_LAYOUTS[name][0])
    return G.mean(axis=1, keepdims=True) if name == "uniform" else G


def lens_per_rad(G):
    """Power per rad of lens strength (W) and the strongest usable lens (rad), as run_lens.py reports them:
    the cheapest drive within 2% of the lens, or the least-squares drive where 2% can't be met."""
    r = fit_lens(G, xe)
    p = least_power_lens(G, xe)
    return (r["power_per_rad"] if np.isnan(p) else p), r["a_limit"]


def max_lens_at_angle(G, theta):
    """Strongest lens (rad) that still leaves Strehl ≥ 0.8 on a beam at theta (bisection, as run_lens.py)."""
    ramp = dphi(theta) * xe / PITCH_UM

    def strehl_at(a):
        target = ramp + a * unit_target(xe)
        P, c = fit_powers(G, target)
        return float(np.abs(np.exp(1j * (G @ P - target - c)).mean()) ** 2)

    if strehl_at(0.0) < STREHL_MIN:
        return np.nan
    lo, hi = 0.0, 60.0
    if strehl_at(hi) >= STREHL_MIN:
        return hi
    for _ in range(30):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if strehl_at(mid) >= STREHL_MIN else (lo, mid)
    return lo


def fig10():
    """The thermal lens: where the heaters go, why a wide heater makes a parabola, what it costs,
    and how much lens is left on top of a steered beam."""
    fig = plt.figure(figsize=(ps.FULL, 5.5))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.3], hspace=0.48)
    row1 = outer[0].subgridspec(1, 2, wspace=0.3)
    row2 = outer[1].subgridspec(1, 2, width_ratios=[1.75, 1], wspace=0.42)
    a, b = fig.add_subplot(row1[0]), fig.add_subplot(row1[1])
    c, d = fig.add_subplot(row2[0]), fig.add_subplot(row2[1])

    # (a) the layouts, seen from above along the array axis
    for i, name in enumerate(LENS_LAYOUTS):
        x_h, color, ls, label = LENS_LAYOUTS[name]
        y = -i
        a.add_patch(Rectangle((-40, y - 0.13), 80, 0.26, color=ps.MAT["oxide"], lw=0))   # membrane
        for left in (-41.2, 40):                                                         # cavity walls
            a.add_patch(Rectangle((left, y - 0.3), 1.2, 0.6, color=ps.MAT["si"], lw=0))
        a.plot(xe, np.full_like(xe, y), "|", color=ps.MAT["core"], ms=4, mew=0.7)       # emitters
        if name == "uniform":
            a.add_patch(Rectangle((-32, y - 0.09), 64, 0.18, color=ps.MAT["heater"], lw=0, zorder=3))
        else:
            for xh in x_h:
                a.add_patch(Rectangle((xh - 0.8, y - 0.09), 1.6, 0.18, color=ps.MAT["heater"], lw=0,
                                      zorder=3))
        a.plot([-40, -33], [y + 0.42] * 2, color=color, ls=ls, lw=1.6)                   # line key
        wires = 1 if name in ("centre", "uniform") else len(x_h)
        a.text(-31, y + 0.42, f"{label}  ·  {wires} wire{'s' if wires > 1 else ''}", fontsize=6.6,
               color=ps.INK, va="center")
    a.set_xlim(-42, 42)
    a.set_ylim(-3.45, 0.75)
    a.set_yticks([])
    a.grid(False)
    a.spines[["left"]].set_visible(False)
    a.set_xticks([-32, -16, 0, 16, 32])
    a.set_xlabel("position (µm); cavity walls at ±40 µm")
    ps.panel_label(a, "a")
    a.set_title("Heater layouts (gold) under the 33 emitters", loc="left", fontsize=8.5, color=ps.INK, pad=10)

    # (b) a wide heater makes a parabola, a point heater a tent: k·t·T'' = -q (vacuum limit)
    x = np.linspace(-40, 40, 801)
    X = xe.max()
    inside = np.abs(x) <= X
    b.plot(x[inside], 1 - (x[inside] / X) ** 2, color=ps.INK, lw=0.9, ls=(0, (3, 2)), zorder=5)
    for name in ("centre", "uniform"):
        x_h, color, ls, _ = LENS_LAYOUTS[name]
        G = lens_matrix(name)
        P, off = fit_powers(G, unit_target(xe))
        cols = np.column_stack([PHASE_COEFF * far(x, v) for v in x_h])
        phi = cols.mean(axis=1) * P.sum() if name == "uniform" else cols @ P
        b.plot(x, phi - off + 1, color=color, ls=ls, lw=1.6)
    for edge in (-X, X):
        b.axvline(edge, color=ps.AXIS, lw=0.7)
    b.annotate("one centre heater:\na tent", (3, 1.12), xytext=(10, 1.2), fontsize=6.6, color=ps.INK2,
               va="center", arrowprops=dict(arrowstyle="-", color=ps.INK2, lw=0.5))
    b.annotate("64 µm heater:\non the parabola", (-20, 1 - (20 / X) ** 2), xytext=(-39, 1.05),
               fontsize=6.6, color=ps.INK2, va="center",
               arrowprops=dict(arrowstyle="-", color=ps.INK2, lw=0.5))
    b.text(0, 0.42, "target lens\n(dashed)", ha="center", fontsize=6.6, color=ps.INK2)
    b.set_xlim(-40, 40)
    b.set_ylim(-0.35, 1.32)
    b.set_xticks([-32, -16, 0, 16, 32])
    b.set_xlabel("position (µm)")
    b.set_ylabel("phase lead ÷ lens strength")
    ps.panel_label(b, "b")
    b.set_title("Why spread heat makes a lens", loc="left", fontsize=8.5, color=ps.INK, pad=10)

    # (c) power against focal length (stronger lens to the right), against a heater on every emitter.
    # Computed in lens strength a, as in run_lens.py, and plotted at f = a·f / a.
    F = lambda strength: LENS_AF / 1000 / np.asarray(strength)          # focal length, mm
    a_grid = np.geomspace(0.5, 40, 400)
    c.set_xscale("log")
    c.set_yscale("log")
    c.set_xlim(F(0.5), F(40))                                             # decreasing: stronger lens →
    c.set_ylim(0.2, 1500)
    p_pi = {"standard": np.pi / (PHASE_COEFF * standard(0.0, 0.0)),
            "undercut": np.pi / (PHASE_COEFF * undercut(0.0, 0.0)), "fang": 0.49e-3}
    # Grey: phases wrapped (a Fresnel lens), so the cost stops growing past a ≈ 2π; labelled on the flat part
    refs = [("standard", "-", "standard stack", (5.5, 700)),
            ("undercut", (0, (4, 2)), "undercut membrane", (2.3, 45)),
            ("fang", (0, (1, 1.4)), "Fang 2011", (3.4, 15.5))]
    for key, ls, txt, (tx, ty) in refs:
        P = [isolated_power(v, xe, p_pi[key], True) * 1e3 for v in a_grid]
        c.plot(F(a_grid), P, color=ps.MUTED, ls=ls, lw=1.2)
        c.text(F(tx), ty, f"{txt}, P$_\\pi$ {p_pi[key] * 1e3:.3g} mW", fontsize=6.2, color=ps.INK2)
    check, fits = [], {}
    for name in ("uniform", "dual", "centre"):
        _, color, ls, label = LENS_LAYOUTS[name]
        per_rad, a_lim = fits[name] = lens_per_rad(lens_matrix(name))
        check.append(f"{name} {A_REF_LENS * per_rad * 1e3:.3f} mW at f = 1 mm")
        aa = a_grid[a_grid <= a_lim]
        c.plot(F(aa), aa * per_rad * 1e3, color=color, ls=ls, lw=1.8)
        if a_lim < 40:
            c.plot(F(a_lim), [a_lim * per_rad * 1e3], "o", color=color, ms=3.5)
    print("   check against lens_comparison.csv (1.161 / 2.410 / 1.782):", "; ".join(check))
    # Ours: power ∝ a, so every one of our lines has the same slope; the two outer ones are labelled
    # along their line, the middle one (dual-use) by a leader into the empty lower right.
    p0, p1 = c.transData.transform([(F(1), 1.0), (F(10), 10.0)])
    slope = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
    on_line = dict(rotation=slope, rotation_mode="anchor", transform_rotates_text=False, fontsize=6.2,
                   color=ps.INK)
    per_u, per_c, a_c = fits["uniform"][0] * 1e3, fits["centre"][0] * 1e3, fits["centre"][1]
    c.text(F(0.58), 0.58 * per_u * 1.25, "ours: one 64 µm heater, 1 wire", va="bottom", **on_line)
    c.text(F(0.85), 0.85 * per_c / 1.2, "ours: one centre heater, 1 wire", va="top", **on_line)
    c.text(F(a_c * 1.1), a_c * per_c / 1.05, f"f = {F(a_c):.2f} mm", va="top", fontsize=6.2,
           color=ps.INK2)
    c.annotate(f"ours: dual-use, 5 wires, also steers;\nusable to f = {LENS_AF / fits['dual'][1]:.0f} µm",
               (F(13.5), 13.5 * fits["dual"][0] * 1e3), xytext=(F(15), 0.24), fontsize=6.2, color=ps.INK,
               ha="center", va="bottom", arrowprops=dict(arrowstyle="-", color=ps.INK2, lw=0.5, shrinkB=1))
    c.axvspan(0.1, F(40), color=ps.GRID, alpha=0.55, lw=0)
    c.text(0.096, 260, "f < 0.1 mm:\nparaxial only", fontsize=6.0, color=ps.INK2, va="top")
    fmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")
    c.set_xticks([4, 2, 1, 0.5, 0.2, 0.1])
    c.xaxis.set_major_formatter(fmt)
    c.yaxis.set_major_formatter(fmt)
    c.minorticks_off()
    c.set_xlabel("focal length f (mm); stronger lens →")
    c.set_ylabel("power to hold the lens (mW)")
    ps.panel_label(c, "c")
    c.set_title("Power: ours against a heater on every emitter (grey)", loc="left", fontsize=8.5,
                color=ps.INK, pad=10)

    # (d) steer and focus at once: the strongest lens left on top of a steered beam
    angles = np.linspace(0, 22.8, 39)
    for name in ("dual", "steer4"):
        _, color, ls, _ = LENS_LAYOUTS[name]
        G = lens_matrix(name)
        vals = np.array([max_lens_at_angle(G, t) for t in angles])
        d.plot(angles, vals, color=color, ls=ls, lw=1.8)
        print(f"   {name}: lens at 0/10/20° = {max_lens_at_angle(G, 0):.1f} / {max_lens_at_angle(G, 10):.1f}"
              f" / {max_lens_at_angle(G, 20):.1f} rad")
    for f_mm in (1.0, 0.2, 0.1):
        d.axhline(LENS_AF / 1000 / f_mm, color=ps.AXIS, lw=0.7, ls=(0, (3, 3)))
    # The same quantity in focal length on the right (a unit conversion, not a second measure)
    to_f = lambda v: LENS_AF / 1000 / np.clip(v, 1e-9, None)
    right = d.secondary_yaxis("right", functions=(to_f, to_f))
    right.set_yticks([1, 0.2, 0.1])
    right.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    right.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    right.set_ylabel("focal length (mm)", fontsize=7.5)
    d.text(0.5, 23.9, "dual-use, 5 wires", fontsize=6.4, color=ps.INK)
    d.text(0.5, 5.9, "steering-best, 4 wires", fontsize=6.4, color=ps.INK, va="top")
    d.text(12.3, 14.4, "lens-best layouts (no heaters\nat the ends) can't steer", ha="center", va="center",
           fontsize=6.2, color=ps.INK2)
    d.set_xlim(0, 23.5)
    d.set_ylim(0, 26.5)
    d.set_xlabel("steering angle (°)")
    d.set_ylabel("strongest lens on top (rad)")
    ps.panel_label(d, "d")
    d.set_title("Steer and focus, same wires", loc="left", fontsize=8.5, color=ps.INK, pad=10)
    ps.save(fig, "fig10_lens")


A_REF_LENS = LENS_AF / 1000.0      # lens strength for f = 1 mm (≈ 2.08 rad)


# ====================================================================== fig11
def fig11():
    """Faults: what one dead heater costs, and how big each slow disturbance is against the beam."""
    faults = pd.read_csv(RESULTS + "faults.csv")

    def get(section, design, case, quantity):
        row = faults[(faults.section == section) & (faults.design == design) & (faults.case == case)
                     & (faults.quantity == quantity)]
        return float(row.value.iloc[0])

    fig = plt.figure(figsize=(ps.FULL, 6.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.42)
    a, b = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # (a) steering range left after one heater dies and the rest are refitted
    intact = {"vacuum": drive("vacuum", 10)[1]["range_deg"], "slot": drive("slot", 10)[1]["range_deg"],
              "dual": steer_one_direction(phase_matrix(far, xe, DUAL), xe, +1)["range_deg"]}
    groups = [
        ("vacuum", "Vacuum: one heater (34 µm)", "vacuum, 1 heater",
         [("heater at 34 µm", "34 µm fails"),
          ("spare heater at 32 µm", "34 µm fails, a spare at 32 µm takes over")]),
        ("slot", "Slot: heaters at 2, 14, 32 µm", "slot, 3 heaters",
         [("heater at 2 µm", "2 µm fails"), ("heater at 14 µm", "14 µm fails"),
          ("heater at 32 µm", "32 µm (edge) fails"),
          ("spare at 34 µm replacing 32 µm", "32 µm fails, a spare at 34 µm takes over")]),
        ("dual", "Dual-use: heaters at 0, ±16, ±33 µm", "dual-use, 5 heaters",
         [("heater at -33 µm", "−33 µm fails"), ("heater at -16 µm", "−16 µm fails"),
          ("heater at 0 µm", "0 µm fails"), ("heater at 16 µm", "+16 µm fails"),
          ("heater at 33 µm", "+33 µm (edge) fails")]),
    ]
    y, ticks, labels, bold = 0.0, [], [], []
    for key, title, design, cases in groups:
        ticks.append(y)
        labels.append(title)
        bold.append(True)
        y -= 1
        for case, text in [(None, "all working")] + cases:
            v = intact[key] if case is None else get("dead heater", design, case, "range_after_deg")
            a.barh(y, v, height=0.62, color=ps.COLOR[key], lw=0)
            a.text(max(v, 22.8) + 0.3, y, f"{v:.1f}°" if v > 0.05 else "0°", va="center", fontsize=6.6,
                   color=ps.INK2)
            ticks.append(y)
            labels.append(text)
            bold.append(False)
            y -= 1
        y -= 0.4
    a.set_yticks(ticks)
    a.set_yticklabels(labels, fontsize=6.9)
    for t, is_bold in zip(a.get_yticklabels(), bold):
        if is_bold:
            t.set_fontweight("bold")
            t.set_color(ps.INK)
    a.tick_params(axis="y", length=0)
    a.axvline(22.8, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
    a.text(22.8, 0.75, "lobe-free limit, 22.8°", ha="right", fontsize=6.5, color=ps.INK2)
    a.set_xlim(0, 25.5)
    a.set_ylim(y + 0.6, 1.2)
    a.grid(axis="y", visible=False)
    a.set_xlabel("steering range after the others are refitted (°, toward +θ)")
    ps.panel_label(a, "a", x=-0.38)
    a.set_title("One heater dies: the edge heater is the single point of failure", loc="left",
                fontsize=8.5, color=ps.INK, pad=10, x=-0.36)

    # (b) slow disturbances: pointing error at 10°, against a tenth of the beam width
    ours, conv = ps.INK2, ps.MUTED
    dn = lambda kind, d: get("drive noise", d, kind, "pointing_rms_deg")
    designs4 = ["vacuum, 1 heater", "slot, 3 heaters", "undercut, 3 heaters", "dual-use, 5 heaters"]
    ind = [dn("independent", d) for d in designs4]
    com = [dn("common", d) for d in designs4]
    amb = [abs(get("drift", "all designs", f"ambient +10 K, b = {s}", "pointing_error_deg"))
           for s in ("0.0e+00", "7.0e-04", "1.4e-03")]
    vac = {p: abs(get("drift", "vacuum, 1 heater", f"pressure {p} Pa", "pointing_error_deg"))
           for p in (10, 30, 100)}
    weather = [abs(get("drift", "slot in air, 3 heaters", f"pressure {p} atm", "pointing_error_deg"))
               for p in ("0.95", "1.05")]
    rows = [  # (label, [(lo, hi, colour)], [(value, colour, marker, note, note above or right)])
        ("Drive noise, 1% rms", None, None),
        ("each heater's driver on its own", [(min(ind), max(ind), ours)],
         [(dn("independent", "conventional"), conv, "D", "conventional, 33 heaters", "above")]),
        ("shared supply, all heaters together", [(min(com), max(com), ours)],
         [(dn("common", "conventional"), conv, "D", "conventional", "above")]),
        ("Ambient +10 K", None, None),
        ("membrane k(T), dn/dT(T); b unverified", [(min(amb), max(amb), ours, "thin")],
         [(v, ours, "o", t, where) for v, t, where in
          zip(amb, ("0", "b = 0.7 (10⁻³/K)", "1.4"), ("above-l", "above", "above-r"))]),
        ("gas conductivity, slot in air", [],
         [(abs(get("drift", "slot in air, 3 heaters", "ambient +10 K, gas only", "pointing_error_deg")),
           ps.COLOR["slot"], "o", "", "")]),
        ("Package pressure, drive set beforehand", None, None),
        ("slot in air, weather ±5%", [(min(weather), max(weather), ps.COLOR["slot"])], []),
        ("vacuum, 1 Pa rises to 10 / 30 / 100 Pa", [],
         [(v, ps.COLOR["vacuum"], "o", f"{p} Pa", "right") for p, v in vac.items()]),
        ("Left out of the model", None, None),
        ("radiation in vacuum (ε = 1, both faces)", [],
         [(abs(get("drift", "vacuum, 1 heater", "radiation, emissivity 1", "pointing_error_deg")),
           ps.COLOR["vacuum"], "o", "", "")]),
    ]
    y, ticks, labels, bold = 0.0, [], [], []
    for label, spans, points in rows:
        ticks.append(y)
        labels.append(label)
        bold.append(spans is None)
        for lo, hi, color, *style in spans or []:
            b.plot([lo, hi], [y, y], color=color, lw=1.0 if style else 4.5, solid_capstyle="round")
        for v, color, marker, note, where in points or []:
            b.plot([v], [y], marker, color=color, ms=5 if marker == "o" else 4.5, mec="white", mew=0.6,
                   zorder=4)
            if where.startswith("above"):
                ha = {"above": "center", "above-l": "right", "above-r": "left"}[where]
                b.text(v, y + 0.32, note, ha=ha, fontsize=6.0, color=ps.INK2)
            elif where == "right":
                b.text(v * 1.13, y, note, va="center", fontsize=6.0, color=ps.INK2)
        y -= 1
    b.text(min(ind) / 1.15, ticks[1] + 0.32, "our 4 designs", ha="right", fontsize=6.0, color=ps.INK2)
    b.axvline(0.12, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
    b.text(0.13, 0.55, "a tenth of the 1.2° beam", fontsize=6.5, color=ps.INK2)
    b.set_xscale("log")
    b.set_xlim(6e-4, 1.0)
    b.set_ylim(y + 0.5, 1.0)
    b.set_yticks(ticks)
    b.set_yticklabels(labels, fontsize=6.9)
    for t, is_bold in zip(b.get_yticklabels(), bold):
        if is_bold:
            t.set_fontweight("bold")
            t.set_color(ps.INK)
    b.tick_params(axis="y", length=0)
    b.grid(axis="y", visible=False)
    b.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.minorticks_off()
    b.set_xlabel("pointing error at 10° (°, log scale)")
    ps.panel_label(b, "b", x=-0.38)
    b.set_title("Slow disturbances: the angle is only as good as the gain", loc="left", fontsize=8.5,
                color=ps.INK, pad=10, x=-0.36)
    ps.save(fig, "fig11_faults")


FIGURES = {f.__name__: f for f in (fig01, fig02, fig03, fig04, fig05, fig06, fig07, fig08, fig09,
                                   fig10, fig11)}

if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGURES):
        print(name)
        FIGURES[name]()
