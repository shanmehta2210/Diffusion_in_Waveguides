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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

from opa import plotstyle as ps
from opa.constants import PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.conventional import ideal_power
from opa.farfield import element_pattern, far_field, near_field_intensity
from opa.lens import fit_lens, isolated_power, least_power_lens, lens_constant, unit_target
from opa.modes import half_width
from opa.steering import STREHL_MIN, fit_powers, steer_one_direction, strehl_vs_step
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
    b.text(edge + 1.5, 0.42, f"one beam only\ninside ±{edge:.1f}°", fontsize=6.5, color=ps.INK2, va="top")
    b.text(-53, 0.72, "copies of\neach beam\n(grating\nlobes)", fontsize=6.5, color=ps.INK2, va="bottom")
    b.set_xlim(-55, 55)
    b.set_ylim(0, 1.1)
    b.set_xlabel("angle (°)")
    b.set_ylabel("far-field intensity")
    ps.panel_label(b, "b")
    # Every beam uses the same three heaters: the 10° powers times sin θ / sin 10° (0, ½, 1, 1.5, 2).
    b.set_title("Same powers × 0, ½, 1, 1.5, 2: beams at 0–20°", fontsize=8.5, color=ps.INK,
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


def design_steps():
    """The design path at 10°: (name, label, plain-step settling ms, hold mW, pre-emphasised settling ms)."""
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
    return steps


def fig05():
    """The design path on power against settling time, beside published heaters."""
    steps = design_steps()

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


def draw_lens_shapes(ax):
    """A wide heater makes a parabola, a point heater a tent: k·t·T'' = -q (vacuum limit).
    The phase each makes for a unit lens, against the target, with the heaters drawn in gold."""
    x = np.linspace(-40, 40, 801)
    X = xe.max()
    inside = np.abs(x) <= X
    ax.plot(x[inside], 1 - (x[inside] / X) ** 2, color=ps.INK, lw=0.9, ls=(0, (3, 2)), zorder=5)
    for name in ("centre", "uniform"):
        x_h, color, ls, _ = LENS_LAYOUTS[name]
        G = lens_matrix(name)
        P, off = fit_powers(G, unit_target(xe))
        cols = np.column_stack([PHASE_COEFF * far(x, v) for v in x_h])
        phi = cols.mean(axis=1) * P.sum() if name == "uniform" else cols @ P
        ax.plot(x, phi - off + 1, color=color, ls=ls, lw=1.6)
    for edge in (-X, X):
        ax.axvline(edge, color=ps.AXIS, lw=0.7)
    # The heaters themselves, under the curves: the 64 µm sheet and the single centre strip
    ax.add_patch(Rectangle((-32, -0.36), 64, 0.07, color=ps.MAT["heater"], lw=0))
    ax.add_patch(Rectangle((-0.8, -0.22), 1.6, 0.07, color=ps.MAT["heater"], lw=0))
    ax.text(-30, -0.27, "64 µm heater", fontsize=6.2, color=ps.INK2, va="bottom")
    ax.text(1.8, -0.185, "centre heater", fontsize=6.2, color=ps.INK2, va="center")
    ax.annotate("one centre heater:\na tent", (3, 1.12), xytext=(10, 1.2), fontsize=6.6, color=ps.INK2,
                va="center", arrowprops=dict(arrowstyle="-", color=ps.INK2, lw=0.5))
    ax.annotate("64 µm heater:\non the parabola", (-20, 1 - (20 / X) ** 2), xytext=(-39, 1.05),
                fontsize=6.6, color=ps.INK2, va="center",
                arrowprops=dict(arrowstyle="-", color=ps.INK2, lw=0.5))
    ax.text(0, 0.42, "target lens\n(dashed)", ha="center", fontsize=6.6, color=ps.INK2)
    ax.set_xlim(-40, 40)
    ax.set_ylim(-0.42, 1.32)
    ax.set_xticks([-32, -16, 0, 16, 32])
    ax.set_xlabel("position (µm); heaters in gold")
    ax.set_ylabel("phase lead ÷ lens strength")


def draw_lens_power(c, compact=False):
    """Power to hold a lens against focal length: ours (coloured) against a heater on every
    emitter with wrapped phases (grey). Shared by fig10 and the report's lens figure; compact=True
    (the narrower report panel) moves the grey labels into a legend and shortens the others."""
    # (c) power against focal length (stronger lens to the right), against a heater on every emitter.
    # Computed in lens strength a, as in run_lens.py, and plotted at f = a·f / a.
    F = lambda strength: LENS_AF / 1000 / np.asarray(strength)          # focal length, mm
    a_grid = np.geomspace(0.5, 40, 400)
    c.set_xscale("log")
    c.set_yscale("log")
    c.set_xlim(F(0.5), F(40))                                             # decreasing: stronger lens →
    c.set_ylim(0.2, 8000 if compact else 1500)            # compact: headroom for the legend
    p_pi = {"standard": np.pi / (PHASE_COEFF * standard(0.0, 0.0)),
            "undercut": np.pi / (PHASE_COEFF * undercut(0.0, 0.0)), "fang": 0.49e-3}
    # Grey: phases wrapped (a Fresnel lens), so the cost stops growing past a ≈ 2π; labelled on the flat part
    refs = [("standard", "-", "standard stack", (5.5, 700)),
            ("undercut", (0, (4, 2)), "undercut membrane", (2.3, 45)),
            ("fang", (0, (1, 1.4)), "Fang 2011", (3.4, 15.5))]
    for key, ls, txt, (tx, ty) in refs:
        P = [isolated_power(v, xe, p_pi[key], True) * 1e3 for v in a_grid]
        label = f"{txt}, P$_\\pi$ {p_pi[key] * 1e3:.3g} mW"
        c.plot(F(a_grid), P, color=ps.MUTED, ls=ls, lw=1.2, label=label)
        if not compact:
            c.text(F(tx), ty, label, fontsize=6.2, color=ps.INK2)
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
    leader = dict(arrowstyle="-", color=ps.INK2, lw=0.5, shrinkB=1)
    if compact:
        c.legend(loc="upper left", fontsize=5.6, handlelength=2.2, frameon=False,
                 title="grey: one heater per emitter, wrapped", title_fontsize=5.6)
        c.annotate("ours: dual-use, 5 wires", (F(13.5), 13.5 * fits["dual"][0] * 1e3), xytext=(F(13), 0.3),
                   fontsize=6.0, color=ps.INK, ha="center", va="bottom", arrowprops=leader)
    else:
        c.text(F(a_c * 1.1), a_c * per_c / 1.05, f"f = {F(a_c):.2f} mm", va="top", fontsize=6.2,
               color=ps.INK2)
        c.annotate(f"ours: dual-use, 5 wires, also steers;\nusable to f = {LENS_AF / fits['dual'][1]:.0f} µm",
                   (F(13.5), 13.5 * fits["dual"][0] * 1e3), xytext=(F(15), 0.24), fontsize=6.2,
                   color=ps.INK, ha="center", va="bottom", arrowprops=leader)
    c.axvspan(0.1, F(40), color=ps.GRID, alpha=0.55, lw=0)
    if not compact:
        c.text(0.096, 260, "f < 0.1 mm:\nparaxial only", fontsize=6.0, color=ps.INK2, va="top")
    fmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")
    c.set_xticks([4, 2, 1, 0.5, 0.2, 0.1])
    c.xaxis.set_major_formatter(fmt)
    c.yaxis.set_major_formatter(fmt)
    c.minorticks_off()
    c.set_xlabel("focal length f (mm); stronger lens →")
    c.set_ylabel("power to hold the lens (mW)")



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

    draw_lens_shapes(b)
    ps.panel_label(b, "b")
    b.set_title("Why spread heat makes a lens", loc="left", fontsize=8.5, color=ps.INK, pad=10)

    draw_lens_power(c)
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


# ====================================================================== fig12 (slide schematic)
def s_curve(x0, x1, y0, y1, n=40):
    """A smooth S-bend from (x0, y0) to (x1, y1), for drawing waveguide fan-outs."""
    t = np.linspace(0, 1, n)
    return x0 + (x1 - x0) * (1 - np.cos(np.pi * t)) / 2, y0 + (y1 - y0) * t


def draw_top_view(ax, kind):
    """Top view of the chip, not to scale: splitter, heated section, gratings, wires and pads."""
    core, wire = ps.MAT["core"], ps.MUTED
    half = 2.2                                   # the 64 µm array, drawn 4.4 units wide
    scale = half / 32                            # drawing units per µm across the array
    narrow = 5 + xe * scale
    y_in, y_split, y_h0, y_h1, y_g0, y_g1 = 0.3, 1.6, 3.7, 6.3, 8.0, 9.2
    ax.plot([5, 5], [y_in, 0.9], color=core, lw=0.9)
    for xn in narrow:                                                       # 1 → 33 splitter
        ax.plot([5, xn], [0.9, y_split], color=core, lw=0.3)
    if kind == "conventional":
        wide = 5 + np.linspace(-3.9, 3.9, 33)                               # heaters spread out
        for xn, xw in zip(narrow, wide):
            ax.plot(*s_curve(xn, xw, y_split, y_h0 - 0.25), color=core, lw=0.3)
            ax.plot([xw, xw], [y_h0 - 0.25, y_h1 + 0.25], color=core, lw=0.3)
            ax.plot(*s_curve(xw, xn, y_h1 + 0.25, y_g0), color=core, lw=0.3)
        heaters, width = wide, 0.1
        for i, xh in enumerate(heaters):                                    # one wire per heater
            left = i < 17
            k = i if left else 32 - i
            pad_x, pad_y = (0.25 if left else 9.75), 2.0 + k * 0.27
            ax.plot([xh, pad_x], [y_h0, pad_y], color=wire, lw=0.35, zorder=2)
            ax.add_patch(Rectangle((pad_x - 0.13, pad_y - 0.08), 0.26, 0.16, color=ps.INK2, lw=0))
    else:
        for xn in narrow:
            ax.plot([xn, xn], [y_split, y_g0], color=core, lw=0.3)
        ax.add_patch(Rectangle((5 - 40 * scale, y_h0 - 0.35), 80 * scale, y_h1 - y_h0 + 0.7,
                               fc=ps.MAT["oxide"], ec=ps.AXIS, lw=0.6, ls=(0, (3, 2)), zorder=0))
        ax.add_patch(Rectangle((5 + 36 * scale, y_h0 - 0.35), 2 * scale, y_h1 - y_h0 + 0.7,
                               fc="white", ec=ps.INK2, lw=0.5, zorder=4))     # the slot
        heaters, width = 5 + np.array(DESIGNS["slot"][1]) * scale, 0.16
        for k, xh in enumerate(heaters):                                    # three wires, out the near end
            y_w = 3.0 - k * 0.3
            ax.plot([xh, xh, 9.75], [y_h0, y_w, y_w], color=wire, lw=0.6, zorder=2)
            ax.add_patch(Rectangle((9.62, y_w - 0.08), 0.26, 0.16, color=ps.INK2, lw=0))
    for xh in heaters:
        ax.add_patch(Rectangle((xh - width / 2, y_h0), width, y_h1 - y_h0, color=ps.MAT["heater"], lw=0,
                               zorder=3))
    ax.add_patch(Rectangle((narrow[0] - 0.12, y_g0), narrow[-1] - narrow[0] + 0.24, y_g1 - y_g0,
                           fc="white", ec=ps.AXIS, lw=0.6, hatch="------", zorder=1))   # gratings
    for xn in narrow:
        ax.plot([xn, xn], [y_g0, y_g1], color=core, lw=0.3, zorder=2)
    ax.text(5, y_g1 + 0.25, "gratings, 2 µm pitch: light leaves upward", ha="center", fontsize=6.8,
            color=ps.INK2)
    ax.text(5.35, 0.45, "1 → 33 splitter", fontsize=6.8, color=ps.INK2)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    return y_h0, y_h1, scale


def fig12():
    """Slide schematic: a conventional thermo-optic OPA against ours, seen from above."""
    fig, (l, r) = plt.subplots(1, 2, figsize=(ps.FULL, 3.5))
    y_h0, y_h1, scale = draw_top_view(l, "conventional")
    l.text(5, 7.3, "33 heaters, spread out\nso they don't heat each other", ha="center", va="center",
           fontsize=6.8, color=ps.INK, bbox=dict(fc="white", ec="none", pad=1.5), zorder=6)
    l.text(0.1, 1.0, "33 wires\nand pads", ha="left", va="center", fontsize=6.8, color=ps.INK)
    l.set_title("Conventional: a heater and a wire per emitter", fontsize=8.5, color=ps.INK, loc="left")
    draw_top_view(r, "ours")
    r.annotate("suspended\nmembrane,\nno substrate", (5 - 40 * scale + 0.05, 5.8), xytext=(0.2, 6.4),
               fontsize=6.8, color=ps.INK, arrowprops=dict(arrowstyle="-", color=ps.INK2, lw=0.5))
    r.annotate("slot", (5 + 37 * scale, y_h1 + 0.2), xytext=(8.4, 7.0), fontsize=6.8, color=ps.INK,
               arrowprops=dict(arrowstyle="-", color=ps.INK2, lw=0.5))
    r.annotate("3 heaters,\n3 wires", (5 + 14 * scale, 4.3), xytext=(0.2, 4.0), fontsize=6.8,
               color=ps.INK, arrowprops=dict(arrowstyle="-", color=ps.INK2, lw=0.5))
    r.set_title("Ours: a few heaters share one membrane", fontsize=8.5, color=ps.INK, loc="left")
    fig.text(0.5, 0.0, "Top view, not to scale: the array is 64 µm wide, the heated section 1 mm long, "
             "and the heaters sit before the gratings.", ha="center", fontsize=6.6, color=ps.INK2)
    fig.subplots_adjust(wspace=0.06, bottom=0.05, top=0.93)
    ps.save(fig, "fig12_topview")


# ====================================================================== fig13 (slide schematic)
def strehl_curve(G, sign=+1, n=600):
    """Strehl toward the target angle against angle, for one unit fit (as steer_one_direction)."""
    r = steer_one_direction(G, xe, sign)
    step = np.linspace(0, np.pi, n)
    theta = np.degrees(np.arcsin(step * WAVELENGTH_UM / (2 * np.pi * PITCH_UM)))
    return theta, strehl_vs_step(r["residual_per_rad"], step), r


def fig13():
    """Slide schematic: the analysis pipeline. Steps 2-5 run once per heater layout, in a loop
    that keeps the layout with the longest steering range, then the least power."""
    kernel, heaters = DESIGNS["slot"]
    G = phase_matrix(kernel, xe, heaters)
    fit = steer_one_direction(G, xe, +1)
    fig = plt.figure(figsize=(8.0, 3.05))
    steps = [("1  COMSOL, 2D", "one heater at a time,\nat every position:\nits temperature profile"),
             ("2  Pick a layout", "K heater positions\n(here 2, 14 and 32 µm)"),
             ("3  Phase matrix G", "profiles add, so\nphase = G · powers\n(33 emitters × K heaters)"),
             ("4  Fit the powers", "closest to the ramp:\nminimise |G·P − ramp − c|²\nP ≥ 0: heaters only heat\nc free: a shared phase"),
             ("5  Score the layout", "beam quality (Strehl)\nagainst angle; range =\nwhere it falls below 0.8")]
    n = len(steps)
    gap, x0, y0, h = 0.024, 0.006, 0.25, 0.73
    w = (1 - 2 * x0 - (n - 1) * gap) / n
    insets, centres = [], []
    for i, (title, text) in enumerate(steps):
        x = x0 + i * (w + gap)
        centres.append(x + w / 2)
        loop = 1 <= i <= 4
        fig.patches.append(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0,rounding_size=0.015",
                                          transform=fig.transFigure, fc="#eef4fc" if loop else "#f7f6f2",
                                          ec=ps.AXIS, lw=0.7, zorder=-1))     # behind the inset axes
        fig.text(x + 0.01, y0 + h - 0.045, title, fontsize=8, fontweight="bold", color=ps.INK, va="top")
        fig.text(x + 0.01, y0 + 0.025, text, fontsize=6.6, color=ps.INK2, va="bottom", linespacing=1.25)
        insets.append(fig.add_axes([x + 0.02, y0 + 0.29, w - 0.04, 0.31]))
        if i < n - 1:
            fig.patches.append(FancyArrowPatch((x + w + 0.003, y0 + h * 0.68), (x + w + gap - 0.003, y0 + h * 0.68),
                                               transform=fig.transFigure, arrowstyle="-|>",
                                               mutation_scale=9, color=ps.INK2, lw=1.0))
    # The loop: from step 5 back to step 2, underneath the boxes
    y_loop = 0.1
    fig.lines.append(plt.Line2D([centres[4], centres[4], centres[1]], [y0 - 0.005, y_loop, y_loop],
                                transform=fig.transFigure, color=ps.COLOR["slot"], lw=1.1))
    fig.patches.append(FancyArrowPatch((centres[1], y_loop), (centres[1], y0 - 0.005), transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=9, color=ps.COLOR["slot"], lw=1.1))
    fig.text((centres[1] + centres[4]) / 2, y_loop - 0.025,
             "repeat for every layout on a 2 µm grid; keep the longest range, then the least power",
             ha="center", va="top", fontsize=7, color=ps.INK)

    k_ax, l_ax, g_ax, f_ax, s_ax = insets
    x = np.linspace(-40, 40, 401)
    for xh, shade in zip((-24, -8, 8, 24, 34), ps.BLUE_RAMP[:5]):
        k_ax.plot(x, kernel(x, xh) * 1e-3, color=shade, lw=1.0)     # K per mW for a 1 mm heater
    k_ax.set_xlabel("x (µm)", fontsize=6, labelpad=1)

    l_ax.add_patch(Rectangle((-40, -0.25), 76, 0.5, color=ps.MAT["oxide"], lw=0))            # membrane
    l_ax.add_patch(Rectangle((36, -0.25), 2, 0.5, fc="white", ec=ps.INK2, lw=0.5))           # slot
    l_ax.plot(xe, np.zeros_like(xe), "|", color=ps.MAT["core"], ms=5, mew=0.6)
    l_ax.plot(heaters, [0.55] * 3, "v", color=ps.MAT["heater"], ms=6, clip_on=False)
    l_ax.set_xlim(-42, 42)
    l_ax.set_ylim(-0.6, 0.8)
    l_ax.set_xlabel("heaters ▼ above 33 emitters", fontsize=6, labelpad=1)
    l_ax.spines["left"].set_visible(False)

    g_ax.imshow(G, aspect="auto", cmap="Blues", interpolation="nearest")
    g_ax.set_xlabel("3 heaters", fontsize=6, labelpad=1)
    g_ax.set_ylabel("33 emitters", fontsize=6, labelpad=1)

    P = fit["powers_per_rad"]
    for j, shade in enumerate(ps.BLUE_RAMP[:3]):
        f_ax.plot(xe, G[:, j] * P[j], color=shade, lw=0.9)
    f_ax.plot(xe, G @ P, color=ps.COLOR["slot"], lw=2.6)
    f_ax.plot(xe, xe / PITCH_UM + fit["offset_per_rad"], color=ps.INK, lw=0.8, ls=(0, (2, 1.5)))
    f_ax.set_xlabel("sum (thick) on the ramp (dashed)", fontsize=6, labelpad=1)

    for hs, shade in (([34], ps.BLUE_RAMP[1]), ([5, 32], ps.BLUE_RAMP[3]), (heaters, ps.COLOR["slot"])):
        th, S, _ = strehl_curve(phase_matrix(kernel, xe, hs))
        s_ax.plot(th, S, color=shade, lw=1.4 if hs is heaters else 0.9)
    s_ax.axhline(STREHL_MIN, color=ps.INK2, lw=0.6, ls=(0, (2, 2)))
    s_ax.set_ylim(0, 1.05)
    s_ax.set_xlabel("angle 0–22.8°; dashed: 0.8", fontsize=6, labelpad=1)
    for axx in insets:
        axx.set_xticks([])
        axx.set_yticks([])
        axx.grid(False)
        axx.set_facecolor("white")
        for side in ("left", "bottom"):
            axx.spines[side].set_color(ps.AXIS)
    ps.save(fig, "fig13_pipeline")


# ====================================================================== fig14 (slide schematic)
def fig14():
    """Slide schematic: the design path in four steps, with what each change bought and cost."""
    steps = {name: (t, p, t_fast) for name, _, t, p, t_fast in design_steps()}
    columns = [
        ("standard", "1  Standard", "a heater on every emitter", "the conventional design", 33),
        ("undercut", "2  Undercut", "substrate etched away", "heat sank straight down", 3),
        ("slot", "3  Slot", "slot past the edge heater", "edge heat escaped to the near wall", 3),
        ("vacuum", "4  Vacuum", "gas pumped out", "air leaked across slot and cavity", 1),
    ]
    fig = plt.figure(figsize=(ps.FULL, 2.45))
    w, gap, x0 = 0.215, 0.032, 0.012
    for i, (name, title, change, reason, wires) in enumerate(columns):
        x = x0 + i * (w + gap)
        fig.patches.append(Rectangle((x, 0.905), 0.03, 0.045, transform=fig.transFigure,
                                     color=ps.COLOR[name], lw=0))
        fig.text(x + 0.04, 0.927, title, fontsize=8.5, fontweight="bold", color=ps.INK, va="center")
        ax = fig.add_axes([x, 0.47, w, 0.4])
        draw_section(ax, name, labels=False)
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
        t, p, t_fast = steps[name]
        settle = f"settles in {t:.2g} ms" + (" (est.)" if name == "standard" else
                                             f", {t_fast:.2g} shaped")
        fig.text(x + 0.005, 0.40, change, fontsize=7.2, color=ps.INK, va="top")
        fig.text(x + 0.005, 0.31, ("why: " if i else "") + reason, fontsize=6.4, color=ps.INK2, va="top",
                 style="italic")
        power = f"{p:.0f}" if p >= 100 else f"{p:.2g}"
        fig.text(x + 0.005, 0.18, f"{power} mW at 10° · {wires} wire{'s' if wires > 1 else ''}",
                 fontsize=7.4, color=ps.INK, va="top", fontweight="bold")
        fig.text(x + 0.005, 0.08, settle, fontsize=6.6, color=ps.INK2, va="top")
        if i < 3:
            fig.patches.append(FancyArrowPatch((x + w + 0.004, 0.67), (x + w + gap - 0.004, 0.67),
                                               transform=fig.transFigure, arrowstyle="-|>",
                                               mutation_scale=9, color=ps.INK2, lw=1.0))
    ps.save(fig, "fig14_design_path")


# ====================================================================== fig15 (slide schematic)
def draw_emitters(ax, x_em, phase, color, label):
    """Emitters along y = 0 and, below them, a bar per emitter for the phase it is given."""
    for x0 in x_em:
        ax.add_patch(Rectangle((x0 - 0.18, -0.12), 0.36, 0.24, color=ps.MAT["core"], lw=0, zorder=4))
    base = -2.9
    for x0, p in zip(x_em, phase):
        ax.add_patch(Rectangle((x0 - 0.22, base), 0.44, 0.25 + 1.5 * p, color=color, lw=0, alpha=0.85))
    ax.plot([x_em[0] - 0.5, x_em[-1] + 0.5], [base, base], color=ps.AXIS, lw=0.6)
    ax.text(0, base - 0.25, label, ha="center", va="top", fontsize=6.8, color=ps.INK2)


def draw_wavelets(ax, x_em, radii, aims_deg, spread=38):
    """Huygens wavelets: an arc of radius r_n around each emitter, centred on its aim direction."""
    for x0, r, aim in zip(x_em, radii, aims_deg):
        t = np.radians(np.linspace(aim - spread, aim + spread, 50))
        ax.plot(x0 + r * np.sin(t), r * np.cos(t), color=ps.MUTED, lw=0.6, alpha=0.6)


def fig15():
    """Slide schematic: how a phased array steers (a straight phase ramp), how it focuses (a parabola),
    and how heat sets the phase. Geometry is schematic (9 emitters); the numbers are the real ones."""
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(ps.FULL * 1.25, 3.0),
                                  gridspec_kw={"width_ratios": [1, 1, 1.05]})
    x_em = np.arange(-4, 5) * 1.0
    blue, yellow = ps.COLOR["slot"], "#c98a00"          # yellow darkened a little for lines on white

    # (a) steer: a ramp of delay tilts the wavefront; the beam leans toward the delayed (hot) side
    th = np.radians(20)
    R0 = 4.6
    r = R0 - np.sin(th) * x_em
    draw_wavelets(a, x_em, r, np.full_like(x_em, 20.0), spread=30)
    tx, ty = x_em + r * np.sin(th), r * np.cos(th)                 # tangent points: all on one line
    ext = 0.8
    a.plot([tx[0] - ext * np.cos(th), tx[-1] + ext * np.cos(th)],
           [ty[0] + ext * np.sin(th), ty[-1] - ext * np.sin(th)], color=blue, lw=2.2)
    mx, my = R0 * np.sin(th), R0 * np.cos(th)
    a.annotate("", xy=(mx + 2.3 * np.sin(th), my + 2.3 * np.cos(th)), xytext=(mx, my),
               arrowprops=dict(arrowstyle="-|>", color=blue, lw=1.6, mutation_scale=10))
    a.plot([mx, mx], [my, my + 2.4], color=ps.INK2, lw=0.6, ls=(0, (2, 2)))
    t = np.radians(np.linspace(0, 20, 20))
    a.plot(mx + 1.7 * np.sin(t), my + 1.7 * np.cos(t), color=ps.INK2, lw=0.6)
    a.text(mx + 0.35, my + 2.05, "θ", fontsize=8, color=ps.INK)
    a.text(-4.7, 7.0, "wavefront", fontsize=6.8, color=blue, va="top")
    draw_emitters(a, x_em, (x_em + 4) / 8, blue, "phase per emitter: a straight ramp,\nthe same step Δφ each time")
    a.text(-4.8, 9.0, "sin θ = λ·Δφ / (2π·d)", fontsize=8, color=ps.INK, fontweight="bold")
    a.set_title("Steer: a phase ramp tilts the beam", loc="left", fontsize=8.5, color=ps.INK, pad=14)

    # (b) focus: the centre delayed most curves the wavefront inward, to a focus at f
    f, rho = 7.6, 5.6
    D = np.hypot(x_em, f)
    r = D - rho
    aims = np.degrees(np.arctan2(-x_em, f))
    draw_wavelets(b, x_em, r, aims, spread=30)
    for x0 in x_em:
        b.plot([x0, 0], [0, f], color=yellow, lw=0.5, alpha=0.6)
    ang = np.arctan2(-f, x_em)                                          # from the focus to each emitter
    arc = np.linspace(ang[-1] - 0.06, ang[0] + 0.06, 60)
    b.plot(rho * np.cos(arc), f + rho * np.sin(arc), color=yellow, lw=2.2)
    b.plot(0, f, "o", color=yellow, ms=5)
    b.text(0.4, f, "focus", fontsize=6.8, color=ps.INK, va="center")
    b.annotate("", xy=(4.8, f), xytext=(4.8, 0), arrowprops=dict(arrowstyle="<->", color=ps.INK2, lw=0.6))
    b.text(5.0, f / 2, "f", fontsize=8, color=ps.INK, va="center")
    draw_emitters(b, x_em, 1 - (x_em / 4) ** 2, yellow, "phase per emitter: a parabola,\ncentre ahead of the edges by a")
    b.text(-4.8, 9.0, "f = k·X² / (2a)", fontsize=8, color=ps.INK, fontweight="bold")
    b.set_title("Focus: a phase parabola is a lens", loc="left", fontsize=8.5, color=ps.INK, pad=14)

    for axx in (a, b):
        axx.set_xlim(-5.2, 6.4)
        axx.set_ylim(-4.2, 9.6)
        axx.set_aspect("equal")
        axx.axis("off")

    # (c) heat sets the phase: phase = 0.81 rad per K for a 1 mm heater, so the temperature
    # profile must have the shape of the phase profile
    x = np.linspace(-32, 32, 200)
    rad_per_K = PHASE_COEFF * 1e-3
    ramp = dphi(10) * (x + 32) / PITCH_UM                 # the 10° ramp, rad
    c.plot(x, ramp / rad_per_K, color=blue, lw=2.0)
    dome = 20 * (1 - (x / 32) ** 2)                       # a 20 rad lens (f ≈ 0.1 mm), rad
    c.plot(x, dome / rad_per_K, color=yellow, lw=2.0)
    c.text(31, 30, f"steer to 10°:\na {ramp[-1]:.0f} rad = {ramp[-1] / rad_per_K:.0f} K ramp",
           ha="right", va="top", fontsize=6.6, color=ps.INK)
    c.text(-18, 25, f"focus, f ≈ {LENS_AF / 20 / 1000:.1f} mm:\na {20 / rad_per_K:.0f} K dome", ha="center",
           va="bottom", fontsize=6.6, color=ps.INK)
    c.set_xlim(-34, 34)
    c.set_ylim(0, 66)
    c.set_xticks([-32, 0, 32])
    c.set_xlabel("x across the array (µm)")
    c.set_ylabel("temperature rise (K)")
    c.text(-32, 63, "Δφ = (2π/λ)·(dn/dT)·L·ΔT\n= 0.81 rad per K for L = 1 mm", fontsize=6.8, color=ps.INK,
           va="top", fontweight="bold")
    c.set_title("Heat sets the phase: shape the heat", loc="left", fontsize=8.5, color=ps.INK, pad=14)
    fig.text(0.01, 0.0, f"Silicon: dn/dT ≈ 2×10⁻⁴ per K.  d = {PITCH_UM:g} µm pitch, λ = {WAVELENGTH_UM} µm: "
             f"only one beam within ±{np.degrees(np.arcsin(WAVELENGTH_UM / (2 * PITCH_UM))):.1f}° (Δφ ≤ π).  "
             "k = 2π/λ;  X = 32 µm, half the array.", fontsize=6.6, color=ps.INK2)
    fig.subplots_adjust(wspace=0.12, left=0.01, right=0.99, bottom=0.17, top=0.88)
    ps.save(fig, "fig15_opa_basics")


# ====================================================================== fig16 (slide / backup)
def fig16():
    """The placement search on the slot membrane: what each extra heater buys, and where it goes."""
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(ps.FULL * 1.25, 2.75),
                                  gridspec_kw={"width_ratios": [1, 1, 1.12]})
    # (a) range = where the Strehl falls below 0.8, for the best layout of 1, 2 and 3 heaters
    best = [([34], "1 heater"), ([5, 32], "2 heaters"), ([2, 14, 32], "3 heaters")]
    shades = [ps.BLUE_RAMP[1], ps.BLUE_RAMP[3], ps.BLUE_RAMP[5]]
    for (hs, label), shade, (tx, ty) in zip(best, shades, [(4.8, 0.35), (12.5, 0.6), (5.5, 1.13)]):
        th, S, r = strehl_curve(phase_matrix(slot, xe, hs))
        a.plot(th, S, color=shade, lw=1.8)
        rng = r["range_deg"]
        a.plot(rng, np.interp(rng, th, S), "o", color=shade, ms=4.5, mec="white", mew=0.6, zorder=4)
        where = ", ".join(f"{v:g}" for v in hs)
        sep = " " if rng > 22.7 else "\n"           # the top label has room for one line
        a.text(tx, ty, f"{label} ({where} µm):{sep}{'full ' if rng > 22.7 else ''}{rng:.1f}°", fontsize=6.4,
               color=ps.INK, va="top")
    a.axhline(STREHL_MIN, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
    a.text(0.4, STREHL_MIN + 0.015, "0.8", fontsize=6.4, color=ps.INK2, va="bottom")
    a.set_xlim(0, 23.5)
    a.set_ylim(0, 1.16)
    a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_xlabel("beam angle (°)")
    a.set_ylabel("beam quality (Strehl)")
    ps.panel_label(a, "a")
    a.set_title("Each heater straightens the ramp", loc="left", fontsize=8.5, color=ps.INK, pad=10)

    # (b) two heaters: the edge one fixed at 32 µm, the second moved across the array
    pos = np.arange(-32, 31)
    rng = np.array([steer_one_direction(phase_matrix(slot, xe, [v, 32]), xe, +1)["range_deg"] for v in pos])
    alone = steer_one_direction(phase_matrix(slot, xe, [32]), xe, +1)["range_deg"]
    b.plot(pos, rng, color=ps.COLOR["slot"], lw=1.8)
    k = int(np.argmax(np.round(rng, 2)))
    b.plot(pos[k], rng[k], "o", color=ps.COLOR["slot"], ms=4.5, mec="white", mew=0.6, zorder=4)
    b.text(pos[k] + 2, rng[k] + 0.4, f"best: {pos[k]} µm, {rng[k]:.1f}°", fontsize=6.4, color=ps.INK)
    b.axhline(alone, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
    b.text(-31, alone + 0.5, f"edge heater alone: {alone:.1f}°", fontsize=6.4, color=ps.INK2)
    b.set_xlim(-33, 31)
    b.set_ylim(0, 23)
    b.set_xticks([-32, -16, 0, 16, 32])
    b.set_xlabel("second heater's position (µm)")
    b.set_ylabel("steering range (°)")
    ps.panel_label(b, "b")
    b.set_title("Two heaters: slide the second one", loc="left", fontsize=8.5, color=ps.INK, pad=10)

    # (c) three heaters: every pair of interior positions, edge heater at 32 µm (1 µm grid)
    grid = np.arange(-32, 31)
    R = np.full((grid.size, grid.size), np.nan)
    M = np.full_like(R, np.nan)
    for i, p in enumerate(grid):
        for j, q in enumerate(grid):
            if q - p >= 2 and 32 - q >= 2:
                r = steer_one_direction(phase_matrix(slot, xe, [p, q, 32]), xe, +1)
                R[j, i], M[j, i] = r["range_deg"], r["mw_per_deg"]
    im = c.imshow(R, origin="lower", cmap="Blues", vmin=0, vmax=22.8, interpolation="nearest",
                  extent=[grid[0] - 0.5, grid[-1] + 0.5, grid[0] - 0.5, grid[-1] + 0.5])
    c.contour(grid, grid, np.nan_to_num(R), levels=[22.75], colors=[ps.INK2], linewidths=0.7)
    full = R >= 22.75
    jj, ii = np.unravel_index(np.nanargmin(np.where(full, M, np.inf)), R.shape)
    c.plot(2, 14, "*", color=ps.HEAT, ms=9, mec="white", mew=0.6, zorder=5)
    c.text(4, 12.5, "chosen: 2, 14", fontsize=6.4, color=ps.INK, va="top")
    c.text(3, -9, f"{full.sum()} of {np.isfinite(R).sum()} layouts\nreach the full 22.8°\n"
           "(inside the line);\nthe cheapest wins", fontsize=6.2, color=ps.INK2, va="top")
    print(f"   fig16: cheapest full-range layout on the 1 µm grid: {grid[ii]}, {grid[jj]}, 32 at "
          f"{M[jj, ii]:.4f} mW/deg; chosen 2, 14, 32 at {M[14 + 32, 2 + 32]:.4f}")
    cb = fig.colorbar(im, ax=c, fraction=0.05, pad=0.03)
    cb.set_label("steering range (°)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    cb.outline.set_visible(False)
    c.set_xticks([-32, -16, 0, 16, 30])
    c.set_yticks([-32, -16, 0, 16, 30])
    c.grid(False)
    c.set_xlabel("first heater (µm)")
    c.set_ylabel("second heater (µm)")
    ps.panel_label(c, "c")
    c.set_title("Three heaters: try every pair", loc="left", fontsize=8.5, color=ps.INK, pad=10)
    fig.tight_layout(w_pad=1.6)
    ps.save(fig, "fig16_placement")


# ====================================================================== fig17 (slide version of fig11 b)
def fault_value(faults, section, design, case, quantity):
    row = faults[(faults.section == section) & (faults.design == design) & (faults.case == case)
                 & (faults.quantity == quantity)]
    return float(row.value.iloc[0])


def fig17():
    """Slide chart: how far a 10° beam drifts after calibration, drive held fixed (linear scale)."""
    faults = pd.read_csv(RESULTS + "faults.csv")
    get = lambda *key: abs(fault_value(faults, *key))
    designs4 = ["vacuum, 1 heater", "slot, 3 heaters", "undercut, 3 heaters", "dual-use, 5 heaters"]
    noise = [get("drive noise", d, "independent", "pointing_rms_deg") for d in designs4]
    amb = [get("drift", "all designs", f"ambient +10 K, b = {s}", "pointing_error_deg")
           for s in ("0.0e+00", "7.0e-04", "1.4e-03")]
    vac = lambda p: get("drift", "vacuum, 1 heater", f"pressure {p} Pa", "pointing_error_deg")
    weather = max(get("drift", "slot in air, 3 heaters", f"pressure {p} atm", "pointing_error_deg")
                  for p in ("0.95", "1.05"))
    ours = ps.INK2
    rows = [  # (label, (low, high) or value, colour, value text); None = a group heading
        ("Driver noise, 1% rms", None),
        ("ours, any design (angle = amount of heat)", (min(noise), max(noise)), ours, f"{min(noise):.2f}–{max(noise):.2f}°"),
        ("conventional, 33 heaters (angle = wrap pattern)",
         get("drive noise", "conventional", "independent", "pointing_rms_deg"), ps.MUTED, None),
        ("Room 10 K warmer", None),
        ("ours: membrane and silicon (dn/dT(T) unverified)", (min(amb), max(amb)), ours,
         f"{min(amb):.2f}–{max(amb):.2f}°"),
        ("slot in air: the air conducts more", get("drift", "slot in air, 3 heaters", "ambient +10 K, gas only",
                                                   "pointing_error_deg"), ps.COLOR["slot"], None),
        ("Gas pressure changes", None),
        ("slot in air: weather, ±5%", weather, ps.COLOR["slot"], None),
        ("vacuum package: 1 Pa → 10 Pa", vac(10), ps.COLOR["vacuum"], None),
        ("vacuum package: 1 Pa → 30 Pa", vac(30), ps.COLOR["vacuum"], None),
        ("vacuum package: 1 Pa → 100 Pa", vac(100), ps.COLOR["vacuum"], None),
    ]
    fig, ax = plt.subplots(figsize=(ps.FULL * 1.1, 3.4))
    y, ticks, labels, heads = 0.0, [], [], []
    for row in rows:
        ticks.append(y)
        labels.append(row[0])
        heads.append(row[1] is None)
        if row[1] is not None:
            _, v, color, text = row
            if isinstance(v, tuple):          # a range: solid to the low end, lighter to the high end
                ax.barh(y, v[0], height=0.58, color=color, lw=0)
                ax.barh(y, v[1] - v[0], left=v[0], height=0.58, color=color, alpha=0.4, lw=0)
                hi = v[1]
            else:
                ax.barh(y, v, height=0.58, color=color, lw=0)
                hi = v
            # A label that would run into the tolerance line starts just past it instead
            ax.text(0.126 if 0.07 < hi < 0.12 else hi + 0.006, y, text or (f"{v:.3f}°" if v < 0.01 else f"{v:.2f}°"), va="center",
                    fontsize=7, color=ps.INK)
        y -= 1
    tol = 0.12
    ax.axvline(tol, color=ps.INK, lw=0.8, ls=(0, (3, 3)))
    ax.text(tol + 0.005, 0.55, "tolerance: a tenth of the beam width (0.12°)", fontsize=7, color=ps.INK)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7.2)
    for t, head in zip(ax.get_yticklabels(), heads):
        if head:
            t.set_fontweight("bold")
            t.set_color(ps.INK)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(y + 0.5, 1.0)
    ax.set_xlabel("how far the beam moves from 10° (°), drive held at its calibrated values")
    ax.set_title("After calibration, the beam is only as steady as the heat-to-phase gain", loc="left",
                 fontsize=8.5, color=ps.INK, pad=14, x=-0.62)
    fig.text(0.01, -0.02, "Gas and pressure rows: 1D membrane model calibrated to COMSOL (est.). "
             "The fix for all of them is one gain number per condition, or a closed loop.",
             fontsize=6.4, color=ps.INK2)
    ps.save(fig, "fig17_drift")


# ====================================================================== fig18 (slide table)
SCOREBOARD = [  # (topic, ours, best published, verdict, verdict detail)
    ("Few-wire steering", "1–3 wires,\n0.5–1.4 mW per degree",
     "12–17 mW per degree*\n[Yaacobi 2014; Poulton 2016]", "win", "10–30× less power"),
    # The † values are not built arrays: one published shifter's P_π, under all 33 emitters of an ideal
    # array (no crosstalk, phases wrapped), as fig04's grey curves. Each line of ours pairs with the same line.
    ("Power at 10°", "5.3 mW, vacuum (1 wire)\n14 mW, slot (3 wires)",
     "1.5 mW: 33 × folded record† [Lu 2015]\n14 mW: 33 × suspended† [Fang 2011]", "mixed",
     "lose 3.6× to record\ntie with suspended"),
    ("Speed, 0 → 10°", "2–9 ms plain step\n0.16–1.4 ms shaped drive",
     "2.4 µs – 0.14 ms\nsingle heaters", "lose", "power × time 1.2–16×\nworse, even shaped"),
    ("Tunable lens", "≈ 1 mW per rad\n1–5 wires",
     "40–110 mW per rad*\n[Liu 2023; Hosseini 2025]", "win", "30–100× less power"),
    ("Calibration", "trim once; then one fit\n+ one gain per angle",
     "a search per beam direction,\n40 s – 1 h [Lian 2024]", "win", "given the trim"),
    ("Per-waveguide errors", "not correctable:\nneeds a one-time trim",
     "fully corrected by\nthe 33 heaters", "lose", "trim required"),
]
VERDICT = {"win": ("▲ Win", "#0ca30c", "white"), "lose": ("▼ Lose", "#d03b3b", "white"),
           "mixed": ("◆ Mixed", "#fab219", ps.INK)}


def fig18():
    """Slide table: where the design wins and loses against the best published devices."""
    fig = plt.figure(figsize=(8.0, 3.7))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    cols = [0.015, 0.215, 0.47, 0.75]
    top, row_h = 0.9, 0.128
    for x, head in zip(cols, ["", "Ours (simulated)", "Best published (measured)", "Verdict"]):
        ax.text(x, top + 0.035, head, fontsize=8.5, fontweight="bold", color=ps.INK, va="center")
    ax.plot([0.01, 0.99], [top, top], color=ps.INK2, lw=0.8)
    for i, (topic, ours, theirs, verdict, detail) in enumerate(SCOREBOARD):
        yc = top - (i + 0.5) * row_h
        if i % 2 == 0:
            ax.add_patch(Rectangle((0.01, yc - row_h / 2), 0.98, row_h, color="#f7f6f2", lw=0, zorder=0))
        ax.text(cols[0], yc, topic, fontsize=8, fontweight="bold", color=ps.INK, va="center")
        ax.text(cols[1], yc, ours, fontsize=7.4, color=ps.INK, va="center", linespacing=1.3)
        ax.text(cols[2], yc, theirs, fontsize=7.4, color=ps.INK2, va="center", linespacing=1.3)
        word, fc, tc = VERDICT[verdict]
        ax.add_patch(FancyBboxPatch((cols[3], yc - 0.033), 0.085, 0.066, boxstyle="round,pad=0,rounding_size=0.012",
                                    fc=fc, ec="none", zorder=2))
        ax.text(cols[3] + 0.0425, yc, word, fontsize=7.4, fontweight="bold", color=tc, ha="center", va="center",
                zorder=3)
        ax.text(cols[3] + 0.1, yc, detail, fontsize=7.4, color=ps.INK, va="center")
    bottom = top - len(SCOREBOARD) * row_h
    ax.plot([0.01, 0.99], [bottom, bottom], color=ps.AXIS, lw=0.6)
    ax.text(0.015, bottom - 0.04, "* the paper's array scaled to 33 emitters (est.).   † our estimate, never "
            "built: one published shifter's P_π, under each of 33 emitters (no crosstalk, phases wrapped).",
            fontsize=6.6, color=ps.INK2, va="center")
    ax.text(0.015, bottom - 0.085, "Ours: 2D simulation.  Theirs: measured devices.  "
            "Shaped drive = pre-emphasis at up to 3× the hold power.", fontsize=6.6, color=ps.INK2, va="center")
    ps.save(fig, "fig18_scoreboard")


# ====================================================================== fig19 (slide version of fig10)
def lens_drive_phase(name, theta, f_um):
    """What a lens layout actually makes when fitted to ramp(θ) + lens(f): (phase, fitted target, powers).
    f_um = None means no lens. The target includes the free offset, so it lines up with the phase."""
    a = 0.0 if f_um is None else LENS_AF / f_um
    target = dphi(theta) * xe / PITCH_UM + a * unit_target(xe)
    G = lens_matrix(name)
    P, c = fit_powers(G, target)
    return G @ P, target + c, P


def fig19(report=False):
    """Slide figure for the lens, from a 2D light simulation above the array (opa/farfield.py):
    (a, b) one small heater against one wide heater asked for the same lens; (c) the yardstick,
    the share of a perfect lens's light that reaches the target focus; (d-f) steering, focusing
    and both at once with the same 5 wires."""
    f_um, theta = 150.0, 10.0
    ep = element_pattern()
    x = np.linspace(-60, 120, 361)
    z = np.linspace(40, 260, 441)
    X, Z = np.meshgrid(x, z)
    on_target = lambda ph, th, f: near_field_intensity(
        ph, xe, np.array([[f * np.cos(np.radians(th)) ** 2 * np.sin(np.radians(th))]]),
        np.array([[f * np.cos(np.radians(th)) ** 3]]), ep)[0, 0]
    perfect = lambda th, f: on_target(dphi(th) * xe / PITCH_UM + LENS_AF / f * unit_target(xe), th, f)
    ref_map = near_field_intensity(LENS_AF / f_um * unit_target(xe), xe, X, Z, ep)
    ref_peak = ref_map[np.abs(X) < 35].max()          # the perfect lens's brightest point: the yardstick

    fig = plt.figure(figsize=(8.4, 6.7))
    grid = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.34)
    norm = matplotlib.colors.PowerNorm(0.5, vmin=0, vmax=1)    # square-root scale: faint light stays visible

    def light_map(cell, name, th, f, title, note, own_scale, first_col):
        phase, target, _ = lens_drive_phase(name, th, f)
        sub = cell.subgridspec(2, 1, height_ratios=[4, 1.15], hspace=0.08)
        m = fig.add_subplot(sub[0])
        s = fig.add_subplot(sub[1], sharex=m)
        I = near_field_intensity(phase, xe, X, Z, ep)
        beam = np.abs(X - Z * np.tan(np.radians(th))) < 35
        m.imshow(I / (I[beam].max() if own_scale else ref_peak), origin="lower", cmap="Blues", norm=norm,
                 extent=[x[0], x[-1], z[0], z[-1]], interpolation="bilinear")
        if f is not None:                  # where a perfect lens would put the focus
            R = f * np.cos(np.radians(th)) ** 2
            m.plot(R * np.sin(np.radians(th)), R * np.cos(np.radians(th)), "+", color=ps.HEAT, ms=9, mew=1.4)
        else:                              # the intended beam direction
            m.plot([0, z[-1] * np.tan(np.radians(th))], [0, z[-1]], color=ps.HEAT, lw=0.8, ls=(0, (3, 3)))
        m.set_xlim(x[0], x[-1])
        m.set_ylim(z[0], z[-1])
        m.grid(False)
        m.tick_params(labelbottom=False)
        m.set_title(title, loc="left", fontsize=8, color=ps.INK, pad=4)
        m.text(-56, z[-1] - 6, note, fontsize=6.3, color=ps.INK, va="top",
               bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
        if first_col:
            m.set_ylabel("height above chip (µm)", fontsize=7)
        # Below the map: the phase the heaters make across the array, against the target
        color = ps.COLOR["dual"] if name == "dual" else ps.COLOR["undercut"]
        s.plot(xe, phase, color=color, lw=1.6)
        s.plot(xe, target, color=ps.INK, lw=0.8, ls=(0, (2, 1.5)))
        lo, hi = min(phase.min(), target.min()), max(phase.max(), target.max())
        base = lo - 0.32 * (hi - lo)
        x_h = LENS_LAYOUTS[name][0]
        if name == "uniform":
            s.add_patch(Rectangle((-32, base), 64, 0.13 * (hi - lo), color=ps.MAT["heater"], lw=0))
        else:
            s.plot(x_h, [base + 0.06 * (hi - lo)] * len(x_h), "^", color=ps.MAT["heater"], ms=4)
        s.set_ylim(base, hi + 0.1 * (hi - lo))
        s.set_yticks([])
        s.grid(False)
        s.spines["left"].set_visible(False)
        s.set_xticks([-32, 0, 32])
        s.tick_params(axis="x", labelsize=6.5)
        s.set_xlabel("across the array (µm)", fontsize=6.5, labelpad=1)
        s.text(40, (lo + hi) / 2, "phase made\n(dashed: ideal)", fontsize=5.8, color=ps.INK2, va="center")
        return phase

    # Row 1: one small heater against one wide heater, both asked for f = 150 µm
    s_centre = on_target(lens_drive_phase("centre", 0, f_um)[0], 0, f_um) / perfect(0, f_um)
    s_wide = on_target(lens_drive_phase("uniform", 0, f_um)[0], 0, f_um) / perfect(0, f_um)
    light_map(grid[0, 0], "centre", 0, f_um, "(a) One small heater",
              f"a tent of phase: a long streak,\nno single focus. Light at the +:\n"
              f"{s_centre:.2f} of a perfect lens", False, True)
    light_map(grid[0, 1], "uniform", 0, f_um, "(b) One 64 µm wide heater",
              f"a parabola: a sharp focus\nat the + (f = {f_um / 1000:.2f} mm).\nLight at the +: {s_wide:.2f}",
              False, False)

    c = fig.add_subplot(grid[0, 2])
    if report:
        draw_lens_power(c, compact=True)
        c.set_title("(c) Power to hold the lens", loc="left", fontsize=8, color=ps.INK, pad=4)
        c.xaxis.label.set_size(7)
        c.yaxis.label.set_size(7)
    # (c, slide version) the yardstick at every focal length: lines from the phase error, dots from
    # the light simulation
    f_line = np.geomspace(2000, 40, 300)
    f_dots = np.array([1000, 500, 300, 150, 100])
    yard = c if not report else None
    for name, label, (tx, ty) in (("uniform", "one 64 µm heater", (1.6, 1.02)),
                                  ("dual", "dual-use, 5 wires", (0.3, 0.82)),
                                  ("centre", "one small heater", (0.9, 0.62))):
        _, color, ls, _ = LENS_LAYOUTS[name]
        G = lens_matrix(name)
        P, off = fit_powers(G, unit_target(xe))
        S = strehl_vs_step(G @ P - unit_target(xe) - off, LENS_AF / f_line)
        dots = [on_target(lens_drive_phase(name, 0, f)[0], 0, f) / perfect(0, f) for f in f_dots]
        s_at = np.interp(f_dots[::-1], f_line[::-1], S[::-1])[::-1]
        print(f"   {name}: light simulation − Strehl at f = {list(f_dots)} µm: "
              + ", ".join(f"{d - t:+.3f}" for d, t in zip(dots, s_at)))
        if yard is None:
            continue
        c.plot(f_line / 1000, S, color=color, ls=ls, lw=1.6)
        c.plot(f_dots / 1000, dots, "o", color=color, ms=3.5, mfc="white", mew=1.0)
        c.text(tx, ty, label, fontsize=6.3, color=ps.INK, va="bottom")
    if yard is not None:                    # slide version: the yardstick panel
        c.axhline(STREHL_MIN, color=ps.INK2, lw=0.7, ls=(0, (3, 3)))
        c.text(1.9, STREHL_MIN - 0.02, "0.8: a good focus", fontsize=6.0, color=ps.INK2, va="top")
        c.axvline(f_um / 1000, color=ps.AXIS, lw=0.8)
        c.text(f_um / 1000 * 0.94, 1.05, "maps", fontsize=6.0, color=ps.INK2, ha="left", va="bottom")
        c.axvspan(0.1, 0.04, color=ps.GRID, alpha=0.6, lw=0)
        c.set_xscale("log")
        c.set_xlim(2, 0.04)
        c.set_xticks([2, 1, 0.5, 0.2, 0.1, 0.05])
        c.minorticks_off()
        c.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        c.set_ylim(0, 1.12)
        c.set_xlabel("focal length asked for (mm); stronger →", fontsize=7)
        c.set_ylabel("light at the target focus ÷ perfect lens", fontsize=7)
        c.set_title("(c) How we judge a lens", loc="left", fontsize=8, color=ps.INK, pad=4)
        c.text(0.03, 0.04, "lines: phase-error formula (Strehl)\ndots: light simulation",
               transform=c.transAxes, ha="left", fontsize=5.8, color=ps.INK2)

    # Row 2: the same 5 wires steer, focus, or both
    ph, tg, _ = lens_drive_phase("dual", theta, None)
    steer_s = float(np.abs(np.exp(1j * (ph - tg)).mean()) ** 2)       # far-field Strehl of the steered beam
    s_both = on_target(lens_drive_phase("dual", theta, f_um)[0], theta, f_um) / perfect(theta, f_um)
    s_focus = on_target(lens_drive_phase("dual", 0, f_um)[0], 0, f_um) / perfect(0, f_um)
    light_map(grid[1, 0], "dual", theta, None, "(d) Steer: a ramp",
              f"beam tilted {theta:.0f}°, not focused\n(Strehl {steer_s:.2f})", True, True)
    light_map(grid[1, 1], "dual", 0, f_um, "(e) Focus: a parabola",
              f"spot at {f_um / 1000:.2f} mm, straight up\nlight at the +: {s_focus:.2f}", True, False)
    light_map(grid[1, 2], "dual", theta, f_um, "(f) Both: ramp + parabola",
              f"spot at {f_um / 1000:.2f} mm, {theta:.0f}° off axis\nlight at the +: {s_both:.2f}", True, False)
    if not report:
        fig.text(0.5, 0.005, "2D light simulation in the plane across the array (square-root colour scale). "
                 "Rows: (a, b) heaters on the undercut membrane; (d–f) the dual-use layout, 5 wires at 0, ±16, "
                 "±33 µm. Faint side streaks are the 2 µm pitch's copies.", ha="center", fontsize=6.0,
                 color=ps.INK2)
    print(f"   fig19: light at the focus / perfect: small {s_centre:.3f}, wide {s_wide:.3f}, "
          f"dual focus {s_focus:.3f}, dual steer+focus {s_both:.3f}; steering Strehl {steer_s:.3f}")
    ps.save(fig, "rep_lens" if report else "fig19_lens_slide")


def rep_lens():
    """Report version of fig19: panel (c) shows the power to hold the lens instead of the yardstick."""
    fig19(report=True)


FIGURES = {f.__name__: f for f in (fig01, fig02, fig03, fig04, fig05, fig06, fig07, fig08, fig09,
                                   fig10, fig11, fig12, fig13, fig14, fig15, fig16, fig17,
                                   fig18, fig19, rep_lens)}

if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGURES):
        print(name)
        FIGURES[name]()
