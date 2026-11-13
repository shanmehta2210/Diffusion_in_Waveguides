"""Validation figures for the report: the COMSOL model checks and the far-field engine's
analytic limits. Writes figures/rep_<name>.pdf and .png.

Run all:   python make_validation_figures.py
Run one:   python make_validation_figures.py rep_mode

Every number drawn or annotated is computed here from the COMSOL exports (the CSVs in the
project folder) or from the same opa/ functions the analyses use. Nothing is typed in from
the notes, so a figure cannot disagree with the data. Each function prints the numbers it
found; they are collected in figures/REPORT_VALIDATION_NOTES.md.

Colour: the four stack colours keep their meaning (grey = standard stack, aqua = undercut).
Checks that belong to no stack use dark ink for the COMSOL result, a grey line for the
analytic prediction, and one purple accent for a second prediction.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Rectangle
from scipy.optimize import brentq

from opa import plotstyle as ps
from opa.comsol_io import read_comsol
from opa.constants import DATA_DIR, DNEFF_DT, K0, N_EMITTERS, PITCH_UM, WAVELENGTH_UM
from opa.farfield import U_GRID, beam_metrics, element_pattern, far_field, strehl
from opa.thermal import StandardKernel, UndercutKernel, emitter_positions

sys.stdout.reconfigure(encoding="utf-8")
ps.apply()

W = 6.3                 # figure width in inches: the report's 160 mm text width
ACCENT = "#7b5ea7"      # purple: a second prediction; not one of the stack colours
BROWN = "#9c6b3f"
SMALL = 7.2             # annotation font size, pt

# Indices used in the COMSOL models (build guide, Models A and E)
N_SI, N_OX = 3.476, 1.444
W_CORE, H_CORE = 0.50, 0.22          # strip cross-section, µm
T0 = 293.15


def shade(color, t):
    """Blend a colour toward white (t < 0) or toward ink (t > 0); t in [-1, 1]."""
    c, target = np.array(to_rgb(color)), np.array(to_rgb("#ffffff" if t < 0 else ps.INK))
    return tuple(c + abs(t) * (target - c))


def note(ax, x, y, text, **kw):
    """Small annotation in axes coordinates."""
    kw.setdefault("fontsize", SMALL)
    kw.setdefault("color", ps.INK2)
    kw.setdefault("transform", ax.transAxes)
    ax.text(x, y, text, **kw)


# ==================================================================== rep_mode
def load_mode(file):
    """Complex Ex and Ey of the TE0 mode on the export grid: arrays [y, x]."""
    _, data = read_comsol(DATA_DIR / file)
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    # The mode is almost purely imaginary (an arbitrary global phase), so keep both parts.
    ex = (data[:, 2] + 1j * data[:, 3]).reshape(len(y), len(x))
    ey = (data[:, 4] + 1j * data[:, 5]).reshape(len(y), len(x))
    return x, y, ex, ey


def wall_jump(x, ex_row):
    """Ratio of |Ex| just outside to just inside the core's side wall at x = +W_CORE/2.

    Ex is normal to that wall, so the normal displacement n²·Ex is continuous and the field
    must jump by (n_Si/n_ox)². Each side is extrapolated to the wall from its own three
    nearest grid points: a straight line inside (the field is a smooth cosine there) and an
    exponential outside (it decays).
    """
    wall = W_CORE / 2
    a = np.abs(ex_row)
    inside = np.where((x > 0) & (x < wall))[0][-3:]
    outside = np.where(x > wall)[0][:3]
    e_in = np.polyval(np.polyfit(x[inside], a[inside], 1), wall)
    e_out = np.exp(np.polyval(np.polyfit(x[outside], np.log(a[outside]), 1), wall))
    return e_out / e_in


def rep_mode():
    x, y, ex, ey = load_mode("mfield_Tmin.csv")
    _, _, ex_hot, ey_hot = load_mode("mfield_Tmax.csv")
    amp = np.abs(ex) / np.abs(ex).max()

    fig = plt.figure(figsize=(W, 3.05))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.75, 1], width_ratios=[1.12, 1, 1.05],
                            hspace=0.14, wspace=0.52, left=0.075, right=0.985, top=0.93, bottom=0.14)

    # (a) |Ex| of the TE0 mode, and a cut through the core centre
    a = fig.add_subplot(grid[0, 0])
    cmap = LinearSegmentedColormap.from_list("ink_blue", ["#ffffff", *ps.BLUE_RAMP])
    a.imshow(amp, extent=[x[0], x[-1], y[0], y[-1]], origin="lower", cmap=cmap, vmin=0, vmax=1,
             interpolation="bilinear", aspect="auto")
    a.add_patch(Rectangle((-W_CORE / 2, -H_CORE / 2), W_CORE, H_CORE, fill=False, ec=ps.INK, lw=0.8))
    a.set_xlim(-1, 1)
    a.set_ylim(-0.62, 0.62)
    a.grid(False)
    a.set_ylabel("y (µm)")
    a.tick_params(labelbottom=False)
    a.text(0, 0.48, "$|E_x|$ of the TE$_0$ mode", ha="center", fontsize=SMALL, color=ps.INK2)
    a.text(0.29, -0.2, "Si core\n500 × 220 nm", fontsize=6.5, color=ps.INK2, va="top")
    ps.panel_label(a, "a", x=-0.12)

    cut = fig.add_subplot(grid[1, 0], sharex=a)
    mid = np.argmin(np.abs(y))                       # the row through the core centre, y = 0
    cut.plot(x, amp[mid], color=ps.INK, lw=1.2)
    for wall in (-W_CORE / 2, W_CORE / 2):
        cut.axvline(wall, color=ps.AXIS, lw=0.8)
    jump = wall_jump(x, ex[mid])
    cut.set_ylim(0, 1.12)
    cut.set_yticks([0, 0.5, 1])
    cut.set_xlabel("x (µm)")
    cut.set_ylabel("$|E_x|$ at y = 0")
    cut.text(-0.97, 1.04, f"wall jump ×{jump:.1f}", fontsize=6.3, color=ps.INK2, va="top")
    cut.text(0.97, 1.04, f"$(n_{{Si}}/n_{{ox}})^2$ = {(N_SI / N_OX) ** 2:.1f}", fontsize=6.3, color=ps.INK2,
             va="top", ha="right")

    # (b) effective index against temperature, and what is left after a straight-line fit
    _, t = read_comsol(DATA_DIR / "neff_vs_T.csv")
    temp, neff = t[:, 0], t[:, 1]
    slope, intercept = np.polyfit(temp, neff, 1)
    resid = neff - (slope * temp + intercept)
    b = fig.add_subplot(grid[0, 1])
    b.plot(temp, slope * temp + intercept, color=ps.MUTED, lw=1.2, zorder=1)
    b.plot(temp, neff, "o", color=ps.INK, ms=3.5, zorder=2)
    b.set_ylabel("effective index $n_{eff}$")
    b.tick_params(labelbottom=False)
    note(b, 0.05, 0.93, f"slope {slope * 1e4:.3f} × 10$^{{-4}}$ K$^{{-1}}$\n"
         f"= {slope / 1.86e-4:.2f} × bulk Si", va="top")
    ps.panel_label(b, "b", x=-0.16)
    r = fig.add_subplot(grid[1, 1], sharex=b)
    r.axhline(0, color=ps.MUTED, lw=0.8)
    r.plot(temp, resid * 1e6, "o-", color=ps.INK, ms=3, lw=0.9)
    r.set_ylim(-5.5, 5.5)
    r.set_xlabel("temperature (K)")
    r.set_ylabel("residual (10$^{-6}$)")
    note(r, 0.5, 0.86, f"largest {np.abs(resid).max() * 1e6:.1f} × 10$^{{-6}}$", ha="center", va="top", fontsize=6.5)

    # Convergence runs: 1 base; 2, 3 finer meshes; 4, 5 larger boxes
    _, conv = read_comsol(DATA_DIR / "n_eff_sweepcheck.csv")
    n_conv = conv[:, -1]
    mesh_spread = np.abs(n_conv[1:3] - n_conv[0]).max()
    box_spread = np.abs(n_conv[3:5] - n_conv[0]).max()

    # (c) single-emitter pattern across the array, cold and 100 K hotter
    c = fig.add_subplot(grid[:, 2])
    ep_cold = element_pattern()
    ep_hot = element_pattern(file="mfield_Tmax.csv")
    deg = np.linspace(0, 60, 601)
    u = np.sin(np.radians(deg))
    edge = np.degrees(np.arcsin(WAVELENGTH_UM / (2 * PITCH_UM)))     # single-beam field of view
    c.axvspan(0, edge, color=ps.GRID, alpha=0.55, lw=0)
    c.plot(deg, ep_cold(u), color=ps.INK, lw=1.5, label="293 K")
    c.plot(deg, ep_hot(u), color=ps.HEAT, lw=1.1, ls=(0, (4, 2)), label="393 K")
    for angle in (20, edge):
        value = float(ep_cold(np.sin(np.radians(angle))))
        c.plot(angle, value, "o", color=ps.INK, ms=3.5)
        c.annotate(f"{value:.2f} at {angle:.3g}°", (angle, value), xytext=(angle + 3, value + (0.045 if angle == 20 else -0.01)),
                   fontsize=SMALL, color=ps.INK2, va="center")
    change = np.abs(ep_hot(u) - ep_cold(u)).max()
    c.set_xlim(0, 60)
    c.set_ylim(0, 1.05)
    c.set_xlabel("angle across the array (deg)")
    c.set_ylabel("single-emitter pattern")
    c.legend(loc="lower left", handlelength=1.8)
    c.text(edge / 2, 0.40, "one-beam\nfield", ha="center", fontsize=6.5, color=ps.INK2)
    note(c, 0.97, 0.24, f"largest change\nover 0–60°: {change * 100:.1f}%", ha="right", va="bottom", fontsize=6.5)
    ps.panel_label(c, "c", x=-0.14, y=1.01)

    # Numbers for the notes
    energy_ex = np.sum(np.abs(ex) ** 2) / np.sum(np.abs(ex) ** 2 + np.abs(ey) ** 2)
    overlap = np.abs(np.vdot(ex, ex_hot)) ** 2 / (np.vdot(ex, ex).real * np.vdot(ex_hot, ex_hot).real)
    i_max = np.argmax(np.abs(ep_hot(u) - ep_cold(u)))
    wide = np.sin(np.radians(np.linspace(0, 90, 901)))
    j_max = np.argmax(np.abs(ep_hot(wide) - ep_cold(wide)))
    print(f"   n_eff(293.15 K) = {neff[0]:.5f}; slope {slope:.4e} /K (constants.py: {DNEFF_DT:.3e}); "
          f"max residual {np.abs(resid).max():.2e}; ratio to bulk {slope / 1.86e-4:.3f}")
    print(f"   phase per K per mm: {K0 * 1e3 * slope:.3f} rad; 2π needs {2 * np.pi / (K0 * 1e3 * slope):.2f} K·mm")
    print(f"   convergence: finer meshes move n_eff by {mesh_spread:.1e}, larger boxes by {box_spread:.1e}")
    print(f"   wall jump {jump:.3f} against (n_Si/n_ox)² = {(N_SI / N_OX) ** 2:.3f}")
    print(f"   EP(20°) = {ep_cold(np.sin(np.radians(20))):.4f}, EP({edge:.2f}°) = {ep_cold(WAVELENGTH_UM / (2 * PITCH_UM)):.4f}")
    print(f"   EP change 293→393 K: {change * 100:.2f}% at {deg[i_max]:.0f}° (0–60°); "
          f"{np.abs(ep_hot(wide) - ep_cold(wide)).max() * 100:.2f}% at {np.degrees(np.arcsin(wide[j_max])):.0f}° (0–90°)")
    print(f"   |Ex|² share of |Ex|²+|Ey|² on the grid: {energy_ex * 100:.2f}%; Ex overlap cold/hot: {overlap:.6f}")
    ps.save(fig, "rep_mode")


# ================================================================ rep_coupling
def rep_coupling():
    _, even = read_comsol(DATA_DIR / "neff_even_vs_d.csv")
    _, odd = read_comsol(DATA_DIR / "neff_odd_vs_d.csv")
    d = even[:, 0]
    # The two supermodes are solved on the same half-model mesh, so the mesh error (about
    # 1e-5 on each index) cancels in the difference.
    dn = even[:, -1] - odd[:, -1]
    slope, intercept = np.polyfit(d, np.log(dn), 1)          # ln Δn falls linearly with spacing
    d_fit = np.linspace(0.7, 3.1, 200)
    dn_fit = np.exp(slope * d_fit + intercept)
    lc = WAVELENGTH_UM / (2 * dn) * 1e-3                     # coupling length, mm
    lc_fit = WAVELENGTH_UM / (2 * dn_fit) * 1e-3

    # Power moved to the neighbour after a length L is sin²(πL/2L_c). Keeping it under 1%
    # over 2 mm of side-by-side waveguide needs L_c above this:
    run_mm, transfer = 2.0, 0.01
    lc_min = np.pi * run_mm / (2 * np.arcsin(np.sqrt(transfer)))
    d_min = (np.log(WAVELENGTH_UM / (2 * lc_min * 1e3)) - intercept) / slope
    lc_design = float(lc[np.isclose(d, PITCH_UM)][0])

    fig, (a, b) = plt.subplots(1, 2, figsize=(W, 2.55))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.92, bottom=0.17, wspace=0.30)

    a.plot(d_fit, dn_fit, color=ps.MUTED, lw=1.2, label=f"fit: decay {-slope:.2f} µm$^{{-1}}$")
    a.plot(d, dn, "o", color=ps.INK, ms=4, label="COMSOL, $n_{even} - n_{odd}$")
    a.set_yscale("log")
    a.set_xlabel("centre spacing d (µm)")
    a.set_ylabel("supermode index split Δn")
    a.set_xlim(0.7, 3.1)
    a.legend(loc="upper right", handlelength=1.6)
    note(a, 0.05, 0.07, f"largest miss of the fit:\n×{np.exp(np.abs(np.log(dn) - (slope * d + intercept)).max()):.2f}",
         fontsize=6.5)
    ps.panel_label(a, "a", x=-0.10)

    b.plot(d_fit, lc_fit, color=ps.MUTED, lw=1.2)
    b.plot(d, lc, "o", color=ps.INK, ms=4)
    b.axhline(lc_min, color=ACCENT, lw=1.0, ls=(0, (4, 2)))
    b.text(3.05, lc_min / 1.8, f"under 1% transferred over 2 mm:\n$L_c$ > {lc_min:.0f} mm, d > {d_min:.2f} µm",
           ha="right", va="top", fontsize=SMALL, color=ACCENT)
    b.plot(PITCH_UM, lc_design, "o", ms=8, mfc="none", mec=ps.INK, mew=1.0)
    b.annotate(f"design pitch 2 µm\n$L_c$ = {lc_design / 1e3:.1f} m", (PITCH_UM, lc_design), xytext=(1.0, lc_design * 20),
               fontsize=SMALL, color=ps.INK2, arrowprops=dict(arrowstyle="-", color=ps.MUTED, lw=0.7,
                                                              shrinkA=1, shrinkB=5))
    b.set_yscale("log")
    b.set_xlabel("centre spacing d (µm)")
    b.set_ylabel("coupling length $L_c$ = λ/2Δn (mm)")
    b.set_xlim(0.7, 3.1)
    ps.panel_label(b, "b", x=-0.10)

    print("   d (µm), Δn, L_c:")
    for di, dni, lci in zip(d, dn, lc):
        print(f"     {di:.2f}  {dni:.3e}  {lci:.4g} mm")
    print(f"   fit: ln Δn slope {slope:.3f} per µm (all {len(d)} points); "
          f"per 0.25 µm Δn falls ×{np.exp(-slope * 0.25):.2f}")
    print(f"   1% over 2 mm needs L_c ≥ {lc_min:.2f} mm, reached at d = {d_min:.3f} µm")
    _, iso = read_comsol(DATA_DIR / "neff_vs_T.csv")
    n_iso = iso[0, 1]                                        # the isolated strip, Model A
    close = d <= 1.75     # beyond this the split is below the mesh-to-mesh scatter of each index
    print(f"   even above and odd below the isolated index {n_iso:.7f} for d ≤ 1.75 µm: "
          f"{bool(np.all(even[close, -1] > n_iso) and np.all(odd[close, -1] < n_iso))}; "
          f"for all d: {bool(np.all(even[:, -1] > n_iso) and np.all(odd[:, -1] < n_iso))}")
    for lo, hi in ((0.8, 2.0), (1.0, 3.0), (1.25, 3.0)):
        sel = (d >= lo) & (d <= hi)
        print(f"   fit over d = {lo}–{hi} µm: decay {-np.polyfit(d[sel], np.log(dn[sel]), 1)[0]:.3f} per µm")
    ps.save(fig, "rep_coupling")


# =========================================================== rep_kernel_checks
def rep_kernel_checks():
    _, data = read_comsol(DATA_DIR / "kernel_checks.csv")
    x = data[:, 0]
    runs = data[:, 2:]                                   # six columns: the base and five variants
    base = runs[:, 0]
    centre = np.argmin(np.abs(x))
    peak = base[centre]
    labels = ["base", "oxide cores", "3 µm buried oxide", "200 µm substrate", "±300 µm domain", "finer mesh"]
    colors = [ps.COLOR["standard"], ps.INK, ACCENT, BROWN, ps.INK, ACCENT]
    styles = ["-", (0, (1, 1.2)), "-", (0, (4, 2)), "-", (0, (4, 2))]
    array = np.abs(x) <= 32.0

    fig = plt.figure(figsize=(W, 4.75))
    outer = fig.add_gridspec(2, 2, height_ratios=[1, 0.95], width_ratios=[1, 1], hspace=0.42, wspace=0.36,
                             left=0.085, right=0.985, top=0.95, bottom=0.095)

    # (a) the six runs, centre to the domain edge
    a = fig.add_subplot(outer[0, 0])
    right = x >= 0
    for k in [2, 3, 1, 4, 5, 0]:
        a.plot(x[right], runs[right, k], color=colors[k], ls=styles[k], lw=1.6 if k == 0 else 1.0,
               label=labels[k], zorder=3 if k == 0 else 2)
    a.set_yscale("log")
    a.set_xlim(0, 150)
    a.set_xlabel("distance from the heater (µm)")
    a.set_ylabel("core-height rise (K per W/m)")
    handles, names = a.get_legend_handles_labels()
    order = [5, 2, 0, 1, 3, 4]
    a.legend([handles[i] for i in order], [names[i] for i in order], loc="upper right", fontsize=6.8,
             handlelength=2.0, labelspacing=0.3)
    ps.panel_label(a, "a", x=-0.10)

    # (b) difference from the base run across the array, in % of the base peak
    inner = outer[0, 1].subgridspec(2, 1, hspace=0.16)
    b1 = fig.add_subplot(inner[0])
    b2 = fig.add_subplot(inner[1], sharex=b1)
    diff = (runs - base[:, None]) / peak * 100
    for k in (1, 2, 3):
        b1.plot(x[array], diff[array, k], color=colors[k], ls=styles[k], lw=1.1)
    b1.set_ylim(-9, 31)
    b1.set_yticks([0, 10, 20, 30])
    b1.tick_params(labelbottom=False)
    b1.text(-31, 25, f"3 µm buried oxide:\n+{diff[centre, 2]:.0f}% at the heater", fontsize=6.5, color=ACCENT, va="center")
    b1.text(-31, 4.6, f"200 µm substrate: +{diff[array, 3].mean():.1f}%", fontsize=6.5, color=BROWN)
    b1.text(31, -6.3, f"oxide cores: +{diff[centre, 1]:.1f}% at the heater", fontsize=6.5, color=ps.INK, ha="right",
            va="center")
    for k in (4, 5):
        b2.plot(x[array], diff[array, k], color=colors[k], ls=styles[k], lw=1.1)
    b2.axhline(0, color=ps.MUTED, lw=0.6)
    b2.set_ylim(-0.075, 0.05)
    b2.set_yticks([-0.05, 0, 0.05])
    b2.text(-31, 0.033, f"finer mesh: at most {np.abs(diff[array, 5]).max():.2f}%", fontsize=6.5, color=ACCENT, va="center")
    b2.text(-31, -0.062, f"±300 µm domain: at most {np.abs(diff[array, 4]).max():.2f}%", fontsize=6.5, color=ps.INK,
            va="center")
    b2.set_xlabel("position across the array (µm)")
    b2.set_xlim(-32, 32)
    b2.set_xticks([-32, -16, 0, 16, 32])
    b1.set_ylabel("change (%)")
    b2.set_ylabel("change (%)")
    fig.align_ylabels([b1, b2])
    ps.panel_label(b1, "b", x=-0.17)

    # (c) undercut stack: the profile changes shape as the heater moves toward the cavity wall
    c = fig.add_subplot(outer[1, 0])
    uc = UndercutKernel()
    xs = np.linspace(-60, 60, 1201)
    c.axvspan(-60, -40, color=ps.GRID, alpha=0.6, lw=0)
    c.axvspan(40, 60, color=ps.GRID, alpha=0.6, lw=0)
    heaters = [0, 8, 16, 24, 32]
    for i, xh in enumerate(heaters):
        col = shade(ps.COLOR["undercut"], -0.45 + 0.9 * i / (len(heaters) - 1))
        c.plot(xs, uc(xs, xh), color=col, lw=1.4)
        c.plot(xh, uc(xh, xh), "v", color=col, ms=3.5)
        c.text(xh + (0 if i else 1.5), uc(xh, xh) + 0.10, f"{xh}" if i else "heater at 0", ha="center" if i else "right",
               va="bottom", fontsize=6.5, color=ps.INK2)
    # What a shift-invariant kernel would give: the centred profile slid to 32 µm
    c.plot(xs, uc(xs - 32, 0), color=ps.MUTED, lw=0.9, ls=(0, (4, 2)))
    c.text(-38.5, 1.9, "dashed: the centred\nprofile slid to 32 µm", fontsize=6.5, color=ps.INK2, ha="left")
    c.text(-50, 3.02, "substrate", ha="center", fontsize=6.5, color=ps.INK2)
    c.text(50, 3.02, "substrate", ha="center", fontsize=6.5, color=ps.INK2)
    c.set_xlim(-60, 60)
    c.set_ylim(0, 3.25)
    c.set_xlabel("position across the array (µm)")
    c.set_ylabel("core-height rise (K per W/m)")
    ps.panel_label(c, "c", x=-0.10)

    # (d) Model C′: mode solved in the real, uneven temperature field, against the shortcut
    # Δn_eff = dn_eff/dT × (core temperature rise from the kernel)
    d = fig.add_subplot(outer[1, 1])
    _, cp = read_comsol(DATA_DIR / "cprime_neff.csv")
    std = StandardKernel()
    n_cold = cp[0, -1]
    offsets, ratios = [], []
    for x_h, q_l, n_hot in cp[1:, [0, 1, -1]]:
        rise = std(0.0, x_h) * q_l                    # K: kernel (K per W/m) × heater power (W/m)
        offsets.append(x_h)
        ratios.append((n_hot - n_cold) / (DNEFF_DT * rise))
    d.axhspan(0.995, 1.005, color=ps.GRID, alpha=0.7, lw=0)
    d.axhline(1, color=ps.MUTED, lw=0.8)
    d.plot(offsets, ratios, "o", color=ps.INK, ms=5)
    for xo, ro in zip(offsets, ratios):
        d.text(xo, ro - 0.0011, f"{ro:.3f}", ha="center", va="top", fontsize=SMALL, color=ps.INK2)
    d.set_xlim(-1, 5)
    d.set_ylim(0.990, 1.010)
    d.set_xticks(offsets)
    d.set_xlabel("heater offset from the core (µm)")
    d.set_ylabel("solved Δ$n_{eff}$ ÷ shortcut")
    note(d, 0.03, 0.95, "coupled heat-and-mode solve, 50 W/m", va="top", fontsize=6.5)
    note(d, 0.97, 0.80, "±0.5% band", ha="right", va="top", fontsize=6.5, color=ps.MUTED)
    ps.panel_label(d, "d", x=-0.17)

    # Numbers for the notes
    def half_width(curve):
        i = centre + np.argmax(curve[centre:] < curve[centre] / 2)
        return np.interp(curve[centre] / 2, [curve[i], curve[i - 1]], [x[i], x[i - 1]])

    print("   run: K(0) [K per W/m], half-width [µm], change at the heater, largest |change| over ±32 µm")
    for k, lab in enumerate(labels):
        print(f"     {lab:18s} {runs[centre, k]:.4f}  {half_width(runs[:, k]):.2f}  {diff[centre, k]:+.2f}%  "
              f"{np.abs(diff[array, k]).max():.2f}%")
    print(f"   200 µm substrate offset over ±32 µm: {diff[array, 3].min():.3f}% to {diff[array, 3].max():.3f}% "
          f"(spread {np.ptp(diff[array, 3]) / 100:.1e} of the peak)")
    print(f"   oxide cores beyond 10 µm: at most {np.abs(diff[(np.abs(x) >= 10) & array, 1]).max():.2f}% of the peak")
    print(f"   fractions of the peak at 2, 4, 6, 8 µm: " + ", ".join(f"{std(v, 0) / peak:.2f}" for v in (2, 4, 6, 8)))
    print(f"   tail at 20–32 µm: {std(32, 0):.4f}–{std(20, 0):.4f} K per W/m")
    _, energy = read_comsol(DATA_DIR / "checks_energy.csv")
    print(f"   energy balance, six check runs: largest |1 − flux| = {np.abs(energy[:, -1] - 1).max():.1e} W/m")
    k00 = uc(0, 0)
    print(f"   undercut K(0;0) = {k00:.3f}; fractions at 2, 8, 16, 32 µm: "
          + ", ".join(f"{uc(v, 0) / k00:.2f}" for v in (2, 8, 16, 32)))
    print(f"   undercut self-heating at x_h = 32 µm: {uc(32, 32) / k00:.2f} of the centre value")
    for shift in (5, 10, 20, 32):
        miss = np.abs(uc(xs, shift) - uc(xs - shift, 0)).max() / k00
        print(f"   slid centred profile misses the true one at x_h = {shift} µm by {miss * 100:.0f}% of the peak")
    print("   C′ solved/shortcut at offsets " + ", ".join(f"{o:g} µm: {r_:.4f}" for o, r_ in zip(offsets, ratios)))
    ps.save(fig, "rep_kernel_checks")


# ================================================================= rep_grating
# Model E parameters (build guide): equivalent 2D core, 20 nm teeth, 50 periods
N_C2D, DNDT_C2D = 3.0604, 2.116e-4
DNDT_OX = 0.95e-5
ETCH, FILL = 0.020, 0.5
L_IN, N_PERIODS = 5.0, 50


def slab_index(n_core, n_clad, thickness, wavelength):
    """Effective index of the fundamental TE mode of a symmetric slab (lengths in µm).

    Inside the core the field is a cosine with transverse wavenumber κ; outside it decays
    at rate γ. Matching them at the two faces gives κ·h = 2·arctan(γ/κ) for the lowest mode.
    """
    k0 = 2 * np.pi / wavelength

    def mismatch(n):
        kappa = k0 * np.sqrt(n_core**2 - n**2)
        gamma = k0 * np.sqrt(n**2 - n_clad**2)
        return kappa * thickness - 2 * np.arctan(gamma / kappa)

    return brentq(mismatch, n_clad + 1e-9, n_core - 1e-9)


def grating_index(wavelength, temp=T0):
    """Mean effective index along the grating: the full and the etched slab, weighted by
    the fill factor. This is what phase matching sees for shallow teeth."""
    n_core = N_C2D + DNDT_C2D * (temp - T0)
    n_clad = N_OX + DNDT_OX * (temp - T0)
    full = slab_index(n_core, n_clad, H_CORE, wavelength)
    etched = slab_index(n_core, n_clad, H_CORE - ETCH, wavelength)
    return FILL * full + (1 - FILL) * etched


def grating_angle(period, wavelength, temp=T0):
    """Grating equation, angle from vertical in air: sin θ = n_eff − λ/Λ (degrees)."""
    return np.degrees(np.arcsin(grating_index(wavelength, temp) - wavelength / period))


def read_grating():
    """Emission angle and leakage of each Model E run, from the field on the monitor line.

    COMSOL uses exp(+iωt), so a wave travelling toward +x has a phase that falls with x:
    sin θ = −(dφ/dx)·λ/2π. The fit runs from x = 8 µm (past the start of the grating at
    5 µm) to 1 µm before the grating ends, where the emitted wave is clean.
    """
    _, power = read_comsol(DATA_DIR / "grating_power.csv")
    _, line = read_comsol(DATA_DIR / "grating_line.csv")
    x = line[:, 0]
    runs = []
    for r, (period, lam_nm, temp, _, refl, trans, up, down) in enumerate(power):
        lam = lam_nm / 1000
        ez = line[:, 2 + 3 * r] + 1j * line[:, 3 + 3 * r]
        flux = line[:, 4 + 3 * r]                                   # upward Poynting flux, W/m²
        window = (x >= 8) & (x <= L_IN + N_PERIODS * period - 1) & np.isfinite(ez.real)
        phase_slope = np.polyfit(x[window], np.unwrap(np.angle(ez[window])), 1)[0]
        decay = -np.polyfit(x[window], np.log(flux[window]), 1)[0]  # power leakage rate, 1/µm
        runs.append(dict(period=period, lam=lam, temp=temp, refl=refl, trans=trans, up=up, down=down,
                         angle=np.degrees(np.arcsin(-phase_slope * lam / (2 * np.pi))),
                         decay=decay, x=x, flux=flux, window=window))
    return runs


def rep_grating():
    runs = read_grating()
    design = runs[0]
    by_period = sorted([design] + runs[1:7], key=lambda q: q["period"])
    by_lam = sorted([design] + runs[7:13], key=lambda q: q["lam"])
    by_temp = [design] + runs[13:17]
    combined = runs[17]

    # Three panels: the emission angle against period, wavelength and temperature. The upward-flux
    # profile is only used for the printed leakage length below, not drawn.
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(W, 2.5))
    fig.subplots_adjust(left=0.075, right=0.99, top=0.92, bottom=0.18, wspace=0.34)

    def rate(xv, yv):
        return np.polyfit(xv, yv, 1)[0]

    # (a) period
    period = np.array([q["period"] for q in by_period])
    angle_p = np.array([q["angle"] for q in by_period])
    fine = np.linspace(0.655, 0.725, 100)
    a.plot(fine * 1000, [grating_angle(p, WAVELENGTH_UM) for p in fine], color=ps.MUTED, lw=1.2,
           label="grating equation")
    a.plot(period * 1000, angle_p, "o", color=ps.INK, ms=4, label="COMSOL (2D)")
    a.set_xlabel("grating period Λ (nm)")
    a.set_ylabel("emission angle in air (deg)")
    a.legend(loc="upper left", handlelength=1.6)
    note(a, 0.97, 0.06, f"{rate(period * 1000, angle_p) * 10:+.2f}° per 10 nm", ha="right")
    ps.panel_label(a, "a", x=-0.10)

    # (b) wavelength: the 2D model, and the same angles corrected to the real strip
    lam = np.array([q["lam"] for q in by_lam])
    angle_l = np.array([q["angle"] for q in by_lam])
    fine = np.linspace(1.495, 1.605, 100)
    b.plot(fine * 1000, [grating_angle(design["period"], w) for w in fine], color=ps.MUTED, lw=1.2,
           label="grating equation, 2D slab")
    b.plot(lam * 1000, angle_l, "o", color=ps.INK, ms=4, label="COMSOL (2D)")
    # The 2D slab stands in for the strip: same index at 1550 nm, but a different change with
    # wavelength. Shift sin θ by the difference between the two index changes.
    _, disp = read_comsol(DATA_DIR / "neff_vs_lambda.csv")
    lam_strip, n_strip = disp[:, 0] / 1000, disp[:, -1]
    ref = np.argmin(np.abs(lam_strip - WAVELENGTH_UM))
    n_group = n_strip[ref] - WAVELENGTH_UM * np.gradient(n_strip, lam_strip)[ref]
    slab_ref = slab_index(N_C2D, N_OX, H_CORE, WAVELENGTH_UM)

    def strip_angle(sin_2d, w):
        strip_change = np.interp(w, lam_strip, n_strip) - n_strip[ref]
        slab_change = slab_index(N_C2D, N_OX, H_CORE, w) - slab_ref
        return np.degrees(np.arcsin(sin_2d + strip_change - slab_change))

    angle_strip = np.array([strip_angle(np.sin(np.radians(t)), w) for t, w in zip(angle_l, lam)])
    b.plot(lam * 1000, angle_strip, "s-", color=ACCENT, ms=3.5, lw=1.0, label="corrected to the real strip")
    b.set_xlabel("wavelength (nm)")
    b.set_ylabel("emission angle in air (deg)")
    b.legend(loc="upper right", handlelength=1.4, fontsize=6.3)
    near = np.abs(lam - WAVELENGTH_UM) < 0.015                       # the three runs around 1550 nm
    note(b, 0.04, 0.17, f"2D: {rate(lam[near] * 1000, angle_l[near]):+.3f}°/nm", va="bottom")
    note(b, 0.04, 0.06, f"strip: {rate(lam[near] * 1000, angle_strip[near]):+.3f}°/nm",
         va="bottom", color=ACCENT)
    ps.panel_label(b, "b", x=-0.10)

    # (c) temperature
    temp = np.array([q["temp"] for q in by_temp])
    angle_t = np.array([q["angle"] for q in by_temp])
    fine = np.linspace(T0 - 3, T0 + 103, 50)
    c.plot(fine - T0, [grating_angle(design["period"], WAVELENGTH_UM, t) for t in fine], color=ps.MUTED, lw=1.2,
           label="grating equation")
    c.plot(temp - T0, angle_t, "o", color=ps.INK, ms=4, label="COMSOL (2D)")
    c.set_xlabel("uniform temperature rise (K)")
    c.set_ylabel("emission angle in air (deg)")
    c.legend(loc="upper left", handlelength=1.6)
    note(c, 0.97, 0.06, f"{rate(temp, angle_t) * 100:+.2f}° per 100 K", ha="right")
    ps.panel_label(c, "c", x=-0.10)

    emitted = design["up"] + design["down"]

    # Numbers for the notes
    print(f"   slab check: 220 nm slab n_eff = {slab_ref:.5f} (strip: {n_strip[ref]:.5f}); "
          f"grating mean index {grating_index(WAVELENGTH_UM):.5f}")
    print(f"   design point: COMSOL {design['angle']:.3f}°, grating equation {grating_angle(design['period'], WAVELENGTH_UM):.3f}°")
    miss = max(abs(q["angle"] - grating_angle(q["period"], q["lam"], q["temp"])) for q in runs)
    print(f"   largest |COMSOL − grating equation| over all 18 runs: {miss:.3f}°")
    print(f"   rates (COMSOL / grating equation): period {rate(period * 1000, angle_p) * 10:+.3f} / "
          f"{rate(period * 1000, [grating_angle(p, WAVELENGTH_UM) for p in period]) * 10:+.3f} ° per 10 nm; "
          f"wavelength (1540–1560 nm) {rate(lam[near] * 1000, angle_l[near]):+.4f} / "
          f"{rate(lam[near] * 1000, [grating_angle(design['period'], w) for w in lam[near]]):+.4f} ° per nm; "
          f"temperature {rate(temp, angle_t) * 100:+.3f} / "
          f"{rate(temp, [grating_angle(design['period'], WAVELENGTH_UM, t) for t in temp]) * 100:+.3f} ° per 100 K")
    print(f"   wavelength rate over the whole 1500–1600 nm: 2D {rate(lam * 1000, angle_l):+.4f}, strip {rate(lam * 1000, angle_strip):+.4f} ° per nm")
    print(f"   strip-corrected angles: " + ", ".join(f"{w * 1000:.0f} nm: {t2:.2f}→{t3:.2f}°" for w, t2, t3 in zip(lam, angle_l, angle_strip)))
    print(f"   strip group index at 1550 nm: {n_group:.3f}; strip rate at 1550: {rate(lam[near] * 1000, angle_strip[near]):+.4f} ° per nm; "
          f"sweep 1500–1600 nm: {angle_strip[0] - angle_strip[-1]:.1f}°")
    slab_disp = (slab_index(N_C2D, N_OX, H_CORE, 1.56) - slab_index(N_C2D, N_OX, H_CORE, 1.54)) / 0.02
    print(f"   slab dn/dλ at 1550 nm: {slab_disp:.3f} per µm; slab group index {slab_ref - WAVELENGTH_UM * slab_disp:.2f}")
    shifts = [by_period[-1]["angle"] - design["angle"], by_lam[-1]["angle"] - design["angle"], by_temp[-1]["angle"] - design["angle"]]
    print(f"   run 18 (Λ 0.72 µm, 1600 nm, +100 K): {combined['angle']:.3f}°; design + separate shifts "
          f"({shifts[0]:+.2f}, {shifts[1]:+.2f}, {shifts[2]:+.2f}) = {design['angle'] + sum(shifts):.3f}°")
    sums = [q["refl"] + q["trans"] + q["up"] + q["down"] for q in runs]
    print(f"   power sum R+T+up+down: {min(sums):.4f}–{max(sums):.4f}; largest R {max(q['refl'] for q in runs):.1e}")
    shares = [q["up"] / (q["up"] + q["down"]) for q in by_lam]
    print(f"   up share: design {design['up'] / emitted:.3f}; across 1500–1600 nm {min(shares):.2f}–{max(shares):.2f}")
    print(f"   power leakage rate at the design point {design['decay']:.5f} per µm (field {design['decay'] / 2:.5f}); "
          f"1/e power length {1 / design['decay']:.1f} µm; across wavelength field rate "
          f"{by_lam[-1]['decay'] / 2:.4f}–{by_lam[0]['decay'] / 2:.4f} per µm")
    print(f"   leakage cross-check: exp(−{N_PERIODS * design['period']:.1f} µm × rate) = "
          f"{np.exp(-N_PERIODS * design['period'] * design['decay']):.3f} against the port's transmission {design['trans']:.3f}")
    ps.save(fig, "rep_grating")


# ========================================================= rep_farfield_checks
def sinc2_constants():
    """Half-power width and first side lobe of sin²(x)/x², the pattern of a uniform aperture.

    For an aperture of width N·d the pattern is sinc²(π·N·d·u/λ), so its full width at half
    maximum in u is (this factor)·λ/(N·d). The factor is the textbook 0.886.
    """
    half = brentq(lambda v: (np.sin(v) / v) ** 2 - 0.5, 0.5, 2.0)
    v = np.linspace(np.pi, 2 * np.pi, 20001)[1:-1]
    return 2 * half / np.pi, -10 * np.log10(((np.sin(v) / v) ** 2).max())


def rep_farfield_checks():
    xe = emitter_positions()
    ep = element_pattern()
    n, d = N_EMITTERS, PITCH_UM
    width_factor, sinc_lobe_db = sinc2_constants()

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(W, 2.6), gridspec_kw=dict(width_ratios=[1, 1.28, 0.95]))
    fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.185, wspace=0.42)

    # (a) all emitters in phase: the engine against the closed-form sum
    engine = far_field(np.zeros(n), xe, ep)
    psi = K0 * d * U_GRID                                  # phase step between neighbours seen from angle u
    in_phase = np.abs(np.sin(psi / 2)) < 1e-9             # broadside and the grating lobes: 0/0
    closed = np.ones_like(psi)                             # the limit where every term adds in phase
    closed[~in_phase] = (np.sin(n * psi[~in_phase] / 2) / (n * np.sin(psi[~in_phase] / 2))) ** 2
    closed *= ep(U_GRID)
    m = beam_metrics(U_GRID, engine, ep)
    deg = np.degrees(np.arcsin(U_GRID))
    floor = 1e-6
    a.plot(deg, 10 * np.log10(np.maximum(engine, floor)), color=ps.INK, lw=1.3, label="far-field engine")
    a.plot(deg, 10 * np.log10(np.maximum(closed, floor)), color=ps.HEAT, lw=0.9, ls=(0, (2, 2)), label="closed form")
    a.axhline(-m["sls_db"], color=ps.MUTED, lw=0.7, ls=(0, (4, 2)))
    a.text(7.8, -m["sls_db"] + 0.8, f"−{m['sls_db']:.2f} dB", ha="right", va="bottom", fontsize=6.5, color=ps.INK2)
    half_deg = m["fwhm_deg"] / 2
    a.plot([-half_deg, half_deg], [-3.01, -3.01], color=ACCENT, lw=1.6)
    a.text(half_deg + 0.5, -3.0, f"FWHM {m['fwhm_deg']:.2f}°", fontsize=6.5, color=ACCENT, va="center")
    a.set_xlim(-8, 8)
    a.set_ylim(-40, 16)
    a.set_yticks([-40, -30, -20, -10, 0])
    a.set_xlabel("angle (deg)")
    a.set_ylabel("intensity (dB)")
    a.legend(loc="upper left", fontsize=6.5, handlelength=1.6, borderaxespad=0.2)
    ps.panel_label(a, "a", x=-0.13)

    # (b) uniform ramps: the beam must land at sin θ = Δφ·λ/(2π·d)
    steps = [0.5, 1.5, np.pi]
    names = ["0.5 rad", "1.5 rad", "π rad"]
    cols = [shade(ps.INK, -0.55), shade(ps.INK, -0.25), ps.INK]
    b.plot(deg, ep(U_GRID), color=ps.MUTED, lw=0.9, ls=(0, (4, 2)))
    b.text(31, 0.90, "single-emitter\npattern", fontsize=6.3, color=ps.INK2, va="center")
    landed = []
    for step, name, col in zip(steps, names, cols):
        pattern = far_field(step * xe / d, xe, ep)
        b.plot(deg, pattern, color=col, lw=1.0, label=name)
        u_pred = step * WAVELENGTH_UM / (2 * np.pi * d)
        metrics = beam_metrics(U_GRID, pattern, ep)
        landed.append((step, metrics["u_peak"], u_pred, metrics["strehl"]))
        for copy in (u_pred, u_pred - WAVELENGTH_UM / d):          # the beam and its grating-lobe copy
            b.plot(np.degrees(np.arcsin(copy)), 1.085, "v", color=col, ms=3.5, clip_on=False)
    b.legend(loc="center", bbox_to_anchor=(0.825, 0.34), fontsize=6.3, handlelength=0.8, handletextpad=0.4,
             labelspacing=0.25, borderpad=0)
    b.set_xlim(-62, 62)
    b.set_ylim(0, 1.06)
    b.set_xticks([-60, -30, 0, 30, 60])
    b.set_xlabel("angle (deg)")
    b.set_ylabel("intensity ÷ broadside peak")
    note(b, 0.5, 1.085, "▼ predicted: sin θ = Δφ·λ/2πd, and its copy at −λ/d", ha="center", va="bottom", fontsize=6.3)
    ps.panel_label(b, "b", x=-0.10)

    # (c) random phase errors: mean on-axis intensity against the Ruze relation
    rng = np.random.default_rng(0)
    sigmas = np.array([0.1, 0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0])
    trials = 4000
    mean = np.array([np.mean([strehl(rng.normal(0, s, n), xe, 0.0) for _ in range(trials)]) for s in sigmas])
    fine = np.linspace(0, 3.1, 200)
    ruze = np.exp(-fine**2)
    # Finite N adds a floor: the N self-terms of |Σ exp(iε)|² never average away.
    finite = ruze + (1 - ruze) / n
    c.plot(fine, ruze, color=ps.MUTED, lw=0.9, ls=(0, (4, 2)), label="$e^{-σ^2}$")
    c.plot(fine, finite, color=ps.MUTED, lw=1.3, label=f"with the 1/N floor")
    c.plot(sigmas, mean, "o", color=ps.INK, ms=3.8, label="Monte Carlo")
    c.set_yscale("log")
    c.set_ylim(1e-3, 1.6)
    c.set_xlim(0, 3.1)
    c.set_xlabel("phase error σ (rad rms)")
    c.set_ylabel("mean Strehl ratio")
    c.legend(loc="lower left", fontsize=6.3, handlelength=1.3, borderaxespad=0.1)
    ps.panel_label(c, "c", x=-0.14)

    # Numbers for the notes
    expected = np.exp(-sigmas**2) + (1 - np.exp(-sigmas**2)) / n
    print(f"   uniform phase: engine − closed form, largest difference {np.abs(engine - closed).max():.1e}")
    print(f"   beam at {m['angle_deg']:.4f}°; FWHM {m['fwhm_deg']:.3f}° = {m['fwhm_u']:.5f} in sin θ; "
          f"{width_factor:.4f}·λ/(N·d) = {width_factor * WAVELENGTH_UM / (n * d):.5f} (ratio {m['fwhm_u'] / (width_factor * WAVELENGTH_UM / (n * d)):.4f})")
    print(f"   first side lobe {m['sls_db']:.2f} dB (continuous aperture: {sinc_lobe_db:.2f} dB); "
          f"grating lobe {m['grating_lobe_db']:.2f} dB below the beam at ±{np.degrees(np.arcsin(WAVELENGTH_UM / d)):.1f}°; Strehl {m['strehl']:.6f}")
    for step, found, pred, s in landed:
        print(f"   ramp {step:.3f} rad: peak at sin θ = {found:+.4f} ({np.degrees(np.arcsin(found)):+.2f}°), "
              f"predicted {pred:+.4f} ({np.degrees(np.arcsin(pred)):+.2f}°), Strehl {s:.4f}")
    pattern = far_field(np.pi * xe / d, xe, ep)
    u_edge = WAVELENGTH_UM / (2 * d)
    print(f"   ramp π: I(+edge)/I(−edge) = {np.interp(u_edge, U_GRID, pattern) / np.interp(-u_edge, U_GRID, pattern):.5f}; "
          f"each beam {np.interp(u_edge, U_GRID, pattern):.3f} of a broadside beam; edge at {np.degrees(np.arcsin(u_edge)):.2f}°")
    print(f"   Ruze, {trials} trials, N = {n}: largest |Monte Carlo − formula| = {np.abs(mean - expected).max():.4f}")
    for s, mc, ex_ in zip(sigmas, mean, expected):
        print(f"     σ = {s:.1f}: {mc:.4f} against {ex_:.4f}")
    ps.save(fig, "rep_farfield_checks")


FIGURES = {f.__name__: f for f in (rep_mode, rep_coupling, rep_kernel_checks, rep_grating, rep_farfield_checks)}

if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGURES):
        print(name)
        FIGURES[name]()
