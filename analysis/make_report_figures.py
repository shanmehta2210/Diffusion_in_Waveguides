"""Figures drawn only for the report (rep_*), one function per figure.

Run all:   python make_report_figures.py
Run some:  python make_report_figures.py rep_ruze rep_transients

Same rules as make_figures.py: every number plotted or annotated is read from the COMSOL
exports or analysis/results/*.csv, or computed here with the opa/ functions the analyses use.
All drives are plain steps. Each function prints the numbers it computed (they are copied
into figures/REPORT_RESULTS_NOTES.md).
"""
import functools
import itertools
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

import make_figures as mf            # kernels, the three designs, published heaters (fig05)
from opa import membrane1d
from opa import plotstyle as ps
from opa.comsol_io import column_parameters, read_comsol
from opa.constants import DATA_DIR, N_EMITTERS, PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.conventional import ideal_power
from opa.fabrication import ruze
from opa.modes import half_width
from opa.nonlinear import (fit_heating, gain_corrected_drive, model_f_runs, pointing_deg,
                           refit_drive, strehl, true_phase, true_rise)
from opa.placement import evaluate_layout
from opa.steering import STREHL_MIN, steer_one_direction
from opa.thermal import phase_matrix
from opa.transient import StepResponses, settle_time

RESULTS = Path(__file__).parent / "results"
xe = mf.xe
dphi = mf.dphi
U_PER_RAD = WAVELENGTH_UM / (2 * np.pi * PITCH_UM)      # sin θ per rad of ramp step
FOV_DEG = np.degrees(np.arcsin(np.pi * U_PER_RAD))       # 22.8°, edge of the single-beam field
LABEL = {"standard": "standard stack", "undercut": "undercut", "slot": "slot in air",
         "vacuum": "vacuum"}
DASH = (0, (3, 3))


def plain(axis):
    """Plain numbers (0.1, 1, 10) on a log axis."""
    axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))


def crossing(times, y, level):
    """First time y(t) reaches level, interpolated in log t."""
    j = int(np.flatnonzero(y >= level)[0])
    w = (level - y[j - 1]) / (y[j] - y[j - 1])
    return float(np.exp(np.log(times[j - 1]) + w * np.log(times[j] / times[j - 1])))


# ====================================================================== rep_transients
def model_d():
    """Model D: step responses at x = 0, 2, 4, 8, 16, 32 µm, as a fraction of the stationary rise.

    The first three output times are a solver start-up glitch and are dropped. The undercut
    run stops at 98-99% of its final value, so both stacks are normalised to the stationary
    kernel (Model C), not to their own last point.
    """
    names, data = read_comsol(DATA_DIR / "transient.csv")
    par = [column_parameters(n) for n in names[2:]]
    t = np.array([p["t"] for p in par])
    uc = np.array([p["uc"] for p in par])
    x = data[:, 0]
    out = {}
    for flag, name, kernel in ((0, "standard", mf.standard), (1, "undercut", mf.undercut)):
        sel = uc == flag
        out[name] = (t[sel][3:], data[:, 2:][:, sel][:, 3:] / kernel(x, 0.0)[:, None])
    return x, out


@functools.lru_cache(maxsize=None)
def r4_designs():
    """The three membrane designs stepped from 0 to 10° (R4 step responses, plain step).

    Layouts are the best one-way layouts on each R4 run's heater grid, as stored in
    results/R4_settling.csv by run_R4.py. The phase is linear in the drive, so
    φ_n(t) = 814·Σ_k P_k·S_k(x_n, t) with P the steady-state fit for 10°.
    """
    table = pd.read_csv(RESULTS / "R4_settling.csv")
    files = {"undercut": ("R4a", StepResponses("transient_R4a.csv", mirror=True)),
             "slot": ("R4b", StepResponses("transient_R4b.csv")),
             "vacuum": ("R4c", StepResponses("transient_R4vac.csv", default_x_h=34.0))}
    out = {}
    for name, (tag, sr) in files.items():
        rows = table[table.case.str.startswith(tag)]
        row = rows[rows.K == rows.K.max()].iloc[0]
        heaters = [float(v) for v in str(row.heaters_um).split(",")]
        combo = [int(np.flatnonzero(np.isclose(sr.x_h, h))[0]) for h in heaters]
        fit = steer_one_direction(PHASE_COEFF * sr.final[:, combo], xe, +1)
        P = dphi(10) * fit["powers_per_rad"]                               # W per heater
        phase = PHASE_COEFF * np.einsum("tkn,k->tn", sr.S[:, combo], P)    # rad, (times, cores)
        s = strehl(phase - dphi(10) * xe / PITCH_UM)                       # Strehl toward 10°
        # A step of size δ leaves the same error pattern scaled by δ/10° (linearity).
        err = phase[-1] - phase
        slew = {d: settle_time(sr.times, strehl(err * dphi(d) / dphi(10)), STREHL_MIN) * 1e3
                for d in (1, 2, 5, 10)}
        out[name] = dict(heaters=heaters, times_ms=sr.times * 1e3, strehl=s, hold_mW=P.sum() * 1e3,
                         settle_ms=settle_time(sr.times, s, STREHL_MIN) * 1e3, slew_ms=slew,
                         range_deg=fit["range_deg"])
    return out


def rep_transients():
    """Speed with plain step drives: one heater's rise, the beam settling, and step size."""
    x, md = model_d()
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(ps.FULL, 2.75),
                                  gridspec_kw=dict(width_ratios=[1.25, 1.25, 0.85]))
    styles = {0: "-", 8: (0, (4, 2)), 32: (0, (1, 1.5))}
    for name in ("standard", "undercut"):
        t, s = md[name]
        for core, ls in styles.items():
            row = s[int(np.argmin(np.abs(x - core)))]
            a.plot(t * 1e3, row, color=ps.COLOR[name], ls=ls, lw=1.4)
        s0 = s[int(np.argmin(np.abs(x)))]
        t63, t90 = crossing(t, s0, 0.63), crossing(t, s0, 0.90)
        print(f"   Model D {name}: heated core 63% / 90% at {t63 * 1e6:.3g} / {t90 * 1e6:.3g} µs; "
              f"core at 32 µm 63% / 90% at {crossing(t, s[-1], 0.63) * 1e6:.3g} / "
              f"{crossing(t, s[-1], 0.90) * 1e6:.3g} µs; last value / stationary {s0[-1]:.4f}")
    for level in (0.63, 0.90):
        a.axhline(level, color=ps.INK2, lw=0.7, ls=DASH)
        a.text(2.6e-5, level + 0.015, f"{level:.0%}", fontsize=6.5, color=ps.INK2, ha="left")
    a.text(1.6e-3, 0.42, "standard\nstack", fontsize=7, color=ps.INK2, ha="right", va="center")
    a.text(4.2, 0.36, "undercut", fontsize=7, color=ps.INK2, ha="right")
    for core, ls in styles.items():
        a.plot([], [], color=ps.INK2, ls=ls, lw=1.2, label=f"core {core} µm from the heater")
    a.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.3), fontsize=6.6, ncol=1)
    a.set_xscale("log")
    a.set_xlim(2e-5, 5)
    a.set_ylim(0, 1.08)
    a.set_xticks([1e-4, 1e-2, 1])
    a.set_xticklabels(["0.1 µs", "10 µs", "1 ms"])
    a.minorticks_off()
    a.set_xlabel("time after the heater switches on")
    a.set_ylabel("rise / stationary rise")
    ps.panel_label(a, "a")

    des = r4_designs()
    for name, d in des.items():
        n_h = len(d["heaters"])
        b.plot(d["times_ms"], d["strehl"], color=ps.COLOR[name], lw=1.6,
               label=f"{LABEL[name]}, {n_h} heater{'s' if n_h > 1 else ''}: {d['settle_ms']:.1f} ms")
        b.plot(d["settle_ms"], STREHL_MIN, "o", color=ps.COLOR[name], ms=4.5, zorder=4)
        c.plot(list(d["slew_ms"]), list(d["slew_ms"].values()), "o-", color=ps.COLOR[name], ms=3.5,
               lw=1.3)
        print(f"   R4 {name}: heaters {d['heaters']} µm, range {d['range_deg']:.1f}°, hold at 10° "
              f"{d['hold_mW']:.2f} mW, settles 0→10° in {d['settle_ms']:.2f} ms; error of a 1/2/5/10° "
              f"step settles in " + " / ".join(f"{v:.2f}" for v in d["slew_ms"].values()) + " ms")
    b.axhline(STREHL_MIN, color=ps.INK2, lw=0.7, ls=DASH)
    b.text(0.012, STREHL_MIN + 0.02, "Strehl 0.8", fontsize=6.5, color=ps.INK2)
    b.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.3), fontsize=6.6, handlelength=1.2,
             title="settles (Strehl ≥ 0.8 from then on) in:", title_fontsize=6.6, alignment="left")
    b.set_xscale("log")
    b.set_xlim(0.01, 32)
    b.set_ylim(0, 1.05)
    plain(b.xaxis)
    b.set_xlabel("time after a 0 → 10° step (ms)")
    b.set_ylabel("Strehl toward 10°")
    ps.panel_label(b, "b")

    c.set_xscale("log")
    c.set_yscale("log")
    c.set_xticks([1, 2, 5, 10])
    c.minorticks_off()
    c.set_xlim(0.8, 12.5)
    c.set_ylim(0.25, 12)
    plain(c.xaxis)
    plain(c.yaxis)
    c.set_xlabel("step size (°)")
    c.set_ylabel("settling time (ms)")
    ps.panel_label(c, "c")
    fig.tight_layout(w_pad=1.2)
    ps.save(fig, "rep_transients")


# ====================================================================== rep_nonlinear
SLOT_HEATERS = [2, 14, 32]
B_BRACKET = (0.0, 0.7e-3, 1.4e-3)     # relative growth of silicon's dn/dT per K (bracket)
HEATER_LENGTH_M = 1e-3


def rep_nonlinear():
    """Model F's heating shortfall, and what it does to the slot design's beam."""
    x, runs, rises = model_f_runs()
    i0 = int(np.argmin(np.abs(x)))                                         # the middle core
    # (a) how far the linear thermal model is off, and the pointwise law that corrects it;
    # (b) the slot design's beam after every heater is refitted to the true (nonlinear) response.
    fig, (a, c) = plt.subplots(1, 2, figsize=(ps.FULL, 2.6))
    b = None

    A = {}
    th = np.linspace(1, 205, 50)
    for uc, name in ((0, "standard"), (1, "undercut")):
        A[name] = fit_heating(x, runs, rises, uc)
        lin = [rises[i0, i] for i, p in enumerate(runs) if int(p["uc"]) == uc and int(p["nl"]) == 0][0]
        pts = [(p["Q_l"] * lin, rises[i0, i] / (p["Q_l"] * lin)) for i, p in enumerate(runs)
               if int(p["uc"]) == uc and int(p["nl"]) == 1 and p["Q_l2"] == 0]
        a.plot(*zip(*pts), "o", color=ps.COLOR[name], ms=4, zorder=3)
        # Kirchhoff law: ΔT = θ - aθ², so the ratio to the linear rise θ is 1 - aθ.
        a.plot(th, true_rise(th, A[name]) / th, color=ps.COLOR[name], lw=1.2,
               label=f"{LABEL[name].split()[0]}: a = {A[name] * 1e4:.1f}×10⁻⁴ /K")
        print(f"   Model F {name}: a = {A[name]:.3e} /K; middle core actual/linear from "
              f"{pts[0][1]:.4f} at {pts[0][0]:.0f} K to {pts[-1][1]:.4f} at {pts[-1][0]:.0f} K (linear rise)")
    a_mem = A["undercut"]                                                   # used for all membranes
    for dT in (30, 80):
        a.axvline(dT, color=ps.INK2, lw=0.7, ls=DASH)
        a.text(dT + 3, 0.931, f"{dT} K:\n−{100 * a_mem * dT:.1f}%", fontsize=6.5, color=ps.INK2, va="bottom")
        print(f"   shortfall at a linear rise of {dT} K (undercut a): {100 * a_mem * dT:.2f}%")
    a.set_xlim(0, 205)
    a.set_ylim(0.93, 1.005)
    a.set_xlabel("linear-model rise at the middle core (K)")
    a.set_ylabel("COMSOL rise / linear model")
    a.legend(loc="upper right", bbox_to_anchor=(1.06, 1.03), fontsize=6.3, handlelength=1.0)
    ps.panel_label(a, "a")

    # The slot design with 1 mm heaters, under three levels of calibration.
    G = phase_matrix(mf.slot, xe, SLOT_HEATERS)
    p_unit = steer_one_direction(G, xe, +1)["powers_per_rad"]
    angles = np.linspace(0.5, 22.5, 45)
    err = np.zeros((3, len(angles)))
    s_lin, s_gain, s_ref = np.zeros_like(err), np.zeros_like(err), np.zeros_like(err)
    s_model, hot = np.zeros(len(angles)), np.zeros(len(angles))
    for j, theta in enumerate(angles):
        ramp = dphi(theta) * xe / PITCH_UM
        P = dphi(theta) * p_unit
        s_model[j] = strehl(G @ P - ramp)[0]
        hot[j] = true_rise((G @ P).max() / (PHASE_COEFF * HEATER_LENGTH_M), a_mem)
        for i, bb in enumerate(B_BRACKET):
            phi = true_phase(G @ P, HEATER_LENGTH_M, a_mem, bb)
            err[i, j] = pointing_deg(phi, xe) - theta
            s_lin[i, j] = strehl(phi - ramp)[0]
            Pg, _ = gain_corrected_drive(G, P, xe, theta, HEATER_LENGTH_M, a_mem, bb)
            s_gain[i, j] = strehl(true_phase(G @ Pg, HEATER_LENGTH_M, a_mem, bb) - ramp)[0]
            Pr = refit_drive(G, ramp, Pg, HEATER_LENGTH_M, a_mem, bb)
            s_ref[i, j] = strehl(true_phase(G @ Pr, HEATER_LENGTH_M, a_mem, bb) - ramp)[0]
    for dT in (30, 80):                       # where the hottest core passes 30 K and 80 K
        t_flag = float(np.interp(dT, hot, angles))
        print(f"   slot, 1 mm heaters: hottest core reaches {dT} K at {t_flag:.1f}°")

    box = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)
    shade = ps.BLUE_RAMP[4]
    c.fill_between(angles, s_ref.min(axis=0), s_ref.max(axis=0), color=shade, alpha=0.45, lw=0,
                   label="with nonlinearity, heaters refitted\n(band: dn/dT growth b = 0 to 1.4×10⁻³ /K)")
    c.plot(angles, s_model, color=ps.INK, lw=1.0, ls=(0, (4, 2)), label="linear model")
    c.legend(loc="lower left", fontsize=6.3, handlelength=1.6)
    c.axhline(STREHL_MIN, color=ps.INK2, lw=0.7, ls=DASH)
    c.set_xlim(0, 22.8)
    c.set_ylim(0.5, 1.02)
    c.set_xlabel("target angle (°)")
    c.set_ylabel("Strehl toward the target")
    ps.panel_label(c, "b")
    for theta in (5, 10, 15, 20):
        j = int(np.argmin(np.abs(angles - theta)))
        print(f"   slot {angles[j]:.1f}°: hottest core {hot[j]:.0f} K, linear model {s_model[j]:.3f}; "
              f"b = 0/0.7/1.4e-3: pointing error " + "/".join(f"{v:+.2f}" for v in err[:, j])
              + "°; Strehl linear fit " + "/".join(f"{v:.3f}" for v in s_lin[:, j])
              + "; gain " + "/".join(f"{v:.3f}" for v in s_gain[:, j])
              + "; refit " + "/".join(f"{v:.3f}" for v in s_ref[:, j]))
    fig.tight_layout(w_pad=1.3)
    ps.save(fig, "rep_nonlinear")


# ====================================================================== rep_ruze
def rep_ruze():
    """Independent per-waveguide errors: the heaters leave the beam on the Ruze curve."""
    f = pd.read_csv(RESULTS / "fab_white.csv")
    f = f[f.angle_deg == 10]
    budget = pd.read_csv(RESULTS / "fab_budget.csv")
    # Left: the path-length budget without trim; right: what the heaters can do about the errors.
    fig, (b, a) = plt.subplots(1, 2, figsize=(ps.FULL, 2.8), gridspec_kw=dict(width_ratios=[1, 1.35]))

    sig = np.linspace(0, 3.05, 200)
    a.plot(sig, ruze(sig, N_EMITTERS), color=ps.INK, lw=1.0, ls=(0, (4, 2)), zorder=5,
           label="Ruze: $e^{-σ^2} + (1 - e^{-σ^2})/N$")
    for name, key in (("vacuum, 1 heater", "vacuum"), ("slot, 3 heaters", "slot"),
                      ("dual-use, 5 heaters", "dual")):
        d = f[f.design == name].sort_values("sigma_rad")
        a.plot(d.sigma_rad, d.strehl_uncorrected, "o", mfc="white", mec=ps.COLOR[key], mew=1.1, ms=4.5)
        a.plot(d.sigma_rad, d.strehl_corrected, "o-", color=ps.COLOR[key], ms=3, lw=1.0, label=name)
        a.plot([0, 0.06], [d.strehl_trimmed.iloc[0]] * 2, color=ps.COLOR[key], lw=2.2)
        gain = (d.strehl_corrected - d.strehl_uncorrected).max()
        print(f"   {name} at 10°: trimmed {d.strehl_trimmed.iloc[0]:.3f}; largest gain from correction "
              f"{gain:.3f}; corrected − Ruze between {(d.strehl_corrected - d.ruze).min():+.3f} and "
              f"{(d.strehl_corrected - d.ruze).max():+.3f}; (K+1)/N = "
              f"{(int(name.split(', ')[1][0]) + 1) / N_EMITTERS:.3f}")
        for s0 in (0.5, 1.0):
            r = d[np.isclose(d.sigma_rad, s0)].iloc[0]
            print(f"      σ = {s0}: uncorrected {r.strehl_uncorrected:.3f}, corrected "
                  f"{r.strehl_corrected:.3f}, Ruze {r.ruze:.3f}")
    a.plot([], [], "o", mfc="white", mec=ps.INK2, mew=1.1, ms=4.5, label="open: nominal drive kept")
    a.plot([], [], "o", color=ps.INK2, ms=3, label="filled: heaters re-fitted")
    a.text(0.1, 1.035, "bars at σ = 0: after a perfect trim", fontsize=6.5, color=ps.INK2)
    a.set_xlim(0, 3.1)
    a.set_ylim(0, 1.1)
    a.set_xlabel("independent phase error per waveguide, σ (rad)")
    a.set_ylabel("mean Strehl at 10°")
    a.legend(loc="upper right", fontsize=6.6)
    ps.panel_label(a, "b")

    # Without trim: errors grow along each waveguide's independent path, σ² = L/L_coh,
    # so the Strehl (Ruze, large N) is exp(-L/L_coh).
    L = np.linspace(0, 2.5, 200)
    for wg, ls in (("strip", "-"), ("rib", (0, (4, 2)))):
        row = budget[(budget.kind == "budget") & (budget.waveguide == wg) & (budget.condition == "S=0.9")]
        l_coh = -row.length_mm.iloc[0] / np.log(0.9)
        b.plot(L, np.exp(-L / l_coh), color=ps.INK2, ls=ls, lw=1.2,
               label=f"{wg} waveguide, $L_{{coh}}$ = {l_coh:.2f} mm")
        h = budget[(budget.kind == "heated") & (budget.waveguide == wg)]
        b.plot(h.length_mm, h.strehl, "o", color=ps.COLOR["slot"], ms=4.5, zorder=4)
        if wg == "strip":
            for _, r in h.iterrows():
                b.annotate(r.condition.replace("dT=", "≤ ").replace("K", " K"), (r.length_mm, r.strehl),
                           xytext=(4, 4), textcoords="offset points", fontsize=6.5, color=ps.INK2)
        print(f"   {wg}: L_coh = {l_coh:.2f} mm; path for Strehl 0.9 = {row.length_mm.iloc[0]:.2f} mm; "
              "heated section alone at 10° (30/60/80 K): "
              + ", ".join(f"{r.length_mm:.2f} mm → {r.strehl:.2f}" for _, r in h.iterrows()))
    b.plot([], [], "o", color=ps.COLOR["slot"], ms=4.5,
           label="slot design's heated section at 10°,\nfor a hottest-core limit")
    b.axhline(0.9, color=ps.INK2, lw=0.7, ls=DASH)
    b.set_xlim(0, 2.5)
    b.set_ylim(0, 1.03)
    b.set_xlabel("independent path per waveguide (mm)")
    b.set_ylabel("Strehl without trim")
    b.legend(loc="lower left", fontsize=6.4)
    ps.panel_label(b, "a")
    fig.tight_layout(w_pad=1.6)
    ps.save(fig, "rep_ruze")


# ====================================================================== rep_range_vs_K
RANGE_CACHE = RESULTS / "report_range_vs_K.csv"
K_EXHAUSTIVE = 3                         # every layout on the 2 µm grid up to this many heaters
K_MAX = {"standard": 33, "undercut": 8, "far": 8, "slot": 8, "vacuum": 8}


def one_way_layouts():
    """Best one-way (+θ) layouts against the number of heaters, for each kernel.

    Ranked by range at Strehl ≥ 0.8 (to 0.01°), then by power per degree. Up to K_EXHAUSTIVE
    heaters every layout on a 2 µm grid is tried; beyond that one heater at a time is added
    to the previous best layout (greedy), which is a lower bound on what a full search finds.
    Cached in results/report_range_vs_K.csv (delete the file to recompute, about a minute).
    """
    if RANGE_CACHE.exists():
        return pd.read_csv(RANGE_CACHE)
    kernels = {"standard": (mf.standard, np.arange(-32, 33, 2.0)),
               "undercut": (mf.undercut, np.arange(-32, 33, 2.0)),     # heaters within ±32 µm
               "far": (mf.far, np.arange(-32, 41, 2.0)),               # heaters out to the wall
               "slot": (mf.slot, np.arange(-32, 35, 2.0)),
               "vacuum": (mf.vacuum, np.arange(-32, 35, 2.0))}
    rows = []
    for name, (kernel, grid) in kernels.items():
        G_all = phase_matrix(kernel, xe, grid)

        def score(combo):
            r = steer_one_direction(G_all[:, list(combo)], xe, +1)
            return (round(r["range_deg"], 2), -r["mw_per_deg"]), r

        best = None
        for K in range(1, K_MAX[name] + 1):
            if K <= K_EXHAUSTIVE:
                trials = itertools.combinations(range(len(grid)), K)
            else:
                trials = [tuple(sorted(best[0] + (j,))) for j in range(len(grid)) if j not in best[0]]
            top = None
            for combo in trials:
                key, r = score(combo)
                if top is None or key > top[1]:
                    top = (combo, key, r)
            best = top
            rows.append(dict(kernel=name, K=K, search="exhaustive" if K <= K_EXHAUSTIVE else "greedy",
                             heaters_um=" ".join(f"{grid[j]:g}" for j in best[0]),
                             range_deg=best[2]["range_deg"], mw_per_deg=best[2]["mw_per_deg"]))
        print(f"   searched {name}")
    table = pd.DataFrame(rows)
    table.to_csv(RANGE_CACHE, index=False)
    return table


def rep_range_vs_K():
    """Design map: range and power against the number of heaters, and against heated-section pitch."""
    table = one_way_layouts()
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(ps.FULL, 2.9),
                                  gridspec_kw=dict(width_ratios=[1.2, 1.2, 1]))
    series = [("standard", "standard", "-", "standard stack"),
              ("undercut", "undercut", "-", "undercut, heaters within ±32 µm"),
              ("far", "undercut", (0, (3, 2)), "undercut, heaters out to the walls"),
              ("slot", "slot", "-", "slot in air"), ("vacuum", "vacuum", "-", "vacuum")]
    for key, colour, ls, label in series:
        d = table[table.kernel == key]
        a.plot(d.K, d.range_deg, "o", ls=ls, color=ps.COLOR[colour], ms=3, lw=1.3, label=label)
        b.plot(d.K, d.mw_per_deg, "o", ls=ls, color=ps.COLOR[colour], ms=3, lw=1.3)
        full = d[d.range_deg > FOV_DEG - 0.01]
        print(f"   one-way {key}: " + "; ".join(
            f"K={r.K}: {r.range_deg:.1f}° at {r.mw_per_deg:.2f} mW/deg [{r.heaters_um}]"
            for r in d[d.K <= 5].itertuples())
            + (f"; full field first at K = {full.K.min()}" if len(full) else
               f"; best {d.range_deg.max():.1f}° at K = {d.K.iloc[int(d.range_deg.argmax())]}"))
    # Two-sided layouts (analysis 2, ranked by +range plus −range): half the span is the ± range.
    pl = pd.read_csv(RESULTS / "placement.csv")
    u = pl[pl["stack"] == "undercut"]
    a.plot(u.K, u.span_deg / 2, "D", mfc="white", mec=ps.COLOR["undercut"], mew=1.1, ms=4, ls="none",
           label="undercut, two-sided: ± range (±32 µm)")
    b.plot(u.K, u.mw_per_deg, "D", mfc="white", mec=ps.COLOR["undercut"], mew=1.1, ms=4, ls="none")
    far4 = evaluate_layout(mf.far, xe, [-32.75, -6.7, 9.4, 32.6])       # layout from run_modelG
    a.plot(4, far4["span_deg"] / 2, "*", mfc="white", mec=ps.COLOR["undercut"], mew=1.1, ms=9,
           label="same, heaters out to ±32.7 µm")
    b.plot(4, far4["mw_per_deg"], "*", mfc="white", mec=ps.COLOR["undercut"], mew=1.1, ms=9)
    print(f"   two-sided undercut (placement.csv): ± range " + ", ".join(
        f"K={r.K}: {r.span_deg / 2:.1f}°" for r in u.itertuples())
        + f"; far kernel K=4 at (−32.75, −6.7, 9.4, 32.6): ±{far4['span_deg'] / 2:.1f}° at "
          f"{far4['mw_per_deg']:.2f} mW/deg")
    a.axhline(FOV_DEG, color=ps.INK2, lw=0.7, ls=DASH)
    a.text(0.9, FOV_DEG + 0.5, f"{FOV_DEG:.1f}°: single-beam field of view", fontsize=6.5, color=ps.INK2,
           ha="left", va="bottom")
    for ax in (a, b):
        ax.set_xscale("log")
        ax.set_xticks([1, 2, 3, 4, 6, 8, 16, 33])
        ax.minorticks_off()
        plain(ax.xaxis)
        ax.set_xlim(0.85, 40)
        ax.set_xlabel("number of heaters, K")
    a.set_ylim(0, 26.5)
    a.set_ylabel("range toward +θ at Strehl ≥ 0.8 (°)")
    b.set_yscale("log")
    plain(b.yaxis)
    b.set_ylabel("steering power (mW/deg)")
    ps.panel_label(a, "a")
    ps.panel_label(b, "b")

    # Standard stack against the heated-section pitch d_h. The heaters sit before the gratings,
    # so d_h can differ from the 2 µm emitter pitch: the phase step per emitter sets the angle
    # (through the emitter pitch), while d_h sets how many pitches the heat spreads over.
    l_th = half_width(mf.standard.x, mf.standard.k)
    pitches = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    for K, shade in zip((4, 8, 16, 33), ps.BLUE_RAMP[1:]):
        rng = []
        for d_h in pitches:
            x_heated = (np.arange(N_EMITTERS) - (N_EMITTERS - 1) / 2) * d_h
            heaters = np.linspace(x_heated[0], x_heated[-1], K)
            r = steer_one_direction(phase_matrix(mf.standard, x_heated, heaters), x_heated, +1, pitch=d_h)
            rng.append(np.degrees(np.arcsin(r["dphi_max"] * U_PER_RAD)))   # angle at the 2 µm emitters
        c.plot(pitches, rng, "o-", color=ps.MUTED if K == 33 else shade, ms=3, lw=1.3)
        c.annotate(f"K = {K}", (pitches[-1], rng[-1]), xytext=(4, 0), textcoords="offset points",
                   fontsize=6.6, color=ps.INK2, va="center")
        print(f"   standard stack, {K} evenly spaced heaters, d_h = {pitches} µm: range "
              + ", ".join(f"{v:.2f}" for v in rng) + f"° (L_th = {l_th:.2f} µm, L_th/d_h = "
              + ", ".join(f"{l_th / p:.2f}" for p in pitches) + ")")
    c.set_yscale("log")
    plain(c.yaxis)
    c.set_xlim(1.4, 4.7)
    c.set_xlabel("heated-section pitch (µm)")
    c.set_ylabel("standard stack: range (°)")
    top = c.secondary_xaxis("top", functions=(lambda d: l_th / np.maximum(d, 1e-9),
                                              lambda r: l_th / np.maximum(r, 1e-9)))
    top.set_xticks([2.2, 1.6, 1.1, 0.84])
    top.set_xlabel("spreading length / pitch", fontsize=7.5)
    top.tick_params(labelsize=7)
    ps.panel_label(c, "c", x=-0.08)
    fig.legend(*a.get_legend_handles_labels(), loc="lower center", ncol=3, fontsize=6.6,
               bbox_to_anchor=(0.5, -0.13))
    fig.tight_layout(w_pad=1.2)
    ps.save(fig, "rep_range_vs_K")


# ====================================================================== rep_power_speed
def rep_power_speed():
    """The design path on hold power against settling time, beside published heaters (plain steps)."""
    des = r4_designs()
    steps = [("standard", "standard stack, 33 heaters", mf.standard_settling_ms(),
              ideal_power(xe, dphi(10), mf.P_PI_STANDARD) * 1e3)]
    for name in ("undercut", "slot", "vacuum"):
        n = len(des[name]["heaters"])
        steps.append((name, f"{LABEL[name]}, {n} heater{'s' if n > 1 else ''}", des[name]["settle_ms"],
                      des[name]["hold_mW"]))
    fig, ax = plt.subplots(figsize=(ps.FULL, 3.6))
    for prod, (lx, ly) in ((1, (1.3e-3, 520)), (10, (0.013, 1300)), (100, (0.13, 1300))):
        t = np.array([1e-3, 30])                           # lines of equal power × time
        ax.plot(t, prod / t, color=ps.AXIS, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.text(lx, ly, f"{prod} mW·ms", fontsize=6.5, color=ps.MUTED, ha="left", va="center")
    for label, t, p, off in mf.PUBLISHED:
        ax.plot(t, p, "o", color=ps.MUTED, ms=4.5, zorder=3)
        ax.annotate(label, (t, p), xytext=off, textcoords="offset points", fontsize=6.5,
                    color=ps.INK2, ha="left" if off[0] > 0 else "right", va="center")
    for (_, _, t0, p0), (_, _, t1, p1) in zip(steps[:-1], steps[1:]):
        ax.annotate("", xy=(t1, p1), xytext=(t0, p0),
                    arrowprops=dict(arrowstyle="-|>", color=ps.INK2, lw=1.0, shrinkA=6, shrinkB=6,
                                    mutation_scale=8), zorder=2)
    for name, label, t, p in steps:
        ax.plot(t, p, "o", color=ps.COLOR[name], ms=7, zorder=4)
        ax.annotate(f"{label}: {p:.3g} mW", (t, p), xytext=(7, 0), textcoords="offset points",
                    fontsize=7, color=ps.INK, va="center", fontweight="bold")
        print(f"   {label}: hold {p:.2f} mW at 10°, settles in {t:.3g} ms, power × time {p * t:.1f} mW·ms")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 30)
    ax.set_ylim(0.8, 2000)
    plain(ax.xaxis)
    plain(ax.yaxis)
    ax.set_xlabel("time to settle a 0 → 10° step (ms)")
    ax.set_ylabel("power to hold 10° (mW)")
    ps.save(fig, "rep_power_speed")


# ====================================================================== rep_membrane_fit
L_HEAL_UM = 63.0                 # healing length fitted to the undercut kernel (beyond_steering.md)
WALL_UM = 40.0                   # cavity walls at ±40 µm
KT_OXIDE = 1.38 * 4.22e-6        # W/K: k_ox × membrane thickness, nothing fitted


def two_wall_membrane(x_um, x_h_um, kt=membrane1d.KT, heal_um=L_HEAL_UM):
    """Rise (K per W/m) on a membrane held at ambient at both walls, with loss to the air below.

    k·t·T'' - (k·t/L²)·T = -δ(x - x_h), T(±40 µm) = 0. Either side of the source the solution
    is a sinh that vanishes at its wall; matching the heat flux at the source gives
        T = L·sinh((x_< + W)/L)·sinh((W - x_>)/L) / (k·t·sinh(2W/L)).
    """
    L = heal_um * 1e-6
    lo = (np.minimum(x_um, x_h_um) + WALL_UM) * 1e-6
    hi = (WALL_UM - np.maximum(x_um, x_h_um)) * 1e-6
    T = L * np.sinh(lo / L) * np.sinh(hi / L) / (kt * np.sinh(2 * WALL_UM * 1e-6 / L))
    return np.where(np.abs(x_um) <= WALL_UM, T, 0.0)


def rep_membrane_fit():
    """The 1D membrane model against the COMSOL kernels."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(ps.FULL, 2.8))
    x = np.linspace(-48, 48, 961)
    peak = float(mf.undercut(0.0, 0.0))
    for xh, shade in zip((0, 16, 32), ("#0d7a52", ps.COLOR["undercut"], "#8fdcc0")):
        a.plot(x, mf.undercut(x, xh), color=shade, lw=1.6)
        a.plot(x, two_wall_membrane(x, xh), color=shade, lw=1.0, ls=(0, (4, 2)))
        err = np.max(np.abs(two_wall_membrane(xe, xh) - mf.undercut(xe, xh))) / peak
        a.annotate(f"heater at {xh} µm", (xh, float(mf.undercut(xh, xh))),
                   xytext={0: (0, 7), 16: (7, 3), 32: (4, 7)}[xh], textcoords="offset points",
                   ha={0: "center", 16: "left", 32: "right"}[xh], fontsize=6.5, color=ps.INK2)
        print(f"   undercut, heater at {xh} µm: worst |1D − COMSOL| over the cores = {100 * err:.1f}% of K(0;0)")
    for w in (-WALL_UM, WALL_UM):
        a.axvline(w, color=ps.INK2, lw=0.7, ls=DASH)
    a.text(-WALL_UM - 1, 0.15, "cavity wall", rotation=90, fontsize=6.5, color=ps.INK2, ha="right")
    a.plot([], [], color=ps.INK2, lw=1.6, label="COMSOL (2D)")
    a.plot([], [], color=ps.INK2, lw=1.0, ls=(0, (4, 2)),
           label=f"1D membrane, L = {L_HEAL_UM:g} µm,\nk·t = {membrane1d.KT * 1e6:.1f} µW/K")
    a.legend(loc="upper left", fontsize=6.5, handlelength=1.8)
    a.set_xlim(-48, 48)
    a.set_ylim(0, 4.4)
    a.set_xlabel("position across the array (µm)")
    a.set_ylabel("rise (K per W/m of heater power)")
    ps.panel_label(a, "a")

    # Slot membrane, heater at 34 µm. Air at 1 atm, and COMSOL's "vacuum" (k = 1e-4 W/m/K).
    xs = np.linspace(membrane1d.X0, membrane1d.X1, 400)
    D, w = 15e-6, 2e-6
    for name, kernel, k_gas in (("slot", mf.slot, membrane1d.k_gas(D, 101325.0)),
                                ("vacuum", mf.vacuum, 1e-4)):
        k_slot = membrane1d.k_gas(w, 101325.0) if name == "slot" else 1e-4
        model = membrane1d.profiles(membrane1d.A_H * k_gas / D,
                                    membrane1d.A_S * k_slot * membrane1d.T_MEM / w, [34.0], x_emit=xs)[:, 0]
        b.plot(x, kernel(x, 34.0), color=ps.COLOR[name], lw=1.6)
        b.plot(xs, model, color=ps.COLOR[name], lw=1.0, ls=(0, (4, 2)))
        at_cores = np.interp(xe, xs, model)
        comsol = kernel(xe, 34.0)
        fit_c = steer_one_direction(phase_matrix(kernel, xe, [34.0]), xe, +1)
        fit_m = steer_one_direction(PHASE_COEFF * at_cores[:, None], xe, +1)
        line = np.polyval(np.polyfit(xe, comsol, 1), xe)
        print(f"   {name}, heater at 34 µm: COMSOL {fit_c['mw_per_deg']:.3f} mW/deg (range "
              f"{fit_c['range_deg']:.1f}°), 1D model {fit_m['mw_per_deg']:.3f} mW/deg; worst |1D − COMSOL| "
              f"over the cores {100 * np.max(np.abs(at_cores - comsol)) / comsol.max():.1f}% of the peak; "
              f"centre core {float(kernel(0.0, 34.0)):.2f} K per W/m; COMSOL rms off a straight line "
              f"{100 * np.sqrt(np.mean((comsol - line) ** 2)) / (line[-1] - line[0]):.3f}% of the rise")
        if name == "vacuum":
            vac_mw = fit_c["mw_per_deg"]
    # Conduction limit: all the heat crosses to the far wall, so dT/dx = q/(k·t), a straight line.
    xl = np.array([-WALL_UM, 34.0])
    b.plot(xl, (xl + WALL_UM) * 1e-6 / KT_OXIDE, color=ps.INK, lw=0.9, ls=(0, (1, 1.5)))
    limit = KT_OXIDE / (PHASE_COEFF * PITCH_UM * 1e-6) / U_PER_RAD * np.pi / 180 * 1e3   # mW/deg
    print(f"   conduction limit k·t/(814·d) = {limit:.3f} mW/deg; COMSOL vacuum {vac_mw:.3f} "
          f"({100 * (vac_mw / limit - 1):.1f}% above)")
    b.text(-46, 15.2, f"vacuum, one heater: {vac_mw:.2f} mW/deg (COMSOL)\n"
           f"conduction limit k·t/(814·d): {limit:.3f} mW/deg", fontsize=6.6, color=ps.INK2, va="top")
    b.text(12, 6.9, "vacuum", fontsize=7, color=ps.INK2, ha="left")
    b.text(20, 2.6, "slot in air", fontsize=7, color=ps.INK2)
    b.plot([], [], color=ps.INK2, lw=1.6, label="COMSOL (2D)")
    b.plot([], [], color=ps.INK2, lw=1.0, ls=(0, (4, 2)), label="1D membrane with gas")
    b.plot([], [], color=ps.INK, lw=0.9, ls=(0, (1, 1.5)), label="straight line, slope q/(k·t)")
    b.legend(loc="upper left", bbox_to_anchor=(0.0, 0.84), fontsize=6.5, handlelength=1.8)
    for w_ in (-WALL_UM, 36.0):
        b.axvline(w_, color=ps.INK2, lw=0.7, ls=DASH)
    b.text(37, 6.0, "slot", rotation=90, fontsize=6.5, color=ps.INK2)
    b.set_xlim(-48, 48)
    b.set_ylim(0, 15.5)
    b.set_xlabel("position across the array (µm)")
    b.set_ylabel("rise (K per W/m), heater at 34 µm")
    ps.panel_label(b, "b")
    fig.tight_layout(w_pad=1.8)
    ps.save(fig, "rep_membrane_fit")


FIGURES = {f.__name__: f for f in (rep_transients, rep_nonlinear, rep_ruze, rep_range_vs_K,
                                   rep_power_speed, rep_membrane_fit)}

if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGURES):
        print(name)
        FIGURES[name]()
