"""Analysis 5: effective degrees of freedom. How many phase patterns can heat make?

A heater under every emitter (33 heaters), offset projected out. Parts:
  1. checks: on a ring the matrix is circulant and its singular values must
     equal the Fourier transform exactly; on an open array they approach it
     (Szegő) when the array is much wider than the kernel;
  2. usable modes against L_th/d: the standard stack with the pitch varied
     from 1.5 to 4 µm (it is shift-invariant, so any pitch can be sampled),
     plus the undercut, far, slot and vacuum kernels at 2 µm. For the standard
     stack, a Gaussian kernel of the same half-width predicts about
     0.80·N·d/L_th usable modes at the 10x threshold (√(2 ln 10)·√(2 ln 2)/π);
  3. what the cheapest modes look like, and what they cost;
  4. a sharp pattern: splitting the beam in two at ±3°, ±5° and ±10°, fitted
     with all 33 heaters on each stack. Two-beam efficiency against the ideal,
     and the heater power.
Writes results/modes.csv, modes_shapes.csv and modes_two_beam.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np

from opa.constants import PITCH_UM
from opa.farfield import array_factor
from opa.modes import half_width, modes, toeplitz_spectrum, two_beam_targets, usable
from opa.steering import fit_powers
from opa.thermal import (StandardKernel, SweptKernel, UndercutKernel, emitter_positions,
                         phase_matrix)

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
N = 33
TAUS = (0.1, 0.01)

# ---------------------------------------------------------------- 1. check
print("1. Checks")
# Exact: on a ring (periodic distances) the matrix is circulant, and its singular values ARE the |DFT|.
xe = emitter_positions(N, 2.0)
dist = np.abs(xe[:, None] - xe[None, :])
ring = np.exp(-0.5 * (np.minimum(dist, N * 2.0 - dist) / 4.0) ** 2)
s_ring = np.sort(np.linalg.svd(ring, compute_uv=False))[::-1]
dft = np.sort(np.abs(np.fft.fft(ring[0])))[::-1]
print(f"  circulant (ring) Gaussian, width 4 µm: largest |singular value - |DFT|| = {np.max(np.abs(s_ring - dft)):.1e}"
      f" (relative to {s_ring[0]:.2f})")
print("  open array (Toeplitz): singular values against the Fourier transform; Szegő holds for arrays much wider than the kernel")
for w in (2.0, 4.0, 8.0):
    xe = emitter_positions(N, 2.0)
    G = np.exp(-0.5 * ((xe[:, None] - xe[None, :]) / w) ** 2)
    s = np.linalg.svd(G, compute_uv=False)
    f = np.sort(toeplitz_spectrum(lambda dx: np.exp(-0.5 * (dx / w) ** 2), 2.0, N))[::-1]
    worst = np.max(np.abs(s[:8] - f[:8]) / f[0])
    print(f"  width {w:g} µm: first 5 singular values {np.round(s[:5] / s[0], 3)};  Fourier {np.round(f[:5] / f[0], 3)};"
          f"  worst gap over the first 8, as a fraction of the largest: {worst:.3f}")

# ---------------------------------------------------------------- 2. usable modes
print("\n2. Usable modes (offset removed) against L_th/d, a heater under every emitter")
std = StandardKernel()
L_std = half_width(std.x, std.k)
cases = []
for d in (1.5, 2.0, 2.5, 3.0, 4.0):
    xe = emitter_positions(N, d)
    cases.append((f"standard, d = {d:g} µm", phase_matrix(std, xe, xe), L_std, d, xe))
xe = emitter_positions()
under = UndercutKernel()
far = SweptKernel("kernel_undercut_far.csv", "far", True, 41)
slot = SweptKernel("kernel_slot.csv", "slot", False, 67)
vac = SweptKernel("kernel_vacuum.csv", "vacuum", False)
vac_heaters = [h for h in vac.x_h if abs(h) <= 32]              # its grid: every 2 µm, the cores
for name, kern, heaters, prof in (
        ("undercut", under, xe, under(under.x, 0.0)),
        ("undercut, far kernel", far, xe, far(far.x, 0.0)),
        ("slot", slot, xe, slot.column(0.0)),
        ("vacuum", vac, vac_heaters, vac.column(0.0))):
    cases.append((name, phase_matrix(kern, xe, heaters), half_width(kern.x, prof), PITCH_UM, xe))

rows, shapes = [], []
print(f"  {'case':24} {'L_th':>6} {'L_th/d':>7}  usable at 10× / 100×   gain of modes 1-6 relative to mode 1")
for name, G, L, d, x in cases:
    s, U, V = modes(G)
    n10, n100 = usable(s, 0.1), usable(s, 0.01)
    rows.append([name, d, L, L / d, n10, n100, *s[:12]])
    guess = f"  (Gaussian estimate {min(N - 1, 0.80 * N * d / L):4.1f})" if name.startswith("standard") else ""
    print(f"  {name:24} {L:5.1f}  {L / d:6.2f}   {n10:5d} / {n100:3d}            "
          + " ".join(f"{v:.3f}" for v in s[:6] / s[0]) + guess)
    if name in ("standard, d = 2 µm", "undercut", "slot", "vacuum"):
        for i in range(4):
            u = U[:, i] * np.sign(U[np.argmax(np.abs(U[:, i])), i])    # sign convention: largest entry positive
            shapes.append([name, i + 1, s[i], *u])
with open(RESULTS / "modes.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["case", "pitch_um", "L_th_um", "L_th_over_d", "usable_10x", "usable_100x",
                *[f"s{i + 1}_rad_per_W" for i in range(12)]])
    w.writerows(rows)
with open(RESULTS / "modes_shapes.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["case", "mode", "gain_rad_per_W", *[f"x{v:g}" for v in xe]])
    w.writerows(shapes)

# ---------------------------------------------------------------- 3. mode shapes and cost
print("\n3. The cheapest modes: shape and the power for 1 rad rms of that pattern (mW, 2-norm of the heater powers)")
for name, G, L, d, x in cases:
    if name not in ("standard, d = 2 µm", "undercut", "slot", "vacuum"):
        continue
    s, U, V = modes(G)
    descr = []
    for i in range(4):
        u = U[:, i]
        # name the shape by its best-matching Legendre polynomial
        t = np.linspace(-1, 1, len(u))
        leg = [np.polynomial.legendre.Legendre.basis(k)(t) for k in range(1, 6)]
        corr = [abs(np.dot(u, p - p.mean())) / np.linalg.norm(p - p.mean()) for p in leg]
        k = int(np.argmax(corr))
        label = ["ramp", "dome", "S-curve", "W", "5th"][k] if max(corr) > 0.9 else "ripple"
        descr.append(f"{label} {np.sqrt(len(u)) / s[i] * 1e3:.1f}")
    print(f"  {name:22} " + ";  ".join(descr))

# ---------------------------------------------------------------- 4. two beams
print("\n4. A sharp pattern: split the beam into two at ±θ, all 33 heaters, least squares")
print("   two-beam efficiency (light in the two target beams) against the ideal pattern; heater power for the sine")
rows = []
for name, G, L, d, x in cases:
    if name not in ("standard, d = 2 µm", "undercut", "slot", "vacuum"):
        continue
    for theta in (3, 5, 10):
        us = np.sin(np.radians(theta))
        parts = []
        for kind, target in two_beam_targets(x, us).items():
            ideal = array_factor(target, x, np.array([us, -us])).sum()
            P, c = fit_powers(G, target)
            eff = array_factor(G @ P, x, np.array([us, -us])).sum()
            rows.append([name, theta, kind, ideal, eff, P.sum() * 1e3])
            power = f", {P.sum() * 1e3:4.0f} mW" if kind == "sine" else ""   # binary: least squares chases the edges
            parts.append(f"{kind}: {eff:.2f} of ideal {ideal:.2f}{power}")
        print(f"  {name:22} ±{theta:2d}°  " + ";   ".join(parts))
with open(RESULTS / "modes_two_beam.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["case", "theta_deg", "target", "ideal_efficiency", "fitted_efficiency", "power_mW"])
    w.writerows(rows)
