"""Analysis 6: fabrication errors. How much can a few smooth heaters correct?

Designs: the final layouts (one-way steering toward +θ):
  vacuum, 1 heater at 34 µm         kernel_vacuum.csv
  slot, 3 heaters at 2, 14, 32 µm   kernel_slot.csv
  undercut, 3 heaters at -8, 10, 32 today's undercut kernel
  dual-use, 5 heaters at 0, ±16, ±33 kernel_undercut_far.csv (steers and focuses)
plus evenly spaced slot layouts with 5, 9 and 17 heaters, to see how many
heaters it would take.

Parts:
  1. checks against Ruze and the least-squares projection limit;
  2. independent per-waveguide errors, σ = 0.1-3 rad: uncorrected, corrected
     by the heaters, trimmed, and Ruze (the brief's analysis 6);
  3. correlated errors: Strehl against correlation length ℓ at σ = 1 rad;
  4. the cost of correcting: extra heater power, ours against conventional;
  5. the path-length budget: how long each waveguide's independent path may be
     without trim, against the heated length the temperature limit asks for.

Strehl is toward the target angle, offset free (as in steering.py), averaged
over Monte Carlo trials. Writes results/fab_white.csv, fab_correlation.csv,
fab_power.csv and fab_budget.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np

from opa.constants import PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.fabrication import (corrected, expected_strehl, phase_errors, ruze, strehl, uncorrected,
                             wrapped_power)
from opa.steering import fit_powers, steer_one_direction
from opa.thermal import SweptKernel, UndercutKernel, emitter_positions, phase_matrix

sys.stdout.reconfigure(encoding="utf-8")
RESULTS = Path(__file__).parent / "results"
xe = emitter_positions()
N = len(xe)
rng = np.random.default_rng(798)
L_COH_MM = {"strip": 4.17, "rib": 1.61}   # Yang 2015 [V]


def dphi(theta_deg):
    """Ramp step per emitter (rad) for steering angle θ."""
    return 2 * np.pi * PITCH_UM * np.sin(np.radians(theta_deg)) / WAVELENGTH_UM


def ramp(theta_deg):
    return dphi(theta_deg) * xe / PITCH_UM


slot = SweptKernel("kernel_slot.csv", "slot", False, 67)
vacuum = SweptKernel("kernel_vacuum.csv", "vacuum", False)
far = SweptKernel("kernel_undercut_far.csv", "far", True, 41)
today = UndercutKernel()


class Design:
    def __init__(self, name, kernel, heaters):
        self.name, self.heaters = name, list(heaters)
        self.G = phase_matrix(kernel, xe, heaters)
        self.p_unit = steer_one_direction(self.G, xe, +1)["powers_per_rad"]

    def drive(self, theta):
        return dphi(theta) * self.p_unit


designs = [
    Design("vacuum, 1 heater", vacuum, [34]),
    Design("slot, 3 heaters", slot, [2, 14, 32]),
    Design("undercut, 3 heaters", today, [-8, 10, 32]),
    Design("dual-use, 5 heaters", far, [-33, -16, 0, 16, 33]),
]
even = [Design(f"slot, {k} even", slot, list(np.round(np.linspace(-32, 32, k)))) for k in (5, 9, 17)]

# ---------------------------------------------------------------- 1. checks
print("1. Checks")
for sigma in (0.1, 0.5, 1.0, 2.0, 3.0):
    e = phase_errors(xe, sigma, 0, 20000, rng)
    mc = strehl(e).mean()
    print(f"  Ruze, σ = {sigma:3.1f} rad: Monte Carlo {mc:.4f}, analytic {ruze(sigma, N):.4f}")

print("  Exact Gaussian expectation at 10° (paired ±ε samples cancel the linear noise term):")
for sigma in (0.1, 0.3):
    for d in designs:
        e = phase_errors(xe, sigma, 0, 2000, rng)
        e = np.vstack([e, -e])
        A = np.hstack([d.G, np.ones((N, 1))])
        Q = np.eye(N) - A @ np.linalg.pinv(A)
        r0 = d.G @ d.drive(10) - ramp(10)
        s_corr, _ = corrected(d.G, ramp(10), e)
        print(f"    σ = {sigma}  {d.name:22} uncorrected {uncorrected(d.G, d.drive(10), ramp(10), e).mean():.4f} "
              f"(exact {expected_strehl(r0, np.eye(N), sigma):.4f});  corrected {s_corr.mean():.4f} "
              f"(exact {expected_strehl(-Q @ ramp(10), Q, sigma):.4f})")
e = phase_errors(xe, 1.0, 1e5, 2000, rng)
d = designs[1]
print(f"  ℓ = 100 mm (a pure offset), σ = 1 rad, slot at 10°: uncorrected "
      f"{uncorrected(d.G, d.drive(10), ramp(10), e).mean():.4f}, trimmed "
      f"{strehl(d.G @ d.drive(10) - ramp(10))[0]:.4f} (should match: an offset costs nothing)")

# ---------------------------------------------------------------- 2. independent errors
print("\n2. Independent per-waveguide errors: mean Strehl (uncorrected / corrected by heaters / trimmed)")
SIGMAS = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
ANGLES = (0, 5, 10)
TRIALS = 400
rows = []
for d in designs:
    for theta in ANGLES:
        s_trim = strehl(d.G @ d.drive(theta) - ramp(theta))[0]
        line = []
        for sigma in SIGMAS:
            e = phase_errors(xe, sigma, 0, TRIALS, rng)
            s_unc = uncorrected(d.G, d.drive(theta), ramp(theta), e).mean()
            s_corr, _ = corrected(d.G, ramp(theta), e)
            rows.append([d.name, theta, sigma, s_unc, s_corr.mean(), s_trim, ruze(sigma, N)])
            line.append(f"{s_unc:.2f}/{s_corr.mean():.2f}")
        print(f"  {d.name:22} {theta:2d}°  trim {s_trim:.3f} | " + "  ".join(
            f"σ={s}: {l}" for s, l in zip(SIGMAS, line) if s in (0.3, 0.5, 1.0, 2.0)))
print("  Ruze:", "  ".join(f"σ={s}: {ruze(s, N):.2f}" for s in (0.3, 0.5, 1.0, 2.0)))
with open(RESULTS / "fab_white.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["design", "angle_deg", "sigma_rad", "strehl_uncorrected", "strehl_corrected",
                "strehl_trimmed", "ruze"])
    w.writerows(rows)

# ---------------------------------------------------------------- 3. correlated errors
print("\n3. Correlated errors, σ = 1 rad: mean Strehl after heater correction, against correlation length ℓ")
ELLS = [0, 2, 4, 8, 16, 32, 64, 128, 1000]
rows = []
print(f"  {'design':22} {'θ':>3}  " + "  ".join(f"ℓ={ell:<4g}" for ell in ELLS))
for d in designs + even:
    for theta in ((0, 10) if d in designs else (10,)):
        s_trim = strehl(d.G @ d.drive(theta) - ramp(theta))[0]   # the design's own steering residual
        vals = []
        for ell in ELLS:
            e = phase_errors(xe, 1.0, ell, TRIALS, rng)
            s_corr, _ = corrected(d.G, ramp(theta), e)
            s_unc = uncorrected(d.G, d.drive(theta), ramp(theta), e).mean()
            rows.append([d.name, theta, ell, s_unc, s_corr.mean(), s_trim])
            vals.append(s_corr.mean())
        print(f"  {d.name:22} {theta:2d}°  " + "  ".join(f"{v:6.2f}" for v in vals) + f"   (trim {s_trim:.2f})")
e_rows = {ell: phase_errors(xe, 1.0, ell, TRIALS, rng) for ell in ELLS}
print(f"  {'uncorrected':22}      " + "  ".join(f"{strehl(e_rows[ell]).mean():6.2f}" for ell in ELLS))
with open(RESULTS / "fab_correlation.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["design", "angle_deg", "ell_um", "strehl_uncorrected", "strehl_corrected", "strehl_trimmed"])
    w.writerows(rows)

# ---------------------------------------------------------------- 4. power cost of correcting
print("\n4. Power cost of correcting errors (mean over trials)")
print("   ours: extra heater power over the nominal drive (mW); conventional: extra power in units of P_π")
rows = []
for sigma, ell in ((0.5, 0), (1.0, 0), (1.0, 64), (1.0, 1000)):
    e = phase_errors(xe, sigma, ell, TRIALS, rng)
    conv = {theta: np.mean([wrapped_power(ramp(theta) - ei, 1.0) - wrapped_power(ramp(theta), 1.0)
                            for ei in e]) for theta in (0, 10)}
    for d in designs:
        for theta in (0, 10):
            s_corr, p = corrected(d.G, ramp(theta), e)
            extra = (p.mean() - d.drive(theta).sum()) * 1e3
            rows.append([d.name, sigma, ell, theta, s_corr.mean(), extra, conv[theta]])
    print(f"  σ = {sigma} rad, ℓ = {ell:g} µm: conventional +{conv[0]:.1f} P_π at 0°, +{conv[10]:.1f} P_π at 10°;  ours: "
          + "; ".join(f"{r[0].split(',')[0]} +{r[5]:.2f} mW @{r[3]}°" for r in rows[-len(designs) * 2:]))
with open(RESULTS / "fab_power.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["design", "sigma_rad", "ell_um", "angle_deg", "strehl_corrected", "extra_mW_ours",
                "extra_conventional_in_P_pi"])
    w.writerows(rows)

# ---------------------------------------------------------------- 5. path-length budget
print("\n5. Path-length budget without trim (independent errors, Ruze)")
rows = []
for target_s in (0.95, 0.9, 0.8):
    for kind, lc in L_COH_MM.items():
        L = lc * np.log(1 / target_s)          # σ² = L/L_coh and S ≈ exp(-σ²)
        rows.append(["budget", f"S={target_s}", kind, L, np.nan])
        print(f"  Strehl {target_s}: independent path ≤ {L:.2f} mm ({kind})")
d = designs[1]
phi_hot = float((d.G @ d.drive(10)).max())    # hottest core phase at 10°, rad
print(f"  The slot design's hottest core carries {phi_hot:.0f} rad at 10°. Heated length for a limit on its rise:")
for dT in (30, 60, 80):
    L_h = phi_hot / (PHASE_COEFF * dT) * 1e3   # mm
    for kind, lc in L_COH_MM.items():
        s = np.exp(-L_h / lc)
        rows.append(["heated", f"dT={dT}K", kind, L_h, s])
    print(f"    ΔT ≤ {dT} K: L_h ≥ {L_h:.2f} mm; that section alone gives Strehl "
          f"{np.exp(-L_h / L_COH_MM['strip']):.2f} (strip), {np.exp(-L_h / L_COH_MM['rib']):.2f} (rib) untrimmed")
with open(RESULTS / "fab_budget.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["kind", "condition", "waveguide", "length_mm", "strehl"])
    w.writerows(rows)
