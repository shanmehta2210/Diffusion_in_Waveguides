"""Check the pipeline before using it: COMSOL numbers, then analytic limits.

Run from the analysis folder:  python check_pipeline.py
Every line prints PASS or FAIL with the value found and the value expected.
"""
import sys

import numpy as np

from opa.comsol_io import read_comsol
from opa.conventional import ideal_power, with_crosstalk
from opa.constants import DATA_DIR, DNEFF_DT, N_EMITTERS, PHASE_COEFF, PITCH_UM, WAVELENGTH_UM
from opa.farfield import U_GRID, beam_metrics, element_pattern, far_field, strehl
from opa.steering import fit_powers, steer_one_direction, strehl_vs_step
from opa.thermal import StandardKernel, UndercutKernel, emitter_positions, phase_matrix

sys.stdout.reconfigure(encoding="utf-8")  # the Windows console defaults to cp1252, which has no λ
failures = 0


def check(label, value, expected, tol):
    global failures
    ok = abs(value - expected) <= tol
    failures += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {label}: {value:.6g} (expected {expected:.6g} ± {tol:.2g})")


print("--- Constants against the COMSOL data")
_, t = read_comsol(DATA_DIR / "neff_vs_T.csv")
check("dn_eff/dT fitted from neff_vs_T [1/K]", np.polyfit(t[:, 0], t[:, 1], 1)[0], DNEFF_DT, 2e-7)
check("phase coefficient [rad K^-1 m^-1]", PHASE_COEFF, 814.4, 0.1)

print("--- Thermal kernels")
std, uc = StandardKernel(), UndercutKernel()
x = std.x
check("standard K(0) [K per W/m]", std(0, 0), 0.1959, 5e-4)
check("standard P_pi [mW]", np.pi / (PHASE_COEFF * std(0, 0)) * 1e3, 19.7, 0.1)
for dx, ratio in [(2, 0.75), (4, 0.41), (6, 0.22), (8, 0.12)]:
    check(f"standard K({dx})/K(0)", std(dx, 0) / std(0, 0), ratio, 0.01)
check("standard symmetry max|K(x)-K(-x)|/K(0)", np.abs(std(x, 0) - std(-x, 0)).max() / std(0, 0), 0, 1e-3)

check("undercut K(0; 0) [K per W/m]", uc(0, 0), 2.85, 0.01)
check("undercut P_pi [mW]", np.pi / (PHASE_COEFF * uc(0, 0)) * 1e3, 1.35, 0.01)
for dx, ratio in [(2, 0.97), (8, 0.81), (16, 0.61), (32, 0.23)]:
    check(f"undercut K({dx}; 0)/K(0; 0)", uc(dx, 0) / uc(0, 0), ratio, 0.01)
# The mirror rule K(x; -x_h) = K(-x; x_h) needs the x_h = 0 column itself to be symmetric.
check("undercut x_h=0 column symmetry (rel.)", np.abs(uc.k[:, 0] - uc.k[::-1, 0]).max() / uc(0, 0), 0, 1e-3)
# Interpolation in x_h: predict the x_h = 5 column from x_h = 4 and 6 (a 2 µm gap,
# twice what the code ever bridges), relative to the peak.
mid = 0.5 * (uc.k[:, 4] + uc.k[:, 6])
check("undercut x_h interpolation error over 2 µm (rel.)", np.abs(mid - uc.k[:, 5]).max() / uc(0, 0), 0, 0.01)

print("--- Coupled heat-and-mode check (Model C')")
_, cp = read_comsol(DATA_DIR / "cprime_neff.csv")
n_cold = cp[0, 3]
for x_h, q, n_eff in cp[1:, [0, 1, 3]]:
    predicted = DNEFF_DT * std(0, x_h) * q
    check(f"C' heater at {x_h:g} µm: measured/predicted dn_eff", (n_eff - n_cold) / predicted, 1, 0.005)

print("--- Heater-to-phase matrix")
xe = emitter_positions()
centre = N_EMITTERS // 2
check("emitter span [µm]", xe[-1] - xe[0], (N_EMITTERS - 1) * PITCH_UM, 1e-9)
check("centre emitter position [µm]", xe[centre], 0, 1e-12)
check("max distance from emitters to thermal-model cores [µm]", np.abs(xe - 2 * np.round(xe / 2)).max(), 0, 1e-12)
g = phase_matrix(std, xe, [xe[centre]])
check("standard: pi / G for heater over its emitter [mW]", np.pi / g[centre, 0] * 1e3, 19.7, 0.1)
g_left = phase_matrix(uc, xe, [-9.5])
g_right = phase_matrix(uc, xe, [9.5])
check("undercut: mirrored heaters give mirrored G", np.abs(g_left[::-1, 0] - g_right[:, 0]).max(), 0, 1e-12)

print("--- Element pattern (Model A mode field)")
ep = element_pattern()
check("EP(20 deg)", ep(np.sin(np.radians(20))), 0.90, 0.01)
check("EP(22.8 deg)", ep(np.sin(np.radians(22.8))), 0.87, 0.01)

print("--- Far field: analytic limits")
nd = N_EMITTERS * PITCH_UM
flat = beam_metrics(U_GRID, far_field(np.zeros(N_EMITTERS), xe, ep), ep)
check("uniform phase: beam angle [deg]", flat["angle_deg"], 0, 0.01)
check("uniform phase: FWHM in u / (0.886 λ/Nd)", flat["fwhm_u"] / (0.886 * WAVELENGTH_UM / nd), 1, 0.01)
check("uniform phase: first side lobe [dB]", flat["sls_db"], 13.26, 0.2)
check("uniform phase: Strehl", flat["strehl"], 1, 1e-6)
# At d = 2 µm > λ the copies at u = ±λ/d are always visible; their level is set by EP alone.
u_gl = WAVELENGTH_UM / PITCH_UM
check("uniform phase: grating lobe [dB] vs -10 log EP(λ/d)", flat["grating_lobe_db"], -10 * np.log10(ep(u_gl)), 0.05)

for dphi in [0.3, 1.0, 2.0, -1.5]:
    phases = dphi * xe / PITCH_UM
    m = beam_metrics(U_GRID, far_field(phases, xe, ep), ep)
    check(f"ramp {dphi:+.1f} rad/emitter: sin(theta)", m["u_peak"], dphi * WAVELENGTH_UM / (2 * np.pi * PITCH_UM), 2e-4)
    check(f"ramp {dphi:+.1f} rad/emitter: Strehl", m["strehl"], 1, 1e-3)

# Δφ = π: two equal beams at ±arcsin(λ/2d) = ±22.8 deg, the edge of the single-beam field of view.
ff = far_field(np.pi * xe / PITCH_UM, xe, ep)
u_edge = WAVELENGTH_UM / (2 * PITCH_UM)
check("ramp pi: beam angle magnitude [deg]", abs(beam_metrics(U_GRID, ff, ep)["angle_deg"]), 22.8, 0.05)
check("ramp pi: I(+edge) / I(-edge)", np.interp(u_edge, U_GRID, ff) / np.interp(-u_edge, U_GRID, ff), 1, 1e-3)
# Just past π the copy on the other side is nearer broadside, so it is the stronger beam.
dphi = np.pi + 0.2
m = beam_metrics(U_GRID, far_field(dphi * xe / PITCH_UM, xe, ep), ep)
check("ramp pi+0.2: beam jumps to (Δφ-2π)λ/(2πd)", m["u_peak"], (dphi - 2 * np.pi) * WAVELENGTH_UM / (2 * np.pi * PITCH_UM), 2e-4)

print("--- Strehl against the Ruze relation")
rng = np.random.default_rng(0)
for sigma in [0.3, 0.7, 1.2]:
    s = np.mean([strehl(rng.normal(0, sigma, N_EMITTERS), xe, 0.0) for _ in range(4000)])
    # Finite N adds a floor: <S> = exp(-σ²) + (1 - exp(-σ²))/N.
    expected = np.exp(-sigma**2) + (1 - np.exp(-sigma**2)) / N_EMITTERS
    check(f"mean Strehl, sigma = {sigma} rad", s, expected, 0.01)

print("--- Least-squares steering")
g4 = phase_matrix(uc, xe, [-30, -10, 10, 30])
ramp = xe / PITCH_UM
p1, c1 = fit_powers(g4, ramp)
p2, c2 = fit_powers(g4, 2.5 * ramp)
check("scaling: P(2.5 x target) / (2.5 P(target))", p2.sum() / (2.5 * p1.sum()), 1, 1e-6)
check("scaling: offset(2.5 x target) / (2.5 offset(target))", c2 / (2.5 * c1), 1, 1e-6)
fit = steer_one_direction(g4, xe, +1)
dphi = 1.3
u_target = dphi * WAVELENGTH_UM / (2 * np.pi * PITCH_UM)
check("Strehl shortcut equals farfield.strehl at the target",
      strehl_vs_step(fit["residual_per_rad"], [dphi])[0], strehl(dphi * (g4 @ fit["powers_per_rad"]), xe, u_target), 1e-9)
# Optimality of the fit: for every heater actually in use (P > 0), changing its power
# slightly cannot reduce the error, so its column is orthogonal to the residual.
check("residual is orthogonal to the fitted columns (unit ramp)",
      np.abs(g4[:, fit["powers_per_rad"] > 0].T @ fit["residual_per_rad"]).max(), 0, 1e-8)

print("--- Conventional benchmark")
p_pi = np.pi / (PHASE_COEFF * std(0, 0))
check("isolated heaters, broadside: power [W]", ideal_power(xe, 0.0, p_pi), 0, 1e-12)
# Δφ = π alternates 0, π: the least-power choice puts π on the 16 emitters at odd x/d.
check("isolated heaters, Δφ = π: power / P_π", ideal_power(xe, np.pi, p_pi) / p_pi, 16, 1e-9)
# With an isolated (diagonal) G the least-squares route must reproduce the exact answer.
g_iso = np.eye(N_EMITTERS) * np.pi / p_pi
power_iso, s_iso, _ = with_crosstalk(g_iso, xe, 1.0, WAVELENGTH_UM / (2 * np.pi * PITCH_UM))
check("diagonal G: fitted power / exact isolated power", power_iso / ideal_power(xe, 1.0, p_pi), 1, 1e-6)
check("diagonal G: Strehl", s_iso, 1, 1e-6)

print(f"\n{failures} check(s) failed" if failures else "\nAll checks passed")
print(f"Note: at broadside the grating lobes at ±{np.degrees(np.arcsin(u_gl)):.1f} deg sit only "
      f"{flat['grating_lobe_db']:.2f} dB below the main beam; main-lobe power fraction {flat['main_lobe_fraction']:.2f}.")
raise SystemExit(failures)
