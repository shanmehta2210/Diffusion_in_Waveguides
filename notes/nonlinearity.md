# Nonlinearity at high temperature (7 Nov 2026)

The brief's analysis 8. Code: `analysis/opa/nonlinear.py`, `analysis/run_nonlinearity.py` (a few seconds). Data: Model F (`nonlinear_line.csv`), Model A (`neff_vs_T.csv`), `heater_Tmax_slot.csv`. Results: `analysis/results/nonlinear_designs.csv`.

## Bottom line

- **The heating nonlinearity follows a simple pointwise law.** The real rise is ΔT = θ − aθ², where θ is the linear model's rise and a = 3.0×10⁻⁴ per K on the membrane. It reproduces every Model F undercut profile to within 0.25% of the peak (up to 187 K), and two heaters together to 0.14%. So it applies to any multi-heater drive: superpose linearly, then transform.
- **On the beam it acts almost entirely as a gain error.** With 1 mm heaters at 10° and our single linear fit:
  - the beam lands 0.22° short (heating alone) to 0.28° long (if silicon's dn/dT grows as fast as the literature suggests);
  - the Strehl toward the target falls from 0.96 to 0.84–0.89 (slot).

  Rescaling the drive by one number per angle (a gain of 0.975–1.022 at 10°) brings it back to 0.95–0.96. Re-optimising every heater adds at most 0.008 at 10°.
- **So the calibration claim gets a small amendment.** It becomes: one linear fit, plus a one-dimensional gain curve against angle. That is still far simpler than per-element calibration.
- **The brief's flags:** with 1 mm heaters the hottest core passes 30 K above about 5° and 80 K above about 13°.
  - Up to 10° the gain correction leaves at most about 0.01 of Strehl on the table.
  - At 15–20° (90–120 K) some distortion of the shape remains. The slot design loses up to 0.07 after gain correction at 20°, recovered partly by a full refit or by longer heaters (with trim, longer heaters are acceptable).
- **Heaters run close to the hottest core in steady state:** 37, 47 and 65 K for the slot design at 10°, against a hottest core of 61 K. The pre-emphasis bursts (3× the hold power in one heater for about 0.2 ms) have a steady-state bound of 160–200 K. The real transient peak is lower, but it was not simulated with temperature-dependent properties.
- **Open:** how fast silicon's dn/dT grows with temperature. I couldn't verify the literature value today, so it is bracketed. The heating part, from our own Model F, is solid.

## Physics

**Heating: why a pointwise law.**
- Oxide conducts better when hot, so a heated membrane spreads heat more easily and the rise falls below the linear prediction.
- For steady conduction with k = k₀(1 + βΔT), the Kirchhoff transform θ = ∫ k/k₀ dT = ΔT + βΔT²/2 satisfies the *linear* heat equation, with the same sources and boundary values. So the linear model predicts θ exactly, and inverting gives ΔT = θ − (β/2)θ² to first order.
- It is exact only for one conducting material with fixed-temperature boundaries. Air under the membrane and the silicon cores break that slightly, which is why the fit is good to 0.25% rather than to zero.
- On the standard stack the profile is not purely local (far from the heater the real rise is slightly *above* linear), but the law still fits to within 0.19%.

**Optics: dn/dT grows with temperature.**
- Our COMSOL index model used a constant dn/dT. Model A's n_eff is linear in T, with curvature 2c₂/c₁ = 2.3×10⁻⁵ per K.
- Real silicon's dn/dT rises with temperature. Write its relative growth as b = (1/(dn/dT))·d(dn/dT)/dT. The index change is then (dn/dT)(ΔT + bΔT²/2).
- My recollection of the Cocorullo et al. 1999 fit (300–550 K at 1523 nm; Appl. Phys. Lett. 74, 3338, DOI 10.1063/1.123337) gives b ≈ 1.4×10⁻³ per K near 300 K.
  - I confirmed the paper exists, but neither search results nor a mirror of Komma et al. 2012 gave the fitted coefficients.
  - So b is bracketed at 0, 0.7 and 1.4×10⁻³ per K. **Before the report, check the coefficients in the paper itself.**

**Together.** The phase goes as φ ≈ φ_lin·(1 + (b/2 − a)·θ) to first order.
- The two effects pull in opposite directions; they nearly cancel at b ≈ 2a = 0.6×10⁻³ per K.
- Across a ramp, θ grows linearly with x, so the excess phase grows as x², a parabola.
- A parabola fitted with a straight line leaves only √(1/180) ≈ 7.5% of its end-to-end height as rms residual. So the beam mostly *moves* (a gain error) and barely *blurs*. That is why one gain per angle fixes almost everything.

## Checks

| Check | Result |
|---|---|
| Kirchhoff law against Model F, undercut, 4–70 W/m (peaks 11–187 K) | worst error 0.05–0.24% of peak; the linear model is off by 0.4–6.5% |
| Same, standard stack, 50–1000 W/m | 0.04–0.19%; linear 0.38–7.05% |
| Two heaters together | plain sum off by 2.84% (undercut) and 1.73% (standard), matching the brief; Kirchhoff law 0.14% and 0.08% |
| Brief's "1% at 30 K, 3% at 80 K" | 0.9% and 2.4% with a = 3.0×10⁻⁴; the brief's 3% was a rounded reading |
| Fitted a | undercut 3.0×10⁻⁴ per K, standard 3.5×10⁻⁴ (brief: about 3.3×10⁻⁴) |

The final designs (slot, vacuum) have no Model F run of their own. They use the undercut's a: the same oxide membrane carries the heat, and the Kirchhoff argument is cleanest for exactly that 1D membrane conduction.

## Results: final designs, 1 mm heaters

Strehl toward the target, under three calibrations:
- **linear fit:** our single linear fit;
- **gain:** the same drive rescaled by one number so the beam points right;
- **refit:** every heater re-optimised.

b is in units of 10⁻³ per K.

| Design | Angle | Hottest core | Linear model | b | Pointing error, linear fit | Linear fit | Gain | Refit |
|---|---|---|---|---|---|---|---|---|
| Slot, 3 heaters | 5° | 31 K | 0.990 | 0 / 1.4 | −0.05° / +0.07° | 0.985 / 0.980 | 0.990 / 0.989 | 0.990 / 0.989 |
| | 10° | 61 K | 0.960 | 0 / 1.4 | −0.22° / +0.26° | 0.887 / 0.837 | 0.961 / 0.950 | 0.965 / 0.953 |
| | 15° | 91 K | 0.912 | 0 / 1.4 | −0.49° / +0.58° | 0.608 / 0.474 | 0.912 / 0.879 | 0.929 / 0.893 |
| | 20° | 119 K | 0.851 | 0 / 1.4 | −0.87° / +1.01° | 0.227 / 0.101 | 0.837 / 0.780 | 0.887 / 0.811 |
| Vacuum, 1 heater | 10° | 63 K | 1.000 | 0 / 1.4 | −0.22° / +0.28° | 0.913 / 0.869 | 0.996 / 0.995 | same as gain (one heater) |
| | 15° | 93 K | 1.000 | 0 / 1.4 | −0.50° / +0.60° | 0.628 / 0.508 | 0.981 / 0.981 | same |
| | 20° | 122 K | 1.000 | 0 / 1.4 | −0.90° / +1.05° | 0.218 / 0.118 | 0.940 / 0.951 | same |
| Undercut, 3 heaters | 10° | 62 K | 0.895 | 0 / 1.4 | −0.23° / +0.26° | 0.848 / 0.774 | 0.897 / 0.884 | 0.897 / 0.892 |

At b = 0.7 the two effects almost cancel: the pointing error is under 0.1° at every angle (`nonlinear_designs.csv`). At 2° nothing changes (13 K).

**Heater length (slot, 10°):**

| Heater length | Hottest core | Linear fit | Gain-corrected |
|---|---|---|---|
| 0.5 mm | 120 K | 0.59–0.69 | 0.94–0.96 |
| 1 mm | 61 K | 0.84–0.89 | 0.95–0.96 |
| 2 mm | 31 K | 0.92–0.94 | 0.96 |

Longer heaters cut the nonlinearity in proportion. With trim (the design requirement from `fabrication.md`), the extra waveguide length costs nothing in fixed errors.

## Heater temperatures (slot design, 1 mm heaters)

`heater_Tmax_slot.csv` gives each heater's hottest point per W/m. The other heaters' rise at its position is added linearly, then the Kirchhoff law is applied.

| Case | Heaters at 2, 14, 32 µm |
|---|---|
| Steady 10° (1.40, 0.91, 11.45 mW) | 37, 47, 65 K (hottest core 61 K) |
| Pre-emphasis burst, 41.3 mW in one heater, steady-state bound | 161, 178, 202 K |

- The bursts last about 0.1–0.2 ms, so the real peak is below the bound.
- Metal heaters tolerate these temperatures.
- The open point is accuracy. The pre-emphasis model is linear. During a burst the heater region is 5% or more into the nonlinear range, so the boost delivers slightly less heat than planned. The robust drive (±5% time-scale band) absorbs some of this, but a nonlinear transient was not simulated.

## What this means for the argument

- **Amend the calibration claim:** "after a one-time trim, one linear fit plus one gain curve against angle covers every direction". Still O(1) parameters per angle, against O(N) for conventional arrays.
- **The nonlinearity doesn't threaten the headline numbers up to 10°.** At 15–20° with 1 mm heaters it costs up to 0.07 of Strehl after gain correction; longer heaters fix that.
- **The vacuum design's single heater can't correct shape distortion** (a refit equals the gain correction). That is why it ends at 0.94–0.95 at 20° while its linear model says 1.00.

## Caveats

- **b is unverified.** Check the Cocorullo 1999 coefficients before quoting the + sign of the pointing error. The heating part (a) comes from our own COMSOL Model F.
- **No Model F run exists for the slot or vacuum geometry,** so the undercut's a is used. Air in the slot conducts more when hot as well, which would add to a slightly.
- **Pre-emphasis transients were not re-simulated with temperature-dependent properties.**
- **Heater temperature uses the core-height kernel for the other heaters' contribution:** a slight underestimate.
