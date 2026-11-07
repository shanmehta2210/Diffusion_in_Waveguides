# Prep for the chosen avenues (2 Nov 2026, night)

Shreyansh picked three avenues. Nothing below has been started; this note is preparation only. Numbers marked *est.* are rough.

---

## 1. The lens: comparisons against conventional approaches

**What we know so far.**
- Undercut stack, one centre heater: up to 5.5 rad edge-to-centre (f ≈ 0.38 mm) at 0.56 mW per rad.
- Five heaters: up to 24 rad (f ≈ 85 µm) at about 2.2 mW per rad.
- The lens is 1D (cylindrical) only. The other axis is set by the grating.

**Caveat to check first: conventional arrays can wrap the lens (a Fresnel lens).** With one heater per channel, the phases only need to cover 0–2π.
- **Weak lenses (a ≲ 2π).** Wrapping saves nothing.
  - Isolated undercut heaters need ≈ 21.3·a/π·P_π ≈ 50 mW at 5.5 rad.
  - Ours needs ≈ 3 mW, so about 16× less (est.).
- **Strong lenses (a = 24 rad).**
  - A wrapped conventional lens costs about N·P_π ≈ 45 mW.
  - Ours costs about 53 mW, so the two are roughly equal (est.).

**Expected conclusion:** our advantage is largest for weak-to-moderate focusing, f ≳ 0.3 mm.

**Competitors to compare against:**

| Competitor | Comparison |
|---|---|
| Isolated heaters, unwrapped and wrapped (standard and undercut) | Power, number of wires, beam quality |
| A Doerr-style fixed-shape lens heater (graded widths, one wire) | Same number of wires, but the lens shape is frozen. Ours can also steer at the same time |
| Wavelength-tuned focusing (DeSantis 2024, spiral OPA) | No phase-shifter power, but needs a tunable laser, and focus is tied to wavelength |
| Off-chip tunable lenses | Optional; needs a literature check |

**Metrics:**
- focal-length range;
- power per unit of optical power (1/f);
- Strehl, or focal-spot width against the diffraction limit;
- number of wires;
- response time;
- steering and focusing at once (a joint non-negative fit).

**Python:** quadratic targets with the existing kernels. A wrapped-quadratic version of `conventional.ideal_power` gives the Fresnel baseline.

---

## 2. Response time and multi-heater pre-emphasis

**What the applications need.**
- FSO tip-tilt: 1–70 Hz loop bandwidth (Frost 2024).
- Holographic tweezers: about 20 Hz (agent, abstract only).
- Retinal stimulation and FSO acquisition: speed not stated.

**What we have now.**
- The undercut rises in 0.58 / 1.63 ms (63% / 90%) at the heated core, and 0.97 / 2.02 ms at 32 µm (Model D).
- The ramp is the slowest pattern: τ ≈ 1.7 ms in air (est.), so f_3dB ≈ 90 Hz. That may already be enough for the slow niches. Pre-emphasis would widen the margin and open faster ones.

**Physics.**
- A temperature pattern with spatial wavenumber q decays at a rate αq² + h/(ρct), so broad patterns are slow.
- Pre-emphasis means driving harder at first, then backing off.
- **Constraint:** heater power is ≥ 0, so cooling cannot be forced. Speeding up cooling therefore needs a static bias, so the system runs around an operating point and can drive both up and down. The cost is static power.

**Plan.**
1. **A time-domain membrane surrogate:** ρc t ∂T/∂t = k t T'' − hT + Σ q_k(t) δ(x − x_k).
   - Fit ρc_eff (and check h) against the undercut columns of `transient.csv`. Those columns only cover a heater at x = 0, at 6 positions.
2. **Input design as convex optimisation.** The system is linear, so time-discretised heater waveforms P_k(t) enter linearly.
   - Find the shortest horizon T for which a feasible drive exists. Feasible means it reaches the target ramp within a set Strehl, with 0 ≤ P_k ≤ P_max (optionally around a bias).
   - This is a linear or quadratic program, bisected over T.
   - Compare: a single step, single-heater overdrive, and K-heater pre-emphasis.
3. **Report:**
   - settling time against peak-power budget and bias power;
   - energy per switch;
   - the beam's path during the transition (angle and Strehl against time).
4. **COMSOL confirmation.** Two options:
   - a time-dependent run with the optimised waveforms on the full 2D model;
   - the frequency-domain run C4, giving K(x; ω). Its transfer functions validate the surrogate dynamics directly.

---

## 3. Slot dissipation, and heaters far out and beside the core (COMSOL build spec)

**Base model:** `modelC_kernel.mph`, undercut (uc = 1, W_uc = 80 µm, D_uc = 15 µm).

| Case | Change | Heater positions | Purpose |
|---|---|---|---|
| S0 | No slot | x_h = 32 … 40 µm, step 1 | "Heaters out far": tests the 16.7° saturation directly |
| S1 | Slot through the whole membrane (cladding + BOX, opening into the cavity) from x = 36 to 38 µm. New parameters x_slot = 37 µm, w_slot = 2 µm | x_h = 33, 34, 35 (outside the emitters, inside the slot) **and** x_h = −32 … 32, step 2, with the slot present | The ramp membrane. The interior columns are needed because the slot changes every heater's kernel. **Mirror symmetry is lost**, so negative x_h must be simulated, not mirrored |
| S2 | As S1, with D_uc = 50 and 100 µm | x_h = 34, plus a coarse interior set (step 4) | Less air loss, straighter ramp |
| S3 | As S1, with k_air → 1e-6 (vacuum limit) | x_h = 34, plus a coarse interior set | Upper bound on linearity and efficiency |
| S4 | Heater beside the core (your geometry), standard and undercut, plus an optional Si slab of 0 / 90 / 150 nm | Your choice | Heater placement in the plane of the core |

**Export.** Same format as `kernel_undercut.csv`:
- temperature at core height, y = 0.11 µm;
- x from −150 to 150 µm in 0.1 µm steps;
- one column per x_h;
- a header naming the case.

**Predictions the run will test** (membrane surrogate, est.):
- **S1:** heater at 35 plus one interior heater gives 17.8° one way at about 0.81 mW/deg.
- **S2/S3:** a single heater at 35 gives the full 22.8° at about 0.55 mW/deg.
- **S0:** no great gain unless the heater sits beyond the last emitter.

**Practical notes.**
- Keep the slot at least 2 µm wide; Parsons 2025 saw damage when trenches were only 1 µm apart.
- The heater (2 µm wide) has to fit between the last core edge (32.25 µm) and the slot.

---

## Not started: heaters over the gratings

Shreyansh wants to study this first. The open questions are in `beyond_steering.md` §2.4, row D:
- uniform against graded heating;
- the coupling between the two axes;
- how an undercut under the grating changes the up/down emission split.
