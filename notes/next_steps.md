# Next steps: what we found and what we can do (2 Nov 2026, due 11 Nov)

This note is the roadmap. Details and references are in two other notes:
- `literature_review.md`: crosstalk handling, heater sizes, prior art.
- `beyond_steering.md`: applications, geometries, simulation list, fabrication errors.

Numbers marked *est.* are quick calculations, not results.

---

## 1. Where we stand

**Done.** All 53 pipeline checks pass, and three analyses are complete:
- steering with evenly spaced heaters;
- optimised heater placement;
- the benchmark against conventional designs.

**Key results:**
- Undercut stack, 4 heaters: ±16.0° at Strehl ≥ 0.8 and about 4.9 mW/deg. The range saturates at ±16.7° for 5 or more heaters.
- Standard stack with few heaters: under 3.4° of span.

**What the literature changed:**

1. **Fewer heaters do not let us pack emitters closer.**
   - Conventional OPAs already separate the two. Their heaters sit 15–100 µm apart and the waveguides then converge to emitter pitches of 0.78–3.5 µm (Phare 2018, Miller 2020, Doylend 2011).
   - Emitter pitch is limited by optical coupling between waveguides, not by heat.
2. **Prior art for few control wires already exists.**
   - Van Acoleyen 2009: one shaped electrode over 16 gratings at 2 µm pitch, steering 2.3°.
   - Yaacobi 2014 and Poulton 2017: phase shifters cascaded on the feed waveguide, steered by 1–3 wires.
   - Hosseini 2025: meandered heaters, 4 wires.
   - Nobody uses diffusion from a few optimised heaters. Our idea is new, but the gap is narrow.
3. **A fair power comparison changes the headline.** The earlier "12× less power" compared us with *non-undercut* conventional heaters.
   - Against *undercut* isolated heaters, we win below about 8° (6.5× at 1°, 1.6× at 5°) and lose above it.
   - A cascade (est.) beats today's design by about 2.5× at every angle. It does so with very short, very hot stages: Yaacobi reports about 500 °C.

**Honest positioning:** small angles, plus phase profiles that are not ramps.

---

## 2. The physics behind every next step

**(a) Heat is a spatial low-pass filter.**
- Steady heat conduction smooths whatever the heaters put in, so we can only make smooth phase profiles: an offset, a ramp, a dome (quadratic), slow ripples.
- We cannot make the sharp 2π resets that conventional arrays use (phase wrapping), nor per-waveguide random patterns.
- Two consequences:
  - The full ramp has to be built physically: 45 rad end to end at 10°. That sets the temperature.
  - Fabrication errors split into a smooth part we can correct and a local part we cannot.

**(b) The undercut is a 1D membrane.**
- About 4.2 µm of oxide spans an 80 µm cavity. It is held at ambient where it meets the cavity walls and loses a little heat into the air gap below.
- The 1D model T'' − T/L² = −q/(k t), with healing length L = √(k t / h_air), fits COMSOL with only two parameters:
  - L = 63 µm (first-principles estimate 58 µm);
  - k t = 5.8×10⁻⁶ W/K (first-principles 5.82×10⁻⁶).
- It matches COMSOL to within 7–9% of the peak and reproduces the one-heater steering result (4.6° against 4.4°).
- This is our cheap design tool for geometries COMSOL hasn't simulated.

**(c) A heater at the centre is a lens.**
- A central heater makes a dome: highest in the middle, falling towards the walls, close to a parabola near the centre.
- Silicon's dn/dT > 0, so the middle gets more optical path, as in a convex glass lens. The lens therefore converges.

**(d) A slot turns the profile into a straight ramp.**
- Between a heater and a heat sink, with no sources or losses in between, the heat equation reduces to T'' = 0: a straight line.
- Today the heater sits inside the cavity, so heat flows both ways. The profile is a tent whose peak lies inside the array.
- Instead, cut a slot through the membrane just past a heater placed *outside* the emitter span. Heat can then only escape by crossing the whole array to the far wall, so every emitter sits on the straight part.
- Air loss bends the line slightly (it becomes a sinh). A second heater, or a longer L, removes that bend.

**(e) Scaling law for the slotted ramp.**
- All the heater power must cross the membrane: P/L_h = k t · dT/dx (Fourier's law), where L_h is the heater length.
- The steering angle fixes the phase gradient: dφ/dx = Δφ/d = 814 · L_h · dT/dx.
- Combining the two: **P = k t (Δφ/d) / 814**. This depends on neither N nor heater length. It gives about 5 mW at 10°.
- What does grow with array size:
  - the hottest-core temperature, ΔT_max = 32Δφ / (814 L_h), which grows with N;
  - the response time, ~ span²/α, which grows as N².

**(f) Broad profiles are the slowest.**
- A profile with spatial wavelength Λ decays at rate α(2π/Λ)² + h/(ρct).
- The ramp is the broadest profile we use, so it is the slowest: about 1.7 ms in air and 3 ms in vacuum (est.).
- Driving several heaters with a timed sequence (pre-emphasis) could speed it up. I found no published multi-heater version.

**(g) Fabrication errors.**
- Variations in waveguide width and thickness change n_eff along each waveguide.
- **The smooth part is correctable.** About 80% of width variance is smooth across a wafer (Xing 2022). Across a 66 µm array it looks like an offset plus a gradient plus some curvature, which our heaters correct. In simulation, with a correlation length ≥ 16 µm, 4–8 heaters bring Strehl to ≥ 0.98.
- **The local part is not.** It is independent from waveguide to waveguide and accumulates like a random walk: σ² ≈ L / L_coh, with L_coh = 4.17 mm for strip waveguides and 1.61 mm for rib (Yang 2015). Our heaters cannot touch it.
- **This creates a trade-off unique to our design.** Without wrapping we need long heaters to keep cool, but long heaters mean long independent waveguides and more local error:

| Maximum angle | Heater length for a 30 K rise | Strehl from local error (strip) |
|---|---|---|
| 5° | 0.9 mm | 0.80 |
| 10° | 1.8 mm | 0.64 |
| 16° | 2.9 mm | 0.50 |

At small angles the trade-off is mild, which is another reason to favour the small-angle niche.

**(h) The membrane bends when heated.**
- A suspended stack of different layers bows when heated. An out-of-plane displacement δz adds a phase of about k₀ δz, so 100 nm gives 0.4 rad. A parabolic bow acts as defocus.
- The photoelastic (stress-optic) effect is about 3.5% of the thermo-optic effect and follows ΔT, so it only rescales dn/dT and calibration absorbs it.

---

## 3. Python work with the data we already have

| # | Task | What it gives the report |
|---|---|---|
| P1 | **Small-angle niche** (≤ ±8°). Power, Strehl, hottest-core temperature against heater length, and fabrication-error Strehl, compared with isolated undercut heaters and the cascade (labelled as an estimate). | The main result: where we win and why |
| P2 | **Tunable lens.** Quadratic targets with the existing kernels; range of focal lengths, power per rad, beam quality; steering and focusing together (a joint non-negative fit). *So far:* one centre heater reaches 5.5 rad at 0.56 mW/rad (f ≈ 0.38 mm); 5 heaters reach 24 rad (f ≈ 85 µm). Isolated heaters need about 9.2 mW/rad. | The second application, and something a cascade cannot do |
| P3 | **Fabrication errors** (the brief's analysis 6, upgraded). Error covariance = smooth part (Xing) + local part (Yang). Monte Carlo in three cases: uncorrected, corrected by the K heaters, trimmed. Include the heater-length trade-off and the comparison with conventional, including how many calibration parameters each needs. | Answers "more, less or equally harmful" |
| P4 | **Wider heaters by superposition.** Heat conduction is linear, so a heater 2m µm wide is exactly m adjacent 2 µm heaters at equal power. Widths of 2–10 µm need no COMSOL. Heater length enters only as a temperature constraint. | Answers the heater-size question |
| P5 | **Heated-section pitch d_h** (the design maps, reinterpreted). The standard kernel is shift-invariant, so we just resample it. For the undercut, use the surrogate, because the cavity scales with the span. | Power-against-speed design map |
| P6 | **Surrogate screening** of ramp-membrane variants before COMSOL. *So far:* with a slot and today's air gap, 2 heaters give 17.8° one way at 0.81 mW/deg (≈ 8 mW at 10°). With L ≥ 150 µm, one heater covers the full 22.8° at 0.55 mW/deg. | Which geometries are worth COMSOL time |
| P7 | **Remaining items from the brief:** singular values of G (5), speed and energy (7), nonlinearity (8). The singular vectors are now useful: they show which shapes (offset, ramp, dome) are cheap. | Completes the brief |

---

## 4. New COMSOL runs, in priority order

| # | Run | Build | Decides |
|---|---|---|---|
| C1 | **Ramp membrane** | A slot through the membrane just outside the last emitter. Heater at x_h = 34–38 µm. Cavity depth 15, 50 and 100 µm, plus a vacuum case (k_air → 0). Export K(x; x_h) in the same format as `kernel_undercut.csv` | Whether the 0.5–0.8 mW/deg prediction survives. Also covers the pending "heaters beyond 32 µm" run |
| C2 | **Heaters beside the core** (your run) + **silicon slab sweep** of 0, 90 and 150 nm, standard stack | — | Whether we get broad spreading at µs speed. Half-width est. 4.4–6.9 µm, against 3.3 µm today |
| C3 | **Thermo-mechanical bow** | Heat Transfer + Solid Mechanics with thermal expansion and residual oxide stress. Export out-of-plane displacement at core height against x | Whether bending is a hidden error |
| C4 | **Frequency-domain heat transfer** | Same 2D model; K(x; ω) for 10 Hz – 10 kHz | The bandwidth of each profile; the basis for pre-emphasis |
| C5 | **3D: finite cavity and heater, with the grating region** | — | Temperature at the first grating period (the Si cores carry heat along the waveguide about as well as the oxide); extra power lost through the cavity ends (est. ~8% for a 1 mm heater); whether 2D is valid (keep the heated region 120–150 µm from the cavity ends) |
| C6 | **Sensitivities** | Cavity width ±2 µm, membrane thickness ±5%, heater offset ±0.2 µm: about 6 runs | δG for the error model |
| C7 | **Mode solver** | dn_eff/dw for waveguides 0.5 µm against 1–2 µm wide | How much wider heated waveguides cut local error |
| C8 | *Optional:* heater over the grating region (uniform and graded), and Model E re-run over a membrane | — | Steering and focusing along the grating at a fixed wavelength |

---

## 5. Design options

Heaters stay **before the gratings**. That frees the heated section's pitch d_h and waveguide width from the emitter pitch:
- d_h trades power against speed: P ∝ 1/d_h and τ ∝ d_h². The lower limit is 1.5 µm.
- Wider heated waveguides lower the fabrication error.

| Option | Verdict |
|---|---|
| Steering heaters over the gratings | **Reject.** It needs about 370 K at 10°, tilts the emitters about 4° along the grating, and smears the beam |
| Ramp section followed by a lens section | **Promising.** Each cavity can be designed for one job |
| Today's design plus 1–2 heaters over the gratings | **Good extension.** Gives 2D steering at a fixed wavelength with 4–6 wires |
| Folded or multi-pass heated section | **Low priority** |

---

## 6. Applications to motivate the report

1. **Tunable 1D focus.** Tweezers, retinal stimulation, wide-beam acquisition followed by narrow-beam tracking in free-space optical (FSO) links.
2. **Tuning arrayed waveguide gratings (AWGs).** A ramp shifts the centre wavelength and a dome sets the dispersion. Precedent: Doerr's thermo-optic lens, in silica at watt-level power.
3. **Slow tip-tilt correction and alignment** in FSO links. Required loop bandwidth is 1–70 Hz (Frost 2024).
4. **2D steering at a fixed wavelength**, by adding a grating heater (Kim 2019).
5. **Microring filter banks:** a gradient corrects the systematic part only.

Poor fits: adaptive optics, coherent beam combining, programmable meshes.

---

## 7. Fabrication errors: the verdict in brief

| Error | Ours compared with conventional |
|---|---|
| Local random | **More harmful**: needs trim, or short and wide heated sections |
| Smooth wafer-level variation | **Equal** |
| Heater gain spread | **Less**: K gains to calibrate instead of N |
| Calibration effort | **Less**: O(K), and one fit covers all angles |
| Ageing at high temperature | **More**: unwrapped phases run hotter |
| Kernel-shape errors | **Untested** (C6) |

Everything is linear and cheap to model. Only δG (C6) and the bow (C3) need COMSOL.

---

## 8. Decisions for you

1. **Framing.** Small-angle niche + lens + scaling-law physics?
2. **Which COMSOL run first:** C1 (ramp membrane) or C2 (your side heaters)?
3. **The cascade:** include it in the report as a clearly labelled estimate?
4. **Static trim:** keep assuming it, or limit the claims to small angles, where working without trim is plausible?
   - **Decided 7 Nov: trim is a stated requirement of the design.** Analysis 6 (`fabrication.md`) shows the heaters can't correct independent errors, and an untrimmed 10° array is limited to Strehl about 0.78 by its heated section alone.

---

## 9. Suggested schedule (the brief's timeline, adjusted)

| Dates | Python | COMSOL (you) |
|---|---|---|
| 2–4 Nov | P1–P4, P6 | C1, C2 |
| 5–6 Nov | Fold in C1/C2; P5, P7 | C3, C4 if time |
| 7–8 Nov | Figures | — |
| 9–11 Nov | Write-up | — |

C5–C8 are stretch goals.
