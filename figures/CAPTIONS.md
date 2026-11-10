# Figures: draft captions

To regenerate everything, run `python make_figures.py` from `analysis/`. To redo one figure, run `python make_figures.py fig05`.

Each figure is saved twice: a `.pdf` for LaTeX (vector, fonts embedded) and a `.png` at 220 dpi for slides. The colours stay the same in every figure:
- grey: the standard stack, or a heater on every emitter;
- aqua: undercut;
- blue: slot in air;
- orange: vacuum;
- yellow: dual-use.

---

**fig01_stacks.** Where the heat leaves, in the four structures studied (cross-sections; heater and core thickness exaggerated). Above each cross-section is the profile it makes at 10°.
- (a) Standard stack: heat drops 2 µm into the substrate, so each emitter needs its own heater, set to its wrapped phase (0–2π).
- (b) Undercut: with the substrate removed, heat spreads sideways through the membrane to both cavity walls. The edge heater, however, sends most of its heat to the nearby wall (≈90%, 1D estimate).
- (c) Slot: a slot through the membrane past the edge heater closes that short exit, so all the heat crosses the array and the temperature falls as one ramp. Air still leaks heat across the slot and into the cavity.
- (d) Vacuum: removing the gas leaves the far wall as the only exit.

Also saved one panel per stack: `fig01_standard`, `fig01_undercut`, `fig01_slot`, `fig01_vacuum`.

**fig02_kernels.** Temperature profile of one 1 mm heater.
- (a) Heater at the array centre. The standard stack confines heat to a half-width of 3.3 µm; the undercut spreads it to 20.6 µm (10 pitches) and is 14× hotter per mW at the peak.
- (b) Heater at the array edge (34 µm), K per mW. Without a slot, little heat reaches the array. A slot in air gives a near-linear ramp. In vacuum the ramp is straight and 2.5× steeper.
- Source: 2D COMSOL kernels (`kernel_checks`, `kernel_undercut`, `_far`, `_slot`, `_vacuum`).

**fig03_ramp_and_beams.** The slot design (heaters at 2, 14 and 32 µm).
- (a) Phase across the array at 10°. The edge heater does most of the work; the two interior heaters straighten the ramp. Their sum matches the ideal ramp (Strehl 0.96).
- (b) Far field steered to 0–20° in 5° steps with the same three heaters, under the single-emitter pattern. Grating lobes appear outside ±22.8°, set by the 2 µm pitch.

**fig04_power_vs_angle.** Power to hold each angle, N = 33 at 2 µm pitch.
- Grey: a heater on every emitter (33 wires), phases wrapped, at three phase-shifter efficiencies. Those curves wiggle because the best wrap pattern changes with angle; that is real for a finite array.
- Ours: 1–3 wires. The undercut layout's range stops at 15.9° because its kernel data ends at 32 µm.
- All 2D simulation. On a real chip, apply the ×0.67–1.5 simulated-to-measured band and ×1.5–1.9 for a real stack (Q1, Q3).

**fig05_power_vs_speed.** Power to hold 10° against the time to settle a 0 → 10° step.
- Coloured dots: our design path. Each step cuts hold power (580 → 49 → 14 → 5.3 mW) and slows the response.
- Rings: the same design with pre-emphasis, chip power ≤ 3× hold (6–12× faster).
- Grey: published single heaters, scaled to 33 isolated wrapped heaters at their measured P_π; time is the slower of rise and fall.
- Dashed lines: equal power × time.
- The standard stack's time is estimated from one heater's step to 92% (`transient.csv`).

**fig06_preemphasis.** Pre-emphasis, cold start 0 → 10°.
- (a, b) Chip power for a plain step and for the optimised drive, capped at 3× hold.
- (c, d) Strehl toward 10°. Settling (Strehl ≥ 0.8 from then on) drops from 3.42 to 0.43 ms for the slot design and from 8.54 to 1.43 ms in vacuum.
- Source: `preemphasis_traces.csv`.

**fig07_modes.** Which phase patterns heat can make.
- (a) Singular values of each phase matrix, with a heater under every emitter and the free offset removed, relative to the largest. Wide spreading leaves only a few usable patterns, those costing at most 10× the cheapest: 15 for the standard stack at 2 µm pitch, 4, 3 and 2 for the membranes.
- (b) The slot design's three cheapest patterns, a ramp, a dome and a ripple, with their relative power cost.

**fig08_gas_pressure.** Gas in the cavity and slot, from the 1D membrane model calibrated to COMSOL (est.).
- (a) Steering power against pressure for four cavity and slot sizes. It approaches the conduction limit k·t·Δφ/(814·d) = 0.506 mW/deg below about 100 Pa. Diamonds are the COMSOL runs R2 (air) and R3a (vacuum).
- (b) Vacuum design, drive set at 1 Pa: how far the beam falls short of 10° if the package pressure rises. It passes a tenth of the beam width at about 26 Pa.

**fig09_fabrication.** Correcting fabrication phase errors (σ = 1 rad rms, Gaussian correlation length ℓ) by refitting the heater powers, beam at 10°. Mean of the Monte Carlo trials.
- Errors smoother than about 16–32 µm are corrected; independent, per-waveguide errors are not, which is why static trim is a design requirement.
- More heaters help only at short ℓ.
- The slot design's ceiling is its own steering residual (0.96).

**fig10_lens.** The thermal lens: one wire focuses, five wires steer and focus. All layouts are on the undercut membrane, with heaters allowed out to the cavity walls (far kernel). Method as in `run_lens.py`.
- (a) Heater layouts (gold) under the 33 emitters. The line styles match (b)–(d).
- (b) Phase across the array for a unit lens, in which the centre leads the edges by 1 rad. In the vacuum limit the membrane obeys k·t·T″ = −q. Heat spread evenly therefore gives a parabola (T″ constant), and a point heater gives a tent (T″ = 0 away from it). The 64 µm heater sits on the target. The centre heater's tent limits it to a lens strength of a ≤ 5.4 rad.
- (c) Power to hold a lens against its focal length, f = 2076 rad·µm / a.
  - Grey: one isolated heater per emitter, phases wrapped (a Fresnel lens), at three P_π: the standard stack (19.7 mW), our undercut membrane (1.35 mW), and the best published suspended heater (Fang 2011, 0.49 mW). Wrapping stops their cost growing above a ≈ 2π, while ours grows in proportion to a.
  - For weak lenses (f ≳ 0.5 mm) ours is 8–16× cheaper than isolated heaters with the same P_π, and 3–6× cheaper than Fang's. The crossover with Fang is at f ≈ 0.14–0.19 mm.
  - Our power is the cheapest drive that keeps every emitter within 2% of the lens. The centre-heater and dual-use layouts can't meet 2%, so theirs is the least-squares drive.
  - Shaded: f < 0.1 mm, where the paraxial lens is only approximate.
- (d) The strongest lens (Strehl ≥ 0.8) that can sit on top of a beam steered to each angle. With 5 wires, the dual-use layout keeps f = 90–124 µm over the full 0–22.8°. The best 4-heater steering layout manages only weak lenses. The best lens layouts have no heater at the array ends, so they can't steer.

**fig11_faults.** Fault tolerance, beam at 10°. Source: `run_faults.py` (`faults.csv`).
- (a) Steering range toward +θ after one heater dies and the remaining heaters are refitted.
  - Losing the edge heater, the one that drives the ramp, leaves under 1° in every design. Losing an interior heater costs range, but the design still steers.
  - A spare edge heater (one more wire) restores nearly all of it: 22.8° in vacuum (spare at 32 µm) and 22.0° in the slot design (spare at 34 µm).
  - The dual-use layout's −33 and −16 µm heaters serve −θ, so losing them leaves +θ untouched.
  - The undercut design behaves like the slot design (`faults.md`).
- (b) Pointing error at 10° from slow disturbances, against a tenth of the 1.2° beam width (dashed).
  - Colours: dark, all four of our designs; grey diamonds, a conventional array with a heater on every emitter; blue, the slot design in air; orange, vacuum.
  - Our angle is carried by the amount of heat, so a 1% gain error moves the beam by 1% of 10°, i.e. 0.1°. That is 20–90× more than the conventional array, whose angle is carried by its wrap pattern.
  - Ambient +10 K changes the gain by (b − 2a)·ΔT, with a = 3×10⁻⁴/K from the conductivity and b the unverified growth of dn/dT. The result is 0.06° at b = 0, 0.01° at 0.7×10⁻³/K and 0.08° at 1.4×10⁻³/K; it passes through zero near b = 6×10⁻⁴/K.
  - A vacuum package must stay below about 25 Pa.
  - Pressure and gas rows: 1D membrane model calibrated to COMSOL (est.).
  - Radiation, left out of the COMSOL vacuum run, moves the beam 0.04° and calibrates out.
