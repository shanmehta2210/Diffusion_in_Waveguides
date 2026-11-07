# Fabrication errors: what a few smooth heaters can correct (7 Nov 2026)

The brief's analysis 6, on the final designs. Code: `analysis/opa/fabrication.py`, `analysis/run_fabrication.py` (about 20 s). Results:
- `analysis/results/fab_white.csv`: independent errors;
- `fab_correlation.csv`: correlated errors;
- `fab_power.csv`: the cost of correcting;
- `fab_budget.csv`: path lengths.

## Bottom line

- **Independent per-waveguide errors can't be corrected by the heaters.**
  - At σ = 0.5 rad the corrected Strehl at 10° is 0.72–0.80, against Ruze's 0.79. The low end is the undercut, whose own steering residual already costs it 10%. At σ = 1 rad it is 0.38–0.43, against 0.39.
  - A least-squares fit onto K heater shapes plus the offset removes on average only (K + 1)/N of independent error variance: 2/33 for the vacuum design, 6/33 for the dual-use one.
- **Smooth errors are corrected** once their correlation length reaches about 16–32 µm, roughly a quarter to a half of the 64 µm aperture.

  | Design at 10° (σ = 1 rad) | ℓ = 16 µm | ℓ = 32 µm |
  |---|---|---|
  | Uncorrected | 0.64 | 0.78 |
  | Vacuum, 1 heater | 0.82 | 0.96 |
  | Slot, 3 heaters | 0.93 | its own steering limit, 0.96 |
  | Dual-use, 5 heaters | 0.98 | 0.99 |

  The corrected Strehl saturates at the design's own steering Strehl ("trim" in the tables). That is the corrected value with the error fully removed.
- **The correction is nearly free**: under 0.2 mW for the slot and vacuum designs. Per-element calibration also corrects everything, but at 0° it costs 10–17 P_π for σ = 0.5–1 rad of independent error.
- **So the design needs a one-time static trim beyond small angles,** and the trim must cover each waveguide's whole independent path, not just the heated section.
  - Without trim, keeping Strehl ≥ 0.9 limits each waveguide's independent path to 0.44 mm (strip waveguide).
  - A 10° heated section kept below 60 K is already 1.04 mm long. That alone limits an untrimmed array to Strehl 0.78.
  - The brief assumes trim. This analysis says the assumption isn't optional, so it should be stated as a requirement of the design.
- **The heaters' own job remains real.** They remove smooth errors (wafer gradients, pattern-density effects, on-chip temperature gradients) and drift after trim, at almost no power. Only K gains need calibrating.

## Error model

ε = local part + smooth part. Both are Gaussian, so the Strehl is averaged over Monte Carlo draws (400 per point; fixed seed 798).
- **Local:** independent between waveguides. Along an independent path of length L it accumulates as σ² ≈ L/L_coh, with L_coh = 4.17 mm for strip and 1.61 mm for rib waveguides [Yang 2015, V]. Mapping L_coh to a phase variance assumes Adar's definition, which is unconfirmed (see `beyond_steering.md` §4.1).
- **Smooth:** a Gaussian process across the array, with correlation exp(−Δx²/2ℓ²) [the smooth/local split follows Xing 2022, V]. ℓ = 0 is the independent case.

How each case is handled:
- **Uncorrected:** keep the nominal steering drive.
- **Corrected:** refit the heaters to target − ε, with non-negative least squares and a free offset (`steering.fit_powers`).
- **Trimmed:** ε removed.
- **Conventional:** one isolated heater per emitter, wrapped. It always corrects fully; its cost is the extra wrapped power.

Designs (all steer one way, toward +θ):
- vacuum, 1 heater at 34 µm;
- slot, 3 heaters at 2, 14 and 32 µm;
- undercut, 3 heaters at −8, 10 and 32 µm;
- dual-use, 5 heaters at 0, ±16 and ±33 µm (steers and focuses);
- for "how many heaters?": evenly spaced slot layouts with 5, 9 and 17 heaters.

## Checks

- **Ruze, finite N:** for independent Gaussian errors with nothing corrected, E[S] = e^{−σ²} + (1 − e^{−σ²})/N. Monte Carlo against analytic:

  | σ (rad) | Monte Carlo | Analytic |
  |---|---|---|
  | 0.1 | 0.9903 | 0.9904 |
  | 0.5 | 0.7849 | 0.7855 |
  | 1 | 0.3882 | 0.3870 |
  | 2 | 0.0485 | 0.0481 |
  | 3 | 0.0307 | 0.0304 |

- **Exact Gaussian expectation, with each design's own residual and the least-squares projection** (`expected_strehl`).
  - Paired ±ε samples cancel the linear noise term.
  - Vacuum, slot and undercut agree to 4 digits, corrected and uncorrected, at σ = 0.1 and 0.3 rad.
  - The dual-use corrected value falls short of the unconstrained formula (0.978 against 0.986 at σ = 0.1). The reason: its −33 µm heater carries zero power when steering one way, so it can't apply negative corrections there. That is the "heaters can only heat" constraint at work, not an error.
  - My first version of this check assumed a zero steering residual and missed by up to 15% for the undercut. It was replaced by the exact form.
- **A pure offset costs nothing:** with ℓ = 100 mm, the uncorrected Strehl equals the trimmed one (0.9596).

## Results

### Independent errors (the brief's σ sweep), 10°: uncorrected / corrected

| Design | Trimmed | σ = 0.1 | 0.3 | 0.5 | 1 | 2 | 3 rad |
|---|---|---|---|---|---|---|---|
| Vacuum, 1 heater | 1.00 | 0.99 / 0.99 | 0.92 / 0.92 | 0.79 / 0.79 | 0.38 / 0.39 | 0.05 / 0.05 | 0.03 / 0.03 |
| Slot, 3 heaters | 0.96 | 0.95 / 0.95 | 0.88 / 0.89 | 0.76 / 0.77 | 0.38 / 0.41 | 0.05 / 0.05 | 0.03 / 0.03 |
| Undercut, 3 heaters | 0.90 | 0.89 / 0.89 | 0.82 / 0.83 | 0.70 / 0.72 | 0.35 / 0.38 | 0.04 / 0.05 | 0.03 / 0.03 |
| Dual-use, 5 heaters | 0.99 | 0.98 / 0.98 | 0.90 / 0.91 | 0.78 / 0.80 | 0.39 / 0.43 | 0.05 / 0.06 | 0.03 / 0.03 |
| Ruze e^{−σ²} + (1 − e^{−σ²})/N | | 0.99 | 0.92 | 0.79 | 0.39 | 0.05 | 0.03 |

- At 0° and 5° the pattern is the same (`fab_white.csv`). The values sit higher where a design's own steering residual is smaller. At 0° the uncorrected values match Ruze to within 0.01, and correction adds at most 0.05.
- The heaters gain at most 0.05: corrected curves sit on the Ruze line.
- Per-element calibration would give Strehl 1 at every σ.

### Correlated errors, σ = 1 rad: Strehl after heater correction

| Design | θ | ℓ = 0 | 2 | 4 | 8 | 16 | 32 | 64 | 128 | 1000 µm | Trimmed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Uncorrected | | 0.38 | 0.41 | 0.46 | 0.52 | 0.64 | 0.78 | 0.92 | 0.98 | 1.00 | |
| Vacuum, 1 heater | 0° | 0.39 | 0.41 | 0.48 | 0.58 | 0.74 | 0.88 | 0.97 | 0.99 | 1.00 | 1.00 |
| | 10° | 0.40 | 0.43 | 0.52 | 0.62 | 0.82 | 0.96 | 1.00 | 1.00 | 1.00 | 1.00 |
| Slot, 3 heaters | 0° | 0.39 | 0.43 | 0.50 | 0.61 | 0.78 | 0.87 | 0.96 | 0.99 | 1.00 | 1.00 |
| | 10° | 0.41 | 0.49 | 0.60 | 0.77 | 0.93 | 0.96 | 0.96 | 0.96 | 0.96 | 0.96 |
| Undercut, 3 heaters | 0° | 0.41 | 0.45 | 0.53 | 0.65 | 0.82 | 0.91 | 0.97 | 0.99 | 1.00 | 1.00 |
| | 10° | 0.38 | 0.44 | 0.56 | 0.74 | 0.87 | 0.89 | 0.89 | 0.89 | 0.89 | 0.89 |
| Dual-use, 5 heaters | 0° | 0.43 | 0.54 | 0.70 | 0.88 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| | 10° | 0.43 | 0.53 | 0.68 | 0.89 | 0.98 | 0.99 | 0.99 | 0.99 | 0.99 | 0.99 |
| Slot, 5 even | 10° | 0.45 | 0.52 | 0.68 | 0.88 | 0.98 | 0.99 | 0.99 | 0.99 | 0.99 | 0.99 |
| Slot, 9 even | 10° | 0.46 | 0.59 | 0.75 | 0.93 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| Slot, 17 even | 10° | 0.47 | 0.59 | 0.79 | 0.95 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

Reading:
- **How short an error the layout can correct depends on how many heaters it has and how they are spread:**

  | Layout | Corrects errors down to about |
  |---|---|
  | 1 heater | ℓ ≈ 32 µm |
  | 3 one-sided heaters | ℓ ≈ 16 µm |
  | 5 two-sided heaters, or 5–9 evenly spaced | ℓ ≈ 8–16 µm |
  | 17 heaters | ℓ ≈ 8 µm |

- **No layout touches ℓ ≲ 4 µm.** That is two pitches, the white-noise limit. Heat simply can't make features that narrow.
- **At 0°, one-sided designs correct less** (slot 0.87 against 0.96 at ℓ = 32 µm). With no steering drive to turn down, a gradient of the wrong sign can't be removed. A standing bias tilt, which the trim can absorb, would fix that.

### Cost of correcting (mean extra power)

| Errors | Vacuum | Slot | Undercut | Dual-use | Conventional (wrapped) |
|---|---|---|---|---|---|
| Independent, σ = 0.5 rad, 0° | 0.01 mW | 0.07 mW | 0.23 mW | 28 mW (for 0.81 instead of 0.79) | +10 P_π |
| Independent, σ = 1 rad, 0° | 0.03 | 0.15 | 0.52 | 58 | +17 P_π |
| Smooth, ℓ = 64 µm, σ = 1 rad, 0° | 0.04 | 0.13 | 0.51 | 2.2 | +3.9 P_π |
| Any of these at 10° | ≤ 0.05 | ≤ 0.02 | ≤ 0.05 | 0.05–22 | about 0 (±3 P_π spread) |

- **The dual-use number is least squares overfitting noise.** It spends tens of mW chasing errors it can barely reduce. In practice the correction would be capped or regularised.
- **At 10° the conventional cost is about zero.** Wrapped phases are spread over 0–2π anyway, so random offsets add nothing on average. The "−2 P_π" in `fab_power.csv` is statistics: the best-offset search picks up a little from the noise itself.
- **At 0° calibration is not free for conventional arrays.** 10–17 P_π is 200–340 mW at the standard stack's 19.7 mW, but under 1 mW at the record 50 µW.

### Path-length budget without trim (independent errors, Ruze)

| Target Strehl | Independent path, strip | Rib |
|---|---|---|
| 0.95 | 0.21 mm | 0.08 mm |
| 0.90 | 0.44 mm | 0.17 mm |
| 0.80 | 0.93 mm | 0.36 mm |

The heated section alone (slot design; its hottest core carries 51 rad at 10°):

| Hottest-core limit at 10° | Heated length | Strehl from that section alone, untrimmed: strip / rib |
|---|---|---|
| 30 K | 2.09 mm | 0.61 / 0.27 |
| 60 K | 1.04 mm | 0.78 / 0.52 |
| 80 K | 0.78 mm | 0.83 / 0.62 |

At small angles the section can be short. Under 60 K it is 0.21 mm at 2° (Strehl 0.95) and 0.52 mm at 5° (0.88) (est.: the hottest phase scales with sin θ). But the splitter tree and routing add their own independent length. Typical splitter trees and routing run to millimetres per channel (est.), so working without trim is plausible only for very compact layouts at small angles.

## Physics

- **Why the heaters can't fix independent errors.** The correction lives in the space spanned by the K heater profiles plus the offset: m = K + 1 dimensions out of N = 33. Independent noise spreads its variance evenly over all N directions, so a projection removes only m/N of it. For the vacuum design that is 6%, and for 17 heaters it would be 55%. This is the same low-pass filtering that limits patterns: heat can't make features narrower than its spreading length.
- **Why smooth errors are easy.** An error correlated over ℓ ≫ pitch has most of its variance in a few broad shapes: offset, gradient, curvature. Those are exactly what the heaters make. The offset is free, every design makes a gradient, and two-sided layouts make curvature. Once ℓ exceeds roughly the heater spacing, nearly all the variance is reachable.
- **Why the one-sided designs are limited at 0°.** Heaters only heat. Removing a negative gradient needs the steering drive to back off, and at 0° there is nothing to back off. The dual-use layout shows the same thing at 10° on its idle left-hand heater.

## Against conventional, error by error

| Error | Conventional | Ours |
|---|---|---|
| Local random path error | Corrected fully. Free at large angles; 10–17 P_π at 0° for σ = 0.5–1 rad | **Not correctable by the heaters: needs trim**, or independent paths ≤ 0.44 mm for Strehl 0.9 |
| Smooth wafer and die variation, gradients (ℓ ≳ 16–32 µm) | Corrected | Corrected, at < 0.2 mW |
| On-chip temperature gradients, slow drift after trim | Corrected | Corrected (smooth) |
| Actuator gain spread | N gains | K gains |
| Calibration effort | Per element, typically per beam direction | One trim, then one linear fit covers every angle |
| Kernel-shape errors (cavity width, thickness, heater position) | Not applicable | **Untested**: needs R8 |

## What this means for the argument

**Decision (Shreyansh, 7 Nov): the report states static trim as a requirement of the design.**

- **Make the static trim an explicit requirement,** not a footnote.
  - The fair trade becomes: "a one-time trim plus K wires", against "N wires plus calibration per beam direction".
  - Trim technologies exist: Ge implant with local anneal [Jayatilleka 2021, V], and Sb₂S₃ phase-change [Chen 2023, V]. We found none shown on an OPA (an agent search on 2 Nov, so not exhaustive), and a hot anneal at 2 µm pitch would also trim the neighbours (`beyond_steering.md` §4.3). That is a real risk to name.
- **Without trim the claim shrinks to small angles in compact layouts:** total independent path ≤ 0.4 mm for Strehl 0.9.
- **The heaters' error-correction role is real but narrow.** They handle smooth and environmental errors and drift, cheaply.

## Caveats

- Gaussian errors with squared-exponential correlation are a modelling choice. Real errors mix scales; Xing 2022 also reports a pattern-density edge effect over about 65 µm, which sits in the correctable range.
- Least-squares correction isn't Strehl-optimal at large σ, but at σ ≳ 1 rad nothing short of per-element control helps anyway.
- Kernel-shape errors (δG) are not included; R8 would give them.
- How much wider waveguides in the heated section would lengthen L_coh is unquantified (mode-solver run C7, not done).
