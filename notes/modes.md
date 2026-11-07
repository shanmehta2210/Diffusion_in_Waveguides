# Effective degrees of freedom: which phase patterns heat can make (7 Nov 2026)

The brief's analysis 5. Code: `analysis/opa/modes.py`, `analysis/run_modes.py` (seconds). Results:
- `analysis/results/modes.csv`: gains and usable-mode counts;
- `modes_shapes.csv`: the first four mode shapes;
- `modes_two_beam.csv`: the beam-splitting test.

## Bottom line

- **Wide spreading leaves only the few smoothest patterns.** Even with a heater under every emitter, the number of phase patterns that cost at most 10× the cheapest one is:

  | Kernel | L_th/d | Usable patterns (10×) | Usable patterns (100×) |
  |---|---|---|---|
  | Standard stack, d = 4 µm | 0.84 | all 32 | all 32 |
  | Standard stack, d = 2 µm | 1.7 | 15 | 32 |
  | Standard stack, d = 1.5 µm | 2.2 | 11 | 27 |
  | Undercut | 10 | 4 | 14 |
  | Slot | 14 | 3 | 10 |
  | Vacuum | 15 | 2 | 7 |

  L_th is the spreading half-width (half-width at half maximum of one heater's profile). The offset is excluded throughout, because it is free.
- **For the standard stack the count follows a simple law:** about 0.8·N·d/L_th. It matches a Gaussian kernel of the same half-width, (√(2 ln 10)·√(2 ln 2)/π)·N·d/L_th, to within one mode at every pitch from 1.5 to 4 µm.
- **On the membranes the cheap patterns are exactly the ones we use.**

  | Membrane | Cheapest patterns, in order |
  |---|---|
  | Slot | ramp, then dome (2.3× dearer) |
  | Vacuum | ramp, then dome (5.4× dearer) |
  | Undercut | dome and ramp at almost equal cost |

  The third pattern already costs 4–14× the cheapest. So the membranes are 2–4-pattern devices: built for steering and focusing, and nothing else.
- **Heaters only heat, and on the vacuum membrane that removes even the dome.** Every vacuum heater's profile rises monotonically toward the slot (its most negative step is 0.02% of the peak) and stays within 5% of flat on the slot side. So any non-negative combination is monotonic. The dome listed as its second mode needs negative heater power. In practice the vacuum design steers and does nothing else. The slot design leaks heat to air on the slot side, so it can make weak bumps; the undercut, with a wall on each side, makes both ramp and dome (as used for the lens).
- **Sharp patterns need the standard stack.** Splitting the beam into two at ±3°, ±5° and ±10°:
  - with a heater under every emitter, the standard stack reaches the ideal efficiency of the same 33-element array, or comes within 0.04 of it (smooth sine split: 0.62–0.76, at 165–670 mW);
  - the membranes reach at most 0.27–0.46 of 0.76 at ±3°, and essentially nothing (≤ 0.03) at ±5° or more.

  It isn't the number of heaters that fails, it is the kernel. A spreading length ten times the pitch can't make features a few pitches wide.
- **The answer to the abstract's question, "the arrangements best suited to sharp patterns":** short spreading (L_th/d ≲ 1–2) with a heater under every emitter. That is the conventional layout, with its crosstalk folded into the fit. No few-heater arrangement makes sharp patterns, so the report should say plainly that this is the trade.

## Method

- **Matrix.** G is N × N with a heater under every emitter, the most a stack can do. For the vacuum kernel the heaters sit on its 2 µm grid, which is the cores.
- **Offset.** The shared offset doesn't move the beam, so each column's mean over emitters is removed first.
- **Modes.** The singular value decomposition G = U·S·Vᵀ gives independent modes. Heater pattern v_i makes phase pattern u_i with gain s_i (rad per W), so the power for a given amount of pattern i scales as 1/s_i.
- **Usable.** A mode is usable if s_i ≥ 0.1·s₁ (at most 10× the power of the cheapest) or, looser, ≥ 0.01·s₁.
- **Kernels:**
  - standard stack at pitches 1.5–4 µm (it only depends on x − x_h, so any pitch can be sampled);
  - undercut, far (R1), slot and vacuum at 2 µm.
- **Sharp-pattern test:** splitting the beam into two at ±θ.
  - Binary 0/π grating: ideally 0.81 of the light in the two beams for an infinite array.
  - Smooth sine grating 1.84·cos(k₀x·sin θ): 0.68.
  - On the finite 33-element array the ideals are 0.77–0.89 (binary) and 0.62–0.76 (sine).
  - Each stack's 33 heaters are least-squares fitted to the target, and the light in the two target beams is read off the array factor.

## Checks

- **Exact.** On a ring (periodic distances) the matrix is circulant, and its singular values must equal the magnitude of the discrete Fourier transform. They do, to 3×10⁻¹⁵.
- **Approximate (Szegő).** On an open array the singular values approach the kernel's Fourier transform when the array is much wider than the kernel. For Gaussians of width 2, 4 and 8 µm on the 64 µm array, the first eight differ from the transform by at most 4%, 7% and 14% of the largest. The gap grows with width relative to the array, as expected from edge effects.
- **Analytic count.** A Gaussian kernel passes the modes with exp(−q²w²/2) ≥ τ, so it has about N·d·√(2 ln(1/τ))/(π·w) of them. In terms of the half-width L_th = w·√(2 ln 2), that is 0.80·N·d/L_th at τ = 0.1. Predicted 11.8, 15.8, 19.7, 23.7 and 31.6 for d = 1.5–4 µm; found 11, 15, 19, 24 and 32.
- **Consistency.** R1's far kernel gives identical modes to the original undercut kernel for heaters inside ±32 µm. It is the same structure, only swept further out.

## Results

### Gains of the first six modes, relative to the first

| Kernel | L_th | L_th/d | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|---|
| Standard, d = 1.5 µm | 3.3 µm | 2.23 | 1 | 0.88 | 0.73 | 0.61 | 0.49 | 0.40 |
| Standard, d = 2 µm | 3.3 | 1.67 | 1 | 0.91 | 0.80 | 0.72 | 0.61 | 0.53 |
| Standard, d = 4 µm | 3.3 | 0.84 | 1 | 0.94 | 0.89 | 0.85 | 0.80 | 0.76 |
| Undercut | 20.6 | 10.3 | 1 | 0.98 | 0.23 | 0.21 | 0.09 | 0.08 |
| Slot | 28.4 | 14.2 | 1 | 0.43 | 0.16 | 0.10 | 0.06 | 0.04 |
| Vacuum | 29.1 | 14.6 | 1 | 0.19 | 0.07 | 0.04 | 0.02 | 0.02 |

### Mode shapes and their cost

The cost is the size of the heater-power vector (its 2-norm) for 1 rad rms of the pattern. Use it to compare modes within one stack.

| Kernel | Mode 1 | Mode 2 | Mode 3 | Mode 4 |
|---|---|---|---|---|
| Undercut | dome, 0.5 mW | ramp, 0.5 mW | ripple, 2.3 mW | W-shape, 2.6 mW |
| Slot | ramp, 0.3 | dome, 0.6 | ripple, 1.7 | ripple, 2.8 |
| Vacuum | ramp, 0.1 | dome, 0.6 | S-curve, 1.5 | ripple, 2.8 |
| Standard, d = 2 µm | 9.4 | 10.3 | 11.7 | 13.1 |

On the standard stack the gains are nearly equal, so the modes mix and have no clean polynomial shapes. Every pattern costs about the same. That is the regime where heaters act as independent phase shifters.

### Splitting the beam in two: efficiency against the ideal for the same array

| Kernel | Split | Binary 0/π (ideal) | Sine (ideal) | Sine, heater power |
|---|---|---|---|---|
| Standard, d = 2 µm | ±3° | 0.88 (0.89) | 0.76 (0.76) | 165 mW |
| | ±5° | 0.77 (0.77) | 0.62 (0.62) | 319 mW |
| | ±10° | 0.85 (0.85) | 0.68 (0.72) | 669 mW |
| Undercut | ±3° | 0.24 | 0.27 | 150 mW |
| | ±5° / ±10° | 0.16 / 0.00 | 0.02 / 0.02 | — |
| Slot | ±3° | 0.36 | 0.46 | 112 mW |
| | ±5° / ±10° | 0.22 / 0.00 | 0.03 / 0.02 | — |
| Vacuum | ±3° to ±10° | ≤ 0.02 | ≤ 0.02 | — |

- Splitting at ±θ needs a phase grating of period λ/sin θ: 30, 18 and 9 µm for 3°, 5° and 10°. Those features are far narrower than the membranes' 20–30 µm spreading.
- The binary pattern's heater power isn't quoted: least squares chases its sharp edges through nearly unreachable modes and spends watts on them. The smooth sine is the fair test.
- The vacuum design can't make the pattern at all. With non-negative power its profiles are all monotonic (see Bottom line), and a grating isn't.

## Physics

- **Heat is a spatial low-pass filter, and the singular values are its transfer function.**
  - For a kernel that only depends on x − x_h, the modes are close to cosines across the array, and each gain is the kernel's Fourier transform at that spatial frequency.
  - A kernel of half-width L_th passes spatial frequencies up to about 1/L_th. An aperture N·d wide holds about N·d/L_th such features.
  - The constant 0.8 comes from the threshold (10×) and the kernel's shape.
- **Membranes are not shift-invariant.** The walls and the slot set the boundary conditions, so the cheap modes are the low-order solutions of the membrane itself: a ramp between a heater and a sink, and a dome between two sinks.
  - The slot design's one-sided boundary makes the ramp cheapest; the undercut's two walls make the dome as cheap as the ramp.
  - In vacuum nothing leaks downward, so the ramp is cheaper still. Every other pattern has to be built against the conduction path, which is why the gains fall off so steeply (0.19 for the dome).
- **This is the same filter that limits error correction** (`fabrication.md`): errors narrower than the spreading length are out of reach. The singular values say how far out of reach.

## What this means for the argument

- **The abstract's "penalises sharp patterns" now has a number.** The membranes give 2–4 usable patterns out of 32, against 15 for the standard stack at the same 2 µm pitch, and all 32 at 4 µm.
- **The design choice is a choice of L_th/d.** Large spreading makes steering and focusing cheap (the ramp and dome are the cheapest modes) and everything else expensive or impossible. Small spreading with a heater per emitter is the conventional design.
- **There is no middle ground that makes sharp patterns with few wires,** because the number of usable patterns can never exceed the number of heaters either.
- **What the membranes can do:** offset (free), ramp (steering), dome (focusing), and at higher cost an S-curve (cubic, a small aberration correction). That matches the applications the comparison doc already lists. Multi-beam, fan-out and holograms are out.

## Caveats

- **The usable count depends on the threshold.** The 10× and 100× numbers are given, and the 12 largest gains are in `modes.csv` for other choices.
- **The singular values ignore that heaters can only heat.** On the standard stack and the undercut, a free offset turns that into a power cost (a bias), not lost patterns. On the one-sided membranes it removes patterns outright: the vacuum design's non-negative combinations are all monotonic. The least-squares beam-splitting fit, which respects it, is the more direct test.
- **The pitch sweep uses the standard kernel only,** as the brief planned. A membrane's spreading scales with its own geometry, not the pitch, so its L_th/d at another pitch would need a new kernel.
