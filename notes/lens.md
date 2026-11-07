# The thermal lens: results (6 Nov 2026)

Code: `analysis/opa/lens.py`, `analysis/run_lens.py`. Results: `analysis/results/lens_layouts.csv`, `lens_comparison.csv`. Kernels: today's undercut (`kernel_undercut.csv`), R1's far kernel (heaters to ±40 µm), R4a's step responses. The slot kernels are one-sided, so they can't make a symmetric dome and are left out.

## Definitions

A 1D (cylindrical) lens of focal length f needs φ(x) = c − k₀x²/(2f) across the array, in the paraxial limit.
- **Lens strength a:** the centre leads the edges by a = k₀X²/(2f), with X = 32 µm. So a·f = 2076 rad·µm.
  - a = 1 rad is f = 2.1 mm.
  - a = 2.08 rad is f = 1 mm.
  - a = 2π is f = 0.33 mm.
- **Sign:** dn/dT > 0, so heating the centre gives a converging lens.
- **Spot size:** the line focus is about λf/D wide (D = 64 µm). That is about 24 µm at f = 1 mm and 8 µm at 0.33 mm.
- **Quality:** the Strehl ratio at the focus, ≥ 0.8.
- **Power:** the cheapest drive that keeps every emitter within 2% of the lens profile (an LP). The least-squares drive is also listed in the CSV.

**Why the power needs care.** On a shared membrane many drives make nearly the same profile, and the least-squares one isn't the cheapest. With one heater per emitter, for example, it costs 1.13 mW/rad against 0.93 for the cheapest.

## Physics

The membrane obeys k·t·T″ − h·T = −q(x), with q the heater power per unit area.

- **A uniform heater makes a parabola.** A parabola needs T″ = constant, which means heat spread evenly across the array. In vacuum (h = 0), one heater the width of the array therefore makes an exact parabola:
  - ΔT(0) − ΔT(X) = q·X²/(2kt);
  - so P/a = 2·k·t·W/(814·X²) = **0.89 mW/rad** for W = 64 µm.
  - Checked: code against the analytic value, to 4 digits.
- **In air** the loss term bends the profile slightly. COMSOL's kernel gives 1.16 mW/rad, and it is still effectively perfect up to a ≥ 60 rad.
- **A single centre heater makes a tent.** T″ = 0 except at the source, so it only resembles a parabola: 0.56 mW/rad, usable to a = 5.4 rad (f ≥ 0.38 mm).
- **Each extra heater gets closer to a uniform source:**

  | Heaters | Usable to |
  |---|---|
  | 3 | 23 rad |
  | 4 | 38 rad |
  | 5 or more | not limited by the heat profile; only the paraxial approximation limits it |

- **Scaling (vacuum-ideal, est.).** Per unit optical power (1/f), the wide heater costs ∝ k·t·X. Isolated heaters cost ∝ N·X² ∝ X³ at fixed pitch, so the advantage grows as X² with aperture, for weak (unwrapped) lenses.

## Our layouts (undercut; the far kernel gives the same, since no layout wants heaters beyond ±22 µm)

| Design | Wires | Usable to | mW/rad | Power at f = 1 mm |
|---|---|---|---|---|
| One centre heater | 1 | 5.4 rad (f ≥ 381 µm) | 0.56 (least squares) | 1.2 mW |
| K = 3 at 0, ±16 | 3 | 23 rad (f ≥ 90 µm) | 0.86 (least squares) | 1.8 mW |
| K = 5 at 0, ±11, ±22 | 5 | ≥ 60 rad | 0.91 | 1.9 mW |
| One uniform 64 µm heater | 1 | ≥ 60 rad | 1.16 | 2.4 mW |
| Best fixed-shape single heater (graded width) | 1 | ≥ 60 rad | 1.13 | 2.3 mW |
| One heater per emitter (shared membrane) | 33 | ≥ 60 rad | 0.93 | 1.9 mW |
| One heater per emitter, standard stack | 33 | ≥ 60 rad | 35.7 | 74 mW |

The cheapest-drive LP can't hold K = 1 and K = 3 within 2%, so their numbers are least-squares drives.

## Against conventional heaters (mW)

| Lens | a = 1 (f 2.1 mm) | 2.08 (1 mm) | π (0.66 mm) | 2π (0.33 mm) | 4π (0.17 mm) | 24 (86 µm) |
|---|---|---|---|---|---|---|
| Isolated, undercut P_π 1.35 mW, unwrapped | 9.2 | 19.1 | 28.8 | 57.7 | 115 | 220 |
| Isolated, same, wrapped (Fresnel) | 9.2 | 19.1 | 28.8 | 34.0 | 34.5 | 35.4 |
| Isolated, best published P_π 0.49 mW [Fang 2011], wrapped | 3.3 | 6.9 | 10.4 | 12.3 | 12.5 | 12.8 |
| Isolated, standard stack, wrapped | 134 | 277 | 420 | 495 | 502 | 515 |
| **Ours, one centre heater** | **0.6** | **1.2** | **1.8** | — | — | — |
| **Ours, K = 5** | **0.9** | **1.9** | **2.9** | **5.7** | **11.5** | **21.9** |
| **Ours, one uniform heater** | **1.2** | **2.4** | **3.6** | **7.3** | **14.6** | **27.9** |

Reading:
- **Weak lenses (f ≳ 0.5 mm):** we are 8–16× cheaper than isolated heaters with the same P_π, and 3–6× cheaper than the best published isolated suspended heater. The range runs from the one-wire uniform heater to the one centre heater.
- **Wrapping caps the conventional cost.** A Fresnel lens's power stops growing above about a = 2π. Our power grows linearly with a. The crossovers:
  - against the same P_π: about a = 38 rad (f ≈ 55 µm), outside the paraxial range;
  - against the best published heater: about a = 14 rad (f ≈ 150 µm).
- **Correction to the prep note.** It said our advantage held only for a ≲ 2π. With optimised layouts it holds up to f ≈ 150 µm even against the best published heater.
- **Wires:** 1–5 against 33. The single-heater versions are fixed-shape lenses with tunable strength, the same as a Doerr-style shaped heater, but on a shared membrane at about 1 mW/rad. Doerr's silica AWG lens used about 3–7 W [Doerr 2002, 2006, V], on a different platform.
- **Core temperature is low.** A lens of a = 2.08 rad (f = 1 mm) needs only 2.5 K of centre-to-edge difference for a 1 mm heater (a/(814·L_h)).

## Steering and focusing at once (far kernel)

Largest usable lens strength at each steering angle, Strehl ≥ 0.8:

| Layout | 0° | 10° | 16° | 20° |
|---|---|---|---|---|
| Lens-best K = 5 (0, ±11, ±22) | ≥ 60 rad | can't steer | can't steer | can't steer |
| Steering-best K = 4 (±32.7, −6.7, 9.4) | 9.6 | 7.4 | 4.8 | 2.3 |
| **Dual-use K = 5 (0, ±16, ±33)** | **23.0** | **21.2** | **19.4** | **18.0** |

The dual-use layout steers the full ±22.8° at 5.4 mW/deg and focuses down to f ≈ 90–115 µm at any angle, with 5 wires. A conventional array needs 33 wires for either job. The lens-best layouts can't steer because they have no heater at the array's ends.

## Speed (R4a step responses, plain step drive)

A lens of f = 0.66 mm reaches Strehl 0.8 in **0.44–0.50 ms** (K = 1, 3, 5). That is about the same as a 1° steering correction on the same membrane (0.36 ms), and 4× faster than a full 0 → 10° slew (2.0 ms).

The physics: how long the beam takes to settle depends mostly on how big the phase change is, not on its shape.
- Strehl 0.8 allows roughly 0.5 rad of error.
- A π lens therefore only needs to reach about 85% of its final profile.
- A 10° ramp spans about 45 rad across the array, so it needs about 99%, which takes several time constants.
- A 1° step spans about 4.5 rad, comparable to the lens, and settles about as fast.

## Caveats

- **Paraxial:** fine for f ≳ 100 µm (NA ≲ 0.3); approximate below that.
- **1D only:** focus along the grating axis would need heaters over the gratings (R10, on hold).
- **The isolated-heater baseline** uses our own membrane's P_π (1.35 mW). That is why the best published value (0.49 mW) is also shown.
- **Not computed, only cited:**
  - wavelength-tuned focusing [DeSantis 2024, A]: no heater power, but it needs a tunable laser and ties focus to wavelength;
  - Kuo 2022 [A], which uses 30 phase shifters.
