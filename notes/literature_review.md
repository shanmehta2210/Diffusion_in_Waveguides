# Literature grounding: applications, crosstalk, heaters (2 Nov 2026)

What this note is for: before posing any design optimisation, check where a few-heater, heat-spreading OPA could be useful, how conventional OPAs handle crosstalk, and how large heaters really are.

**Verification key** (for every reference listed at the end):
- [V] title, authors and year checked against OpenAlex/arXiv, or the full text was read;
- [A] reported by a search agent from an abstract or excerpt, not independently checked.

Numbers marked "est." are my calculations with this project's constants.

---

## 1. Applications that need little steering

How our design translates for comparison (est.):
- **Bandwidth:** τ = 0.6–2 ms gives a bandwidth of about 80–265 Hz.
- **Beam width:** the 66 µm aperture gives a beam about 1.2° wide.
  - ±12° then holds about 20 resolvable spots, and ±16.7° about 28.
  - A telescope divides the angle but keeps the spot count, so for long links the spot count is the real range metric.

| Application | Range needed | Speed needed | Fit |
|---|---|---|---|
| Ground or terrestrial FSO tip-tilt (0.75–10.6 km) | 5–75 µrad stroke at the telescope [Frost 2024, V] | 1–70 Hz [Frost 2024, V] | **Good** |
| Satellite FSO fine pointing | about ±1° behind a telescope; OPA PAT demo 5° × 5° [A] | kHz loops, 46–60 µs rise [A] | Marginal: too slow for jitter, could be a bias or point-ahead stage |
| Optical tweezers | about ±6° (est. from >1 mm travel at 5 mm) [Sneh 2024, V] | holographic tweezers about 20 Hz [A] | **Good** |
| Retinal or neural optical stimulation | ±5.1° [Hosseini 2025, V] | not stated | **Good**, and it is direct prior art (see §2) |
| Chip-to-chip or fibre alignment drift | about 1° [weak] | s to min (unsourced) | Good on paper, weak evidence |
| Data-centre FSO, FireFly | about ±20° | 20 ms acceptable [A] | Marginal: range too large |
| Automotive LiDAR, foveated LiDAR, ion addressing, OCT, AR displays | wide, or µs speed | — | Poor |

**Strongest niches:**
- slow tip-tilt or bias pointing for FSO terminals, where size, weight and power are tight;
- focused-spot biophotonics (tweezers, stimulation), where ms response is enough and heat near tissue is a cost;
- static or slow alignment.

---

## 2. How conventional OPAs handle thermal crosstalk, and close prior art

### Crosstalk is avoided by spacing, trenches or calibration, not tolerated

**Heaters on a wide pitch, then a fan-in to the dense emitters:**

| Work | Phase-shifter spacing | Emitter pitch |
|---|---|---|
| Phare 2018 | 20 µm | 775 nm |
| Miller 2020 | 20 µm | 1.3 µm |
| Doylend 2011 | 100 µm (79 µm trench) | 3.5 µm |

[Phare, Miller: V; Doylend: A]

- **Rule of thumb:** phase shifters more than about 15 µm apart keep crosstalk under 10% [Voyant patent, A].
- **Undercut alone:** P_π 60 → 1.2 mW [van Niekerk 2022, A]; P_π 43 → 1.07 mW, 1.85 ms rise [Parsons 2025, V].
- **Calibration:** crosstalk-matrix compensation [Milanizadeh 2019, A]; genetic or lookup-table phase search [Phare 2018; Doylend 2011].

**Consequence for "pack the emitters closer": no.**
- Emitter pitch is already decoupled from heater pitch, either by fan-in or by putting the phase shifters on a bus.
- The pitch limit is optical coupling. Superlattices or alternating widths reach 775–780 nm pitch [Phare 2018; Song 2015, A].
- I found no 1D OPA where heater crosstalk sets the emitter pitch.

**Weak advantages that do survive:**
- No fan-in region or fan-out wiring. Over 90% of Phare 2018's chip area was fan-out wiring [A].
- K pads and DACs instead of N.

### Prior art for "few control signals"

| Work | Mechanism | Signals | Power scaling |
|---|---|---|---|
| Van Acoleyen 2009, Opt. Lett. 34, 1477 [V metadata] | one shaped Ti electrode over 16 gratings at 2 µm pitch; 2.3° steering | 1 | N² (Yaacobi 2014 says so, citing it) |
| Yaacobi 2014, Opt. Lett. 39, 4575 [V, full text read] | cascaded doped-Si shifters on the feed bus, 16 antennas at 2 µm pitch | 1 | N (equal power per antenna) |
| Poulton 2017, Opt. Lett. 42, 4091 [V metadata] | cascaded, 50 antennas | 3 groups, to correct fabrication phase noise | N |
| Hosseini 2024/2025 (Witzens group), SiN at 525 nm [V abstract] | meandered heaters, ±5.1° transverse | 4 (with focusing) | — |
| Li 2022; binary-tree OPAs [A] | log₂N electrode trees | 4 for 16 channels | — |

**Yaacobi 2014, as read in the full text:**
- One voltage on all cascaded shifters gives a constant Δφ between antennas.
- 51° range, about 18 mW per antenna, 100 kHz bandwidth.
- The heated paths are only about 1.5 µm long and reach about 500 °C (about 3.2% index change).
- They criticise Van Acoleyen's single electrode because "the power needed for steering increases as the square of the number of antennas".

**No paper found that deliberately uses lateral heat diffusion from K = 3–8 optimised heaters.** The novelty is real but narrow: Van Acoleyen and the US Navy triangular-heater patent [A] already build the ramp from heater *geometry*, and Hosseini from heater *layout*.

### Like-for-like power comparison

Steering power in mW, N = 33, toward +θ:

| Angle | Conventional, isolated, standard | Conventional, isolated, undercut (P_π 1.35 mW) | Cascade, standard (est.) | Cascade, undercut (est.) | **Ours, undercut, K = 4** |
|---|---|---|---|---|---|
| 1° | 468 | 32 | 28 | 2.0 | **4.9** |
| 5° | 566 | 39 | 142 | 9.7 | **24.7** |
| 10° | 580 | 40 | 282 | 19.4 | **49.1** |
| 16° | 597 | 41 | 448 | 30.8 | **78.0** |

How the cascade estimate is built (est.):
- Power = (N − 1)·Δφ·P_π/π, with the same P_π as our stacks.
- This assumes each stage is long enough to stay cool. Yaacobi's stages are about 1.5 µm and run at about 500 °C. At 10° the linear model would put a 1.5 µm stage about 1150 K above ambient. Long cascade stages need bus meanders between taps; how Poulton 2017 sized them is unchecked.

**Reading:**
- Against parallel-fed isolated undercut heaters, we win below about 8°: 6.5× at 1°, 1.6× at 5°.
- Against a cascade with the same P_π, we lose by about 2.5× at every angle.
- What we keep:
  - parallel feed (equal amplitudes from a splitter tree, no tap engineering, no loss or error accumulating along a bus);
  - long, cool heaters (core ΔT about 30 K at 16° with 3 mm heaters);
  - no fan-in;
  - 3–4 signals.

---

## 3. Heater dimensions

| Quantity | Typical | Extremes | Source |
|---|---|---|---|
| Width | 1–3 µm | 0.8 µm; ~20 µm meander (simulated) | Jacques 2019 [V]: optimum TiN 2.5 µm vs 2.0 µm minimum; others [A] |
| Length (straight) | 50–125 µm | ~1–1.5 mm | Harris 2014: 61.6 µm; Fang 2011: 1 mm suspended arm [V]; Hosseini meander 1.25–1.55 mm |
| Vertical gap to core | 1–2 µm | 0.75 µm lateral (doped) | ours (2 µm) is the conservative end [A] |
| P_π, standard stack | 20–30 mW | 12.7 mW | ours 19.7 mW, in line [A] |
| P_π, undercut | 0.5–1.5 mW | 0.05–0.14 mW folded | Fang 2011: 0.49 mW, 266 µs [V]; Parsons 2025: 1.07 mW [V] |
| τ, standard / undercut | 2–40 µs / 0.15–2 ms | — | ours 8–21 µs / 0.6–2 ms, in line |
| TiN failure | I ≈ 15 mA, P ≈ 136–179 mW (2–3 µm strips) | — | Limongi 2026 [A] |
| Core ΔT in practice | 60–80 K routine (est. from Harris and the AIM PDK) | ~500 °C (Yaacobi) | — |

**Placement beside the core** (doped Si in the slab): similar P_π, about 2.5× faster [Jacques 2019].

**Undercut size:** I found no precedent for one 80 µm × 15 µm cavity under 33 waveguides. Existing undercuts cover one phase shifter or a few folds. A cavity this wide would need a row of release holes, and suspended films carry fatigue risk.

**What this means for the model:**
- **Width:** the typical range (1–3 µm) is narrow, and P_π barely depends on it.
- **Wider heaters:** these change the spreading kernel, and so they matter for us. Because heat conduction is linear, a heater of width 2m µm is exactly m adjacent 2 µm heaters at equal power, assuming uniform current density. Its kernel can therefore be built from the existing CSVs with no new COMSOL run:
  - standard stack: any multiple of 2 µm;
  - undercut: positions on the 1 µm column grid.
- **Length:** in 2D only the temperature depends on length, ΔT = φ/(814·L). Length enters as a constraint (core ΔT ≤ 30 or 80 K gives a minimum L), not as a variable to optimise.
- **Heaters beside the core in the silicon slab:** this is a new geometry and needs COMSOL.

---

## References

Verification: V = checked against OpenAlex/arXiv or read in full text; A = from a search agent, not checked independently.

1. **[V]** Van Acoleyen K., Bogaerts W., Jágerská J., Le Thomas N., Houdré R., Baets R., Opt. Lett. 34, 1477 (2009). doi:10.1364/OL.34.001477
2. **[V]** Yaacobi A. et al., "Integrated phased array for wide-angle beam steering," Opt. Lett. 39, 4575 (2014). https://pmg.mit.edu/documents/ol-39-15-4575.pdf
3. **[V]** Poulton C.V. et al., "Coherent solid-state LIDAR with silicon photonic optical phased arrays," Opt. Lett. 42, 4091 (2017). doi:10.1364/OL.42.004091
4. **[V]** Hosseini P. et al., arXiv:2507.15018 (2025); ECIO 2024, arXiv:2406.14406
5. **[V]** Phare C.T., Shin M.C., Miller S.A., Stern B., Lipson M., arXiv:1802.04624 (2018)
6. **[V]** Miller S.A. et al., Optica 7, 3 (2020). doi:10.1364/OPTICA.7.000003
7. **[V]** Jacques M. et al., Opt. Express 27, 10456 (2019). doi:10.1364/OE.27.010456
8. **[V]** Fang Q. et al., IEEE PTL 23, 525 (2011). doi:10.1109/LPT.2011.2114336
9. **[V]** Parsons R. et al., arXiv:2503.09009 (2025)
10. **[V]** Frost A.M. et al., arXiv:2407.04745 (2024)
11. **[V]** Sneh T. et al., Nat. Commun. (2024). doi:10.1038/s41467-024-52273-x
12. **[A]** Doylend J.K. et al., Opt. Express 19, 21595 (2011)
13. **[A]** Milanizadeh M. et al., J. Lightwave Technol. 37, 1325 (2019)
14. **[A]** Song W. et al., Nat. Commun. 6, 7027 (2015). doi:10.1038/ncomms8027
15. **[A]** Harris N.C. et al., Opt. Express 22, 10487 (2014)
16. **[A]** van Niekerk M. et al., CLEO 2022, JTh3B.24
17. **[A]** Li et al., Opt. Express 30, 44029 (2022). doi:10.1364/OE.475350
18. **[A]** Pruessner, Rabinovich et al. (US Navy), US 9,989,831 (2018)
19. **[A]** Sacher W.D. et al., Opt. Lett. 47, 1073 (2022). doi:10.1364/OL.441609
20. **[A]** Limongi et al., J. Micromech. Microeng. 36, 085005 (2026)
