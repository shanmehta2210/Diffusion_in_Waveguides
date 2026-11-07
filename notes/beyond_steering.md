# Beyond steering: applications, geometries, simulations, fabrication errors (2 Nov 2026)

This follows `literature_review.md`.

**Tags used throughout:**
- **[V]**: the reference was checked against OpenAlex or arXiv, or read in full text.
- **[A]**: the reference was reported by a search agent, mostly from abstracts.
- **"est."**: my own calculation. Calculations marked **(surrogate)** use the 1D membrane model in §2.1.

The scratch scripts are in the session scratchpad:
- `membrane_and_errors.py`
- `ramp_membrane.py`
- `other_profiles.py`

## 0. Headline findings

1. **The heat kernel is naturally a lens.** In today's undercut stack, a single heater at the array centre makes a quadratic, centre-hot phase.
   - It reaches 5.5 rad from edge to centre at overlap ≥ 0.8, at 0.56 mW per rad. That is a 1D focal length of about 0.38 mm.
   - With 5 heaters it reaches 24 rad, a focal length of about 85 µm.
   - Isolated undercut heaters need about 9.2 mW per rad for the same lens, roughly 16× more (est.).
   - A cascade cannot make this profile at all.
2. **A "ramp membrane" geometry could beat every competitor on steering power** (surrogate). Cut a release slot just past an edge heater, so that heat can only leave by crossing the array to the far wall.
   - With today's 15 µm air gap, 2 heaters reach 17.8° one way at 0.81 mW/deg, about 8 mW at 10°.
   - With a healing length ≥ 150 µm (a deeper cavity or lower pressure), a **single heater** covers the full 22.8° at 0.55 mW/deg.
   - For comparison at 10°: today's design needs 49 mW, a cascade about 19 mW (est.), and isolated undercut heaters 40 mW.
   - Scaling law: P = k·t·(dφ/dx)/(k₀·dn_eff/dT) **does not depend on N or on heater length**. The hottest-core temperature grows ∝ N, and the response time grows ∝ span².
3. **Fabrication errors split cleanly into two kinds.**
   - Smooth errors are wafer-level width and thickness variation, about 80% of the width variance [Xing 2022, V]. They are corrected by 4–8 heaters: Strehl ≥ 0.98 once the correlation length is ≥ 16 µm.
   - Local random errors are uncorrectable and grow with the independent waveguide length: σ² ≈ L/L_coh, with L_coh = 4.17 mm for strip waveguides [Yang 2015, V].
   - Unwrapped ramps need long heaters to stay cool. That sets up a **temperature-versus-fabrication-error trade-off**, which is mild at small angles and severe at 16° (§4.2).

---

## 1. Applications beyond simple steering, ranked

| # | Application | Phase needed | Fit | Why / precedent |
|---|---|---|---|---|
| 1 | **Tunable 1D focus and divergence**: tweezers, retinal or optogenetic stimulation, FSO wide-beam acquisition then narrow-beam tracking | offset + quadratic | **Strong** | Heat already makes the shape (headline 1). Today's tunable-focus OPAs use 30 individual phase shifters [Kuo 2022, A] or wavelength tuning [DeSantis 2024, A]. Hosseini 2025 [V] uses meandered heaters. |
| 2 | **AWG / phased-array demux tuning**: a ramp shifts the centre wavelength, a quadratic sets the dispersion | ramp + quadratic | **Strong fit, low novelty** | Doerr's thermo-optic lens (silica, about 3–7 W, about 2 ms) [Doerr 2002, 2006, V]. An undercut SOI version would be new, at mW power (est.). The random per-arm errors still need trimming. |
| 3 | **Slow FSO tip-tilt and bias pointing; alignment drift** | ramp | Good | From the first note: 1–70 Hz loop bandwidth [Frost 2024, V]. |
| 4 | **2D steering at a fixed wavelength**: heat all gratings evenly to steer along the grating axis | ramp + common mode | Good add-on | Kim 2019 [V]: 0.016°/mW, 10°. The US Navy patent [A] claims a triangular heater plus a uniform grating heater. Our heaters sit before the gratings, so the free offset c does *not* do this; it needs one extra heater over the grating region. |
| 5 | **Gradient tuning of microring filter banks**: common shift plus channel spacing | offset + ramp | Partial | It corrects only the systematic part; random spacing errors (8 GHz σ) still need per-ring trim [Holzwarth 2008, V]. |
| 6 | **Tunable fan-out**, a sinusoidal phase grating | sinusoid | Weak | Only the standard stack can make ripples. Five heaters reach 2.8 rad amplitude at about 78 mW/rad, so an equal 3-way split costs about 110 mW (est.). The undercut is too smooth. |
| — | Adaptive-optics receivers, coherent beam combining, programmable meshes, structured illumination, ghost imaging | random or full-rank | Poor | These need independent per-element phases [A]. |

**Lens sign.** dn/dT > 0, so heating the centre adds optical path there, as a convex glass lens does. That makes it converging, which is the useful sign for tweezers and stimulation.

---

## 2. Geometries worth exploring

### 2.1 Membrane surrogate (validated)

Model: a 1D oxide membrane with T(±40 µm) = 0 and air loss below, which gives a healing length L.
- **Fitted to COMSOL:** L = 63 µm and k·t = 5.80×10⁻⁶ W/K. From first principles these are 58 µm and k_ox × 4.22 µm = 5.82×10⁻⁶.
- **Fit error:** 7–9% of K(0;0) at x_h = 0–32 µm.
- **Steering check:** for one heater at +32 µm it predicts 4.6° at 5.2 mW/deg; COMSOL gives 4.4° at 4.7 mW/deg.

It is good for ranking geometries, not for final numbers.

### 2.2 Candidates

| ID | Geometry | Expected effect | Source of the expectation |
|---|---|---|---|
| G1 | **Ramp membrane**: release slot next to a heater placed outside the emitter span (x_h ≈ 34–38 µm) | Linear temperature over the array; numbers in headline 2 | surrogate; no precedent found [A] |
| G2 | **Deeper cavity or reduced pressure** | Mainly *straighter* ramps; power improves only about 1.8×. Air conduction stays near its full value down to about 450 Pa for a 15 µm gap and needs ≲ 50 Pa to drop 10× (est., agent) | surrogate + gas kinetics |
| G3 | **Thinner membrane** or thinner cladding over the waveguides (the heater is no longer above the light in G1) | P ∝ t at the same span-limited speed; the healing length L ∝ √t shrinks, so the ramp curves more | est. |
| G4 | **Centre-heater lens** (today's stack, anchored at both walls) | Headline 1 | computed from the COMSOL kernel |
| G5 | **Heaters beside the core** (your run) + **a silicon slab as a heat spreader** in the standard stack | A continuous slab raises lateral conductance (k_Si = 148 against oxide 1.38). L_th ≈ √(k_Si·t_slab·t_box/k_ox) ≈ 4.4 µm for a 90 nm slab and 6.9 µm for 220 nm, against 3.3 µm today. That tunes L_th/d while keeping µs speed | est. |
| G6 | **Bidirectional ramp membrane**: (a) a built-in static tilt plus one heater sweeping 0 → 2θ_max, or (b) heaters at both ends on narrow tethers | (a) is simplest but doubles the peak ramp; (b) is a thermal-resistance design problem | est. |
| G7 | **2D array on one membrane** | In vacuum, T = ax + by is a valid steady state, but it needs linearly graded *boundary* temperatures, not point heaters. Bus-fed 2D arrays integrate temperature along each feed path | agent, physics |
| G8 | **Thermal shunts**: Si pillars or metal vias pinning temperature lines inside the cavity | Kernel shaping, or splitting one cavity into independent sub-arrays | [Sysak 2011, A], [Ning 2024, A] |

### 2.3 Design freedom from heaters sitting before the gratings

Because the heaters sit before the gratings, the **heated section and the emitting section can have different geometry**, joined by a short taper or fan-in.

- **The emitter pitch d sets the field of view** (±22.8° at 2 µm). Only the emitter index n matters for the phase ramp.
- **The heated-section pitch d_h sets the thermal regime** (L_th/d_h).
  - For a ramp membrane, P = k·t·(Δφ/d_h)/(k₀·dn/dT) ∝ 1/d_h, and the response time ∝ span² ∝ d_h². The pitch therefore trades power against speed directly.
  - The hottest-core temperature is set by 32Δφ and the heater length, not by d_h.
  - The floor is d_h ≥ 1.5 µm (Model B: under 1% coupling over 2 mm).
- **The heated-section waveguide width can be larger** (1–2 µm, which needs a larger d_h). This lowers sensitivity to width errors and lengthens L_coh, which eases the §4.2 trade-off.

### 2.4 Other places to put the heaters

| Layout | Verdict | Reason (est.) |
|---|---|---|
| **A. Transverse heaters before the gratings** (today) | **Keep as the main line** | Cold gratings; heater length (mm) and heated-section pitch and width are free |
| **B. Transverse heaters over the gratings** | **Reject for steering** | Three problems: (1) The grating emits within about 110 µm (1/e), so the heater can only be about 150 µm long. 10° then needs a core rise of 45/(814 × 150 µm) ≈ 370 K. (2) That rise tilts the hottest emitter by about 4° along the grating, against a beam about 0.8° wide there. (3) Phase builds up *while* the grating emits, so the transverse ramp grows along y. Each slice of the grating then steers differently, smearing the beam |
| **C. Two functional sections in series**: a slotted ramp membrane (steering) followed by a lens membrane anchored at both walls (focus) | **Promising** | Phases add along each waveguide, so each section's cavity can be designed for one profile: slot for the ramp, two anchors for the dome. Cost: more independent length (§4.2) and more footprint |
| **D. A + one or two heaters over the grating region** (longitudinal steering or focus) | **Good extension** | Uniform heating steers along the grating axis (+1.16° per 100 K, Model E); a gradient along y chirps the grating, giving a longitudinal focus. That gives 2D steering at a fixed wavelength with 4–6 wires. On the standard stack about 0.7 W per 1.16°; on an undercut about 13 mW per 1.16° (rough). Caveats: non-uniform temperature across x makes each emitter steer differently along y, and an undercut under the grating changes the substrate reflection, so the up/down split changes. Needs simulation 5b and Model E re-runs |
| **E. Multi-pass or folded heated section** | Low priority | Folding cuts the temperature at fixed footprint, but in a 33-guide parallel array folding is awkward. It does not help fabrication errors, because the error grows with path length, not footprint |

This also gives the brief's pitch sweep (analysis 4) a cleaner meaning: it sweeps d_h, the heated-section pitch, not the emitter pitch. The cost is a short taper or fan-in whose independent length adds local phase error, so keep it short.

---

## 3. COMSOL simulations, in priority order

| P | Simulation | Type | What it decides |
|---|---|---|---|
| 1 | **Ramp membrane (G1)**: slot + heater at 34–38 µm; cavity depth 15 / 50 / 100 µm; k_air → 0 (vacuum) | geometry + validation | Confirms or kills the 0.5–0.8 mW/deg prediction. It also covers the pending "heaters beyond 32 µm" run |
| 2 | **Heaters beside the core (your run)**, plus a Si-slab thickness sweep (0, 90, 150 nm) | geometry | Whether the standard stack can get broad spreading at µs speed |
| 3 | **Thermo-mechanical bow of the heated undercut membrane**: Solid Mechanics + thermal expansion + residual oxide stress | **new physics / error** | An out-of-plane bow δz adds about k₀·δz of phase, so 100 nm gives 0.4 rad, larger than the photoelastic term. A parabolic bow acts as defocus. Photoelastic Δn is about 3.5% of the thermo-optic term and only rescales dn/dT [Dickmann 2023, V]. Buckling from heat alone is unlikely (critical stress about 650 MPa, est.). Asymmetric stacks will bow |
| 4 | **Frequency-domain heat transfer** on the same 2D model, giving K(x; ω) | **new physics**, cheap | Each spatial mode decays at αq² + h/(ρct): the ramp is the slowest mode, with τ ≈ 1.7 ms in air and 3 ms in vacuum (est.). This is the basis for multi-heater pre-emphasis. Single-heater differential drive gave about 12× speed-up [Geis 2004, V]; no multi-heater dynamic control was found [A] |
| 5 | **Heat leaking from the heated section into the gratings** (merged into simulation 6). The heaters sit *before* the gratings | check | The heat leaking along the waveguides must leave the gratings cold. A ramp that reaches the gratings would tilt each emitter along the grating by up to 1.16° per 100 K, against a beam about 0.8° wide there. Note that the Si cores carry about as much heat along y as the oxide does: per core, k·A ≈ 1.6e-11 W·m/K, against 1.2e-11 for the oxide in one 2 µm pitch cell (est.). *Optional extension:* a separate heater along the grating gives a tunable focus in the other axis |
| 6 | **3D finite cavity and heater, including the grating region** | validation | 2D is exact for the integrated phase if no heat leaves through the cavity ends. End effects decay over about 38–51 µm, so keep the heated region about 120–150 µm from the ends (est.). Also check the temperature at the first grating period, and the extra power lost through the ends, roughly 2ℓ/L_h ≈ 8% for a 1 mm heater (est.) |
| 7 | **Sensitivities for the error model**: cavity width ±2 µm, membrane thickness ±5%, heater offset ±0.2 µm | error | Finite differences of G (about 6 runs). This feeds the linear error model in §4 |
| 8 | **Mode solver: dn_eff/dw for 0.5 versus 1–2 µm wide waveguides** | error, cheap | How much a wide heated section extends L_coh; Li 2022 used multimode widths for calibration-free arrays [A] |
| 9 | Thermal-wave drive: AC heaters with phase offsets | speculative | Heat goes as I², so this needs DC + AC drive. The wave decays e⁻¹ per radian, about 3.4× across the array at 100 Hz. Frame it as spatio-temporal control, not a travelling-wave scanner |

---

## 4. Fabrication errors: how ours compares, and how to model it

### 4.1 Error model

Everything is linear: φ = (G + δG)·P + c + e_path + e_env.

- **e_path** splits into a smooth part and a local part.
  - Smooth: wafer-level and die-level width and thickness variation, about 80% of the width variance; a small circuit sees "common-mode variation… with maybe an additional linear gradient" [Xing 2022, V].
  - Local random: σ_w ≈ 1.37 nm, uncorrelated between sites [Xing 2022, V]. Its accumulated phase follows σ² ≈ L/L_coh, with L_coh = 4.17 mm (strip) or 1.61 mm (rib) [Yang 2015, V]. That I_coh-to-variance mapping assumes Adar's definition, which is unconfirmed.
- **δG**: heater resistance and P_π spread (undercut CV about 1% on a 300 mm wafer [Parsons 2025, V]); cavity size; heater position.
- **e_env**: ambient offsets and on-chip gradients; these are smooth.

### 4.2 Results so far (est.)

**Correctable fraction** (σ = 1 rad, Gaussian correlation length ℓ; mean Strehl after least-squares correction):

| ℓ (µm) | none | undercut K = 4 | undercut K = 8 | standard K = 8 |
|---|---|---|---|---|
| 0 | 0.39 | 0.44 | 0.49 | 0.49 |
| 4 | 0.45 | 0.64 | 0.84 | 0.85 |
| 16 | 0.65 | 0.98 | 1.00 | 0.99 |
| 64 | 0.93 | 1.00 | 1.00 | 1.00 |

**Temperature versus fabrication-error trade-off.** The heater length is set by a core rise of 30 K, and the uncorrectable local error by S ≈ exp(−L/L_coh):

| Max angle | End-emitter phase | L for 30 K | L for 80 K | S_fab (strip, 30 K length) | S_fab (rib) |
|---|---|---|---|---|---|
| 2° | 9 rad | 0.37 mm | 0.14 mm | 0.91 | 0.79 |
| 5° | 23 rad | 0.93 mm | 0.35 mm | 0.80 | 0.56 |
| 10° | 45 rad | 1.84 mm | 0.69 mm | 0.64 | 0.32 |
| 16° | 72 rad | 2.93 mm | 1.10 mm | 0.50 | 0.16 |

The 30 K limit is conservative. The linear-model error at 80 K (3%) is deterministic and can be calibrated away; the real ceiling is ageing. Waveguides inside suspended heaters have been extrapolated to drift by more than 1 rad over years at 150–200 °C [Xue 2025, A].

**Mitigations:**
- strip rather than rib waveguides;
- wide waveguides in the heated section (simulation 8);
- short heaters, which favour the small-angle niche;
- static trim.

### 4.3 Verdict by error type

| Error | Conventional (N heaters) | Ours | Harm to ours |
|---|---|---|---|
| Local random path error | Fully corrected, at static power that rises to ~N·P_π when σ ≳ π (est.) | Not correctable; grows with heater length | **More** (needs trim or short and wide heated sections) |
| Smooth wafer and die variation | Corrected | Corrected by K heaters (ℓ ≥ 16 µm) | **Equal** |
| Pattern-density edge effect (~65 µm radius [Xing, V]) | Corrected | Partly, because it sits at the array edges | **Slightly more**; add dummy waveguides |
| Actuator gain (P_π, R) | N gains to calibrate | K gains; G1 is close to one gain | **Less** |
| Kernel shape (cavity, thickness, heater position) | Not applicable to isolated heaters | Today's K = 4 depends on positions; G1's linear profile should only rescale to first order | **To test** (simulation 7) |
| Ambient and on-chip gradients | Corrected | Corrected with the same heaters | **Equal** |
| High-temperature ageing | Wrapped phases, at most 2π per channel | Unwrapped, so hotter | **More** |
| Calibration effort | O(N) per angle: 3–10·N iterations, 40 s to 1 h per beam [Lian 2024, V; Kim 2023, A] | One O(N) trim; then O(K), and one unit fit covers every angle | **Less** |
| Trim technology | Optional | Required beyond small angles. Options: Ge implant (0.01π precision, about 700 °C local anneal, one-way) [Jayatilleka 2021, V]; Sb₂S₃ (5-bit, 77-day retention) [Chen 2023, V]. **None has been shown on an OPA**, and hot anneals at 2 µm pitch would also trim the neighbours | Risk |

**Can it be modelled efficiently? Yes.**
- The error model is linear algebra.
- The local part is analytic (exp(−L/L_coh)).
- The smooth part is a Monte Carlo over a covariance built from Xing's decomposition, which takes seconds.
- Only δG needs COMSOL (simulation 7, about 6 runs), plus 1–2 coupled runs for the membrane bow (simulation 3).

---

## References added in this note

**Verified [V]:**
- Doerr C.R. et al., OFC 2002. doi:10.1109/OFC.2002.1036755
- Doerr C.R. et al., IEEE PTL 18 (2006). doi:10.1109/LPT.2006.874676
- Kim S.-H. et al., Opt. Lett. 44, 411 (2019). doi:10.1364/OL.44.000411
- Guerber S. et al., Opt. Express 30, 35246 (2022). doi:10.1364/OE.464540
- Holzwarth C.W. et al., JVST B 26, 2164 (2008). doi:10.1116/1.3021389
- Yong Z. et al., Opt. Express 30, 7225 (2022; OpenAlex lists 2021). doi:10.1364/OE.448614
- Ceccarelli F. et al., arXiv:2406.16770 (2024)
- Pertsch T. et al., Appl. Phys. Lett. 80, 3247 (2002). doi:10.1063/1.1476720
- Dickmann J. et al., Sci. Rep. 13, 19455 (2023). doi:10.1038/s41598-023-46819-0
- Geis M.W. et al., IEEE PTL 16, 2514 (2004). doi:10.1109/LPT.2004.835194
- Xing Y., Dong J., Khan U., Bogaerts W., ACS Photonics 10, 928 (2023; OpenAlex lists 2022). doi:10.1021/acsphotonics.2c01194
- Yang Y. et al., "Phase coherence length in silicon photonic platform," Opt. Express 23, 16890 (2015). doi:10.1364/OE.23.016890
- Jayatilleka H. et al., J. Lightwave Technol. 39, 5083 (2021). doi:10.1109/JLT.2021.3079801
- Chen R. et al., Nat. Commun. 14, 3465 (2023). doi:10.1038/s41467-023-39180-3
- Lian D. et al., Nanophotonics 13, 29 (2024). doi:10.1515/nanoph-2023-0519
- Parsons R. et al., arXiv:2503.09009 (2025)

**Agent-reported [A]:**
- Kuo S.-I. et al., Opt. Lett. 47, 3676 (2022)
- DeSantis D.M. et al., Opt. Express 32, 44567 (2024)
- Notaros J. et al., J. Lightwave Technol. 36, 5912 (2018)
- Pruessner M.W. et al., US 9,989,831 (2018)
- Xue et al., Opt. Express 33, 43542 (2025)
- Kim J.-Y. et al., Sci. Rep. 13, 19929 (2023)
- Li et al., Opt. Express 30, 44029 (2022)
- Sysak M.N. et al., OFC 2011 OWZ6
- Ning N. et al., APL Photonics (2024)
- Martinez A. et al., Sci. Rep. 14, 20178 (2024)
