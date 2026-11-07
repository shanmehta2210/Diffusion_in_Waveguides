# COMSOL run sheet (6 Nov 2026)

**Deadlines:** presentation Thu 12 Nov; report 18 Nov.

**Before any run.** Save each run as a copy of an existing model, so the originals stay untouched. Put the CSV exports in the project root, next to the existing ones.

**Export format.** Use the same as `kernel_undercut.csv`:
- temperature rise T−T0 at core height y = 0.11 µm;
- x from −150 to 150 µm in 0.1 µm steps (3001 points, the existing export node);
- one column per heater position;
- headers left as COMSOL writes them (`T-T0 (K) @ ... x_h=...`).

The Python reader already parses that.

**Geometry coordinates** (from `modelCparams.txt`):

| Layer | y range |
|---|---|
| core | 0 to 0.22 µm |
| BOX | −2 to 0 µm |
| cladding | 0.22 to 2.22 µm |
| heater | 2.22 to 2.34 µm |
| undercut cavity | x = −40 to 40 µm, y = −17 to −2 µm |

---

## Needed for Thursday

### R1. Heaters far out, no slot (quick: 8 solves)

- **Model:** a copy of the model that wrote `kernel_undercut.csv`, with uc = 1.
- **Sweep:** x_h = 33, 34, …, 40.
  - x_h = 39 puts the heater edge on the cavity wall.
  - x_h = 40 straddles the wall.
- **Export:** `kernel_undercut_far.csv`
- **Tests:**
  - whether the ±16.7° saturation was caused by clipping heaters at 32 µm;
  - the surrogate's prediction that, without a slot, moving out gains little.

### R2. Slot "ramp membrane", today's 15 µm cavity (headline run)

**Model:** `modelC_slot.mph`, a copy with uc = 1. Build it as follows.

1. **Parameters:**
   - `x_slot 37[um]` "Slot centre"
   - `w_slot 2[um]` "Slot width"
2. **Geometry.** Add a Rectangle:
   - width `w_slot`;
   - height `t_box + h_wg + t_clad`;
   - corner at (`x_slot - w_slot/2`, `-t_box`).
   
   It runs from the cavity roof (y = −2) to the cladding top (y = 2.22), cutting the whole membrane. Use Form Union, so it becomes its own domain.
3. **Material:** assign the slot domain the same air as the cavity.
   - If the model has an air domain above the cladding, the slot opens into it.
   - If the top surface is a boundary condition instead, give the slot's top edge the same condition.
4. **Heater sweep:** `x_h = range(-32,1,34)`, 67 solves.
   - **Do not mirror.** The slot breaks the left–right symmetry, so negative positions must be simulated.
   - Keep x_h ≤ 34 so the 2 µm heater (33–35 µm) stays clear of the slot (36–38 µm). The last core edge is at 32.25 µm.
   - Keep the slot ≥ 2 µm wide: Parsons 2025 saw damage with 1 µm trenches.
5. **Also evaluate:** the maximum temperature in the heater domain for each x_h (a Global Evaluation, exported to a second small CSV). This is for reliability.
6. **Export:** `kernel_slot_D15.csv` (and `kernel_slot_D15_heaterTmax.csv`).

**Prediction to test** (membrane surrogate, still only an estimate): a heater at 34–35 µm plus one interior heater gives about 17.8° one way at about 0.8 mW/deg, against 4.9 mW/deg today.

### R4. Transients exported at every core (for response time and pre-emphasis)

- **Model:** a copy of `modelD_transient.mph`.
- **Why:** the system is linear, so the step response from each heater position to each core is the *entire* dynamic model. Pre-emphasis can then be optimised in Python.
- **Export points:** all 33 cores, x = −32, −30, …, 32 at y = 0.11. The old run had only 6 points.
- **Output times:** `10^{range(-7,0.1,-1.5)}`, i.e. 1e-7 s to 3.2e-2 s, 56 points.
  - The slotted ramp is slower, roughly 2–3 ms, and the old run reached only 98–99% by 3.16 ms.
- **Solver:**
  - relative tolerance **1e-4**;
  - time stepping "Intermediate" or "Strict".
  
  The old standard-stack run overshot by up to 3%, a numerical artefact.
- **Input:** the 1 W/m step as before, or a smoothed step with about 10 ns transition, which removes the start-up glitch.
- **Cases:**
  - **(a)** Undercut, no slot: x_h = 0, 4, …, 32 (9 runs). Mirroring is fine here.
  - **(b)** Slot model from R2: x_h = −32, −28, …, 32, plus 34 (18 runs).
- **Export:** `transient_cores_uc.csv` and `transient_cores_slot.csv`. Times and x_h will be in the headers.

### R3. Slot with a deeper cavity, and vacuum (if time allows before Thursday)

- **Copies of R2 with:**
  - D_uc = 50 µm;
  - D_uc = 100 µm. **Set t_sub = 200 µm**, otherwise the cavity reaches the bottom of the simulated substrate.
  - Vacuum limit: a material with k = 1e-6 W/(m·K) (ρ and Cp arbitrary) on the cavity, the slot, and any air above. Radiation is ignored, so this is an upper bound on how straight the ramp can get.
- **Sweep:** x_h = `range(-32,2,34)`, 34 solves each.
- **Export:** `kernel_slot_D50.csv`, `kernel_slot_D100.csv`, `kernel_slot_vac.csv`.
- **Prediction:** a single heater at 34–35 µm covers the full 22.8° at about 0.55 mW/deg.

---

## For the report (by 18 Nov)

### R5. Heaters beside the core (your design), plus a silicon slab

- Standard stack, and undercut if time.
- Si slab of 0, 90 and 150 nm, continuous across the array.
- Sweep heater positions as for `kernel_checks.csv`.
- **Export:** `kernel_side_<stack>_slab<t>.csv`
- **Prediction:** the heat half-width grows from 3.3 µm to about 4.4–6.9 µm at µs speed.

### R6. Thermo-mechanical bow of the heated membrane

**Physics.** Add Solid Mechanics with the Thermal Expansion multiphysics coupling.
- Plane strain, because the membrane is long along the waveguides.
- Geometric nonlinearity on.
- Bottom of the substrate fixed.

**Materials:**

| Material | E | ν | α |
|---|---|---|---|
| Si | 170 GPa | 0.28 | 2.6e-6 /K |
| SiO₂ | 70 GPa | 0.17 | 0.5e-6 /K |
| TiN | ~500 GPa | — | 9.4e-6 /K |

**Residual oxide stress.** Use Initial Stress on the BOX and cladding, and sweep it: 0, −100 and −300 MPa (compressive). Treat this as a parameter, since the real value depends on the process.

**Cases:**
- heater off (residual stress only);
- heater on, at a power giving about 30 K and about 80 K at the hottest core;
- for both the undercut and the slot geometry.

**Export:**
- vertical displacement v(x) at core height;
- optionally σxx, σyy and σzz in the cores, for the photoelastic estimate.

**What it decides.** A vertical shift δz adds a phase of about k₀·δz, so 100 nm gives 0.4 rad. A parabolic bow acts as defocus.

### R7. 3D validation, including the grating region

- Thermal model only.
- Finite cavity, about 400 µm long along the waveguides.
- Heater about 200 µm long, centred.
- Waveguides continue past the cavity end onto the standard stack, where the gratings start.
- Use mirror symmetry at the heater's mid-length to halve the model.
- **Export:**
  - ∫ΔT dy along each core, compared with 2D × heater length;
  - the temperature along a core from the cavity end into the grating region.
- **What it decides:**
  - whether the 2D model is valid (keep the heated region about 120–150 µm from the cavity ends);
  - whether the gratings stay cold.

### R8. Sensitivities for the error model (about 8 solves on the R2 model)

- Cavity width W_uc = 78 and 82 µm.
- t_box and t_clad ±5% (1.9 and 2.1 µm).
- Heater offset ±0.2 µm, for the chosen design's heaters.
- **Export:** `kernel_slot_sens.csv`

### Optional

- **R9. Frequency-domain heat transfer:** K(x; ω) from 10 Hz to 10 kHz. This cross-checks R4.
- **R10. Heater over the grating region.** On hold until you've looked into it.
