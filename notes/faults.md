# Fault tolerance (8 Nov 2026)

Code: `analysis/run_faults.py` (about 5 s), which uses the 1D membrane model with gas now shared in `analysis/opa/membrane1d.py` (moved out of `run_gas_loss.py`; `gas_loss.csv` reproduces byte for byte). Results: `analysis/results/faults.csv`. Written up as Q5 in the Q&A doc.

## Bottom line

- **Our pointing accuracy is our gain accuracy.** The angle is carried by the amount of heat, so anything that scales heat into phase moves the beam by the same fraction: 1% of drive is 0.1° at 10°, about a tenth of the 1.2° beam. A conventional array carries the angle in its wrap pattern and moves 0.005° for the same noise.
- **The edge heater is a single point of failure** in every design. A spare edge heater (one more wire) restores the full range: in vacuum a spare at 32 µm alone gives 22.8° at 0.533 mW/deg; in the slot design a spare at 34 µm gives 22.0° at 1.405 mW/deg.
- **A vacuum package needs a pressure spec.** Calibrated at 1 Pa with the drive held fixed, the beam falls 0.12° short (a tenth of its width) at about 26 Pa, and 0.45° short at 100 Pa. Q2's "1 mbar is enough" holds for power, not for drift.
- **Radiation, left out of the COMSOL vacuum run, is negligible:** +0.2–0.4% power.
- **Requirements this adds:** a spare edge heater; a constant-power driver stable to about 0.1%; a stable package (below about 25 Pa, held steady, or a heavy gas at 1 atm); a temperature sensor and one gain coefficient, or closed loop.

## Results (beam at 10°)

### One heater dies

| Design | Heater lost | Beam before recalibration | Refitted range |
|---|---|---|---|
| Vacuum, 1 heater | 34 µm | back to rest (0°) | none |
| Slot, 3 heaters | 2 / 14 / 32 µm | 9.54° (Strehl 0.57) / 9.42° (0.87) / 1.02° (0.33) | 10.2° / 16.2° / 0.9° |
| Undercut, 3 heaters | −8 / 10 / 32 µm | 10.15° (0.67) / 9.78° (0.66) / 0.06° (0.33) | 9.7° / 8.8° / 0.4° |
| Dual-use, 5 heaters | −33 / −16 / 0 / 16 / 33 µm | 10.00° (0.99) / 10.16° (0.94) / 10.00° (0.70) / 9.68° (0.82) / 0.16° (0.27) | 22.8° / 22.8° / 15.6° / 15.0° / 0.5° |
| Conventional, 33 heaters | any one | 10° (Strehl 0.93 mean, 0.85 worst) | full: the free offset can make the dead emitter's target zero |

### Drive noise, 1% rms (4000 trials)

| Design | Independent: jitter, Strehl | Shared (supply): jitter |
|---|---|---|
| Vacuum / slot / undercut / dual-use | 0.100° / 0.090° / 0.098° / 0.097°; Strehl unchanged | about 0.101° for all |
| Conventional | 0.0047°, Strehl 0.9988 | 0.0012° |

### Drift after calibration, drive fixed

| Case | Pointing change at 10° |
|---|---|
| Vacuum, 1 Pa → 10 / 30 / 100 / 300 / 1000 Pa | −0.044° / −0.14° / −0.45° / −1.14° / −2.5° |
| Vacuum → 10 kPa / 1 atm | −5.1° / −6.2° |
| Slot in air, 0.75 / 0.95 / 1.05 atm | +0.14° / +0.023° / −0.021° |
| Ambient +10 K, membrane and silicon (b = 0, 0.7, 1.4 ×10⁻³/K) | −0.060° / +0.010° / +0.080° |
| Ambient +10 K, gas only (slot in air) | −0.14° |
| Radiation, ε = 0.5 / 1, both faces | −0.021° / −0.042°; refitted 0.505 / 0.506 mW/deg against 0.504 |

## Method

- **Designs:** as in the fabrication and nonlinearity analyses: vacuum [34], slot [2, 14, 32], undercut [−8, 10, 32] on the original undercut kernel (range 15.9°), dual-use [−33, −16, 0, 16, 33] on the far kernel. All steer toward +θ.
- **Conventional array:** a heater under every emitter, isolated, phases wrapped into [0, 2π). Its error is added to the unwrapped ramp, which points exactly at 10°.
- **Beam:** pointing is the least-squares tilt of the phase across the array; Strehl is |mean exp(i r)|² of what is left after the tilt.
- **Dead heater:** before recalibration the others keep their 10° drive; after it, the remaining heaters are refitted with the usual non-negative fit.
- **Drive noise:** each heater's power is multiplied by (1 + σε), with ε independent per heater or shared, σ = 1%.
- **Drift:** the 1D membrane model (k·t = 5.8e-6 W/K; gas factors a_h = 1.29, a_s = 1.55 fitted to R2 and R3b). The drive is fitted at the calibration condition, then held while pressure or gas conductivity changes.
- **Ambient temperature:** conductivity k₀(1 + 2aΔT) with a = 3.0e-4/K (Model F) and dn/dT growth b, so the phase per watt changes by (b − 2a) per K. For the gas, k_air ∝ T^0.8 and the mean free path ∝ T at fixed pressure.
- **Radiation:** both faces at h = 4εσT³ each, added as a loss under the membrane in the vacuum model.

## Caveats

- **The drift numbers come from the 1D model (est.).** It is calibrated to COMSOL within 3.5% in air, but the pressure dependence itself has not been simulated in COMSOL.
- **b is unverified** (see `nonlinearity.md`), so the ambient term's sign is open.
- **The 1D model's slot design has Strehl 0.93 at 10°,** against COMSOL's 0.96. Only the change in pointing is used from it.
- **Channel failure rates aren't modelled.** Fewer channels means fewer failures, but each matters more.
- **Heater temperatures** are from `nonlinearity.md`: 65 K at 10°, about 150 K at 22.8° (scaled), and pre-emphasis bursts bounded by 202 K.
