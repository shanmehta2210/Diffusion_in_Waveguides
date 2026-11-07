# Pre-emphasis: shaping the heater drives in time (7 Nov 2026)

Code: `analysis/opa/preemphasis.py`, `analysis/run_preemphasis.py`. Data: the R4 step responses (`transient_R4a/R4b/R4vac.csv`). Results:
- `analysis/results/preemphasis.csv`: sections A–D;
- `preemphasis_robust.csv`: E–F;
- `preemphasis_traces.csv`: beam path against time.

A full run takes over an hour, so sections can be run alone: `python run_preemphasis.py E F`.

## Bottom line

- **Cold start to 10°, with the chip allowed 3× its hold power for a moment:**

  | Design | Plain step | Pre-emphasis | Speed-up |
  |---|---|---|---|
  | Slot | 3.42 ms | 0.43 ms | 8× |
  | Vacuum | 8.54 ms | 1.43 ms | 6× |
  | Undercut | 1.94 ms | 0.16 ms | 12× |

  At 2× power the speed-ups are 4.9×, 4.2× and 8.6×.
- **Power × settling time at 10°** falls from 45–96 mW·ms to 6–8 mW·ms. Published heaters driven with plain steps score 0.5–5 mW·ms. The 10–100× gap shrinks to about 1.2–16×.
- **Independent waveforms per heater are what make it work.** One shared overdrive waveform gives only 3–4×, and heats the hottest core 1.5–2.3× above its final temperature. With independent waveforms that overshoot is 1.08× (slot), and capping it at 5% costs nothing.
- **Two heaters get most of it.** At m = 3, two heaters give 7.3×, three give 8.0× and all seventeen give 10.7×. The few-wire design keeps nearly all of the multi-heater speed-up.
- **The catch is model accuracy.** Pre-emphasis without feedback needs the thermal time scale right to 1–2%. If the chip's heat flow is 2% slower than modelled, the slot's 8× drops to 3×. A drive designed to work across a ±5% band keeps 6.8× (slot) and 5.6× (vacuum).
- **Small corrections gain most and don't care about the model.** A 1° correction from a held angle settles 16× (slot) and 10× (vacuum) faster. A ±10% time-scale error costs under 20% of that.
- **Falling to 0° can't be sped up**, because heaters can't cool and they are all already off at 0°. Falling to a lower angle that isn't zero speeds up 3–7×. If fast falling slews matter, define 0° at a small standing bias.

## Method

**Model.** The R4 step responses S_k(x_n, t) are the whole linear dynamics. For heater drives that are constant on each time interval, the phase is a sum of rectangular-pulse responses (`PulseModel`). Between output times the steps are interpolated in log t (PCHIP); after 31.6 ms they are flat.

**"Settled".** The Strehl toward the target stays ≥ 0.8 from time T on, the same test as R4's "arrive". For the optimisation we use a stricter rms test, which keeps it convex:

  σ(t) = rms of the phase error about its mean ≤ σ* = √(2(1 − √0.8)) = 0.4595 rad.

This guarantees Strehl ≥ 0.8, because |mean e^{iε}| ≥ mean cos(ε − ε̄) ≥ 1 − σ²/2. It is nearly tight: at σ* the Strehl is 0.80–0.81, and rms settling comes out only 1–2% later than Strehl settling. All quoted times are the actual Strehl times, recomputed from the optimised drive on a fine grid.

**Optimisation.** For a given T, the drive is a linear program:
- 0 ≤ P_k(t);
- Σ_k P_k(t) ≤ m·P10 at every moment (the chip's power limit, with P10 the hold power at 10°);
- the drive ends on the plain step's hold drive, so all methods end in the same steady state;
- σ(t) ≤ σ* from T on.

The rms limit is a cone constraint. We impose it with cutting planes: tangent cuts u·ε ≤ σ* are added wherever a solution breaks the limit, on a dense grid of check times running out to where every pulse has died away. Bisection on T finds the shortest settling time.

Among the drives that settle by T, the one returned departs least from the hold drive (Σ dt·|P − P_hold|). We switched to this objective after the least-energy version had two problems:
- left to choose its own hold drive, it ended at Strehl 0.81 and 9.65° instead of the step's 1.00 and 10°;
- with the hold drive fixed, it still lingered at the tolerance edge for several T.

**Drive grid.**
- 30 intervals before T (14 µs each for the slot at m = 3, so the driver needs about a 70 kHz update rate).
- 15 intervals from T to 3T, then hold.

## Checks

- **Single thermal mode, analytic.** The optimal overdrive settles in τ·ln(m/(m − 1 + ε)). The LP matches it to 0.3% at m = 1 and 0.1% at m = 2, 5 and 10.
- **Plain step.** The pulse model, sampled at the data times, gives R4's settling times exactly (3.507, 8.769, 1.960 ms). On the fine grid it gives 2.5% less (3.42, 8.54, 1.94 ms), because R4 interpolated between samples 26% apart. The fine-grid values are used here.
- **One heater, m = 1.** The LP gives back the plain step (8.63 against 8.59 ms, within the bisection tolerance), as it must.
- **Grid convergence.** For the slot at m = 3:
  - doubling the drive and check grids leaves T unchanged at 0.240 ms (an earlier per-heater-cap version of the problem; the final settings are 30/15);
  - halving them gives 0.242 ms;
  - doubling the free-drive window (3T to 6T) gives no change.
- **Monotonic in power.** Settling time falls monotonically with m in every design, once the check grid was made dense. An earlier sparse grid let error bumps hide between checks.
- **Final state.** Every shaped drive ends at the plain step's final Strehl: 0.966 (slot), 1.000 (vacuum).

## Results

### A. Cold start, 0 → 10°, against peak power (m = peak chip power / 10° hold power)

| Design | Drive | m = 1 | 1.5 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|---|---|
| Vacuum (1 heater) | overdrive | 8.54 | 2.87 | 2.01 | 1.43 | 1.09 | 0.88 ms |
| Slot (0, 12, 32) | shared waveform | 3.42 | 1.40 | 1.09 | 0.87 | 0.73 | 0.64 ms |
| Slot | independent | 2.78 | 1.06 | 0.69 | **0.43** | 0.26 | 0.17 ms |
| Undercut (−12, 8, 32) | shared waveform | 1.94 | 0.92 | 0.75 | 0.63 | 0.55 | 0.50 ms |
| Undercut | independent | 0.61 | 0.30 | 0.23 | **0.16** | 0.11 | 0.088 ms |

Hottest core during the transient, relative to its final temperature, at m = 3:

| Design | Shared waveform | Independent / overdrive |
|---|---|---|
| Vacuum | — | 1.25 |
| Slot | 1.54 | 1.08 |
| Undercut | 2.27 | 1.41 |

At m = 10 the shared waveform reaches 2.5× (slot) and 5.0× (undercut). For scale, the final hottest core is about 62 K above ambient for 1 mm heaters.

**Extra energy per switch** (beyond the hold power): 9–13 µJ for 0 → 10° at m = 2–3. That equals 0.2 ms (undercut) to 2 ms (vacuum) of hold power.

### B. Overshoot capped at 5% (m = 3)

| Design | Uncapped | Capped |
|---|---|---|
| Slot | 0.43 ms | 0.43 ms |
| Undercut | 0.16 ms | 0.17 ms |
| Vacuum | 1.43 ms | 1.58 ms |

So the multi-heater speed-up doesn't depend on overheating anything. (A cap at exactly 1.0 makes the LP degenerate, so 5% is used.)

### C. Slews from a held angle (m = 2), where heaters can also be turned down

| Slew | Vacuum | Slot | Undercut |
|---|---|---|---|
| 0 → 1° | 2.37 → 0.23 ms (10×) | 0.80 → 0.042 (19×) | 0.36 → 0.037 (10×) |
| 5 → 6° | 2.36 → 0.22 (11×) | 0.83 → 0.051 (16×) | 0.41 → 0.036 (11×) |
| 5 → 4° | 2.36 → 0.35 (7×) | 0.80 → 0.15 (5×) | 0.37 → 0.12 (3×) |
| 5 → 10° | 6.72 → 1.19 (6×) | 2.68 → 0.38 (7×) | 1.56 → 0.13 (12×) |
| 10 → 5° | 6.72 → 1.72 (4×) | 2.58 → 0.74 (3.5×) | 1.34 → 0.50 (2.7×) |
| 10 → 0° | 8.54 → 8.54 (1×) | 3.32 → 3.32 (1×) | 1.73 → 1.61 (1.1×) |

- **Small rising corrections gain most.** The chip's power limit (2 × P10) is 20× the 1° hold power, so there is a lot of headroom to overdrive with.
- **Falling slews gain 3–7×,** because heaters can drop to zero, below their final power.
- **Falling to exactly 0° gains nothing.** Every heater is already at zero. Speeding it up would mean heating the cold side to flatten the ramp, but the slot and vacuum designs have no heater there. The undercut's −12 µm heater buys only 10%.

**Extra energy for a 1° correction** is about 0.7 µJ. At a 70 Hz loop that is about 50 µW, negligible.

### D. More heaters (slot, cold start to 10°)

| Heaters | 10° hold | Plain step | m = 2 | m = 3 |
|---|---|---|---|---|
| 1 shared waveform over (0, 12, 32) | 13.8 mW | 3.42 ms | 1.09 ms | 0.87 ms |
| K = 2 (4, 32) | 13.8 mW | 3.60 ms | 0.77 ms (4.7×) | 0.49 ms (7.3×) |
| K = 3 (0, 12, 32) | 13.8 mW | 3.42 ms | 0.69 ms (4.9×) | 0.43 ms (8.0×) |
| K = 17 (every 4 µm, −32 to 32) | 14.5 mW | 3.36 ms | not run | 0.31 ms (10.7×) |

**Most of the gain comes from a second independent channel.**
- Two heaters already reach 7.3× at m = 3, against 3.9× for one shared waveform.
- A third heater adds about 10%.
- Seventeen heaters add only 1.4× more.

So the few-wire design keeps nearly all of the multi-heater speed-up. (The 17-heater case was run at m = 3 only, because it takes about 15 minutes per case.)

### E. Sensitivity: the drive designed on the nominal model, run on a chip whose thermal time scale is off

Settling time (ms) when the real time scale is a factor f × the model's:

| Case | f = 0.90 | 0.95 | 0.98 | 1.00 | 1.02 | 1.05 | 1.10 |
|---|---|---|---|---|---|---|---|
| Slot 0 → 10°, m = 3 | 0.73 | 0.57 | 0.43 | **0.43** | 1.11 | 1.41 | 1.74 |
| Vacuum 0 → 10°, m = 3 | 1.57 | 1.44 | 1.43 | **1.43** | 3.04 | 3.74 | 4.68 |
| Undercut 0 → 10°, m = 3 | 0.33 | 0.25 | 0.20 | **0.16** | 0.72 | 0.87 | 1.03 |
| Slot 5 → 6°, m = 2 | 0.049 | 0.050 | 0.050 | **0.051** | 0.052 | 0.053 | 0.060 |
| Vacuum 5 → 6°, m = 2 | 0.21 | 0.22 | 0.22 | **0.22** | 0.23 | 0.23 | 0.26 |

A plain step simply rescales: f × its nominal time.

**It is lopsided.**
- If the chip is faster than the model, the boost lands early and its excess decays quickly: little harm.
- If the chip is slower, the planned cut-back comes before the heat has arrived. The leftover error then sits in the slowest mode and takes most of a plain step's tail to decay.

### F. Robust design: one drive that must settle on models at f = 1 − b, 1 and 1 + b

Worst settling time over the band, with nine factors checked, including between the design points:

| Case | Band | Robust drive | Nominal drive | Plain step |
|---|---|---|---|---|
| Slot 0 → 10°, m = 3 | ±2% | 0.44 ms | 1.11 ms | 3.49 ms |
| | ±5% | **0.51** | 1.41 | 3.59 |
| | ±10% | 0.75 | 1.74 | 3.76 |
| Vacuum 0 → 10°, m = 3 | ±2% | 1.46 | 3.04 | 8.71 |
| | ±5% | **1.54** | 3.74 | 8.97 |
| | ±10% | 2.75 | 4.68 | 9.40 |
| Slot 5 → 6°, m = 2 | ±10% | 0.056 | 0.060 | 0.91 |
| Vacuum 5 → 6°, m = 2 | ±10% | 0.25 | 0.26 | 2.59 |

Designing for ±5% uncertainty costs only 8–18% against the nominal optimum, so it's the sensible default.

## Physics: why it works

**Why a plain step is slow at 10°.** The ramp spans about 45 rad (rms 13.4 rad about its mean). σ* = 0.46 rad means the profile must be within 3.4% of its final shape: several time constants of the slowest mode. A 1° step only has to get within 34%, so it settles about 3× faster in units of τ.

**One heater.** A single mode overdriven at m settles in τ·ln(m/(m − 1 + ε)) instead of τ·ln(1/ε). With ε = 0.034 that is an ideal 5.1× at m = 2, 8.7× at 3 and 33× at 10.

The vacuum design gets 4.2×, 6.0× and 9.7×. It falls short at high m because the membrane has several modes with different speeds:
- overdriving long enough for the slow mode overshoots the fast ones;
- the beam swings past the target (to 13.5° at m = 3) before settling;
- the hottest core reaches 1.25× (m = 3) to 2× (m = 10) its final temperature.

**One shared waveform for several heaters** has the same limit, only worse. It controls one combination of modes. Overdriving the steady pattern overheats the region near the strong end heater while the far side still lags, so the speed-up saturates at 3–5× however much power is allowed.

**Independent waveforms.** For the slot at m = 3, the solution is:
1. the middle heater (12 µm) at full power for the first 0.2 ms;
2. then the centre heater (0 µm);
3. the end heater (32 µm) boosting last, about 0.32–0.40 ms;
4. then all three settle onto the hold drive.

Heat that has far to travel goes in first. Heat that shows up quickly at the nearby cores goes in last. Each heater only has to fill its own stretch of membrane (about 20 µm), and diffusion time grows with distance squared. So the ramp can form roughly (64/20)² ≈ 10× faster than when it all comes from one end (est.), without overheating any core.

**Bias.**
- From a held angle, a heater can drop below its final power, so falling slews speed up as well.
- At 0° there is nothing left to turn down. The remedy is a standing bias, which the static trim can absorb as a fixed tilt. Its cost is static power.

## Against the literature

| | Plain step | With pre-emphasis (m = 3) | Robust ±5% | Published |
|---|---|---|---|---|
| 10° slew, power × settling time | 45–96 mW·ms | 6–8 mW·ms | 7.0 (slot), 8.1 (vacuum) | 0.5–5 mW·ms, plain step (comparison doc) |
| 1° correction, per-degree power × settling time | 1.2–2.0 mW·ms | 0.07–0.18 mW·ms | about the same | 0.4–1.7 mW·ms, best isolated heaters, plain step (est.) |
| Speed-up from pre-emphasis | — | 6–12× (10° slew); 11–16× (1° correction) | 5.6–6.8× | 3.7–12× at device level [Geis 2004, Atabaki 2010, Lin 2022, Gondo 2021] |

Reading:
- **Against published devices driven with plain steps:**
  - 10° slew: the 10–100× gap shrinks to about 1.2–16×, which is level with slow suspended heaters. The comparison doc's target for a 5–10× gain was 0.35–0.7 ms (5–10 mW·ms) for the slot; we reach 0.43 ms (5.9 mW·ms), or 0.51 ms (7.0 mW·ms) with the robust drive.
  - Small corrections: we move ahead.
- **Against devices that use pre-emphasis too:** their device-level gains (3.7–12×) are similar to ours, so the gap stays roughly where it was. An isolated heater is a single-input, nearly single-mode system, and single-mode overdrive gets close to its ideal.
- **What is new here** is the shaping. On a shared membrane, independent waveforms beat a shared waveform by 2–4×, with far less overshoot. No thermally steered OPA uses pre-emphasis yet (comparison doc).

## Caveats and next steps

- **Model accuracy.** Pre-emphasis without feedback assumes a known thermal time scale. On a real chip, the step responses would be measured in place (each heater stepped, the far field or a tap monitor recorded), giving the same data as R4. Then a robust drive (±5%) is designed. Closed-loop tip-tilt (the main application) corrects residual error anyway, and small corrections are insensitive (section E).
- **Heater rating.**
  - The slot drive puts the whole 3× budget (41.5 mW) into one heater at a time, 38× that heater's own hold power. The cores stay within 8% of their final temperature, but the heater element itself must take 41.5 mW for about 0.2 ms.
  - For the vacuum design the peak is 15.9 mW.
- **Linearity.** Superposition assumes temperature-independent material properties. The boost raises local temperatures briefly, so any temperature-dependent conductivity would show up there first.
  - **Suggested COMSOL confirmation:** feed the optimised slot waveform (m = 3) into the R4 time-dependent model as an interpolation function, and compare the core temperatures with the pulse model's prediction. This tests interpolation, superposition and any nonlinearity in one run.
- **2D model:** the heaters are treated as 1 mm long; end effects are ignored, as in R4.
- **Model scope.** Only a uniform time-scale error was tested (a heat capacity or diffusivity error). Errors in kernel shape are not covered.
