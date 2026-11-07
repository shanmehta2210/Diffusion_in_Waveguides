"""Pre-emphasis: heater drives that settle the beam faster than a plain step.

A plain step switches each heater straight to its final power. The phase then
creeps toward the target at the pace of heat diffusion, and the slowest thermal
mode sets the settling time. Pre-emphasis overdrives first and backs off later.
With several heaters it can also reshape the temperature profile on the way,
and from a held angle it can cut a heater's power below its final value (which
a cold start cannot: heaters only heat).

Model (opa/transient.py): the step responses S_k(x_n, t) are the whole linear
dynamics. A drive that is constant on each interval [t_j, t_{j+1}), with
P_k(t) = P_kj, gives the phase
    φ_n(t) = φ_free,n(t) + 814·Σ_kj P_kj·[S_k(x_n, t - t_j) - S_k(x_n, t - t_{j+1})],
one rectangular-pulse response per interval and heater. The last interval runs
to infinity, so its powers are the final hold drive. φ_free is what the
starting state does on its own once its drive is removed at t = 0: zero from
cold, and φ_a - 814·Σ_k P_a,k·S_k(x_n, t) from a beam held steady with drive P_a.

Beam quality: the rms phase error about its mean (the offset is free),
    σ(t) = rms_n(ε_n(t) - mean ε),   ε_n = φ_n - r_n,
must stay ≤ σ* = √(2(1 - √0.8)) = 0.4595 rad from time T on. That guarantees
Strehl ≥ 0.8, because |mean e^{iε}| ≥ mean cos(ε_n - mean ε) ≥ 1 - σ²/2, and it
is nearly tight (a Strehl of 0.80-0.81 at σ*). σ(t) ≤ σ* is a convex (second-
order cone) constraint on the drive; we impose it with cutting planes on a
linear program: start from cuts along a few smooth error shapes, solve, add a
tangent cut u·ε ≤ σ* (u the unit direction of the error) wherever the solution
breaks the limit, repeat. Every cut is a valid outer bound, so "infeasible" is
final, and the loop stops once all errors are within limit. The shortest T is
found by bisection. verify() recomputes the actual Strehl settling time on a
fine grid; that is the number we quote.

Limits on the drive: P_kj ≥ 0 (heaters can only heat) and Σ_k P_kj ≤ P_cap at
every moment (the chip's power limit, so one heater and three are compared on
the same budget). The drive can be made to end on a given hold drive, normally
the plain step's, so both end in the same steady state; among the drives that
settle by T it then returns the one that departs least from that hold drive.
An optional cap on every core's phase limits the temperature overshoot, and the
drive can be required to work on several models at once (robust design).
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import linprog

from .constants import PHASE_COEFF
from .steering import STREHL_MIN
from .transient import settle_time

SIGMA_STREHL = float(np.sqrt(2 * (1 - np.sqrt(STREHL_MIN))))  # 0.4595 rad: rms error that guarantees Strehl 0.8


def strehl(phase_error):
    """|mean exp(iε)|² over emitters, for each row of phase errors."""
    return np.abs(np.exp(1j * np.atleast_2d(phase_error)).mean(axis=1)) ** 2


def rms_error(phase_error):
    """rms over emitters of the phase error about its mean (the shared offset is free), per row."""
    e = np.atleast_2d(phase_error)
    return np.sqrt(np.mean((e - e.mean(axis=1, keepdims=True)) ** 2, axis=1))


class PulseModel:
    """Phase response of a chosen set of heaters to piecewise-constant drives.

    sr is a StepResponses (or anything with times, x_h and S[i_t, k, n]).
    Between output times the step responses are interpolated (PCHIP in log t);
    after the last output (31.6 ms in R4) they are taken as flat.
    """

    def __init__(self, sr, heaters):
        self.x_h = np.asarray(heaters, float)
        cols = [int(np.flatnonzero(np.isclose(sr.x_h, xh))[0]) for xh in self.x_h]
        S = sr.S[:, cols, :]                      # (n_t, K, N), K per W/m
        self.t_data = np.asarray(sr.times)
        self.S_final = S[-1]                      # (K, N)
        self._interp = PchipInterpolator(np.log(self.t_data), S, axis=0)
        self.K, self.N = S.shape[1], S.shape[2]
        self.G = PHASE_COEFF * self.S_final.T     # steady kernel (N, K), rad per W

    def step(self, tau):
        """S_k(x_n, τ) for an array of delays τ (s): shape τ.shape + (K, N). Zero for τ ≤ 0."""
        tau = np.asarray(tau, float)
        out = np.zeros(tau.shape + (self.K, self.N))
        t0, t1 = self.t_data[0], self.t_data[-1]
        mid = (tau >= t0) & (tau < t1)
        out[mid] = self._interp(np.log(tau[mid]))
        out[tau >= t1] = self.S_final
        early = (tau > 0) & (tau < t0)            # before the first output: straight line from 0
        out[early] = self.S_first[None] * (tau[early] / t0)[:, None, None]
        return out

    @property
    def S_first(self):
        return self._interp(np.log(self.t_data[0]))

    def pulse_matrix(self, t_check, t_switch):
        """Phase (rad) at each check time per W in each (heater, interval): shape (n_check, N, K, J).

        Interval j is [t_switch[j], t_switch[j+1]); the last one never ends.
        """
        edges = np.append(t_switch, np.inf)
        on = self.step(t_check[:, None] - edges[None, :-1])     # (n_check, J, K, N)
        off = self.step(t_check[:, None] - edges[None, 1:])
        return PHASE_COEFF * np.transpose(on - off, (0, 3, 2, 1))

    def free_response(self, t, P_start):
        """Phase (n_t, N) from a steady start with drive P_start (W) when that drive is removed at t = 0."""
        P_start = np.asarray(P_start, float)
        return self.G @ P_start - PHASE_COEFF * np.einsum("tkn,k->tn", self.step(t), P_start)

    def phase(self, t, t_switch, P, P_start):
        """Phase (n_t, N) for the drive P (K, J) on the switch grid, from a steady start at P_start."""
        out = self.free_response(t, P_start)
        for i in range(0, len(t), 500):          # in chunks, to keep the pulse matrix small
            A = self.pulse_matrix(t[i:i + 500], t_switch)
            out[i:i + 500] += np.einsum("tnkj,kj->tn", A, P)
        return out

    def step_drive_phase(self, t, P_final, P_start):
        """Phase for the plain step: every heater jumps from P_start to P_final at t = 0."""
        return self.free_response(t, P_start) + PHASE_COEFF * np.einsum(
            "tkn,k->tn", self.step(t), np.asarray(P_final, float))


def grids(T, n_fast=30, n_slow=15, free_factor=3.0, per_interval=10, n_tail=60, t_flat=0.0316):
    """Switch times and check times for 'settled from T on'.

    The drive changes freely on [0, T_free], T_free = free_factor·T (n_fast
    intervals before T, n_slow after), then holds. Errors are checked from T to
    T_free at per_interval points per drive interval (dense, so that no error
    bump can hide between checks), then at log-spaced times until every pulse
    has died away (T_free + t_flat, where the step responses are flat), and
    once more at the end, which stands for steady state.
    """
    T_free = free_factor * T
    t_switch = np.concatenate([np.linspace(0, T, n_fast, endpoint=False),
                               np.linspace(T, T_free, n_slow, endpoint=False), [T_free]])
    t_check = np.concatenate([np.linspace(T, T_free, per_interval * n_slow + 1),
                              T_free + np.logspace(-6, np.log10(t_flat), n_tail),
                              [T_free + 2 * t_flat]])
    t_before = np.linspace(0, T, 3 * n_fast, endpoint=False)[1:]   # for the optional overshoot cap
    return t_switch, t_check, t_before


def _seed_directions(N, degree=6):
    """Unit, mean-free error shapes to seed the cuts: Legendre polynomials of degree 1..degree, both signs."""
    V = np.polynomial.legendre.legvander(np.linspace(-1, 1, N), degree)[:, 1:]
    V = V - V.mean(axis=0)
    V = V / np.linalg.norm(V, axis=0)
    return np.hstack([V, -V]).T                                          # (2·degree, N)


def design_drive(model, T, P_start, target, cap, P_final=None, sigma=SIGMA_STREHL, phase_cap=None,
                 tol=0.005, max_iter=300, **grid_kw):
    """A drive keeping the rms error ≤ sigma from T on, or None if there is none.

    target: the ramp r_n (rad). cap: largest total power of all heaters together
    at any moment (W), the chip's supply limit; it also caps each heater.
    P_final: the hold drive the waveform must end on (W per heater). Pass the
    plain step's steady drive, so both end in the same state. Among the drives
    that settle by T, the one returned departs least from that hold drive,
    Σ dt·|P - P_final|: it boosts only as much as it must. Without P_final it is
    the least-energy drive, which also chooses its own hold drive; that hold is
    cheaper but its beam only just meets the tolerance (worse Strehl and pointing
    than the plain step's final state), so pass P_final for fair comparisons.
    phase_cap: if given, no core's phase may exceed it at any time (overshoot limit).
    model may be a list of PulseModels (same heaters): the drive must then work on
    every one of them, e.g. the nominal thermal model and versions with its time
    scale a few per cent off (robust design). The constraints simply stack.
    The cuts use the radius sigma/(1 + tol) and the loop accepts errors up to
    sigma, so it ends within a finite number of rounds and the result keeps σ ≤ sigma.
    """
    models = list(model) if isinstance(model, (list, tuple)) else [model]
    t_switch, t_check, t_before = grids(T, **grid_kw)
    K, N, J, C = models[0].K, models[0].N, len(t_switch), len(t_check)
    A, b = [], []
    for M in models:                                                     # one block of check times per model
        A.append(M.pulse_matrix(t_check, t_switch).reshape(C, N, K * J))
        b.append(M.free_response(t_check, P_start) - target)
    A, b = np.concatenate(A), np.concatenate(b)
    A = A - A.mean(axis=1, keepdims=True)                                # mean-free: offset drops out
    b = b - b.mean(axis=1, keepdims=True)                                # error_i = A_i·x + b_i
    radius = sigma / (1 + tol)
    U = _seed_directions(N)
    base = np.unique(np.append(np.arange(0, C, 4), C - 1))             # seed cuts on every 4th check time
    seed = np.concatenate([i * C + base for i in range(len(models))])
    rows = [np.einsum("dn,cnv->cdv", U, A[seed]).reshape(-1, K * J)]
    rhs = [(radius * np.sqrt(N) - np.einsum("dn,cn->cd", U, b[seed])).reshape(-1)]
    if phase_cap is not None:   # every core's phase ≤ phase_cap, also during the run-up
        tb = np.concatenate([t_before, t_check])
        for M in models:
            rows.append(M.pulse_matrix(tb, t_switch).reshape(len(tb) * N, K * J))
            rhs.append(phase_cap - M.free_response(tb, P_start).reshape(-1))
    rows.append(np.tile(np.eye(J), K))                                  # Σ_k P_kj ≤ cap, every interval
    rhs.append(np.full(J, float(cap)))
    bounds = [(0, cap)] * (K * J)
    dt = np.append(np.diff(t_switch), t_switch[-1])
    if P_final is None:
        # Energy: each interval's power × its length; the open-ended hold counts for T_free.
        cost, extra_rows, extra_rhs, n_slack = np.tile(dt, K), np.zeros((0, K * J)), np.zeros(0), 0
    else:
        P_final = np.broadcast_to(np.asarray(P_final, float), (K,))
        for k in range(K):                                               # the open-ended last interval
            bounds[k * J + J - 1] = (P_final[k], P_final[k])
        # Slacks s_kj ≥ |P_kj - P_final,k| on the finite intervals; minimise Σ dt·s.
        pick = np.eye(K * J)[[k * J + j for k in range(K) for j in range(J - 1)]]
        n_slack = len(pick)
        extra_rows = np.vstack([np.hstack([pick, -np.eye(n_slack)]), np.hstack([-pick, -np.eye(n_slack)])])
        extra_rhs = np.concatenate([np.repeat(P_final, J - 1), -np.repeat(P_final, J - 1)])
        cost = np.concatenate([np.zeros(K * J), np.tile(dt[:-1], K)])
        bounds = bounds + [(0, None)] * n_slack
    for _ in range(max_iter):
        A_ub = np.vstack([np.hstack([np.vstack(rows), np.zeros((sum(len(r) for r in rows), n_slack))]),
                          extra_rows.reshape(-1, K * J + n_slack)])
        b_ub = np.concatenate([np.concatenate(rhs), extra_rhs])
        res = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if res.status not in (0, 2):                                     # numerical trouble: try interior point
            res = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs-ipm")
        if res.status == 2:
            return None
        if res.status != 0:
            raise RuntimeError(f"LP failed: {res.message}")
        x = res.x[:K * J]
        e = A @ x + b                                                    # (models·C, N)
        norms = np.linalg.norm(e, axis=1) / np.sqrt(N)                   # rms per check time
        bad = np.flatnonzero(norms > sigma)
        if bad.size == 0:
            return {"T": T, "t_switch": t_switch, "P": x.reshape(K, J), "objective": res.fun,
                    "rms_max": float(norms.max())}
        u = e[bad] / np.linalg.norm(e[bad], axis=1, keepdims=True)       # tangent cut u·e ≤ radius·√N
        rows.append(np.einsum("cn,cnv->cv", u, A[bad]))
        rhs.append(radius * np.sqrt(N) - np.einsum("cn,cn->c", u, b[bad]))
    raise RuntimeError("cutting planes did not converge")


def fastest_drive(model, P_start, target, cap, T_hi, rel_tol=0.01, **kw):
    """Shortest settling time T (to rel_tol) for which design_drive finds a drive, and that drive.

    The grid is rescaled with T, so feasibility is monotone in T only up to
    grid resolution (checked in notes/preemphasis.md: settling time falls
    monotonically with the power cap).
    """
    best = design_drive(model, T_hi, P_start, target, cap, **kw)
    while best is None:
        T_hi *= 1.5
        if T_hi > 1.0:
            raise RuntimeError("no feasible drive below 1 s: the target can't be held at all")
        best = design_drive(model, T_hi, P_start, target, cap, **kw)
    lo = 0.0
    while T_hi - lo > rel_tol * T_hi:
        mid = 0.5 * (lo + T_hi)
        sol = design_drive(model, mid, P_start, target, cap, **kw)
        if sol is None:
            lo = mid
        else:
            T_hi, best = mid, sol
    return best


def fine_times(t_end, n=3000):
    """Fine evaluation grid: uniform to t_end, then log-spaced until the step responses are flat."""
    return np.concatenate([np.linspace(1e-7, t_end, n), t_end + np.logspace(-6, np.log10(0.0316), 200)])


def verify(model, drive, P_start, target):
    """Strehl settling time, the worst rms error after T, and the peak core phase, on a fine grid."""
    T_free = drive["t_switch"][-1]
    t = fine_times(T_free)
    phi = model.phase(t, drive["t_switch"], drive["P"], P_start)
    s = strehl(phi - target)
    after = t >= drive["T"]
    return {
        "t_strehl": settle_time(t, s, STREHL_MIN),
        "worst_rms_after_T": float(rms_error(phi[after] - target).max()),
        "peak_phase": float(phi.max()),
        "final_peak_phase": float(phi[-1].max()),
    }
