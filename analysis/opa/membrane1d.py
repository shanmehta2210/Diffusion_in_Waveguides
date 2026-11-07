"""1D membrane model of the slot design, with heat lost to gas (est., calibrated to COMSOL).

The membrane runs from the far cavity wall (x = -40 µm, the anchor and heat sink) to the
slot edge (x = +36 µm). Along it,
    k·t·T'' - h·T = -q·δ(x - x_h)
  - far wall: T = 0;
  - slot edge: the heat crossing the slot's gas gap, -k·t·T' = G_s·T, with G_s = a_s·k_gas·t/w;
  - below: gas in the cavity of depth D, h = a_h·k_gas/D.
k·t = 5.8e-6 W/K is the membrane fit from the undercut kernel (beyond_steering.md). The two
factors a_s and a_h absorb the 2D shape of the gas paths; run_gas_loss.py fits them to the
COMSOL slot runs R2 (2 µm slot, air) and R3b (4 µm slot, air), within 3.5%, and the vacuum run
R3a checks the result.

Gas at low pressure (temperature-jump model, parallel plates): once the mean free path λ is
comparable to the gap g, conduction falls as k_eff = k_gas·g/(g + 4.0·λ), with λ = 68 nm for
air at 1 atm and 300 K, inversely proportional to pressure (accommodation 0.9).
"""
import numpy as np

from .thermal import emitter_positions

KT = 5.8e-6            # membrane k·t, W/K
T_MEM = 4.22e-6        # membrane thickness, m
K_AIR = 0.026          # W/m/K at 1 atm, 300 K
LAMBDA_ATM = 68e-9     # mean free path of air at 1 atm and 300 K, m
X0, X1 = -40.0, 36.0   # far wall and slot edge, µm
A_H, A_S = 1.29, 1.55  # gas-path factors (cavity, slot) as fitted by run_gas_loss.py


def k_gas(gap_m, pressure_pa, k_bulk=K_AIR, lambda_atm=LAMBDA_ATM):
    """Effective gas conductivity across a gap at a given pressure (temperature-jump model)."""
    lam = lambda_atm * 101325.0 / pressure_pa
    return k_bulk * gap_m / (gap_m + 4.0 * lam)


def profiles(h, G_s, heaters, x_emit=None, dx=0.1):
    """Temperature rise (K per W/m) at the emitters, one column per heater position (µm).

    Finite differences on a 0.1 µm grid; each heater is a 1 W/m line source.
    """
    if x_emit is None:
        x_emit = emitter_positions()
    x = np.arange(X0, X1 + dx / 2, dx)
    n = len(x)
    d = dx * 1e-6
    A = (np.diag(np.full(n, -2 * KT / d**2 - h))
         + np.diag(np.full(n - 1, KT / d**2), 1) + np.diag(np.full(n - 1, KT / d**2), -1))
    A[0, :] = 0
    A[0, 0] = 1                                        # far wall: T = 0
    A[-1, -2] = 2 * KT / d**2                          # slot edge: ghost node with -k t T' = G_s T
    A[-1, -1] = -2 * KT / d**2 - h - 2 * G_s / d
    cols = []
    for xh in heaters:
        b = np.zeros(n)
        b[int(round((xh - X0) / dx))] = -1.0 / d       # 1 W/m point source
        b[0] = 0
        cols.append(np.interp(x_emit, x, np.linalg.solve(A, b)))
    return np.array(cols).T


def model(a_h, a_s, D_um, w_um, heaters, pressure=101325.0, k_bulk=K_AIR, lambda_atm=LAMBDA_ATM):
    """Profiles for a cavity of depth D_um and a slot of width w_um, filled with gas at a pressure."""
    h = a_h * k_gas(D_um * 1e-6, pressure, k_bulk, lambda_atm) / (D_um * 1e-6)
    G_s = a_s * k_gas(w_um * 1e-6, pressure, k_bulk, lambda_atm) * T_MEM / (w_um * 1e-6)
    return profiles(h, G_s, heaters)
