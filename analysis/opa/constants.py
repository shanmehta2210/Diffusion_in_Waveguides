"""Physical constants and fixed design choices (values from the project brief)."""
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parents[2]  # the Project folder holding the CSVs

WAVELENGTH_UM = 1.55
K0 = 2 * np.pi / WAVELENGTH_UM  # free-space wavenumber, rad/µm

N_EFF = 2.44546      # strip effective index at 293.15 K (Model A)
DNEFF_DT = 2.009e-4  # 1/K (Model A)
N_GROUP = 4.18

# Phase per (temperature rise x length): (2π/λ)·dn_eff/dT with λ in metres, so
# about 814.4 rad K^-1 m^-1. Multiplied by a kernel value in K per (W/m) and a
# heater power in W, it gives radians.
PHASE_COEFF = 2 * np.pi / (WAVELENGTH_UM * 1e-6) * DNEFF_DT

# Odd count: one emitter at x = 0, the rest every PITCH_UM either side. 33 at
# 2 µm puts them at 0, ±2, ..., ±32 µm, exactly on the thermal model's cores.
N_EMITTERS = 33
PITCH_UM = 2.0
MIN_PITCH_UM = 1.5  # below this, neighbouring waveguides couple optically (Model B)
