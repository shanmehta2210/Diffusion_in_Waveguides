"""Reading COMSOL CSV exports.

COMSOL writes header lines starting with '%'; the last one names the columns.
Values are comma-separated. Some come out as complex strings such as
'2.43+4e-16i' (the imaginary part is solver noise, so we keep the real part),
and cells are empty where a geometry ended early.
"""
import csv

import numpy as np


def _to_float(s):
    s = s.strip()
    if not s:
        return np.nan
    if s.endswith("i"):
        return complex(s[:-1] + "j").real
    return float(s)


def read_comsol(path):
    """Return (column_names, data) for a COMSOL CSV; data is a 2D float array."""
    names, rows = [], []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("%"):
                names = next(csv.reader([line[1:].strip()]))
            elif line.strip():
                rows.append([_to_float(s) for s in line.split(",")])
    return names, np.array(rows)


def column_parameters(name):
    """Parse the sweep suffix of a column name into a dict.

    'T-T0 (K) @ uc=1; x_h=5'                  -> {'uc': 1.0, 'x_h': 5.0}
    'T-T0 (K) @ 3: uc=0; nl=1; Q_l=50 W/m; ...' -> {'uc': 0.0, 'nl': 1.0, 'Q_l': 50.0, ...}
    """
    if "@" not in name:
        return {}
    suffix = name.split("@", 1)[1]
    if ":" in suffix:  # drop the run number some exports add
        suffix = suffix.split(":", 1)[1]
    params = {}
    for item in suffix.split(";"):
        if "=" not in item:  # e.g. mode-field exports label columns by n_eff only
            continue
        key, value = item.split("=")
        params[key.strip()] = float(value.split()[0])
    return params
