"""Choosing where to put K heaters.

A layout is ranked first by its steering span (range toward +θ plus range
toward -θ, each at Strehl ≥ 0.8), then by its power per degree. Spans are
compared to 0.01°, so two layouts that steer equally far are separated by power.

The search is a pattern search: from a starting layout, try moving each heater
by ±step, keep any move that improves the ranking, and halve the step once no
move helps. It is run from several starting layouts because the ranking can
have more than one local optimum. For K = 1 and 2 an exhaustive grid
(grid_search) checks that it finds the global one.
"""
import itertools

import numpy as np

from .steering import steer_both_directions
from .thermal import phase_matrix

HEATER_WIDTH_UM = 2.0  # heaters closer than this would overlap


def rank_key(result):
    """Larger is better: longest span first, then least power per degree."""
    return (round(result["span_deg"], 2), -result["mw_per_deg"])


def evaluate_layout(kernel, x_emit, x_heat):
    return steer_both_directions(phase_matrix(kernel, x_emit, x_heat), x_emit)


def is_feasible(x_heat, bounds, min_gap=HEATER_WIDTH_UM):
    x = np.sort(x_heat)
    return x[0] >= bounds[0] and x[-1] <= bounds[1] and np.all(np.diff(x) >= min_gap - 1e-9)


def pattern_search(kernel, x_emit, x_start, bounds, steps=(4, 2, 1, 0.5, 0.25)):
    """Improve one starting layout; returns (positions, result)."""
    best_x = np.sort(np.asarray(x_start, dtype=float))
    best = evaluate_layout(kernel, x_emit, best_x)
    for step in steps:
        improved = True
        while improved:
            improved = False
            for k in range(len(best_x)):
                for move in (+step, -step):
                    trial = best_x.copy()
                    # Clip rather than reject, so a heater can land exactly on the
                    # boundary (the outermost heaters usually want to be there).
                    trial[k] = np.clip(trial[k] + move, *bounds)
                    if trial[k] == best_x[k] or not is_feasible(trial, bounds):
                        continue
                    result = evaluate_layout(kernel, x_emit, trial)
                    if rank_key(result) > rank_key(best):
                        best_x, best, improved = np.sort(trial), result, True
    return best_x, best


def random_layout(rng, n_heaters, bounds):
    while True:
        x = np.sort(rng.uniform(*bounds, n_heaters))
        if is_feasible(x, bounds):
            return x


def optimise_layout(kernel, x_emit, n_heaters, bounds, n_random_starts=5, seed=0):
    """Best layout of n_heaters found from an evenly spaced start plus random starts."""
    rng = np.random.default_rng(seed)
    starts = [np.linspace(*bounds, n_heaters) if n_heaters > 1 else np.array([bounds[1]])]
    starts += [random_layout(rng, n_heaters, bounds) for _ in range(n_random_starts)]
    found = [pattern_search(kernel, x_emit, s, bounds) for s in starts]
    return max(found, key=lambda pair: rank_key(pair[1]))


def grid_search(kernel, x_emit, n_heaters, bounds, step=1.0):
    """Try every layout on a grid of heater positions (only practical for K ≤ 2)."""
    grid = np.arange(bounds[0], bounds[1] + step / 2, step)
    best_x, best = None, None
    for combo in itertools.combinations(grid, n_heaters):
        x = np.array(combo)
        if not is_feasible(x, bounds):
            continue
        result = evaluate_layout(kernel, x_emit, x)
        if best is None or rank_key(result) > rank_key(best):
            best_x, best = x, result
    return best_x, best
