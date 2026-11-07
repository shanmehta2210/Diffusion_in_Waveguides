"""Which phase patterns heat can make: singular values of the heater-to-phase matrix.

Put a heater under every emitter, the most any layout of this stack could do,
and write φ = G·P. The singular value decomposition G = U·S·Vᵀ splits this into
independent modes: heater pattern v_i makes phase pattern u_i with gain s_i
(rad per W). A pattern with a small s_i needs a lot of power. A shared offset
costs nothing (it doesn't move the beam), so it is projected out first: the
mean is removed from every column of G.

Heat diffusion is a spatial low-pass filter, so the modes come out ordered by
spatial frequency: offset (removed), ramp, dome, then faster and faster ripples,
each with a smaller gain. For a kernel that only depends on x - x_h (the
standard stack), the gains approach the kernel's Fourier transform at the
mode's spatial frequency (Szegő's theorem for Toeplitz matrices). So a kernel
of half-width L_th passes about N·d/L_th usable modes: wide spreading (large
L_th/d) leaves only the few smoothest patterns.

The usable count is the number of modes whose gain is at least a fraction tau
of the largest: tau = 0.1 means the mode costs at most 10× the power per rad of
the cheapest one.
"""
import numpy as np

from .constants import PHASE_COEFF


def mean_free(G):
    """Remove each column's mean over emitters: the shared offset is free."""
    return G - G.mean(axis=0, keepdims=True)


def modes(G):
    """Singular values (rad per W) and the phase patterns u_i (columns), offset removed."""
    U, s, Vt = np.linalg.svd(mean_free(G), full_matrices=False)
    return s, U, Vt.T


def usable(s, tau):
    """How many modes have gain at least tau times the largest."""
    return int(np.sum(s >= tau * s[0]))


def half_width(x, profile):
    """Half-width at half maximum of a peaked profile (µm), interpolated."""
    peak = int(np.argmax(profile))
    half = profile[peak] / 2
    right = peak + int(np.flatnonzero(profile[peak:] < half)[0])
    left = peak - int(np.flatnonzero(profile[:peak + 1][::-1] < half)[0])
    xr = np.interp(half, [profile[right], profile[right - 1]], [x[right], x[right - 1]])
    xl = np.interp(half, [profile[left], profile[left + 1]], [x[left], x[left + 1]])
    return 0.5 * (xr - xl)


def toeplitz_spectrum(kernel_values, pitch, n):
    """Fourier transform of a sampled shift-invariant kernel at the n cosine frequencies πm/(n·d).

    kernel_values(dx) gives G for a separation dx. For a large array the
    singular values of the Toeplitz matrix approach these numbers (Szegő).
    """
    j = np.arange(-(n - 1), n)
    k = kernel_values(j * pitch)
    q = np.pi * np.arange(n) / (n * pitch)
    return np.abs(np.cos(np.outer(q, j * pitch)) @ k)


def two_beam_targets(x, u_split):
    """Phase patterns that split the beam into two at u = ±u_split.

    binary: 0/π grating, each first order gets (2/π)² ≈ 0.405 of the light (0.81 total);
    sine:   1.84·cos(k0 x u), each first order J1(1.84)² ≈ 0.339 (0.68 total), and smooth.
    """
    from .constants import K0
    c = np.cos(K0 * u_split * x)
    return {"binary": np.where(c < 0, np.pi, 0.0), "sine": 1.84 * c}
