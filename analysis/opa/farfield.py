"""Far field across the array, and the beam metrics.

Angles are handled as u = sin(theta), measured across the array (in the plane
perpendicular to the waveguides). The transverse wavevector k0·u is conserved
through the flat oxide/air interface, so u is directly the angle in air.

With emitter amplitudes a_n at positions x_n and phases φ_n,
    I(u) = EP(u) · |Σ a_n exp(i(φ_n - k0 x_n u))|² / (Σ|a_n|)²,
normalised so that a perfect beam with EP = 1 has peak 1. A uniform ramp
φ_n = Δφ·x_n/d puts the beam at u = Δφ·λ/(2πd).
"""
import numpy as np

from .comsol_io import read_comsol
from .constants import DATA_DIR, K0, PITCH_UM, WAVELENGTH_UM

# Step 1e-4 in sin(theta); the beamwidth at N = 32, d = 2 µm is about 0.02.
U_GRID = np.linspace(-1, 1, 20001)


def element_pattern(data_dir=DATA_DIR, file="mfield_Tmin.csv"):
    """Single-emitter intensity pattern EP(u), normalised to 1 at u = 0.

    Approximation (stated in the report): the emitted field across the array
    has the waveguide mode's x-profile, so EP is the |Fourier transform|² of
    Ex summed over y. Returns a function of u.
    """
    _, data = read_comsol(data_dir / file)
    x = np.unique(data[:, 0])
    assert np.allclose(data[: len(x), 0], x), "expected x to vary fastest in the export"
    # The mode is almost purely imaginary (arbitrary global phase), so keep both parts.
    ex = (data[:, 2] + 1j * data[:, 3]).reshape(-1, len(x))  # rows: y, columns: x
    line = ex.sum(axis=0)

    u = np.linspace(-1, 1, 2001)
    transform = np.exp(-1j * K0 * np.outer(u, x)) @ line
    ep = np.abs(transform) ** 2 / np.abs(line.sum()) ** 2

    def ep_of_u(uu):
        return np.interp(uu, u, ep)

    return ep_of_u


def array_factor(phases, x_emit, u=U_GRID, amplitudes=None):
    """|Σ a_n exp(i(φ_n - k0 x_n u))|² / (Σ|a_n|)²: 1 where all emitters add in phase."""
    a = np.ones(len(x_emit)) if amplitudes is None else np.asarray(amplitudes)
    field = np.exp(1j * (np.asarray(phases)[None, :] - K0 * np.outer(u, x_emit))) @ a
    return np.abs(field) ** 2 / np.sum(np.abs(a)) ** 2


def far_field(phases, x_emit, ep, u=U_GRID, amplitudes=None):
    """Intensity I(u) = EP(u) · array factor, sampled on u."""
    return ep(u) * array_factor(phases, x_emit, u, amplitudes)


def strehl(phases, x_emit, u0, amplitudes=None):
    """Intensity toward u0 relative to a perfect ramp aimed at u0 (EP cancels).

    For random phase errors of spread σ its average is about exp(-σ²), the Ruze relation.
    """
    return array_factor(phases, x_emit, np.atleast_1d(u0), amplitudes)[0]


def _half_max_width(u, intensity, lo, i_pk, hi):
    """Width in u between the half-maximum crossings either side of the peak."""
    half = intensity[i_pk] / 2
    left = np.arange(lo, i_pk + 1)
    right = np.arange(i_pk, hi + 1)
    below_left = left[intensity[left] < half]
    below_right = right[intensity[right] < half]
    if below_left.size == 0 or below_right.size == 0:
        return np.nan, np.nan, np.nan  # the main lobe never drops to half: no clean beam
    j, k = below_left[-1], below_right[0]
    u_left = np.interp(half, intensity[j : j + 2], u[j : j + 2])
    u_right = np.interp(half, intensity[k - 1 : k + 1][::-1], u[k - 1 : k + 1][::-1])
    return u_right - u_left, u_left, u_right


def beam_metrics(u, intensity, ep, pitch=PITCH_UM):
    """Beam angle, width, side lobes, main-lobe power and Strehl ratio of a far field.

    Side lobes are reported two ways, because at d = 2 µm > λ grating lobes
    (copies of the beam λ/d away in u) are always visible:
      sls_db          main peak over the largest lobe within ±λ/(2d) of it, the
                      window holding one copy of the pattern (-13.3 dB for a uniform array)
      grating_lobe_db main peak over the strongest copy outside that window
      sls_all_db      main peak over the largest lobe anywhere (the brief's definition)
    All three are positive dB; larger is better.
    """
    i_pk = int(np.argmax(intensity))
    u_pk, i_max = u[i_pk], intensity[i_pk]

    # Main lobe: from the peak out to the first minimum on each side.
    lo = i_pk
    while lo > 0 and intensity[lo - 1] < intensity[lo]:
        lo -= 1
    hi = i_pk
    while hi < len(u) - 1 and intensity[hi + 1] < intensity[hi]:
        hi += 1
    outside = np.ones(len(u), dtype=bool)
    outside[lo : hi + 1] = False

    width_u, u_left, u_right = _half_max_width(u, intensity, lo, i_pk, hi)
    fwhm_deg = np.degrees(np.arcsin(u_right) - np.arcsin(u_left))

    window = np.abs(u - u_pk) < WAVELENGTH_UM / (2 * pitch)

    def suppression_db(mask):
        return 10 * np.log10(i_max / intensity[mask].max()) if mask.any() else np.inf

    # In-plane power fraction: integrate over angle, theta = arcsin(u).
    theta = np.arcsin(u)
    main_power = np.trapezoid(intensity[lo : hi + 1], theta[lo : hi + 1])

    return {
        "angle_deg": np.degrees(np.arcsin(u_pk)),
        "u_peak": u_pk,
        "fwhm_u": width_u,
        "fwhm_deg": fwhm_deg,
        "sls_db": suppression_db(outside & window),
        "grating_lobe_db": suppression_db(~window),
        "sls_all_db": suppression_db(outside),
        "main_lobe_fraction": main_power / np.trapezoid(intensity, theta),
        "strehl": i_max / ep(u_pk),
    }
