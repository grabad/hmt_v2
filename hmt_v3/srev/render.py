"""
PSF rendering, camera noise, and TIFF stack writing for the SR-EV pipeline.

PSF is an astigmatic 2D Gaussian, `sigma_{x,y}(z)` following the same
functional form ThunderSTORM's own Astigmatism3D fitter assumes and
calibrates -- using it here means ground-truth z produces widths via the
exact formula ThunderSTORM will later try to invert, so a downstream
ThunderSTORM run tests its real fitting bias/uncertainty rather than a
PSF-model mismatch. Each emitter is rasterized as a pixel-INTEGRATED
Gaussian (via `scipy.special.erf`) over a small bounded window, not
point-sampled and not FFT'd -- the direct fix for the legacy MATLAB
pipeline's fresh 2D FFT of a binary pupil per emitter per frame.

Noise is applied in real ADU space (photons -> QE -> Poisson shot noise ->
EM-gain Gamma stage -> read noise -> conversion gain -> baseline -> clip),
not the photon-equivalent excess-noise-factor shortcut the legacy pipeline
used -- the output needs to look like literal raw sensor counts for a user
to configure ThunderSTORM's own camera-setup dialog against it.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import erf

try:
    import tifffile
except ImportError:  # pragma: no cover - surfaced at call time, not import time
    tifffile = None


@dataclass
class RenderConfig:
    """
    Astigmatic PSF calibration:
        sigma_{x,y}(z) = sigma0_nm * sqrt(1 + u^2 + astig_a*u^3 + astig_b*u^4),
        u = (z - c_{x,y}) / depth_scale_nm,  c_x = +focal_offset_nm, c_y = -focal_offset_nm

    astig_a/astig_b default to 0 (a pure symmetric quadratic defocus curve
    in each axis) -- the cubic/quartic terms exist to match a REAL
    ThunderSTORM astigmatism calibration if one is available; there is no
    universal default for them.
    """
    sigma0_nm: float = 150.0
    focal_offset_nm: float = 250.0
    depth_scale_nm: float = 400.0
    astig_a: float = 0.0
    astig_b: float = 0.0
    window_sigma_multiple: float = 4.0


@dataclass
class CameraConfig:
    """
    Field names deliberately match ThunderSTORM's own camera-setup dialog
    so a user can copy values across directly when telling ThunderSTORM
    how to convert this stack's ADU counts back to photons.
    """
    pixel_size_nm: float = 100.0
    frame_shape: tuple = None  # (height, width) in pixels; None = not yet sized, see fit_frame_to_points
    quantum_efficiency: float = 0.9
    em_gain: float = 1.0
    e_per_adu: float = 1.0
    read_noise_e: float = 1.0
    base_level_adu: float = 100.0
    background_photons_per_pixel: float = 5.0
    bit_depth_max_adu: int = 65535


def fit_frame_to_points(x_nm, y_nm, pixel_size_nm, margin_nm=300.0):
    """
    Computes a (height, width) frame size in pixels that comfortably
    contains every given point plus a margin, along with the
    (x_offset_nm, y_offset_nm) that shifts those points' own coordinate
    system so its minimum (inclusive of margin) lands at pixel (0, 0).

    render_photon_frame/render_stack_to_tiff take x_nm/y_nm as already
    frame-relative (pixel (0, 0) == nm (0, 0)) with no offset concept of
    their own -- config-N.txt's coordinates are centered on an arbitrary
    origin (frequently not at zero, and often negative), so this offset
    must be applied to every coordinate column that shares that
    coordinate system (monomers, primaries, secondaries, fluorophores
    alike) BEFORE calling any rendering function, all by the same amount,
    so they stay mutually consistent. See pipeline.run_pipeline for where
    this is actually applied.

    Args:
        x_nm, y_nm:    emitter (or any reference point set's) coordinates
                       in their original, arbitrary-origin coordinate
                       system
        pixel_size_nm: camera pixel size
        margin_nm:     padding added on every side, so a PSF centered
                       near the point set's edge isn't immediately
                       clipped by the frame boundary

    Returns:
        (frame_shape, x_offset_nm, y_offset_nm): frame_shape is
        (height, width) in pixels; adding (x_offset_nm, y_offset_nm) to a
        point's (x_nm, y_nm) moves it into this frame's coordinate system
    """
    x_nm = np.asarray(x_nm, dtype=float)
    y_nm = np.asarray(y_nm, dtype=float)
    if len(x_nm) == 0:
        raise ValueError("fit_frame_to_points needs at least one point")

    x_offset_nm = margin_nm - x_nm.min()
    y_offset_nm = margin_nm - y_nm.min()
    width_px = int(np.ceil((x_nm.max() + x_offset_nm + margin_nm) / pixel_size_nm))
    height_px = int(np.ceil((y_nm.max() + y_offset_nm + margin_nm) / pixel_size_nm))
    return (height_px, width_px), x_offset_nm, y_offset_nm


def astigmatic_sigma_nm(z_nm, config):
    """
    Returns (sigma_x_nm, sigma_y_nm) for the given z (nm), per RenderConfig's
    astigmatic calibration. Vectorized over an array of z values.
    """
    z_nm = np.asarray(z_nm, dtype=float)

    def _sigma(z, c):
        u = (z - c) / config.depth_scale_nm
        val = 1.0 + u**2 + config.astig_a * u**3 + config.astig_b * u**4
        return config.sigma0_nm * np.sqrt(np.maximum(val, 1e-6))

    sigma_x = _sigma(z_nm, config.focal_offset_nm)
    sigma_y = _sigma(z_nm, -config.focal_offset_nm)
    return sigma_x, sigma_y


def _pixel_integrated_gaussian(cx_px, cy_px, sigma_x_px, sigma_y_px, half_window_px):
    """
    A local (2*half_window_px+1)^2 patch of a 2D Gaussian integrated over
    each pixel's area (the standard separable erf-difference form), not
    point-sampled -- point-sampling is a known systematic bias when the
    PSF is only a few pixels wide, as it typically is here. Normalized to
    sum to ~1.0 (short of the window's truncation loss); the caller scales
    by the emitter's photon count.

    Returns (y0, x0, patch): (y0, x0) is the integer pixel index of the
    patch's top-left corner, for slicing directly into a frame buffer.
    """
    x_center_px = int(np.floor(cx_px))
    y_center_px = int(np.floor(cy_px))
    x0 = x_center_px - half_window_px
    y0 = y_center_px - half_window_px
    xs = np.arange(x0, x0 + 2 * half_window_px + 1, dtype=float)
    ys = np.arange(y0, y0 + 2 * half_window_px + 1, dtype=float)

    def _axis_profile(pix_coords, center, sigma):
        sigma = max(sigma, 1e-3)
        edges_lo = pix_coords - 0.5
        edges_hi = pix_coords + 0.5
        return 0.5 * (
            erf((edges_hi - center) / (np.sqrt(2) * sigma))
            - erf((edges_lo - center) / (np.sqrt(2) * sigma))
        )

    profile_x = _axis_profile(xs, cx_px, sigma_x_px)
    profile_y = _axis_profile(ys, cy_px, sigma_y_px)
    patch = np.outer(profile_y, profile_x)
    return y0, x0, patch


def render_photon_frame(x_nm, y_nm, z_nm, photons, frame_shape, camera_config, render_config):
    """
    Renders one frame's expected-photon image (no noise, no background)
    from a handful of active emitters -- a Python loop over emitters, not
    over pixels or frames, so its cost scales with typical SMLM sparsity
    (tens to hundreds of active emitters) rather than image size or movie
    length.

    Args:
        x_nm, y_nm, z_nm: (n,) true emitter positions for THIS frame only
        photons:          (n,) expected photon count for this frame
        frame_shape:      (height, width) in pixels
        camera_config:    CameraConfig (uses pixel_size_nm)
        render_config:    RenderConfig (PSF calibration + window size)

    Returns:
        (H, W) float64 array of expected photons per pixel
    """
    H, W = frame_shape
    buf = np.zeros((H, W), dtype=np.float64)
    if len(x_nm) == 0:
        return buf

    sigma_x_nm, sigma_y_nm = astigmatic_sigma_nm(z_nm, render_config)
    sigma_x_px = sigma_x_nm / camera_config.pixel_size_nm
    sigma_y_px = sigma_y_nm / camera_config.pixel_size_nm
    cx_px = np.asarray(x_nm, dtype=float) / camera_config.pixel_size_nm
    cy_px = np.asarray(y_nm, dtype=float) / camera_config.pixel_size_nm

    for i in range(len(x_nm)):
        half_window = int(np.ceil(
            render_config.window_sigma_multiple * max(sigma_x_px[i], sigma_y_px[i])
        )) + 1
        y0, x0, patch = _pixel_integrated_gaussian(
            cx_px[i], cy_px[i], sigma_x_px[i], sigma_y_px[i], half_window
        )
        y1, x1 = y0 + patch.shape[0], x0 + patch.shape[1]

        cy0, cy1 = max(y0, 0), min(y1, H)
        cx0, cx1 = max(x0, 0), min(x1, W)
        if cy0 >= cy1 or cx0 >= cx1:
            continue  # emitter's whole window falls outside the frame

        py0, py1 = cy0 - y0, patch.shape[0] - (y1 - cy1)
        px0, px1 = cx0 - x0, patch.shape[1] - (x1 - cx1)
        buf[cy0:cy1, cx0:cx1] += patch[py0:py1, px0:px1] * photons[i]

    return buf


def apply_camera_noise(photon_img, camera_config, rng):
    """
    Full noise chain in real ADU space:
      1. mu_e = (photons + background) * QE
      2. shot noise:  N_e ~ Poisson(mu_e)             [applied ONCE, to signal+background combined]
      3. EM gain:     N_out_e ~ Gamma(shape=N_e, scale=em_gain)
                       -- the standard continuous approximation of an EM
                       register's stochastic avalanche multiplication,
                       degenerating to (approximately) N_e itself when
                       em_gain == 1.0 (skipped entirely in that case, both
                       for speed and because Gamma(shape=0, ...) is
                       undefined and must be masked out either way)
      4. read noise:  + Normal(0, read_noise_e), added AFTER gain -- real
                       EMCCD/sCMOS read noise is injected at the output
                       amplifier, post-gain-register
      5. ADU = N_out_e / e_per_adu + base_level_adu, clipped to
         [0, bit_depth_max_adu] and quantized to uint16

    Args:
        photon_img:    (H, W) expected photons per pixel (signal only;
                       background is added inside this function)
        camera_config: CameraConfig
        rng:           np.random.Generator

    Returns:
        (H, W) uint16 array -- literal raw-camera-style pixel counts
    """
    mu_e = (photon_img + camera_config.background_photons_per_pixel) * camera_config.quantum_efficiency
    n_e = rng.poisson(np.maximum(mu_e, 0.0))

    if camera_config.em_gain > 1.0:
        n_out_e = np.zeros_like(n_e, dtype=np.float64)
        nonzero = n_e > 0
        n_out_e[nonzero] = rng.gamma(shape=n_e[nonzero], scale=camera_config.em_gain)
    else:
        n_out_e = n_e.astype(np.float64)

    n_out_e += rng.normal(0.0, camera_config.read_noise_e, size=photon_img.shape)

    adu = n_out_e / camera_config.e_per_adu + camera_config.base_level_adu
    adu = np.clip(adu, 0, camera_config.bit_depth_max_adu)
    return adu.astype(np.uint16)


def expand_active_frames(fluorophores_df, blink_events_df):
    """
    Explodes blink_events (frame_start/frame_end/total_photons/n_frames_on,
    one row per blink) into one row per (frame, fluorophore) for every
    frame that emitter was actually ON, joined to its true static x/y/z --
    the sparse, event-driven table render_stack_to_tiff loops over, sized
    by total active-emitter-frames rather than n_frames * n_fluorophores.

    Args:
        fluorophores_df: output of labeling.place_fluorophores
        blink_events_df: output of photophysics.simulate_blinking

    Returns:
        DataFrame: frame, fluorophore_id, blink_id, photon_rate (constant
                   per-frame photon count for that blink, i.e.
                   total_photons / n_frames_on), x/y/z [nm]
    """
    n_on = blink_events_df["n_frames_on"].to_numpy()
    frame_start = blink_events_df["frame_start"].to_numpy()
    photon_rate = (blink_events_df["total_photons"] / blink_events_df["n_frames_on"]).to_numpy()

    frame = np.concatenate([np.arange(s, s + n) for s, n in zip(frame_start, n_on)]) if len(n_on) else np.array([], dtype=np.int64)
    rep_fluor = np.repeat(blink_events_df["fluorophore_id"].to_numpy(), n_on)
    rep_rate = np.repeat(photon_rate, n_on)
    rep_blink = np.repeat(blink_events_df["blink_id"].to_numpy(), n_on)

    active = pd.DataFrame({
        "frame": frame,
        "fluorophore_id": rep_fluor,
        "blink_id": rep_blink,
        "photon_rate": rep_rate,
    })
    return active.merge(
        fluorophores_df[["fluorophore_id", "x [nm]", "y [nm]", "z [nm]"]],
        on="fluorophore_id", how="left",
    )


def render_stack_to_tiff(fluorophores_df, blink_events_df, n_frames, camera_config,
                          render_config, output_path, rng=None):
    """
    Renders the full movie and streams it to a (Big)TIFF, one frame at a
    time -- never holding the whole stack in memory. `bigtiff=True` is not
    optional: a realistic stack (e.g. 512x512x25000 frames, uint16) is
    ~13 GB, over the classic 4 GB TIFF limit.

    Args:
        fluorophores_df: output of labeling.place_fluorophores
        blink_events_df: output of photophysics.simulate_blinking
        n_frames:        total frames to write (independent of the max
                         frame actually reached by blink_events_df, so a
                         run can pad trailing background-only frames)
        camera_config:   CameraConfig
        render_config:   RenderConfig
        output_path:     destination .tif/.tiff path
        rng:             np.random.Generator (default_rng() if None)
    """
    if tifffile is None:
        raise ImportError(
            "render_stack_to_tiff requires tifffile (pip install tifffile / "
            "add to environment.yml) -- not installed in this environment."
        )
    rng = rng if rng is not None else np.random.default_rng()

    active = expand_active_frames(fluorophores_df, blink_events_df)
    grouped = active.groupby("frame") if len(active) else None
    active_frames = set(grouped.groups.keys()) if grouped is not None else set()
    frame_shape = camera_config.frame_shape

    with tifffile.TiffWriter(output_path, bigtiff=True) as tif:
        for frame in range(n_frames):
            if frame in active_frames:
                g = grouped.get_group(frame)
                photon_img = render_photon_frame(
                    g["x [nm]"].to_numpy(), g["y [nm]"].to_numpy(), g["z [nm]"].to_numpy(),
                    g["photon_rate"].to_numpy(), frame_shape, camera_config, render_config,
                )
            else:
                photon_img = np.zeros(frame_shape, dtype=np.float64)

            adu_img = apply_camera_noise(photon_img, camera_config, rng)
            tif.write(adu_img, contiguous=True)
