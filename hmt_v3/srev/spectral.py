"""
First-order (spectrally dispersed) rendering for dual-dye spectroscopic
SMLM (sSMLM), so ONE acquisition carries BOTH dyes blinking together and
the dyes are told apart afterwards by their emission spectra (e.g. with
RainbowSTORM) rather than by being simulated as separate movies.

This is the v3 counterpart of the legacy MATLAB chain in hmt_v1/Simulation:
    blinking_sim_unrealistic.m  -> zeroth-order (spatial) stack, both dyes
    make_first_order.m          -> each dye's emission spectrum resampled
                                   onto the spectral pixel grid
    make_real_fist.m            -> first-order stack: the PSF shifted along
                                   x by the dye's dispersion and convolved
                                   with its spectrum (conv2(h, first_ord'))
    add_noise.m                 -> noise on both stacks, first-order
                                   background 1.5x the zeroth-order one
    splitmeac.m                 -> (post RainbowSTORM) split the localization
                                   list on spectral centroid: 655-688 nm ->
                                   AF647 channel, 705-737 nm -> CF680 channel
Those scripts depended on three files that never made it into the repo
(`new_cal.mat` for the pixel->wavelength calibration, `count_dist.mat` /
`fwhm.mat` for the blink and width statistics, plus the two dye spectrum
CSVs). Everything they encoded is a declared parameter here instead:

  * DyeSpectrum        -- a dye's emission spectrum: either a parametric
                          two-sided Gaussian (sharp blue edge, longer red
                          tail, the generic shape of a cyanine dye) or a
                          tabulated (wavelength, intensity) table loaded
                          from a vendor CSV. The AF647 / CF680 presets are
                          tuned so their centroid over the default
                          detection window matches what RainbowSTORM
                          actually measures for those dyes on the FOIL
                          instrument (hmt_v1/Simulation/blink_stats:
                          AF647 ~670.5 nm, CF680 ~700 nm).
  * SpectralConfig     -- the instrument: the grating's dispersion
                          (nm per camera pixel) and where a reference
                          wavelength lands relative to the zeroth-order
                          spot (the calibration RainbowSTORM needs), the
                          first-order PSF (wider than the zeroth-order one:
                          ~200 nm vs ~145 nm measured in blink_stats), the
                          zeroth/first photon split, and the first-order
                          background ratio.

Forward model of the spectral (first-order) image, per emitter, per frame:
    I(x, y) = N * f1 * sum_k w_k * G(x - x0 - off0 - shift/D - k, y - y0; sigma1(z))
with N the emitter's photons that frame, f1 the first-order fraction,
w_k the dye's emission spectrum binned into integer pixel offsets k
(dispersed_kernel), off0 the reference offset, D the dispersion, shift the
emitter's own spectral heterogeneity in nm, and G the pixel-integrated
Gaussian PSF of the first-order path. The k-sum is a 1D convolution of the
PSF's x-profile with the binned spectrum -- the same PSF-times-spectrum
construction as make_real_fist.m, done separably instead of as a full 2D
conv2 per emitter. Emitters at the same y whose dispersed spectra overlap
in x corrupt each other's centroids, exactly as in a real acquisition;
that crosstalk is the point of simulating both dyes in one stack.

Calibration: the wavelength axis of the spectral image is
    offset_px(lambda) = reference_offset_px + (lambda - reference_wavelength_nm) / dispersion_nm_per_px
relative to the emitter's zeroth-order x position (wavelength_to_offset_px).
RainbowSTORM's spectral calibration step fits exactly this relation from a
calibration image with lines at known wavelengths -- render_calibration_pair
produces such a pair (pinholes in the spatial frame, monochromatic lines
in the spectral frame) from the same forward model, and
pipeline.run_multichannel_pipeline records the parameters and a
wavelength->offset table in run_config.json for entering them by hand.
"""
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from . import render as srev_render


@dataclass(frozen=True)
class DyeSpectrum:
    """
    One dye's emission spectrum. Frozen so the module-level presets can't
    be mutated in place -- use dataclasses.replace to derive a variant.

    Parametric form (when wavelengths_nm is None): a two-sided Gaussian
    peaking at peak_nm, with the blue-side and red-side sigmas in the
    ratio 1:red_tail_ratio and their sum fixed by fwhm_nm
    (FWHM = sqrt(2 ln 2) * (sigma_blue + sigma_red)). red_tail_ratio > 1
    gives the asymmetric red tail every cyanine dye has, which is also
    why a dye's centroid sits a few nm above its peak.

    Tabulated form (wavelengths_nm / intensities set, e.g. via from_csv):
    linearly interpolated, zero outside the table. peak_nm is then
    informational only.

    heterogeneity_nm: per-FLUOROPHORE random spectral shift (stdev, nm),
    drawn once per dye molecule (pipeline.simulate_channel) and applied as
    a rigid shift of its whole spectrum -- the single-molecule spectral
    heterogeneity the legacy pipeline declared as `SH=10 nm` but never
    wired up. 0 = every molecule of a dye has the identical spectrum.
    """
    name: str
    peak_nm: float = None
    fwhm_nm: float = 42.0
    red_tail_ratio: float = 1.4
    heterogeneity_nm: float = 0.0
    wavelengths_nm: tuple = None
    intensities: tuple = None

    def emission(self, wavelengths_nm):
        """Unnormalized emission intensity (>= 0) at each wavelength (nm)."""
        lam = np.asarray(wavelengths_nm, dtype=float)
        if self.wavelengths_nm is not None:
            return np.clip(
                np.interp(lam, np.asarray(self.wavelengths_nm, dtype=float),
                          np.asarray(self.intensities, dtype=float), left=0.0, right=0.0),
                0.0, None,
            )
        if self.peak_nm is None:
            raise ValueError(f"DyeSpectrum {self.name!r} needs peak_nm or a tabulated spectrum")
        sigma_blue = self.fwhm_nm / (np.sqrt(2.0 * np.log(2.0)) * (1.0 + self.red_tail_ratio))
        sigma_red = sigma_blue * self.red_tail_ratio
        sigma = np.where(lam < self.peak_nm, sigma_blue, sigma_red)
        return np.exp(-0.5 * ((lam - self.peak_nm) / sigma) ** 2)

    @classmethod
    def from_csv(cls, name, path, wavelength_col=0, intensity_col=1, heterogeneity_nm=0.0, **read_csv_kwargs):
        """
        Builds a tabulated DyeSpectrum from a vendor spectrum CSV (e.g. the
        `Alexa_Fluor_647.csv` / `CF 680.csv` files make_first_order.m read,
        which kept wavelength in column 0 and emission in column 2 ->
        intensity_col=2). Columns are positional so header naming doesn't
        matter; pass header=None via read_csv_kwargs for a headerless file.
        """
        table = pd.read_csv(path, **read_csv_kwargs)
        lam = table.iloc[:, wavelength_col].to_numpy(dtype=float)
        inten = table.iloc[:, intensity_col].to_numpy(dtype=float)
        keep = np.isfinite(lam) & np.isfinite(inten)
        order = np.argsort(lam[keep])
        lam, inten = lam[keep][order], inten[keep][order]
        return cls(
            name=name, peak_nm=float(lam[np.argmax(inten)]),
            heterogeneity_nm=heterogeneity_nm,
            wavelengths_nm=tuple(lam.tolist()), intensities=tuple(inten.tolist()),
        )


# Presets. Peak / width / tail are chosen so the centroid over the default
# SpectralConfig window reproduces the measured RainbowSTORM centroids in
# hmt_v1/Simulation/blink_stats (AF647_*.csv: 670.2-670.5 nm; cf680_*.csv:
# 699.0-700.8 nm), i.e. what the real instrument reports for each dye alone,
# and sit inside the legacy splitmeac.m gates (655.6-688.4 / 704.7-737.4 nm
# were that instrument's gates; the CF680 gate there is centred ~20 nm
# redder than these measurements, so re-derive gates from a single-dye
# control rather than reusing those numbers).
AF647 = DyeSpectrum("AF647", peak_nm=665.5, fwhm_nm=42.0, red_tail_ratio=1.4)
CF680 = DyeSpectrum("CF680", peak_nm=695.0, fwhm_nm=45.0, red_tail_ratio=1.4)


@dataclass
class SpectralConfig:
    """
    The sSMLM instrument's spectral path. One instance is shared by every
    channel of an acquisition (it is hardware, not a per-dye property).

    Calibration (see module docstring):
        dispersion_nm_per_px:    grating dispersion on the camera. 6 nm/px
                                 is the value the legacy master2.m header
                                 declared (`SD=6`).
        reference_wavelength_nm: the wavelength that lands
                                 reference_offset_px to the right of the
                                 emitter's zeroth-order x position.
        reference_offset_px:     that landing offset. 0 (as in the legacy
                                 `dif-950` construction) puts the reference
                                 wavelength directly over the spatial
                                 position, which keeps a shared spatial /
                                 spectral frame compact; a real split-camera
                                 layout is a large constant here.
        wavelength_window_nm:    band over which spectra are rendered and
                                 ground-truth centroids are computed --
                                 set it to the range RainbowSTORM will be
                                 told to integrate over, so the recorded
                                 `emission_centroid_nm` is the quantity it
                                 should recover.
        kernel_tail_cutoff:      fraction of a dye's spectrum trimmed from
                                 each end of its dispersed kernel (pure
                                 speed / frame-size economy).

    First-order PSF: isotropic Gaussian with a symmetric defocus curve
    sigma1(z) = psf_sigma0_nm * sqrt(1 + ((z - psf_focal_offset_nm) / psf_depth_scale_nm)^2).
    No astigmatism: the cylindrical lens sits in the zeroth-order path.
    psf_sigma0_nm defaults to the ~200 nm first-order sigma ThunderSTORM
    fits on the real instrument (blink_stats/*_first.csv, mean 198 nm)
    vs ~145 nm for the zeroth order. psf_focal_offset_nm != 0 models a
    biplane-style focal offset between the two paths.

    Photon split: an emitter's photons that frame times
    zeroth_order_fraction go to the spatial image and times
    first_order_fraction to the spectral one. Defaults (0.45 / 0.55) are
    the measured median zeroth:first intensity ratio in blink_stats
    (456:539 photons). The legacy pipeline used 1.0 / 1.0 (full brightness
    in both orders); set both to 1.0 to reproduce that. Note the
    single-channel pipeline.run_pipeline (no grating) renders the full
    photon count and is unaffected by these.

    background_ratio: first-order background per pixel relative to the
    zeroth-order CameraConfig.background_photons_per_pixel (legacy
    add_noise.m: 1.5; measured offset ratio in blink_stats: ~1.3).
    """
    dispersion_nm_per_px: float = 6.0
    reference_wavelength_nm: float = 670.0
    reference_offset_px: float = 0.0
    wavelength_window_nm: tuple = (620.0, 800.0)
    kernel_tail_cutoff: float = 1e-4
    psf_sigma0_nm: float = 200.0
    psf_depth_scale_nm: float = 400.0
    psf_focal_offset_nm: float = 0.0
    zeroth_order_fraction: float = 0.45
    first_order_fraction: float = 0.55
    background_ratio: float = 1.3


def wavelength_to_offset_px(wavelength_nm, config):
    """x offset (px) of `wavelength_nm` in the spectral image, relative to the emitter's zeroth-order x. Vectorized."""
    lam = np.asarray(wavelength_nm, dtype=float)
    return config.reference_offset_px + (lam - config.reference_wavelength_nm) / config.dispersion_nm_per_px


def offset_px_to_wavelength(offset_px, config):
    """Inverse of wavelength_to_offset_px -- what a spectral-image column at that offset means in nm. Vectorized."""
    off = np.asarray(offset_px, dtype=float)
    return config.reference_wavelength_nm + (off - config.reference_offset_px) * config.dispersion_nm_per_px


def sample_emission(dye, config, step_nm=0.1):
    """
    The dye's emission spectrum sampled on a fine wavelength grid over
    config.wavelength_window_nm, normalized so sum(S) * step_nm == 1
    (i.e. S is a probability density over the window).

    Returns:
        (wavelengths_nm, density): two (n,) float arrays
    """
    lo, hi = config.wavelength_window_nm
    lam = np.arange(lo, hi + 0.5 * step_nm, step_nm)
    s = dye.emission(lam)
    total = s.sum() * step_nm
    if total <= 0:
        raise ValueError(
            f"dye {dye.name!r} has no emission inside wavelength_window_nm={config.wavelength_window_nm}"
        )
    return lam, s / total


def emission_centroid_nm(dye, config):
    """
    Intensity-weighted mean wavelength of the dye over
    config.wavelength_window_nm -- the ground-truth spectral centroid a
    RainbowSTORM analysis integrating that same window should recover for
    an unshifted molecule of this dye (before the emitter's own
    heterogeneity shift, which adds to it).
    """
    lam, dens = sample_emission(dye, config)
    return float((lam * dens).sum() / dens.sum())


def dispersed_kernel(dye, config):
    """
    The dye's emission spectrum binned onto the spectral image's pixel
    grid: weight w_k = fraction of the dye's photons whose dispersed
    position falls in the pixel bin at integer offset k, where the offset
    is measured from the emitter's zeroth-order x PLUS reference_offset_px
    (that reference shift, and any per-molecule spectral shift, are applied
    as a fractional shift of the PSF center at render time rather than
    baked into the kernel, so they stay exact -- only the spectrum's own
    shape is quantized to whole pixels, and that is then smoothed by the
    PSF). This is make_first_order.m's `sptimg` (the spectrum resampled
    onto the dispersion axis, normalized to sum 1), minus its missing
    calibration file.

    Returns:
        (offsets, weights): contiguous int array k_min..k_max and matching
        float weights summing to 1, trimmed of a config.kernel_tail_cutoff
        fraction at each end
    """
    lam, dens = sample_emission(dye, config)
    u = (lam - config.reference_wavelength_nm) / config.dispersion_nm_per_px
    k = np.rint(u).astype(np.int64)
    k_min, k_max = int(k.min()), int(k.max())
    weights = np.bincount(k - k_min, weights=dens, minlength=k_max - k_min + 1)
    weights = weights / weights.sum()

    cum = np.cumsum(weights)
    lo = int(np.searchsorted(cum, config.kernel_tail_cutoff))
    hi = int(np.searchsorted(cum, 1.0 - config.kernel_tail_cutoff)) + 1
    weights = weights[lo:hi]
    offsets = np.arange(k_min + lo, k_min + hi, dtype=np.int64)
    return offsets, weights / weights.sum()


def build_kernels(dyes, config):
    """{dye.name: dispersed_kernel(dye, config)} for every dye in `dyes`."""
    return {dye.name: dispersed_kernel(dye, config) for dye in dyes}


def kernel_extent_px(kernels, config):
    """
    (min, max) x offset in px, relative to an emitter's zeroth-order
    position, that any of the given kernels' spectra can reach --
    including reference_offset_px. Used to size a frame that fits the
    spectra (fit_frame_to_spatial_and_spectral) and to define the column
    window measure_spectral_centroids integrates over.
    """
    lo = min(int(off[0]) for off, _ in kernels.values()) + config.reference_offset_px
    hi = max(int(off[-1]) for off, _ in kernels.values()) + config.reference_offset_px
    return lo, hi


def first_order_sigma_nm(z_nm, config):
    """First-order PSF sigma (nm) at each z, per SpectralConfig's symmetric defocus curve. Vectorized."""
    u = (np.asarray(z_nm, dtype=float) - config.psf_focal_offset_nm) / config.psf_depth_scale_nm
    return config.psf_sigma0_nm * np.sqrt(1.0 + u**2)


def fit_frame_to_spatial_and_spectral(x_nm, y_nm, kernels, camera_config, spectral_config,
                                      render_config, margin_nm=300.0):
    """
    render.fit_frame_to_points over the UNION of the emitters' zeroth-order
    positions and the full x extent their dispersed spectra reach (kernel
    extent plus the first-order PSF's own window), so ONE frame shape
    serves both the spatial and the spectral stack with nothing clipped.
    Both stacks having identical dimensions mirrors the legacy pipeline
    (make_real_fist's output was cropped to the zeroth-order size) and a
    split-camera acquisition; the spatial stack simply has empty
    background where the other stack carries spectra.

    Returns:
        (frame_shape, x_offset_nm, y_offset_nm) exactly as
        render.fit_frame_to_points
    """
    x_nm = np.asarray(x_nm, dtype=float)
    y_nm = np.asarray(y_nm, dtype=float)
    px = camera_config.pixel_size_nm
    off_lo, off_hi = kernel_extent_px(kernels, spectral_config)
    psf_reach_nm = render_config.window_sigma_multiple * spectral_config.psf_sigma0_nm
    x_all = np.concatenate([x_nm, x_nm + off_lo * px - psf_reach_nm, x_nm + off_hi * px + psf_reach_nm])
    y_all = np.concatenate([y_nm, y_nm, y_nm])
    return srev_render.fit_frame_to_points(x_all, y_all, px, margin_nm)


def render_spectral_frame(x_nm, y_nm, z_nm, photons, dye_names, spectral_shift_nm, kernels,
                          frame_shape, camera_config, spectral_config, render_config):
    """
    One frame's expected-photon first-order image (no noise, no
    background): each active emitter's first-order PSF, dispersed along +x
    by its dye's binned spectrum. Per emitter this is one pixel-integrated
    Gaussian x-profile (render._pixel_integrated_axis_profile) convolved
    with the dye's kernel, outer-producted with the y-profile -- a
    separable version of make_real_fist.m's conv2(h, first_ord'), costing
    O(window + kernel) instead of a 2D convolution over the whole frame per
    emitter.

    Args:
        x_nm, y_nm, z_nm:   (n,) true emitter positions, frame-relative
        photons:            (n,) photons this frame ALREADY multiplied by
                            spectral_config.first_order_fraction (the
                            caller owns the split)
        dye_names:          (n,) dye name per emitter, keys of `kernels`
        spectral_shift_nm:  (n,) per-emitter spectral heterogeneity shift
        kernels:            build_kernels output
        frame_shape:        (height, width) px
        camera_config:      CameraConfig (pixel_size_nm)
        spectral_config:    SpectralConfig
        render_config:      RenderConfig (window_sigma_multiple only)

    Returns:
        (H, W) float64 expected photons per pixel
    """
    H, W = frame_shape
    buf = np.zeros((H, W), dtype=np.float64)
    if len(x_nm) == 0:
        return buf

    px = camera_config.pixel_size_nm
    sigma_px = first_order_sigma_nm(z_nm, spectral_config) / px
    cx_px = (np.asarray(x_nm, dtype=float) / px
             + spectral_config.reference_offset_px
             + np.asarray(spectral_shift_nm, dtype=float) / spectral_config.dispersion_nm_per_px)
    cy_px = np.asarray(y_nm, dtype=float) / px
    dye_names = np.asarray(dye_names)

    for i in range(len(x_nm)):
        offsets, weights = kernels[dye_names[i]]
        half_window = int(np.ceil(render_config.window_sigma_multiple * sigma_px[i])) + 1
        x_center, y_center = int(np.floor(cx_px[i])), int(np.floor(cy_px[i]))
        x0, y0 = x_center - half_window, y_center - half_window
        xs = np.arange(x0, x0 + 2 * half_window + 1, dtype=float)
        ys = np.arange(y0, y0 + 2 * half_window + 1, dtype=float)
        profile_x = srev_render._pixel_integrated_axis_profile(xs, cx_px[i], sigma_px[i])
        profile_y = srev_render._pixel_integrated_axis_profile(ys, cy_px[i], sigma_px[i])
        # full convolution: element m sits at pixel x0 + offsets[0] + m
        dispersed_x = np.convolve(profile_x, weights)
        patch = np.outer(profile_y, dispersed_x) * photons[i]
        srev_render._add_patch(buf, y0, x0 + int(offsets[0]), patch)

    return buf


def render_dual_stack_to_tiff(fluorophores_df, blink_events_df, n_frames, camera_config,
                              render_config, spectral_config, dyes, spatial_path, spectral_path,
                              rng=None):
    """
    Renders the zeroth-order (spatial) and first-order (spectral) movies of
    ONE mixed-dye emitter set and streams each to its own BigTIFF, frame
    by frame, in lockstep -- frame i of spatial_path and frame i of
    spectral_path are the same exposure, which is what RainbowSTORM pairs
    up. Both stacks use camera_config.frame_shape (size it with
    fit_frame_to_spatial_and_spectral so the spectra fit).

    The spatial frame is render.render_photon_frame (astigmatic PSF, so a
    ThunderSTORM 3D fit of it still works) at photons *
    zeroth_order_fraction; the spectral frame is render_spectral_frame at
    photons * first_order_fraction. Each gets render.apply_camera_noise
    independently (two detectors / two camera halves see independent shot
    and read noise), the spectral one with its background scaled by
    background_ratio -- add_noise.m's `Ibg = Ibgp * 1.5`.

    Args:
        fluorophores_df:  merged fluorophores with "dye" and
                          "spectral_shift_nm" columns (pipeline.merge_channels)
        blink_events_df:  merged blink events keyed to fluorophores_df's
                          fluorophore_id
        n_frames:         frames to write (both stacks)
        camera_config:    CameraConfig with frame_shape set
        render_config:    RenderConfig (spatial PSF)
        spectral_config:  SpectralConfig
        dyes:             the DyeSpectrum objects present in fluorophores_df
        spatial_path, spectral_path: output .tif paths
        rng:              np.random.Generator (default_rng() if None)

    Returns:
        the kernels dict used (build_kernels output), for QC
    """
    if srev_render.tifffile is None:
        raise ImportError("render_dual_stack_to_tiff requires tifffile (pip install tifffile)")
    rng = rng if rng is not None else np.random.default_rng()

    kernels = build_kernels(dyes, spectral_config)
    missing = set(fluorophores_df["dye"].unique()) - set(kernels)
    if missing:
        raise ValueError(f"fluorophores_df has dyes with no DyeSpectrum passed in `dyes`: {sorted(missing)}")

    active = srev_render.expand_active_frames(
        fluorophores_df, blink_events_df, extra_columns=("dye", "spectral_shift_nm")
    )
    grouped = active.groupby("frame") if len(active) else None
    active_frames = set(grouped.groups.keys()) if grouped is not None else set()
    frame_shape = camera_config.frame_shape
    spectral_camera = replace(
        camera_config,
        background_photons_per_pixel=camera_config.background_photons_per_pixel * spectral_config.background_ratio,
    )

    tifffile = srev_render.tifffile
    with tifffile.TiffWriter(spatial_path, bigtiff=True) as tif_spatial, \
         tifffile.TiffWriter(spectral_path, bigtiff=True) as tif_spectral:
        for frame in range(n_frames):
            if frame in active_frames:
                g = grouped.get_group(frame)
                x, y, z = g["x [nm]"].to_numpy(), g["y [nm]"].to_numpy(), g["z [nm]"].to_numpy()
                rate = g["photon_rate"].to_numpy()
                spatial_img = srev_render.render_photon_frame(
                    x, y, z, rate * spectral_config.zeroth_order_fraction,
                    frame_shape, camera_config, render_config,
                )
                spectral_img = render_spectral_frame(
                    x, y, z, rate * spectral_config.first_order_fraction,
                    g["dye"].to_numpy(), g["spectral_shift_nm"].to_numpy(), kernels,
                    frame_shape, camera_config, spectral_config, render_config,
                )
            else:
                spatial_img = np.zeros(frame_shape, dtype=np.float64)
                spectral_img = np.zeros(frame_shape, dtype=np.float64)

            tif_spatial.write(srev_render.apply_camera_noise(spatial_img, camera_config, rng), contiguous=True)
            tif_spectral.write(srev_render.apply_camera_noise(spectral_img, spectral_camera, rng), contiguous=True)

    return kernels


def render_calibration_pair(x_nm, y_nm, line_wavelengths_nm, photons_per_line, frame_shape,
                            camera_config, render_config, spectral_config, rng=None, z_nm=0.0,
                            output_dir=None):
    """
    A single-frame spatial / spectral image pair of a spectral CALIBRATION
    target: point sources (pinholes) at (x_nm, y_nm) in the spatial frame,
    and for each of them one monochromatic line per wavelength in
    line_wavelengths_nm in the spectral frame -- what a pinhole array lit
    by a line lamp or a set of narrow bandpass filters looks like on an
    sSMLM system, and the input RainbowSTORM's spectral-calibration step
    fits its pixel->wavelength relation from. Rendered with the SAME
    forward model as the data (wavelength_to_offset_px, first-order PSF,
    apply_camera_noise), so a calibration fitted from this pair is valid
    for stacks from render_dual_stack_to_tiff with the same SpectralConfig.

    Args:
        x_nm, y_nm:          (n,) pinhole positions, frame-relative
        line_wavelengths_nm: the known wavelengths, e.g. (620, 660, 700, 740, 780)
        photons_per_line:    photons in each line (and in each pinhole's
                             zeroth-order spot)
        frame_shape:         (height, width) px, ideally the data stack's
        camera_config, render_config, spectral_config: as for the data
        rng:                 np.random.Generator (default_rng() if None)
        z_nm:                pinhole z (0 = in focus)
        output_dir:          if given, also writes calibration_spatial.tif
                             and calibration_spectral.tif there

    Returns:
        (spatial_adu, spectral_adu): two (H, W) uint16 frames
    """
    rng = rng if rng is not None else np.random.default_rng()
    x_nm = np.asarray(x_nm, dtype=float)
    y_nm = np.asarray(y_nm, dtype=float)
    z = np.full(len(x_nm), float(z_nm))
    photons = np.full(len(x_nm), float(photons_per_line))

    spatial_img = srev_render.render_photon_frame(x_nm, y_nm, z, photons, frame_shape, camera_config, render_config)

    # one delta-function "dye" per line: kernel weight 1 at offset 0, and the
    # line's wavelength enters as a spectral shift relative to the reference
    line_wavelengths_nm = np.asarray(line_wavelengths_nm, dtype=float)
    kernels = {"line": (np.array([0], dtype=np.int64), np.array([1.0]))}
    x_rep = np.tile(x_nm, len(line_wavelengths_nm))
    y_rep = np.tile(y_nm, len(line_wavelengths_nm))
    z_rep = np.tile(z, len(line_wavelengths_nm))
    shift_rep = np.repeat(line_wavelengths_nm - spectral_config.reference_wavelength_nm, len(x_nm))
    spectral_img = render_spectral_frame(
        x_rep, y_rep, z_rep, np.full(len(x_rep), float(photons_per_line)),
        np.full(len(x_rep), "line"), shift_rep, kernels,
        frame_shape, camera_config, spectral_config, render_config,
    )

    spectral_camera = replace(
        camera_config,
        background_photons_per_pixel=camera_config.background_photons_per_pixel * spectral_config.background_ratio,
    )
    spatial_adu = srev_render.apply_camera_noise(spatial_img, camera_config, rng)
    spectral_adu = srev_render.apply_camera_noise(spectral_img, spectral_camera, rng)

    if output_dir is not None:
        if srev_render.tifffile is None:
            raise ImportError("writing calibration TIFFs requires tifffile (pip install tifffile)")
        from pathlib import Path
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        srev_render.tifffile.imwrite(output_dir / "calibration_spatial.tif", spatial_adu)
        srev_render.tifffile.imwrite(output_dir / "calibration_spectral.tif", spectral_adu)

    return spatial_adu, spectral_adu


def measure_spectral_centroids(spectral_tiff_path, frames, x_nm, y_nm, camera_config, spectral_config,
                               row_half_width_px=2, column_window_px=None, background="frame_median"):
    """
    A deliberately minimal stand-in for RainbowSTORM's centroid step, for
    QC of a rendered stack against its ground truth WITHOUT leaving Python
    -- not a replacement for running RainbowSTORM (no spatial fitting, no
    spectral background model, no crosstalk handling). For each given
    localization it sums the spectral frame over rows y +/- row_half_width_px
    and the columns the detection window maps to, subtracts a background
    level, and returns the intensity-weighted mean wavelength.

    Args:
        spectral_tiff_path: render_dual_stack_to_tiff's spectral_path
        frames:             (n,) frame index per localization
        x_nm, y_nm:         (n,) localization positions, frame-relative
                            (use the ground-truth x/y, or ThunderSTORM's
                            spatial-stack fit)
        camera_config:      CameraConfig used to render
        spectral_config:    SpectralConfig used to render
        row_half_width_px:  rows summed on each side of y
        column_window_px:   (lo, hi) x-offset range in px relative to x
                            (kernel_extent_px output); default = the full
                            wavelength_window_nm mapped through
                            wavelength_to_offset_px
        background:         "frame_median" (the frame's median ADU, robust
                            for sparse frames but biased once streaks
                            cover much of a small frame), "camera" (the
                            exact expected background ADU from
                            camera_config / spectral_config, i.e. what
                            was rendered), or a number in ADU

    Returns:
        (centroid_nm, integrated_adu): two (n,) float arrays; NaN centroid
        where the window fell entirely outside the frame or had no signal
        above background
    """
    if srev_render.tifffile is None:
        raise ImportError("measure_spectral_centroids requires tifffile (pip install tifffile)")
    frames = np.asarray(frames, dtype=np.int64)
    x_px = np.asarray(x_nm, dtype=float) / camera_config.pixel_size_nm
    y_px = np.asarray(y_nm, dtype=float) / camera_config.pixel_size_nm
    if column_window_px is None:
        lo_nm, hi_nm = spectral_config.wavelength_window_nm
        column_window_px = (float(wavelength_to_offset_px(lo_nm, spectral_config)),
                            float(wavelength_to_offset_px(hi_nm, spectral_config)))
    off_lo, off_hi = column_window_px
    if background == "camera":
        expected_bg_e = (camera_config.background_photons_per_pixel * spectral_config.background_ratio
                         * camera_config.quantum_efficiency * max(camera_config.em_gain, 1.0))
        background = camera_config.base_level_adu + expected_bg_e / camera_config.e_per_adu

    centroid = np.full(len(frames), np.nan)
    integrated = np.full(len(frames), np.nan)
    order = np.argsort(frames, kind="stable")

    with srev_render.tifffile.TiffFile(spectral_tiff_path) as tif:
        n_pages = len(tif.pages)
        current_frame, img, bg = -1, None, 0.0
        for i in order:
            f = int(frames[i])
            if f != current_frame:
                if f < 0 or f >= n_pages:
                    continue
                img = tif.pages[f].asarray().astype(np.float64)
                bg = float(np.median(img)) if background == "frame_median" else float(background)
                current_frame = f
            H, W = img.shape
            r0 = max(int(np.floor(y_px[i])) - row_half_width_px, 0)
            r1 = min(int(np.floor(y_px[i])) + row_half_width_px + 1, H)
            c0 = max(int(np.floor(x_px[i] + off_lo)), 0)
            c1 = min(int(np.ceil(x_px[i] + off_hi)) + 1, W)
            if r0 >= r1 or c0 >= c1:
                continue
            # No clipping at zero: negative noise residuals must be kept so they
            # cancel their positive counterparts -- clipping leaves a positive
            # pedestal across the window that pulls every centroid toward the
            # window's center wavelength (measured: +5 nm on AF647 with the
            # default 620-800 nm window).
            spectrum = img[r0:r1, c0:c1].sum(axis=0) - bg * (r1 - r0)
            total = spectrum.sum()
            if total <= 0:
                continue
            cols = np.arange(c0, c1, dtype=float)
            lam = offset_px_to_wavelength(cols - x_px[i], spectral_config)
            centroid[i] = float((lam * spectrum).sum() / total)
            integrated[i] = float(total)

    return centroid, integrated
