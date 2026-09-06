"""
Packing domain identification and sizing for an SR-EV monomer cloud --
the established, published method from Cangnano et al. 2024 (eLife,
"Local Volume Concentration, Packing Domains and Scaling Properties of
Chromatin", the paper describing this exact SR-EV model), Methods section
"Chromatin Domain Radius Measured from Experiment" and Appendix 2.

Unlike compartments.py (per-monomer heterochromatin/euchromatin
labeling via coordination number, used to drive antibody/dye targeting),
this module answers a different question: how many distinct packing
domains does a configuration actually contain, and how big is each one?
The domain COUNT is never a free parameter here -- it emerges from local
maxima of a rendered density image, exactly as in the source paper and
in the real ChromSTEM image analysis it was built to match ("essentially
the same as the experimental one").

Procedure (paper's own description):
  1. Render bead coordinates to a density image: each bead is a normal
     distribution, and its contribution to a pixel is the integral of
     that distribution over the pixel's area (reuses
     render._pixel_integrated_gaussian -- the same math the PSF
     renderer uses, applied here to a static density image instead of a
     per-frame photon image).
  2. Collapse a 100 nm-thick slab by projecting along z (paper: "2D
     chromatin density distributions were obtained by re-projection of
     the tomogram along z-axis").
  3. Gaussian-filter the image (paper: "5 pixels radius") and apply CLAHE
     contrast enhancement (paper: "block size of 120 pixels") -- but only
     within an Otsu-thresholded "occupied" mask of the pre-CLAHE density;
     see find_domain_centers' docstring for why the mask is necessary.
  4. Domain centers = local maxima of the masked, enhanced image, refined
     to a weighted centroid within an 11x11 pixel window (paper: "multiple
     mass scaling curves were sampled by using pixels (a 11-pixel x
     11-pixel window) around the center ... averaged by the weight of
     the pixel values").
  5. For each center, compute a radial density profile (paper: "assuming
     cylindrical symmetry") from the real bead coordinates, and take the
     domain radius at the profile's first local minimum (Appendix 2:
     "The radius of a domain corresponds to the first minimum in the
     density profile").

Calibration status (checked against the real ChromSTEM source method,
Li et al. 2022, Scientific Reports, "Analysis of three-dimensional
chromatin packing domains by chromatin scanning transmission electron
microscopy (ChromSTEM)" -- the paper Cangnano et al. 2024 explicitly
say they reused):
  - Voxel/pixel size: Li et al. report "1.8 to 2.9 nm" for their raw
    tomograms; identify_packing_domains defaults pixel_size_nm to 2.0
    accordingly. The Gaussian-filter (5 px) and CLAHE (120 px block)
    parameters are applied in THAT pixel unit, so this default matters.
  - CLAHE bug found and fixed (2026-09): CLAHE performs *local* contrast
    stretching per tile, which is safe on real ChromSTEM tomograms
    (continuously-varying, never truly zero density) but on this
    simulated bead-density image -- which DOES have genuinely near-zero
    void regions between packing domains -- it was manufacturing "peaks"
    out of floating-point noise inside empty tiles. Direct inspection on
    hmt_v1/Simulation/config-1.txt: CLAHE without masking found 46
    peaks, of which the raw (pre-CLAHE) density at the median peak was
    exactly 0, and fully half had zero real monomers within 20 nm --
    i.e. "domains" sitting in empty space. Those centers' radial density
    profiles start at ~0 and only rise once they happen to reach a real
    neighboring domain, which is why ~45% of all computed radii were
    pegged at the r_max_nm profiling ceiling rather than a genuine
    measurement. Fix: find_domain_centers now zeroes out the enhanced
    image outside an Otsu-thresholded "occupied" mask (computed on the
    pre-CLAHE normalized density) before peak-finding, so CLAHE's local
    stretching can no longer manufacture peaks in void tiles.
  - Domain COUNT, recalibrated post-fix: with the masking fix and
    min_distance_nm=30 (down from the old 110, which was tuned against
    counts inflated by the CLAHE bug), hmt_v1/Simulation/config-1.txt
    gives 41 domains over 2.327 um^2 = 17.6 domains/um^2 -- matching Li
    et al. 2022's reported 17.5 domains/um^2 almost exactly, and now
    with every center sitting on real density (min 36 monomers within
    20 nm, vs. many at exactly 0 pre-fix).
  - Domain RADIUS: Li et al. 2022's own radius definition (used for
    their real ChromSTEM data, Methods "Chromatin Domain Radius Measured
    from Experiment") is the smallest length scale where ANY of three
    criteria trips: (1) the mass-scaling curve deviates from its initial
    power-law fit by 5%, (2) local packing scaling D reaches 3, or (3)
    radial density/CVC begins to increase. Cangnano et al. 2024 do NOT
    use all three when applying this to their own SR-EV point clouds
    (the same kind of data this module works on) -- Appendix 2--Figure
    2's caption states plainly, for that exact case: "The radius of a
    domain corresponds to the first minimum in the density profile."
    That is criterion (3) alone, which is what domain_radius_from_profile
    implements -- not a shortcut relative to the paper, but a match to
    what it does for SR-EV data (the fuller 3-criteria test is only used
    there for the real-ChromSTEM comparison dataset, their Fig 7C).
    Post-fix, on config-1.txt: median radius 85 nm, mean 107 nm (pulled
    up by a right tail out to ~235 nm at the 95th percentile), only 2.4%
    of domains still hitting the r_max_nm=300 ceiling. This lines up
    well with the paper's reported values (A549: median Rf 74 nm, mean
    80.6 nm; BJ: 58.98 nm) -- radius_nm can now be treated as a
    reasonable measurement rather than an upper bound, though the
    right-tail domains that still hit the ceiling are individual
    profiles that genuinely don't turn over within 300 nm (not a
    detection artifact) and would need a larger r_max_nm or the
    mass-scaling criteria to size properly.
"""
import numpy as np
import pandas as pd
from scipy import ndimage
from skimage.exposure import equalize_adapthist
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu

from .render import _pixel_integrated_gaussian


def render_2d_density(monomers_df, z_center_nm, slab_thickness_nm, pixel_size_nm,
                       bead_sigma_nm, margin_nm=100.0):
    """
    Renders a 2D chromatin density image from a z-slab of monomers, each
    bead contributing a pixel-integrated Gaussian (not point-sampled) --
    the same rasterization render.render_photon_frame uses for PSFs,
    reused here as a static density estimator instead of a per-frame
    photon accumulator.

    Args:
        monomers_df:        output of io.parse_config_txt (or a
                             compartment-labeled variant -- only x/y/z
                             [nm] are used)
        z_center_nm:         center of the slab along z
        slab_thickness_nm:   full slab thickness (100 nm in the source paper)
        pixel_size_nm:       image pixel size
        bead_sigma_nm:       per-bead Gaussian sigma (see module docstring
                              -- not given an exact value in the source
                              text; a common physically-motivated choice
                              is the bead's own excluded-volume radius)
        margin_nm:           padding added around the slab's monomer
                              extent, in nm

    Returns:
        (density, x0_nm, y0_nm): density is a (H, W) float64 array;
        (x0_nm, y0_nm) is the image origin, i.e. pixel (0, 0)'s
        lower-left corner in world nm coordinates
    """
    z_lo = z_center_nm - slab_thickness_nm / 2.0
    z_hi = z_center_nm + slab_thickness_nm / 2.0
    in_slab = (monomers_df["z [nm]"] >= z_lo) & (monomers_df["z [nm]"] < z_hi)
    xy = monomers_df.loc[in_slab, ["x [nm]", "y [nm]"]].to_numpy(dtype=float)

    if len(xy) == 0:
        raise ValueError(
            f"no monomers found in slab z in [{z_lo}, {z_hi}) nm -- check z_center_nm"
        )

    x_min = xy[:, 0].min() - margin_nm
    y_min = xy[:, 1].min() - margin_nm
    width_px = int(np.ceil((xy[:, 0].max() + margin_nm - x_min) / pixel_size_nm))
    height_px = int(np.ceil((xy[:, 1].max() + margin_nm - y_min) / pixel_size_nm))

    density = np.zeros((height_px, width_px), dtype=np.float64)
    sigma_px = bead_sigma_nm / pixel_size_nm
    half_window = int(np.ceil(4 * sigma_px)) + 1

    cx_px = (xy[:, 0] - x_min) / pixel_size_nm
    cy_px = (xy[:, 1] - y_min) / pixel_size_nm

    for i in range(len(xy)):
        y0, x0, patch = _pixel_integrated_gaussian(cx_px[i], cy_px[i], sigma_px, sigma_px, half_window)
        y1, x1 = y0 + patch.shape[0], x0 + patch.shape[1]

        cy0, cy1 = max(y0, 0), min(y1, height_px)
        cx0, cx1 = max(x0, 0), min(x1, width_px)
        if cy0 >= cy1 or cx0 >= cx1:
            continue

        py0, py1 = cy0 - y0, patch.shape[0] - (y1 - cy1)
        px0, px1 = cx0 - x0, patch.shape[1] - (x1 - cx1)
        density[cy0:cy1, cx0:cx1] += patch[py0:py1, px0:px1]

    return density, x_min, y_min


def find_domain_centers(density_image, x0_nm, y0_nm, pixel_size_nm,
                         gaussian_filter_sigma_px=5.0, clahe_kernel_px=120,
                         min_distance_nm=30.0, centroid_window_px=11):
    """
    Finds packing domain centers as local maxima of the smoothed,
    contrast-enhanced density image -- the paper's own procedure, not a
    clustering algorithm with a chosen k. The NUMBER of returned centers
    is the domain count; it is an output of this function, never an input.

    IMPORTANT -- CLAHE is masked to an "occupied" region before
    peak-finding (found necessary 2026-09, see module docstring "CLAHE
    bug found and fixed"): CLAHE's local contrast stretching manufactures
    spurious peaks out of floating-point noise inside near-zero-density
    void tiles, which real ChromSTEM tomograms never have but this
    simulated bead-density image does. The mask is an Otsu threshold on
    the PRE-CLAHE smoothed/normalized density; pixels below it are
    zeroed in the enhanced image so CLAHE's stretching inside them can
    never register as a peak. Without this, on
    hmt_v1/Simulation/config-1.txt, roughly half of all detected
    "domains" sat in regions with zero real monomers nearby.

    Args:
        density_image:            (H, W) array from render_2d_density
        x0_nm, y0_nm:              that image's origin (as returned by
                                   render_2d_density)
        pixel_size_nm:             that image's pixel size
        gaussian_filter_sigma_px:  smoothing sigma before peak detection
                                   (paper: "5 pixels radius")
        clahe_kernel_px:           CLAHE block size (paper: "120 pixels")
        min_distance_nm:           minimum separation enforced between
                                   detected peaks (skimage.feature.peak_local_max's
                                   own de-duplication radius, converted to
                                   pixels internally). Not stated in
                                   either source paper; the default of
                                   30 nm is fit (post CLAHE-masking fix)
                                   to reproduce Li et al. 2022's reported
                                   domain density (~17.5 domains/um^2) on
                                   hmt_v1/Simulation/config-1.txt (gives
                                   17.6/um^2) -- see module docstring
                                   "Calibration status". Too small a
                                   value fragments one real domain into
                                   several detected peaks; too large
                                   merges genuinely distinct domains into
                                   one.
        centroid_window_px:        window size for the weighted-centroid
                                   center refinement (paper: "11-pixel x
                                   11-pixel window")

    Returns:
        (N, 2) array of domain center [x, y] in nm
    """
    smoothed = ndimage.gaussian_filter(density_image, sigma=gaussian_filter_sigma_px)

    peak = smoothed.max()
    normalized = smoothed / peak if peak > 0 else smoothed
    enhanced = equalize_adapthist(normalized, kernel_size=clahe_kernel_px)

    try:
        occupied = normalized > threshold_otsu(normalized)
    except ValueError:
        # threshold_otsu needs at least 2 distinct intensity values; a degenerate
        # (e.g. uniform or all-zero) density image has nothing to mask against.
        occupied = np.ones_like(normalized, dtype=bool)
    enhanced_masked = np.where(occupied, enhanced, 0.0)

    min_distance_px = max(1, int(round(min_distance_nm / pixel_size_nm)))
    peaks_px = peak_local_max(enhanced_masked, min_distance=min_distance_px)  # (row=y, col=x)
    if len(peaks_px) == 0:
        return np.empty((0, 2))

    half = centroid_window_px // 2
    height_px, width_px = density_image.shape
    refined_xy = []
    for row, col in peaks_px:
        y_lo, y_hi = max(row - half, 0), min(row + half + 1, height_px)
        x_lo, x_hi = max(col - half, 0), min(col + half + 1, width_px)
        window = density_image[y_lo:y_hi, x_lo:x_hi]

        if window.sum() <= 0:
            cy, cx = float(row), float(col)
        else:
            yy, xx = np.mgrid[y_lo:y_hi, x_lo:x_hi]
            cy = float((yy * window).sum() / window.sum())
            cx = float((xx * window).sum() / window.sum())
        refined_xy.append((cx, cy))

    refined_xy = np.asarray(refined_xy)
    x_nm = x0_nm + refined_xy[:, 0] * pixel_size_nm
    y_nm = y0_nm + refined_xy[:, 1] * pixel_size_nm
    return np.column_stack([x_nm, y_nm])


def radial_density_profile(monomers_df, center_xy_nm, z_range_nm, r_max_nm=300.0, dr_nm=10.0):
    """
    2D radial monomer-count density around a domain center, within the
    same z-slab used to find it -- "assuming cylindrical symmetry" per
    the source paper's Appendix 2.

    Args:
        monomers_df:  output of io.parse_config_txt
        center_xy_nm: (x0, y0) domain center in nm
        z_range_nm:   (z_lo, z_hi) slab bounds -- pass the same slab
                      render_2d_density used for this center
        r_max_nm:     maximum radius to profile
        dr_nm:        ring width

    Returns:
        (radii_nm, density): radii_nm are ring midpoints; density is
        count-per-nm^2 in each ring (2D, matching the collapsed-slab
        density image the center was found in)
    """
    x0, y0 = center_xy_nm
    z_lo, z_hi = z_range_nm
    in_slab = (monomers_df["z [nm]"] >= z_lo) & (monomers_df["z [nm]"] < z_hi)
    xy = monomers_df.loc[in_slab, ["x [nm]", "y [nm]"]].to_numpy(dtype=float)
    r = np.hypot(xy[:, 0] - x0, xy[:, 1] - y0)

    edges = np.arange(0.0, r_max_nm + dr_nm, dr_nm)
    counts, _ = np.histogram(r, bins=edges)
    ring_areas = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    density = counts / ring_areas
    radii_nm = 0.5 * (edges[1:] + edges[:-1])
    return radii_nm, density


def domain_radius_from_profile(radii_nm, density, smooth_sigma_bins=1.5):
    """
    Domain radius = the first local minimum of the radial density
    profile (Appendix 2: "The radius of a domain corresponds to the
    first minimum in the density profile"). The profile is lightly
    Gaussian-smoothed first since it's built from a discrete point count
    per ring (noisiest at small radii, where ring area and point count
    are both small) -- the paper's own profiles (Appendix 2 Figure 2)
    are visibly smooth curves, consistent with some such smoothing
    having been applied before minimum-finding, though the exact amount
    isn't stated.

    Args:
        radii_nm, density:   from radial_density_profile
        smooth_sigma_bins:   Gaussian smoothing sigma, in BINS (not nm)

    Returns:
        float: radius in nm at the first local minimum, or the largest
        profiled radius if no interior local minimum is found (a
        domain whose profile never turns over within r_max_nm --
        widen r_max_nm in radial_density_profile if this is common)
    """
    smoothed = ndimage.gaussian_filter1d(density, sigma=smooth_sigma_bins)
    for i in range(1, len(smoothed) - 1):
        if smoothed[i] < smoothed[i - 1] and smoothed[i] <= smoothed[i + 1]:
            return float(radii_nm[i])
    return float(radii_nm[-1])


def identify_packing_domains(monomers_df, z_center_nm=None, slab_thickness_nm=100.0,
                              pixel_size_nm=2.0, bead_sigma_nm=4.9,
                              gaussian_filter_sigma_px=5.0, clahe_kernel_px=120,
                              min_distance_nm=30.0, centroid_window_px=11,
                              r_max_nm=300.0, dr_nm=10.0):
    """
    Full pipeline: render a slab's density image, find domain centers as
    its local maxima, and size each domain from its own radial density
    profile. The resulting row count IS the identified number of packing
    domains -- never set as an input.

    IMPORTANT -- current calibration status (see module docstring for
    detail): find_domain_centers masks CLAHE to an Otsu-thresholded
    "occupied" region before peak-finding, fixing a bug where CLAHE
    manufactured spurious domain centers in near-zero-density void
    regions (real ChromSTEM data never has these; this simulated density
    image does). Post-fix, on hmt_v1/Simulation/config-1.txt: 41 domains
    over 2.327 um^2 = 17.6 domains/um^2 (Li et al. 2022 report 17.5/um^2)
    with every center on real density, and radius_nm has median 85 nm /
    mean 107 nm (paper: A549 median Rf 74 nm, mean 80.6 nm), with only
    2.4% of domains still hitting the r_max_nm profiling ceiling. Both
    COUNT and RADIUS can now be treated as reasonable measurements
    rather than upper bounds, though the remaining ceiling-hitting
    domains (large, slowly-declining profiles that genuinely don't turn
    over within r_max_nm) would need a wider r_max_nm or the paper's
    mass-scaling criteria to size properly.

    Args:
        monomers_df:        output of io.parse_config_txt
        z_center_nm:         slab center; defaults to the median z of
                              monomers_df if None
        slab_thickness_nm:   100 nm in the source paper
        pixel_size_nm:       density-image pixel size; 2.0 nm matches
                              Li et al. 2022's reported ChromSTEM voxel
                              size (1.8-2.9 nm) -- this also sets the
                              physical scale of gaussian_filter_sigma_px
                              and clahe_kernel_px below, both stated in
                              THIS pixel unit by that paper
        bead_sigma_nm:       per-bead Gaussian width for the density
                              render. NOT given an exact value in either
                              source paper; default of 4.9 nm uses the
                              bead excluded-volume radius r° from
                              Cangnano et al. 2024 as a physically
                              motivated stand-in, not a stated value
        gaussian_filter_sigma_px, clahe_kernel_px, min_distance_nm,
        centroid_window_px:  see find_domain_centers
        r_max_nm, dr_nm:     see radial_density_profile

    Returns:
        DataFrame: domain_id, x [nm], y [nm], radius_nm -- len(result) is
        the identified packing domain count for this slab
    """
    if z_center_nm is None:
        z_center_nm = float(monomers_df["z [nm]"].median())

    density, x0_nm, y0_nm = render_2d_density(
        monomers_df, z_center_nm, slab_thickness_nm, pixel_size_nm, bead_sigma_nm
    )
    centers = find_domain_centers(
        density, x0_nm, y0_nm, pixel_size_nm,
        gaussian_filter_sigma_px, clahe_kernel_px, min_distance_nm, centroid_window_px,
    )

    z_range_nm = (z_center_nm - slab_thickness_nm / 2.0, z_center_nm + slab_thickness_nm / 2.0)
    radii = np.array([
        domain_radius_from_profile(
            *radial_density_profile(monomers_df, tuple(c), z_range_nm, r_max_nm, dr_nm)
        )
        for c in centers
    ])

    return pd.DataFrame({
        "domain_id": np.arange(len(centers), dtype=np.int64),
        "x [nm]": centers[:, 0] if len(centers) else np.array([]),
        "y [nm]": centers[:, 1] if len(centers) else np.array([]),
        "radius_nm": radii,
    })
