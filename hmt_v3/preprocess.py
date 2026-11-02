"""
Pre-processing of two-mark sSMLM localizations: filtering, nucleus masking,
radial contour bands, and density subsampling.

Input is one DataFrame per histone mark (e.g. H3K27me3 and H3K27ac), already
split by spectral classification, with columns "x [nm]", "y [nm]", "z [nm]".

Mask / contour grid convention (shared with measure, simulate and calibrate):
    * binarize_nucleus histograms the combined, PRE-mask localizations on a
      bin_size grid whose origin is their minimum x and y -- get it with
      mask_origin(me3_df, ac_df) and pass that same (x_min, y_min) everywhere
      a function maps localizations onto the mask.
    * Arrays are indexed [x_bin, y_bin] (transpose before imshow).
    * create_radial_contours numbers bands 0 = innermost core ...
      num_bands - 1 = outer periphery, and -1 = outside the nucleus.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage import binary_fill_holes, distance_transform_edt, gaussian_filter, label


# --- Localization filters ------------------------------------------------------

def filter_quantile(df, column, percent=0.05):
    """
    Keeps rows whose `column` lies within its [percent, 1 - percent] quantile
    range (default: the 5th-95th percentile), trimming outliers in e.g.
    photon count or localization uncertainty.

    Args:
        df:      localization DataFrame
        column:  name of the column to trim on, e.g. "uncertainty [nm]"
        percent: fraction trimmed from EACH tail (0.05 -> keep 5th-95th pct)

    Returns:
        pd.DataFrame: the rows of df inside the quantile range (a view, not a copy)
    """
    lower_bound = df[column].quantile(percent)
    upper_bound = df[column].quantile(1 - percent)
    
    filtered_df = df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]
    
    return filtered_df

def filter_axial(df, padding=1):
    """
    Drops localizations piled up at the edges of the axial fitting range.

    Astigmatic z fitting clamps out-of-range emitters to the minimum / maximum
    calibrated z, so the extreme z values hold a spike of mis-localized points.
    Keeps only points strictly inside (min z + padding, max z - padding).

    Args:
        df:      localization DataFrame with "z [nm]"
        padding: extra margin (nm) trimmed inside each edge of the z range

    Returns:
        pd.DataFrame: the rows of df strictly inside the padded z range
    """
    lower_bound = df["z [nm]"].min() + padding
    upper_bound = df["z [nm]"].max() - padding

    filtered_df = df[(df["z [nm]"] > lower_bound) & (df["z [nm]"] < upper_bound)]

    return filtered_df

def mask_origin(*dfs):
    """
    Returns the (x_min, y_min) origin in nm of the pixel grid that binarize_nucleus
    and create_radial_contours build.

    binarize_nucleus starts its 2D histogram at the minimum x and y of the
    dataframes it is given (before any masking). Every downstream function that maps
    localizations onto the mask or the contour bands — extract_radial_density_profile,
    assign_contour_bands, place_seeds / place_seeds_matched, add_noise_locs,
    measure_pair_correlation, generate_nucleus — must use this same origin, or the
    localizations are shifted into the wrong bands (a small nm offset can be many
    pixels, which empties the core bands).

    IMPORTANT: pass the SAME dataframes you gave binarize_nucleus (the pre-mask
    channels), not the masked outputs it returns. Masking drops outlying points, so
    the masked minimum is larger than the grid origin.

    Args:
        *dfs: the pre-mask dataframes passed to binarize_nucleus, each with
              columns [x [nm], y [nm]]

    Returns:
        (x_min, y_min) in nm
    """
    comb = pd.concat(dfs)
    return float(comb["x [nm]"].min()), float(comb["y [nm]"].min())


# --- Nucleus mask and radial contours ------------------------------------------

def auto_mask_threshold(me3_df, ac_df, bin_size=50, sigma=4.0, factor=0.5, floor=0.05):
    """
    Data-driven `thresh` for binarize_nucleus: `factor` x the Otsu threshold of the
    smoothed density map (over pixels above `floor`, i.e. ignoring empty field of view).

    A fixed thresh is in absolute smoothed locs-per-pixel, so a value tuned on one
    cell fails on a cell labelled several-fold more sparsely: only the dense
    heterochromatin rim clears it, the ring doesn't close, binary_fill_holes can't
    fill the interior, and most of the nucleus is thrown away. Otsu scales with the
    cell's own density. Otsu alone sits between "dense rim" and "everything else",
    so it still drops sparse interior; half of it (factor=0.5) lands between
    "nucleus" and "background" across the cells this was checked on, and on dense
    cells agrees with the previously hand-tuned thresh=2 to within a few % of area.

    Args:
        me3_df, ac_df: the same (pre-mask) dataframes you will pass to binarize_nucleus
        bin_size, sigma: must match the binarize_nucleus call
        factor:        multiplier on the Otsu threshold
        floor:         smoothed-density cutoff below which pixels are treated as empty
                       field of view and left out of the Otsu histogram

    Returns:
        float: thresh to pass to binarize_nucleus
    """
    comb_df = pd.concat([me3_df, ac_df])
    x = comb_df["x [nm]"].values
    y = comb_df["y [nm]"].values
    x_bins = np.arange(x.min(), x.max() + bin_size, bin_size)
    y_bins = np.arange(y.min(), y.max() + bin_size, bin_size)
    density_map, _, _ = np.histogram2d(x, y, bins=(x_bins, y_bins))
    values = gaussian_filter(density_map, sigma=sigma)
    values = values[values > floor]

    hist, edges = np.histogram(values, bins=256)
    centers = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    cum_mass = np.cumsum(hist * centers)
    m0 = cum_mass / np.maximum(w0, 1)
    m1 = (cum_mass[-1] - cum_mass) / np.maximum(w1, 1)
    otsu = centers[np.argmax(w0 * w1 * (m0 - m1) ** 2)]
    return float(factor * otsu)


def binarize_nucleus(me3_df, ac_df, thresh, bin_size=50, sigma=4.0, num_nuclei=1, show_plots=False):
    """
    Builds a binary nucleus mask from the combined localization density of
    both marks, and keeps only the localizations that fall inside it.

    Both channels are pooled, histogrammed on a bin_size grid starting at their
    minimum x / y (see mask_origin), Gaussian-smoothed, and thresholded at
    `thresh`. Holes are filled (so nucleoli and other signal-free interiors
    count as inside), and only the `num_nuclei` largest connected components
    are kept.

    Args:
        me3_df, ac_df: per-mark DataFrames with "x [nm]", "y [nm]"
        thresh:        smoothed locs-per-pixel threshold for "inside"; see
                       auto_mask_threshold for a density-adaptive value
        bin_size:      pixel size of the density grid (nm)
        sigma:         Gaussian smoothing width, in pixels
        num_nuclei:    number of largest connected components to keep
        show_plots:    if True, show density map, mask, and raw vs. masked points

    Returns:
        tuple: (binary_mask, me3_masked, ac_masked)
            binary_mask: 2D bool array indexed [x_bin, y_bin]
            me3_masked, ac_masked: the rows of each input inside the mask
    """
    comb_df = pd.concat([me3_df, ac_df])
    x = comb_df["x [nm]"].values
    y = comb_df["y [nm]"].values
    
    x_bins = np.arange(x.min(), x.max() + bin_size, bin_size)
    y_bins = np.arange(y.min(), y.max() + bin_size, bin_size)
    density_map, _, _ = np.histogram2d(x, y, bins=(x_bins, y_bins))

    smoothed_map = gaussian_filter(density_map, sigma=sigma)
    binary_mask = smoothed_map > thresh
    binary_mask = binary_fill_holes(binary_mask)
    
    labeled_mask, num_features = label(binary_mask)
    if num_features > num_nuclei:
        component_sizes = np.bincount(labeled_mask.ravel())
        component_sizes[0] = 0  # Set background size to 0 so it's ignored
        
        largest_labels = np.argsort(component_sizes)[-num_nuclei:]
        binary_mask = np.isin(labeled_mask, largest_labels)
    
    # Map original coordinates to the binary mask grid to filter localizations
    x_idx = np.clip(np.digitize(x, x_bins) - 1, 0, binary_mask.shape[0] - 1)
    y_idx = np.clip(np.digitize(y, y_bins) - 1, 0, binary_mask.shape[1] - 1)
    in_mask = binary_mask[x_idx, y_idx]
    
    me3_masked = me3_df[in_mask[:len(me3_df)]]
    ac_masked = ac_df[in_mask[len(me3_df):]]

    if show_plots:
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        axes = axes.flatten()

        extent = [x_bins.min(), x_bins.max(), y_bins.min(), y_bins.max()]
        
        vmax = np.percentile(density_map, 99)
        im = axes[0].imshow(density_map.T, origin="lower", extent=extent, vmax=vmax)
        
        divider = make_axes_locatable(axes[0])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        axes[0].set_title("Density Map")
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        axes[0].set_aspect('equal', adjustable='box')

        axes[1].imshow(binary_mask.T, cmap="gray", origin="lower", extent=extent)
        axes[1].set_title("Binary Mask")
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        axes[1].set_aspect('equal', adjustable='box')

        axes[2].scatter(x, y, s=0.05, alpha=0.3, color='blue')
        axes[2].set_title("Raw Localizations")
        axes[2].set_xlim([x_bins.min(), x_bins.max()])
        axes[2].set_ylim([y_bins.min(), y_bins.max()])
        axes[2].set_xticks([])
        axes[2].set_yticks([])
        axes[2].set_aspect('equal', adjustable='box')

        axes[3].scatter(x[in_mask], y[in_mask], s=0.05, alpha=0.3, color='blue')
        axes[3].set_title("Masked Localizations")
        axes[3].set_xlim([x_bins.min(), x_bins.max()])
        axes[3].set_ylim([y_bins.min(), y_bins.max()])
        axes[3].set_xticks([])
        axes[3].set_yticks([])
        axes[3].set_aspect('equal', adjustable='box')

        plt.tight_layout()
        plt.show()

    return binary_mask, me3_masked, ac_masked


def create_radial_contours(binary_mask, num_bands=100, show_plots=False):
    """
    Splits an arbitrarily shaped nucleus mask into concentric radial bands
    using the Euclidean distance transform, so localizations can be binned by
    depth into the nucleus (periphery vs. interior) regardless of its shape.

    Distance to the nearest background pixel is binned into num_bands equal
    steps from 0 to its maximum. Bands are numbered 0 = innermost core up to
    num_bands - 1 = outermost periphery; pixels outside the mask are -1.
    Because bands are equal distance steps, a mask whose deepest pixel is
    fewer than num_bands pixels from the edge leaves some bands empty (the edge
    pixels, at distance 1, then land a few bands short of num_bands - 1).

    Args:
        binary_mask: 2D bool mask indexed [x_bin, y_bin], from binarize_nucleus
        num_bands:   number of radial bands
        show_plots:  if True, show the distance map and the band map

    Returns:
        tuple: (distance_map, contour_bands)
            distance_map:  2D float array, distance (pixels) to the mask edge
            contour_bands: 2D int array, band index per pixel (-1 outside)
    """
    # Calculate the distance of each pixel inside the mask to the nearest background pixel
    distance_map = distance_transform_edt(binary_mask)
    
    # Create discrete contour bands
    contour_bands = np.zeros_like(distance_map, dtype=int)
    max_dist = distance_map.max()
    
    if max_dist > 0:
        # Bin the distances into the requested number of bands, then flip so
        # band 0 is the innermost core and band num_bands - 1 the outer periphery
        bins = np.linspace(0, max_dist, num_bands + 1)
        
        # Only assign bands to pixels inside the mask
        inside_mask = binary_mask > 0
        outside_mask = binary_mask <= 0
        
        contour_bands[inside_mask] = np.clip(np.digitize(distance_map[inside_mask], bins), 1, num_bands)
        contour_bands[inside_mask] = num_bands - contour_bands[inside_mask]
        contour_bands[outside_mask] = contour_bands[outside_mask] - 1
        
    if show_plots:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        
        im1 = axes[0].imshow(distance_map.T, origin='lower', cmap='viridis')
        axes[0].set_title('Distance Transform Map')
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
        
        cmap = plt.get_cmap('tab20_r').copy()
        cmap.set_under('white')
        im2 = axes[1].imshow(contour_bands.T, origin='lower', cmap=cmap, vmax=num_bands, vmin=0)
        axes[1].set_title(f'Radial Contours ({num_bands} Bands)')
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
        
        plt.tight_layout()
        plt.show()
        
    return distance_map, contour_bands


def assign_contour_bands(df, contour_bands, bin_size=50, x_min=None, y_min=None):
    """
    Looks up the radial contour band (from create_radial_contours) of every
    localization and adds it as a "contour_band" column (-1 = outside nucleus).

    Pass x_min / y_min = mask_origin(me3_df, ac_df) of the pre-mask inputs.
    If omitted they default to df's own minimum x / y, which is only the mask
    origin when df is the full pre-mask data; for masked or single-channel
    data it shifts every point into the wrong band.

    Args:
        df:            DataFrame with "x [nm]", "y [nm]"
        contour_bands: 2D int array from create_radial_contours
        bin_size:      pixel size in nm (must match contour_bands)
        x_min, y_min:  mask grid origin in nm; None = df's own minimum

    Returns:
        pd.DataFrame: a copy of df with a "contour_band" column
    """
    df = df.copy()
    coords = df[["x [nm]", "y [nm]"]].to_numpy()

    if x_min is None:
        x_min = coords[:, 0].min()
    if y_min is None:
        y_min = coords[:, 1].min()

    x_idx = np.clip(((coords[:, 0] - x_min) // bin_size).astype(int), 0, contour_bands.shape[0] - 1)
    y_idx = np.clip(((coords[:, 1] - y_min) // bin_size).astype(int), 0, contour_bands.shape[1] - 1)

    df["contour_band"] = contour_bands[x_idx, y_idx]
    return df


# --- Density subsampling -------------------------------------------------------

def sample_lower_densities(me3_df, ac_df, num_samples=10, random_state=42, show_plots=False, num_plots=4):
    """
    Randomly thins both channels to a ladder of lower labelling densities, to
    mimic reduced staining efficiency (independent thinning).

    Fractions are num_samples evenly spaced values from 1.0 down to
    1 / num_samples; fraction 1.0 returns unsampled copies.

    Args:
        me3_df, ac_df: per-mark localization DataFrames
        num_samples:   number of density fractions
        random_state:  seed for DataFrame.sample (ac uses random_state + 1)
        show_plots:    if True, scatter-plot a few of the subsampled fields
        num_plots:     how many fractions to show when show_plots is True

    Returns:
        tuple: (sampled_data, fractions)
            sampled_data: list of (me3_sampled, ac_sampled, fraction) tuples,
                          densest first
            fractions:    1D array of the fractions used
    """
    fractions = np.linspace(1.0, 1.0/num_samples, num_samples)
    
    sampled_data = []
    for frac in fractions:
        me3_samp = me3_df.sample(frac=frac, random_state=random_state) if frac < 1.0 else me3_df.copy()
        ac_samp = ac_df.sample(frac=frac, random_state=random_state+1) if frac < 1.0 else ac_df.copy() 
        sampled_data.append((me3_samp, ac_samp, frac))
        
    if show_plots:
        num_plots = min(num_plots, len(sampled_data))
        # Grab evenly spaced indices
        indices_to_plot = np.linspace(0, len(sampled_data) - 1, num_plots, dtype=int)
        
        cols = min(4, num_plots)
        rows = int(np.ceil(num_plots / cols))
        
        fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
        if rows * cols == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
            
        for i, ax in enumerate(axes):
            if i < num_plots:
                idx = indices_to_plot[i]
                me3_samp, ac_samp, frac = sampled_data[idx]
                
                ax.scatter(me3_samp["x [nm]"], me3_samp["y [nm]"], s=0.05, alpha=0.3, color='blue')
                ax.scatter(ac_samp["x [nm]"], ac_samp["y [nm]"], s=0.05, alpha=0.3, color='blue')
                
                ax.set_title(f"Density Fraction: {frac:.2f}")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_aspect('equal', adjustable='box')
            else:
                ax.axis('off') # Hide empty subplots
            
        plt.tight_layout()
        plt.show()

    return sampled_data, fractions
