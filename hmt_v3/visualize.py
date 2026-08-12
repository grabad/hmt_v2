import napari
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import napari
from sklearn.cluster import DBSCAN
from scipy.ndimage import gaussian_filter, binary_fill_holes, label, distance_transform_edt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.qhull import QhullError
from scipy.optimize import minimize_scalar



def plot_clusters_napari(me3_df, ac_df, scaling_factor=10):
    me3_df = me3_df.copy()
    ac_df = ac_df.copy()
    
    me3_df["plot_size"] = me3_df['sigmax [nm]'].to_numpy()/110*scaling_factor
    ac_df["plot_size"] = ac_df['sigmax [nm]'].to_numpy()/110*scaling_factor

    coords = ["x [nm]", "y [nm]", "z [nm]"]
    me3_border = [0, 0.7, 0, 1]
    ac_border = [0.7, 0, 0, 1]

    me3_nano = me3_df[me3_df["cluster_label"] != -1]
    me3_noise = me3_df[me3_df["cluster_label"] == -1]
    ac_nano = ac_df[ac_df["cluster_label"] != -1]
    ac_noise = ac_df[ac_df["cluster_label"] == -1]

    viewer = napari.Viewer()
    
    if not ac_noise.empty:
        viewer.add_points(
            data=ac_noise[coords],
            features=ac_noise,
            face_color=ac_noise["cluster_color"].tolist(),
            border_color=ac_border, 
            border_width=0.01, 
            size=ac_noise["plot_size"],    
            name="H3K27ac Noise",
            out_of_slice_display=True
        )
    
    if not me3_noise.empty:   
        viewer.add_points(
            data=me3_noise[coords],
            features=me3_noise,
            face_color=me3_noise["cluster_color"].tolist(), 
            border_color=me3_border, 
            border_width=0.01, 
            size=me3_noise["plot_size"],   
            name="H3K27me3 Noise",
            out_of_slice_display=True
        )
    
    if not ac_nano.empty:                     
        viewer.add_points(
            data=ac_nano[coords],
            features=ac_nano,
            face_color=ac_nano["cluster_color"].tolist(), 
            border_color=ac_border, 
            border_width=0.01, 
            size=ac_nano["plot_size"],    
            name="H3K27ac Nanodomains",
            out_of_slice_display=True
        )
        
    if not me3_nano.empty:           
        viewer.add_points(
            data=me3_nano[coords],
            features=me3_nano,
            face_color=me3_nano["cluster_color"].tolist(),
            border_color=me3_border, 
            border_width=0.01, 
            size=me3_nano["plot_size"],   
            name="H3K27me3 Nanodomains",
            out_of_slice_display=True
        )

    viewer.dims.order = (2, 1, 0)
    viewer.dims.ndisplay = 3

def plot_seeds_napari(me3_seeds=None, ac_seeds=None, seed_size=15):
    """
    Visualizes raw (N, 3) numpy matrices of ground-truth seeds in 3D.
    
    Args:
        me3_seeds (np.ndarray, optional): (N, 3) matrix of H3K27me3 seeds. Plotted in green.
        ac_seeds (np.ndarray, optional): (N, 3) matrix of H3K27ac seeds. Plotted in red.
        seed_size (int): Visual size of the points.
    """
    viewer = napari.Viewer()
    
    # Exact colors from your clustering script
    me3_color = [0, 0.7, 0, 1]  # Green
    ac_color = [0.7, 0, 0, 1]   # Red
    
    # Plot H3K27me3 Seeds
    if me3_seeds is not None and len(me3_seeds) > 0:
        viewer.add_points(
            data=me3_seeds,
            face_color=me3_color,
            border_width=0, 
            size=seed_size,    
            name="H3K27me3 Seeds",
            out_of_slice_display=True
        )
        
    # Plot H3K27ac Seeds
    if ac_seeds is not None and len(ac_seeds) > 0:
        viewer.add_points(
            data=ac_seeds,
            face_color=ac_color,
            border_width=0, 
            size=seed_size,    
            name="H3K27ac Seeds",
            out_of_slice_display=True
        )

    # Maintain the exact spatial orientation from your cluster viewer
    viewer.dims.order = (2, 1, 0)
    viewer.dims.ndisplay = 3
    
def plot_nanodomain_2d(df, seeds=None, title="Nanodomains", point_size=10, alpha=0.5, plot_seeds=True, ax=None):
    """
    2D scatter plot (XY projection) of localizations from spawn_nanodomains.
    Each cluster_label is drawn in a distinct color; cluster_label == -1 is drawn in gray.

    Args:
        df: DataFrame with columns [x [nm], y [nm], z [nm], cluster_label]
        seeds: optional (N, 3) array of seed coordinates — plotted as red stars
        title: plot title
        point_size: marker size for localizations
        alpha: opacity of localization points
        ax: optional existing Axes to draw into; if None a new figure is created
    """
    _standalone = ax is None
    if _standalone:
        fig = plt.figure(figsize=(7, 6), label=' ')
        ax = fig.add_subplot(111)

    noise = df[df['cluster_label'] == -1]
    if not noise.empty:
        ax.scatter(noise['x [nm]'], noise['y [nm]'],
                   s=point_size, alpha=alpha * 0.6, color='gray', label='Noise')

    cmap = plt.get_cmap('tab10')
    for cluster_id, group in df[df['cluster_label'] != -1].groupby('cluster_label'):
        color = cmap(int(cluster_id) % 10)
        ax.scatter(group['x [nm]'], group['y [nm]'],
                   s=point_size, alpha=alpha, color=color)

    if seeds is not None and plot_seeds:
        seeds = np.asarray(seeds)
        ax.scatter(seeds[:, 0], seeds[:, 1],
                   s=80, color='red', marker='*', zorder=5, label='Seeds')
        ax.legend(loc='upper left', fontsize=8)

    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_title(title)
    ax.set_aspect('equal')

    if _standalone:
        plt.tight_layout()
        plt.show()

def plot_nanodomain_3d(df, seeds=None, title="Nanodomains", point_size=10, alpha=0.5, plot_seeds=False, plot_hulls=True, hull_alpha=0.15):
    """
    3D scatter plot of localizations from spawn_nanodomains.
    Each cluster_label is drawn in a distinct color.

    Args:
        df: DataFrame with columns [x [nm], y [nm], z [nm], cluster_label]
        seeds: optional (N, 3) array of seed coordinates — plotted as red stars
        title: plot title
        point_size: marker size for localizations
        alpha: opacity of localization points
        plot_hulls: if True, draw the convex hull of each cluster
        hull_alpha: opacity of convex hull faces
    """
    fig = plt.figure(figsize=(7, 6), label=' ')
    ax = fig.add_subplot(111, projection='3d')

    noise = df[df['cluster_label'] == -1]
    if not noise.empty:
        pts = noise[['x [nm]', 'y [nm]', 'z [nm]']].to_numpy()
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
                   s=point_size, alpha=alpha * 0.6, color='gray', label='Noise')

    cmap = plt.get_cmap('tab10')
    for cluster_id, group in df[df['cluster_label'] != -1].groupby('cluster_label'):
        color = cmap(int(cluster_id) % 10)
        pts = group[['x [nm]', 'y [nm]', 'z [nm]']].to_numpy()
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
                   s=point_size, alpha=alpha, color=color)

        if plot_hulls and len(pts) >= 4:
            try:
                hull = ConvexHull(pts)
                faces = [pts[simplex] for simplex in hull.simplices]
                poly = Poly3DCollection(faces, alpha=hull_alpha,
                                        facecolor=color, edgecolor=color, linewidth=0.7)
                ax.add_collection3d(poly)
            except QhullError:
                pass

    if seeds is not None and plot_seeds:
        seeds = np.asarray(seeds)
        ax.scatter(seeds[:, 0], seeds[:, 1], seeds[:, 2],
                   s=80, color='red', marker='*', zorder=5, label='Seeds')
        ax.legend(loc='upper left', fontsize=8)

    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_zlabel('Z (nm)')
    ax.set_title(title)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.show()


def plot_rdf_adf(me3_rdf, me3_adf, ac_rdf, ac_adf, step=10, rdf_title=None, adf_title=None):
    """
    Plots the RDF and ADF for H3K27me3 and H3K27ac on shared axes for direct comparison.

    Args:
        me3_rdf: 1D RDF array for H3K27me3
        me3_adf: 1D ADF array for H3K27me3
        ac_rdf:  1D RDF array for H3K27ac
        ac_adf:  1D ADF array for H3K27ac
        step:    bin width in nm (must match extraction)
    """
    me3_rdf = np.asarray(me3_rdf, dtype=float)
    me3_adf = np.asarray(me3_adf, dtype=float)
    ac_rdf  = np.asarray(ac_rdf,  dtype=float)
    ac_adf  = np.asarray(ac_adf,  dtype=float)

    r_centers   = (np.arange(len(me3_rdf)) + 0.5) * step
    z_centers   = (np.arange(len(me3_adf)) + 0.5) * step
    bar_width   = step * 0.45

    _, axes = plt.subplots(1, 2, figsize=(10, 4), label=' ')

    # RDF
    axes[0].bar(r_centers - bar_width / 2, me3_rdf*1000, width=bar_width, color='#77DD76', label='H3K27me3')
    axes[0].bar(r_centers + bar_width / 2, ac_rdf*1000,  width=bar_width, color='#FF6962', label='H3K27ac')
    axes[0].set_xlabel('Radial distance (nm)')
    axes[0].set_ylabel('Local Density (millilocs / nm²)')
    if rdf_title is not None:
        axes[0].set_title(rdf_title)
    else:
        axes[0].set_title('Radial Distribution Function')
    axes[0].legend()

    # ADF — normalize so both channels are directly comparable
    me3_adf_norm = me3_adf / me3_adf.sum()
    ac_adf_norm  = ac_adf  / ac_adf.sum()

    axes[1].bar(z_centers - bar_width / 2, me3_adf_norm, width=bar_width, color='#77DD76', label='H3K27me3')
    axes[1].bar(z_centers + bar_width / 2, ac_adf_norm,  width=bar_width, color='#FF6962', label='H3K27ac')
    axes[1].set_xlabel('|Z offset| (nm)')
    axes[1].set_ylabel('Probability')
    if adf_title is not None:
        axes[1].set_title(adf_title)
    else:
        axes[1].set_title('Axial Distribution Function')
    axes[1].legend()

    plt.tight_layout()
    plt.show()


def plot_seed_spacing(empirical_nnd, sim_seeds, bins=30, label_real='Real domains',
                      label_sim='Simulated seeds', ax=None):
    """
    Compares the inter-domain nearest-neighbour distance (NND) distribution of the
    real nucleus against the simulated seeds. Use this to verify that a
    spacing-matched placement (place_seeds_matched / generate_nucleus) reproduces
    the empirical domain arrangement rather than an independent-Poisson one.

    Args:
        empirical_nnd: 1D array of real NNDs (from simulate.measure_seed_spacing)
        sim_seeds:     (N, 3) or (N, 2) array of simulated seed coordinates
        bins:          number of histogram bins (shared across both curves)
        label_real:    legend label for the empirical distribution
        label_sim:     legend label for the simulated distribution
        ax:            optional existing Axes; a new figure is created if None
    """
    empirical_nnd = np.asarray(empirical_nnd, dtype=float)
    sim_seeds = np.asarray(sim_seeds, dtype=float)

    tree = cKDTree(sim_seeds[:, :2])
    d, _ = tree.query(sim_seeds[:, :2], k=2)
    sim_nnd = d[:, 1]

    _standalone = ax is None
    if _standalone:
        fig = plt.figure(figsize=(7, 5), label=' ')
        ax = fig.add_subplot(111)

    lo = 0.0
    hi = max(empirical_nnd.max(), sim_nnd.max())
    edges = np.linspace(lo, hi, bins + 1)

    ax.hist(empirical_nnd, bins=edges, density=True, alpha=0.5,
            color='#4C72B0', label=f'{label_real} (median {np.median(empirical_nnd):.0f} nm)')
    ax.hist(sim_nnd, bins=edges, density=True, alpha=0.5,
            color='#DD8452', label=f'{label_sim} (median {np.median(sim_nnd):.0f} nm)')

    ax.set_xlabel('Nearest-neighbour distance between domains (nm)')
    ax.set_ylabel('Probability density')
    ax.set_title('Inter-domain spacing: real vs. simulated')
    ax.legend()

    if _standalone:
        plt.tight_layout()
        plt.show()


def plot_domain_packing(result, ax=None):
    """
    Plots the domain overlap-margin distribution from
    simulate.diagnose_domain_packing: nearest-neighbour centroid distance minus
    the sum of the two domains' radii, for every true domain. Negative values
    mean that domain's nearest neighbour geometrically overlaps it — a limit no
    DBSCAN parameter choice can get past, since the two domains' localizations
    are genuinely interspersed in the underlying point cloud.

    Args:
        result: dict from simulate.diagnose_domain_packing
        ax:     optional existing Axes; a new figure is created if None
    """
    margin = np.asarray(result['overlap_margin'], dtype=float)
    if len(margin) == 0:
        print("Fewer than 2 domains — nothing to plot.")
        return

    _standalone = ax is None
    if _standalone:
        fig = plt.figure(figsize=(7, 5), label=' ')
        ax = fig.add_subplot(111)

    ax.hist(margin, bins=40, color='#4C72B0', alpha=0.8)
    ax.axvline(0, color='black', linestyle='--', linewidth=1.2,
              label=f'{result["frac_overlapping"]:.0%} of domains overlap their nearest neighbour')

    ax.set_xlabel('Nearest-neighbour overlap margin (nm)\n[NN centroid distance − (radius₁ + radius₂)]')
    ax.set_ylabel('Number of domains')
    ax.set_title('Domain packing: negative = geometrically overlapping')
    ax.legend()

    if _standalone:
        plt.tight_layout()
        plt.show()


def plot_pair_correlation(pcf_real, pcf_sim=None, r_domain=None,
                          label_real='Real', label_sim='Simulated', color=None):
    """
    Plots Ripley's L(r) - r and the pair correlation g(r) for the real cell and,
    optionally, the simulated reconstruction, on shared axes. Overlap means the
    simulation reproduces the density-independent clustering structure measured by
    simulate.measure_pair_correlation.

    Reference lines: L(r) - r = 0 and g(r) = 1 are complete spatial randomness.
    L(r) - r > 0 (g > 1) is clustering; L(r) - r < 0 (g < 1) is regularity.

    Args:
        pcf_real:   dict from simulate.measure_pair_correlation (real cell)
        pcf_sim:    optional dict from measure_pair_correlation (simulation)
        r_domain:   optional domain-scale radius (nm) to mark with a vertical line
        label_real: legend label for the real curve
        label_sim:  legend label for the simulated curve
    """
    _, axes = plt.subplots(1, 2, figsize=(11, 4.5), label=' ')

    if color == 'green':
        line_color = '#77DD76'

    elif color == 'red':
        line_color = '#FF6962'
    else:
        line_color = 'black'

    lw = 2
    axes[0].plot(pcf_real['r'], pcf_real['L'] - pcf_real['r'], color=line_color, label=label_real, lw=lw)
    if pcf_sim is not None:
        axes[0].plot(pcf_sim['r'], pcf_sim['L'] - pcf_sim['r'], color=line_color, label=label_sim, ls='--', lw=lw)
    axes[0].set_xlabel('r (nm)')
    axes[0].set_ylabel('L(r) - r (nm)')
    axes[0].set_title("Clustering Similarity via Ripley's L")
    axes[0].legend(loc='upper right')

    axes[1].plot(pcf_real['r'], pcf_real['g'], color=line_color, label=label_real, lw=lw)
    if pcf_sim is not None:
        axes[1].plot(pcf_sim['r'], pcf_sim['g'], color=line_color, label=label_sim, ls='--', lw=lw)
    axes[1].set_xlabel('r (nm)')
    axes[1].set_ylabel('g(r)')
    axes[1].set_title('Clustering Similarity via Pairwise Correlation')
    axes[1].legend(loc='upper right')

    # if r_domain is not None:
    #     for ax in axes:
    #         ax.axvline(r_domain, color='green', lw=0.8, ls=':')

    plt.tight_layout()
    plt.show()


def plot_dbscan_calibration(me3_result, ac_result, title='DBSCAN Calibration: H3K27me3 vs H3K27ac'):
    """
    Plots the H3K27me3 and H3K27ac DBSCAN calibration curves on separate axes
    (one row per mark): eps fit as a power law of density, min_samples fit as
    a linear function of density. Reads the 'densities' / 'best_epsilons' /
    'eps_coef' / 'eps_exp' / 'eps_r_squared' / 'best_min_samples' /
    'min_samples_slope' / 'min_samples_intercept' / 'min_samples_r_squared'
    keys returned by simulate.optimize_dbscan_size_matching.

    Args:
        me3_result: dict returned by optimize_dbscan_size_matching for H3K27me3
        ac_result:  dict returned by optimize_dbscan_size_matching for H3K27ac
        title:      figure suptitle
    """
    me3_color = '#77DD76'
    ac_color = '#FF6962'

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), label=' ')

    for row, (result, color, label) in enumerate(
        [(me3_result, me3_color, 'H3K27me3'), (ac_result, ac_color, 'H3K27ac')]
    ):
        ax_eps, ax_ms = axes[row, 0], axes[row, 1]

        dens = np.asarray(result['densities'], dtype=float)
        if len(dens) == 0:
            continue
        dens_smooth = np.linspace(dens.min(), dens.max(), 200)

        eps = np.asarray(result['best_epsilons'], dtype=float)
        eps_coef, eps_exp = result['eps_coef'], result['eps_exp']
        ax_eps.scatter(dens, eps, color=color, label=f'R²={result["eps_r_squared"]:.3f}')
        ax_eps.plot(dens_smooth, eps_coef * dens_smooth ** eps_exp, color=color, linestyle='--')
        ax_eps.set_xlabel('Density [locs/nm²]')
        ax_eps.set_ylabel('Epsilon [nm]')
        ax_eps.set_title(f'{label} — Epsilon Calibration')

        min_samples = np.asarray(result['best_min_samples'], dtype=float)
        ms_slope, ms_intercept = result['min_samples_slope'], result['min_samples_intercept']
        ax_ms.scatter(dens, min_samples, color=color, label=f'R²={result["min_samples_r_squared"]:.3f}')
        ax_ms.plot(dens_smooth, ms_slope * dens_smooth + ms_intercept, color=color, linestyle='--')
        ax_ms.set_xlabel('Density [locs/nm²]')
        ax_ms.set_ylabel('min_samples')
        ax_ms.set_title(f'{label} — min_samples Calibration')

    for ax in axes.ravel():
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.suptitle(title)
    plt.tight_layout()
    plt.show()


def plot_size_vs_density(density_fractions, arbitrary_sizes, corrected_sizes, true_diameter_nm,
                         arbitrary_color='#6B8FD6', corrected_color='#F5A623',
                         title='Measured Nanodomain Size vs. Localization Density', text_loc=(0.2, 0.15)):
    """
    Line plot of a single simulated nanodomain's measured size (geometric_2D
    diameter) across a labelling-density gradient, comparing DBSCAN clustered
    with a fixed ("arbitrary") eps/min_samples pair against DBSCAN clustered with
    density-corrected eps/min_samples (cluster.cluster_dbscan_density_corrected).

    A horizontal dashed line marks the domain's true diameter — known exactly,
    since simulate.simulate_single_domain_scene generates a single isolated
    domain rather than a dense field of possibly-overlapping ones. The closer a
    curve tracks that line across density, the less its size measurement is
    biased by how densely the sample happened to be labelled; an "arbitrary"
    curve that drifts away from the true line as density drops (while the
    corrected curve stays flat) is the artifact density correction exists to fix.

    Color here encodes the arbitrary/corrected condition directly (rather than a
    mark identity, as plot_dbscan_calibration's me3/ac colors do) since that
    condition is the entity being compared in this figure.

    Args:
        density_fractions: x-values, the labelling-density fraction (1.0 = full
                           density) each point was subsampled to
        arbitrary_sizes:    measured diameter (nm) at each fraction using fixed
                           eps/min_samples
        corrected_sizes:    measured diameter (nm) at each fraction using
                           density-corrected eps/min_samples
        true_diameter_nm:   the domain's known true diameter (nm)
        arbitrary_color:    line color for the arbitrary-params series
        corrected_color:    line color for the density-corrected-params series
        title:              figure title
    """
    fig, ax = plt.subplots(figsize=(6.5, 4.8))

    ax.axhline(true_diameter_nm, color='0.3', linestyle=':', linewidth=2,
               label=f'True diameter ({true_diameter_nm:.0f} nm)')

    ax.plot(density_fractions, arbitrary_sizes, 'o-', color=arbitrary_color,
            label='Constant parameters')
    ax.plot(density_fractions, corrected_sizes, 's-', color=corrected_color,
            label='Density-corrected parameters')

    ax.set_xlabel('Localization density fraction')
    ax.set_ylabel('Measured diameter [nm]')
    ax.set_title(title)
    ax.legend(loc='lower right', fontsize=12)
    ax.grid(True, alpha=0.3)

    if any(np.isnan(s) for s in arbitrary_sizes) or any(np.isnan(s) for s in corrected_sizes):
        text_loc_x = text_loc[0]
        text_loc_y = text_loc[1]
        ax.text(text_loc_x, text_loc_y, 'Failed to\nidentify nanodomain', transform=ax.transAxes,
                ha='center', va='bottom', color='red', fontsize=16)

    plt.tight_layout()
    plt.show()