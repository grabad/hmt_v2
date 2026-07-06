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


def plot_rdf_adf(me3_rdf, me3_adf, ac_rdf, ac_adf, step=10):
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

    _, axes = plt.subplots(1, 2, figsize=(9, 4), label=' ')

    # RDF
    axes[0].bar(r_centers - bar_width / 2, me3_rdf, width=bar_width, color='#77DD76', label='H3K27me3')
    axes[0].bar(r_centers + bar_width / 2, ac_rdf,  width=bar_width, color='#FF6962', label='H3K27ac')
    axes[0].set_xlabel('Radial distance (nm)')
    axes[0].set_ylabel('Neighbor density (locs / nm²)')
    axes[0].set_title('Radial Distribution Function')
    axes[0].legend()

    # ADF — normalize so both channels are directly comparable
    me3_adf_norm = me3_adf / me3_adf.sum()
    ac_adf_norm  = ac_adf  / ac_adf.sum()

    axes[1].bar(z_centers - bar_width / 2, me3_adf_norm, width=bar_width, color='#77DD76', label='H3K27me3')
    axes[1].bar(z_centers + bar_width / 2, ac_adf_norm,  width=bar_width, color='#FF6962', label='H3K27ac')
    axes[1].set_xlabel('|Z offset| (nm)')
    axes[1].set_ylabel('Probability')
    axes[1].set_title('Axial Distribution Function')
    axes[1].legend()

    plt.tight_layout()
    plt.show()