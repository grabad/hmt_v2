import napari
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import napari
from sklearn.cluster import DBSCAN
from scipy.ndimage import gaussian_filter, binary_fill_holes, label, distance_transform_edt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.qhull import QhullError
from scipy.optimize import minimize_scalar

def extract_empirical_parameters(real_df, raw_me3_df, raw_ac_df, contour_bands, bin_size=50, sdis=200, step=10):
    """
    Extracts the 1D Radial Distribution Function (RDF) from experimental data.
    This captures the spatial architecture (tight vs. sprawling) of the domains.
    """
    comb_df = pd.concat([raw_me3_df, raw_ac_df])
    x_min = comb_df["x [nm]"].min()
    y_min = comb_df["y [nm]"].min()
    
    coords_2d = real_df[["x [nm]", "y [nm]"]].to_numpy()
    
    # --- 1. Extract Radial Density Profile (for Seed Placement) ---
    x_idx = np.clip(((coords_2d[:, 0] - x_min) // bin_size).astype(int), 0, contour_bands.shape[0] - 1)
    y_idx = np.clip(((coords_2d[:, 1] - y_min) // bin_size).astype(int), 0, contour_bands.shape[1] - 1)
    loc_bands = contour_bands[x_idx, y_idx]
    
    num_bands = contour_bands.max()
    radial_density_profile = np.zeros(num_bands)
    for i in range(1, num_bands + 1):
        band_area = np.sum(contour_bands == i)
        locs_in_band = np.sum(loc_bands == i)
        if band_area > 0:
            radial_density_profile[i-1] = locs_in_band / band_area
            
    # --- 2. Extract Radial Distribution Function (RDF) ---
    tree = cKDTree(coords_2d)
    radii = np.arange(step, sdis + step, step)
    
    # Calculate average density of the whole nucleus to normalize the RDF
    # This prevents 'scaling' issues between different cells
    nuc_area_nm2 = np.sum(contour_bands > 0) * (bin_size**2)
    avg_density = len(coords_2d) / nuc_area_nm2
    
    rdf_profile = []
    prev_counts = np.zeros(len(coords_2d))
    
    for r in radii:
        current_counts = tree.query_ball_point(coords_2d, r, return_length=True)
        ring_counts = current_counts - prev_counts
        prev_counts = current_counts
        
        # Mean neighbors in this ring across all points
        mean_neighbors = np.mean(ring_counts)
        
        # Normalize by the area of the ring to get local density
        ring_area = np.pi * (r**2 - (r-step)**2)
        local_density = mean_neighbors / ring_area
        
        # RDF = Local Density / Global Average Density
        # Values > 1 indicate clustering/enrichment at this radius
        rdf_profile.append(local_density / avg_density)
            
    z_degradation = np.full(len(radii), 1.5)
            
    return radial_density_profile, np.array(rdf_profile), z_degradation, x_min, y_min