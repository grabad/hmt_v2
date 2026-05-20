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