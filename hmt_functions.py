import napari
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import napari
from sklearn.cluster import DBSCAN

def filter_quantile(df, percent=0.05, column="z [nm]"):
    lower_bound = df[column].quantile(percent)
    upper_bound = df[column].quantile(1 - percent)
    
    filtered_df = df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]
    
    return filtered_df


def cluster_dbscan(loc_df, eps=50, min_samples=3):
    loc_np = loc_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy()
    loc_clust = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(loc_np)    
    loc_df["cluster_label"] = loc_clust
    
    
    rng = np.random.default_rng(0)  # reproducible colors
    alpha = 0.5

    loc_unique_clusters = np.unique(loc_clust[loc_clust != -1])
    loc_cluster_colors = {c: np.append(rng.random(3), 1) for c in loc_unique_clusters}
    loc_colors = []

    for c in loc_clust:
        if c == -1:
            loc_colors.append((0.7, 0.7, 0.7, alpha))
        else:
            loc_colors.append(loc_cluster_colors[c])
            
    loc_df["cluster_color"] = loc_colors

def plot_clusters_napari(me3_df, ac_df, scaling_factor=10):
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
    viewer.add_points(
        data=ac_noise[coords],
        features=ac_noise,
        face_color=ac_noise["cluster_color"].tolist(),
        border_color=ac_border, 
        border_width=0.01, 
        size=ac_noise["plot_size"],    
        name="H3K27ac Noise"
    )           
    viewer.add_points(
        data=me3_noise[coords],
        features=me3_noise,
        face_color=me3_noise["cluster_color"].tolist(), 
        border_color=me3_border, 
        border_width=0.01, 
        size=me3_noise["plot_size"],   
        name="H3K27me3 Noise"
    )                     
    viewer.add_points(
        data=ac_nano[coords],
        features=ac_nano,
        face_color=ac_nano["cluster_color"].tolist(), 
        border_color=ac_border, 
        border_width=0.01, 
        size=ac_nano["plot_size"],    
        name="H3K27ac Nanodomains"
    )           
    viewer.add_points(
        data=me3_nano[coords],
        features=me3_nano,
        face_color=me3_nano["cluster_color"].tolist(),
        border_color=me3_border, 
        border_width=0.01, 
        size=me3_nano["plot_size"],   
        name="H3K27me3 Nanodomains"
)