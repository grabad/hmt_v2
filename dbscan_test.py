import numpy as np
import pandas as pd
import locan as lc
import matplotlib.pyplot as plt
import napari
import inspect

if __name__=="__main__":
    lc.configuration.N_JOBS = -1

    me3_locdata = lc.locan_io.load_thunderstorm_file("test_data/k27_k27_thaw009_me3.csv")
    ac_locdata = lc.locan_io.load_thunderstorm_file("test_data/k27_k27_thaw009_ac.csv")

    me3_noise, me3_clusters = lc.cluster_dbscan(me3_locdata, eps=100, min_samples=3)
    ac_noise, ac_clusters = lc.cluster_dbscan(ac_locdata, eps=100, min_samples=3)
    
    print(me3_locdata)