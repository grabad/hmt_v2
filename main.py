import numpy as np
import pandas as pd
import locan as lc
import matplotlib.pyplot as plt
import napari

if __name__=="__main__":
    me3_locdata = lc.locan_io.load_thunderstorm_file("test_data/k27_k27_thaw009_me3.csv")
    ac_locdata = lc.locan_io.load_thunderstorm_file("test_data/k27_k27_thaw009_ac.csv")

    me3_size = me3_locdata.data['sigmax [nm]'].to_numpy()/110
    ac_size = ac_locdata.data['sigmax [nm]'].to_numpy()/110

    viewer = napari.Viewer()
    viewer.add_points(me3_locdata.coordinates, name="H3K27me3", face_color="green", border_color="green", border_width=0, size=me3_size, blending="additive")
    viewer.add_points(ac_locdata.coordinates, name="H3K27ac", face_color="red", border_color="red", border_width=0, size=ac_size, blending="additive")
    viewer.show()
    napari.run()
