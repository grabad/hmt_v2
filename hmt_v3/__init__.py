"""
hmt_v3: quantify the spatial distribution of histone post-translational
modifications (e.g. H3K27me3 vs H3K27ac) from two-colour 3D sSMLM localizations.

Modules, in the order the analysis uses them (see README.md for the full guide):

    preprocess   filter localizations, mask the nucleus, radial contour bands
    measure      density-independent statistics of real data (radial profile,
                 RDF/ADF, pair correlation g(r) / Ripley's K, noise fraction)
    simulate     synthetic nuclei with known nanodomains, built from `measure`
    calibrate    fit DBSCAN eps / min_samples vs. density on simulated ground truth
    cluster      density-corrected DBSCAN and per-nanodomain size/count
    postprocess  two-mark spatial metrics (SNCR, cross g(r), Ripley's cross-K,
                 Delaunay / MST / mesh mixing)
    sweep        toy two-mark centroid model for validating the postprocess metrics
    export       write simulated localizations as ThunderSTORM CSVs
    visualize    matplotlib / napari plots
    srev         SR-EV chromatin-polymer -> labelled, rendered sSMLM movie simulator

Every DataFrame of localizations uses the RainbowSTORM / ThunderSTORM column
names "x [nm]", "y [nm]", "z [nm]".
"""
from . import preprocess
from . import measure
from . import simulate
from . import calibrate
from . import cluster
from . import postprocess
from . import sweep
from . import export
from . import visualize
from . import srev
