"""
SR-EV: raw config-N.txt -> antibody/fluorophore labeling -> photophysics ->
rendered TIFF stack (single channel, srev.pipeline.run_pipeline) or a
dual-dye spatial + spectral sSMLM stack pair for RainbowSTORM
(srev.pipeline.run_multichannel_pipeline, srev.spectral), plus supporting
compartment/domain analysis. See srev.pipeline for the orchestrators.
"""
from . import io
from . import labeling
from . import photophysics
from . import compartments
from . import domains
from . import render
from . import spectral
from . import ground_truth
from . import pipeline
from . import visualize
