"""
SR-EV: raw config-N.txt -> antibody/fluorophore labeling -> photophysics ->
rendered TIFF stack, plus supporting compartment/domain analysis. See
srev.pipeline for the top-level orchestrator.
"""
from . import io
from . import labeling
from . import photophysics
from . import compartments
from . import domains
from . import render
from . import ground_truth
from . import pipeline
