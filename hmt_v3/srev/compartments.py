"""
Per-monomer heterochromatin/euchromatin classification for an SR-EV
monomer cloud, via coordination number (CN).

This replaces an earlier version of this module that tried to run
k-means with a cluster count k derived from an SR-EV generator parameter
(first a hardcoded 10, matching the legacy MATLAB vis_SREV_int.m; then a
guess that the run's "L" folder-name parameter set the count). Both were
wrong. The actual method — confirmed against Cangnano et al. 2024 (eLife,
"Local Volume Concentration, Packing Domains and Scaling Properties of
Chromatin"), the paper describing this exact SR-EV model — needs no
cluster count at all for per-monomer classification: each nucleosome's
coordination number, the count of other nucleosomes within an 11.5 nm
radius, directly measures how deeply it sits inside a packed domain (CN
0 = isolated, CN 12 = fully packed domain interior — their Figure 3 uses
this to color monomers heterochromatin (high CN, red) vs euchromatin (low
CN, blue)). Verified empirically against hmt_v1/Simulation/config-1.txt
(with length_unit_to_nm=10, see io.parse_config_txt): CN over the
whole configuration spans exactly 0-12 with a broad distribution, matching
the paper's stated range.

The paper separately describes a fully different procedure for counting
and sizing "packing domains" as a population-level statistic (rendering a
Gaussian-bead density image, finding its local maxima, and measuring each
domain's radius from a radial density profile — their Appendix 2). That
is a domain-counting/validation tool, not a per-monomer labeling step, and
is not implemented here; see the module docstring note in
pipeline.py if it's needed later.
"""
import numpy as np
from scipy.spatial import cKDTree


def compute_coordination_number(monomers_df, coordination_radius_nm=11.5):
    """
    Counts, for every monomer, how many OTHER monomers lie within
    coordination_radius_nm — the CN measure from Cangnano et al. 2024.
    A monomer at the interior of a packing domain has many close
    neighbors (CN up to 12 in the source paper); an isolated monomer in
    an open region has few or none (CN down to 0).

    Args:
        monomers_df:            output of io.parse_config_txt
        coordination_radius_nm: neighbor-counting radius (11.5 nm in the
                                 source paper — do not change without a
                                 reason, since it's what makes CN values
                                 comparable to their reported 0-12 range)

    Returns:
        np.ndarray of shape (len(monomers_df),): integer coordination
        number per monomer, in monomers_df's row order
    """
    xyz = monomers_df[["x [nm]", "y [nm]", "z [nm]"]].to_numpy(dtype=float)
    tree = cKDTree(xyz)
    counts = tree.query_ball_point(xyz, coordination_radius_nm, return_length=True)
    return counts - 1  # exclude each point counting itself


def classify_compartments(monomers_df, coordination_radius_nm=11.5, heterochromatin_cn_threshold=6):
    """
    Labels each monomer "heterochromatin" (densely coordinated, CN >=
    threshold) or "euchromatin" (openly coordinated, CN < threshold),
    using compute_coordination_number.

    The source paper shows CN as a continuous red-blue colormap (Figure 3)
    rather than stating a single binary cutoff; it does describe CN~6 as
    characteristic of the periphery/surface of a packing domain (the
    transition zone), which is why 6 — the midpoint of the reported 0-12
    range — is used here as the default split. Revisit this threshold
    once real per-mark labeling density data is available to calibrate
    against, rather than treating 6 as an established biological cutoff.

    Args:
        monomers_df:                  output of io.parse_config_txt
        coordination_radius_nm:       passed through to
                                       compute_coordination_number
        heterochromatin_cn_threshold: monomers with CN >= this are
                                       labeled heterochromatin

    Returns:
        monomers_df, unmodified except for two added columns:
            coordination_number: int, from compute_coordination_number
            compartment: "heterochromatin" | "euchromatin"
    """
    monomers_df = monomers_df.copy()
    cn = compute_coordination_number(monomers_df, coordination_radius_nm)
    monomers_df["coordination_number"] = cn
    monomers_df["compartment"] = np.where(
        cn >= heterochromatin_cn_threshold, "heterochromatin", "euchromatin"
    )
    return monomers_df
