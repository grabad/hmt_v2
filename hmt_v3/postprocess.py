import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial import cKDTree, Delaunay
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree


# --- Toy centroid generation --------------------------------------------------

def simulate_toy_centroids(n_domains=150, integration_level=0.5, field_size_nm=10000.0,
                            pair_jitter_nm=60.0, n_domains_ac=None, rng=None):
    """
    Two synthetic sets of nanodomain centroids (me3, ac) inside a square field,
    with a single knob controlling how spatially "integrated" (co-localized)
    the two channels are.

    integration_level in [0, 1]:
        0 -> me3 and ac centroids are two fully independent CSR (uniform
             random) point processes -- no relationship between channels.
        1 -> every centroid in the smaller-count channel is paired to one in
             the other, offset by a small Gaussian jitter (pair_jitter_nm) --
             the channels sit on top of each other as much as their relative
             abundance allows.
        0 < x < 1 -> that fraction of the smaller channel's domains are
             paired; the rest are independent, unpaired CSR centroids in
             each channel.

    Args:
        n_domains:        number of me3 domains
        integration_level: fraction of domains that are spatially paired
        field_size_nm:     side length of the square field (nm)
        pair_jitter_nm:    stdev of the 2D offset between a paired me3/ac centroid
        n_domains_ac:      number of ac domains; None (default) makes it equal
                           to n_domains. Set this different from n_domains to
                           test an abundance imbalance between channels --
                           colocalization_fraction / cross_pair_correlation /
                           self_nonself_contact_ratio's *_norm curves are
                           designed to stay comparable across that imbalance,
                           while the raw self_nonself_contact_ratio curves are
                           NOT (see that function's docstring) -- this is the
                           case that makes the raw-vs-normalized distinction
                           actually visible, unlike the balanced default.
        rng:               np.random.default_rng() instance; None makes a fresh one

    Returns:
        me3_df, ac_df: DataFrames [x [nm], y [nm], pair_id]
            pair_id >= 0 links a paired me3/ac centroid across channels;
            pair_id == -1 marks an independent (unpaired) centroid.
    """
    rng = rng or np.random.default_rng()
    integration_level = float(np.clip(integration_level, 0.0, 1.0))
    n_domains_ac = n_domains if n_domains_ac is None else n_domains_ac

    n_paired = int(round(integration_level * min(n_domains, n_domains_ac)))
    n_me3_indep = n_domains - n_paired
    n_ac_indep = n_domains_ac - n_paired

    parent_xy = rng.uniform(0, field_size_nm, size=(n_paired, 2))
    me3_paired = parent_xy + rng.normal(0, pair_jitter_nm, size=(n_paired, 2))
    ac_paired = parent_xy + rng.normal(0, pair_jitter_nm, size=(n_paired, 2))
    pair_ids = np.arange(n_paired)

    me3_indep = rng.uniform(0, field_size_nm, size=(n_me3_indep, 2))
    ac_indep = rng.uniform(0, field_size_nm, size=(n_ac_indep, 2))

    me3_xy = np.vstack([me3_paired, me3_indep])
    ac_xy = np.vstack([ac_paired, ac_indep])
    me3_pair_id = np.concatenate([pair_ids, np.full(n_me3_indep, -1)])
    ac_pair_id = np.concatenate([pair_ids, np.full(n_ac_indep, -1)])

    me3_df = pd.DataFrame(np.column_stack([me3_xy, me3_pair_id]), columns=["x [nm]", "y [nm]", "pair_id"])
    ac_df = pd.DataFrame(np.column_stack([ac_xy, ac_pair_id]), columns=["x [nm]", "y [nm]", "pair_id"])
    return me3_df, ac_df


def generate_integration_sweep(integration_levels=(0.0, 0.25, 0.5, 0.75, 1.0), n_domains=150,
                               n_domains_ac=None, field_size_nm=10000.0, pair_jitter_nm=60.0, rng=None):
    """
    Generates me3/ac toy centroids at each integration_level in one call, so
    every downstream metric (colocalization_fraction, cross_pair_correlation,
    self_nonself_contact_ratio, or a custom one) is computed on and plotted
    against the SAME generated centroids instead of each caller re-simulating
    its own draw. The balanced case (equal me3/ac counts) and the imbalanced
    case (n_domains_ac != n_domains) both go through this one function --
    just pass n_domains_ac to get the imbalanced distribution.

    Args:
        integration_levels: iterable of integration_level values to sweep
        n_domains, n_domains_ac, field_size_nm, pair_jitter_nm: forwarded to
            simulate_toy_centroids at each level (see its docstring)
        rng: np.random.default_rng() instance, reused (not reset) across
             levels so the whole sweep is reproducible from one seed while
             each level still draws fresh points; None makes a fresh one

    Returns:
        list of dicts, one per level, each with keys:
            'integration_level', 'me3_seeds', 'ac_seeds'
    """
    rng = rng or np.random.default_rng()
    sweep = []
    for level in integration_levels:
        me3_seeds, ac_seeds = simulate_toy_centroids(
            n_domains=n_domains, n_domains_ac=n_domains_ac, integration_level=level,
            field_size_nm=field_size_nm, pair_jitter_nm=pair_jitter_nm, rng=rng)
        sweep.append({"integration_level": level, "me3_seeds": me3_seeds, "ac_seeds": ac_seeds})
    return sweep


# --- Centroid-based spatial-distribution metrics -------------------------------

def cross_nn_distances(centroids_a, centroids_b):
    """For each point in A, distance to its nearest neighbour in B (nm)."""
    tree_b = cKDTree(centroids_b)
    dist, _ = tree_b.query(centroids_a, k=1)
    return dist


def colocalization_fraction(centroids_a, centroids_b, radius_nm):
    """
    Fraction of A centroids that have at least one B centroid within radius_nm.
    A single-number summary of cross-channel integration at one length scale.
    """
    return float(np.mean(cross_nn_distances(centroids_a, centroids_b) <= radius_nm))


def cross_pair_correlation(centroids_a, centroids_b, field_size_nm, r_max=800.0, dr=25.0):
    """
    Cross pair-correlation function g_ab(r): density of B centroids in an
    annulus at distance r from an A centroid, relative to the density expected
    under complete spatial randomness (g=1 -> no association, g>1 -> attraction/
    co-clustering, g<1 -> exclusion/segregation). Quantifies integration as a
    function of length scale instead of a single threshold.

    Edge correction: only A centroids at least r_max from the field boundary
    are used as references, so every ring drawn around a reference stays
    inside the field (simple border exclusion -- fine for a toy field, not a
    substitute for the translation-corrected estimator in
    hmt.simulate.measure_pair_correlation).

    Returns:
        r: (M,) bin centers (nm)
        g: (M,) g_ab(r); NaN where no reference points had enough margin
    """
    a = np.asarray(centroids_a, dtype=float)
    b = np.asarray(centroids_b, dtype=float)

    margin = (a[:, 0] > r_max) & (a[:, 0] < field_size_nm - r_max) & \
             (a[:, 1] > r_max) & (a[:, 1] < field_size_nm - r_max)
    refs = a[margin]

    radii = np.arange(dr, r_max + dr, dr)
    r_centers = radii - dr / 2
    g = np.full(len(radii), np.nan)
    if len(refs) == 0:
        return r_centers, g

    rho_b = len(b) / field_size_nm ** 2
    tree_b = cKDTree(b)

    prev_counts = np.zeros(len(refs))
    for i, r in enumerate(radii):
        counts = np.array(tree_b.query_ball_point(refs, r, return_length=True), dtype=float)
        ring_counts = counts - prev_counts
        prev_counts = counts
        ring_area = np.pi * (r ** 2 - (r - dr) ** 2)
        expected = rho_b * ring_area
        g[i] = np.mean(ring_counts) / expected if expected > 0 else np.nan

    return r_centers, g


def self_nonself_contact_ratio(me3_centroids, ac_centroids, field_size_nm, r_max=800.0, dr=25.0):
    """
    For each point, counts same-mark ("self") and other-mark ("nonself")
    neighbours within a growing disk of radius r (cumulative, not a ring).
    Those counts are summed across all reference points of a given type
    BEFORE taking the ratio -- sum(self_count(r)) / sum(nonself_count(r)) --
    rather than averaging each point's own ratio. This matters: at small r
    many points have zero nonself neighbours, and averaging per-point ratios
    would have to drop those (0/0 undefined), which silently changes WHICH
    points are represented in the average as r grows -- points that happen to
    already have a close nonself neighbour (e.g. paired/integrated points)
    dominate the small-r average, then get diluted as more points pick up
    their first nonself neighbour at larger r, producing a spurious upward
    drift that has nothing to do with the real spatial structure. Summing
    counts first avoids that: every reference point contributes at every r
    (contributing 0 where appropriate), so the curve reflects only the actual
    change in local self/nonself density with r.

    Traces this out as a curve over r, computed separately over me3-centered
    and ac-centered points -- the two marks need not cluster the same way
    (e.g. one mark tightly self-clustered while the other is more dispersed),
    so the two curves can differ.

    Returns both the raw count ratio and a density-normalized version:
    raw ratio(r) > 1 means self neighbours outnumber nonself neighbours within
    r; but the raw curve is also confounded by which mark simply has more
    points overall -- e.g. 2x as many me3 points as ac pushes the me3 curve up
    and the ac curve down even with zero real spatial structure, converging to
    each mark's total count ratio as r covers the whole field. The normalized
    curve divides by the field-wide self/nonself DENSITY ratio (self count
    excludes the reference point itself) so that normalized ratio(r) = 1 is
    the CSR/chance baseline regardless of that abundance imbalance -- values
    above/below 1 there reflect only clustering structure, not raw counts.

    Args:
        me3_centroids, ac_centroids: (N, 2) arrays of [x, y] in nm
        field_size_nm: side length of the square field (nm)
        r_max, dr:     max disk radius and step (nm) for the r sweep

    Returns:
        dict with keys:
            'r':                (M,) disk radius at each step (nm)
            'me3_ratio_raw':    (M,) summed self/nonself count ratio, me3-centered
            'me3_ratio_norm':   (M,) me3_ratio_raw / me3_expected_ratio
            'ac_ratio_raw', 'ac_ratio_norm': same, ac-centered
            'me3_expected_ratio', 'ac_expected_ratio': the raw-ratio CSR/chance
                                baseline for each curve (used to normalize, and
                                to draw a reference line on the raw plot)
    """
    me3 = np.asarray(me3_centroids, dtype=float)
    ac = np.asarray(ac_centroids, dtype=float)
    radii = np.arange(dr, r_max + dr, dr)

    def _curve(self_pts, nonself_pts):
        margin = (self_pts[:, 0] > r_max) & (self_pts[:, 0] < field_size_nm - r_max) & \
                 (self_pts[:, 1] > r_max) & (self_pts[:, 1] < field_size_nm - r_max)
        refs = self_pts[margin]
        if len(refs) == 0 or len(nonself_pts) == 0:
            return np.full(len(radii), np.nan), np.nan

        tree_self = cKDTree(self_pts)
        tree_nonself = cKDTree(nonself_pts)

        rho_self = (len(self_pts) - 1) / field_size_nm ** 2
        rho_nonself = len(nonself_pts) / field_size_nm ** 2
        expected_ratio = rho_self / rho_nonself

        raw_curve = np.empty(len(radii))
        for i, r in enumerate(radii):
            self_counts = np.array(
                tree_self.query_ball_point(refs, r, return_length=True), dtype=float) - 1
            nonself_counts = np.array(
                tree_nonself.query_ball_point(refs, r, return_length=True), dtype=float)

            total_nonself = nonself_counts.sum()
            raw_curve[i] = self_counts.sum() / total_nonself if total_nonself > 0 else np.nan

        return raw_curve, expected_ratio

    me3_raw, me3_expected = _curve(me3, ac)
    ac_raw, ac_expected = _curve(ac, me3)

    return {
        'r': radii,
        'me3_ratio_raw': me3_raw, 'me3_ratio_norm': me3_raw / me3_expected,
        'ac_ratio_raw': ac_raw, 'ac_ratio_norm': ac_raw / ac_expected,
        'me3_expected_ratio': me3_expected, 'ac_expected_ratio': ac_expected,
    }


def ripleys_k_cross(centroids_a, centroids_b, field_size_nm, r_max=800.0, dr=25.0):
    """
    Bivariate (cross-type) Ripley's K function K_ab(r): for a reference point
    of type A, the (border-corrected) mean count of type B points within a
    disk of radius r, scaled by area/n_b so that K_ab(r) = pi*r^2 under
    complete spatial randomness (CSR) between the two types -- independent of
    how many points each type has. This is the disk-cumulative sibling of
    cross_pair_correlation (K_ab is its running integral:
    K_ab(r) = 2*pi * integral_0^r g_ab(s) s ds), and the standard tool in
    spatial statistics for testing/quantifying attraction or repulsion
    between two point types -- summing counts before any division (like
    self_nonself_contact_ratio) rather than averaging per-point statistics
    keeps it well-behaved at low counts.

    Interpretation, via the variance-stabilized L-function
    (L_ab(r) = sqrt(K_ab(r)/pi), CSR = r):
        L_ab(r) - r > 0  -> more B neighbours near A than expected by chance
                            (attraction / integration at scale r)
        L_ab(r) - r < 0  -> fewer B neighbours near A than expected
                            (repulsion / segregation at scale r)
        L_ab(r) - r ~= 0 -> no cross-type spatial structure at scale r

    Computed in both directions (A-centered and B-centered): the
    border-exclusion correction below uses a different reference-point set
    for each direction, so the two agree only in expectation for a
    stationary process and can differ slightly at finite sample size --
    both are reported, mirroring self_nonself_contact_ratio's
    me3-centered/ac-centered convention.

    Edge correction: only reference points at least r_max from the field
    boundary are used, so every disk drawn around a reference stays inside
    the field (same border-exclusion approach as cross_pair_correlation /
    self_nonself_contact_ratio -- simpler than the translation-corrected
    estimator in hmt.simulate.measure_pair_correlation, which needs an
    explicit nucleus mask this toy field doesn't have).

    Args:
        centroids_a, centroids_b: (N, 2) arrays of [x, y] in nm
        field_size_nm: side length of the square field (nm)
        r_max, dr:     max radius and step (nm) for the r sweep

    Returns:
        dict with keys:
            'r':            (M,) radii (nm)
            'K_ab', 'L_ab': (M,) A-centered cross-K/L (mean B-neighbour count
                             around A points, scaled)
            'K_ba', 'L_ba': (M,) B-centered cross-K/L (mean A-neighbour count
                             around B points, scaled)
    """
    a = np.asarray(centroids_a, dtype=float)
    b = np.asarray(centroids_b, dtype=float)
    radii = np.arange(dr, r_max + dr, dr)
    area = field_size_nm ** 2

    def _K(ref_pts, other_pts):
        margin = (ref_pts[:, 0] > r_max) & (ref_pts[:, 0] < field_size_nm - r_max) & \
                 (ref_pts[:, 1] > r_max) & (ref_pts[:, 1] < field_size_nm - r_max)
        refs = ref_pts[margin]
        if len(refs) == 0 or len(other_pts) == 0:
            return np.full(len(radii), np.nan)

        tree_other = cKDTree(other_pts)
        K = np.empty(len(radii))
        for i, r in enumerate(radii):
            counts = np.array(tree_other.query_ball_point(refs, r, return_length=True), dtype=float)
            K[i] = (area / len(other_pts)) * np.mean(counts)
        return K

    K_ab = _K(a, b)
    K_ba = _K(b, a)
    L_ab = np.sqrt(np.clip(K_ab, 0, None) / np.pi)
    L_ba = np.sqrt(np.clip(K_ba, 0, None) / np.pi)

    return {'r': radii, 'K_ab': K_ab, 'L_ab': L_ab, 'K_ba': K_ba, 'L_ba': L_ba}


# --- Graph-theory metrics ---------------------------------------------------

def _delaunay_edges(points):
    """Unique undirected edges (index pairs into `points`) of the Delaunay triangulation over `points`."""
    tri = Delaunay(points)
    edge_set = set()
    for simplex in tri.simplices:
        for i in range(3):
            a, b = int(simplex[i]), int(simplex[(i + 1) % 3])
            edge_set.add((a, b) if a < b else (b, a))
    return np.array(sorted(edge_set))


def delaunay_channel_mixing(me3_centroids, ac_centroids):
    """
    Builds a Delaunay triangulation over the POOLED me3+ac centroids -- a
    parameter-free "who is whose spatial neighbour" graph, unlike every other
    metric in this module, which all need a radius r -- and asks how often an
    edge connects two points of the SAME channel (homotypic) vs DIFFERENT
    channels (heterotypic).

    Reports two related numbers:
      heterotypic_fraction: raw fraction of triangulation edges that are
          heterotypic. Confounded by the me3/ac abundance ratio -- same
          caveat as self_nonself_contact_ratio's raw curves: with very
          unequal channel counts, heterotypic edges are rarer/more common
          just from there being fewer/more of one type, regardless of real
          mixing.
      assortativity: Newman's categorical assortativity coefficient on the
          channel label, which corrects for that abundance imbalance -- the
          graph-theory analog of self_nonself_contact_ratio's *_norm curves
          or Ripley's L(r) - r.

    assortativity in [-1, 1]:
        +1 -> every edge connects same-channel points (fully segregated)
         0 -> edges connect channels no more/less than expected from their
              relative abundance alone (no mixing preference)
        -1 -> every edge connects different-channel points (maximally mixed,
              "checkerboard" alternation)

    Args:
        me3_centroids, ac_centroids: (N, 2) arrays of [x, y] in nm

    Returns:
        dict with keys:
            'points', 'labels': pooled [x, y] array (me3 then ac) and their
                                channel labels (0=me3, 1=ac) -- together with
                                'edges', enough to redraw the graph
            'edges':            (E, 2) int array of point indices (into
                                'points') for each triangulation edge
            'n_edges', 'n_heterotypic': edge counts
            'heterotypic_fraction': raw fraction of edges that are heterotypic
            'assortativity':        Newman's categorical assortativity coefficient
    """
    me3 = np.asarray(me3_centroids, dtype=float)
    ac = np.asarray(ac_centroids, dtype=float)
    points = np.vstack([me3, ac])
    labels = np.concatenate([np.zeros(len(me3), dtype=int), np.ones(len(ac), dtype=int)])

    edges = _delaunay_edges(points)
    n_edges = len(edges)
    edge_labels = labels[edges]
    heterotypic = edge_labels[:, 0] != edge_labels[:, 1]
    n_heterotypic = int(heterotypic.sum())
    heterotypic_fraction = n_heterotypic / n_edges if n_edges > 0 else np.nan

    # Newman's categorical assortativity: symmetric 2x2 mixing matrix over
    # edge endpoint types, each edge contributing to both (i, j) and (j, i)
    # so the matrix is symmetric regardless of edge storage order.
    M = np.zeros((2, 2))
    for l0, l1 in edge_labels:
        M[l0, l1] += 1
        M[l1, l0] += 1
    e = M / M.sum() if M.sum() > 0 else M
    a_marg = e.sum(axis=1)
    denom = 1.0 - np.sum(a_marg ** 2)
    assortativity = (np.trace(e) - np.sum(a_marg ** 2)) / denom if denom > 0 else np.nan

    return {
        'points': points, 'labels': labels, 'edges': edges,
        'n_edges': n_edges, 'n_heterotypic': n_heterotypic,
        'heterotypic_fraction': heterotypic_fraction, 'assortativity': assortativity,
    }


def friedman_rafsky_test(me3_centroids, ac_centroids, n_permutations=500, rng=None):
    """
    Friedman-Rafsky two-sample test: builds a minimum spanning tree (MST)
    over the POOLED me3+ac centroids and counts how many MST edges connect a
    me3 point to an ac point ("cross-type" edges). Under the null hypothesis
    that me3 and ac centroids are drawn from the same spatial distribution
    (no channel-specific clustering), the channel LABELS are exchangeable
    across the fixed set of pooled positions -- so the null distribution of
    the cross-type edge count is obtained by permuting labels over the SAME
    MST (positions and tree topology don't change, only which point is
    called me3 vs ac). This is a Monte Carlo approximation to Friedman &
    Rafsky's (1979) exact permutation null -- equivalent to their closed-form
    normal approximation, without needing to re-derive its variance formula.

    Few cross-type edges (observed << null) -> the two channels form
    separate, self-clustered branches of the tree (segregated). Many
    cross-type edges (observed >> null) -> me3 and ac points are interleaved
    throughout the tree (integrated/mixed) -- possibly even more mixed than
    chance (e.g. tightly paired points that are each other's nearest
    neighbour, so the tree keeps stitching them together).

    Args:
        me3_centroids, ac_centroids: (N, 2) arrays of [x, y] in nm
        n_permutations: number of label permutations for the null distribution
        rng:            np.random.default_rng() instance; None makes a fresh one

    Returns:
        dict with keys:
            'points', 'labels': as in delaunay_channel_mixing
            'edges':            (N-1, 2) int array of point indices (into
                                'points') for each MST edge
            'observed_cross_edges', 'total_edges': counts
            'observed_fraction':   observed_cross_edges / total_edges
            'null_mean_fraction', 'null_std_fraction': permutation-null moments,
                                same units as observed_fraction
            'z_score':             (observed - null_mean) / null_std
            'p_value':             two-sided empirical permutation p-value
    """
    rng = rng or np.random.default_rng()
    me3 = np.asarray(me3_centroids, dtype=float)
    ac = np.asarray(ac_centroids, dtype=float)
    points = np.vstack([me3, ac])
    labels = np.concatenate([np.zeros(len(me3), dtype=int), np.ones(len(ac), dtype=int)])

    dist = squareform(pdist(points))
    mst = minimum_spanning_tree(dist).tocoo()
    edges = np.column_stack([mst.row, mst.col])
    total_edges = len(edges)

    def _cross_count(lab):
        return int(np.sum(lab[edges[:, 0]] != lab[edges[:, 1]]))

    observed = _cross_count(labels)

    null_counts = np.empty(n_permutations, dtype=int)
    perm_labels = labels.copy()
    for i in range(n_permutations):
        rng.shuffle(perm_labels)
        null_counts[i] = _cross_count(perm_labels)

    null_mean = float(null_counts.mean())
    null_std = float(null_counts.std(ddof=1))
    z_score = (observed - null_mean) / null_std if null_std > 0 else np.nan
    p_value = float(np.mean(np.abs(null_counts - null_mean) >= abs(observed - null_mean)))

    return {
        'points': points, 'labels': labels, 'edges': edges,
        'observed_cross_edges': observed, 'total_edges': total_edges,
        'observed_fraction': observed / total_edges if total_edges > 0 else np.nan,
        'null_mean_fraction': null_mean / total_edges if total_edges > 0 else np.nan,
        'null_std_fraction': null_std / total_edges if total_edges > 0 else np.nan,
        'z_score': z_score, 'p_value': p_value,
    }


def mst_merge_curve(me3_centroids, ac_centroids, r_max=800.0):
    """
    H0 persistent homology ("merge curve") over the pooled me3+ac centroids:
    the same minimum spanning tree used by friedman_rafsky_test, replayed in
    increasing edge-length order. Sorting an MST's own edges by length
    exactly reproduces the order components merge as a Vietoris-Rips
    filtration radius grows -- Kruskal's algorithm run on ONLY the MST's own
    edges can never hit a cycle, since they're already a tree, so "add edges
    shortest-first" and "the MST's edges in length order" are the same
    sequence. This is the standard fact that a point cloud's 0-dimensional
    persistent homology equals its single-linkage merge tree, so no separate
    TDA library (or hand-rolled union-find) is needed -- the MST already
    computed for friedman_rafsky_test contains the full H0 barcode.

    Unlike friedman_rafsky_test, which collapses the whole tree into one
    number (total cross-type edge count vs. permutation null), this keeps
    the RADIUS each merge happened at -- answering "at what length scale
    does cross-channel merging start dominating" rather than "how much
    cross-channel merging is there overall". A merge event is heterotypic if
    its edge directly connects a me3 point to an ac point.

    Args:
        me3_centroids, ac_centroids: (N, 2) arrays of [x, y] in nm
        r_max: merge events beyond this radius are dropped -- the last few
               merges of any MST span most of the field and are dominated by
               whichever two components happened to be left, not
               nanodomain-scale structure (matches the r_max convention used
               by every radius-based metric elsewhere in this module)

    Returns:
        dict with keys:
            'radius':               (M,) merge radius (nm), increasing, <= r_max
            'heterotypic':          (M,) bool, whether that merge event's
                                    edge directly connects a me3 point to an
                                    ac point
            'cum_heterotypic_frac': (M,) cumulative fraction of merges up to
                                    and including this radius that were
                                    heterotypic
            'chance_fraction':      the label-permutation-null EXPECTED
                                    heterotypic fraction (same quantity as
                                    friedman_rafsky_test's null_mean_fraction,
                                    computed here in closed form since it is
                                    the same for every edge by exchangeability,
                                    so it doesn't vary with radius):
                                    2*n_me3*n_ac / (N*(N-1))
    """
    me3 = np.asarray(me3_centroids, dtype=float)
    ac = np.asarray(ac_centroids, dtype=float)
    points = np.vstack([me3, ac])
    labels = np.concatenate([np.zeros(len(me3), dtype=int), np.ones(len(ac), dtype=int)])
    n = len(points)

    dist = squareform(pdist(points))
    mst = minimum_spanning_tree(dist).tocoo()
    edges = np.column_stack([mst.row, mst.col])
    lengths = dist[edges[:, 0], edges[:, 1]]

    order = np.argsort(lengths)
    edges, lengths = edges[order], lengths[order]

    keep = lengths <= r_max
    edges, lengths = edges[keep], lengths[keep]

    heterotypic = labels[edges[:, 0]] != labels[edges[:, 1]]
    cum_heterotypic_frac = np.cumsum(heterotypic) / np.arange(1, len(heterotypic) + 1)

    n_me3, n_ac = len(me3), len(ac)
    chance_fraction = 2 * n_me3 * n_ac / (n * (n - 1))

    return {
        'radius': lengths, 'heterotypic': heterotypic,
        'cum_heterotypic_frac': cum_heterotypic_frac, 'chance_fraction': chance_fraction,
    }


def _segment_crossing_matrix(points_a, edges_a, points_b, edges_b):
    """
    Boolean (E_a, E_b) matrix: whether edge i of mesh A geometrically crosses
    edge j of mesh B (a proper intersection between the two line segments,
    via the standard pairwise cross-product orientation test, fully
    vectorised over all edge pairs at once). Collinear/touching edge cases
    are ignored -- they have probability zero for continuous random point
    positions, so not worth the extra special-casing here.
    """
    def cross2(o, p, q):
        return (p[..., 0] - o[..., 0]) * (q[..., 1] - o[..., 1]) - (p[..., 1] - o[..., 1]) * (q[..., 0] - o[..., 0])

    P1 = points_a[edges_a[:, 0]][:, None, :]
    P2 = points_a[edges_a[:, 1]][:, None, :]
    Q1 = points_b[edges_b[:, 0]][None, :, :]
    Q2 = points_b[edges_b[:, 1]][None, :, :]

    d1 = cross2(Q1, Q2, P1)
    d2 = cross2(Q1, Q2, P2)
    d3 = cross2(P1, P2, Q1)
    d4 = cross2(P1, P2, Q2)

    return ((d1 > 0) != (d2 > 0)) & ((d3 > 0) != (d4 > 0))


def mesh_overlap(me3_centroids, ac_centroids, n_permutations=0, rng=None):
    """
    Builds two SEPARATE Delaunay triangulations -- one over me3 centroids
    only, one over ac centroids only (unlike delaunay_channel_mixing, which
    pools both channels into a single mesh) -- and asks where the two
    independent meshes physically cross in space: a me3 edge and an ac edge
    "overlap" if the two line segments actually intersect, not just if their
    endpoints happen to be close.

    This is a genuinely different notion of integration than every other
    metric in this module: not "is a me3 domain's nearest neighbour an ac
    domain" (colocalization_fraction, delaunay_channel_mixing), but "does
    the local CONNECTIVITY STRUCTURE (the mesh skeleton) of one channel
    spatially thread through the other's". Two channels can have their
    domains interleaved at the point level yet have meshes that barely
    cross if both are locally sparse, or vice versa -- not simply a
    rescaling of the point-based metrics.

    Args:
        me3_centroids, ac_centroids: (N, 2) arrays of [x, y] in nm
        n_permutations: if > 0, also runs a permutation-null significance
                        test (see "Null model" below); 0 (default) skips it
                        and returns crossing counts/fractions only. Each
                        permutation rebuilds BOTH triangulations from
                        scratch (unlike friedman_rafsky_test's null, which
                        reuses one fixed tree and only reshuffles labels),
                        so this is noticeably slower per permutation --
                        kept off by default for that reason.
        rng:            np.random.default_rng() instance; None makes a
                        fresh one (only used if n_permutations > 0)

    Null model (when n_permutations > 0): pool all me3+ac points, randomly
    relabel them into two groups of the same sizes (n_me3, n_ac), rebuild
    BOTH triangulations from that fake split, and recount crossings. Repeat
    n_permutations times to get the crossing count expected under "same
    positions, no real channel-specific structure", then compare the real
    split's crossing count against that null distribution (z-score,
    p-value) -- the same label-permutation logic as friedman_rafsky_test,
    just applied to two rebuilt triangulations instead of one fixed tree.

    Returns:
        dict with keys:
            'me3_points', 'ac_points':   the two input arrays, unchanged
            'me3_edges', 'ac_edges':     (E1, 2) / (E2, 2) int arrays of
                                         point indices (into their own
                                         me3_points / ac_points) per mesh edge
            'crossing_points':           (n_crossings, 2) [x, y] coordinates
                                         where a me3 edge and an ac edge
                                         actually cross, for plotting
            'n_crossings':               total number of (me3 edge, ac edge)
                                         pairs that physically cross
            'me3_crossing_fraction':     fraction of me3 edges that cross at
                                         least one ac edge
            'ac_crossing_fraction':      fraction of ac edges that cross at
                                         least one me3 edge
            'null_mean_crossings', 'null_std_crossings', 'z_score', 'p_value':
                                         only present if n_permutations > 0
    """
    me3 = np.asarray(me3_centroids, dtype=float)
    ac = np.asarray(ac_centroids, dtype=float)

    me3_edges = _delaunay_edges(me3)
    ac_edges = _delaunay_edges(ac)
    crossing = _segment_crossing_matrix(me3, me3_edges, ac, ac_edges)
    n_crossings = int(crossing.sum())

    crossing_idx = np.argwhere(crossing)
    if len(crossing_idx) > 0:
        p1 = me3[me3_edges[crossing_idx[:, 0], 0]]
        p2 = me3[me3_edges[crossing_idx[:, 0], 1]]
        q1 = ac[ac_edges[crossing_idx[:, 1], 0]]
        q2 = ac[ac_edges[crossing_idx[:, 1], 1]]
        dP = p2 - p1
        dQ = q2 - q1
        denom = dP[:, 0] * dQ[:, 1] - dP[:, 1] * dQ[:, 0]
        t = ((q1[:, 0] - p1[:, 0]) * dQ[:, 1] - (q1[:, 1] - p1[:, 1]) * dQ[:, 0]) / denom
        crossing_points = p1 + t[:, None] * dP
    else:
        crossing_points = np.empty((0, 2))

    result = {
        'me3_points': me3, 'ac_points': ac,
        'me3_edges': me3_edges, 'ac_edges': ac_edges,
        'crossing_points': crossing_points, 'n_crossings': n_crossings,
        'me3_crossing_fraction': float(crossing.any(axis=1).mean()) if len(me3_edges) > 0 else np.nan,
        'ac_crossing_fraction': float(crossing.any(axis=0).mean()) if len(ac_edges) > 0 else np.nan,
    }

    if n_permutations > 0:
        rng = rng or np.random.default_rng()
        points = np.vstack([me3, ac])
        n_me3, n = len(me3), len(points)
        idx = np.arange(n)

        null_counts = np.empty(n_permutations, dtype=int)
        for i in range(n_permutations):
            rng.shuffle(idx)
            fake_me3, fake_ac = points[idx[:n_me3]], points[idx[n_me3:]]
            fake_crossing = _segment_crossing_matrix(
                fake_me3, _delaunay_edges(fake_me3), fake_ac, _delaunay_edges(fake_ac))
            null_counts[i] = fake_crossing.sum()

        null_mean = float(null_counts.mean())
        null_std = float(null_counts.std(ddof=1))
        result['null_mean_crossings'] = null_mean
        result['null_std_crossings'] = null_std
        result['z_score'] = (n_crossings - null_mean) / null_std if null_std > 0 else np.nan
        result['p_value'] = float(np.mean(np.abs(null_counts - null_mean) >= abs(n_crossings - null_mean)))

    return result


# --- Plotting --------------------------------------------------------------

def plot_sncr_sweep(sweep, field_size_nm, r_max=800.0, dr=25.0, title=None):
    """
    Computes self_nonself_contact_ratio for every level in `sweep` (as returned
    by generate_integration_sweep) and plots the raw and density-normalized
    me3-centered/ac-centered curves in a 2x2 grid, one color per integration_level.

    Args:
        sweep:      list of dicts from generate_integration_sweep, each with
                    'integration_level', 'me3_seeds', 'ac_seeds'
        field_size_nm, r_max, dr: forwarded to self_nonself_contact_ratio
        title:      figure suptitle; a sensible default is used if None

    Returns:
        (fig, axes): the created matplotlib Figure and 2x2 Axes array
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    cmap = plt.get_cmap("viridis")
    me3_expected = ac_expected = None

    for row in sweep:
        level = row["integration_level"]
        ratios = self_nonself_contact_ratio(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            field_size_nm=field_size_nm, r_max=r_max, dr=dr)
        me3_expected, ac_expected = ratios["me3_expected_ratio"], ratios["ac_expected_ratio"]

        color = cmap(level)
        axes[0, 0].plot(ratios["r"], ratios["me3_ratio_raw"], marker="o", markersize=3, color=color, label=f"level={level}")
        axes[0, 1].plot(ratios["r"], ratios["ac_ratio_raw"], marker="o", markersize=3, color=color, label=f"level={level}")
        axes[1, 0].plot(ratios["r"], ratios["me3_ratio_norm"], marker="o", markersize=3, color=color, label=f"level={level}")
        axes[1, 1].plot(ratios["r"], ratios["ac_ratio_norm"], marker="o", markersize=3, color=color, label=f"level={level}")

    axes[0, 0].axhline(me3_expected, color="gray", linestyle="--", label=f"CSR (chance, raw={me3_expected:.2f})")
    axes[0, 1].axhline(ac_expected, color="gray", linestyle="--", label=f"CSR (chance, raw={ac_expected:.2f})")
    axes[1, 0].axhline(1.0, color="gray", linestyle="--", label="CSR (chance)")
    axes[1, 1].axhline(1.0, color="gray", linestyle="--", label="CSR (chance)")

    axes[0, 0].set_title("me3-centered (raw)")
    axes[0, 1].set_title("ac-centered (raw)")
    axes[1, 0].set_title("me3-centered (normalized)")
    axes[1, 1].set_title("ac-centered (normalized)")

    axes[0, 0].set_ylabel("self/nonself count ratio")
    axes[1, 0].set_ylabel("normalized self/nonself ratio")
    for ax in axes[1, :]:
        ax.set_xlabel("r (nm)")
    for ax in axes.ravel():
        ax.legend(fontsize=7)

    plt.suptitle(title or "Self-nonself contact ratio")
    plt.tight_layout()
    plt.show()
    return fig, axes


def plot_ripley_k_sweep(sweep, field_size_nm, r_max=800.0, dr=25.0, title=None):
    """
    Computes ripleys_k_cross for every level in `sweep` (as returned by
    generate_integration_sweep) and plots the detrended L_ab(r) - r /
    L_ba(r) - r curves, one color per integration_level.

    Args:
        sweep:      list of dicts from generate_integration_sweep, each with
                    'integration_level', 'me3_seeds', 'ac_seeds'
        field_size_nm, r_max, dr: forwarded to ripleys_k_cross
        title:      figure suptitle; a sensible default is used if None

    Returns:
        (fig, axes): the created matplotlib Figure and 1x2 Axes array
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    cmap = plt.get_cmap("viridis")

    for row in sweep:
        level = row["integration_level"]
        k = ripleys_k_cross(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            field_size_nm=field_size_nm, r_max=r_max, dr=dr)

        color = cmap(level)
        axes[0].plot(k["r"], k["L_ab"] - k["r"], marker="o", markersize=3, color=color, label=f"level={level}")
        axes[1].plot(k["r"], k["L_ba"] - k["r"], marker="o", markersize=3, color=color, label=f"level={level}")

    axes[0].axhline(0.0, color="gray", linestyle="--", label="CSR (no association)")
    axes[1].axhline(0.0, color="gray", linestyle="--", label="CSR (no association)")

    axes[0].set_title("me3-centered")
    axes[1].set_title("ac-centered")
    axes[0].set_ylabel("L(r) - r (nm)")
    for ax in axes:
        ax.set_xlabel("r (nm)")
        ax.legend(fontsize=7)

    plt.suptitle(title or "Ripley's K/L")
    plt.tight_layout()
    plt.show()
    return fig, axes


def _draw_graph(ax, points, edges, labels, heterotypic_color="darkorange", homotypic_color="lightgray",
                heterotypic_width=1.2, homotypic_width=0.6):
    """
    Draws a point+edge graph (Delaunay triangulation or MST) on `ax`: nodes
    colored by channel (green=me3, red=ac), edges colored by whether they
    connect the same channel (homotypic, faint gray) or different channels
    (heterotypic, highlighted). Shared by plot_delaunay_sweep and
    plot_mst_sweep since both draw the same kind of point+edge graph.
    """
    heterotypic = labels[edges[:, 0]] != labels[edges[:, 1]]
    segments = points[edges]  # (E, 2, 2): E edges, 2 endpoints, [x, y] each

    ax.add_collection(LineCollection(segments[~heterotypic], colors=homotypic_color,
                                     linewidths=homotypic_width, zorder=1))
    ax.add_collection(LineCollection(segments[heterotypic], colors=heterotypic_color,
                                     linewidths=heterotypic_width, zorder=2))
    ax.scatter(points[labels == 0, 0], points[labels == 0, 1], s=15, color="green", zorder=3, label="me3")
    ax.scatter(points[labels == 1, 0], points[labels == 1, 1], s=15, color="red", zorder=3, label="ac")
    ax.set_aspect("equal")


def plot_delaunay_sweep(sweep, title=None):
    """
    For every level in `sweep` (from generate_integration_sweep), draws the
    Delaunay triangulation over the pooled me3+ac centroids side by side
    across integration_level (heterotypic edges highlighted, homotypic edges
    faint), then plots heterotypic_fraction (raw) and assortativity
    (chance-corrected) vs. integration_level underneath -- see
    delaunay_channel_mixing for what each measures.

    Args:
        sweep: list of dicts from generate_integration_sweep
        title: suptitle for the graph row; a sensible default is used if None

    Returns:
        list of dicts from delaunay_channel_mixing, one per level (sweep order)
    """
    results = [
        delaunay_channel_mixing(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy())
        for row in sweep
    ]

    fig, axes = plt.subplots(1, len(sweep), figsize=(4 * len(sweep), 4), sharey=True)
    axes = np.atleast_1d(axes)
    for row, result, ax in zip(sweep, results, axes):
        _draw_graph(ax, result["points"], result["edges"], result["labels"])
        ax.set_title(f"level={row['integration_level']}\nheterotypic={result['heterotypic_fraction']:.2f}")

    axes[0].legend(loc="upper right", fontsize=9)
    plt.suptitle(title or "Delaunay triangulation, channel mixing")
    plt.tight_layout()
    plt.show()

    levels = [row["integration_level"] for row in sweep]
    fig, axes2 = plt.subplots(1, 2, figsize=(11, 4.5))

    axes2[0].plot(levels, [r["heterotypic_fraction"] for r in results], "o-", color="darkorange")
    axes2[0].set_title("heterotypic edge fraction (raw)")
    axes2[0].set_xlabel("integration_level")
    axes2[0].set_ylim(-0.05, 1.05)

    axes2[1].plot(levels, [r["assortativity"] for r in results], "o-", color="steelblue")
    axes2[1].axhline(0.0, color="gray", linestyle="--", label="no mixing preference")
    axes2[1].set_title("assortativity (chance-corrected)")
    axes2[1].set_xlabel("integration_level")
    axes2[1].set_ylim(-1.05, 1.05)
    axes2[1].legend(fontsize=8)

    plt.tight_layout()
    plt.show()
    return results


def plot_mst_sweep(sweep, n_permutations=500, rng=None, title=None):
    """
    For every level in `sweep` (from generate_integration_sweep), draws the
    minimum spanning tree over the pooled me3+ac centroids side by side
    across integration_level (cross-type edges highlighted), then plots the
    Friedman-Rafsky cross-edge fraction (observed vs. permutation null) and
    z-score vs. integration_level underneath -- see friedman_rafsky_test for
    what each measures.

    Args:
        sweep:          list of dicts from generate_integration_sweep
        n_permutations: forwarded to friedman_rafsky_test
        rng:            np.random.default_rng() instance, reused across
                        levels; None makes a fresh one
        title:          suptitle for the tree row; a sensible default is
                        used if None

    Returns:
        list of dicts from friedman_rafsky_test, one per level (sweep order)
    """
    rng = rng or np.random.default_rng()
    results = [
        friedman_rafsky_test(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            n_permutations=n_permutations, rng=rng)
        for row in sweep
    ]

    fig, axes = plt.subplots(1, len(sweep), figsize=(4 * len(sweep), 4), sharey=True)
    axes = np.atleast_1d(axes)
    for row, result, ax in zip(sweep, results, axes):
        _draw_graph(ax, result["points"], result["edges"], result["labels"],
                    heterotypic_width=1.5, homotypic_width=0.8)
        ax.set_title(f"level={row['integration_level']}\nz={result['z_score']:.1f}")

    axes[0].legend(loc="upper right", fontsize=9)
    plt.suptitle(title or "Minimum spanning tree, channel mixing")
    plt.tight_layout()
    plt.show()

    levels = [row["integration_level"] for row in sweep]
    fig, axes2 = plt.subplots(1, 2, figsize=(11, 4.5))

    axes2[0].plot(levels, [r["observed_fraction"] for r in results], "o-", color="darkorange", label="observed")
    axes2[0].plot(levels, [r["null_mean_fraction"] for r in results], "o--", color="gray", label="permutation null")
    axes2[0].set_title("cross-type edge fraction")
    axes2[0].set_xlabel("integration_level")
    axes2[0].set_ylim(-0.05, 1.05)
    axes2[0].legend(fontsize=8)

    axes2[1].plot(levels, [r["z_score"] for r in results], "o-", color="steelblue")
    axes2[1].axhline(0.0, color="gray", linestyle="--")
    axes2[1].axhline(1.96, color="gray", linestyle=":", linewidth=1, label="|z| = 1.96")
    axes2[1].axhline(-1.96, color="gray", linestyle=":", linewidth=1)
    axes2[1].set_title("Friedman-Rafsky z-score")
    axes2[1].set_xlabel("integration_level")
    axes2[1].legend(fontsize=8)

    plt.tight_layout()
    plt.show()
    return results


def plot_merge_curve_sweep(sweep, r_max=800.0, title=None):
    """
    Computes mst_merge_curve for every level in `sweep` (from
    generate_integration_sweep) and plots the cumulative heterotypic-merge
    fraction vs. merge radius, one color per integration_level -- how EARLY
    (at what length scale) cross-channel merging starts dominating, as a
    complement to plot_mst_sweep's single aggregate number per level.

    Args:
        sweep: list of dicts from generate_integration_sweep
        r_max: forwarded to mst_merge_curve
        title: figure title; a sensible default is used if None

    Returns:
        list of dicts from mst_merge_curve, one per level (sweep order)
    """
    results = []
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("viridis")
    chance = None

    for row in sweep:
        level = row["integration_level"]
        curve = mst_merge_curve(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            r_max=r_max)
        results.append(curve)
        chance = curve["chance_fraction"]

        ax.plot(curve["radius"], curve["cum_heterotypic_frac"], color=cmap(level), label=f"level={level}")

    ax.axhline(chance, color="gray", linestyle="--", label=f"chance ({chance:.2f})")
    ax.set_xlabel("merge radius (nm)")
    ax.set_ylabel("cumulative heterotypic-merge fraction")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=8)
    ax.set_title(title or "MST merge curve (H0 persistence)")
    plt.tight_layout()
    plt.show()
    return results


def plot_mesh_overlap_sweep(sweep, n_permutations=100, rng=None, title=None):
    """
    For every level in `sweep` (from generate_integration_sweep), draws the
    SEPARATE me3 and ac Delaunay meshes overlaid (me3 edges green, ac edges
    red, actual crossing points marked with an orange x) side by side across
    integration_level, then plots the crossing fraction (raw, both
    directions) and -- if n_permutations > 0 -- the permutation z-score vs.
    integration_level underneath. See mesh_overlap for what each measures.

    Args:
        sweep:          list of dicts from generate_integration_sweep
        n_permutations: forwarded to mesh_overlap; each permutation rebuilds
                        BOTH triangulations from scratch, so this is pricier
                        per-call than friedman_rafsky_test's permutations --
                        kept modest by default. 0 skips the null entirely
                        (crossing counts/fractions only, no z-score panel)
        rng:            np.random.default_rng() instance, reused across
                        levels; None makes a fresh one
        title:          suptitle for the mesh row; a sensible default is
                        used if None

    Returns:
        list of dicts from mesh_overlap, one per level (sweep order)
    """
    rng = rng or np.random.default_rng()
    results = [
        mesh_overlap(
            row["me3_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            row["ac_seeds"][["x [nm]", "y [nm]"]].to_numpy(),
            n_permutations=n_permutations, rng=rng)
        for row in sweep
    ]
    has_null = "z_score" in results[0]

    fig, axes = plt.subplots(1, len(sweep), figsize=(4 * len(sweep), 4), sharey=True)
    axes = np.atleast_1d(axes)
    for row, result, ax in zip(sweep, results, axes):
        me3_pts, ac_pts = result["me3_points"], result["ac_points"]
        me3_segs = me3_pts[result["me3_edges"]]
        ac_segs = ac_pts[result["ac_edges"]]

        ax.add_collection(LineCollection(me3_segs, colors="green", linewidths=0.7, alpha=0.6, zorder=1))
        ax.add_collection(LineCollection(ac_segs, colors="red", linewidths=0.7, alpha=0.6, zorder=1))
        ax.scatter(me3_pts[:, 0], me3_pts[:, 1], s=12, color="green", zorder=3, label="me3")
        ax.scatter(ac_pts[:, 0], ac_pts[:, 1], s=12, color="red", zorder=3, label="ac")
        if len(result["crossing_points"]) > 0:
            ax.scatter(result["crossing_points"][:, 0], result["crossing_points"][:, 1],
                      s=4, color="darkorange", marker="x", linewidths=0.7, zorder=4, label="crossing")

        z_label = f"  z={result['z_score']:.1f}" if has_null else ""
        ax.set_title(f"level={row['integration_level']}\ncrossings={result['n_crossings']}{z_label}")
        ax.set_aspect("equal")

    mesh_title = title or "Delaunay mesh overlap, me3 vs. ac"
    axes[0].legend(loc="upper right", fontsize=8)
    plt.suptitle(mesh_title)
    plt.tight_layout()
    plt.show()

    levels = [row["integration_level"] for row in sweep]
    fig, axes2 = plt.subplots(1, 3 if has_null else 1, figsize=(5.5 * (3 if has_null else 1), 4.5))
    axes2 = np.atleast_1d(axes2)

    axes2[0].plot(levels, [r["me3_crossing_fraction"] for r in results], "o-", color="green", label="me3 edges")
    axes2[0].plot(levels, [r["ac_crossing_fraction"] for r in results], "o-", color="red", label="ac edges")
    axes2[0].set_title("fraction of edges crossing the other mesh")
    axes2[0].set_xlabel("integration_level")
    axes2[0].set_ylim(-0.05, 1.05)
    axes2[0].legend(fontsize=8)

    if has_null:
        axes2[1].plot(levels, [r["n_crossings"] for r in results], "o-", color="darkorange", label="observed")
        axes2[1].plot(levels, [r["null_mean_crossings"] for r in results], "o--", color="gray", label="permutation null")
        axes2[1].set_title("total edge crossings")
        axes2[1].set_xlabel("integration_level")
        axes2[1].legend(fontsize=8)

        axes2[2].plot(levels, [r["z_score"] for r in results], "o-", color="steelblue")
        axes2[2].axhline(0.0, color="gray", linestyle="--")
        axes2[2].axhline(1.96, color="gray", linestyle=":", linewidth=1, label="|z| = 1.96")
        axes2[2].axhline(-1.96, color="gray", linestyle=":", linewidth=1)
        axes2[2].set_title("mesh-overlap z-score")
        axes2[2].set_xlabel("integration_level")
        axes2[2].legend(fontsize=8)

    plt.suptitle(f"{mesh_title} -- crossing statistics")
    plt.tight_layout()
    plt.show()
    return results
