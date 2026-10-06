"""
Visualization of SR-EV structures -- the monomer cloud, the antibody /
fluorophore placement chain, and identified packing domains. Replaces the
legacy pipeline's OVITO-based inspection step, and feeds OVITO again rather
than a molecular viewer: this is particle data, and OVITO is a particle
visualizer.

Three views of the same structure, all fed by one shared assembly step
(build_particle_table), so every view agrees on species, colors, radii and
units:

  view_napari(...)        Interactive 3D, driven entirely from Python, no
                          separate application. Default choice while
                          working -- napari is already this codebase's
                          viewer (see hmt_v3/visualize.py) and stays
                          responsive at 10^5-10^6 points. Needs a real
                          OpenGL context, so it will not run headless.
  plot_structure_3d(...)  Static matplotlib figure that renders INLINE in a
                          notebook. Good for a quick look or a figure
                          panel; auto-thins large clouds.
  export_for_ovito(...)   Writes extended XYZ for OVITO, for large scenes,
                          publication rendering, and OVITO's own analysis
                          modifiers.

WHY EXTENDED XYZ: a raw config-N.txt is already plain XYZ (count / comment /
"type x y z"), so OVITO can open one directly -- but only as untyped points
in the generator's own length units, with no compartment or coordination
information. This writer emits the same format with an extended `Properties`
header, which is what lets a single file carry nanometre positions, per-
particle color, physical radii, and the scalar fields worth coloring by.
Everything below is a per-particle property; there is no molecular topology
being faked (no residues, chains, or B-factor columns standing in for
physical quantities), and no 99999-atom or coordinate-field limit.

Properties written, and what OVITO does with each:
    species        particle TYPE name -- becomes a named type in OVITO's
                   type list, so each population can be colored, sized or
                   hidden independently. The two dye populations get
                   distinct types (DYE_HET / DYE_EUC), which is what makes
                   the two fluorophore species separable in one scene.
    pos            position, in NANOMETRES (see UNITS below)
    Color          per-particle RGB, so the scene is already correctly
                   colored on open with no modifier set up. NOTE that
                   OVITO gives an explicit Color property priority over
                   its per-particle-type colors, so while this is present
                   the type list's colors have no visible effect -- pass
                   include_color=False to omit it and let OVITO color by
                   type instead (a Color Coding modifier overrides it
                   either way, since it runs later in the pipeline)
    Radius         per-particle display radius in nm -- monomers draw at
                   r° = 4.9 nm, domains at their measured radius_nm
    Transparency   used to make domain spheres see-through, so the beads
                   inside them stay visible
    Coordination   integer CN, mapping onto the same property OVITO's own
                   Coordination Analysis modifier produces -- so Color
                   Coding lists it exactly as if OVITO had computed it, and
                   OVITO's version can be used to independently check
                   compartments.compute_coordination_number
    compartment    0 unknown / 1 heterochromatin / 2 euchromatin, numeric
                   so it can drive Color Coding
    layer          0 monomers / 1 primaries / 2 secondaries /
                   3 fluorophores / 4 domains

UNITS: positions are written in nanometres, matching every other number
this package reports, so distances measured in OVITO read directly in nm.
(Opening a raw config-N.txt instead gives the generator's native Umin units,
which are 10x smaller -- see io.parse_config_txt.)

NOTE on pipeline outputs: pipeline.run_pipeline shifts every coordinate
table into frame-relative space before rendering, so tables taken from its
return value are in that shifted frame (relative geometry is unchanged;
absolute positions differ from the source config-N.txt by
run_config.json's "coordinate_offset_nm").
"""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import Normalize

_XYZ_COLS = ["x [nm]", "y [nm]", "z [nm]"]

_LAYER_ORDER = ("monomers", "primaries", "secondaries", "fluorophores", "domains")
_LAYER_CODE = {name: i for i, name in enumerate(_LAYER_ORDER)}

_COMPARTMENT_CODE = {"none": 0, "heterochromatin": 1, "euchromatin": 2}
_COMPARTMENT_SHORT = {"heterochromatin": "HET", "euchromatin": "EUC"}

# Base species names; a channel suffix is appended for the label layers when
# the channel is known, so a two-channel scene has separately selectable
# DYE_HET / DYE_EUC types rather than one undifferentiated dye population.
_LAYER_SPECIES = {
    "monomers": "MONO",
    "primaries": "PRI",
    "secondaries": "SEC",
    "fluorophores": "DYE",
    "domains": "DOMAIN",
}

_LAYER_DISPLAY = {
    "monomers":     {"name": "Monomers",             "color": (0.60, 0.60, 0.60, 1.00)},
    "primaries":    {"name": "Primary antibodies",   "color": (0.50, 0.50, 0.50, 1.00)},
    # Blue rather than orange: orange sits too close to the euchromatin red
    # below to tell apart in a dense scene, and blue is otherwise unused.
    "secondaries":  {"name": "Secondary antibodies", "color": (0.35, 0.55, 0.90, 1.00)},
    "fluorophores": {"name": "Fluorophores",         "color": (1.00, 0.00, 1.00, 1.00)},
    "domains":      {"name": "Packing domains",      "color": (1.00, 0.90, 0.20, 1.00)},
}

# Green heterochromatin / red euchromatin, matching the convention already used
# for the localization track in hmt_v3/visualize.py (me3_border = [0, 0.7, 0],
# ac_border = [0.7, 0, 0]) -- H3K27me3 marks heterochromatin, H3K27ac marks
# euchromatin, so the same mark keeps the same color across both halves of the
# toolbox. Monomers read as muted context; dyes are brighter and more saturated
# versions of the same hue, so a dye is obviously distinct from the beads
# underneath it while still reading as belonging to that compartment.
_COMPARTMENT_COLORS = {
    "heterochromatin": (0.00, 0.70, 0.00, 1.00),  # green
    "euchromatin":     (0.70, 0.00, 0.00, 1.00),  # red
    "none":            (0.60, 0.60, 0.60, 1.00),
}
_DYE_COLORS = {
    "heterochromatin": (0.20, 1.00, 0.40, 1.00),  # bright spring green
    "euchromatin":     (1.00, 0.25, 0.45, 1.00),  # bright pink-red
    "none":            (1.00, 0.00, 1.00, 1.00),
}

_DOMAIN_TRANSPARENCY = 0.7


def resolve_fluorophore_channel(fluorophores_df, secondaries_df, primaries_df, monomers_df):
    """
    Walks the fluorophore -> secondary -> primary -> monomer foreign-key
    chain to recover which compartment each dye is actually reporting on.

    A fluorophores_df carries no compartment of its own -- it only knows its
    secondary_id -- so in a single-run labeling of a compartment-classified
    monomer cloud this is the only way to tell the two dye species apart.
    (In the two-channel pipeline workflow the channels come from separate
    run_pipeline calls and are already known, in which case pass
    `channel=` to build_particle_table instead of using this.)

    Args:
        fluorophores_df: output of labeling.place_fluorophores
        secondaries_df:  output of labeling.place_secondaries
        primaries_df:    output of labeling.place_primaries
        monomers_df:     the monomer table the labeling ran on, already
                         through compartments.classify_compartments

    Returns:
        np.ndarray of compartment strings ("heterochromatin",
        "euchromatin", or "none" where the chain could not be resolved),
        in fluorophores_df's row order

    Raises:
        ValueError: if monomers_df has no "compartment" column, or a
            required key column is missing from any table in the chain
    """
    if "compartment" not in monomers_df.columns:
        raise ValueError(
            "monomers_df has no 'compartment' column -- run "
            "compartments.classify_compartments first, or pass channel= explicitly"
        )
    for label, df, col in (
        ("fluorophores_df", fluorophores_df, "secondary_id"),
        ("secondaries_df", secondaries_df, "primary_id"),
        ("primaries_df", primaries_df, "monomer_id"),
    ):
        if col not in df.columns:
            raise ValueError(f"{label} is missing its '{col}' foreign key")

    resolved = (
        fluorophores_df[["secondary_id"]]
        .merge(secondaries_df[["secondary_id", "primary_id"]], on="secondary_id", how="left")
        .merge(primaries_df[["primary_id", "monomer_id"]], on="primary_id", how="left")
        .merge(monomers_df[["monomer_id", "compartment"]], on="monomer_id", how="left")
    )
    return resolved["compartment"].fillna("none").to_numpy()


def build_particle_table(monomers_df=None, primaries_df=None, secondaries_df=None,
                         fluorophores_df=None, domains_df=None, domain_z_nm=None,
                         channel=None, stride=1, layers=None, bead_radius_nm=4.9,
                         label_radius_nm=4.9, dye_radius_nm=4.9):
    """
    Flattens whichever srev tables are supplied into one tidy particle
    table, tagging each row with the species, channel, radius and scalar
    fields the exporter and the viewers all read from.

    Handles all three shapes of input:
      * raw SR-EV          -- monomers_df only, no compartment column:
                              everything is species MONO, coordination 0
      * het/eu labeled     -- monomers_df through
                              compartments.classify_compartments: monomers
                              become species HET / EUC and carry their CN
      * labeled + dyes     -- add the labeling tables. Each dye's channel is
                              taken from `channel` if given, otherwise
                              resolved per-dye through the FK chain (see
                              resolve_fluorophore_channel), producing
                              distinct DYE_HET / DYE_EUC species with
                              distinct colors

    For a two-channel acquisition (two pipeline.run_pipeline calls), build
    one table per channel with an explicit `channel=` and concatenate:

        het = build_particle_table(monomers_df=mono, fluorophores_df=f_het,
                                   channel="heterochromatin")
        eu  = build_particle_table(fluorophores_df=f_eu,
                                   channel="euchromatin")
        particles = pd.concat([het, eu], ignore_index=True)

    Args:
        monomers_df:     output of io.parse_config_txt, optionally already
                          through compartments.classify_compartments
        primaries_df, secondaries_df, fluorophores_df:
                          outputs of the corresponding labeling.place_* calls
        domains_df:       output of domains.identify_packing_domains --
                          radius_nm becomes the display radius, so domains
                          draw at their measured size
        domain_z_nm:      z to place domain centers at. REQUIRED with
                          domains_df, because identify_packing_domains works
                          on a single projected z-slab and returns x/y only;
                          pass the slab center used
        channel:          "heterochromatin" / "euchromatin" / None. Labels
                          this call's antibody and dye layers. None means
                          resolve each dye's channel from the FK chain when
                          the whole chain plus a compartment column is
                          available, and leave it unknown otherwise
        stride:           keep every Nth monomer only (1 = all). Applies to
                          monomers alone -- the label chain and domains are
                          orders of magnitude smaller and are never thinned,
                          so a strided export still shows every dye
        layers:           which layers to EMIT, e.g.
                          ("monomers", "fluorophores"). None emits every
                          layer that was supplied. Layers left out are
                          still used for their side effects -- notably,
                          passing primaries_df/secondaries_df while
                          emitting only monomers and fluorophores keeps the
                          dye-channel FK walk working (see
                          resolve_fluorophore_channel) while keeping the
                          antibodies out of the picture
        bead_radius_nm, label_radius_nm, dye_radius_nm:
                          display radii, all defaulting to r° = 4.9 nm (the
                          bead excluded-volume radius from Cangnano et al.
                          2024) so every point draws the same size. This
                          view is about WHERE things are, not how big they
                          are, and at true relative scale a dye is ~1/5 the
                          radius of a nucleosome and effectively invisible
                          next to it. For a volume-accurate picture instead,
                          pass roughly label_radius_nm=5.0 (IgG hydrodynamic
                          radius) and dye_radius_nm=0.7 (a small organic
                          dye). None of these affect the simulation -- they
                          are display only; placement is governed by the
                          linker lengths and exclusion_radius_nm in
                          labeling.LabelingConfig

    Returns:
        DataFrame with columns layer, species, compartment, x, y, z (nm),
        radius_nm, coordination_number, transparency. `compartment` is the
        per-row channel: a monomer's own compartment, and for the label
        layers whatever `channel` said or the FK walk resolved.

    Raises:
        ValueError: if domains_df is given without domain_z_nm, if a
            supplied table lacks its x/y/z [nm] columns, or if `channel`
            is not a recognized compartment name
    """
    if channel is not None and channel not in _COMPARTMENT_SHORT:
        raise ValueError(
            f"channel must be 'heterochromatin', 'euchromatin', or None, got {channel!r}"
        )

    frames = []

    def _add(layer, df, *, compartment, coordination, radius, z=None):
        n = len(df)
        compartment = np.asarray(compartment, dtype=object)
        species = np.array([
            _species_name(layer, c) for c in compartment
        ], dtype=object)
        frames.append(pd.DataFrame({
            "layer": layer,
            "species": species,
            "compartment": compartment,
            "x": df["x [nm]"].to_numpy(dtype=float),
            "y": df["y [nm]"].to_numpy(dtype=float),
            "z": (df["z [nm]"].to_numpy(dtype=float) if z is None
                  else np.full(n, float(z))),
            "radius_nm": np.broadcast_to(np.asarray(radius, dtype=float), (n,)).copy(),
            "coordination_number": np.asarray(coordination, dtype=np.int64),
            "transparency": _DOMAIN_TRANSPARENCY if layer == "domains" else 0.0,
        }))

    def _check(label, df, need_z=True):
        cols = _XYZ_COLS if need_z else _XYZ_COLS[:2]
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise ValueError(f"{label}: missing required column(s) {missing}")

    if monomers_df is not None and len(monomers_df):
        _check("monomers_df", monomers_df)
        mono = monomers_df.iloc[::stride] if stride > 1 else monomers_df
        compartment = (mono["compartment"].to_numpy() if "compartment" in mono.columns
                       else np.full(len(mono), "none", dtype=object))
        cn = (mono["coordination_number"].to_numpy() if "coordination_number" in mono.columns
              else np.zeros(len(mono), dtype=np.int64))
        _add("monomers", mono, compartment=compartment, coordination=cn, radius=bead_radius_nm)

    # The dye channel is what distinguishes the two fluorophore species, so it
    # is resolved once here and reused; antibodies inherit whatever `channel`
    # says, since walking the chain for them adds nothing a caller doing a
    # single-compartment run does not already know.
    dye_compartment = None
    if fluorophores_df is not None and len(fluorophores_df):
        if channel is not None:
            dye_compartment = np.full(len(fluorophores_df), channel, dtype=object)
        elif (secondaries_df is not None and primaries_df is not None
                and monomers_df is not None and "compartment" in monomers_df.columns):
            dye_compartment = resolve_fluorophore_channel(
                fluorophores_df, secondaries_df, primaries_df, monomers_df
            )
        else:
            dye_compartment = np.full(len(fluorophores_df), "none", dtype=object)

    for layer, df, radius in (
        ("primaries", primaries_df, label_radius_nm),
        ("secondaries", secondaries_df, label_radius_nm),
    ):
        if df is not None and len(df):
            _check(f"{layer}_df", df)
            comp = np.full(len(df), channel if channel is not None else "none", dtype=object)
            _add(layer, df, compartment=comp, coordination=np.zeros(len(df)), radius=radius)

    if fluorophores_df is not None and len(fluorophores_df):
        _check("fluorophores_df", fluorophores_df)
        _add("fluorophores", fluorophores_df, compartment=dye_compartment,
             coordination=np.zeros(len(fluorophores_df)), radius=dye_radius_nm)

    if domains_df is not None and len(domains_df):
        _check("domains_df", domains_df, need_z=False)
        if domain_z_nm is None:
            raise ValueError(
                "domain_z_nm is required with domains_df -- identify_packing_domains "
                "returns x/y only (it works on one projected z-slab), so the centers "
                "have no z of their own. Pass the slab center you used."
            )
        comp = np.full(len(domains_df), channel if channel is not None else "none", dtype=object)
        _add("domains", domains_df, compartment=comp,
             coordination=np.zeros(len(domains_df)),
             radius=domains_df["radius_nm"].to_numpy(dtype=float), z=domain_z_nm)

    if not frames:
        raise ValueError("no tables supplied -- nothing to export")

    particles = pd.concat(frames, ignore_index=True)

    if layers is not None:
        unknown = set(layers) - set(_LAYER_ORDER)
        if unknown:
            raise ValueError(
                f"unknown layer name(s) {sorted(unknown)} -- valid layers are {list(_LAYER_ORDER)}"
            )
        # Filtered only at the end, so layers left out have already done their
        # work: primaries/secondaries can resolve the dye channel above and
        # then be dropped from the output.
        particles = particles[particles["layer"].isin(set(layers))].reset_index(drop=True)
        if len(particles) == 0:
            raise ValueError(
                f"layers={tuple(layers)} selected nothing -- none of those layers "
                "were among the supplied tables"
            )

    return particles


def _species_name(layer, compartment):
    """Species (OVITO particle type) name for one row: base name, suffixed by compartment where known."""
    if layer == "monomers":
        return _COMPARTMENT_SHORT.get(compartment, _LAYER_SPECIES["monomers"])
    base = _LAYER_SPECIES[layer]
    if layer == "domains":
        return base
    short = _COMPARTMENT_SHORT.get(compartment)
    return f"{base}_{short}" if short else base


def particle_colors(particles, color_by="species", cn_range=(0, 12), colormap="coolwarm"):
    """
    Per-particle RGBA, computed here with matplotlib rather than delegated
    to each viewer's own colormap handling, so the OVITO export, the napari
    view and the inline matplotlib figure color a structure identically.

    Fluorophores are ALWAYS colored by their channel, in every mode -- the
    two dye species staying visually distinct is the point of separating
    them, so it is not something a color_by choice should be able to
    switch off.

    Args:
        particles: DataFrame from build_particle_table
        color_by:  "species" (default)     -- each particle type its own
                   color: monomers green (heterochromatin) / red
                   (euchromatin), antibodies grey and blue, domains yellow.
                   "compartment" is accepted as a synonym, since a
                   monomer's species IS its compartment.
                   "coordination_number"  -- monomers shaded across
                   cn_range instead; every other layer is unaffected
        cn_range:  (low, high) for the CN shading; 0-12 is the range
                   Cangnano et al. 2024 report
        colormap:  any matplotlib colormap name

    Returns:
        (N, 4) float array of RGBA in particles' row order
    """
    if color_by not in ("species", "compartment", "coordination_number"):
        raise ValueError(
            "color_by must be 'species', 'compartment', or 'coordination_number', "
            f"got {color_by!r}"
        )

    rgba = np.zeros((len(particles), 4), dtype=float)
    layers = particles["layer"].to_numpy()

    for layer, spec in _LAYER_DISPLAY.items():
        mask = layers == layer
        if mask.any():
            rgba[mask] = spec["color"]

    mono = layers == "monomers"
    if mono.any():
        if color_by == "coordination_number":
            lo, hi = cn_range
            values = particles.loc[mono, "coordination_number"].to_numpy(dtype=float)
            rgba[mono] = colormaps[colormap](Normalize(vmin=lo, vmax=hi)(values))
        else:
            # A monomer's species IS its compartment (HET / EUC), so coloring
            # by species means coloring by compartment -- without this, both
            # types would come out the same flat grey and "color by species"
            # would not distinguish anything.
            comp = particles.loc[mono, "compartment"].to_numpy()
            rgba[mono] = [_COMPARTMENT_COLORS.get(c, _COMPARTMENT_COLORS["none"]) for c in comp]

    dyes = layers == "fluorophores"
    if dyes.any():
        comp = particles.loc[dyes, "compartment"].to_numpy()
        rgba[dyes] = [_DYE_COLORS.get(c, _DYE_COLORS["none"]) for c in comp]

    return rgba


def write_extended_xyz(path, particles, color_by="species", cn_range=(0, 12),
                       colormap="coolwarm", include_color=True):
    """
    Writes a particle table as extended XYZ for OVITO. See the module
    docstring for the property list and what OVITO does with each.

    Args:
        path:      destination .xyz path
        particles: DataFrame from build_particle_table
        color_by, cn_range, colormap: see particle_colors -- these decide
                   the baked-in Color property, i.e. how the scene looks
                   the moment it opens. Ignored when include_color=False.
        include_color: whether to write the per-particle Color property.

                   TRUE (default) bakes colors in, so the file opens
                   already looking right -- but OVITO gives an explicit
                   Color property priority over its per-particle-type
                   colors, so the type colors in OVITO's type list will
                   have no visible effect. A Color Coding modifier still
                   works, because it runs later in the pipeline and
                   overwrites Color.

                   FALSE omits the column, leaving OVITO to color by
                   particle TYPE. Use this when you want to drive the
                   colors from OVITO's type list yourself (i.e. color by
                   species rather than by whatever color_by baked in).
                   Every numeric field is written either way, so Color
                   Coding on Coordination/compartment/layer is available
                   in both cases.

    Returns:
        the written path, as a string
    """
    compartment_code = np.array(
        [_COMPARTMENT_CODE.get(c, 0) for c in particles["compartment"].to_numpy()],
        dtype=np.int64,
    )

    columns = {
        "species": particles["species"].to_numpy(),
        "x": particles["x"].to_numpy(dtype=float),
        "y": particles["y"].to_numpy(dtype=float),
        "z": particles["z"].to_numpy(dtype=float),
    }
    if include_color:
        rgba = particle_colors(particles, color_by, cn_range, colormap)
        columns.update({"r": rgba[:, 0], "g": rgba[:, 1], "b": rgba[:, 2]})
    columns.update({
        "radius": particles["radius_nm"].to_numpy(dtype=float),
        "transparency": particles["transparency"].to_numpy(dtype=float),
        "coordination": particles["coordination_number"].to_numpy(dtype=np.int64),
        "compartment": compartment_code,
        "layer": np.array([_LAYER_CODE[l] for l in particles["layer"].to_numpy()], dtype=np.int64),
    })
    out = pd.DataFrame(columns)

    properties = (
        "Properties=species:S:1:pos:R:3:"
        + ("Color:R:3:" if include_color else "")
        + "Radius:R:1:Transparency:R:1:Coordination:I:1:compartment:I:1:layer:I:1"
    )
    comment = (
        f'{properties} units="nm" '
        f'compartment_codes="0=none,1=heterochromatin,2=euchromatin" '
        f'layer_codes="0=monomers,1=primaries,2=secondaries,3=fluorophores,4=domains"'
    )

    path = Path(path)
    with open(path, "w", newline="") as f:
        f.write(f"{len(out)}\n")
        f.write(comment + "\n")
        out.to_csv(f, sep=" ", header=False, index=False,
                   float_format="%.3f", lineterminator="\n")
    return str(path)


def export_for_ovito(output_dir, particles=None, basename="srev",
                     color_by="species", cn_range=(0, 12), colormap="coolwarm",
                     include_color=True, **table_kwargs):
    """
    One-call export: builds the particle table if one isn't supplied,
    writes the extended XYZ, and returns the path. Open the result with
    File > Load File in OVITO -- it arrives already colored and sized, with
    each population a separately selectable particle type.

    Args:
        output_dir:    directory to write into (created if missing)
        particles:     a prebuilt table from build_particle_table. Pass this
                       for a multi-channel scene (build one table per
                       channel and pd.concat them -- see
                       build_particle_table). If None, one is built from
                       **table_kwargs
        basename:      stem for the written file
        color_by, cn_range, colormap, include_color:
                       see write_extended_xyz. Note include_color=True
                       (the default) means OVITO shows THESE colors and
                       ignores its own per-type colors -- pass
                       include_color=False to color by particle type in
                       OVITO instead.
        **table_kwargs: forwarded to build_particle_table when `particles`
                       is None (monomers_df, primaries_df, ..., channel,
                       stride, radii)

    Returns:
        dict with "structure" (path), "n_particles", and "species" (the
        particle type names present, i.e. what OVITO's type list will show)
    """
    if particles is None:
        particles = build_particle_table(**table_kwargs)
    elif table_kwargs:
        raise ValueError(
            "pass either `particles` or the table arguments, not both "
            f"(got both, with table arguments {sorted(table_kwargs)})"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{basename}.xyz"
    write_extended_xyz(path, particles, color_by, cn_range, colormap, include_color)

    return {
        "structure": str(path),
        "n_particles": len(particles),
        "species": sorted(particles["species"].unique().tolist()),
    }


def structure_layers(particles, color_by="species", cn_range=(0, 12),
                     colormap="coolwarm"):
    """
    Splits a particle table into one point-layer description per species,
    ready to hand to napari. Pure data assembly -- imports no viewer, so it
    can be unit-tested without a display.

    Splitting by SPECIES rather than by layer is what gives the two dye
    populations their own napari layers, so they can be toggled and
    recolored independently, matching how they appear as separate types in
    OVITO.

    Point sizes are DIAMETERS IN NANOMETRES (2 * radius_nm), which is
    directly meaningful because napari sizes points in data coordinates and
    this package's coordinates are nm: a monomer really is drawn 9.8 nm
    across and a packing domain at its measured diameter.

    Args:
        particles:                     DataFrame from build_particle_table
        color_by, cn_range, colormap:  see particle_colors

    Returns:
        list of dicts, each with keys data ((N, 3) xyz array), face_color,
        size, name, opacity -- in draw order, monomers first
    """
    rgba = particle_colors(particles, color_by, cn_range, colormap)

    layers = []
    for layer in _LAYER_ORDER:
        in_layer = (particles["layer"] == layer).to_numpy()
        if not in_layer.any():
            continue
        for species in pd.unique(particles.loc[in_layer, "species"]):
            mask = in_layer & (particles["species"] == species).to_numpy()
            subset = particles.loc[mask]
            name = _LAYER_DISPLAY[layer]["name"]
            compartment = subset["compartment"].iloc[0]
            if compartment in _COMPARTMENT_SHORT:
                # Every split species needs a distinct layer name, monomers
                # included -- napari silently renames duplicates, which would
                # make the het/eu split unreadable in the layer list.
                name = f"{name} ({compartment})"
            layers.append({
                "data": subset[["x", "y", "z"]].to_numpy(dtype=float),
                "face_color": rgba[mask],
                "size": 2.0 * subset["radius_nm"].to_numpy(dtype=float),
                "name": name,
                # Domains are hollow context around the beads they contain, so
                # they are drawn faint enough to see through.
                "opacity": 1.0 - _DOMAIN_TRANSPARENCY if layer == "domains" else 1.0,
            })

    return layers


def view_napari(particles=None, color_by="species", cn_range=(0, 12),
                colormap="coolwarm", viewer=None, **table_kwargs):
    """
    Opens the structure in napari -- interactive 3D, driven from Python, no
    separate application. Recommended for looking at a structure while
    working; use export_for_ovito for large scenes and rendering.

    Follows the same conventions as hmt_v3/visualize.py's existing napari
    functions (3D display, out_of_slice_display, (2, 1, 0) dims order) so
    both halves of the toolbox behave the same way.

    Args:
        particles:      a prebuilt table from build_particle_table, or None
                        to build one from **table_kwargs
        color_by, cn_range, colormap: see particle_colors
        viewer:         an existing napari.Viewer to add layers to; a new
                        one is created if None
        **table_kwargs: forwarded to build_particle_table when `particles`
                        is None

    Returns:
        the napari.Viewer, so further layers or settings can be added

    Raises:
        ImportError: if napari isn't installed (it is a GUI dependency; the
            OVITO export path deliberately does not need it)
    """
    try:
        import napari
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise ImportError(
            "view_napari requires napari (already in environment.yml; "
            "`conda install -c conda-forge napari pyqt` if missing). For a "
            "viewer-free look, use plot_structure_3d or export_for_ovito."
        ) from exc

    if particles is None:
        particles = build_particle_table(**table_kwargs)
    elif table_kwargs:
        raise ValueError("pass either `particles` or the table arguments, not both")

    viewer = viewer if viewer is not None else napari.Viewer()
    for layer in structure_layers(particles, color_by, cn_range, colormap):
        viewer.add_points(
            data=layer["data"],
            face_color=layer["face_color"],
            size=5,
            name=layer["name"],
            opacity=0.35,
            border_width=0,
            border_color=layer["face_color"],
            out_of_slice_display=True,
        )

    viewer.dims.order = (2, 1, 0)
    viewer.dims.ndisplay = 3
    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = "nm"
    return viewer


def plot_structure_3d(particles=None, color_by="species", cn_range=(0, 12),
                      colormap="coolwarm", max_points=20000, point_size=2.0,
                      ax=None, title=None, **table_kwargs):
    """
    Static 3D scatter that renders INLINE in a notebook -- no separate
    window, no viewer install, and unlike napari it works headless. Best
    for a quick sanity look or a figure panel; use view_napari to actually
    explore a structure, or export_for_ovito to render one.

    matplotlib's 3D scatter slows to a crawl well before 10^5 points and
    cannot size markers in data units, so this thins the cloud to
    `max_points` (reporting that it did) and takes marker size in points.

    Args:
        particles:      a prebuilt table from build_particle_table, or None
                        to build one from **table_kwargs
        color_by, cn_range, colormap: see particle_colors
        max_points:     thinning cap applied to the assembled table (after
                        any `stride`); None disables it
        point_size:     matplotlib marker size (points^2), not nanometres
        ax:             existing 3D axes to draw into; one is created if None
        title:          axes title
        **table_kwargs: forwarded to build_particle_table when `particles`
                        is None

    Returns:
        (fig, ax)
    """
    import matplotlib.pyplot as plt  # lazy: keeps a backend off the import path

    if particles is None:
        particles = build_particle_table(**table_kwargs)
    elif table_kwargs:
        raise ValueError("pass either `particles` or the table arguments, not both")

    if max_points is not None and len(particles) > max_points:
        # integer stride, so a table only just over the cap rounds to 1 and would
        # be sliced (and announced) as a no-op -- skip it entirely in that case
        keep = len(particles) // max_points
        if keep > 1:
            particles = particles.iloc[::keep]
            print(f"plot_structure_3d: thinned to {len(particles)} points (every {keep}th) "
                  f"for matplotlib -- raise max_points, or use view_napari / "
                  f"export_for_ovito for the full cloud")

    rgba = particle_colors(particles, color_by, cn_range, colormap)

    if ax is None:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = ax.figure

    ax.scatter(particles["x"], particles["y"], particles["z"], c=rgba, s=point_size,
               linewidths=0, depthshade=False)

    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_zlabel("z [nm]")
    if title:
        ax.set_title(title)
    ax.set_box_aspect((1, 1, 1))
    return fig, ax
