"""
Parser for externally-generated SR-EV chromatin polymer configurations.

An SR-EV run (a polymer-physics simulation of a chromatin fiber, generated
outside this repository) is handed off as a single config-N.txt file: a
LAMMPS-style dump with a 2-line header followed by one line per monomer.
This module turns that file into the monomer coordinate table the rest of
the srev package (labeling, photophysics, rendering) builds on.
"""
import gzip
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

_TIMESTEP_RE = re.compile(r"Timestep:\s*(-?\d+)")
_SREV_PARAM_RE = re.compile(r"(?:^|-)([A-Za-z]+)=([0-9.]+)")


def parse_config_txt(path, length_unit_to_nm=1.0):
    """
    Parses a config-N.txt SR-EV dump into a monomer coordinate table.

    File format (confirmed against hmt_v1/Simulation/config-1.txt):
        line 1: declared point count, e.g. "186742"
        line 2: a free-text line embedding the LAMMPS timestep, e.g.
                " Atoms. Timestep: 1300"
        remaining lines: one monomer per line, whitespace-separated. The
                LAST 3 tokens on each line are taken as x, y, z in the
                polymer simulation's native length units; any tokens
                BEFORE those are kept as metadata columns (type_0, type_1,
                ... in left-to-right order) rather than assumed to mean
                anything specific. The one sample file inspected carries a
                single leading column that is a constant "1" for every
                row (almost certainly a LAMMPS atom `type`, uniform
                because that run has no compartment/species labeling
                embedded yet) — this isn't required to hold for every
                config-N.txt, so the column count is read from the file,
                not hard-coded.

    Coordinates are in the polymer simulation's native units on input;
    `length_unit_to_nm` rescales them to nanometers. For the SR-EV model
    (Cangnano et al. 2024, eLife, "Local Volume Concentration, Packing
    Domains and Scaling Properties of Chromatin"), coordinates are
    generated in units of Umin, the minimum bond length, which the paper
    sets to 10 nm — i.e. length_unit_to_nm=10 for those runs. Confirmed
    empirically against hmt_v1/Simulation/config-1.txt: under a scale of
    10, the resulting nearest-neighbor spacing comes out to ~9.8-10.2 nm,
    matching the paper's stated 9.8 nm bead diameter almost exactly. The
    legacy MATLAB pipeline instead hardcoded `xpos*4.9` — 4.9 nm is
    r° = 0.49*Umin, the bead's EXCLUDED-VOLUME RADIUS in that same paper,
    not the coordinate length unit; the two constants coincide in the
    "SRRW-MARCELO-...-r=0.49-ru=10" folder-naming convention (see
    parse_srev_folder_name below) but mean different things. Using 4.9 as
    the coordinate scale understates every spatial measurement by roughly
    a factor of 2. This is a per-run parameter regardless — never assumed
    by this function, always supplied by the caller.

    Transparently reads gzip-compressed files: a path ending in ".gz" is
    opened with gzip.open in text mode instead of the builtin open, with
    no other change to the parsing logic. This is the expected case for
    a full SR-EV batch -- raw config-N.txt files compress roughly 4-6x
    (numeric text), so a batch that would be tens to a hundred-plus GB
    uncompressed fits in a few GB as .txt.gz, without ever needing an
    uncompressed copy on disk.

    Args:
        path:               path to a config-N.txt or config-N.txt.gz file
        length_unit_to_nm:  multiplier converting native coordinate units
                             to nanometers (10.0 for SR-EV runs following
                             the Cangnano et al. 2024 convention — i.e.
                             the run's "ru" value, NOT "r")

    Returns:
        (monomers_df, metadata):
            monomers_df: DataFrame with columns
                monomer_id, x [nm], y [nm], z [nm], x_raw, y_raw, z_raw,
                plus any leading columns as type_0, type_1, ... (cast to
                int64 where every value is integral, else left as parsed)
            metadata: dict with keys
                n_points_declared, n_points_parsed, timestep
                (None if "Timestep: <int>" wasn't found in line 2),
                header_line2 (raw text), source_path, length_unit_to_nm

    Raises:
        ValueError: if line 1 isn't an integer, if a data line has fewer
            than 3 columns, or if the parsed row count doesn't match the
            declared point count (a truncated file or a misidentified
            header would otherwise fail silently downstream instead).
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        line1 = f.readline()
        line2 = f.readline()

        try:
            n_points_declared = int(line1.strip())
        except ValueError as exc:
            raise ValueError(
                f"{path}: expected an integer point count on line 1, got {line1!r}"
            ) from exc

        timestep_match = _TIMESTEP_RE.search(line2)
        timestep = int(timestep_match.group(1)) if timestep_match else None

        data = pd.read_csv(f, sep=r"\s+", header=None)

    if len(data) != n_points_declared:
        raise ValueError(
            f"{path}: header declares {n_points_declared} points but "
            f"{len(data)} data lines were parsed"
        )

    n_cols = data.shape[1]
    if n_cols < 3:
        raise ValueError(
            f"{path}: each data line needs at least 3 columns (x y z), found {n_cols}"
        )

    coords_raw = data.iloc[:, -3:].to_numpy(dtype=float)
    lead_cols = data.iloc[:, : n_cols - 3]

    monomers_df = pd.DataFrame({
        "monomer_id": np.arange(len(data), dtype=np.int64),
        "x [nm]": coords_raw[:, 0] * length_unit_to_nm,
        "y [nm]": coords_raw[:, 1] * length_unit_to_nm,
        "z [nm]": coords_raw[:, 2] * length_unit_to_nm,
        "x_raw": coords_raw[:, 0],
        "y_raw": coords_raw[:, 1],
        "z_raw": coords_raw[:, 2],
    })

    for i, col in enumerate(lead_cols.columns):
        values = lead_cols[col]
        try:
            as_int = values.astype(np.int64)
        except (ValueError, TypeError):
            monomers_df[f"type_{i}"] = values
        else:
            monomers_df[f"type_{i}"] = as_int if np.all(as_int == values) else values

    metadata = {
        "n_points_declared": n_points_declared,
        "n_points_parsed": len(monomers_df),
        "timestep": timestep,
        "header_line2": line2.rstrip("\n"),
        "source_path": str(path),
        "length_unit_to_nm": length_unit_to_nm,
    }

    return monomers_df, metadata


def parse_srev_folder_name(name):
    """
    Parses the key=value parameters out of an SR-EV run's folder name, e.g.
    "SRRW-MARCELO-L=130-N=186741-phi=0.08-alpha=1.2-bmax=20-r=0.49-ru=10".

    Known fields, per Cangnano et al. 2024 (eLife, "Local Volume
    Concentration, Packing Domains and Scaling Properties of Chromatin")
    and this naming convention:
        L:     present in every run inspected so far as a constant 130
               regardless of N/phi/alpha, so it is NOT a per-run domain
               or cluster count to feed into compartment classification
               (an earlier version of this docstring assumed it was --
               that was wrong; see compartments.py, which uses a
               per-monomer coordination number instead and needs no
               domain count at all). Its exact meaning in the generator
               is not established here; returned as int in case a future
               run varies it meaningfully.
        N:     monomer count the run was generated with (int) -- matches
               phi = N*(r_deg_nm/Rc_nm)^3 with r_deg_nm = 4.9 nm and
               Rc_nm = 650 nm fixed across all runs in the paper (verified:
               186741*(4.9/650)**3 = 0.080). Cross-check against a parsed
               config-N.txt's own declared point count; a small off-by-one
               difference against the file's header is possible and not
               necessarily an error.
        phi:   volume fraction / packing density (float)
        alpha: polymer scaling exponent (float) -- matches the legacy
               alph1pX filename suffix, e.g. alpha=1.2 -> "alph1p2"
        r:     the bead EXCLUDED-VOLUME RADIUS r° in units of ru (float)
               -- r=0.49 means r° = 0.49*ru nm (4.9 nm when ru=10), the
               non-overlap radius between beads. This is NOT a coordinate
               scale factor.
        ru:    Umin, the minimum bond length / coordinate length unit, in
               nm (float) -- THIS is length_unit_to_nm for
               parse_config_txt (10.0 in every run inspected so far), not
               r and not r*ru. An earlier version of this function
               returned "length_unit_to_nm": r*ru, which silently
               reproduced the legacy MATLAB pipeline's scaling bug
               (see parse_config_txt) -- fixed below.
        bmax:  present in the naming convention; meaning not established
               from available documentation, returned as a float but not
               otherwise interpreted by this module

    Any other "key=value" segment present in the name is still returned
    (as a float, or as int for a name matched case-insensitively to "l"
    or "n") rather than silently dropped, so a run with additional
    parameters doesn't lose them here.

    Args:
        name: folder name, or a path (only the final path component is
              used), containing "-key=value" segments

    Returns:
        dict of {key: int|float}. Includes "length_unit_to_nm": ru
        (copied, not r*ru) whenever ru is present in the name.
    """
    base = os.path.basename(str(name).rstrip("/\\"))
    params = {}
    for key, value in _SREV_PARAM_RE.findall(base):
        key_lower = key.lower()
        if key_lower == "l":
            params["L"] = int(float(value))
        elif key_lower == "n":
            params["N"] = int(float(value))
        else:
            params[key_lower] = float(value)

    if "ru" in params:
        params["length_unit_to_nm"] = params["ru"]

    return params


def parse_copy_parameters_txt(path):
    """
    Parses a copy-parameters.txt file sitting alongside config-N.txt in an
    SRRW-* run folder: a plain list of numbers, one per line, in a fixed
    order, with no keys or header at all.

    Line order (1-indexed, as written in the file): N, L, alpha, r, bmax,
    param_6, ru, param_8. Confirmed against six example files spanning
    N in {186741, 280112, 373483} and alpha in {1.1, 1.15, 1.2}: N and
    alpha are the only values that ever change; L (130), r (0.49), bmax
    (20), and ru (10) are identical in every file and match the
    corresponding SRRW-* folder-name parameters (see
    parse_srev_folder_name) exactly. param_6 and param_8 are NOT part of
    the folder-naming convention and are also constant (0.1 and 150.0
    respectively) across every file checked -- likely fixed generator
    constants (param_8 = 150.0 is plausibly the "local cutoff Umax", the
    maximum bond length Cangnano et al. 2024's Appendix 3 describes only
    symbolically and never gives a numeric value for) but that mapping
    is a guess, not confirmed, so both are returned under generic names
    rather than a guessed one that could be silently wrong downstream.

    Because this file is a direct, machine-written record of the
    generator's own inputs rather than a human/script-composed folder
    label, prefer it over parse_srev_folder_name when both are available
    -- see resolve_srev_run_params, which does exactly that and
    cross-checks the two against each other.

    Args:
        path: path to a copy-parameters.txt file

    Returns:
        dict: N (int), L (int), alpha, r, bmax, ru, param_6, param_8
        (all float except N/L), and length_unit_to_nm (= ru)

    Raises:
        ValueError: if the file has fewer than 8 non-blank numeric lines
    """
    with open(path, "r") as f:
        values = [float(line.strip()) for line in f if line.strip()]

    if len(values) < 8:
        raise ValueError(f"{path}: expected 8 numeric lines, found {len(values)}")

    n, l, alpha, r, bmax, param_6, ru, param_8 = values[:8]

    return {
        "N": int(round(n)),
        "L": int(round(l)),
        "alpha": alpha,
        "r": r,
        "bmax": bmax,
        "ru": ru,
        "param_6": param_6,
        "param_8": param_8,
        "length_unit_to_nm": ru,
    }


def resolve_srev_run_params(config_path, rtol=1e-3):
    """
    Resolves an SRRW-* run's parameters for a given config-N.txt, preferring
    a sibling copy-parameters.txt (a direct machine-written record of the
    generator's inputs) over the folder name (a human/script-composed
    label that could in principle be mistyped, truncated, or renamed) --
    falling back to the folder name alone if no copy-parameters.txt exists
    next to this config file.

    When both are available, cross-checks every field the two sources
    share (N, L, alpha, r, ru; bmax is checked too if the folder name has
    it) and raises if any disagree by more than `rtol` -- silently
    trusting one source over a real mismatch would be worse than erroring.

    Args:
        config_path: path to a config-N.txt file (copy-parameters.txt is
                     looked for in the same directory)
        rtol:        relative tolerance for the cross-check

    Returns:
        dict of run parameters (see parse_copy_parameters_txt if a
        copy-parameters.txt was found, else parse_srev_folder_name's
        output for this config's parent folder), plus "source": either
        "copy-parameters.txt" or "folder-name"

    Raises:
        ValueError: if both sources are available but disagree beyond
            rtol on a shared field
    """
    config_path = Path(config_path)
    folder_params = parse_srev_folder_name(config_path.parent.name)

    params_path = config_path.parent / "copy-parameters.txt"
    if not params_path.exists():
        folder_params["source"] = "folder-name"
        return folder_params

    file_params = parse_copy_parameters_txt(params_path)
    for key in ("N", "L", "alpha", "r", "ru", "bmax"):
        if key in folder_params and key in file_params:
            a, b = folder_params[key], file_params[key]
            if abs(a - b) > rtol * max(abs(a), abs(b), 1e-9):
                raise ValueError(
                    f"{config_path}: folder name says {key}={a} but "
                    f"copy-parameters.txt says {key}={b}"
                )

    file_params["source"] = "copy-parameters.txt"
    return file_params
