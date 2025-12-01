#!/usr/bin/env python

"""
Build FRB host JSON files for the web frontend.

This script reads the generic host template:

    frb/data/Galaxies/frb_host_template.json

and, for each FRB with an existing host file

    FRBxxxxxx_host.json

plus the corresponding FRB file

    frb/data/FRBs/FRBxxxxxx.json

produces a web-ready version

    FRBxxxxxx_host_web.json

in the same FRB subdirectory.

Usage:

    python frb/builds/build_hosts_web.py FRBs --frb FRB20180916B
    python frb/builds/build_hosts_web.py FRBs              # build all
"""

import argparse
import copy
import json
from pathlib import Path

from pkg_resources import resource_filename


# --------------------------------------------------------------------
# Locate BASE_DIRs and template robustly
# --------------------------------------------------------------------
def find_galaxies_base_dir() -> Path:
    """
    Locate frb/data/Galaxies/frb_host_template.json.
    """
    here = Path(__file__).resolve()

    # 1) Look for a repo-style layout relative to this file
    for parent in [here] + list(here.parents):
        candidate_base = parent / "frb" / "data" / "Galaxies"
        candidate_template = candidate_base / "frb_host_template.json"
        if candidate_template.exists():
            return candidate_base

    # 2) Fall back to package data
    try:
        pkg_dir = Path(resource_filename("frb", "data/Galaxies")).resolve()
        if (pkg_dir / "frb_host_template.json").exists():
            return pkg_dir
    except Exception:
        pass

    raise FileNotFoundError(
        "Could not find frb_host_template.json either in a repo-style "
        "checkout or in the installed 'frb' package."
    )


def find_frbs_base_dir() -> Path:
    """
    Locate frb/data/FRBs.
    """
    here = Path(__file__).resolve()

    # 1) Repo-style layout
    for parent in [here] + list(here.parents):
        candidate_base = parent / "frb" / "data" / "FRBs"
        if candidate_base.exists():
            return candidate_base

    # 2) Installed package data
    try:
        pkg_dir = Path(resource_filename("frb", "data/FRBs")).resolve()
        if pkg_dir.exists():
            return pkg_dir
    except Exception:
        pass

    raise FileNotFoundError(
        "Could not find frb/data/FRBs either in a repo-style "
        "checkout or in the installed 'frb' package."
    )


# --------------------------------------------------------------------
# JSON helpers
# --------------------------------------------------------------------
def load_json(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def save_json(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    # Don't sort keys so metadata/data grouping stays as we construct it.
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


# --------------------------------------------------------------------
# Utility: recursive key search
# --------------------------------------------------------------------
def find_key_recursive(data, key):
    """
    Recursively search for `key` anywhere inside `data` (dicts/lists).

    Returns the first value found, or None if not present.
    """
    if isinstance(data, dict):
        if key in data:
            return data[key]
        for v in data.values():
            found = find_key_recursive(v, key)
            if found is not None:
                return found
    elif isinstance(data, list):
        for item in data:
            found = find_key_recursive(item, key)
            if found is not None:
                return found
    return None


def get_scalar_from_source(source_data, key):
    """
    Find `key` in source_data (recursively) and return a scalar:

    - If the value is a dict with a 'value' field, return that.
    - Otherwise return the raw value.
    - If not found, return None.
    """
    if not source_data:
        return None
    raw_val = find_key_recursive(source_data, key)
    if raw_val is None:
        return None
    if isinstance(raw_val, dict) and "value" in raw_val:
        return raw_val["value"]
    return raw_val

def override_in_tree(tree, key, new_val):
    """
    Recursively walk `tree` (dicts/lists) and replace any occurrence
    of `key` with `new_val`. Does not create new keys, only overrides
    existing ones.
    """
    if isinstance(tree, dict):
        if key in tree:
            tree[key] = new_val
        for v in tree.values():
            override_in_tree(v, key, new_val)
    elif isinstance(tree, list):
        for item in tree:
            override_in_tree(item, key, new_val)


# --------------------------------------------------------------------
# Data extraction based on template (metadata)
# --------------------------------------------------------------------
def extract_data_by_template(template, source_data):
    """
    Build a 'data' dict from a metadata template and a source_data dict
    (merged host + FRB JSON).

    Rules:
      * Only include fields that have corresponding metadata entries
        (i.e., keys present in the template tree).
      * 'metadata' controls the shape (grouping) of the output, but we
        do NOT try to mirror its internal keys (description, unit, etc.).
      * For a leaf metadata block (e.g., template['DM'] is a dict whose
        values are not all dicts), we treat the key ('DM') as the data
        key and grab the value of that key from source_data, searching
        recursively if needed.
      * If the found value is a dict with a 'value' key (e.g. {"unit":..., "value":...}),
        we store only the numeric part (val['value']) in the data block.
    """
    if not isinstance(template, dict):
        return {}

    result = {}

    for key, t_val in template.items():
        # Case 1: metadata value is a dict
        if isinstance(t_val, dict):
            # Heuristic: if ANY of the values in t_val are non-dicts,
            # then this is a *leaf metadata block* describing a single
            # quantity (e.g., DM, ra_frb, etc.), not a pure grouping node.
            any_non_dict = any(not isinstance(v, dict) for v in t_val.values())

            if any_non_dict:
                # Leaf metadata: look for this key anywhere in source_data
                raw_val = find_key_recursive(source_data, key)
                if raw_val is not None:
                    # If the value is a dict with a 'value' field, keep only that
                    if isinstance(raw_val, dict) and "value" in raw_val:
                        result[key] = raw_val["value"]
                    else:
                        result[key] = raw_val
            else:
                # Grouping node: recurse into children, but still search
                # over the full source_data, not restricted by key path.
                sub_data = extract_data_by_template(t_val, source_data)
                if sub_data:
                    result[key] = sub_data

        # Case 2: metadata value is NOT a dict (string, number, etc.)
        else:
            # Treat this as metadata for a single quantity with name=key.
            raw_val = find_key_recursive(source_data, key)
            if raw_val is not None:
                if isinstance(raw_val, dict) and "value" in raw_val:
                    result[key] = raw_val["value"]
                else:
                    result[key] = raw_val

    return result



# --------------------------------------------------------------------
# Core build logic
# --------------------------------------------------------------------
def build_single(
    frb_name: str,
    galaxies_base: Path,
    frbs_base: Path,
    template: dict,
    overwrite: bool = False,
):
    """
    Build FRB<...>_host_web.json for a single FRB.

    Parameters
    ----------
    frb_name : str
        Name like 'FRB20180916B'.
    galaxies_base : Path
        Path to frb/data/Galaxies.
    frbs_base : Path
        Path to frb/data/FRBs.
    template : dict
        Parsed template JSON (treated entirely as metadata).
    overwrite : bool
        If True, overwrite existing *_host_web.json.
    """
    # Galaxy/host JSON: e.g. FRB20180916B_host.json
    host_matches = list(galaxies_base.rglob(f"{frb_name}_host.json"))
    if len(host_matches) == 0:
        print(f"[WARN] No *_host.json found for {frb_name}")
        return

    host_path = host_matches[0]

    # FRB JSON: e.g. frb/data/FRBs/FRB20180916B.json
    frb_json_path = frbs_base / f"{frb_name}.json"

    print(f"[INFO] Building web host JSON for {frb_name}")
    print(f"[INFO]   Host JSON: {host_path}")
    if frb_json_path.exists():
        print(f"[INFO]   FRB JSON:  {frb_json_path}")
    else:
        print(f"[WARN]   No FRB JSON found at {frb_json_path}; using host JSON only")

    # Load host data
    host_data = load_json(host_path)

    # Load FRB data if available
    frb_data = {}
    if frb_json_path.exists():
        frb_data = load_json(frb_json_path)

    # Merge data sources:
    # start with host_data, then overlay frb_data so FRB entries
    # (e.g., DM, ra_frb) are available and can override if needed.
    merged_data = {}
    merged_data.update(host_data)
    merged_data.update(frb_data)

    # Metadata: just a copy of the template (we do NOT overwrite it with values)
    metadata = copy.deepcopy(template)

    # Data: only fields with corresponding metadata, using merged_data.
    # If something like DM is nested, find_key_recursive() will locate it.
    data = extract_data_by_template(template, merged_data)

    # ------------------------------------------------------------------
    # Enforce specific sources for coordinates:
    #   - 'ra' and 'dec' from host_data (galaxy/host file)
    #   - 'ra_frb' and 'dec_frb' from frb_data (FRB file)
    # ------------------------------------------------------------------

    # Host (galaxy) coordinates
    ra_host = get_scalar_from_source(host_data, "ra")
    dec_host = get_scalar_from_source(host_data, "dec")
    if ra_host is not None:
        data["ra"] = ra_host
    if dec_host is not None:
        data["dec"] = dec_host

    # FRB coordinates: primary keys 'ra_frb' / 'dec_frb';
    # if for some reason the FRB JSON only has 'ra'/'dec', fall back to those.
    ra_frb = get_scalar_from_source(frb_data, "ra_frb")
    if ra_frb is None:
        ra_frb = get_scalar_from_source(frb_data, "ra")
    dec_frb = get_scalar_from_source(frb_data, "dec_frb")
    if dec_frb is None:
        dec_frb = get_scalar_from_source(frb_data, "dec")

    if ra_frb is not None:
        data["ra_frb"] = ra_frb
    if dec_frb is not None:
        data["dec_frb"] = dec_frb

    # ------------------------------------------------------------------
    # eellipse: 'a' and 'b' from FRB JSON, convert deg -> arcsec
    # ------------------------------------------------------------------
    eellipse_src = frb_data.get("eellipse", {})
    if isinstance(eellipse_src, dict):
        # Values in FRB JSON are assumed to be in degrees
        a_deg = get_scalar_from_source(eellipse_src, "a")
        b_deg = get_scalar_from_source(eellipse_src, "b")

        # Convert to arcseconds
        a_arcsec = a_deg * 3600.0 if a_deg is not None else None
        b_arcsec = b_deg * 3600.0 if b_deg is not None else None

        if a_arcsec is not None or b_arcsec is not None:
            # Ensure we have an 'eellipse' dict in the data block
            if "eellipse" not in data or not isinstance(data["eellipse"], dict):
                data["eellipse"] = {}

            if a_arcsec is not None:
                data["eellipse"]["a"] = a_arcsec
            if b_arcsec is not None:
                data["eellipse"]["b"] = b_arcsec

        # ------------------------------------------------------------------
    # r-band magnitude:
    #   1) Prefer DECaL_r / DECaL_r_err
    #   2) Otherwise Pan-STARRS_r / Pan-STARRS_r_err
    #   r_mag_ref = "DECaL" or "Pan-STARRS" accordingly.
    # ------------------------------------------------------------------
    r_mag = None
    r_mag_err = None
    r_mag_ref = None

    # First try DECaL
    decal_r     = get_scalar_from_source(host_data, "DECaL_r")
    decal_r_err = get_scalar_from_source(host_data, "DECaL_r_err")

    if decal_r is not None:
        r_mag = decal_r
        r_mag_err = decal_r_err  # may be None
        r_mag_ref = "DECaL"
    else:
        # Fallback to Pan-STARRS
        ps_r     = get_scalar_from_source(host_data, "Pan-STARRS_r")
        ps_r_err = get_scalar_from_source(host_data, "Pan-STARRS_r_err")

        if ps_r is not None:
            r_mag = ps_r
            r_mag_err = ps_r_err
            r_mag_ref = "Pan-STARRS"

    # Write into data block if we found something
    if r_mag is not None:
        data["r_mag"] = r_mag
        if r_mag_err is not None:
            data["r_mag_err"] = r_mag_err
        if r_mag_ref is not None:
            data["r_mag_ref"] = r_mag_ref


    out_path = host_path.parent / f"{frb_name}_host_web.json"

    if out_path.exists() and not overwrite:
        print(f"[INFO]   -> Skipping existing file {out_path}")
        return

    combined = {
        "metadata": metadata,
        "data": data,
    }

    save_json(out_path, combined)

    print(f"[INFO]   -> Wrote {out_path}")


def build_all(galaxies_base: Path, frbs_base: Path, template: dict, overwrite: bool = False):
    """
    Build *_host_web.json for all FRBs with *_host.json under galaxies_base.
    """
    host_files = sorted(galaxies_base.rglob("*_host.json"))

    if not host_files:
        print(f"[WARN] No *_host.json files found under {galaxies_base}")
        return

    for host_path in host_files:
        if host_path.name == "frb_host_template.json":
            continue

        frb_name = host_path.stem.replace("_host", "")
        build_single(frb_name, galaxies_base, frbs_base, template, overwrite=overwrite)


# --------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------
def parse_args(options=None):
    """
    Parse command-line arguments.

    We mirror the style of frb_build_web by requiring a first positional
    argument specifying what to build (here: 'FRBs').
    """
    parser = argparse.ArgumentParser(
        description=(
            "Build FRB *_host_web.json files with separate 'metadata' and 'data' "
            "blocks, where 'data' comes from the merged host+FRB JSON and only "
            "contains fields that have corresponding metadata."
        )
    )

    parser.add_argument(
        "what",
        help="What to build (must be 'FRBs' to mirror frb_build_web).",
        choices=["FRBs"],
    )

    parser.add_argument(
        "--frb",
        type=str,
        default=None,
        help="Specific FRB name (e.g. FRB20180916B). If omitted, build all.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing *_host_web.json files.",
    )

    return parser.parse_args(options)


def main(options=None):
    """
    Entry point for script / console usage.
    """
    args = parse_args(options)

    galaxies_base = find_galaxies_base_dir()
    frbs_base = find_frbs_base_dir()
    template_path = galaxies_base / "frb_host_template.json"

    print(f"[INFO] Using Galaxies directory: {galaxies_base}")
    print(f"[INFO] Using FRBs directory:     {frbs_base}")
    print(f"[INFO] Using template:           {template_path}")

    if not template_path.exists():
        raise FileNotFoundError(f"Template JSON not found at {template_path}")

    template = load_json(template_path)

    if args.frb:
        build_single(args.frb, galaxies_base, frbs_base, template, overwrite=args.overwrite)
    else:
        build_all(galaxies_base, frbs_base, template, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
