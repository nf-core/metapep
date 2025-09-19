#!/usr/bin/env python3
"""
Prepare input files for MHC binding prediction tools.

This script:
- parses and normalizes allele names,
- validates peptide sequences,
- filters peptides by MHC class–specific length windows,
- optionally filters alleles against a supported-alleles JSON,
- and writes tool-specific input files in the required formats.

Author: Jonas Scheid (adapted)
"""

import argparse
import json
import logging
from enum import Enum
from pathlib import Path
import re
import pandas as pd
import mhcgnomes

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# --------------------------------------------------------------------------- #
# Tool-specific global length defaults
# These are conservative defaults used by additional per-tool filters below.
# The primary class-specific filtering is driven by CLI arguments (Class I/II).
# --------------------------------------------------------------------------- #
class MinLength(Enum):
    MHCFLURRY = 5
    MHCNUGGETS = 5
    NETMHCPAN = 8
    NETMHCIIPAN = 9

class MaxLength(Enum):
    MHCFLURRY = 15
    MHCNUGGETS_CLASSI = 15
    MHCNUGGETS_CLASSII = 30
    NETMHCPAN = 14
    NETMHCIIPAN = 50

# Max count guardrails (extend as needed)
class MaxNumberOfAlleles(Enum):
    NETMHCPAN = 50


# --------------------------- Utility helpers --------------------------- #
def parse_int_or_none(x: "str | None") -> "int | None":
    """Return int(x) or None if x is empty/None-ish."""
    if x is None:
        return None
    s = str(x).strip().lower()
    if s in {"", "none", "null", "na"}:
        return None
    return int(x)

def normalize_tools(tools_val) -> list[str]:
    """
    Normalize --tools input into a non-empty list of tool names.
    Accepts comma-separated string or a list; drops empty items.
    """
    if isinstance(tools_val, list):
        return [t for t in tools_val if str(t).strip()]
    return [t.strip() for t in str(tools_val or "").split(",") if t.strip()]

def parse_supported_json(path_str: "str | None") -> dict:
    """
    Robustly read a supported-alleles JSON mapping (tool -> [alleles]).
    Returns {} (empty mapping) if the file is missing, empty, malformed, or not a dict.
    """
    if not path_str:
        logging.warning("No --supported-alleles-json provided; continuing with empty mapping.")
        return {}
    p = Path(path_str)
    if not p.exists():
        logging.warning(f"supported_alleles_json not found: {p} (continuing with empty mapping).")
        return {}
    if p.stat().st_size == 0:
        logging.warning(f"supported_alleles_json is empty: {p} (continuing with empty mapping).")
        return {}
    try:
        with p.open() as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            logging.warning(f"supported_alleles_json is not a dict: {p} (continuing with empty mapping).")
            return {}
        return data
    except Exception as e:
        logging.warning(f"Could not parse supported_alleles_json ({p}): {e} (continuing with empty mapping).")
        return {}

# -------------------- Allele formatting helpers -------------------- #
def to_mouse_star_format(a: str) -> str:
    """
    Normalize common mouse allele spellings to star format, e.g.:
      H2-Db  -> H2-D*b
      H2-Kb  -> H2-K*b
      H2-Ld  -> H2-L*d
    Non-mouse alleles are returned unchanged.
    """
    s = re.sub(r'^H-?2-', 'H2-', a)  # normalize H-2- -> H2-
    return re.sub(r'^(H2-)([DKL])([a-z])$', r'\1\2*\3', s)

def to_human_star_format(a: str) -> str:
    """
    Normalize human allele names to 'HLA-<LOCUS>*<XX>:<YY(:ZZ)?>' style.

    Examples (→ normalized):
      A0201, A*0201, HLA-A0201, A*02:01           → HLA-A*02:01
      DRB1*1501                                   → HLA-DRB1*15:01
      DQB1*0302                                   → HLA-DQB1*03:02

    Already-correct strings are returned unchanged.
    """
    s = a.strip()

    # already correct
    if re.match(r'^HLA-[A-Z0-9]+[*][0-9]{2}:[0-9]{2,3}(?::[0-9]{2,3})*$', s):
        return s

    # Class I compact numeric (A/B/C)
    m = re.match(r'^(?:HLA-)?([ABC])\*?([0-9]{2})([0-9]{2,3})(?::([0-9]{2,3}))?$',
                 s, flags=re.IGNORECASE)
    if m:
        locus = m.group(1).upper(); g1, g2 = m.group(2), m.group(3)
        rest = f":{m.group(4)}" if m.group(4) else ""
        return f"HLA-{locus}*{g1}:{g2}{rest}"

    # Class I with colon but missing 'HLA-'
    m = re.match(r'^([ABC])\*([0-9]{2}):([0-9]{2,3})(?::([0-9]{2,3}))?$',
                 s, flags=re.IGNORECASE)
    if m:
        locus = m.group(1).upper(); g1, g2 = m.group(2), m.group(3)
        rest = f":{m.group(4)}" if m.group(4) else ""
        return f"HLA-{locus}*{g1}:{g2}{rest}"

    # Class II compact (DPA1/DPB1/DQA1/DQB1/DRB1–DRB5)
    m = re.match(r'^(?:HLA-)?(DP[AB]1|DQ[AB]1|DRB[1-5])\*?([0-9]{2})([0-9]{2,3})(?::([0-9]{2,3}))?$',
                 s, flags=re.IGNORECASE)
    if m:
        locus = m.group(1).upper(); g1, g2 = m.group(2), m.group(3)
        rest = f":{m.group(4)}" if m.group(4) else ""
        return f"HLA-{locus}*{g1}:{g2}{rest}"

    # Class II with colon but missing 'HLA-'
    m = re.match(r'^(DP[AB]1|DQ[AB]1|DRB[1-5])\*([0-9]{2}):([0-9]{2,3})(?::([0-9]{2,3}))?$',
                 s, flags=re.IGNORECASE)
    if m:
        locus = m.group(1).upper(); g1, g2 = m.group(2), m.group(3)
        rest = f":{m.group(4)}" if m.group(4) else ""
        return f"HLA-{locus}*{g1}:{g2}{rest}"

    # Minimal forms like A*01 → just add HLA- prefix
    m = re.match(r'^([ABC])\*([0-9]{2})$', s, flags=re.IGNORECASE)
    if m:
        return f"HLA-{m.group(1).upper()}*{m.group(2)}"

    return s


def parse_alleles(allele_str: "str | None") -> list[str]:
    """
    Parse alleles from:
      - a file path (.txt/.csv/.tsv) with one allele per line, OR
      - a string like "A*02:01;B*07:02" (semicolons or whitespace separated).

    Steps:
      - drop empty/None-ish tokens,
      - normalize via mhcgnomes (robust parsing),
      - convert to preferred human/mouse star formats,
      - de-duplicate while preserving order.
    """
    raw: list[str] = []
    if allele_str and str(allele_str).lower() not in {"null", "none"}:
        s = str(allele_str).strip()
        if s.endswith((".txt", ".csv", ".tsv")) and Path(s).exists():
            with open(s, "r") as fh:
                raw = [ln.strip() for ln in fh]
        else:
            raw = [x.strip() for x in re.split(r'[;\s]+', str(allele_str or '').strip()) if x.strip()]

    bad = {"", "null", "none", "na", "n/a"}
    cleaned = [a for a in raw if a and a.lower() not in bad]

    # Pre-normalization tweak: some mouse inputs want H-2- instead of H2-
    def to_mhcgnomes_mouse(a: str) -> str:
        return re.sub(r'^H2-', 'H-2-', a)

    normalized: list[str] = []
    for a in cleaned:
        try:
            normalized.append(mhcgnomes.parse(to_mhcgnomes_mouse(a)).to_string())
        except Exception:
            # Ignore invalid alleles silently (downstream filters may warn)
            pass

    # Apply preferred star-style formatting for readability/consistency
    normalized = [to_human_star_format(to_mouse_star_format(x)) for x in normalized]

    # Stable de-duplication
    seen = set()
    deduped = []
    for a in normalized:
        if a not in seen:
            seen.add(a)
            deduped.append(a)
    return deduped

def keep_supported_alleles(alleles: list[str], tool: str, supported: list[str]) -> list[str]:
    """
    If a supported-alleles list is available for `tool`, drop alleles not in it.
    Otherwise, return the input list as-is and let the tool validate downstream.
    """
    if not supported:
        return alleles
    tool_ok = [a for a in alleles if a in supported]
    dropped = set(alleles) - set(tool_ok)
    if dropped:
        logging.warning(f"Ignoring not supported alleles for {tool}: {sorted(dropped)}")
    if not tool_ok:
        logging.warning(f"No supported alleles for {tool} remain after filtering.")
    return tool_ok

def has_valid_aas(peptide: str) -> bool:
    """Return True if `peptide` contains only canonical 20 AA letters."""
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    return all(aa in valid for aa in peptide)

def filter_by_length(df: pd.DataFrame, min_len: int, max_len: int, col: str) -> pd.DataFrame:
    """Filter dataframe rows by peptide length in column `col` (inclusive bounds)."""
    return df[df[col].str.len().between(min_len, max_len)]


# --------------------------- CLI --------------------------- #
def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True,
                   help="TSV/CSV with peptides. Default expectation: tab-separated with a 'sequence' column.")
    p.add_argument("--supported-alleles-json", dest="supported_alleles_json", default=None,
                   help="Optional JSON mapping {tool: [alleles,...]}; used to filter alleles per tool.")
    p.add_argument("--prefix", required=True, help="Output file prefix.")
    p.add_argument("--mhc-class", dest="mhc_class", choices=["I", "II"], required=True,
                   help="Controls primary peptide length window: Class I or Class II.")
    p.add_argument("--alleles", default=None,
                   help="Either a path (.txt/.csv/.tsv) with one allele per line, or a semicolon/whitespace-separated list.")
    p.add_argument("--tools", default="mhcflurry",
                   help="Comma-separated list of tools to prepare for (e.g., 'mhcflurry').")
    p.add_argument("--peptide-col-name", dest="peptide_col_name", default="sequence",
                   help="Column name holding the peptide sequences (default: 'sequence').")

    # Lengths (None allowed; defaults are set for Class I, Class II must be explicit)
    p.add_argument("--min-peptide-length-classI",  dest="minI",  type=parse_int_or_none, default=None)
    p.add_argument("--max-peptide-length-classI",  dest="maxI",  type=parse_int_or_none, default=None)
    p.add_argument("--min-peptide-length-classII", dest="minII", type=parse_int_or_none, default=None)
    p.add_argument("--max-peptide-length-classII", dest="maxII", type=parse_int_or_none, default=None)
    return p


# --------------------------- main --------------------------- #
def main():
    args = build_arg_parser().parse_args()

    tools = normalize_tools(args.tools)
    if not tools:
        tools = ["mhcflurry"]
    logging.info(f"Tools: {tools}")

    # Read supported-alleles mapping robustly
    supported_map = parse_supported_json(args.supported_alleles_json)

    # Parse / normalize alleles, or fallback to supported ones per tool
    alleles_norm = parse_alleles(args.alleles)
    if not alleles_norm:
        # Fallback: gather supported alleles from JSON for the selected tools
        fallback = []
        for t in tools:
            supp = supported_map.get(t, [])
            fallback.extend(supp if isinstance(supp, list) else [])
        # Stable de-duplication
        seen = set()
        alleles_norm = [a for a in fallback if not (a in seen or seen.add(a))]

        if alleles_norm:
            logging.warning(
                f"No valid alleles provided; falling back to supported alleles for tools={tools} "
                f"({len(alleles_norm)} total)."
            )
        else:
            raise ValueError(
                "No valid alleles provided and supported_alleles_json is empty for the selected tools. "
                "Please provide --alleles or a non-empty supported_alleles.json."
            )

    # Per-tool allele filtering against supported lists
    tool_allele_input = {
        t: ";".join(keep_supported_alleles(alleles_norm, t, supported_map.get(t, [])))
        for t in tools
    }
    with open(f"{args.prefix}_allele_input.json", "w") as fh:
        json.dump(tool_allele_input, fh)

    # Read peptides
    # Default: tab-separated with a 'sequence' column (nf-core style). Adjust if needed.
    df = pd.read_csv(args.input, sep="\t")
    if args.peptide_col_name not in df.columns:
        raise ValueError(
            f"Peptide column '{args.peptide_col_name}' not found in {args.input}. "
            f"Columns present: {list(df.columns)}"
        )
    logging.info(f"Read {len(df)} rows from {args.input}")

    # Keep only canonical 20-AA sequences
    df = df[df[args.peptide_col_name].apply(has_valid_aas)]

    # Primary class-based length window
    if args.mhc_class == "I":
        min_len = args.minI if args.minI is not None else 9
        max_len = args.maxI if args.maxI is not None else 11
    else:  # Class II
        if args.minII is None or args.maxII is None:
            raise ValueError("Class II selected but min/max peptide length for Class II is not set.")
        min_len = args.minII
        max_len = args.maxII

    df_len = filter_by_length(df, min_len, max_len, args.peptide_col_name)
    if df_len.empty:
        raise ValueError("No peptides left after applying MHC class length filters.")
    logging.info(
        f"{len(df_len)} peptides remain after length filtering "
        f"(Class {args.mhc_class}, {min_len}-{max_len})."
    )

    # Tool-specific secondary filters and outputs
    tool_cfg = {
        "mhcflurry":    {"min": MinLength.MHCFLURRY.value,   "max": MaxLength.MHCFLURRY.value,          "suffix": "mhcflurry_input.csv",    "mhc_class": "I"},
        "mhcnuggets":   {"min": MinLength.MHCNUGGETS.value,  "max": MaxLength.MHCNUGGETS_CLASSI.value,  "suffix": "mhcnuggets_input.tsv",   "mhc_class": "I"},
        "mhcnuggetsii": {"min": MinLength.MHCNUGGETS.value,  "max": MaxLength.MHCNUGGETS_CLASSII.value, "suffix": "mhcnuggetsii_input.tsv", "mhc_class": "II"},
        "netmhcpan":    {"min": MinLength.NETMHCPAN.value,   "max": MaxLength.NETMHCPAN.value,          "suffix": "netmhcpan_input.tsv",    "mhc_class": "I"},
        "netmhciipan":  {"min": MinLength.NETMHCIIPAN.value, "max": MaxLength.NETMHCIIPAN.value,        "suffix": "netmhciipan_input.tsv",  "mhc_class": "II"},
    }

    for t in tools:
        cfg = tool_cfg.get(t)
        if not cfg or cfg["mhc_class"] != args.mhc_class:
            continue

        df_tool = filter_by_length(df_len, cfg["min"], cfg["max"], args.peptide_col_name)
        if df_tool.empty:
            logging.info(f"No peptides match tool-specific length for {t}; skipping.")
            continue

        logging.info(f"Preparing {len(df_tool)} peptides for {t}...")

        # Basic sanity guard for netMHCpan family
        if t in {"netmhcpan", "netmhciipan"} and len(alleles_norm) > MaxNumberOfAlleles.NETMHCPAN.value:
            raise ValueError(f"Number of alleles {len(alleles_norm)} exceeds NetMHCpan limit of {MaxNumberOfAlleles.NETMHCPAN.value}.")

        if t == "mhcflurry":
            alleles_str = tool_allele_input.get(t, "")
            if not alleles_str:
                logging.warning("No supported alleles for mhcflurry after filtering; skipping mhcflurry output.")
                continue
            # MHCflurry expects a 2-column CSV: peptide,allele
            df_out = df_tool[[args.peptide_col_name]].copy()
            df_out["allele"] = [alleles_str.split(";")] * len(df_out)
            df_out = df_out.explode("allele").reset_index(drop=True)
            df_out.rename(columns={args.peptide_col_name: "peptide"}, inplace=True)
            df_out[["peptide", "allele"]].to_csv(f"{args.prefix}_{cfg['suffix']}", index=False)
        else:
            # Default for others here: one-column TSV of sequences (named 'sequence')
            df_tool[[args.peptide_col_name]] \
                .rename(columns={args.peptide_col_name: "sequence"}) \
                .to_csv(f"{args.prefix}_{cfg['suffix']}", sep="\t", header=False, index=False)

    # Minimal versions file (nf-core style)
    versions = {
        "prepare_prediction_input.py": {
            "pandas": pd.__version__,
            "mhcgnomes": getattr(mhcgnomes, "__version__", "unknown"),
        }
    }
    with open("versions.yml", "w") as fh:
        for k, v in versions.items():
            fh.write(f"{k}:\n")
            for kk, vv in v.items():
                fh.write(f"  {kk}: {vv}\n")


if __name__ == "__main__":
    main()
