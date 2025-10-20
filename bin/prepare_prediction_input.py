#!/usr/bin/env python3
"""
Prepare input file for ONE MHC binding prediction method(original-like).

- Normalizes alleles via mhcgnomes
- Filters alleles against supported_alleles.json for the selected pred_method
- Global length filter, then pred_method specific filter
- Writes:
    <prefix>_allele_supported.txt  (semicolon list of normalized alleles)
    <prefix>_{pred_method}_input.tsv/csv  (pred_method-specific peptide input)
    <prefix>_idmap.csv             (ONLY for Mhcnuggets; peptide ↔ allele with allele_id)
"""
import argparse
import json
import logging
from enum import Enum
from pathlib import Path

import pandas as pd 
import mhcgnomes

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

class MinLength(Enum):
    MHCFLURRY = 5
    MHCNUGGETS = 5
    

class MaxLength(Enum):
    MHCFLURRY = 15
    MHCNUGGETS_CLASSI = 15
    MHCNUGGETS_CLASSII = 30


def has_valid_aas(s: str) -> bool:
    return all(c in "ACDEFGHIKLMNPQRSTVWY" for c in s)

def filter_by_length(df: pd.DataFrame, lo: int, hi: int, col: str) -> pd.DataFrame:
    return df[df[col].str.len().between(lo, hi)]


class Version:
    @staticmethod
    def get_versions(mods: list) -> dict:
        return {m.__name__: getattr(m,"__version__","NA") for m in mods}
    @staticmethod
    def format_yaml_like(data: dict, indent: int = 0) -> str:
        s = ""
        for k,v in data.items():
            pad = "  "*indent
            if isinstance(v, dict):
                s += f"{pad}{k}:\n{Version.format_yaml_like(v, indent+1)}"
            else:
                s += f"{pad}{k}: {v}\n"
        return s

def parse_arguments():
    parser = argparse.ArgumentParser(description="Prepare input for MHC binding prediction")
    
    parser.add_argument("--input", required=True, help="Input TSV file with peptides")
    parser.add_argument("--prefix", required=True, help="Output file prefix")
    parser.add_argument("--pred_method", required=True, help="Prediction method")
    parser.add_argument("--allele_id", default="", help="Allele ID")
    parser.add_argument("--allele_name", default="", help="Allele name")
    parser.add_argument("--mhc_class", default="", help="MHC class")
    parser.add_argument("--min_pep_len", type=int, required=True, help="Minimum peptide length")
    parser.add_argument("--max_pep_len", type=int, required=True, help="Maximum peptide length")
    parser.add_argument("--alleles", required=True, help="Alleles file")
    parser.add_argument("--supported_alleles_json", required=True, help="Supported alleles JSON file")
    parser.add_argument("--peptide_col_name", default="sequence", help="Name of peptide column")
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    pred_method = (args.pred_method or "").strip().lower()
    # Leite Klasse automatisch ab, falls nicht gesetzt
    mhc_class = (args.mhc_class or ("II" if pred_method == "mhcnuggets-class-2" else "I")).strip()
    logging.info("pred_method=%s -> pred_method=%s, class=%s", args.pred_method, pred_method, mhc_class)

    # Sanitize arguments
    args.allele_id = int(args.allele_id)
    


    # Load supported_alleles.json
    sup = Path(args.supported_alleles_json)
    with sup.open() as fh:
       supported_map = json.load(fh)
   

    supported_for_pred_method = supported_map.get(pred_method)
    if supported_for_pred_method is None:
       raise KeyError(f"No supported allele list for pred_method '{pred_method}' in supported_alleles.json.")


    # --------------------------
    # Read allele table and normalize like original
    # --------------------------

    alleles_tbl = Path(args.alleles)
    df_al = pd.read_csv(alleles_tbl, sep=None, engine="python", comment="#")

    # First column is ID, second is name
    df_al.columns = ["allele_id", "allele_name"]
    cols = {"allele_id": "allele_id", "allele_name": "allele_name"}

    # normalize via mhcgnomes
    df_al["allele_name_norm"] = (
        df_al[cols["allele_name"]]
        .astype(str).str.strip()
        .apply(lambda a: mhcgnomes.parse(a).to_string())
    )

    # Identify and warn about unsupported alleles
    unsupported_alleles = df_al[~df_al["allele_name_norm"].isin(supported_for_pred_method)]["allele_name"].tolist()
    if unsupported_alleles:
        print(f"Warning: The following alleles are not supported by {pred_method}: {', '.join(unsupported_alleles)}")

    # Filter to keep only supported alleles  
    df_al = df_al[df_al["allele_name_norm"].isin(supported_for_pred_method)]
    if df_al.empty:
       raise ValueError(f"None of the provided alleles are supported by {pred_method}.")

    # Select exactly one allele for every chunk
    sel_row = None
    if args.allele_id is not None:
        # allele_id ist immer int → direkt numerisch vergleichen
         sel_row = df_al[df_al[cols["allele_id"]] == int(args.allele_id)]
    
    allele_name_norm = sel_row["allele_name_norm"].iloc[0]
    allele_id_val    = int(sel_row[cols["allele_id"]].iloc[0])

    # per-chunk allele artifacts for downstream
    with open(f"{args.prefix}_allele_supported.txt", "w") as fh:
        fh.write(allele_name_norm + "\n")
    


    # -----------------------------
    # Read peptides TSV 
    # -----------------------------

    df = pd.read_csv(args.input, sep=None, engine="python", comment="#")
    df = df.rename(columns={"peptide_sequence": "sequence"})

    # Keep only valid amino acids
    df = df[df["sequence"].apply(has_valid_aas)] #inwiefern wurde das vorher schon geprüft?

    # Global length filter 
    df_len = filter_by_length(df, args.min_pep_len, args.max_pep_len, "sequence")
    if df_len.empty:
        raise ValueError("No peptides left after global length filtering.")

    #  pred_method-specific extra filter & I/O 
    pred_method_cfg = {
        "mhcflurry":     {"min": MinLength.MHCFLURRY.value,   "max": MaxLength.MHCFLURRY.value,          "suffix":"mhcflurry_input.csv",    "mhc_class":"I"},
        "mhcnuggets-class-1":    {"min": MinLength.MHCNUGGETS.value,  "max": MaxLength.MHCNUGGETS_CLASSI.value,  "suffix":"mhcnuggets_input.tsv",   "mhc_class":"I"},
        "mhcnuggets-class-2":  {"min": MinLength.MHCNUGGETS.value,  "max": MaxLength.MHCNUGGETS_CLASSII.value, "suffix":"mhcnuggets_input.tsv", "mhc_class":"II"},
    }

    cfg = pred_method_cfg.get(pred_method)
    if not cfg:
        raise ValueError(f"Unsupported pred_method: {pred_method}")
    if cfg["mhc_class"] != mhc_class:
        logging.warning("Selected pred_method '%s' is for class %s, but mhc_class is '%s'.",
                        pred_method, cfg["mhc_class"], mhc_class)

    df_pred_method = filter_by_length(df_len, cfg["min"], cfg["max"], "sequence")
    if df_pred_method.empty:
        raise ValueError(f"No peptides for {pred_method} after pred_method-specific length filter.")

    logging.info("Preparing %d peptides for %s (allele_id=%s, allele=%s) ...",
                 len(df_pred_method), pred_method, allele_id_val, allele_name_norm)

    # Write per pred_method inputs (single allele for this chunk) 
    if pred_method == "mhcflurry":
        # Build table with peptide+allele and attach stable IDs (needed by POSTPROCESS_PREDICTION)
        if "peptide_id" in df_pred_method.columns:
            prep = df_pred_method[["peptide_id", "sequence"]].copy()
            prep["peptide_id"] = prep["peptide_id"].astype(str)
        else:
            tmp = df_pred_method[["sequence"]].reset_index().rename(columns={"index": "peptide_id"})
            tmp["peptide_id"] = tmp["peptide_id"].astype(str)
            prep = tmp[["peptide_id", "sequence"]]

        prep = prep.rename(columns={"sequence": "peptide"})
        prep["allele"]    = allele_name_norm
        prep["allele_id"] = allele_id_val

        # This file is used both by the mhcflurry (it reads peptide+allele)
        # and by POSTPROCESS_PREDICTIONS (needs peptide_id/allele_id to reduce)
        prep[["peptide", "allele", "peptide_id", "allele_id"]].to_csv(
            f"{args.prefix}_{cfg['suffix']}", index=False
        )

    else: # mhcnuggets-class-1 / mhcnuggets-class-2
        # one-column 'sequence' file for predictor
        df_pred_method[["sequence"]].to_csv(f"{args.prefix}_{cfg['suffix']}", sep="\t", header=False, index=False)

        # idmap is only created for mhcnuggets (peptide, allele, peptide_id, allele_id)
        if "peptide_id" in df_pred_method.columns:
            idmap = df_pred_method[["peptide_id", "sequence"]].copy()
            idmap["peptide_id"] = idmap["peptide_id"].astype(str)
        else:
            tmp = df_pred_method[["sequence"]].reset_index().rename(columns={"index": "peptide_id"})
            tmp["peptide_id"] = tmp["peptide_id"].astype(str)
            idmap = tmp[["peptide_id", "sequence"]]

        idmap = idmap.rename(columns={"sequence": "peptide"})
        idmap["allele"]    = allele_name_norm
        idmap["allele_id"] = allele_id_val
        idmap[["peptide", "allele", "peptide_id", "allele_id"]].to_csv(
            f"{args.prefix}_idmap.csv", index=False
        )

    
    # versions
    versions_this_module = {}
    versions_this_module[f"PREPARE_PREDICTION_INPUT"] = Version.get_versions([pd, mhcgnomes])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))

if __name__ == "__main__":
    main()
