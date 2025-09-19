#!/usr/bin/env python
"""
Runs MHCnuggets predictions for specified alleles and merges results.

Author: Jonas Scheid
License: MIT
"""
import argparse
import shlex
import logging

import numpy as np
import pandas as pd
from mhcnuggets.src.predict import predict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


class Arguments:
    """
    Parses the arguments, including the ones coming from $task.ext.args.
    """
    def __init__(self) -> None:
        self.input = "$tsv"
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        self.mhc_class = "$meta.mhc_class"
        # Versuche zuerst meta.alleles, sonst meta.alleles_supported
        _als = "$meta.alleles" if "$meta.alleles" != "null" else "$meta.alleles_supported"
        self.alleles = _als.split(";") if _als != "null" else []
        self.parse_ext_args("$task.ext.args")

    def parse_ext_args(self, args_string: str) -> None:
        # skip when there are no extended arguments
        if args_string == "null":
            args_string = ""
        args_list = shlex.split(args_string)
        parser = argparse.ArgumentParser()
        args = parser.parse_args(args_list)
        for attr in vars(args):
            setattr(self, attr, getattr(args, attr))


class Version:
    """
    Parse the versions of the modules used in the script.
    """

    @staticmethod
    def get_versions(modules: list) -> dict:
        """
        This function takes a list of modules and returns a dictionary with the
        versions of each module.
        """
        return {module.__name__: module.__version__ for module in modules}

    @staticmethod
    def format_yaml_like(data: dict, indent: int = 0) -> str:
        """
        Formats a dictionary to a YAML-like string.

        Args:
            data (dict): The dictionary to format.
            indent (int): The current indentation level.

        Returns:
            yaml_str: A string formatted as YAML.
        """
        yaml_str = ""
        for key, value in data.items():
            spaces = "  " * indent
            if isinstance(value, dict):
                yaml_str += f"{spaces}{key}:\\n{Version.format_yaml_like(value, indent + 1)}"
            else:
                yaml_str += f"{spaces}{key}: {value}\\n"
        return yaml_str
def main():
    args = Arguments()

    # --- Eingabe normalisieren: CSV -> eine Peptid-Sequenz pro Zeile ---
    peptides_path = args.input
    try:
        df = pd.read_csv(args.input, sep=None, engine="python")  # autodetect , ; \t
        lower = {c.lower(): c for c in df.columns}
        pep_col = lower.get("peptide") or lower.get("peptide_sequence")
        if pep_col:
            # nur nicht-leere Sequenzen, Whitespace trimmen
            s = df[pep_col].dropna().astype(str).str.strip()
            if not s.empty:
                s.to_csv("peptides.txt", index=False, header=False)
                peptides_path = "peptides.txt"
    except Exception:
        # Wenn das Lesen fehlschlägt, ist es vermutlich schon eine Plain-List
        pass

    if not args.alleles:
        raise SystemExit("No alleles provided to MHCnuggets (meta.alleles / meta.alleles_supported missing)")

    # Vorbereiten / Vorhersagen je Allel
    frames = []
    for allele in args.alleles:
        mhcnuggets_allele = allele.replace('*', '').replace('H2', 'H-2')
        # Ranks nur für Nicht-Maus
        compute_rank = 'H-2' not in mhcnuggets_allele

        predict(
            class_=args.mhc_class,
            peptides_path=peptides_path,   # <— WICHTIG
            mhc=mhcnuggets_allele,
            output=f"{args.prefix}_{allele}.csv",
            rank_output=compute_rank
        )


        if compute_rank:
            tmp_df = pd.read_csv(f"{args.prefix}_{allele}_ranks.csv")
        else:
            tmp_df = pd.read_csv(f"{args.prefix}_{allele}.csv")
            tmp_df["rank"] = np.nan

        tmp_df["allele"] = allele
        frames.append(tmp_df)

    if not frames:
        raise SystemExit("No predictions produced")

    predicted = pd.concat(frames, ignore_index=True)
    outfile = (
        f"{args.prefix}_predicted_mhcnuggets.csv"
        if args.mhc_class == "I"
        else f"{args.prefix}_predicted_mhcnuggetsii.csv"
    )
    predicted.to_csv(outfile, index=False)

    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))
        # No __version__ dunder or similar available, need to hardcode version
        f.write('mhcnuggets: 2.4.0')


if __name__ == "__main__":
    main()
