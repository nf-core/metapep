#!/usr/bin/env python
"""
Runs MHCnuggets predictions for specified alleles and merges results.

Author: Jonas Scheid (adapted)
"""
import argparse
import shlex
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from mhcnuggets.src.predict import predict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

class Arguments:

    def __init__(self) -> None:
        # output PREPARE_PREDICTION_INPUT: *_mhcnuggets_input.tsv (Single-column list presenting only sequences; without header)
        self.input = "$tsv"
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        self.mhc_class = "$meta.mhc_class"

        raw_supported = "$meta.alleles_supported"   # from {meta.id}_allele_supported.txt
        if not raw_supported or raw_supported == "null":
             raise SystemExit("No alleles provided in meta.alleles_supported")

        # Semicolon-separated list (one entry per chunk)
        self.alleles = [a.strip() for a in raw_supported.split(";") if a.strip()]

    def parse_ext_args(self, args_string: str) -> None:
        if args_string == "null":
            args_string = ""
        args_list = shlex.split(args_string)
        parser = argparse.ArgumentParser()
        args = parser.parse_args(args_list)
        for attr in vars(args):
            setattr(self, attr, getattr(args, attr))


class Version:
    @staticmethod
    def get_versions(modules: list) -> dict:
        """
        This function takes a list of modules and returns a dictionary with the
        versions of each module.
        """
        #return {module.__name__: getattr(module, "__version__", "NA") for module in modules}
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

    peptides_path = args.input

    logging.info(f"Running MHCnuggets class={args.mhc_class} on {peptides_path}")
    logging.info("First alleles (raw)   : %s", args.alleles[:5])

    predicted_df = []
    for allele in args.alleles:
        # Format for MHCnuggets: remove asterisk, mouse 'H2' -> 'H-2'
        mhcnuggets_allele = allele.replace("*", "").replace("H2", "H-2")
        # MHCnuggets cannot compute ranks for mouse alleles
        compute_rank = "H-2" not in mhcnuggets_allele

        predict(
            class_=args.mhc_class,
            peptides_path=peptides_path,
            mhc=mhcnuggets_allele,
            output=f"{args.prefix}_{allele}.csv",
            rank_output=compute_rank,
        )

        if compute_rank:
            df = pd.read_csv(f"{args.prefix}_{allele}_ranks.csv")
        else:
            df = pd.read_csv(f"{args.prefix}_{allele}.csv")
            df["rank"] = np.nan  # rank-Spalte nachziehen

        # Keep the allele as provided by mhcgnomes for the downstream join
        df["allele"] = allele
        predicted_df.append(df)

    predicted_df = pd.concat(predicted_df, ignore_index=True)
    filename_out = (
        f"{args.prefix}_predicted_mhcnuggets.csv"
        if args.mhc_class == "I"
        else f"{args.prefix}_predicted_mhcnuggetsii.csv"
    )
    predicted_df.to_csv(filename_out, index=False)

    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, pd])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))
        f.write("mhcnuggets: 2.4.0")


if __name__ == "__main__":
    main()
