#!/usr/bin/env python3
#
# Author: Sabrina Krakau
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################

import argparse
import sys
import os
import csv

import pandas as pd


def parse_args(args=None):
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Collect some stats.")

    # INPUT FILES
    parser.add_argument(
        "-ppo",
        "--protein-peptide-occ",
        help="Path to the protein peptide occurrences input file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-epo",
        "--entities-proteins-occ",
        help="Path to the entity protein occurrences input file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-meo",
        "--microbiomes-entities-occ",
        help="Path to the microbiome entity occurrences input file",
        type=str,
        required=True,
    )
    parser.add_argument("-c", "--conditions", help="Path to the conditions input file", type=str, required=True)

    # OUTPUT FILES
    parser.add_argument("-o", "--outfile", help="Path to the output file", type=argparse.FileType("w"), required=True)

    return parser.parse_args()


def main(args=None):
    # Parse command line arguments
    args = parse_args(args)

    # Read input files
    protein_peptide_occs = pd.read_csv(args.protein_peptide_occ, sep="\t")
    entities_proteins_occs = pd.read_csv(args.entities_proteins_occ, sep="\t")
    microbiomes_entities_occs = pd.read_csv(args.microbiomes_entities_occ, sep="\t").drop(columns="entity_weight")
    conditions = pd.read_csv(args.conditions, sep="\t").drop(columns="condition_id")

    print("Joining input data...", flush=True)

    conditions_proteins = (
        conditions.merge(microbiomes_entities_occs)
        .drop(columns="microbiome_id")
        .merge(entities_proteins_occs)
        .drop(columns="entity_id")
    )

    # condition_name, unique_proteins
    unique_protein_counts = (
        conditions_proteins[["condition_name", "protein_id"]]
        .drop_duplicates()
        .groupby("condition_name")
        .size()
        .reset_index(name="unique_protein_count")
    )

    print("Unique protein counts:", file=args.outfile, sep="\t", flush=True)
    print("condition_name", "unique_protein_count", file=args.outfile, sep="\t", flush=True)
    for index, row in unique_protein_counts.iterrows():
        print(row["condition_name"], row["unique_protein_count"], file=args.outfile, sep="\t", flush=True)

    # condition_name, peptide_id, condition_peptide_count
    conditions_peptides = (
        conditions_proteins.merge(protein_peptide_occs)
        .drop(columns="protein_id")
        .groupby(["condition_name", "peptide_id"])["count"]
        .sum()
        .reset_index(name="condition_peptide_count")
    )

    # condition_name, total_peptide_count
    total_peptide_counts = (
        conditions_peptides.groupby("condition_name")["condition_peptide_count"]
        .sum()
        .reset_index(name="total_peptide_count")
    )

    print(file=args.outfile)
    print("Total peptide counts:", file=args.outfile, sep="\t", flush=True)
    print("condition_name", "total_peptide_count", file=args.outfile, sep="\t", flush=True)
    for index, row in total_peptide_counts.iterrows():
        print(row["condition_name"], row["total_peptide_count"], file=args.outfile, sep="\t", flush=True)

    # condition_name, unique_peptide_count
    unqiue_peptide_counts = (
        conditions_peptides.groupby("condition_name").size().reset_index(name="unique_peptide_count")
    )

    print(file=args.outfile)
    print("Unique peptide counts:", file=args.outfile, sep="\t", flush=True)
    print("condition_name", "unique_peptide_count", file=args.outfile, sep="\t", flush=True)
    for index, row in unqiue_peptide_counts.iterrows():
        print(row["condition_name"], row["unique_peptide_count"], file=args.outfile, sep="\t", flush=True)

    # unique peptides across all conditions
    all_conditions_unqiue_peptide_counts = len(conditions_peptides["peptide_id"].drop_duplicates())

    print(file=args.outfile)
    print(
        "Unique peptides across all conditions:",
        all_conditions_unqiue_peptide_counts,
        file=args.outfile,
        sep="\t",
        flush=True,
    )

    print("Done!", flush=True)


if __name__ == "__main__":
    sys.exit(main())
