#!/usr/bin/env python3
# Written by Sabrina Krakau, Leon Kuchenbecker, and Till Englert under the MIT license

import argparse
import sys

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
    protein_peptide_occs["protein_id"] = pd.to_numeric(protein_peptide_occs["protein_id"], downcast="unsigned")
    protein_peptide_occs["peptide_id"] = pd.to_numeric(protein_peptide_occs["peptide_id"], downcast="unsigned")
    protein_peptide_occs["count"] = pd.to_numeric(protein_peptide_occs["count"], downcast="unsigned")

    entities_proteins_occs = pd.read_csv(args.entities_proteins_occ, sep="\t")
    entities_proteins_occs["entity_id"] = pd.to_numeric(entities_proteins_occs["entity_id"], downcast="unsigned")
    entities_proteins_occs["protein_id"] = pd.to_numeric(entities_proteins_occs["protein_id"], downcast="unsigned")

    microbiomes_entities_occs = pd.read_csv(
        args.microbiomes_entities_occ, usecols=["microbiome_id", "entity_id"], sep="\t"
    )
    microbiomes_entities_occs["microbiome_id"] = pd.to_numeric(
        microbiomes_entities_occs["microbiome_id"], downcast="unsigned"
    )
    microbiomes_entities_occs["entity_id"] = pd.to_numeric(microbiomes_entities_occs["entity_id"], downcast="unsigned")

    conditions = pd.read_csv(args.conditions, usecols=["condition_name", "microbiome_id"], sep="\t")
    # TODO "condition_name" as " "category"?
    # (first try: memory went up, check how to use properly)
    conditions["microbiome_id"] = pd.to_numeric(conditions["microbiome_id"], downcast="unsigned")

    # Process data condition-wise to reduce memory usage
    for condition_name in conditions["condition_name"]:
        print("Process condition:", condition_name, flush=True)
        print("Condition name:", condition_name, file=args.outfile, sep="\t", flush=True)

        conditions_proteins = (
            conditions[conditions.condition_name == condition_name]
            .merge(microbiomes_entities_occs)
            .drop(columns="microbiome_id")
            .merge(entities_proteins_occs)
            .drop(columns="entity_id")
        )

        # condition_name, unique_proteins
        unique_protein_count = len(conditions_proteins[["condition_name", "protein_id"]].drop_duplicates())
        print("Unique proteins:", unique_protein_count, file=args.outfile, sep="\t", flush=True)

        # condition_name, peptide_id, condition_peptide_count
        conditions_peptides = (
            conditions_proteins.merge(protein_peptide_occs)
            .drop(columns="protein_id")
            .groupby(["condition_name", "peptide_id"])["count"]
            .sum()
            .reset_index(name="condition_peptide_count")
        )

        # condition_name, total_peptide_count
        total_peptide_count = sum(conditions_peptides["condition_peptide_count"])
        print("Total peptides:", total_peptide_count, file=args.outfile, sep="\t", flush=True)

        # condition_name, unique_peptide_count
        unique_peptide_count = len(conditions_peptides)
        print("Unique peptides:", unique_peptide_count, file=args.outfile, sep="\t", flush=True)
        print(file=args.outfile)

    # unique peptides across all conditions
    all_conditions_unqiue_peptide_counts = len(protein_peptide_occs["peptide_id"].drop_duplicates())
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
