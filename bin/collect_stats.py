#!/usr/bin/env python3

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
    parser.add_argument(
        "-c",
        "--conditions",
        help="Path to the conditions input file",
        type=str,
        required=True)
    parser.add_argument(
        "-a",
        "--alleles",
        help="Path to the allele input file",
        type=str,
        required=True)
    parser.add_argument(
        "-pr",
        "--predictions",
        help="Path to the predictions input file",
        type=str,
        required=True)
    parser.add_argument(
        "-bt",
        "--binder_threshold",
        help="Score threshold to call binder",
        type=float,
        default=0.5,
        required=True)

    # OUTPUT FILES
    parser.add_argument("-o", "--outfile", help="Path to the output file", type=str, required=True)

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
    conditions["microbiome_id"] = pd.to_numeric(conditions["microbiome_id"], downcast="unsigned")

    predictions = pd.read_csv(args.predictions, usecols=["peptide_id", "prediction_score", "allele_id"], sep="\t")
    predictions["peptide_id"] = pd.to_numeric(predictions["peptide_id"], downcast="unsigned")
    predictions["allele_id"] = pd.to_numeric(predictions["allele_id"], downcast="unsigned")
    predictions["prediction_score"] = pd.to_numeric(predictions["prediction_score"], downcast="float")
    predictions = predictions[predictions["prediction_score"]>=args.binder_threshold]

    alleles = pd.read_csv(args.alleles, usecols=["allele_id", "allele_name"], sep="\t").set_index("allele_id")
    allele_names = alleles["allele_name"].to_dict()

    predictions = pd.concat([df.set_index("peptide_id").rename(columns={"prediction_score":f"prediction_score_allele_{allele_id}"}).drop("allele_id", axis=1) for allele_id, df in predictions.groupby("allele_id")], join="outer", axis=1)

    stat_table = []

    # Process data condition-wise to reduce memory usage
    for condition_name in conditions["condition_name"]:

        # TODO solve same entity between conditions -> duplication of data within merge
        conditions_proteins = (
            conditions[conditions.condition_name == condition_name]
            .merge(microbiomes_entities_occs)
            .drop(columns="microbiome_id")
            .merge(entities_proteins_occs)
            .drop(columns="entity_id")
        )

        # condition_name, unique_proteins
        unique_protein_count = len(conditions_proteins[["condition_name", "protein_id"]].drop_duplicates())

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

        # condition_name, unique_peptide_count
        unique_peptide_count = len(conditions_peptides)

        # peptide_id, prediction_score_allele_0, prediction_score_allele_1,...prediction_score_allele_n
        conditions_peptides = conditions_peptides.set_index("peptide_id").drop(["condition_name","condition_peptide_count"], axis=1).join(predictions, how="inner")

        # number of best prediction allele binder
        num_binder = len(conditions_peptides)

        # collect all info into table row per condition
        row_temp = {"Condition name":condition_name, "Unique proteins":unique_protein_count, "Total peptides":total_peptide_count, "Unique peptides":unique_peptide_count, "# Binder (any allele)": num_binder}
        for col in conditions_peptides.columns:
            allele = allele_names[int(col.split('_')[-1])]
            row_temp[f"# Binders for allele {allele}"] = len(conditions_peptides[col].dropna())

        stat_table.append(row_temp)

    stat_table = pd.DataFrame(stat_table)

    # collect info across all conditions
    all_conditions_unique_protein_count = len(protein_peptide_occs["protein_id"].drop_duplicates())
    all_conditions_total_peptide_counts = stat_table["Total peptides"].sum()
    all_conditions_unique_peptide_counts = len(protein_peptide_occs["peptide_id"].drop_duplicates())
    all_conditions_unique_binder_counts = len(predictions)

    row_total = {"Condition name":"total", "Unique proteins":all_conditions_unique_protein_count, "Total peptides":all_conditions_total_peptide_counts, "Unique peptides":all_conditions_unique_peptide_counts, "# Binder (any allele)": all_conditions_unique_binder_counts}
    for col in conditions_peptides.columns:
        allele = allele_names[int(col.split('_')[-1])]
        row_total[f"# Binders for allele {allele}"] = len(conditions_peptides[col].dropna())

    stat_table = stat_table.append(row_total, ignore_index=True)
    stat_table.to_csv(args.outfile, sep="\t", index=False)

if __name__ == "__main__":
    sys.exit(main())
