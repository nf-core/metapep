#!/usr/bin/env python3
# Written by Sabrina Krakau, Leon Kuchenbecker, and Till Englert under the MIT license

import argparse
import datetime
import os
import sys

import pandas as pd

####################################################################################################


def parse_args(args=None):
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Prepare entity binding rates for plotting.")

    # INPUT FILES
    parser.add_argument("-p", "--predictions", help="Path to the predictions input file", type=str, required=True)
    parser.add_argument(
        "-ppo",
        "--protein-peptide-occ",
        help="Path to the protein peptide occurences input file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-epo",
        "--entities-proteins-occ",
        help="Path to the entity protein occurences input file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-meo",
        "--microbiomes-entities-occ",
        help="Path to the microbiome entity occurences input file",
        type=str,
        required=True,
    )
    parser.add_argument("-c", "--conditions", help="Path to the conditions input file", type=str, required=True)
    parser.add_argument(
        "-cam", "--condition-allele-map", help="Path to the condition allele map input file", type=str, required=True
    )
    parser.add_argument("-a", "--alleles", help="Path to the allele input file", type=str, required=True)
    parser.add_argument("-m", "--method", help="Used epitope prediction method", type=str, required=True)

    # OUTPUT FILES
    parser.add_argument("-o", "--outdir", help="Path to the output directory", type=str, required=True)

    # PARAMETERS
    parser.add_argument(
        "-pc",
        "--chunk-size",
        help=(
            "Chunk size with respect to peptide_ids used for internal processing to limit memory usage. Default: 500000"
        ),
        type=int,
        default=500000,
    )
    parser.add_argument(
        "-sst",
        "--syfpeithi_score_threshold",
        help=("Threshold for binder/non-binder calling when using SYFPEITHI epitope prediction method. Default: 0.5"),
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "-mst",
        "--mhcf_mhcn_score_threshold",
        help=(
            "Threshold for binder/non-binder calling when using MHCflurry or MHCnuggets epitope prediction methods. Default: 0.426"
        ),
        type=float,
        default=0.426,
    )
    parser.add_argument(
        "-mlld",
        "--mem_log_level_deep",
        help="Enable 'deep' option for pandas memory usage output ('deep' enables accurate usage values, but increases runtime). Default: None. ",
        default=False,
        action="store_true",
    )

    return parser.parse_args()


def call_binder(score, method, syfpeithi_score_threshold, mhcfn_score_threshold):
    """
    Scoring threshold is based on the nf-core/epitopeprediction pipeline.
    For SYFPEITHI the scoring threshold is a "half of maximum score". After
    normalization the highest achievable score is 1. For MHCflurry and
    MHCnuggets the score is an 0 to 1 scoring base on the
    affinity score (IC50) and is calculated by: 1-log_50000(affinity_score)
    in this scale the old threshold of 500 is: 0.426 and the higher the better.
    """
    if method == "syfpeithi":
        return score >= syfpeithi_score_threshold
    else:
        return score >= mhcfn_score_threshold


def main(args=None):
    args = parse_args(args)

    now = datetime.datetime.now()
    print("Start date and time : ")
    print(now.strftime("%Y-%m-%d %H:%M:%S"))

    # Read input files
    predictions = pd.read_csv(args.predictions, sep="\t", index_col="peptide_id").sort_index()
    predictions["allele_id"] = pd.to_numeric(predictions["allele_id"], downcast="unsigned")
    predictions["prediction_score"] = pd.to_numeric(predictions["prediction_score"], downcast="float")
    # NOTE could be read in chunk-wise if this gets bottleneck.
    # E.g. sort and split before by peptide_ids into chunks, read in as file list, process chunk-wise (protein_peptide_occs accordingly)
    # (would decrease code readability though)
    protein_peptide_occs = pd.read_csv(args.protein_peptide_occ, sep="\t", index_col="peptide_id").sort_index()
    protein_peptide_occs["protein_id"] = pd.to_numeric(protein_peptide_occs["protein_id"], downcast="unsigned")
    protein_peptide_occs["count"] = pd.to_numeric(protein_peptide_occs["count"], downcast="unsigned")

    entities_proteins_occs = pd.read_csv(args.entities_proteins_occ, sep="\t")
    entities_proteins_occs["entity_id"] = pd.to_numeric(entities_proteins_occs["entity_id"], downcast="unsigned")
    entities_proteins_occs["protein_id"] = pd.to_numeric(entities_proteins_occs["protein_id"], downcast="unsigned")

    microbiomes_entities_occs = pd.read_csv(args.microbiomes_entities_occ, sep="\t")
    microbiomes_entities_occs["microbiome_id"] = pd.to_numeric(
        microbiomes_entities_occs["microbiome_id"], downcast="unsigned"
    )
    microbiomes_entities_occs["entity_id"] = pd.to_numeric(microbiomes_entities_occs["entity_id"], downcast="unsigned")
    microbiomes_entities_occs["entity_weight"] = pd.to_numeric(
        microbiomes_entities_occs["entity_weight"], downcast="float"
    )

    conditions = pd.read_csv(args.conditions, sep="\t")
    conditions["condition_id"] = pd.to_numeric(conditions["condition_id"], downcast="unsigned")
    conditions["microbiome_id"] = pd.to_numeric(conditions["microbiome_id"], downcast="unsigned")

    condition_allele_map = pd.read_csv(args.condition_allele_map, sep="\t")
    condition_allele_map["condition_id"] = pd.to_numeric(condition_allele_map["condition_id"], downcast="unsigned")
    condition_allele_map["allele_id"] = pd.to_numeric(condition_allele_map["allele_id"], downcast="unsigned")

    alleles = pd.read_csv(args.alleles, sep="\t")
    alleles["allele_id"] = pd.to_numeric(alleles["allele_id"], downcast="unsigned")

    if args.mem_log_level_deep:
        print_mem = "deep"
    else:
        print_mem = None

    print("\nInfo: predictions", flush=True)
    predictions.info(verbose=False, memory_usage=print_mem)
    print("\nInfo: protein_peptide_occs", flush=True)
    protein_peptide_occs.info(verbose=False, memory_usage=print_mem)

    # Create output directory if it doesn't exist
    if os.path.exists(args.outdir) and not os.path.isdir(args.outdir):
        print("ERROR - The target path is not a directory", file=sys.stderr)
        sys.exit(2)
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("Prepare df with protein info ...", file=sys.stderr, flush=True, end="")
    # Prepare df for joining protein information
    protein_info = (
        microbiomes_entities_occs.merge(conditions)
        .drop(columns="microbiome_id")
        .merge(condition_allele_map)
        .drop(columns="condition_id")
        .merge(entities_proteins_occs)
    )
    # -> protein_id, entity_id, entity_weight, condition_name, allele_id
    # merged against condition_allele_map to keep only entities, and thus proteins, for which a prediction is requested for the current allele
    print("\nInfo: protein_info", flush=True)
    protein_info.info(verbose=False, memory_usage=print_mem)

    # Process predictions chunk-wise based on peptide_ids
    # (predictions chunk can contain more than chunk_size rows due to multiple alleles)
    entity_results = pd.DataFrame()
    max_peptide_id = predictions.index.max()
    for i in range(0, max_peptide_id, args.chunk_size):
        print("\nChunk peptide_ids: ", i, " - ", i + args.chunk_size - 1)

        now = datetime.datetime.now()
        print("Time: ...")
        print(now.strftime("%Y-%m-%d %H:%M:%S"))

        # Join predictions with protein_ids and further protein info
        # NOTE would be faster to process directly respective chunks (see comment above)
        data = (
            predictions[(predictions.index >= i) & (predictions.index < i + args.chunk_size)]
            .join(
                protein_peptide_occs[
                    (protein_peptide_occs.index >= i) & (protein_peptide_occs.index < i + args.chunk_size)
                ]
            )
            .merge(protein_info)
        )  # based on protein_id, allele_id
        # -> index, prediction_score, allele_id, protein_id, count, entity_id, entity_weight, condition_name
        # (protein_id could be dropped, but no big impact here)

        # Call binder based on prediction_score
        data["binder"] = data["prediction_score"].apply(
            call_binder,
            method=args.method,
            syfpeithi_score_threshold=args.syfpeithi_score_threshold,
            mhcfn_score_threshold=args.mhcf_mhcn_score_threshold,
        )
        data.drop(columns="prediction_score", inplace=True)

        # Count total number of peptides and number of binders for each entity, allele and condition (including multiple counts within proteins)
        # Create extra column containing only counts for peptides that are classified as 'binder' to allow efficient aggregation
        data["binder_count"] = 0
        data.loc[data["binder"], "binder_count"] = data["count"]
        data = (
            data.groupby(["entity_id", "allele_id", "condition_name", "entity_weight"], group_keys=False)
            .agg(
                count_binders=pd.NamedAgg(column="binder_count", aggfunc="sum"),
                count_peptides=pd.NamedAgg(column="count", aggfunc="sum"),
            )
            .reset_index()
        )
        # -> index, entity_id, allele_id, condition_name, entity_weight, count_binders, count_peptides

        # Append to entity results
        entity_results = pd.concat([entity_results, data], ignore_index=True)
        print("\nInfo: entity_results", flush=True)
        entity_results.info(verbose=False, memory_usage=print_mem)

    # Combine entity results from different chunks
    entity_results = (
        entity_results.groupby(["entity_id", "allele_id", "condition_name", "entity_weight"], group_keys=False)[
            ["count_binders", "count_peptides"]
        ]
        .sum()
        .reset_index()
    )
    entity_results["binding_rate"] = entity_results["count_binders"] / entity_results["count_peptides"]
    data.drop(columns=["count_binders", "count_peptides"], inplace=True)

    # Write out results for each allele
    for allele_id in alleles.allele_id:
        with open(os.path.join(args.outdir, "entity_binding_ratios.allele_" + str(allele_id) + ".tsv"), "w") as outfile:
            entity_results[entity_results.allele_id == allele_id][
                ["condition_name", "binding_rate", "entity_weight"]
            ].to_csv(outfile, sep="\t", index=False, header=True)
    print("Done!", flush=True)


if __name__ == "__main__":
    sys.exit(main())
