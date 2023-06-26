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
import numpy as np
import datetime

####################################################################################################


def parse_args(args=None):
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Prepare prediction score distribution for plotting.")

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
        "-mlld",
        "--mem_log_level_deep",
        help="Enable 'deep' option for pandas memory usage output ('deep' enables accurate usage values, but increases runtime). Default: None. ",
        default=False,
        action="store_true",
    )

    return parser.parse_args()


def main(args=None):
    args = parse_args(args)
    if args.mem_log_level_deep:
        print_mem = "deep"
    else:
        print_mem = None

    now = datetime.datetime.now()
    print("Start date and time : ")
    print(now.strftime("%Y-%m-%d %H:%M:%S"))

    # Read input files
    predictions = pd.read_csv(args.predictions, sep="\t", index_col="peptide_id").sort_index()
    predictions["allele_id"] = pd.to_numeric(predictions["allele_id"], downcast="unsigned")
    predictions["prediction_score"] = pd.to_numeric(predictions["prediction_score"], downcast="float")

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

    # Create prediction_score bin intervals to reduce number of data points for plotting
    bin_results, bin_edges = pd.cut(predictions["prediction_score"], bins=1000, retbins=True)

    print("Prepare df with protein info ...", flush=True)
    # Prepare df for joining protein information
    # (get condition_name, entity_weights and filter for condition entities)
    protein_info = (
        microbiomes_entities_occs.merge(conditions)
        .drop(columns="microbiome_id")
        .merge(condition_allele_map)
        .drop(columns="condition_id")
        .merge(entities_proteins_occs)
        .drop(columns="entity_id")
    )
    # -> protein_id, entity_weight, condition_name, allele_id
    # merged against condition_allele_map to keep only entities, and thus proteins, for which a prediction is requested for the current allele
    print("\nInfo: protein_info", flush=True)
    protein_info.info(verbose=False, memory_usage=print_mem)

    # Prepare output files for each allele
    outfile_dict = {}
    for allele_id in alleles.allele_id:
        outfile = open(os.path.join(args.outdir, "prediction_scores.allele_" + str(allele_id) + ".tsv"), "w")
        outfile_dict[allele_id] = outfile
    print_header = ["prediction_score", "condition_name", "weight_sum"]

    try:
        # Process predictions chunk-wise based on peptide_ids
        # (predictions chunk can contain more than chunk_size rows due to multiple alleles)
        max_peptide_id = predictions.index.max()
        for chunk_start in range(0, max_peptide_id, args.chunk_size):
            chunk_end = min(chunk_start + args.chunk_size - 1, max_peptide_id)
            print("\nChunk peptide_ids: ", chunk_start, " - ", chunk_end)

            now = datetime.datetime.now()
            print("Time: ...")
            print(now.strftime("%Y-%m-%d %H:%M:%S"))

            data = (
                predictions[(predictions.index >= chunk_start) & (predictions.index <= chunk_end)]
                .join(
                    protein_peptide_occs[
                        (protein_peptide_occs.index >= chunk_start) & (protein_peptide_occs.index <= chunk_end)
                    ]
                )
                .reset_index(names="peptide_id")
                .merge(protein_info)
                .drop(columns="protein_id")
            )
            # (merge() might change index, so reset_index(..) before)
            # -> index, peptide_id, prediction_score, allele_id, count, entity_weight, condition_name

            # Add for each prediction, protein and entity aggregated weights containing entity_weights * counts
            data["agg_weight"] = data["entity_weight"] * data["count"]
            data.drop(columns=["entity_weight", "count"], inplace=True)

            print("\nInfo: data after merging protein_info", flush=True)
            data.info(verbose=False, memory_usage=print_mem)

            data = (
                data.groupby(["peptide_id", "prediction_score", "condition_name", "allele_id"])["agg_weight"]
                .sum()
                .reset_index(name="weight_sum")
                .drop(columns="peptide_id")
            )
            # -> index, prediction_score, condition_name, allele_id, weight_sum

            # NOTE
            # for each peptide in a condition the weight is computed as follows:
            # - the sum of all weights of the corresponding entity_weights, each weighted by the number of proteins in which the peptide occurs
            # - multiple occurrences of the peptide within one protein are counted
            # if the no entity weights are provided, weight_sum corresponds to the number of peptide occurrences within the condition

            print("\nInfo: final data", flush=True)
            data.info(verbose=False, memory_usage=print_mem)

            # Assign prediction_scores to precomputed bins, group and sum up weights
            data["prediction_score_bin"] = pd.cut(data["prediction_score"], bins=bin_edges)
            data = (
                data.groupby(["prediction_score_bin", "condition_name", "allele_id"])
                .agg(
                    agg_weight=("weight_sum", "sum"),
                )
                .reset_index()
            )
            # For the sake of simplicity, just use bin center for plotting (anyway smoothed for violin plots)
            # (todo: think about more elegant way for visualization if necessary)
            data["prediction_score_bin_center"] = data["prediction_score_bin"].apply(
                lambda x: np.mean([x.left, x.right])
            )

            # Write out results for each allele
            for allele_id in outfile_dict.keys():
                data[data.allele_id == allele_id][
                    ["prediction_score_bin_center", "condition_name", "agg_weight"]
                ].to_csv(outfile_dict[allele_id], mode="a", sep="\t", index=False, header=print_header)
            print_header = False

        print("Done!", flush=True)

    finally:
        for outfile in outfile_dict.values():
            outfile.close()


if __name__ == "__main__":
    sys.exit(main())
