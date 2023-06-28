#!/usr/bin/env python

"""
Provides a command line tool to validate and transform the samplesheet and create the data tables needed for nf-core/metapep.
"""


import os
import sys
import errno
import argparse
import pandas as pd
from epytope.Core import Allele
from epytope.EpitopePrediction import EpitopePredictorFactory
import epytope


def parse_args(args=None):
    Description = "Reformat nf-core/metapep samplesheet file, check its contents and create the data tables."
    Epilog = "Example usage: python check_samplesheet_create_tables.py -i <FILE_IN> -m <MICROBIOMES_OUT> -c <CONDITIONS_OUT> -a <ALLELES_OUT> -ca <CONDITION_ALLELES_OUT> -pm <prediction_method> -pmv <prediction_method_version>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="FILE",
        type=argparse.FileType("r"),
        help="Input samplesheet file containing: condition, type, microbiome_path, alleles, weights_path.",
    )
    parser.add_argument(
        "-m",
        "--microbiomes",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: microbiome_id, microbiome_path, microbiome_type, weights_path.",
    )
    parser.add_argument(
        "-c",
        "--conditions",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: condition_id, condition_name, microbiome_id.",
    )
    parser.add_argument(
        "-a",
        "--alleles",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: allele_id, allele_name.",
    )
    parser.add_argument(
        "-ca",
        "--conditions_alleles",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: condition_id, allele_id.",
    )
    parser.add_argument(
        "-pm",
        "--prediction_method",
        required=True,
        metavar="STRING",
        type=str,
        help="Chosen method for epitope prediction",
    )
    parser.add_argument(
        "-pmv",
        "--pred_method_version",
        required=True,
        metavar="STRING",
        type=str,
        help="Chosen version of epitope prediction method",
    )
    parser.add_argument(
        "-pl",
        "--peptide_lengths",
        required=True,
        metavar="TUPLE",
        nargs="+",
        help="Peptide lengths as given in parameters (min max)",
    )
    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# TODO for 'proteins' type
# - allow weight input for type 'proteins' as well! (for now use equal weight ?)


def process_samplesheet(args):
    """
    Check that the tabular samplesheet has the structure expected by nf-core/metapep and create the data tables.

    Header structure:
    condition,type,microbiome_path,alleles,weights_path

    For an example see:
    https://github.com/nf-core/metapep/raw/dev/assets/samplesheet.csv
    """
    input_table = pd.read_csv(args.input)
    input_table_cp = input_table.copy()

    # check if microbiome_path file extensions are valid
    for type, fname in zip(input_table["type"], input_table["microbiome_path"]):
        if type == "proteins":
            print_error(
                "Invalid type '"
                + type
                + "' specified in "
                + args.input.name
                + ". Type 'proteins' is not yet supported! Valid types are 'taxa', 'bins' and 'assembly'."
            )
        if type not in ["taxa", "assembly", "bins"]:
            print_error(
                "Invalid type '"
                + type
                + "' specified in "
                + args.input.name
                + ". Valid types are 'taxa', 'bins' and 'assembly'."
            )
        if type == "taxa" and not fname.lower().endswith((".tsv")):
            print_error(
                "In "
                + args.input.name
                + " specified file "
                + fname
                + " of type 'taxa' has invalid file extension. Valid extensions is '.tsv'."
            )
        if type == "proteins" and not fname.lower().endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            print_error(
                "In "
                + args.input.name
                + " specified file "
                + fname
                + " of type 'proteins' has invalid file extension. Valid extensions are '.fa', '.fa.gz', '.fasta' and"
                " '.fasta.gz'."
            )
        if type == "assembly" and not fname.lower().endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            print_error(
                "In "
                + args.input.name
                + " specified file "
                + fname
                + " of type 'assembly' has invalid file extension. Valid extensions are '.fa', '.fa.gz', '.fasta' and"
                " '.fasta.gz'."
            )

    # check if microbiome_path files exist
    # for fname in input_table["microbiome_path"]:
    #     if not os.path.isfile(fname):
    #         sys.exit("In " + args.input.name + " specified file " + fname + " does not exist!")
    # NOTE not possible for urls, will be checked afterwards during channel creation for microbiome files

    # check if condition names unique
    if len(input_table["condition"]) != len(input_table["condition"].drop_duplicates()):
        sys.exit("Input file " + args.input.name + " contains duplicated conditions! Please specify unique conditions.")

    # check if weight_path is valid
    for type, weights_path in zip(input_table["type"], input_table["weights_path"]):
        if not (type == "assembly" or type == "bins") and not pd.isnull(weights_path):
            sys.exit(
                "Input file "
                + args.input.name
                + " contains 'weights_path' specified for type '"
                + type
                + "'! Currently input weights are only supported for type 'assembly' or 'bins."
            )
        if not pd.isnull(weights_path) and not weights_path.lower().endswith(".tsv"):
            sys.exit(
                "In "
                + args.input.name
                + " specified 'weights_path' "
                + weights_path
                + " has invalid file extension. The extension must be '.tsv'."
            )

    # microbiome_id - microbiome_path - microbiome_type
    microbiomes = (
        input_table[["microbiome_path", "type", "weights_path"]]
        .drop_duplicates()
        .rename({"type": "microbiome_type"}, axis=1)
    )
    microbiomes["microbiome_id"] = range(len(microbiomes))

    # Create bare id for each microbiome path to reduce redundancy in protein generation
    mcrb_uni = {}
    x = 0
    for path in microbiomes["microbiome_path"]:
        if path not in mcrb_uni:
            mcrb_uni[path] = x
            x += 1
        else:
            continue

    microbiomes["microbiome_bare_id"] = [mcrb_uni[path] for path in microbiomes["microbiome_path"]]

    microbiomes[["microbiome_id", "microbiome_path", "microbiome_type", "weights_path", "microbiome_bare_id"]].to_csv(
        args.microbiomes, sep="\t", index=False
    )

    # condition id - condition name - microbiome id
    conditions = input_table.merge(microbiomes)[["condition", "microbiome_id"]].rename(
        {"condition": "condition_name"}, axis=1
    )  # conditions unique (checked in nextflow)
    conditions["condition_id"] = range(len(conditions))

    conditions[["condition_id", "condition_name", "microbiome_id"]].to_csv(args.conditions, sep="\t", index=False)

    # allele id - allele name
    unique_alleles = {allele for allele_list in input_table["alleles"] for allele in allele_list.split(" ")}

    # Check if alleles are supported by chosen predictor and check if tool generally supports chosen lengths

    # TODO Parameter for prediction method version
    predictor = EpitopePredictorFactory(args.prediction_method, version=args.pred_method_version)
    for allele in unique_alleles:
        if not Allele(allele) in predictor.supportedAlleles:
            sys.exit(
                "\n\n\n\nThe chosen allele: "
                + allele
                + " is not available for the chosen prediction method: "
                + args.prediction_method
                + ":"
                + args.pred_method_version
                + "\n\nFurther information on which allele is supported for "
                + "which prediction method can be found when running: "
                + "'nextflow run metapep -profile <YOURPROFILE> --outdir <OUTDIR> --show_supported_models"
            )

    # Check if peptide lengths are listed as supported for the allele model by epytope:
    # Note in some cases for individual alleles certain lengths are not supported, which is not captured here
    peptide_lengths = range(int(args.peptide_lengths[0]), int(args.peptide_lengths[1]) + 1)
    checked_pep_lens = set(peptide_lengths)
    for pep_len in peptide_lengths:
        if pep_len in predictor.supportedLength:
            checked_pep_lens.remove(pep_len)
        if not pep_len in predictor.supportedLength:
            sys.exit(
                "\n\n\n\nThe chosen lengths: "
                + ", ".join([str(i) for i in checked_pep_lens])
                + " are not available for the chosen prediction method: "
                + str(args.prediction_method)
                + ":"
                + str(args.pred_method_version)
                + "\n\nFurther information on which peptide lengths are supported for "
                + "which prediction method can be found when running: "
                + "'nextflow run metapep -profile <YOURPROFILE> --outdir <OUTDIR> --show_supported_models"
            )

    alleles = pd.DataFrame({"allele_name": list(unique_alleles)})
    alleles["allele_id"] = range(len(alleles))
    alleles[["allele_id", "allele_name"]].to_csv(args.alleles, sep="\t", index=False)

    # condition id - allele id
    conditions_alleles = pd.DataFrame(
        [
            (row["condition"], allele_name)
            for _, row in input_table.iterrows()
            for allele_name in row["alleles"].split(" ")
        ],
        columns=["condition_name", "allele_name"],
    )
    conditions_alleles = conditions_alleles.merge(conditions).merge(alleles)[["condition_id", "allele_id"]]
    conditions_alleles.to_csv(args.conditions_alleles, sep="\t", index=False)

    input_table_cp.to_csv("samplesheet.valid.csv", index=False)
    print("Done!")


def main(args=None):
    args = parse_args(args)
    process_samplesheet(args)


if __name__ == "__main__":
    sys.exit(main())
