#!/usr/bin/env python

"""Provide a command line tool to validate and transform tabular samplesheets."""


import os
import sys
import errno
import argparse
import pandas as pd
import numpy as np


def parse_args(args=None):
    Description = "Reformat nf-core/metapep samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-i', "--input", required=True, metavar='FILE', type=argparse.FileType('r'), help="Input samplesheet file containing: condition, type, microbiome_path, alleles, weights_id.")
    parser.add_argument('-m', "--microbiomes", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: microbiome_id, microbiome_path, microbiome_type, 'conditions', 'alleles', 'weights_ids.")
    parser.add_argument('-c', "--conditions", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, condition_name.")
    parser.add_argument('-a', "--alleles", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: allele_id, allele_name.")
    parser.add_argument('-w', "--weights", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: weights_id, weights_path.")
    parser.add_argument('-cm', "--conditions_microbiomes", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, microbiome_id.")
    parser.add_argument('-ca', "--conditions_alleles", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, allele_id.")
    parser.add_argument('-cw', "--conditions_weights", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, weights_id.")
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

def check_samplesheet(args):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    condition,type,microbiome_path,alleles,weights_path

    For an example see:
    https://github.com/nf-core/metapep/raw/dev/assets/samplesheet.csv
    """
    input_table = pd.read_csv(args.input)
    input_table_cp = input_table.copy()

    # check if microbiome_path file extensions are valid
    for type, fname in zip(input_table["type"], input_table["microbiome_path"]):
        if type == "proteins":
            print_error("Invalid type '" + type + "' specified in " + args.input.name + ". Type 'proteins' is not yet supported! Valid types are 'taxa', 'bins' and 'assembly'.")
        if type not in ["taxa", "assembly", "bins"]:
            print_error("Invalid type '" + type + "' specified in " + args.input.name + ". Valid types are 'taxa', 'bins' and 'assembly'.")
        if type == "taxa" and not fname.lower().endswith(('.txt', '.tsv')):
            print_error("In " + args.input.name + " specified file " + fname + " of type 'taxa' has invalid file extension. Valid extensions are '.txt' and '.tsv'.")
        if type == "proteins" and not fname.lower().endswith(('.fa', '.fa.gz', '.fasta', '.fasta.gz')):
            print_error("In " + args.input.name + " specified file " + fname + " of type 'proteins' has invalid file extension. Valid extensions are '.fa', '.fa.gz', '.fasta' and '.fasta.gz'.")
        if type == "assembly" and not fname.lower().endswith(('.fa', '.fa.gz', '.fasta', '.fasta.gz')):
            print_error("In " + args.input.name + " specified file " + fname + " of type 'assembly' has invalid file extension. Valid extensions are '.fa', '.fa.gz', '.fasta' and '.fasta.gz'.")

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
        if not type == 'assembly' and not pd.isnull(weights_path):
            sys.exit("Input file " + args.input.name + " contains 'weights_path' specified for type '" + type + "'! Currently input weights are only supported for type 'assembly'.")
        if not pd.isnull(weights_path) and not weights_path.lower().endswith('.tsv'):
            sys.exit("In " + args.input.name + " specified 'weights_path' " + weights_path + " has invalid file extension. The extension must be '.tsv'.")

    db_table = input_table.copy()

    # weights id - weights path
    weights_group = db_table.groupby(["weights_path"], sort=False)

    db_table["weights_id"] = pd.Series([n if n >= 0 else np.nan for n in weights_group.ngroup()], dtype="Int64")
    weights_group.agg('first').reset_index()[["weights_id", "weights_path"]].to_csv(args.weights, sep="\t", index=False)

    # microbiome_id - microbiome_path - microbiome_type
    microbiome_group = db_table.groupby(["microbiome_path", "type"], sort=False)

    if microbiome_group.ngroups != len(db_table["microbiome_path"].drop_duplicates()):
        sys.exit("Conflicting types were specified for the same microbiome path!")

    db_table["microbiome_id"] = microbiome_group.ngroup()
    microbiome_group \
        .agg({'microbiome_id': 'first', 'condition': lambda x: ';'.join(list(x)), 'alleles': lambda x: ';'.join(list(x)), 'weights_id': lambda x: ';'.join([str(y) if not pd.isna(y) else "" for y in list(x)])}) \
        .reset_index() \
        .rename({"type":"microbiome_type", "condition":"conditions", "weights_id":"weights_ids"}, axis=1)[['microbiome_id', 'microbiome_path', 'microbiome_type', 'conditions', 'alleles', 'weights_ids']] \
        .to_csv(args.microbiomes, sep="\t", index=False)

    # condition id - condition name - microbiome id
    db_table = db_table.reset_index().rename({'index':'condition_id', 'condition':'condition_name'}, axis=1)
    db_table[["condition_id", "condition_name"]].to_csv(args.conditions, sep="\t", index=False)

    # condition id - microbiome id
    db_table[["condition_id", "microbiome_id"]].to_csv(args.conditions_microbiomes, sep="\t", index=False)

    # allele id - allele name
    db_table["alleles"] = db_table.apply(lambda row: list(row.alleles.split(' ')), axis=1)
    alleles_exploded = db_table[["condition_id", "alleles"]].explode("alleles")
    allele_group = alleles_exploded.groupby(["alleles"], sort=False)
    alleles_exploded["allele_id"] = allele_group.ngroup()
    allele_group.agg('first').reset_index().rename({"alleles":"allele_name"}, axis=1)[["allele_id", "allele_name"]].to_csv(args.alleles, sep="\t", index=False)

    # condition id - allele id
    alleles_exploded[["condition_id", "allele_id"]].to_csv(args.conditions_alleles, sep="\t", index=False)



    # condition id - weights id
    db_table[["condition_id", "weights_id"]].dropna().to_csv(args.conditions_weights, sep="\t", index=False)

    input_table_cp.to_csv("samplesheet.valid.csv", index=False)
    print("Done!")

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args)


if __name__ == "__main__":
    sys.exit(main())
