#!/usr/bin/env python

import os
import sys
import errno
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Reformat nf-core/metapep samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-i', "--input", required=True, metavar='FILE', type=argparse.FileType('r'), help="Input samplesheet file containing: condition, type, microbiome_path, alleles, weights_path.")
    parser.add_argument('-m', "--microbiomes", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: microbiome_id, microbiome_path, microbiome_type, weights_path.")
    parser.add_argument('-c', "--conditions", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, condition_name, microbiome_id.")
    parser.add_argument('-a', "--alleles", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: allele_id, allele_name.")
    parser.add_argument('-ca', "--conditions_alleles", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, allele_id.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# TODO nf-core: Update the check_samplesheet function
# TODO for 'proteins' type
# - allow weight input for type 'proteins' as well! (for now use equal weight ?)

def check_samplesheet(args):
    """
    This function checks that the samplesheet follows the following structure:

    condition,type,microbiome_path,alleles,weights_path

    For an example see:
    TODO
    """
    input_table = pd.read_csv(args.input)

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

    # microbiome_id - microbiome_path - microbiome_type
    microbiomes = input_table[["microbiome_path", "type", "weights_path"]].drop_duplicates().rename({"type":"microbiome_type"}, axis=1)
    microbiomes["microbiome_id"] = range(len(microbiomes))

    if len(microbiomes) != len(microbiomes["microbiome_path"].drop_duplicates()):
        sys.exit("Conflicting types or weights were specified for the same microbiome path!")

    microbiomes[["microbiome_id", "microbiome_path", "microbiome_type", "weights_path"]].to_csv(args.microbiomes, sep="\t", index=False)

    # condition id - condition name - microbiome id
    conditions = input_table.merge(microbiomes)[["condition", "microbiome_id"]].rename({"condition":"condition_name"}, axis=1)    # conditions unique (checked in nextflow)
    conditions["condition_id"] = range(len(conditions))

    conditions[["condition_id", "condition_name", "microbiome_id"]].to_csv(args.conditions, sep="\t", index=False)

    # allele id - allele name
    unique_alleles = { allele for allele_list in input_table["alleles"] for allele in allele_list.split(' ') }

    alleles = pd.DataFrame({"allele_name":list(unique_alleles)})
    alleles["allele_id"] = range(len(alleles))
    alleles[["allele_id", "allele_name"]].to_csv(args.alleles, sep="\t", index=False)

    # condition id - allele id
    conditions_alleles = pd.DataFrame([ (row["condition"], allele_name) for _, row in input_table.iterrows() for allele_name in row["alleles"].split(' ') ], columns = ["condition_name", "allele_name"])
    conditions_alleles = conditions_alleles.merge(conditions).merge(alleles)[["condition_id", "allele_id"]]
    conditions_alleles.to_csv(args.conditions_alleles, sep="\t", index=False)

    print("Done!")

    # sample_mapping_dict = {}
    # with open(file_in, "r") as fin:

    #     ## Check header
    #     MIN_COLS = 2
    #     # TODO nf-core: Update the column names for the input samplesheet
    #     HEADER = ["sample", "fastq_1", "fastq_2"]
    #     header = [x.strip('"') for x in fin.readline().strip().split(",")]
    #     if header[: len(HEADER)] != HEADER:
    #         print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
    #         sys.exit(1)

    #     ## Check sample entries
    #     for line in fin:
    #         lspl = [x.strip().strip('"') for x in line.strip().split(",")]

    #         # Check valid number of columns per row
    #         if len(lspl) < len(HEADER):
    #             print_error(
    #                 "Invalid number of columns (minimum = {})!".format(len(HEADER)),
    #                 "Line",
    #                 line,
    #             )
    #         num_cols = len([x for x in lspl if x])
    #         if num_cols < MIN_COLS:
    #             print_error(
    #                 "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
    #                 "Line",
    #                 line,
    #             )

    #         ## Check sample name entries
    #         sample, fastq_1, fastq_2 = lspl[: len(HEADER)]
    #         sample = sample.replace(" ", "_")
    #         if not sample:
    #             print_error("Sample entry has not been specified!", "Line", line)

    #         ## Check FastQ file extension
    #         for fastq in [fastq_1, fastq_2]:
    #             if fastq:
    #                 if fastq.find(" ") != -1:
    #                     print_error("FastQ file contains spaces!", "Line", line)
    #                 if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
    #                     print_error(
    #                         "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
    #                         "Line",
    #                         line,
    #                     )

    #         ## Auto-detect paired-end/single-end
    #         sample_info = []  ## [single_end, fastq_1, fastq_2]
    #         if sample and fastq_1 and fastq_2:  ## Paired-end short reads
    #             sample_info = ["0", fastq_1, fastq_2]
    #         elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
    #             sample_info = ["1", fastq_1, fastq_2]
    #         else:
    #             print_error("Invalid combination of columns provided!", "Line", line)

    #         ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2 ] }
    #         if sample not in sample_mapping_dict:
    #             sample_mapping_dict[sample] = [sample_info]
    #         else:
    #             if sample_info in sample_mapping_dict[sample]:
    #                 print_error("Samplesheet contains duplicate rows!", "Line", line)
    #             else:
    #                 sample_mapping_dict[sample].append(sample_info)

    # ## Write validated samplesheet with appropriate columns
    # if len(sample_mapping_dict) > 0:
    #     out_dir = os.path.dirname(file_out)
    #     make_dir(out_dir)
    #     with open(file_out, "w") as fout:
    #         fout.write(",".join(["sample", "single_end", "fastq_1", "fastq_2"]) + "\n")
    #         for sample in sorted(sample_mapping_dict.keys()):

    #             ## Check that multiple runs of the same sample are of the same datatype
    #             if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
    #                 print_error("Multiple runs of a sample must be of the same datatype!", "Sample: {}".format(sample))

    #             for idx, val in enumerate(sample_mapping_dict[sample]):
    #                 fout.write(",".join(["{}_T{}".format(sample, idx + 1)] + val) + "\n")
    # else:
    #     print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args)


if __name__ == "__main__":
    sys.exit(main())
