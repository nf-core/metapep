#!/usr/bin/env python

"""Provide a command line tool to validate and transform tabular samplesheets."""


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


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="fastq_1",
        second_col="fastq_2",
        single_col="single_end",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._single_col = single_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_pair(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0, "Sample input is required."
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        assert len(row[self._first_col]) > 0, "At least the first FASTQ file is required."
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._second_col]) > 0:
            self._validate_fastq_format(row[self._second_col])

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            row[self._single_col] = False
            assert (
                Path(row[self._first_col]).suffixes == Path(row[self._second_col]).suffixes
            ), "FASTQ pairs must have the same file extensions."
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FORMATS), (
            f"The FASTQ file has an unrecognized extension: {filename}\n"
            f"It should be one of: {', '.join(self.VALID_FORMATS)}"
        )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename the sample if more than one sample,
        FASTQ file combination exists.

        """
        assert len(self._seen) == len(self.modified), "The pair of sample name and FASTQ must be unique."
        if len({pair[0] for pair in self._seen}) < len(self._seen):
            counts = Counter(pair[0] for pair in self._seen)
            seen = Counter()
            for row in self.modified:
                sample = row[self._sample_col]
                seen[sample] += 1
                if counts[sample] > 1:
                    row[self._sample_col] = f"{sample}_T{seen[sample]}"


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = handle.read(2048)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical(f"The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    handle.seek(0)
    return dialect

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

    input_table_cp.to_csv("samplesheet.valid.csv", index=False)
    print("Done!")

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args)


if __name__ == "__main__":
    sys.exit(main())
