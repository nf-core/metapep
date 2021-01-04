#!/usr/bin/env python3
####################################################################################################
#
# Author: Sabrina Krakau, Leon Kuchenbecker
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

import sys
import csv
import pandas as pd
import io
import argparse


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input", required=True, metavar='FILE', type=argparse.FileType('r'), help="Input file containing: condition_name, microbiome_type, microbiome_path and allele_name(s).")
    parser.add_argument('-m', "--microbiomes", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: microbiome_id, microbiome_path, microbiome_type.")
    parser.add_argument('-c', "--conditions", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, condition_name, microbiome_id.")
    parser.add_argument('-a', "--alleles", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: allele_id, allele_name.")
    parser.add_argument('-ca', "--conditions_alleles", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, allele_id.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    input_table = pd.read_csv(args.input, sep='\t')

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
    unique_alleles = { allele for allele_list in input_table["alleles"] for allele in allele_list.split(',') }

    alleles = pd.DataFrame({"allele_name":list(unique_alleles)})
    alleles["allele_id"] = range(len(alleles))
    alleles[["allele_id", "allele_name"]].to_csv(args.alleles, sep="\t", index=False)

    # condition id - allele id
    conditions_alleles = pd.DataFrame([ (row["condition"], allele_name) for _, row in input_table.iterrows() for allele_name in row["alleles"].split(',') ], columns = ["condition_name", "allele_name"])
    conditions_alleles = conditions_alleles.merge(conditions).merge(alleles)[["condition_id", "allele_id"]]
    conditions_alleles.to_csv(args.conditions_alleles, sep="\t", index=False)

    print("Done!")


if __name__ == "__main__":
    sys.exit(main())
