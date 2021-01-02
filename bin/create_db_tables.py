#!/usr/bin/env python3
####################################################################################################
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

import sys
import csv
import io
import argparse


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input", required=True, metavar='FILE', type=argparse.FileType('r'), help="Input file containing: condition name, input data type, input data path and alleles.")
    parser.add_argument('-m', "--microbiomes", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: microbiome id, path.")
    parser.add_argument('-cm', "--condition_microbiome", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition id, microbiome id.")
    parser.add_argument('-a', "--alleles", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: allele id, allele name.")
    parser.add_argument('-ca', "--condition_allele", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition id, allele.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    input_table = list(csv.reader(args.input, delimiter='\t'))

    # microbiome_id - path
    microbiome_paths = set([ row[2] for row in input_table ])
    dict_microbiomes = {}
    print("microbiome_id", "path", sep='\t', file=args.microbiomes)
    for microbiome_id, microbiome_path in enumerate(microbiome_paths):
        dict_microbiomes[microbiome_path] = microbiome_id
        print(microbiome_id, microbiome_path, sep='\t', file=args.microbiomes, flush=True)

    # condition id - condition name - microbiome id
    dict_condition_path = { row[0]:row[2] for row in input_table }
    print("condition_id", "condition_name", "microbiome_id", sep='\t', file=args.condition_microbiome)
    for condition_id, condition_name in enumerate(dict_condition_path):
        microbiome_id = dict_microbiomes[dict_condition_path[condition_name]]
        print(condition_id, condition_name, microbiome_id, sep='\t', file=args.condition_microbiome, flush=True)

    # allele id - allele name
    allele_names = set([ allele for row in input_table for allele in row[3].split(',') ])
    dict_alleles = {}
    print("allele_id", "allele_name", sep='\t', file=args.alleles)
    for allele_id, allele_name in enumerate(allele_names):
        dict_alleles[allele_name] = allele_id
        print(allele_id, allele_name, sep='\t', file=args.alleles, flush=True)

    # condition id - allele id
    dict_condition_alleles = { row[0]:row[3].split(',') for row in input_table }
    print("condition_id", "allele_id", sep='\t', file=args.condition_allele)
    for condition_id, condition_name in enumerate(dict_condition_alleles):
        for allele_name in dict_condition_alleles[condition_name]:
            allele_id = dict_alleles[allele_name]
            print(condition_id, allele_id, sep='\t', file=args.condition_allele, flush=True)

    print("Done!")


if __name__ == "__main__":
    sys.exit(main())
