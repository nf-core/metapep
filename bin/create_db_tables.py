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
    parser.add_argument('-i', "--input", required=True, metavar='FILE', type=argparse.FileType('r'), help="Input file containing: condition_name, microbiome_type, microbiome_path and allele_name(s).")
    parser.add_argument('-m', "--microbiomes", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: microbiome_id, microbiome_path, microbiome_type.")
    parser.add_argument('-c', "--conditions", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, condition_name, microbiome_id.")
    parser.add_argument('-a', "--alleles", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: allele_id, allele_name.")
    parser.add_argument('-ca', "--condition_allele", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: condition_id, allele_id.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    input_table = list(csv.reader(args.input, delimiter='\t'))

    # microbiome_id - microbiome_path - microbiome_type
    dict_microbiome_path_type = {}
    for row in input_table:
        if row[2] in dict_microbiome_path_type:
            if dict_microbiome_path_type[row[2]] != row[1]:
                print("ERROR: for microbiome path ", row[2], " different types '", dict_microbiome_path_type[row[2]], "' and '", row[1], "' were specified.", sep='')
                sys.exit("Conflicting types were specified for microbiome file path!")
        else:
            dict_microbiome_path_type[row[2]] = row[1]

    dict_microbiomes_path_id = {}
    print("microbiome_id", "microbiome_path", "microbiome_type", sep='\t', file=args.microbiomes)
    for microbiome_id, microbiome_path in enumerate(dict_microbiome_path_type):
        microbiome_type = dict_microbiome_path_type[microbiome_path]
        print(microbiome_id, microbiome_path, microbiome_type, sep='\t', file=args.microbiomes, flush=True)
        dict_microbiomes_path_id[microbiome_path] = microbiome_id

    # condition id - condition name - microbiome id
    dict_condition_path = { row[0]:row[2] for row in input_table }
    print("condition_id", "condition_name", "microbiome_id", sep='\t', file=args.conditions)
    for condition_id, condition_name in enumerate(dict_condition_path):
        microbiome_id = dict_microbiomes_path_id[dict_condition_path[condition_name]]
        print(condition_id, condition_name, microbiome_id, sep='\t', file=args.conditions, flush=True)

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
