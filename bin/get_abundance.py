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
import gzip
import csv
import io
import argparse

from Bio import SeqIO

import sys


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', "--proteins", required=True, metavar='FILE', help="Input file containing: protein_tmp_id, protein_sequence.")
    parser.add_argument('-d', "--depths", required=True, metavar='FILE', type=argparse.FileType('r'), help="File containing contig depths.")
    parser.add_argument('-m', "--microbiome_id", required=True, help="Corresponding microbiome_id.")
    parser.add_argument('-o', "--output", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: protein_tmp_id, protein_weight, microbiome_id.")
    return parser.parse_args(args)

def get_prot_depth(proteinId, dict_contig_depths):
    # get contig id the protein is originating from ("contigid_*")
    contig_id = proteinId.rsplit('_', 1)[0]
    return dict_contig_depths[contig_id]


def main(args=None):
    args = parse_args(args)

    with gzip.open(args.proteins, "rt") as handle:
        proteinIds = [ row[0] for row in csv.reader(handle, delimiter='\t') ]

    dict_contig_depths = { row[0]:row[1] for row in csv.reader(args.depths, delimiter='\t') }
    print(dict_contig_depths)

    dict_protein_depths = { id:get_prot_depth(id, dict_contig_depths) for id in proteinIds }
    print("protein_tmp_id", "protein_weight", "microbiome_id", sep='\t', file=args.output)
    for prot in dict_protein_depths:
        print(prot, dict_protein_depths[prot], args.microbiome_id, sep='\t', file=args.output, flush=True)


if __name__ == "__main__":
    sys.exit(main())