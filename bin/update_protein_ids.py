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
# NOTE
# entrez proteins of all microbiome input files already within one file (proteins.entrez.tsv.gz)
# proteins from different assemblies have non-unique names originating from enumerated contigs, should get new ids assigned separately for each file (proteins.pred_${microbiome_id}.tsv.gz)
# proteins from 'proteins' input type: not known if unique or not, handle separately for now (in case of unique ids this causes unnecessary redundancy; could add parameter for this in future)

import sys
import gzip
import csv
import io
import argparse

from Bio import SeqIO

import sys


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-ip', "--in_proteins", required=True, nargs="+", metavar='FILE', help="Input file containing: protein_tmp_id, protein_sequence.")
    parser.add_argument('-ipm', "--in_proteins_microbiomes", required=True, nargs="+", metavar='FILE', type=argparse.FileType('r'), help="Input file containing: protein_tmp_id, protein_weight, microbiome_id.")
    parser.add_argument('-op', "--out_proteins", required=True, metavar='FILE', help="output file containing: protein_id, protein_sequence.")
    parser.add_argument('-opm', "--out_proteins_microbiomes", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: protein_id, protein_weight, microbiome_id.")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    id = 0
    # prepare outputs
    with gzip.open(args.out_proteins, "wt") as out_handle:
        print("protein_id", "protein_sequence", sep='\t', file=out_handle)
        print("protein_id", "protein_weight", "microbiome_id", sep='\t', file=args.out_proteins_microbiomes)

        # for each pair in input file lists
        for proteins, proteins_microbiomes in zip(args.in_proteins, args.in_proteins_microbiomes):
            dict_oldId_newId = {}
            # write new 'proteins.tsv.gz'
            with gzip.open(proteins, "rt") as in_handle:
                reader = csv.reader(in_handle, delimiter='\t')
                header = next(reader)
                for row in reader:
                    print(id, row[1], sep='\t', file=out_handle, flush=True)
                    if row[0] in dict_oldId_newId:
                        print("ERROR: protein ID " + row[0] + " was already read from " +  proteins + " file!")
                        sys.exit("Protein ID non-unique for proteins file!")
                    dict_oldId_newId[row[0]] = id
                    id += 1
            # write new 'proteins_microbiomes.tsv'
            reader = csv.reader(proteins_microbiomes, delimiter='\t')
            header = next(reader)
            for row in reader:
                if row[0] not in dict_oldId_newId:
                    print(row)
                    print("ERROR: protein ID " + row[0] + " from " + proteins_microbiomes.name + " not in dictionary for " +  proteins + " file!")
                    sys.exit("Protein ID not in dictionary!")
                newId = dict_oldId_newId[row[0]]
                print(newId, row[1], row[2], sep='\t', file=args.out_proteins_microbiomes, flush=True)

    print("Processed " + str(id) + " proteins.")
    print("Done!")

if __name__ == "__main__":
    sys.exit(main())