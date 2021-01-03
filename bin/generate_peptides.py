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
import gzip
import csv
import pandas as pd
#import tqdm
import io
import argparse
import time

from Bio import SeqIO
from pprint import pprint
from datetime import datetime
from collections import Counter

import sys


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', "--proteins", required=True, metavar='FILE', help="Compressed TSV file containing: protein_id, protein_sequence.")
    parser.add_argument('-mn', "--min_len", required=True, metavar='N', type=int, help="Min. peptide length.")
    parser.add_argument('-mx', "--max_len", required=True, metavar='N', type=int, help="Max. peptide length.")
    parser.add_argument('-pp', '--peptides', required=True, metavar='FILE', help='Output file containing peptides.') # use str type to allow compression of output
    parser.add_argument('-l', '--prot_lengths', required=True, metavar='FILE', type=argparse.FileType('w'), help='Output file containing protein lengths.')
    return parser.parse_args(args)

def gen_peptides(prot_seq, k):
    return [ prot_seq[i:(i+k)] for i in range(len(prot_seq)-k) ]


def main(args=None):
    args = parse_args(args)

    with gzip.open(args.proteins, "rt") as handle:
        # get protein id, sequence, len
        protid_protseq_protlen = pd.DataFrame(
            [ (row[0], row[1], len(row[1])) for row in csv.reader(handle, delimiter='\t') ],
            columns = ['protein','sequence', 'length']
            )
    print("# proteins: ", len(protid_protseq_protlen))

    # write out protein lengths
    protid_protseq_protlen[['protein', 'length']].to_csv(args.prot_lengths, sep="\t", index=False)

    ####################
    # generate peptides
    # write header
    out = pd.DataFrame([], columns = ['sequence','id','proteins','counts'])
    out.to_csv(args.peptides, sep="\t", index=False, compression='gzip')

    # for each k
    for k in range(args.min_len, args.max_len + 1):
        print("Generate peptides of length ", k, " ...", flush=True)
        # for each protein generate all peptides of length k
        prot_peptides = pd.DataFrame(
            [ (it.protein, pep) for it in protid_protseq_protlen.itertuples() for pep in gen_peptides(it.sequence, k) ],
            columns = ['protein','peptides']
            )

        print("format results ...", flush=True)
        # count occurences of one peptide in one protein
        prot_peptides = prot_peptides.groupby(['protein','peptides']).size().reset_index(name='count')

        # aggregate for each peptide: -> pep_id     pep_seq    'prot1','prot2',..    3,6,0,..
        results = prot_peptides.groupby('peptides').agg(list)
        results = results.reset_index()
        results = results.assign(id=["P_k" + str(k) + "_" + str(id) for id in results.index])
        # rename column names
        results.columns = ['sequence', 'proteins', 'counts', 'id']
        # convert to string and then joint to get rid of brackets and quotes
        results["proteins"] = results["proteins"].str.join(",") 
        results["counts"] = results["counts"].apply(lambda x : ','.join([ str(e) for e in  x]))
        print(results.head(5), flush=True)
        results[['sequence','id','proteins','counts']].to_csv(args.peptides, mode='a', sep="\t", index=False, header=False, compression='gzip')
        print("# peptides of length ", k, ", (non-unique across proteins): ", len(results))

    print("Done!", flush=True)


if __name__ == "__main__":
    sys.exit(main())
