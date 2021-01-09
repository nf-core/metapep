#!/usr/bin/env python3
#
# Author: Leon Kuchenbecker <leon.kuchenbecker@uni-tuebingen.de>
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
# epytope_predict.py
#
# This program provides a command line interface to the epitope prediction module of the 'epytope'
# python framework.
####################################################################################################

import argparse
import sys
import os
from tqdm import tqdm

import pandas as pd

####################################################################################################

def parse_args():
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Generate chunks of peptides that are to be predicted against a specified allele.")

    # INPUT FILES
    parser.add_argument("-p"     , "--peptides"                 , help="Path to the peptides input file"                        , type=str   , required=True)
    parser.add_argument("-ppo"   , "--protein-peptide-occ"      , help="Path to the protein peptide occurences input file"      , type=str   , required=True)
    parser.add_argument("-mpo"   , "--microbiome-protein-occ"   , help="Path to the microbiome protein occurences input file"   , type=str   , required=True)
    parser.add_argument("-c"     , "--conditions"               , help="Path to the conditions input file"                      , type=str   , required=True)
    parser.add_argument("-cam"   , "--condition-allele-map"     , help="Path to the condition allele map input file"            , type=str   , required=True)
    parser.add_argument("-a"     , "--alleles"                  , help="Path to the allele input file"                          , type=str   , required=True)

    # OUTPUT FILES
    parser.add_argument("-o"     , "--outdir"                   , help="Path to the output directory"                           , type=str                      , required=True)

    # PARAMETERS
    parser.add_argument("-mc"    , "--max-chunk-size"           , help="Maximum chunk size. Default: 5000"                      , type=int                      , default=5000)

    return parser.parse_args()

def write_chunks(data, pbar=None):
    """Takes data in form of a table of peptide_id, peptide_sequence and
    identical allele_name values. The data is partitioned into chunks and
    written into individual output files, prepended with a comment line (#)
    indicating the allele name."""
    global cur_chunk
    for start in range(0, len(data), args.max_chunk_size):
        with open(os.path.join(args.outdir, "peptides_" + str(cur_chunk).rjust(5,"0") + ".txt"), 'w') as outfile:
            print(f"#{data.iloc[0].allele_name}#{data.iloc[0].allele_id}", file = outfile)
            write = data.iloc[start:start+args.max_chunk_size]
            if pbar:
                pbar.update(len(write))
            write[["peptide_id", "peptide_sequence"]].to_csv(outfile, sep='\t', index=False)
            cur_chunk = cur_chunk + 1

####################################################################################################

try:
    # Parse command line arguments
    args = parse_args()

    # Read input files
    peptides                  = pd.read_csv(args.peptides, sep='\t')
    protein_peptide_occs      = pd.read_csv(args.protein_peptide_occ, sep='\t').drop(columns="count")
    microbiome_protein_occs   = pd.read_csv(args.microbiome_protein_occ, sep='\t').drop(columns="protein_weight")
    conditions                = pd.read_csv(args.conditions, sep='\t').drop(columns="condition_name")
    condition_allele_map      = pd.read_csv(args.condition_allele_map, sep='\t')
    alleles                   = pd.read_csv(args.alleles, sep='\t')

    # Create output directory if it doesn't exist
    if os.path.exists(args.outdir) and not os.path.isdir(args.outdir):
        print("ERROR - The target path is not a directory", file = sys.stderr)
        sys.exit(2)
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("Joining input data...", file = sys.stderr, flush=True, end='')

    # Identify which predictions have to be computed
    to_predict = peptides\
            .merge(protein_peptide_occs)\
            .merge(microbiome_protein_occs)\
            .drop(columns="protein_id")\
            .merge(conditions)\
            .drop(columns="microbiome_id")\
            .merge(condition_allele_map)\
            .drop(columns="condition_id")\
            .merge(alleles)\
            .drop_duplicates()

    print(" done.", file = sys.stderr, flush=True)
    print("Writing output files...", file = sys.stderr, flush=True)

    # Write the necessary predictions into chunks of peptide lists
    cur_chunk = 0
    with tqdm(total=len(to_predict), ascii=True, unit=" peptides") as pbar:
        to_predict.groupby("allele_id").apply(lambda x : write_chunks(x, pbar=pbar))

    # We're happy if we got here
    print(f"All done. Written {len(to_predict)} peptide prediction requests into {cur_chunk} chunks.")
    sys.exit(0)
except KeyboardInterrupt:
    print("\nUser aborted.", file = sys.stderr)
    sys.exit(1)
