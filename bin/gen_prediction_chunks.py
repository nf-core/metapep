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

def write_chunks(data):
    """Takes data in form of a table of peptide_id, peptide_sequence and
    identical allele_name values. The data is partitioned into chunks and
    written into individual output files, prepended with a comment line (#)
    indicating the allele name."""
    global cur_chunk
    for start in range(0, len(data), args.max_chunk_size):
        with open(os.path.join(args.outdir, "peptides_" + str(cur_chunk).rjust(5,"0") + ".txt"), 'w') as outfile:
            print(f"#{data.iloc[0].allele_name}#{data.iloc[0].allele_id}", file = outfile)
            data.iloc[start:start+args.max_chunk_size][["peptide_id", "peptide_sequence"]].to_csv(outfile, sep='\t', index=False)
            cur_chunk = cur_chunk + 1

####################################################################################################

try:
    # Parse command line arguments
    args = parse_args()

    # Read input files
    peptides                  = pd.read_csv(args.peptides, sep='\t')
    protein_peptide_occs      = pd.read_csv(args.protein_peptide_occ, sep='\t')
    microbiome_protein_occs   = pd.read_csv(args.microbiome_protein_occ, sep='\t')
    conditions                = pd.read_csv(args.conditions, sep='\t')
    condition_allele_map      = pd.read_csv(args.condition_allele_map, sep='\t')
    alleles                   = pd.read_csv(args.alleles, sep='\t')

    # Create output directory if it doesn't exist
    if os.path.exists(args.outdir) and not os.path.isdir(args.outdir):
        print("ERROR - The target path is not a directory", file = sys.stderr)
        sys.exit(2)
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Identify which predictions have to be computed
    to_predict = peptides\
            .merge(protein_peptide_occs)\
            .merge(microbiome_protein_occs)\
            .merge(conditions)\
            .merge(condition_allele_map)\
            .merge(alleles)[['peptide_id','peptide_sequence', 'allele_id', 'allele_name']]\
            .drop_duplicates()

    # Write the necessary predictions into chunks of peptide lists
    cur_chunk = 0
    to_predict.groupby("allele_id").apply(write_chunks)

    # We're happy if we got here
    print(f"All done. Written {len(to_predict)} peptide prediction requests into {cur_chunk} chunks.")
    sys.exit(0)
except KeyboardInterrupt:
    print("\nUser aborted.", file = sys.stderr)
    sys.exit(1)

