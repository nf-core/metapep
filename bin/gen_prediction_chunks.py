#!/usr/bin/env python3
#
# Author: Leon Kuchenbecker <leon.kuchenbecker@uni-tuebingen.de>, Sabrina Krakau <sabrina.krakau@uni-tuebingen.de>
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
import math

####################################################################################################

def parse_args():
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Generate chunks of peptides that are to be predicted against a specified allele.")

    # INPUT FILES
    parser.add_argument("-p"     , "--peptides"                 , help="Path to the peptides input file"                                                  , type=str, required=True)
    parser.add_argument("-ppo"   , "--protein-peptide-occ"      , help="Path to the protein peptide occurences input file"                                , type=str, required=True)
    parser.add_argument("-epo"   , "--entities-proteins-occ"    , help="Path to the entity protein occurences input file"                                 , type=str, required=True)
    parser.add_argument("-meo"   , "--microbiomes-entities-occ" , help="Path to the microbiome entity occurences input file"                              , type=str, required=True)
    parser.add_argument("-c"     , "--conditions"               , help="Path to the conditions input file"                                                , type=str, required=True)
    parser.add_argument("-cam"   , "--condition-allele-map"     , help="Path to the condition allele map input file"                                      , type=str, required=True)
    parser.add_argument("-a"     , "--alleles"                  , help="Path to the allele input file"                                                    , type=str, required=True)

    # OUTPUT FILES
    parser.add_argument("-o"     , "--outdir"                   , help="Path to the output directory"                                                     , type=str, required=True)

    # PARAMETERS
    parser.add_argument("-mc"    , "--max-chunk-size"           , help="Maximum chunk size used for final output files. Default: 5000"                    , type=int               , default=5000)
    parser.add_argument("-pc"    , "--proc-chunk-size"          , help="Chunk size used for internal processing to limit memory usage. Default: 500000"   , type=int               , default=500000)
    # parser.add_argument("-sn"    , "--sample_n"                 , help="Number of peptides to subsample for each condition"     , type=int)
    # NOTE subsampling option currently not working in pipeline!
    return parser.parse_args()

def write_chunks(data, alleles, remainder=False, pbar=None):
    """Takes data in form of a table of peptide_id, peptide_sequence and
    identical allele_name values. The data is partitioned into chunks and
    written into individual output files, prepended with a comment line (#)
    indicating the allele name."""
    global cur_chunk

    if remainder and len(data) > args.max_chunk_size:
        print("ERROR: Something went wrong!", file = sys.stderr)
        sys.exit(1)

    allele_name = alleles[alleles["allele_id"] == data.iloc[0].allele_id]["allele_name"].iloc[0]
    written = pd.Index([])
    for start in range(0, len(data), args.max_chunk_size):
        # if not handling remainder: only write out full chunks here
        if remainder or len(data) - start >= args.max_chunk_size:
            with open(os.path.join(args.outdir, "peptides_" + str(cur_chunk).rjust(5,"0") + ".txt"), 'w') as outfile:
                print(f"#{allele_name}#{data.iloc[0].allele_id}", file = outfile)
                write = data.iloc[start:start+args.max_chunk_size]
                written = written.append(data.index[start:start+args.max_chunk_size])
                if pbar:
                    pbar.update(len(write))
                write[["peptide_id", "peptide_sequence"]].to_csv(outfile, sep='\t', index=False)
                cur_chunk = cur_chunk + 1

    # delete chunks that were written out already
    data.drop(written, inplace=True)
    return data

####################################################################################################

try:
    # Parse command line arguments
    args = parse_args()

    # Read input files
    # NOTE try out if datatable package can be used and would be faster or more memory efficient
    # peptides                  = pd.read_csv(args.peptides, sep='\t', index_col="peptide_id").sort_index() -> read in chunk-wise
    protein_peptide_occs      = pd.read_csv(args.protein_peptide_occ, usecols=['protein_id', 'peptide_id'], sep='\t').set_index('peptide_id').sort_index()  # NOTE could this be handled more efficiently somehow (easily)?
    entities_proteins_occs    = pd.read_csv(args.entities_proteins_occ, sep='\t')
    microbiomes_entities_occs = pd.read_csv(args.microbiomes_entities_occ, usecols=['microbiome_id', 'entity_id'], sep='\t')
    conditions                = pd.read_csv(args.conditions, usecols=['condition_id', 'microbiome_id'], sep='\t')
    condition_allele_map      = pd.read_csv(args.condition_allele_map, sep='\t')
    alleles                   = pd.read_csv(args.alleles, sep='\t')

    # Print info including memory usage for input DataFrames
    # (Pandas df ~ 3x the size of input data)
    print_mem = 'deep'    # 'deep' (extra computational costs) or None
    print("\nInfo: protein_peptide_occs", flush=True)
    protein_peptide_occs.info(verbose = False, memory_usage=print_mem)
    print("\nInfo: entities_proteins_occs", flush=True)
    entities_proteins_occs.info(verbose = False, memory_usage=print_mem)
    print("\nInfo: microbiomes_entities_occs", flush=True)
    microbiomes_entities_occs.info(verbose = False, memory_usage=print_mem)

    # Create output directory if it doesn't exist
    if os.path.exists(args.outdir) and not os.path.isdir(args.outdir):
        print("ERROR - The target path is not a directory", file = sys.stderr)
        sys.exit(2)
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("\nJoining protein ids with allele info...", flush=True)
    # Prepare df with allele info for proteins for downstream merging against all peptides
    proteins_allele_info = entities_proteins_occs\
        .merge(microbiomes_entities_occs)\
        .drop(columns="entity_id")\
        .merge(conditions)\
        .drop(columns="microbiome_id")\
        .merge(condition_allele_map)\
        .drop(columns="condition_id")\
        .drop_duplicates()\
        .set_index('protein_id')\
        .sort_index()
    # -> protein_id, allele_id

    print("\nInfo: proteins_allele_info", flush=True)
    proteins_allele_info.info(verbose = False, memory_usage=print_mem)

    cur_chunk = 0
    requests = 0
    keep = pd.DataFrame()
    # Process peptides chunk-wise, write out into files with max_chunk_size or keep for next chunk if remaining requests, i.e. < max_chunk_size for one allele
    with pd.read_csv(args.peptides, sep='\t', index_col="peptide_id", chunksize=args.proc_chunk_size) as reader:
        for c, peptides in enumerate(reader):
            print("\nChunk: ", c)
            print("Info: peptides", flush=True)
            peptides.info(verbose = False, memory_usage=print_mem)

            # Identify which predictions have to be computed: join peptides with allele info
            to_predict = peptides\
                    .sort_index()\
                    .join(protein_peptide_occs)\
                    .reset_index()\
                    .set_index('protein_id')\
                    .sort_index()\
                    .join(proteins_allele_info)\
                    .drop_duplicates()\
                    .reset_index(drop=True)
            # -> index, peptide_id, peptide_sequence, allele_id

            # TODO peptides can be deleted?

            print("Info: to_predict", flush=True)
            to_predict.info(verbose = False, memory_usage=print_mem)

            requests += len(to_predict)

            # Concat remaining peptides from previous round with newly processed ones
            # write the required predictions into chunks of peptide lists
            keep = pd.concat([keep, to_predict], ignore_index=True)\
                .groupby("allele_id", group_keys=False)\
                .apply(lambda x : write_chunks(x, alleles))
            # use group_keys=False to avoid generation of extra index with "allele_id"

            print("Info: keep", flush=True)
            keep.info(verbose = False, memory_usage=print_mem)

    # Write out remaining peptides
    keep.groupby("allele_id", group_keys=False).apply(lambda x : write_chunks(x, alleles, remainder=True))

    # We're happy if we got here
    print(f"All done. Written {requests} peptide prediction requests into {cur_chunk} chunks.")
    sys.exit(0)
except KeyboardInterrupt:
    print("\nUser aborted.", file = sys.stderr)
    sys.exit(1)
