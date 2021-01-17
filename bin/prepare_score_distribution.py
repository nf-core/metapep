#!/usr/bin/env python3
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

import argparse
import sys
import os
import csv

import pandas as pd

####################################################################################################

def parse_args(args=None):
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Prepare prediction score distribution for plotting.")

    # INPUT FILES
    parser.add_argument("-p"     , "--predictions"              , help="Path to the predictions input file"                     , type=str   , required=True)
    parser.add_argument("-ppo"   , "--protein-peptide-occ"      , help="Path to the protein peptide occurences input file"      , type=str   , required=True)
    parser.add_argument("-mpo"   , "--microbiome-protein-occ"   , help="Path to the microbiome protein occurences input file"   , type=str   , required=True)
    parser.add_argument("-c"     , "--conditions"               , help="Path to the conditions input file"                      , type=str   , required=True)
    parser.add_argument("-cam"   , "--condition-allele-map"     , help="Path to the condition allele map input file"            , type=str   , required=True)
    parser.add_argument("-a"     , "--alleles"                  , help="Path to the allele input file"                          , type=str   , required=True)

    # OUTPUT FILES
    parser.add_argument("-o"     , "--outdir"                   , help="Path to the output directory"                           , type=str   , required=True)

    # PARAMETERS
    return parser.parse_args()



def main(args=None):
    args = parse_args(args)

    # Read input files
    predictions               = pd.read_csv(args.predictions, sep='\t')
    protein_peptide_occs      = pd.read_csv(args.protein_peptide_occ, sep='\t').drop(columns="count")
    microbiome_protein_occs   = pd.read_csv(args.microbiome_protein_occ, sep='\t')
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

    # for each allele separately (to save mem)
    for allele_id in alleles.allele_id:
        print("Process allele: ", allele_id, flush=True)

        data = predictions[predictions.allele_id == allele_id] \
                .merge(protein_peptide_occs) \
                .merge(microbiome_protein_occs) \
                .drop(columns="protein_id") \
                .merge(conditions) \
                .merge(condition_allele_map) \
                .drop(columns=["allele_id", "condition_id"]) \
                .groupby(["peptide_id", "prediction_score", "microbiome_id"])["protein_weight"] \
                .sum() \
                .reset_index(name="weight_sum") \
                .drop(columns="peptide_id")

        with open(os.path.join(args.outdir, "prediction_scores.allele_" + str(allele_id) + ".tsv"), 'w') as outfile:
            data[["prediction_score", "microbiome_id", "weight_sum"]].to_csv(outfile, sep="\t", index=False, header=True)

        microbiome_ids = data.microbiome_id.drop_duplicates()
        for microbiome_id in microbiome_ids:
            print("Wrote out ", len(data[data.microbiome_id == microbiome_id]), " prediction scores for microbiome_id ", microbiome_id, ".", flush=True)


if __name__ == "__main__":
    sys.exit(main())
