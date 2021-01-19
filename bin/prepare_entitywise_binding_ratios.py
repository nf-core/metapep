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
    parser.add_argument("-epo"   , "--entities-proteins-occ"    , help="Path to the entity protein occurences input file"       , type=str   , required=True)
    parser.add_argument("-meo"   , "--microbiomes-entities-occ" , help="Path to the microbiome entity occurences input file"    , type=str   , required=True)
    parser.add_argument("-c"     , "--conditions"               , help="Path to the conditions input file"                      , type=str   , required=True)
    parser.add_argument("-cam"   , "--condition-allele-map"     , help="Path to the condition allele map input file"            , type=str   , required=True)
    parser.add_argument("-a"     , "--alleles"                  , help="Path to the allele input file"                          , type=str   , required=True)
    parser.add_argument("-m"     , "--method"                   , help="Used epitope prediction method"                         , type=str   , required=True)

    # OUTPUT FILES
    parser.add_argument("-o"     , "--outdir"                   , help="Path to the output directory"                           , type=str   , required=True)

    # PARAMETERS
    return parser.parse_args()


def get_binder(score):
    if args.method == "syfpeithi":
        if score >= 0.50:
            return True
        else
            return False
    else:
        if score <= 500:
            return True
        else
            return False


def main(args=None):
    args = parse_args(args)

    # Read input files
    predictions               = pd.read_csv(args.predictions, sep='\t')
    protein_peptide_occs      = pd.read_csv(args.protein_peptide_occ, sep='\t').drop(columns="count")
    entities_proteins_occs    = pd.read_csv(args.entities_proteins_occ, sep='\t')
    microbiomes_entities_occs = pd.read_csv(args.microbiomes_entities_occ, sep='\t')
    conditions                = pd.read_csv(args.conditions, sep='\t')
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

        # entity-wise: do not account for entity-weight?
        data = predictions[predictions.allele_id == allele_id] \
                .merge(protein_peptide_occs) \
                .merge(entities_proteins_occs) \
                .drop(columns="protein_id") \
                .merge(microbiomes_entities_occs) \
                .merge(conditions) \
                .drop(columns="microbiome_id") \
                .merge(condition_allele_map) \
                .drop(columns=["allele_id", "condition_id"]) \
                .groupby(["entity_name", "peptide_id", "condition_name", "entity_weight"])["prediction_score"].apply(lambda s: (get_binder(s)).sum()/len(s)) \
                .reset_index(name="binding_rate") \
                .drop(columns=["entity_name", "peptide_id"])

        # TODO maybe keep "entity_name" to explore potential outliers?

        with open(os.path.join(args.outdir, "entity_binding_rates.allele_" + str(allele_id) + ".tsv"), 'w') as outfile:
            data[["condition_name", "binding_rate", "entity_weight"]].to_csv(outfile, sep="\t", index=False, header=True)

        # Sanity check, remove
        for condition_name in data.condition_name.drop_duplicates():
            print("Wrote out ", len(data[data.condition_name == condition_name]), " entity binding rates for condition_name ", condition_name, ".", flush=True)


if __name__ == "__main__":
    sys.exit(main())
