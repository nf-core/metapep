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


def call_binder(score, method):
    if method == "syfpeithi":
        return score >= 0.50
    else:
        return score <= 500

def get_binding_ratio(df):
    df_new = df \
            .groupby(["condition_name", "entity_weight"])["binder"] \
            .sum() \
            .reset_index(name="binding_rate")
    df_new["binding_rate"] = df_new["binding_rate"]/float(len(df))
    # -> index, condition_name, entity_weight, binding_rate
    return df_new


def main(args=None):
    args = parse_args(args)

    # Read input files and downcast (critical) df columns that will not be used as indices to save mem usage
    # (loading predictions allele-wise with skiprows doesn't seem to work, could still be considered to be filtered afterwards if this becomes bottleneck)
    predictions                     = pd.read_csv(args.predictions, index_col=['peptide_id', 'allele_id'], sep='\t').sort_index()
    predictions["prediction_score"] = pd.to_numeric(predictions["prediction_score"], downcast="float")

    # (binding ratio: peptide occurrences within multiple proteins of an entity are counted, while occurrences within the same protein are not considered currently)
    protein_peptide_occs         = pd.read_csv(args.protein_peptide_occ, usecols=['protein_id', 'peptide_id'], index_col="protein_id", sep='\t').sort_index()
    entities_proteins_occs       = pd.read_csv(args.entities_proteins_occ, sep='\t')
    microbiomes_entities_occs    = pd.read_csv(args.microbiomes_entities_occ, sep='\t')
    conditions                   = pd.read_csv(args.conditions, sep='\t')
    # convert condition_name column to categorical datatype to reduce memory usage (condition_id column is anyway not used for bigger dfs)
    conditions["condition_name"] = conditions["condition_name"].astype("category")
    condition_allele_map         = pd.read_csv(args.condition_allele_map, sep='\t')
    alleles                      = pd.read_csv(args.alleles, sep='\t')

    print_mem = 'deep'      # 'deep' (extra computational costs) or None
    print("\nInfo: predictions", flush=True)
    predictions.info(verbose = False, memory_usage=print_mem)
    print("\nInfo: protein_peptide_occs", flush=True)
    protein_peptide_occs.info(verbose = False, memory_usage=print_mem)

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

        write_header = True
        with open(os.path.join(args.outdir, "entity_binding_ratios.allele_" + str(allele_id) + ".tsv"), 'w') as outfile:
            # Process predictions entity_wise to reduce memory usage
            for entity_id in microbiomes_entities_occs["entity_id"].drop_duplicates():
                print("Entity_id: ", entity_id, flush=True)

                # Prepare df for joining protein information
                protein_info = microbiomes_entities_occs[microbiomes_entities_occs.entity_id == entity_id] \
                    .merge(conditions) \
                    .drop(columns="microbiome_id") \
                    .merge(condition_allele_map[condition_allele_map.allele_id == allele_id]) \
                    .drop(columns="condition_id") \
                    .merge(entities_proteins_occs) \
                    .drop(columns="entity_id") \
                    .set_index('protein_id') \
                    .sort_index()
                # -> entity_weight, condition_name, allele_id, protein_id
                # merged against condition_allele_map to keep only entities, and thus proteins, for which a prediction is requested for the current allele

                if (len(protein_info) != len(protein_info.reset_index().drop_duplicates())):
                    print("ERROR - protein_info dataframe contains duplicates!", file = sys.stderr)
                    sys.exit(1)

                # Prepare data for current entity_id
                # (avoid copying full peptide dfs -> start with entity subsets)
                data = protein_info \
                        .join(protein_peptide_occs) \
                        .reset_index(drop=True) \
                        .set_index(['peptide_id', 'allele_id']) \
                        .sort_index() \
                        .join(predictions) \
                        .reset_index(drop=True)
                # -> index, entity_weight, condition_name, prediction_score (multiple rows for occurrences in multiple proteins within entity_id)

                data["binder"] = data["prediction_score"].apply(call_binder, method=args.method)
                # print("\nInfo: data 1", flush=True)
                # data.info(verbose = False, memory_usage=print_mem)

                data = data \
                        .drop(columns="prediction_score") \
                        .groupby(["condition_name", "entity_weight"], group_keys=False).apply(lambda x : get_binding_ratio(x)) \
                        .reset_index(drop=True)
                # print("\nInfo: data 2", flush=True)
                # data.info(verbose = False, memory_usage=print_mem)

                if not data.empty:
                    data[["condition_name", "binding_rate", "entity_weight"]].to_csv(outfile, sep="\t", index=False, mode='a', header=write_header)
                    write_header = False

        # Sanity check, remove
        for condition_name in data.condition_name.drop_duplicates():
            print("Wrote out ", len(data[data.condition_name == condition_name]), " entity binding rates for condition_name ", condition_name, ".", flush=True)


if __name__ == "__main__":
    sys.exit(main())
