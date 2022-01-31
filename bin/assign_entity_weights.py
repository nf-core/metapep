#!/usr/bin/env python3
####################################################################################################
#
# Author: Leon Kuchenbecker
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
import argparse

import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-w',   "--weights",            required=True, metavar="PATH",  type=str,                       help="Path to the weights input file containing: weights_id, weights_path.")
    parser.add_argument('-cw',  "--conditions_weights", required=True, metavar="PATH",  type=str,                       help="Path to the condition weights map input file containing: condition_id, weights_id.")
    parser.add_argument('-o',   "--out",                required=True,                  type=argparse.FileType('w'),    help="Output TSV file (condition_id, entity_name, entity_weight)")
    return parser.parse_args(args)

args = parse_args()

weights = pd.read_csv(args.weights, sep="\t")
conditions_weights = pd.read_csv(args.conditions_weights, sep="\t")

# Read all input files and rename columns accoring to data model
dfs = []
for w in weights.itertuples():
    df = pd.read_csv(w.weights_path, sep="\t").rename(columns={
        'contig_name' : 'entity_name',
        'bin_basename' : 'entity_name',
        'weight' : 'entity_weight'
        })
    df["weights_id"] = w.weights_id
    dfs.append(df)

column_names = ["condition_id", "entity_name", "entity_weight"]

# Check if we processed any files and either concatenate them and merge with conditions table or write an empty table with only a header
if dfs:
    print(f"{len(dfs)} input tables provided, writing concatenated table.", file = sys.stderr)
    conditions_weights.merge(pd.concat(dfs))[column_names].to_csv(args.out, sep='\t', index=False, header=True)
else:
    print("No input files provided, writing empty table.", file = sys.stderr)
    print('\t'.join(column_names), file = args.out)
