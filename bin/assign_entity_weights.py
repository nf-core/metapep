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
# for each strain: select largest assembly (for now)

import sys
import argparse

import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', "--microbiome-ids", metavar="ID", nargs="*", type=int, help="List of microbiome IDs")
    parser.add_argument('-w', "--weights-files", metavar="PATH", nargs="*", type=str, help="List of corresponding weights files (entity_name, weight)")
    parser.add_argument('-o', "--out", required=True, type=argparse.FileType('w'), help="Output TSV file (entity_name, microbiome_id, weight)")
    return parser.parse_args(args)

args = parse_args()

if len(args.microbiome_ids) != len(args.weights_files):
    sys.exit("The number of microbiome IDs and weights files has to be identical!")

dfs = []
for mb_id, w_path in zip(args.microbiome_ids, args.weights_files):
    data = pd.read_csv(w_path, sep='\t').rename(columns={
        'contig_name' : 'entity_name',
        'bin_basename' : 'entity_name'
        })
    data['microbiome_id'] = mb_id
    dfs.append(data)

column_names = ["entity_name", "microbiome_id", "weight"]

if dfs:
    print(f"{args.microbiome_ids} input tables provided, writing concatenated table.", file = sys.stderr)
    pd.concat(dfs)[column_names].to_csv(args.out, sep='\t', index=False, header=True)
else:
    print("No input files provided, writing empty table.", file = sys.stderr)
    print('\t'.join(column_names), file = args.out)

