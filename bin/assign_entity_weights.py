#!/usr/bin/env python3
# Written by Sabrina Krakau, Leon Kuchenbecker, and Till Englert under the MIT license

import argparse
import sys

import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--microbiome-ids", metavar="ID", nargs="*", type=int, help="List of microbiome IDs")
    parser.add_argument(
        "-w",
        "--weights-files",
        metavar="PATH",
        nargs="*",
        type=str,
        help="List of corresponding weights files (entity_name, weight)",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        type=argparse.FileType("w"),
        help="Output TSV file (entity_name, microbiome_id, entity_weight)",
    )
    return parser.parse_args(args)


args = parse_args()

# Check user arguments for validity
if len(args.microbiome_ids) != len(args.weights_files):
    sys.exit("The number of microbiome IDs and weights files has to be identical!")

# Read all input files, rename columns accoring to data model and add microbiome_id
dfs = []
for mb_id, w_path in zip(args.microbiome_ids, args.weights_files):
    data = pd.read_csv(w_path, sep="\t").rename(
        columns={"contig_name": "entity_name", "bin_basename": "entity_name", "weight": "entity_weight"}
    )
    data["microbiome_id"] = mb_id
    dfs.append(data)

column_names = ["entity_name", "microbiome_id", "entity_weight"]

# Check if we processed any files and either concatenate them or write an empty table with only a header
if dfs:
    print(f"{len(args.microbiome_ids)} input tables provided, writing concatenated table.", file=sys.stderr)
    pd.concat(dfs)[column_names].to_csv(args.out, sep="\t", index=False, header=True)
else:
    print("No input files provided, writing empty table.", file=sys.stderr)
    print("\t".join(column_names), file=args.out)
