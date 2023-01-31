#!/usr/bin/env python3
#
# Author: Sabrina Krakau <sabrina.krakau@uni-tuebingen.de>
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

import pandas as pd
import math

####################################################################################################


def parse_args(args=None):
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(
        description="Concatenate TSV files."
    )

    # INPUT FILES
    parser.add_argument("-i", "--input", help="Path to input files.", type=str, required=True, nargs="+")

    # OUTPUT FILES
    parser.add_argument("-o", "--output", help="Path to output file.", type=str, required=True)

    # PARAMETERS
    parser.add_argument(
        "-c",
        "--chunk-size",
        help="Chunk size used to process TSV files. Default: 10000",
        type=int,
        default=10000,
    )

    return parser.parse_args()


def main(args=None):
    args = parse_args(args)

    first_header = pd.DataFrame().columns
    for i, filename in enumerate(args.input):
        print("Processing file: ", filename, flush=True)

        # Read input file chunk-wise
        with pd.read_csv(filename, sep="\t", chunksize=args.chunk_size) as reader:
            for j, tsv_chunk in enumerate(reader):
                print(" Chunk: ", j, flush=True)
                if (i == 0 and j == 0):
                    first_header = tsv_chunk.columns
                    print("Header: ", first_header.tolist(), flush=True)
                    tsv_chunk.to_csv(
                        args.output, mode="w", sep="\t", index=False, header=True
                    )
                else:
                    if (j == 0):
                        # Check if header of subsequent input files match header of first input file
                        # (column order must be the same)
                        if (tsv_chunk.columns.tolist() != first_header.tolist()):
                            print("ERROR - header of input file", filename, "does not match the header of the first input file!", file=sys.stderr)
                            sys.exit(1)

                    tsv_chunk.to_csv(
                        args.output, mode="a", sep="\t", index=False, header=False
                    )


if __name__ == "__main__":
    sys.exit(main())
