#!/usr/bin/env python3

import argparse
import sys
import pandas as pd


def parse_args(args=None):
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Concatenate TSV files.")

    # INPUT FILES
    parser.add_argument("-i", "--input", help="Path to input files.", type=str, required=True, nargs="+")

    # OUTPUT FILES
    parser.add_argument("-o", "--output", help="Path to output file.", type=str, required=True)

    # PARAMETERS
    parser.add_argument(
        "-c",
        "--chunk-size",
        help="Chunk size used to process TSV files. Chunk size refers to number of processed lines within tsv file. Default: 10000",
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
                if i == 0 and j == 0:
                    first_header = tsv_chunk.columns
                    print("Header: ", first_header.tolist(), flush=True)
                    tsv_chunk.to_csv(args.output, mode="w", sep="\t", index=False, header=True)
                else:
                    if j == 0:
                        # Check if header of subsequent input files match header of first input file
                        # (column order must be the same)
                        if tsv_chunk.columns.tolist() != first_header.tolist():
                            print(
                                "ERROR - header of input file",
                                filename,
                                "does not match the header of the first input file!",
                                file=sys.stderr,
                            )
                            sys.exit(1)

                    tsv_chunk.to_csv(args.output, mode="a", sep="\t", index=False, header=False)


if __name__ == "__main__":
    sys.exit(main())
