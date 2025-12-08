#!/usr/bin/env python3

import argparse
import sys

import pandas as pd

# Required columns in raw prediction files
raw_prediction_cols = {"sequence", "allele", "predictor", "BA", "rank", "binder"}

# Target schema after normalization
target_cols = ["peptide_id", "allele_id", "prediction_score", "rank"]


def parse_args(args=None):
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Concatenate CSV files into a normalized TSV file.")

    # INPUT FILES
    parser.add_argument("-i", "--input", help="Path to input files.", type=str, required=True, nargs="+")

    # OUTPUT FILE
    parser.add_argument("-o", "--output", help="Path to output TSV file.", type=str, required=True)

    # PARAMETERS
    parser.add_argument(
        "-c",
        "--chunk-size",
        help="Chunk size used to process CSV files. Default: 10000",
        type=int,
        default=10000,
    )

    # MAPPING FILES
    parser.add_argument(
        "--pepmap",
        help="Path to peptides.tsv.gz for mapping sequence to peptide_id.",
        type=str,
        required=False,
        default=None,
    )

    parser.add_argument(
        "--allelemap",
        help="Path to alleles.tsv for mapping allele name to allele_id.",
        type=str,
        required=False,
        default=None,
    )

    return parser.parse_args(args)


def normalize_chunk(df, pepmap=None, allelemap=None, src_name=None):
    """Normalize prediction chunk to target schema -> columns: peptide_id, allele_id, sequence, allele, rank, prediction_score, binder, predictor"""

    # Check if already normalized (has both peptide_id and allele_id columns with values)
    already_normalized = "peptide_id" in df.columns and "allele_id" in df.columns

    if already_normalized:
        # Check if allele_id actually has values (not all empty)
        has_allele_values = df["allele_id"].notna().any()
        if has_allele_values:
            print(f" File already normalized with IDs, skipping mapping", flush=True)
            # Ensure all target columns exist
            for c in target_cols:
                if c not in df.columns:
                    df[c] = pd.NA
            return df[target_cols]
        else:
            print(f" File has ID columns but they are empty, will remap", flush=True)

    # Peptide ID mapping
    if pepmap is not None and ("peptide_id" not in df.columns or df["peptide_id"].isna().all()):
        if "sequence" not in df.columns:
            print(f"ERROR: 'sequence' column not found in {src_name}. Available columns: {df.columns.tolist()}", file=sys.stderr)
            sys.exit(1)

        seq_col = "peptide_sequence" if "peptide_sequence" in pepmap.columns else "sequence"

        # Drop peptide_id if it exists but is empty
        if "peptide_id" in df.columns:
            df = df.drop(columns=["peptide_id"])

        df = df.merge(
            pepmap[["peptide_id", seq_col]],
            left_on="sequence",
            right_on=seq_col,
            how="left"
        ).drop(columns=[seq_col], errors="ignore")

    # Allele ID mapping
    if allelemap is not None and ("allele_id" not in df.columns or df["allele_id"].isna().all()):
        # Find the allele name column in the mapping file
        allele_col_in_map = None
        for possible_col in ["allele_name", "allele"]:
            if possible_col in allelemap.columns:
                allele_col_in_map = possible_col
                break

        if allele_col_in_map is None:
            print(f"WARNING: No allele name column found in allelemap. Available columns: {allelemap.columns.tolist()}", file=sys.stderr)
        else:
            df_allele_col = "allele" if "allele" in df.columns else None

            if df_allele_col:
                print(f" Mapping alleles", file=sys.stderr)

                # Drop allele_id if it exists but is empty
                if "allele_id" in df.columns:
                    df = df.drop(columns=["allele_id"])

                # Direct mapping
                df = df.merge(
                    allelemap[["allele_id", allele_col_in_map]],
                    left_on=df_allele_col,
                    right_on=allele_col_in_map,
                    how="left"
                ).drop(columns=[allele_col_in_map], errors="ignore")

                # Check for unmapped alleles
                unmapped = df[df["allele_id"].isna() & df[df_allele_col].notna()]
                if len(unmapped) > 0:
                    unique_unmapped = unmapped[df_allele_col].unique()
                    print(f"WARNING: {len(unmapped)} rows have unmapped alleles.", file=sys.stderr)
                    print(f"  Unmapped alleles: {unique_unmapped.tolist()}", file=sys.stderr)
                    print(f"  Available in mapping: {allelemap[allele_col_in_map].unique().tolist()}", file=sys.stderr)

    # Rename BA to prediction_score for downstream analysis
    if "BA" in df.columns:
        df = df.rename(columns={"BA": "prediction_score"})

    # Ensure all target columns exist
    for c in target_cols:
        if c not in df.columns:
            df[c] = pd.NA

    # for safety reasons drop any unwanted columns
    df = df.drop(columns=["predictor", "binder", "sequence", "allele"], errors="ignore")

    return df[target_cols]


def main(args=None):
    args = parse_args(args)

    # Load mapping files
    pepmap = pd.read_csv(args.pepmap, sep="\t", compression="infer") if args.pepmap else None
    allelemap = pd.read_csv(args.allelemap, sep="\t", compression="infer") if args.allelemap else None

    if allelemap is not None:
        print(f"Loaded allele mapping with {len(allelemap)} entries", flush=True)
        allele_col = "allele_name" if "allele_name" in allelemap.columns else "allele"
        print(f"Alleles in mapping: {allelemap[allele_col].unique().tolist()}", flush=True)

    first_header = pd.DataFrame().columns
    for i, filename in enumerate(args.input):
        print("Processing file: ", filename, flush=True)

        # Read input file chunk-wise
        with pd.read_csv(filename, sep=None, engine='python', chunksize=args.chunk_size) as reader:
            for j, csv_chunk in enumerate(reader):
                print(" Chunk: ", j, flush=True)

                # Normalize chunk
                csv_chunk = normalize_chunk(csv_chunk, pepmap, allelemap, filename)
                if i == 0 and j == 0:
                    first_header = csv_chunk.columns
                    print("Header: ", first_header.tolist(), flush=True)
                    csv_chunk.to_csv(args.output, mode="w", sep="\t", index=False, header=True)
                else:
                    if j == 0 :
                        # Check if header of subsequent input files match header of first input file
                        # (column order must be the same)

                        if csv_chunk.columns.tolist() != first_header.tolist():
                            print(
                                "ERROR - header of input file",
                                filename,
                                "does not match the header of the first input file!",
                                file=sys.stderr,
                            )
                            sys.exit(1)

                    csv_chunk.to_csv(args.output, mode="a", sep="\t", index=False, header=False)


if __name__ == "__main__":
    sys.exit(main())
