#!/usr/bin/env python3
import argparse, sys
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True) #reduced prediction tsv files
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    #predcitions.tsv.gz should have columns:"peptide_id", "allele_id", "prediction_score", "rank_score"
    dfs = []
    for reduced_tsv in args.inputs:
        df = pd.read_csv(reduced_tsv, sep="\t")
        dfs.append(df)

    cols = ["peptide_id", "allele_id", "prediction_score", "rank_score"]
    merged = pd.concat([d[cols] for d in dfs], ignore_index=True)
    merged.to_csv(args.out, sep="\t", index=False, compression="gzip")
    print(f"Successfully created: {args.out} :)", file=sys.stderr)
    
if __name__ == "__main__":
    sys.exit(main())
