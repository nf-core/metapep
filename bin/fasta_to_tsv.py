#!/usr/bin/env python3

import argparse
import sys
import gzip
from Bio import SeqIO

parser = argparse.ArgumentParser("Converts FASTA to TSV")
parser.add_argument("-i", "--input", help="Input FASTA file", default=sys.stdin)
parser.add_argument("-o", "--output", type=argparse.FileType("o"), help="Output TSV file", default=sys.stdout)
parser.add_argument(
    "-ra", "--remove-asterisk", action="store_true", help="Remove trailing asterisks produced by prodigal"
)
args = parser.parse_args()

with gzip.open(args.input, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if args.remove_asterisk and record.seq[-1] == "*":
            print(f"{record.id}\t{record.seq[:-1]}")
        else:
            print(f"{record.id}\t{record.seq}")
