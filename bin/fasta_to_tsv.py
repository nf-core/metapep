#!/usr/bin/env python3
# Written by Sabrina Krakau, Leon Kuchenbecker, and Till Englert under the MIT license

import argparse
import gzip
import sys

from Bio import SeqIO

parser = argparse.ArgumentParser("Converts FASTA to TSV")
parser.add_argument("-i", "--input", help="Input FASTA file", default=sys.stdin)
parser.add_argument("-o", "--output", type=argparse.FileType("o"), help="Output TSV file", default=sys.stdout)
parser.add_argument(
    "-ra", "--remove-asterisk", action="store_true", help="Remove trailing asterisks produced by prodigal"
)
args = parser.parse_args()

records_out = []
with gzip.open(args.input, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if args.remove_asterisk and record.seq[-1] == "*":
            records_out.append([str(record.id),"\t",str(record.seq[:-1]),"\n"])
        else:
            records_out.append([str(record.id),"\t",str(record.seq),"\n"])
# Two dimensional array to enable sorting
records_out = sorted(records_out, key=lambda x: x[0])
print("".join(["".join(rec) for rec in records_out]))
