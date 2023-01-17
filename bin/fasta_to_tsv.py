#!/usr/bin/env python3
#
# Author: Leon Kuchenbecker <leon.kuchenbecker@uni-tuebingen.de>
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

import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser("Converts FASTA to TSV")
parser.add_argument("-i", "--input", type=argparse.FileType("r"), help="Input FASTA file", default=sys.stdin)
parser.add_argument("-o", "--output", type=argparse.FileType("o"), help="Output TSV file", default=sys.stdout)
parser.add_argument(
    "-ra", "--remove-asterisk", action="store_true", help="Remove trailing asterisks produced by prodigal"
)
args = parser.parse_args()

for record in SeqIO.parse(args.input, "fasta"):
    if args.remove_asterisk and record.seq[-1] == "*":
        print(f"{record.id}\t{record.seq[:-1]}")
    else:
        print(f"{record.id}\t{record.seq}")
