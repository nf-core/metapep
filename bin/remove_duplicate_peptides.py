#!/usr/bin/env python3
####################################################################################################
#
# Author: Antonia Schuster
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
from heapq import merge
from itertools import count, islice
from contextlib import ExitStack, contextmanager
import csv

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input", required=True, nargs='+', metavar='FILE', help="List of peptide files")
    parser.add_argument('-e', "--entities", required=True, nargs='+', metavar='FILE', help="List of entities")
    parser.add_argument('-o', "--output", required=True, metavar='FILE', help="Output")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    def myopen(file, entity):
        try:
            for line in csv.reader(file, delimiter='\t'):
                yield line, entity
        except Exception as e:
            print(e.__traceback__)
    # with open(args.input) as csvfile:
    #     file_list = list(csv.reader(csvfile, delimiter='\t'))
    # print([f for f, s in file_list])
    with ExitStack() as stack, open(args.output, 'w') as output_file:
        files = [myopen(stack.enter_context(open(f)), e) for f, e in zip(args.input, args.entities)]
        header = [next(f)[0] for f in files][0] # sequence, id, protein_ids, counts
        output_file.write('\t'.join(header) + "\n")
        merged_files = merge(*files)
        (peptide_seq, peptide_id, protein_ids, counts), entity = next(merged_files)
        buffer = [[peptide_seq, peptide_id, f"{entity}_{protein_ids}", counts]]
        last = buffer[0][0]
        for line in merged_files:
            (peptide_seq, peptide_id, protein_ids, counts), entity = line
            if peptide_seq == last:  # remove duplicate peptides
                buffer.append([peptide_seq, peptide_id, f"{entity}_{protein_ids}", counts])
            else:
                # empty buffer
                new_peptides, new_id, new_proteins, new_counts = list(map(list, zip(*buffer)))
                output_file.write('\t'.join([new_peptides[0], ';'.join(new_id), ';'.join(new_proteins), ';'.join(new_counts)]) + "\n")
                last = peptide_seq
                buffer = [[peptide_seq, peptide_id, f"{entity}_{protein_ids}", counts]]
        # empty buffer one last time
        new_peptides, new_id, new_proteins, new_counts = list(map(list, zip(*buffer)))
        output_file.write('\t'.join([new_peptides[0], ';'.join(new_id), ';'.join(new_proteins), ';'.join(new_counts)]))



if __name__ == "__main__":
    sys.exit(main())
