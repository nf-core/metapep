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
from functools import total_ordering

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input", required=True, nargs='+', metavar='FILE', help="List of peptide files")
    parser.add_argument('-s', "--samples", required=True, nargs='+', metavar='FILE', help="List of samples")
    parser.add_argument('-c', "--conditions", required=True, nargs='+', metavar='FILE', help="List of conditions")
    parser.add_argument('-t', "--type", required=True, nargs='+', metavar='FILE', help="List of types")
    parser.add_argument('-w', "--weights", required=True, nargs='+', metavar='FILE', help="List of weights ids")
    parser.add_argument('-b', "--bin_basename", required=True, nargs='+', metavar='FILE', help="List of bin basenames")
    parser.add_argument('-o', "--output", required=True, metavar='FILE', help="Output")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)
    print(args.samples)
    print(args.input)
    print(args.conditions)
    print(args.type)
    print(args.weights)
    print(args.bin_basename)

    @total_ordering
    class Peptide:
        """
        This class represents a peptide with metadata
        """
        def __init__(self, line, sample, conditions, input_type, bin_basename):
            line_split =        line.strip().split('\t')
            self.seq =          line_split[0]
            self.id =           line_split[1]
            self.protein_ids =  line_split[2]
            self.counts =       line_split[3]
            self.conditions =   conditions
            self.type =         input_type
            self.entity =       get_entity(self.protein_ids, self.type, bin_basename, sample)

        def __eq__(self, other):
            return self.seq == other.seq

        def __lt__(self, other):
            return self.seq < other.seq

        def __str__(self):
            return self.seq

        def __repr__(self):
            return self.seq

    def gen_peptides(file, sample, conditions, input_type, bin_basename):
        try:
            for line in file:
                yield Peptide(line, sample, conditions, input_type, bin_basename)
        except Exception as e:
            print(e.__traceback__)

    def empty_buffer(buffer, output_file, peptide=None):
        buffer_peptide = buffer[0].seq
        meta = map(';'.join, zip(*[[p.id, p.protein_ids, p.counts, p.conditions, p.entity, p.type] for p in buffer]))
        output_file.write('\t'.join([buffer_peptide, *meta]) + "\n")
        return [peptide]

    def get_entity(protein_ids, input_type, bin_basename, sample):
        if input_type == "assembly":
            return ';'.join(map(lambda x:"_".join(x.split("_")[:-1]), protein_ids.split(';')))
        elif input_type == "bins":
            return bin_basename
        elif input_type == "taxa":
            return sample
        return ""

    with ExitStack() as stack, open(args.output, 'w') as output_file:
        files = [stack.enter_context(open(f)) for f in args.input]
        [next(f) for f in files] # first line of every file is header
        peptides = [gen_peptides(f, s, c, t, b) for f, s, c, t, b in zip(files, args.samples, args.conditions, args.type, args.bin_basename)]
        header = ['sequence','id','protein_ids','counts', 'conditions', 'entities', 'type']
        output_file.write('\t'.join(header) + "\n")
        merged_peptides = merge(*peptides)
        buffer = [next(merged_peptides)]
        last = buffer[0]
        for peptide in merged_peptides:
            if peptide == last:  # remove duplicate peptides
                buffer.append(peptide)
            else:
                # empty buffer
                buffer = empty_buffer(buffer, output_file, peptide)
                last = buffer[0]
        # empty buffer one last time
        empty_buffer(buffer, output_file)

if __name__ == "__main__":
    sys.exit(main())
