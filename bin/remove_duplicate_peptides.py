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
from contextlib import ExitStack

import urllib.request
import os
import gzip

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   "--input", required=True, nargs='+', metavar='FILE', help="List of peptide files")
    parser.add_argument('-s',   "--samples", required=True, nargs='+', metavar='STR', help="List of samples")
    parser.add_argument('-c',   "--conditions", required=True, nargs='+', metavar='STR', help="List of conditions")
    parser.add_argument('-t',   "--type", required=True, nargs='+', metavar='STR', help="List of types")
    parser.add_argument('-wi',  "--weights_ids", required=True, nargs='+', metavar='INT', help="List of weights ids")
    parser.add_argument('-w',   "--weights_table", required=True, metavar='FILE', help="Table containing weights_id, weights_path")
    parser.add_argument('-b',   "--bin_basename", required=True, nargs='+', metavar='STR', help="List of bin basenames")
    parser.add_argument('-f',   "--filter_by_weights", default=False, action='store_true', help="Discard peptides with weight 0")
    parser.add_argument('-bs',  "--buffer_size", default=1, type=int, metavar='INT', help="Size of external merge sort buffers")
    parser.add_argument('-o',   "--output", required=True, metavar='FILE', help="Output")
    return parser.parse_args(args)

class Peptide:
    """
    This class represents a peptide with metadata
    """
    def __init__(self, seq, peptide_id, protein_ids, counts, sample, condition, input_type, bin_basename, weights_id, weights):
        self.seq =                      seq
        self.id =                       peptide_id
        self.protein_ids =              protein_ids
        self.counts =                   counts
        self.condition =                condition
        self.type =                     input_type
        self.entities, self.weights =   get_entities_and_weights(self.protein_ids, input_type, bin_basename, sample, weights_id, weights)
        # TODO: Currently input weights are only supported for type 'assembly', implement for other types

    def __eq__(self, other):
        return self.seq == other.seq

    def __lt__(self, other):
        return self.seq < other.seq

    def __le__(self, other):
        return self.seq <= other.seq

    def __gt__(self, other):
        return self.seq > other.seq

    def __ge__(self, other):
        return self.seq >= other.seq

    def __str__(self):
        return self.seq

    def __repr__(self):
        return self.seq

def gen_peptides(peptide_handle, sample, conditions, input_type, bin_basename, weights_ids, weights, buffer_size):
    try:
        buffer = list(peptide_handle.readlines(buffer_size))
        while buffer:
            for line in buffer:
                line_split =   line.strip().split('\t')
                seq =          line_split[0]
                peptide_id =   line_split[1]
                protein_ids =  line_split[2]
                counts =       line_split[3]
                # for each condition yield one peptide (weights can differ between conditions)
                for condition, weights_id in zip(conditions.split(';'), weights_ids.split(';')):
                    yield Peptide(seq, peptide_id, protein_ids, counts, sample, condition, input_type, bin_basename, weights_id, weights)
            buffer = list(peptide_handle.readlines(buffer_size))
    except Exception as e:
        print(e)

def get_entities_and_weights(protein_ids, input_type, bin_basename, sample, weights_id, weights):
    if input_type == "assembly":
        collect_entities = []
        collect_weights = []
        for protein_id in protein_ids.split(','):
            try:
                entity = "_".join(protein_id.split("_")[:-1])
                weight = weights[int(weights_id)][entity] if weights_id and int(weights_id) in weights.keys() else 1
            except KeyError as e:
                print(f"Incomplete weights: entity {e} missing!")
            if float(weight > 0):
                collect_entities.append(entity)
                collect_weights.append(str(weight))
        return ",".join(collect_entities), ",".join(collect_weights)
    elif input_type == "bins":
        return bin_basename, "1"
    elif input_type == "taxa":
        return sample, "1"
    return ""

def create_weights_dictionary(weights_table):
    d = {}
    with open(weights_table, "r") as weights_entry:
        weights_entry.readline() # ignore header
        for line in weights_entry:
            weights_id, weights_path = line.strip().split("\t")
            # open weights_path, can be local or url
            if os.path.isfile(weights_path):
                with open(weights_path, "r") as entity_weights:
                    next(entity_weights) # ignore header
                    d[int(weights_id)] = {str(contig): float(weight) for contig, weight in (entity_weight.strip().split('\t') for entity_weight in entity_weights)}
            else:
                try:
                    with urllib.request.urlopen(weights_path) as entity_weights:
                        next(entity_weights) # ignore header
                        d[int(weights_id)] = {str(contig): float(weight) for contig, weight in (entity_weight.decode('utf-8').strip().split('\t') for entity_weight in entity_weights)}
                except urllib.error.URLError as e:
                    print(e.reason)
    return d

def empty_buffer(buffer, output_file, duplicate_peptides, peptide=None, print_buffer=False):
    buffer_peptide = buffer[0].seq
    duplicate_peptides = duplicate_peptides +  len(buffer) - 1 # count duplicate peptides
    meta = map(';'.join, zip(*[[p.id, p.protein_ids, p.counts, p.condition, p.entities, p.weights, p.type] for p in buffer if p.weights]))
    meta_list = list(meta)
    if len(meta_list) > 0: # if there are no weights don't write peptide
        output_file.write('\t'.join([buffer_peptide, *meta_list]) + "\n")
        if print_buffer:
            print('\t'.join([buffer_peptide, *meta_list]))
    return [peptide], duplicate_peptides



def main(args=None):
    args = parse_args(args)

    # read weights into dictionary for lookup
    weights = create_weights_dictionary(args.weights_table)

    with ExitStack() as stack, gzip.open(args.output, 'wt') as output_file:
        files = [stack.enter_context(open(f)) for f in args.input]
        [next(f) for f in files] # ignore header

        # generate peptides
        peptides = [gen_peptides(f, s, c, t, b, w, weights, args.buffer_size) for f, s, c, t, b, w in zip(files, args.samples, args.conditions, args.type, args.bin_basename, args.weights_ids)]

        # write header to output file
        header = ['sequence','ids','protein_ids','counts', 'conditions', 'entities', 'weights', 'types']
        output_file.write('\t'.join(header) + "\n")

        # merge peptides (peptides given as input have to be sorted)
        merged_peptides = merge(*peptides)

        # initialize buffer and set peptide as last
        buffer = [next(merged_peptides)]
        last = buffer[0]

        # initialize counters
        n = 1
        duplicate_peptides = 0
        print_buffer = False

        # filter duplicate peptides
        for peptide in merged_peptides:
            n = n+1
            n_step = (n % 1000 == 0)
            if n_step:
                # print to stdout to keep track of how many peptides are done
                print(f"peptide no.: {n}", flush=True)
                print(f"duplicates: {duplicate_peptides}", flush=True)
                print_buffer = True
            if peptide == last:  # duplicate peptides
                buffer.append(peptide)
            else:
                # empty buffer
                buffer, duplicate_peptides = empty_buffer(buffer, output_file, duplicate_peptides, peptide, print_buffer)
                last = buffer[0]
                print_buffer = False

        # empty buffer one last time
        buffer, duplicate_peptides = empty_buffer(buffer, output_file, duplicate_peptides)

        # print stats
        print(f"peptide no.: {n}", flush=True)
        print(f"duplicates: {duplicate_peptides}", flush=True)

if __name__ == "__main__":
    sys.exit(main())
