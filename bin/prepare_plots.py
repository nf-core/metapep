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
import csv
import re
from collections import defaultdict

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   "--input", required=True, metavar='FILE', help="List of epitope prediction files")
    parser.add_argument('-a',   "--allele", required=True, metavar='STR', help="Allele name")
    return parser.parse_args(args)

def call_binder(score, method):
    if method == "syfpeithi":
        return float(score) >= 0.50
    else:
        return float(score) <= 500

def get_binding_ratio(no_binder, n, condition, entity):
    return str(no_binder[condition][entity]/n[condition])

def prepare_plots(reader, pred_scores_file, entity_ratios_file):
    n = defaultdict(int)
    no_binder = defaultdict(lambda: defaultdict(int))
    weight = defaultdict(lambda: defaultdict(float))
    for line in reader:
        binder = bool(line["binder"])
        score = float(line["score"])
        score_dist = []
        for condition, counts, entities, entity_weights in zip(line["conditions"].split(";"), line["counts"].split(";"), line["entities"].split(";"), line["weights"].split(";")):
            # prepare score distribution
            score_dist.append([str(score), condition, str(sum([float(entity_weight) for entity_weight in entity_weights.split(",")]))])
            # prepare entity binding ratio
            n[condition] += sum([int(c) for c in counts.split(",")])
            if binder:
                for entity, entity_weight in zip(entities.split(","), entity_weights.split(",")):
                    no_binder[condition][entity] += sum([int(c) for c in counts.split(",")])
                    weight[condition][entity] = float(entity_weight)
        for entry in score_dist:
            print("\t".join(entry), file=pred_scores_file)
    for condition in n.keys():
        for entity in no_binder[condition].keys():
            entry = [condition, get_binding_ratio(no_binder, n, condition, entity), str(weight[condition][entity])]
            print("\t".join(entry), file=entity_ratios_file)


def main(args=None):
    args = parse_args(args)
    with open(args.input, "r", newline='') as file, open(f"prediction_scores.{args.allele}.tsv", "w") as pred_scores_file, open(f"entity_binding_ratios.{args.allele}.tsv", "w") as entity_ratios_file:
        header = file.readline().strip().split("\t")
        header = [re.sub(r'.*affinity', 'score', h) for h in header]
        header = [re.sub(r'.*binder', 'binder', h) for h in header]
        reader = csv.DictReader(file, delimiter='\t', fieldnames=header)
        header_entity_ratios = ["condition_name", "binding_rate","entity_weight"]
        entity_ratios_file.write("\t".join(header_entity_ratios) + "\n")
        header_pred_score = ["prediction_score", "condition_name", "weight_sum"]
        pred_scores_file.write("\t".join(header_pred_score) + "\n")
        prepare_plots(reader, pred_scores_file, entity_ratios_file)

if __name__ == "__main__":
    sys.exit(main())
