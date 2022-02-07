#!/usr/bin/env python3
####################################################################################################
#
# Author: Leon Kuchenbecker
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

import pandas as pd

####################################################################################################

class PartialWeightsError(RuntimeError):
    pass

####################################################################################################

def parse_args():
    """Parse the command line arguments"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-ece",     "--entrez-conditions-entities",     required=True, metavar="PATH", type=str,                   nargs="?", help="Conditions entities map from Entrez (condition_id, entity_name, entity_weight)")
    parser.add_argument("-nce",     "--nucl-conditions-entities",       required=True, metavar="PATH", type=str,                   nargs="?", help="Conditions entities map from nucleotide methods (condition_id, entity_name, entity_weight)")
    parser.add_argument("-menw",    "--microbiomes-entities",           required=True, metavar="PATH", type=str,                              help="Microbiome entity map (microbiome_id, entity_id).")
    parser.add_argument("-ent",     "--entities",                       required=True, metavar="PATH", type=str,                              help="Entity map (entity_id, entity_name)")
    parser.add_argument("-cond",    "--conditions",                     required=True, metavar="PATH", type=str,                              help="Conditions - microbiomes map (condition_id, condition_name, microbiome_id)")
    parser.add_argument("-o",       "--output",                         required=True, metavar="PATH", type=argparse.FileType('w'),           help="Output file (condition_id, entity_id, entity_weight)")

    return parser.parse_args()


def process_weights(subset):
    """Check the integrity of the weights for a subset of the data representing
    one microbiome. If all weights are present, do nothing. If no weights are
    present, assign uniform weights. Otherwise, raise an error."""

    if subset['entity_weight'].isnull().all():
        subset['entity_weight'] = 1
    elif not subset['entity_weight'].isnull().any():
        pass
    else:
        raise PartialWeightsError(subset['condition_id'].iloc[0])
    return subset

####################################################################################################

args = parse_args()

if not args.entrez_conditions_entities and not args.nucl_conditions_entities:
    sys.exit("Neither --entrez-conditions-entities nor --nucl-conditions-entities were specified. Aborting.")

# Read and join the tables that provide microbiome_id, entity_id and entity_name
entity_microbiome = pd.read_csv(args.microbiomes_entities, sep='\t')
entity            = pd.read_csv(args.entities, sep='\t')

conditions        = pd.read_csv(args.conditions, sep="\t")

entity_condition = entity_microbiome.merge(entity).merge(conditions)

# Read the tables that provide the weights and concatenate them
input_data = pd.concat([ pd.read_csv(e, sep='\t', dtype={"entity_name":str}) for e in [args.entrez_conditions_entities, args.nucl_conditions_entities] if e ])

# Join the weights against the entity ids table, which contains all entities
# that we have observed in upstream processes. Thus, per condition, we expect
# to find weights either for all of them or for none of them.
result = entity_condition.merge(input_data, how="left").drop(columns="entity_name")

# For each condition, we now check whether this assumption is true. If we find
# no weights for a condition, we add uniform weights.
try:
    result = result.groupby("condition_id")\
            .apply(process_weights)
    result[["condition_id", "entity_id", "entity_weight"]].to_csv(args.output, sep='\t', index=False, header=True)
except PartialWeightsError as e:
    sys.exit(f"Inconsist weight specifications. Weights were specified for only a subset of entities in condition with condition ID {e}.")
