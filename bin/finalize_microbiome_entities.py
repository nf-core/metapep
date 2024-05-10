#!/usr/bin/env python3

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
    parser.add_argument(
        "-eme",
        "--entrez-microbiomes-entities",
        required=True,
        metavar="PATH",
        type=str,
        nargs="?",
        help="Microbiomes entities map from Entrez (microbiome_id, entity_name, entity_weight)",
    )
    parser.add_argument(
        "-nme",
        "--nucl-microbiomes-entities",
        required=True,
        metavar="PATH",
        type=str,
        nargs="?",
        help="Microbiomes entities map from nucleotide methods (microbiome_id, entity_name, entity_weight)",
    )
    parser.add_argument(
        "-menw",
        "--microbiomes-entities-noweights",
        required=True,
        metavar="PATH",
        type=str,
        help="Preliminary microbiome entity map (microbiome_id, entity_id) w/o weights.",
    )
    parser.add_argument(
        "-ent", "--entities", required=True, metavar="PATH", type=str, help="Entity map (entity_id, entity_name)"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        metavar="PATH",
        type=argparse.FileType("w"),
        help="Output file (microbiome_id, entity_id, entity_weight)",
    )

    return parser.parse_args()


def process_weights(subset):
    """Check the integrity of the weights for a subset of the data representing
    one microbiome. If all weights are present, do nothing. If no weights are
    present, assign uniform weights. Otherwise, raise an error."""

    if subset["entity_weight"].isnull().all():
        subset["entity_weight"] = 1
    elif not subset["entity_weight"].isnull().any():
        pass
    else:
        raise PartialWeightsError(subset["microbiome_id"].iloc[0])
    return subset


####################################################################################################

args = parse_args()

if not args.entrez_microbiomes_entities and not args.nucl_microbiomes_entities:
    sys.exit("Neither --entrez-microbiome-entities nor --nucl-microbiome-entities were specified. Aborting.")

# Read and join the tables that provide microbiome_id, entity_id and entity_name
entity_microbiome = pd.read_csv(args.microbiomes_entities_noweights, sep="\t")
entity = pd.read_csv(args.entities, sep="\t")

entity_microbiome = entity_microbiome.merge(entity)

# Read the tables that provide the weights and concatenate them
input_data = pd.concat(
    [pd.read_csv(e, sep="\t") for e in [args.entrez_microbiomes_entities, args.nucl_microbiomes_entities] if e]
)

# Join the weights against the entity ids table, which contains all entities
# that we have observed in upstream processes. Thus, per microbiome, we expect
# to find weights either for all of them or for none of them.
result = entity_microbiome.merge(input_data, how="left").drop(columns="entity_name")

# For each microbiome, we now check whether this assumption is true. If we find
# no weights for a microbiome, we add uniform weights.
try:
    result = result.groupby("microbiome_id").apply(process_weights)
    result.to_csv(args.output, sep="\t", index=False, header=True)
except PartialWeightsError as e:
    sys.exit(
        "Inconsist weight specifications. Weights were specified for only a subset of entities in microbiome with"
        f" microbiome ID {e}."
    )
