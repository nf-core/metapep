#!/usr/bin/env python3
####################################################################################################
#
# Author: Sabrina Krakau
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
# NOTE
# entrez proteins of all microbiome input files already within one file (proteins.entrez.tsv.gz)
# proteins from different assemblies have non-unique names originating from enumerated contigs, should get new ids assigned separately for each file (proteins.pred_${microbiome_id}.tsv.gz)
# proteins from 'proteins' input type: not known if unique or not, handle separately for now (in case of unique ids this causes unnecessary redundancy; could add parameter for this in future)

import sys
import gzip
import csv
import io
import argparse

from Bio import SeqIO
import pandas as pd

import sys


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    # input microbiomes
    parser.add_argument(
        "-m", "--microbiomes", type=str, help="Input microbiomes tsv containing microbiome_id and microbiome_bare_id"
    )
    # Predicted Proteins
    parser.add_argument(
        "-pp", "--predicted-proteins", type=str, nargs="*", help="Protein TSV files with predicted proteins"
    )
    parser.add_argument(
        "-ppm",
        "--predicted-proteins-microbiome-ids",
        type=int,
        nargs="*",
        help="Microbiome bare ids of the predicted protein TSV files in corresponding order",
    )
    parser.add_argument(
        "-ppb",
        "--predicted-proteins-bin-basenames",
        type=str,
        nargs="*",
        help="Bin basenames of the predicted protein TSV files in corresponding order, false for type 'assembly'.",
    )
    # Entrez Proteins
    parser.add_argument(
        "-ep", "--entrez-proteins", type=str, nargs="?", help="Protein TSV file with entrez downloaded proteins"
    )
    parser.add_argument(
        "-eep",
        "--entrez-entities-proteins",
        nargs="?",
        required=True,
        type=str,
        help="TSV file associating entrez downloaded proteins with entity_name (taxon_id).",
    )
    parser.add_argument(
        "-eme",
        "--entrez-microbiomes-entities",
        nargs="?",
        required=True,
        type=str,
        help="TSV file associating entrez entity_name with microbiome_id and entity_weight.",
    )
    # Bare Proteins
    parser.add_argument(
        "-bp", "--bare-proteins", type=str, nargs="*", help="Protein TSV file with user-provided proteins"
    )
    parser.add_argument(
        "-bpm",
        "--bare-proteins-microbiome-ids",
        type=int,
        nargs="*",
        help="Microbiome ids of the user provided protein TSV files in corresponding order",
    )

    # Output
    parser.add_argument(
        "-op", "--out-proteins", type=str, required=True, help="Outputh path for the global protein TSV file."
    )
    parser.add_argument(
        "-oep",
        "--out-entities-proteins",
        type=str,
        required=True,
        help="Outputh path for the global entity protein association TSV file.",
    )
    parser.add_argument(
        "-oe", "--out-entities", type=str, required=True, help="Outputh path for the global entities TSV file."
    )
    parser.add_argument(
        "-ome",
        "--out-microbiomes-entities",
        type=str,
        required=True,
        help="Outputh path for the global microbiome entity association TSV file.",
    )

    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    next_protein_id = 0
    next_entity_id = 0

    proteins_columns = ["protein_id", "protein_orig_id", "protein_sequence"]
    entities_proteins_columns = ["entity_id", "protein_id"]
    entities_columns = ["entity_id", "entity_name"]
    microbiomes_entities_columns = ["microbiome_id", "entity_id"]
    with gzip.open(args.out_proteins, "wt") as outfile_proteins, open(
        args.out_entities_proteins, "w"
    ) as outfile_entities_proteins, open(args.out_entities, "w") as outfile_entities, open(
        args.out_microbiomes_entities, "w"
    ) as outfile_microbiomes_entities:
        # HEADERS
        print("\t".join(proteins_columns), file=outfile_proteins)
        print("\t".join(entities_proteins_columns), file=outfile_entities_proteins)
        print("\t".join(entities_columns), file=outfile_entities)
        print("\t".join(microbiomes_entities_columns), file=outfile_microbiomes_entities)

        #
        # PREDICTED PROTEINS
        #
        if args.predicted_proteins:
            # Check validity of runtime arguments
            if not args.predicted_proteins_microbiome_ids or len(args.predicted_proteins_microbiome_ids) != len(
                args.predicted_proteins
            ):
                sys.exit(
                    "An equal number of arguments has to be passed to --predicted-proteins and"
                    " --predicted_proteins_microbiome_ids"
                )
            if not args.predicted_proteins_bin_basenames or len(args.predicted_proteins_bin_basenames) != len(
                args.predicted_proteins
            ):
                sys.exit(
                    "An equal number of arguments has to be passed to --predicted-proteins and"
                    " --predicted_proteins_bin_basenames"
                )
            # Read the microbiomes table:
            microbiomes = pd.read_csv(args.microbiomes, sep="\t")
            # Read all provided files while checking in each microbiome_bare_id
            check_in_microbiome_bare_id = set()
            for microbiome_bare_id, bin_basename, inpath in zip(
                args.predicted_proteins_microbiome_ids, args.predicted_proteins_bin_basenames, args.predicted_proteins
            ):
                if microbiome_bare_id not in check_in_microbiome_bare_id:
                    check_in_microbiome_bare_id.add(microbiome_bare_id)
                    # Read and annotate proteins
                    proteins = pd.read_csv(inpath, sep="\t")
                    if bin_basename == "__ISASSEMBLY__":
                        # retrieve 'entity_name' from 'protein_tmp_id' prefix
                        proteins["entity_name"] = proteins["protein_tmp_id"].map(lambda x: "_".join(x.split("_")[:-1]))
                    else:
                        proteins["entity_name"] = bin_basename

                    proteins["protein_id"] = range(next_protein_id, next_protein_id + len(proteins))
                    next_protein_id += len(proteins)

                    # Check if microbiome is coassembly
                    if len(microbiomes.groupby("microbiome_bare_id").get_group(microbiome_bare_id)) != 1:
                        all_entities = []
                        # Iterate over microbiome_ids associated to current co-assembly
                        # (i.e. microbiome_bare_id) and assign the corresponding microbiome_id to the entities
                        for microbiome_id in microbiomes.groupby("microbiome_bare_id").get_group(microbiome_bare_id)[
                            "microbiome_id"
                        ]:
                            entities = pd.DataFrame()
                            entities = proteins[["entity_name"]].drop_duplicates()
                            entities["entity_id"] = range(next_entity_id, next_entity_id + len(entities))
                            # Instead of microbiome_bare_id append microbiome_id
                            entities["microbiome_id"] = microbiome_id
                            all_entities.append(entities)

                        next_entity_id += len(entities)
                        entities = pd.concat(all_entities)

                    else:
                        entities = proteins[["entity_name"]].drop_duplicates()
                        entities["entity_id"] = range(next_entity_id, next_entity_id + len(entities))
                        next_entity_id += len(entities)
                        # If not coassembled microbiome_id = microbiome_bare_id
                        entities["microbiome_id"] = microbiome_bare_id

                # Write proteins
                proteins.rename(columns={"protein_tmp_id": "protein_orig_id"}, inplace=True)
                proteins[proteins_columns].to_csv(outfile_proteins, sep="\t", header=False, index=False)
                # Write entities_proteins
                proteins.merge(entities)[["entity_id", "protein_id"]].drop_duplicates().to_csv(
                    outfile_entities_proteins, sep="\t", header=False, index=False
                )
                # Write entities
                entities[entities_columns].drop_duplicates().to_csv(
                    outfile_entities, sep="\t", header=False, index=False
                )
                # Write microbiomes - entities
                entities[microbiomes_entities_columns].to_csv(
                    outfile_microbiomes_entities, sep="\t", index=False, header=False
                )

        #
        # ENTREZ PROTEINS
        #
        if args.entrez_proteins:
            # Check validity of runtime arguments
            if not args.entrez_entities_proteins or not args.entrez_microbiomes_entities:
                sys.exit("The --entrez-* flags have to be specified together")
            # Read proteins and associations
            proteins = pd.read_csv(args.entrez_proteins, sep="\t")
            entities_proteins = pd.read_csv(
                args.entrez_entities_proteins, "\t"
            )  # protein_tmp_id (accessionVersion), entity_name (taxon_id)
            microbiomes_entities = pd.read_csv(
                args.entrez_microbiomes_entities, "\t"
            )  # entity_name, microbiome_id, entity_weight

            # Assign protein_id
            proteins["protein_id"] = range(next_protein_id, next_protein_id + len(proteins))
            next_protein_id += len(proteins)

            # Assign entity_id
            entities = microbiomes_entities[["entity_name"]].drop_duplicates()
            entities["entity_id"] = range(next_entity_id, next_entity_id + len(entities))
            next_entity_id += len(entities)

            # Write proteins
            proteins.rename(columns={"protein_tmp_id": "protein_orig_id"})[proteins_columns].to_csv(
                outfile_proteins, sep="\t", header=False, index=False
            )

            entities_microbiomes_proteins = (
                entities_proteins.merge(proteins)
                .merge(entities)
                .merge(microbiomes_entities)[["entity_id", "protein_id", "microbiome_id", "entity_weight"]]
            )

            # Write entities_proteins: 'entity_id', 'protein_id'
            entities_microbiomes_proteins[entities_proteins_columns].to_csv(
                outfile_entities_proteins, sep="\t", header=False, index=False
            )
            # Write entities: 'entity_id', 'entity_name'
            entities[entities_columns].to_csv(outfile_entities, sep="\t", header=False, index=False)
            # Write microbiomes - entities: 'microbiome_id', 'entity_id'
            entities_microbiomes_proteins[microbiomes_entities_columns].drop_duplicates().to_csv(
                outfile_microbiomes_entities, sep="\t", header=False, index=False
            )

        #
        # BARE PROTEINS
        #


if __name__ == "__main__":
    sys.exit(main())
