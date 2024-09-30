#!/usr/bin/env python3

import argparse
import os
import sys

import pandas as pd

####################################################################################################


def parse_args():
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(
        description="Generate chunks of peptides that are to be predicted against a specified allele."
    )

    # INPUT FILES
    parser.add_argument("-p", "--peptides", help="Path to the peptides input file", type=str, required=True)
    parser.add_argument(
        "-ppo",
        "--protein-peptide-occ",
        help="Path to the protein peptide occurences input file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-epo",
        "--entities-proteins-occ",
        help="Path to the entity protein occurences input file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-meo",
        "--microbiomes-entities-occ",
        help="Path to the microbiome entity occurences input file",
        type=str,
        required=True,
    )
    parser.add_argument("-c", "--conditions", help="Path to the conditions input file", type=str, required=True)
    parser.add_argument(
        "-cam", "--condition-allele-map", help="Path to the condition allele map input file", type=str, required=True
    )
    parser.add_argument("-a", "--alleles", help="Path to the allele input file", type=str, required=True)

    # OUTPUT FILES
    parser.add_argument("-o", "--outdir", help="Path to the output directory", type=str, required=True)

    # PARAMETERS
    parser.add_argument(
        "-mc",
        "--max-chunk-size",
        help="Maximum chunk size used for final output files. Default: 5000",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "-pc",
        "--proc-chunk-size",
        help="Chunk size used for internal processing to limit memory usage. Default: 500000",
        type=int,
        default=500000,
    )
    parser.add_argument(
        "-mlld",
        "--mem_log_level_deep",
        help="Enable 'deep' option for pandas memory usage output ('deep' enables accurate usage values, but increases runtime). Default: None. ",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-mcn",
        "--maximum_chunk_number",
        help="Maximum number of generated chunks (resulting prediction processes). Default: 1000",
        type=int,
        default=1000,
    )

    return parser.parse_args()


def write_chunks(data, alleles, max_task_per_allele, remainder=False, pbar=None):
    """Takes data in form of a table of peptide_id, peptide_sequence and
    identical allele_name values. The data is partitioned into chunks and
    written into individual output files, prepended with a comment line (#)
    indicating the allele name."""
    global cur_chunk

    max_chunk_size = args.max_chunk_size

    # Dynamically increase the chunk size dependent on the maximum number of allowed processes.
    if len(data)/max_chunk_size > max_task_per_allele:
        print("WARN: Chunk size is too small and too many chunks are generated. Chunksize is increased to match the maximum number of chunks.")
        max_chunk_size = int(len(data)/max_task_per_allele)+1 # Make sure that all peptides end up in chunks

    if remainder and len(data) > max_chunk_size:
        print("ERROR: Something went wrong!", file=sys.stderr)
        sys.exit(1)

    allele_name = alleles[alleles["allele_id"] == data.iloc[0].allele_id]["allele_name"].iloc[0]
    written = pd.Index([])
    for start in range(0, len(data), max_chunk_size):
        # if not handling remainder: only write out full chunks here
        if remainder or len(data) - start >= max_chunk_size:
            with open(os.path.join(args.outdir, "peptides_" + str(cur_chunk).rjust(5, "0") + ".txt"), "w") as outfile:
                print(f"#{allele_name}#{data.iloc[0].allele_id}", file=outfile)
                write = data.iloc[start : start + max_chunk_size]
                written = written.append(data.index[start : start + max_chunk_size])
                if pbar:
                    pbar.update(len(write))
                write[["peptide_id", "peptide_sequence"]].to_csv(outfile, sep="\t", index=False)
                cur_chunk = cur_chunk + 1

    # delete chunks that were written out already
    data.drop(written, inplace=True)
    return data

####################################################################################################

try:
    # Parse command line arguments
    args = parse_args()

    # Read input files
    # NOTE try out if datatable package can be used and would be faster or more memory efficient
    # downcast df columns that will not be used as indices to save mem usage
    # (skip index columns to avoid upcasting with set_index() and increased runtime)
    protein_peptide_occs = (
        pd.read_csv(args.protein_peptide_occ, usecols=["protein_id", "peptide_id"], sep="\t")
        .set_index("peptide_id")
        .sort_index()
    )  # NOTE could this be handled more efficiently somehow (easily)?

    entities_proteins_occs = pd.read_csv(args.entities_proteins_occ, sep="\t")
    entities_proteins_occs["entity_id"] = pd.to_numeric(entities_proteins_occs["entity_id"], downcast="unsigned")

    microbiomes_entities_occs = pd.read_csv(
        args.microbiomes_entities_occ, usecols=["microbiome_id", "entity_id"], sep="\t"
    )
    microbiomes_entities_occs["microbiome_id"] = pd.to_numeric(
        microbiomes_entities_occs["microbiome_id"], downcast="unsigned"
    )
    microbiomes_entities_occs["entity_id"] = pd.to_numeric(microbiomes_entities_occs["entity_id"], downcast="unsigned")

    conditions = pd.read_csv(args.conditions, usecols=["condition_id", "microbiome_id"], sep="\t")
    conditions["condition_id"] = pd.to_numeric(conditions["condition_id"], downcast="unsigned")
    conditions["microbiome_id"] = pd.to_numeric(conditions["microbiome_id"], downcast="unsigned")

    condition_allele_map = pd.read_csv(args.condition_allele_map, sep="\t")
    condition_allele_map["condition_id"] = pd.to_numeric(condition_allele_map["condition_id"], downcast="unsigned")
    condition_allele_map["allele_id"] = pd.to_numeric(condition_allele_map["allele_id"], downcast="unsigned")

    alleles = pd.read_csv(args.alleles, sep="\t")
    alleles["allele_id"] = pd.to_numeric(alleles["allele_id"], downcast="unsigned")
    # not converting allele_name to categorical type, since it is anyway not used for larger dfs

    # Print info including memory usage for input DataFrames
    # (Pandas df ~ 3x the size of input data)
    if args.mem_log_level_deep:
        print_mem = "deep"
    else:
        print_mem = None
    print("\nInfo: protein_peptide_occs", flush=True)
    protein_peptide_occs.info(verbose=False, memory_usage=print_mem)
    print("\nInfo: entities_proteins_occs", flush=True)
    entities_proteins_occs.info(verbose=False, memory_usage=print_mem)
    print("\nInfo: microbiomes_entities_occs", flush=True)
    microbiomes_entities_occs.info(verbose=False, memory_usage=print_mem)

    # Create output directory if it doesn't exist
    if os.path.exists(args.outdir) and not os.path.isdir(args.outdir):
        print("ERROR - The target path is not a directory", file=sys.stderr)
        sys.exit(2)
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("\nJoining protein ids with allele info...", flush=True)
    # Prepare df with allele info for proteins for downstream merging against all peptides
    proteins_allele_info = (
        entities_proteins_occs.merge(microbiomes_entities_occs)
        .drop(columns="entity_id")
        .merge(conditions)
        .drop(columns="microbiome_id")
        .merge(condition_allele_map)
        .drop(columns="condition_id")
        .drop_duplicates()
        .set_index("protein_id")
        .sort_index()
    )
    # -> protein_id, allele_id

    print("\nInfo: proteins_allele_info", flush=True)
    proteins_allele_info.info(verbose=False, memory_usage=print_mem)

    cur_chunk = 0
    requests = 0
    keep = pd.DataFrame()

    # Define how many chunks may be created per allele
    max_task_per_allele = int(args.maximum_chunk_number/len(alleles["allele_id"])) # cut instead of round to ensure being lower than maximum

    # Process peptides chunk-wise, write out into files with max_chunk_size or keep for next chunk if remaining requests, i.e. < max_chunk_size for one allele
    with pd.read_csv(args.peptides, sep="\t", index_col="peptide_id", chunksize=args.proc_chunk_size) as reader:
        for c, peptides in enumerate(reader):
            print("\nChunk: ", c)
            print("Info: peptides", flush=True)
            peptides.info(verbose=False, memory_usage=print_mem)

            # Identify which predictions have to be computed: join peptides with allele info
            to_predict = (
                peptides.sort_index()
                .join(protein_peptide_occs)
                .reset_index()
                .set_index("protein_id")
                .sort_index()
                .join(proteins_allele_info)
                .drop_duplicates()
                .reset_index(drop=True)
            )
            # -> index, peptide_id, peptide_sequence, allele_id

            # TODO peptides can be deleted?

            print("Info: to_predict", flush=True)
            to_predict.info(verbose=False, memory_usage=print_mem)

            requests += len(to_predict)

            # Concat remaining peptides from previous round with newly processed ones
            # write the required predictions into chunks of peptide lists
            keep = (
                pd.concat([keep, to_predict], ignore_index=True)
                .groupby("allele_id", group_keys=False)
                .apply(lambda x: write_chunks(x, alleles, max_task_per_allele))
            )
            # use group_keys=False to avoid generation of extra index with "allele_id"

            print("Info: keep", flush=True)
            keep.info(verbose=False, memory_usage=print_mem)

    # Write out remaining peptides
    keep.groupby("allele_id", group_keys=False).apply(lambda x: write_chunks(x, alleles, max_task_per_allele, remainder=True))

    # We're happy if we got here
    print(f"All done. Written {requests} peptide prediction requests into {cur_chunk} chunks.")
    sys.exit(0)
except KeyboardInterrupt:
    print("\nUser aborted.", file=sys.stderr)
    sys.exit(1)
