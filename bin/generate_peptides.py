#!/usr/bin/env python3

import argparse
import gzip
import sys

import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--proteins",
        required=True,
        metavar="FILE",
        help="Compressed TSV file containing: protein_id, protein_sequence.",
    )
    parser.add_argument(
        "-p", "--peptides", required=True, metavar="FILE", help="Output file containing: peptide_id, peptide_sequence."
    )  # use str type to allow compression of output
    parser.add_argument(
        "-pp",
        "--proteins_peptides",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: protein_id, peptide_id, count.",
    )
    parser.add_argument(
        "-l",
        "--proteins_lengths",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: protein_id, protein_length.",
    )
    parser.add_argument(
        "-mlld",
        "--mem_log_level_deep",
        help="Enable 'deep' option for pandas memory usage output ('deep' enables accurate usage values, but increases runtime). Default: None. ",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "-pll", "--peptide_lengths", required=True, metavar="LIST", nargs="*", help="Peptide lengths as list."
    )
    return parser.parse_args(args)


# Validate letters of input protein sequences to avoid unnoticed loss of input k-mers
def validate_letters(string, alphabet):
    for letter in string:
        if letter not in alphabet:
            print("ERROR: invalid input letter ", letter, ". The supported alphabet is ", alphabet, ".")
            sys.exit(1)


def gen_peptides(prot_seq, k, prefix):
    return [prot_seq[i : (i + k)] for i in range(len(prot_seq) - k + 1) if prot_seq[i] == prefix]


def main(args=None):
    args = parse_args(args)
    if args.mem_log_level_deep:
        print_mem = "deep"
    else:
        print_mem = None

    # Generate peptides in chunks based on initial AA to reduce memory usage
    # valid amino acid codes:
    # 20 standard ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
    # and extended codes ('B', 'J', 'O', 'U', 'X', 'Z')
    aa_list = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
        "B",
        "J",
        "O",
        "U",
        "X",
        "Z",
    ]

    protid_protseq_protlen = pd.read_csv(args.proteins, sep="\t")
    # downcast df columns where possible (i.e. that will not be used as index for downstream joining)
    protid_protseq_protlen["protein_id"] = pd.to_numeric(protid_protseq_protlen["protein_id"], downcast="unsigned")
    # validate input AAs
    protid_protseq_protlen["protein_sequence"] = protid_protseq_protlen["protein_sequence"].str.upper()
    protid_protseq_protlen["protein_sequence"].apply(validate_letters, alphabet=aa_list)
    # get protein lengths
    protid_protseq_protlen["protein_length"] = protid_protseq_protlen["protein_sequence"].apply(len)
    protid_protseq_protlen["protein_length"] = pd.to_numeric(
        protid_protseq_protlen["protein_length"], downcast="unsigned"
    )

    print("\nInfo: protid_protseq_protlen", flush=True)
    protid_protseq_protlen.info(verbose=False, memory_usage=print_mem)

    print("# proteins: ", len(protid_protseq_protlen))

    # write out protein lengths
    protid_protseq_protlen[["protein_id", "protein_length"]].to_csv(args.proteins_lengths, sep="\t", index=False)

    # Parse peptide_lengths input list to int
    peptide_lengths_int = [int(p_len) for p_len in args.peptide_lengths]

    ####################
    # generate peptides
    with gzip.open(args.peptides, "wt") as pep_handle:
        print_header = True
        id_counter = 0

        # for each k
        for k in peptide_lengths_int:
            print("Generate peptides of length ", k, " ...", flush=True)

            # Note: could be done with prefixes instead of single first letters if this remains bottleneck
            for prefix in aa_list:
                print("with prefix ", prefix, flush=True)
                # for each protein generate all peptides of length k with current prefix (to reduce peak mem usage)
                # (the AA-wise processing causes multiple iterations over the same protein sequences,
                # but the increase of run time is negligible in this context)
                results = pd.DataFrame(
                    [
                        (it.protein_id, pep)
                        for it in protid_protseq_protlen.itertuples()
                        for pep in gen_peptides(it.protein_sequence, k, prefix)
                    ],
                    columns=["protein_id", "peptide_sequence"],
                )

                print("\nInfo: results (['protein_id','peptide_sequence'])", flush=True)
                results.info(verbose=False, memory_usage=print_mem)

                print("format results ...", flush=True)
                # count occurrences of one peptide in one protein
                results = results.groupby(["protein_id", "peptide_sequence"]).size().reset_index(name="count")
                # -> protein_id, peptide_sequence, count
                results["count"] = pd.to_numeric(results["count"], downcast="unsigned")
                # prepare df for joining
                results.set_index("peptide_sequence", inplace=True)
                results.sort_index(inplace=True, kind="stable")

                unique_peptides = pd.DataFrame(index=results.index.drop_duplicates())
                unique_peptides["peptide_id"] = range(id_counter, id_counter + len(unique_peptides))
                id_counter += len(unique_peptides)
                # -> peptide_sequence, peptide_id
                unique_peptides.to_csv(pep_handle, mode="a", sep="\t", index=True, header=print_header)

                results = results.join(unique_peptides)
                # -> protein_id, peptide_sequence, count, peptide_id

                print("\nInfo: results (['protein_id','peptide_sequence','peptide_id','count'])", flush=True)
                results.info(verbose=False, memory_usage=print_mem)

                results[["protein_id", "peptide_id", "count"]].to_csv(
                    args.proteins_peptides, mode="a", sep="\t", index=False, header=print_header
                )

                print("# peptides of length ", k, ", (non-unique across proteins): ", len(results))
                print_header = False

    print("Done!", flush=True)


if __name__ == "__main__":
    sys.exit(main())
