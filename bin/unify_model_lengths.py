#!/usr/bin/env python3


import argparse
import pandas as pd
from epytope.Core import Allele
from epytope.EpitopePrediction import EpitopePredictorFactory


class AlleleParseException(RuntimeError):
    """Represents a failure to parse an allele string"""

    pass


def parse_args():
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(
        description="Checks if the given allele and length is supported by the given peptide prediction tool and unifies to a common denominator"
    )
    parser.add_argument(
        "-i", "--input", required=True, help=("Path to the input samplesheet file, containing the allele names.")
    )
    parser.add_argument(
        "-m", "--method", help="Prediction method to use. Default: syfpeithi.", type=str, default="syfpeithi"
    )
    parser.add_argument(
        "-pll", "--peptide_min_len", help="Minimum length of peptide used for prediction", type=int, default=8
    )
    parser.add_argument(
        "-plh", "--peptide_max_len", help="Maximum length of peptide used for prediction", type=int, default=13
    )
    parser.add_argument(
        "-log_s",
        "--output_log_suffix",
        required=True,
        help="Suffix for the output log file.",
        type=str,
        default="_unify_peptide_lengths",
    )
    return parser.parse_args()


def allele_from_string(allele_string):
    """Obtains an epytope Allele object from a string, handling errors"""
    try:
        return Allele(allele_string)
    except:
        raise AlleleParseException(allele_string)


def check_model_availability(model_name, prediction_method):
    try:
        __import__("epytope.Data.pssms." + prediction_method + ".mat." + model_name, fromlist=[model_name])
        return True
    except ImportError:
        return False


def main():
    args = parse_args()
    log_str = ""

    log_str += "### Unify Peptide Lengths ###\n\n"

    samplesheet = pd.read_csv(args.input)

    # Retrieve unique list of alleles
    alleles_s = {allele for allele_list in samplesheet["alleles"] for allele in allele_list.split(" ")}
    log_str += f"Found the following alleles: {', '.join(alleles_s)}\n\n"
    # Parse alleles to epytope convention
    predictor = EpitopePredictorFactory(args.method)
    alleles = [allele_from_string(allele) for allele in alleles_s]
    conv_alleles = predictor.convert_alleles(alleles)
    # Check if a model is available at given lengths
    input_lengths = [i for i in range(args.peptide_min_len, args.peptide_max_len + 1)]
    log_str += f"Check if models are available at given lengths: {', '.join(map(str, input_lengths))}\n\n"
    allele_availability = []
    for conv_allele, allele_s in zip(conv_alleles, alleles_s):
        for pep_len in input_lengths:
            model_name = f"{conv_allele}_{pep_len}"
            availability = check_model_availability(model_name, args.method)
            if availability:
                log_str += f"Found model for allele {allele_s} with length {pep_len}\n"
            else:
                log_str += f"No model found for allele {allele_s} with length {pep_len}\n"
            allele_availability.append([allele_s, pep_len, model_name, availability])

    allele_availability = pd.DataFrame(
        allele_availability, columns=["Allele", "Peptide_Length", "Allele_Model", "Availability"]
    )

    # Drop all non available models
    allele_availability = allele_availability[allele_availability["Availability"]]
    allele_availability.drop("Availability", axis=1, inplace=True)

    # get intersection of peptide lengths
    len_sets = [set(allele_models["Peptide_Length"]) for allele, allele_models in allele_availability.groupby("Allele")]
    len_intersect = set.intersection(*len_sets)

    log_str += "\nReducing the used peptide lengths to the common denominator\n"
    log_str += f"Following lengths are used for the epitope prediction on the alleles {', '.join(alleles_s)}: {', '.join(map(str, len_intersect))}\n"
    log_str += "All other peptide lengths are omitted from further analysis."
    # Remove all non fitting lengths
    allele_availability = allele_availability[allele_availability["Peptide_Length"].isin(len_intersect)]

    if len_intersect == set(input_lengths):
        log_prefix = "SUCCEEDED"
    elif len(len_intersect) == 0:
        log_prefix = "ERROR"
    else:
        log_prefix = "WARNING"
    log_fname = log_prefix + args.output_log_suffix + ".log"

    with open(log_fname, "w") as log:
        log.write(log_str)

    # Output the unified lengths to the stdout, so it can be used as output for the process
    print(",".join(str(l) for l in len_intersect))


if __name__ == "__main__":
    main()
