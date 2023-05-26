#!/usr/bin/env python

# Written by Sabrina Krakau, Christopher Mohr and released under the MIT license.
# This script originates from the nf-core/epitopeprediction pipeline and is modified for use in nf-core/metapep

import argparse

from epytope.EpitopePrediction import EpitopePredictorFactory


def convert_allele_back(allele):
    if str(allele).startswith("H-2-"):
        # convert internal Epytope representation back to the nf-core/metapep input allele format
        return allele.replace("H-2-", "H2-")
    elif allele.startswith("HLA-"):
        return allele.replace("HLA-", "")
    else:
        raise ValueError(
            "Allele type unknown: " + allele + ". Currently expects allele to start either with 'HLA-' or 'H-2-'."
        )


def __main__():
    parser = argparse.ArgumentParser(
        "Write out information about supported models by Epytope for available prediction tool versions."
    )
    parser.add_argument(
        "-m",
        "--pred_methods",
        help="List of prediction models (sorted like list of --pred_method_versions)",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--pred_method_versions",
        help="List of prediction method versions (sorted like --pred_methods)",
        nargs="+",
        required=True,
    )
    args = parser.parse_args()

    for method, version in zip(args.pred_methods, args.pred_method_versions):
        if version not in EpitopePredictorFactory.available_methods()[method]:
            raise ValueError("The specified version " + version + " for " + method + " is not supported by Epytope.")

        predictor = EpitopePredictorFactory(method, version=version)
        with open(method + ".v" + str(version) + ".supported_alleles.txt", "w") as output:
            for a in sorted(predictor.supportedAlleles):
                output.write(convert_allele_back(a) + "\n")
        with open(method + ".v" + str(version) + ".supported_lengths.txt", "w") as output:
            for l in sorted(predictor.supportedLength):
                output.write(str(l) + "\n")


if __name__ == "__main__":
    __main__()
