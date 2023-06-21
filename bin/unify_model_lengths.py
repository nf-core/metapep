#!/usr/bin/env python3
#
# Author: Till Englert <till.englert@qbic.uni-tuebingen.de>
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

import sys
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
        "-o",
        "--output",
        required=True,
        help=(
            "Path to the output file. The output file will be a "
            "plain text tsv file with a column containing the allele name, a column containing the length and a column containing the corresponding allele models."
        ),
        type=argparse.FileType("w"),
        default="unified_allele_models.tsv",
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

    samplesheet = pd.read_csv(args.input)

    # Retrieve unique list of alleles
    alleles_s = {allele for allele_list in samplesheet["alleles"] for allele in allele_list.split(" ")}
    # Parse alleles to epytope convention
    predictor = EpitopePredictorFactory(args.method)
    alleles = [allele_from_string(allele) for allele in alleles_s]
    conv_alleles = predictor.convert_alleles(alleles)
    # Check if a model is available at given lengths
    allele_availability = []
    for conv_allele, allele_s in zip(conv_alleles, alleles_s):
        for pep_len in range(args.peptide_min_len, args.peptide_max_len + 1):
            model_name = f"{conv_allele}_{pep_len}"
            availibility = check_model_availability(model_name, args.method)
            allele_availability.append([allele_s, pep_len, model_name, availibility])

    allele_availability = pd.DataFrame(
        allele_availability, columns=["Allele", "Peptide_Length", "Allele_Model", "Availability"]
    )

    # Drop all non available models
    allele_availability = allele_availability[allele_availability["Availability"]]
    allele_availability.drop("Availability", axis=1, inplace=True)

    # get intersection of peptide lengths
    len_sets = [set(allele_models["Peptide_Length"]) for allele, allele_models in allele_availability.groupby("Allele")]
    len_intersect = set.intersection(*len_sets)

    # Remove all non fitting lengths
    allele_availability = allele_availability[allele_availability["Peptide_Length"].isin(len_intersect)]

    allele_availability.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
