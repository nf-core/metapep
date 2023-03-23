#!/usr/bin/env python3
#
# Author: Leon Kuchenbecker <leon.kuchenbecker@uni-tuebingen.de>
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
# epytope_predict.py
#
# This program provides a command line interface to the epitope prediction module of the 'epytope'
# python framework.
####################################################################################################

import argparse
import sys
import csv
import warnings
import contextlib
import io

####################################################################################################

import logging

logging.basicConfig(format="%(levelname)s - %(message)s", level=logging.WARNING)

####################################################################################################

import pandas as pd
import numpy as np
from epytope.Core import Allele, Peptide
from epytope.EpitopePrediction import EpitopePredictorFactory

####################################################################################################


class AlleleParseException(RuntimeError):
    """Represents a failure to parse an allele string"""

    pass


class PredictorCreationException(RuntimeError):
    """Represents a failure to create an epitope predictor"""

    pass


class PeptidesParseException(RuntimeError):
    """Represents a failure to parse the provided peptides file"""

    pass


####################################################################################################


@contextlib.contextmanager
def capture_stdout(target):
    """Captures stdout into specified target file object"""
    old = sys.stdout
    try:
        sys.stdout = target
        yield
    finally:
        sys.stdout = old


def list_available_methods(out=sys.stdout):
    """Print a list of available epitope prediction methods"""
    print("The following methods and method versions are available:", file=out)
    for name, version in EpitopePredictorFactory.available_methods().items():
        print(name.ljust(50, "."), "[" + ", ".join(version) + "]", sep="", file=out)


def get_predictor(method, version):
    """Tries to obtain a predictor based on a method and version descriptor, wraps errors"""
    try:
        return EpitopePredictorFactory(method, version=version) if version else EpitopePredictorFactory(method)
    except Exception as e:
        raise PredictorCreationException(str(e))


def validate_method(method, version):
    """Confirm that the specified method / version combination is available and exit accordingly"""
    try:
        get_predictor(method, version)
    except PredictorCreationException:
        sys.exit(1)
    sys.exit(0)


def allele_from_string(s):
    """Obtains an epytope Allele object from a string, handling errors"""
    try:
        return Allele(s)
    except:
        raise AlleleParseException(s)


def fail(s, exitcode):
    """Log an error and exit the program with the specified exit code"""
    logging.error(s)
    sys.exit(exitcode)


def parse_args():
    """Parses the command line arguments specified by the user."""
    parser = argparse.ArgumentParser(description="Perform epitope prediction using the epytope tool bindings.")
    parser.add_argument(
        "-p",
        "--peptides",
        help=(
            "Path to file with input peptides. Two input "
            "formats are supported, a plain text file with one peptide per line or a plain text tsv "
            "file with two named columns, 'peptide_id' and 'peptide_sequence'. Default: stdin."
        ),
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Path to the output file. The output file will be a "
            "plain text tsv file with a column containing either the peptide id or the peptide "
            "sequence (depending on whether ids were provided in the input), a column with the "
            "method name and a column for each allele containing the corresponding score. "
            "Default: stdout."
        ),
        type=argparse.FileType("w"),
        default=sys.stdout,
    )
    parser.add_argument(
        "-m", "--method", help="Prediction method to use. Default: syfpeithi.", type=str, default="syfpeithi"
    )
    parser.add_argument(
        "-mv",
        "--method_version",
        help="Prediction method version to use. Default: Use method-specific default version.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "-sn", "--syfpeithi-norm", help="When using the SYFPEITHI method, normalize the scores", action="store_true"
    )
    parser.add_argument(
        "-lm", "--list_methods", help="List available methods and versions and exit.", action="store_true"
    )
    parser.add_argument(
        "-vm",
        "--validate_method",
        help="Confirm that the requested method / version is available.",
        action="store_true",
    )

    parser.add_argument(
        "allele",
        nargs="*",
        help="Space separated list of MHC alleles to predict binding scores for.",
        type=str,
        default=["A*01:01"],
    )

    return parser.parse_args()


def read_peptides(f):
    """Reads the peptides from the provided file. The provided file can either be a plain text file
    with one line per peptide or a plain text tsv file with two columns, `peptide_id` and
    `peptide_sequence`."""
    try:
        peptides = pd.read_csv(f, sep="\t", header=None, comment="#")

        # Guess whether this data comes with IDs
        nrows, ncols = peptides.shape
        if ncols == 1:
            return None, [Peptide(x) for x in peptides.iloc[:, 0]]
        else:
            peptides.columns = peptides.iloc[0]
            peptides = peptides.iloc[1:]
            return list(peptides["peptide_id"]), [Peptide(x) for x in peptides["peptide_sequence"]]
    except KeyError as e:
        raise PeptidesParseException(
            f"Missing column: {str(e)}. Please provide either a plain list or a table containing named columns"
            " 'peptide_id' and 'peptide_sequence'."
        )
    except pd.errors.EmptyDataError:
        raise PeptidesParseException("The file appears to be empty.")


def write_results(outfile, ids, peptides, predictions):
    """Write the prediction results to the specified output file. If the input peptide sequences
    were annotated with ids, the ids are written into the output table instead of the sequences."""
    # Remove the index and rename the columns of the prediction results.
    predictions.rename({"Seq": "peptide_sequence", "Method": "method"}, axis=1, inplace=True)
    # Check if we have ids from the provided input data and write the results accordingly.
    if ids:
        id_map = pd.DataFrame({"peptide_id": ids, "peptide_sequence": peptides})
        id_map.merge(predictions).drop("peptide_sequence", axis=1).to_csv(outfile, sep="\t", index=False, na_rep="NA")
    else:
        predictions.to_csv(outfile, sep="\t", index=False, na_rep="NA")


def matrix_max(matrix):  # SYFPEITHI NORMALIZATION
    """Returns the maximum attainable score for a pssm"""
    return sum([max(value.values()) for _, value in matrix.items()])


def get_allele_model_max_value(allele, length):  # SYFPEITHI NORMALIZATION
    """Returns the SYFPEITHI pssm max value for a given allele"""
    allele_model = "%s_%i" % (allele, length)
    try:
        return matrix_max(
            getattr(
                __import__("epytope.Data.pssms.syfpeithi" + ".mat." + allele_model, fromlist=[allele_model]), allele_model
            )
        )
    except ImportError:
        return None


def syfpeithi_normalize(predictions):  # SYFPEITHI NORMALIZATION
    """Normalizes syfpeithi prediction scores by dividing by the maximum
    attainable score for a particular model"""
    predictor = EpitopePredictorFactory("syfpeithi")
    alleles = [cname for cname in predictions.columns if cname not in ["Seq", "Method"]]
    conv_alleles = predictor.convert_alleles(alleles)

    lengths = predictions.Seq.apply(len)
    for allele, conv_allele in zip(alleles, conv_alleles):
        scored_lengths = {l for l, score in zip(lengths, predictions[allele]) if not pd.isna(score)}
        max_vals = {l: get_allele_model_max_value(conv_allele, l) for l in scored_lengths}
        new_scores = [
            np.nan if pd.isna(score) else score / max_vals[l] for score, l in zip(predictions[allele], lengths)
        ]
        predictions[allele] = new_scores

    return predictions


####################################################################################################

try:
    # Parse command line arguments
    args = parse_args()

    # Validate specified method if requested
    if args.validate_method:
        validate_method(method=args.method, version=args.method_version)

    # Print methods if requested
    if args.list_methods:
        list_available_methods()
        sys.exit(0)

    # Parse allele names
    alleles = [allele_from_string(allele) for allele in args.allele]

    # Read peptides from provided input file
    ids, peptides = read_peptides(args.peptides)

    # Create predictor
    predictor = get_predictor(method=args.method, version=args.method_version)

    # Run predictor
    predictor_stdout = io.StringIO()
    query_peptides_index = pd.DataFrame({"Seq": peptides}).Seq
    with warnings.catch_warnings(record=True) as warnings:
        try:
            # Redirect stdout output of predictor code to stderr to use stdout
            # exclusively for the results table if writing to stdout was
            # specified by the user.
            with capture_stdout(sys.stderr):
                predictions = (pd.concat(
                    [pred[["Peptides", "Method", 0]].rename(columns={"Peptides":"Seq", 0:allel}) for allel, pred in predictor.predict(peptides, alleles=alleles).unstack().reset_index().groupby("Allele")[["Peptides", "Method", 0]]])
                )
            predictions["Method"].fillna(args.method, inplace=True)
        except ValueError as e:
            predictions = pd.DataFrame({"Seq": query_peptides_index, "Method": args.method})
        for message in {w.message for w in warnings}:
            logging.warning(f"PREDICTION ({predictor.name} {predictor.version}) - {str(message)}")

    # Add missing alleles as NA values
    for missing_allele in [str(allele) for allele in alleles if str(allele) not in predictions.columns]:
        logging.warning(
            f"PREDICTION ({predictor.name} {predictor.version}) yielded no results for allele {missing_allele}"
        )
        predictions[missing_allele] = np.NaN

    # Normalize syfpeithi scores
    if args.method == "syfpeithi" and args.syfpeithi_norm:
        predictions = syfpeithi_normalize(predictions)

    # Write results
    write_results(args.output, ids, peptides, predictions)

    sys.exit(0)
except KeyboardInterrupt:
    fail(f"User requested shutdown.", 1)
except AlleleParseException as e:
    fail(f"Failed to parse allele specification '{e}'.", 2)
except PredictorCreationException as e:
    fail(f"Failed to find a predictor based on your specification ('{args.method}' version '{args.method_version}')", 3)
except PeptidesParseException as e:
    fail(f"The provided input file could not be parsed: {str(e)}", 4)
except Exception as e:
    fail(f"UNEXPECTED - {e}", 999)
