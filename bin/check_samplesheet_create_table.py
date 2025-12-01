#!/usr/bin/env python
# Written by Sabrina Krakau, Leon Kuchenbecker, and Till Englert under the MIT license

import argparse
import sys

import pandas as pd
import mhcgnomes

def parse_args(args=None):
    description = "Reformat nf-core/metapep samplesheet file, check its contents and create the data tables."
    epilog = "Example usage: python check_samplesheet_create_tables.py -i <FILE_IN> -m <MICROBIOMES_OUT> -c <CONDITIONS_OUT> -a <ALLELES_OUT> -ca <CONDITION_ALLELES_OUT> -pm <prediction_method> -pmv <prediction_method_version>"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="FILE",
        type=argparse.FileType("r"),
        help="Input samplesheet file containing: condition, type, microbiome_path, alleles, weights_path.",
    )
    parser.add_argument(
        "-m",
        "--microbiomes",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: microbiome_id, microbiome_path, microbiome_type, weights_path.",
    )
    parser.add_argument(
        "-c",
        "--conditions",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: condition_id, condition_name, microbiome_id.",
    )
    parser.add_argument(
        "-a",
        "--alleles",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: allele_id, allele_name.",
    )
    parser.add_argument(
        "-ca",
        "--conditions_alleles",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: condition_id, allele_id.",
    )
    parser.add_argument(
        "-pm",
        "--prediction_method",
        required=True,
        metavar="STRING",
        type=str,
        help="Chosen method for epitope prediction",
    )
    parser.add_argument(
        "-pl",
        "--peptide_lengths",
        required=True,
        metavar="TUPLE",
        nargs="+",
        help="Peptide lengths as given in parameters (min max)",
    )
    #parser.add_argument(
    #    "-sa",
     #   "--supported_alleles",
      #  required=True,
       # metavar="FILE",
        #type=str,
        #help="Path to supported_alleles.json file",
    #)
    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)

#alleles format check; after this function all db_tables use the same allele names
# alleles format check; after this function all db_tables use the same allele names
def normalize_allele(allele: str) -> str:
    """Return canonical allele string using mhcgnomes; exit with error if parsing fails."""
    allele = (allele or "").strip()

    if not allele:
        sys.exit("ERROR: Empty allele value encountered — please check your input samplesheet.")
        
    parsed = mhcgnomes.parse(allele)

    if parsed:
        return parsed.to_string()

    sys.exit(
        "ERROR: Could not parse allele "
        f"'{allele}' — invalid or unrecognized allele format.\n"
        "Further information on supported alleles and peptide lengths:\n"
        "  nextflow run metapep -profile <YOURPROFILE> --outdir <OUTDIR> --show_supported_models"
    )

#format check after mhcgenomes and check, if allele is supported by pred_method
#same as process as prepare_prediction_input
def validate_alleles_and_lengths(unique_alleles, prediction_method, peptide_lengths, supported_alleles_json):
    """
    Validate that alleles and peptide lengths are supported by the prediction method.
    """
    # Load supported alleles
    with open(supported_alleles_json, 'r') as f:
        supported_alleles_dict = json.load(f)
    
    supported = set(supported_alleles_dict[prediction_method])
    
    # Check alleles
    unsupported_alleles = []
    for allele in unique_alleles:
        if allele not in supported:
            unsupported_alleles.append(allele)
    
    if unsupported_alleles:
        sys.exit(
            f"\n\nERROR: The following alleles are not supported by {prediction_method}:\n"
            f"  {', '.join(unsupported_alleles)}\n\n"
            f"Supported alleles for {prediction_method}: {len(supported)} alleles\n"
            f"Run with '--show_supported_models' to see all supported alleles.\n"
        )
    
    # Validate peptide lengths
    min_len, max_len = peptide_lengths
    
    # Define tool-specific length constraints
    length_constraints = {
        'mhcflurry': (8, 15),
        'mhcnuggets': (8, 15),
        'mhcnuggetsii': (9, 30),
    }
    
    if prediction_method in length_constraints:
        tool_min, tool_max = length_constraints[prediction_method]
        if min_len < tool_min or max_len > tool_max:
            sys.exit(
                f"\n\nERROR: {prediction_method} supports peptide lengths {tool_min}-{tool_max}, "
                f"but you specified {min_len}-{max_len}\n"
            )
    
    print(f" All {len(unique_alleles)} alleles are suported by {prediction_method}")
    print(f" Peptide lengths {min_len}-{max_len} are valid for {prediction_method}")


def process_samplesheet(args):
    """
    Check that the tabular samplesheet has the structure expected by nf-core/metapep and create the data tables.

    Header structure:
    condition,type,microbiome_path,alleles,weights_path

    For an example see:
    https://github.com/nf-core/metapep/raw/dev/assets/samplesheet.csv
    """
    input_table = pd.read_csv(args.input)
    input_table_cp = input_table.copy()

    # check if microbiome_path file extensions are valid
    for type, fname in zip(input_table["type"], input_table["microbiome_path"]):
        if type == "proteins":
            print_error(
                "Invalid type '"
                + type
                + "' specified in "
                + args.input.name
                + ". Type 'proteins' is not yet supported! Valid types are 'taxa', 'bins' and 'assembly'."
            )
        if type not in ["taxa", "assembly", "bins"]:
            print_error(
                "Invalid type '"
                + type
                + "' specified in "
                + args.input.name
                + ". Valid types are 'taxa', 'bins' and 'assembly'."
            )
        if type == "taxa" and not fname.lower().endswith(".tsv"):
            print_error(
                "In "
                + args.input.name
                + " specified file "
                + fname
                + " of type 'taxa' has invalid file extension. Valid extensions is '.tsv'."
            )
        if type == "proteins" and not fname.lower().endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            print_error(
                "In "
                + args.input.name
                + " specified file "
                + fname
                + " of type 'proteins' has invalid file extension. Valid extensions are '.fa', '.fa.gz', '.fasta' and"
                " '.fasta.gz'."
            )
        if type == "assembly" and not fname.lower().endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            print_error(
                "In "
                + args.input.name
                + " specified file "
                + fname
                + " of type 'assembly' has invalid file extension. Valid extensions are '.fa', '.fa.gz', '.fasta' and"
                " '.fasta.gz'."
            )

    # check if condition names unique
    if len(input_table["condition"]) != len(input_table["condition"].drop_duplicates()):
        sys.exit("Input file " + args.input.name + " contains duplicated conditions! Please specify unique conditions.")

    # check if weight_path is valid
    for type, weights_path in zip(input_table["type"], input_table["weights_path"]):
        if not (type == "assembly" or type == "bins") and not pd.isnull(weights_path):
            sys.exit(
                "Input file "
                + args.input.name
                + " contains 'weights_path' specified for type '"
                + type
                + "'! Currently input weights are only supported for type 'assembly' or 'bins."
            )
        if not pd.isnull(weights_path) and not weights_path.lower().endswith(".tsv"):
            sys.exit(
                "In "
                + args.input.name
                + " specified 'weights_path' "
                + weights_path
                + " has invalid file extension. The extension must be '.tsv'."
            )

    # microbiome_id - microbiome_path - microbiome_type
    microbiomes = (
        input_table[["microbiome_path", "type", "weights_path"]]
        .drop_duplicates()
        .rename({"type": "microbiome_type"}, axis=1)
    )
    microbiomes["microbiome_id"] = range(len(microbiomes))

    # Create bare id for each microbiome path to reduce redundancy in protein generation
    mcrb_uni = {}
    x = 0
    for path in microbiomes["microbiome_path"]:
        if path not in mcrb_uni:
            mcrb_uni[path] = x
            x += 1
        else:
            continue

    microbiomes["microbiome_bare_id"] = [mcrb_uni[path] for path in microbiomes["microbiome_path"]]

    microbiomes[["microbiome_id", "microbiome_path", "microbiome_type", "weights_path", "microbiome_bare_id"]].to_csv(
        args.microbiomes, sep="\t", index=False
    )

    # condition id - condition name - microbiome id
    conditions = input_table.merge(microbiomes)[["condition", "microbiome_id"]].rename(
        {"condition": "condition_name"}, axis=1
    )  # conditions unique (checked in nextflow)
    conditions["condition_id"] = range(len(conditions))

    conditions[["condition_id", "condition_name", "microbiome_id"]].to_csv(args.conditions, sep="\t", index=False)

    # allele id - allele name
    raw_alleles = [
        allele
        for allele_list in input_table["alleles"].astype(str)
        for allele in allele_list.split()
    ]
    normalized_alleles = [normalize_allele(a) for a in raw_alleles if a]
    unique_alleles = sorted(set(normalized_alleles))


    # validate alleles and peptide lengths
   # validate_alleles_and_lengths(
    #    unique_alleles,
     #   args.prediction_method,
      #  args.peptide_lengths,
       # args.supported_alleles
    #)


    alleles = pd.DataFrame({"allele_name": sorted(list(unique_alleles))})
    alleles["allele_id"] = range(len(alleles))
    alleles[["allele_id", "allele_name"]].to_csv(args.alleles, sep="\t", index=False)

    # condition id - allele id
    conditions_alleles = pd.DataFrame(
        [
            (row["condition"], normalize_allele(allele_name))
            for _, row in input_table.iterrows()
            for allele_name in row["alleles"].split(" ")
        ],
        columns=["condition_name", "allele_name"],
    )
    conditions_alleles = conditions_alleles.merge(conditions).merge(alleles)[["condition_id", "allele_id"]]
    conditions_alleles.to_csv(args.conditions_alleles, sep="\t", index=False)

    input_table_cp.to_csv("samplesheet.valid.csv", index=False)
    print("Done!")


def main(args=None):
    args = parse_args(args)
    process_samplesheet(args)


if __name__ == "__main__":
    sys.exit(main())
