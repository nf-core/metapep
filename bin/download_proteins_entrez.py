#!/usr/bin/env python3
# Written by Sabrina Krakau, Leon Kuchenbecker, and Till Englert under the MIT license

# for each strain: select largest assembly or given specific (for now)

import argparse
import csv
import gzip
import sys
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from urllib.error import HTTPError

from Bio import Entrez, SeqIO

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t",
        "--taxid_input",
        required=True,
        nargs="+",
        metavar="FILE",
        type=argparse.FileType("r"),
        help="List of microbiome files containing: taxon_id [,assembly_id, abundance].",
    )
    parser.add_argument(
        "-m", "--microbiome_ids", required=True, nargs="+", help="List of corresponding microbiome IDs."
    )
    parser.add_argument("-e", "--email", required=True, help="Email address to use for NCBI access.")
    parser.add_argument("-k", "--key", required=True, help="NCBI key to allow faster access.")
    parser.add_argument(
        "-p",
        "--proteins",
        required=True,
        metavar="FILE",
        help="Output file (compressed) containing: protein_tmp_id, protein_sequence.",
    )
    parser.add_argument(
        "-ta",
        "--taxa_assemblies",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: taxon_id, assembly_id.",
    )
    parser.add_argument(
        "-ep",
        "--entities_proteins",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: protein_tmp_id, entity_name (taxon_id).",
    )
    parser.add_argument(
        "-me",
        "--microbiomes_entities",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing: entity_name (taxon_id), microbiome_id, entity_weight.",
    )
    return parser.parse_args(args)


# get assembly length ("total_length") from entrez
def get_assembly_length(assemblyId):
    success = False
    for attempt in range(3):
        try:
            with Entrez.efetch(db="assembly", rettype="docsum", id=assemblyId) as entrez_handle:
                assembly_stats = Entrez.read(entrez_handle, validate=False)
            time.sleep(1)  # avoid getting blocked by ncbi
            success = True
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(10)
            else:
                raise
    if not success:
        sys.exit("Entrez efetch download failed!")

    root = ET.fromstring("<root>" + str(assembly_stats["DocumentSummarySet"]["DocumentSummary"][0]["Meta"]) + "</root>")
    return int(root.find("./Stats/Stat[@category='total_length'][@sequence_tag='all']").text)


def main(args=None):
    args = parse_args(args)

    # setup entrez email
    Entrez.api_key = args.key
    Entrez.email = args.email

    # read taxonomic ids for download (together with abundances) and write 'microbiomes_entities' output
    taxIds = []
    input_taxids_assemblyids = {}
    print("entity_name", "microbiome_id", "entity_weight", sep="\t", file=args.microbiomes_entities)
    for taxid_input, microbiomeId in zip(args.taxid_input, args.microbiome_ids):
        reader = csv.DictReader(taxid_input, delimiter="\t")
        for row in reader:
            # If the abundance is not defined in the taxid_input it is assigned as 1.
            if row.keys() == set(["taxon_id"]):
                taxIds.append(row["taxon_id"])
                print(row["taxon_id"], microbiomeId, 1, sep="\t", file=args.microbiomes_entities, flush=True)

            elif row.keys() == set(["taxon_id", "abundance"]):
                taxIds.append(row["taxon_id"])
                print(
                    row["taxon_id"],
                    microbiomeId,
                    row["abundance"],
                    sep="\t",
                    file=args.microbiomes_entities,
                    flush=True,
                )

            elif row.keys() == set(["taxon_id", "assembly_id", "abundance"]):
                taxIds.append(row["taxon_id"])
                if row["assembly_id"] != "":
                    input_taxids_assemblyids[row["taxon_id"]] = row["assembly_id"]
                print(
                    row["taxon_id"],
                    microbiomeId,
                    row["abundance"],
                    sep="\t",
                    file=args.microbiomes_entities,
                    flush=True,
                )

            elif row.keys() == set(["taxon_id", "assembly_id"]):
                taxIds.append(row["taxon_id"])
                if row["assembly_id"] != "":
                    input_taxids_assemblyids[row["taxon_id"]] = row["assembly_id"]
                print(row["taxon_id"], microbiomeId, 1, sep="\t", file=args.microbiomes_entities, flush=True)

            else:
                sys.exit(
                    f"The format of the input file '{taxid_input.name}' is invalid!"
                    + "It needs to be a tsv file containing 'taxon_id' and/or optionally 'assembly_id' and 'abundance."
                )

    taxIds = sorted(list(set(taxIds)))
    ####################################################################################################
    # Process TaxIDs

    print("Processing the following taxonomy IDs:")
    print(taxIds)

    ####################################################################################################
    # 0) Check if the taxids link to a strain level organism
    print("Check if taxonomy IDs are on strain rank:")

    success = False
    for attempt in range(3):
        try:
            with Entrez.efetch(db="taxonomy", id=taxIds) as entrez_handle:
                record = Entrez.read(entrez_handle)
                if len(record) == len(taxIds):
                    for taxid, rec in zip(taxIds, record):
                        rank = rec["Rank"]
                        if rank != "strain":
                            sys.exit(f"Strain level check failed for taxid: {taxid}. {taxid} linked to rank: {rank}")
                    time.sleep(1)  # avoid getting blocked by ncbi
                    success = True
                    break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(10)
            else:
                raise
    if not success:
        sys.exit("Entrez efetch download failed!")

    print("Taxids succeeded strain level check.")

    ####################################################################################################
    # 1) for each taxId -> get all assembly IDs // skip if assemblyID is given in input
    print("# taxa: ", len(taxIds))
    print("# taxa without assemblyId: ", len(set(taxIds)-set(input_taxids_assemblyids.keys())))
    if len(set(taxIds)) != len(set(input_taxids_assemblyids.keys())) and len(set(input_taxids_assemblyids.keys())) != 0:
        print("WARNING: A mix of specific assembly Ids and unspecific taxids was chosen!")
    print("for each taxon without assemblyID retrieve assembly IDs ...")

    taxIds_without_assemblyID = set(taxIds)-set(input_taxids_assemblyids.keys())
    success = False
    for attempt in range(3):
        try:
            with Entrez.elink(
                dbfrom="taxonomy", db="assembly", LinkName="taxonomy_assembly", id=taxIds_without_assemblyID
            ) as entrez_handle:
                assembly_results = Entrez.read(entrez_handle)
            time.sleep(1)  # avoid getting blocked by ncbi
            success = True
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(10)
            else:
                raise
    if not success:
        sys.exit("Entrez elink download failed!")

    ####################################################################################################
    # 2) for each taxon -> select one assembly (largest for now) and merge with input assembly ids
    #    Selecting the "largest" assembly is not reproducible, there are taxon_ids (e.g. 83333
    #    (E.coli K12)) which contains many duplicate assembly lengths. This results in different
    #    numbers of downloaded proteins for each assembly, in combination with other assemblies
    #    (as the overlap may be different between assemblies even though having the same size)

    print("get assembly lengths and select largest assembly for each taxon ...")
    dict_taxId_assemblyId = {}
    for tax_record in assembly_results:
        taxId = tax_record["IdList"][0]
        if len(tax_record["LinkSetDb"]) > 0:
            # get all assembly ids
            ids = [assembly_record["Id"] for assembly_record in tax_record["LinkSetDb"][0]["Link"]]
            # get corresponding lengths
            lengths = [get_assembly_length(id) for id in ids]
            # get id for largest assembly
            selected_assemblyId = ids[lengths.index(max(lengths))]
            dict_taxId_assemblyId[taxId] = selected_assemblyId
    # Merge input assembly ids with fetched assembly ids for taxids
    dict_taxId_assemblyId = dict_taxId_assemblyId | input_taxids_assemblyids
    # write taxId - assemblyId out
    print("taxon_id", "assembly_id", sep="\t", file=args.taxa_assemblies, flush=True)
    for taxId in dict_taxId_assemblyId.keys():
        print(taxId, dict_taxId_assemblyId[taxId], sep="\t", file=args.taxa_assemblies, flush=True)

    ####################################################################################################
    # 3) (selected) assembly -> nucleotide sequences

    assemblyIds = dict_taxId_assemblyId.values()
    print("# selected assemblies: ", len(assemblyIds))
    print("for each assembly get nucloetide sequence IDs...")

    success = False
    for attempt in range(3):
        try:
            with Entrez.elink(
                dbfrom="assembly", db="nuccore", LinkName="assembly_nuccore_refseq", id=assemblyIds
            ) as entrez_handle:
                nucleotide_results = Entrez.read(entrez_handle)
            time.sleep(1)  # avoid getting blocked by ncbi
            success = True
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(10)
            else:
                raise
    if not success:
        sys.exit("Entrez elink download failed!")

    ### for each assembly get list of sequence ids
    dict_seqId_assemblyIds = defaultdict(lambda: [])
    for assembly_record in nucleotide_results:
        assemblyId = assembly_record["IdList"][0]
        for record in assembly_record["LinkSetDb"][0]["Link"]:
            dict_seqId_assemblyIds[record["Id"]].append(assemblyId)

    print("# nucleotide sequences (unique): ", len(dict_seqId_assemblyIds.keys()))
    # -> # contigs

    ####################################################################################################
    # 4) nucelotide sequences -> proteins
    print("for each nucleotide sequence get proteins ...")

    success = False
    for attempt in range(3):
        try:
            with Entrez.elink(
                dbfrom="nuccore", db="protein", LinkName="nuccore_protein", id=list(dict_seqId_assemblyIds.keys())
            ) as entrez_handle:
                protein_results = Entrez.read(entrez_handle)
            time.sleep(1)  # avoid getting blocked by ncbi
            success = True
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(10)
            else:
                raise
    if not success:
        sys.exit("Entrez elink download failed!")

    ### for each nucleotide sequence get list of protein ids
    dict_proteinId_assemblyIds = defaultdict(lambda: set())
    for nucleotide_record in protein_results:
        seqId = nucleotide_record["IdList"][0]
        assemblyIds = dict_seqId_assemblyIds[seqId]
        if len(nucleotide_record["LinkSetDb"]) > 0:
            for protein_record in nucleotide_record["LinkSetDb"][0]["Link"]:
                dict_proteinId_assemblyIds[protein_record["Id"]].update(assemblyIds)

    # NOTE:
    # some proteins, such as 487413233, occur within multiple sequences of the assembly!
    # -> assembly only listed once!

    proteinIds = sorted(list(dict_proteinId_assemblyIds.keys()))

    print("# proteins (unique): ", len(proteinIds))
    # -> # proteins with refseq source (<= # IPG proteins)

    ####################################################################################################
    # 5) download protein FASTAs, convert to TSV
    print("    download proteins ...")

    # Chunking into retmax (maximum defined by E-Utilities 9999)
    prot_id_chunks = []
    retmax = 9999
    if len(proteinIds) > retmax:
        for n_chunk in range(1, int(len(proteinIds)/retmax)+1):
            chunk_start = n_chunk*retmax-retmax
            prot_id_chunks.append(proteinIds[chunk_start:chunk_start+retmax])
        prot_id_chunks.append(proteinIds[(n_chunk+1)*retmax-retmax:])
    else:
        prot_id_chunks.append(proteinIds)

    print(f"Generated {len(prot_id_chunks)} chunks and start processing:")

    protein_summaries = []
    # first retrieve mapping for protein UIDs and accession versions
    for prot_id_chunk in prot_id_chunks:
        print(f"Process chunk with {len(prot_id_chunk)} protein_ids:")
        success = False
        for attempt in range(3):
            try:
                with Entrez.esummary(
                    db="protein", id=",".join(prot_id_chunk)
                ) as entrez_handle:  # esummary doesn't work on python lists somehow
                    protein_summaries.extend(Entrez.read(entrez_handle))
                time.sleep(1)  # avoid getting blocked by ncbi
                success = True
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    time.sleep(10)
                else:
                    raise
        if not success:
            sys.exit("Entrez esummary download failed!")

        # download actual fasta records and write out
        success = False
        for attempt in range(3):
            try:
                with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=prot_id_chunk) as entrez_handle:
                    with gzip.open(args.proteins, "wt") as out_handle:
                        print("protein_tmp_id", "protein_sequence", sep="\t", file=out_handle)
                        for record in SeqIO.parse(entrez_handle, "fasta"):
                            print(record.id, record.seq, sep="\t", file=out_handle, flush=True)
                time.sleep(1)  # avoid getting blocked by ncbi
                success = True
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    time.sleep(10)
                else:
                    raise
        if not success:
            sys.exit("Entrez efetch download failed!")

        print(f"Downloaded a total of {len(protein_summaries)} protein summaries and their respective sequences.")

    dict_protein_uid_acc = {}
    for protein_summary in protein_summaries:
        dict_protein_uid_acc[protein_summary["Id"]] = protein_summary["AccessionVersion"]

    if len(proteinIds) != len(dict_protein_uid_acc.keys()):
        sys.exit(f"Unmatched size of downloaded protein ids ({len(dict_protein_uid_acc.keys())}) and input protein ids ({len(proteinIds)})")

    ####################################################################################################
    # 6) write out 'entities_proteins.entrez.tsv'
    print("protein_tmp_id", "entity_name", sep="\t", file=args.entities_proteins)
    dict_assemblyId_taxId = {v: k for k, v in dict_taxId_assemblyId.items()}
    if len(dict_assemblyId_taxId) != len(dict_taxId_assemblyId):
        sys.exit("Creation of dict_assemblyId_taxId failed!")

    for proteinId in proteinIds:
        accVersion = dict_protein_uid_acc[proteinId]
        # write out protein_tmp_id, entity_name (taxon_id)
        for assemblyId in sorted(dict_proteinId_assemblyIds[proteinId]):
            taxId = dict_assemblyId_taxId[assemblyId]
            print(accVersion, taxId, sep="\t", file=args.entities_proteins, flush=True)

    # NOTE for proteins the NCBI accession version ids are written out to enable a mapping to the fasta/tsv output
    # in contrast, for assemblies UIDs are written out

    print("Done!")

if __name__ == "__main__":
    sys.exit(main())
