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
# for each strain: select largest assembly (for now)

import sys
import gzip
import csv
import xml.etree.ElementTree as ET
import tqdm
import io
import argparse
import time

from Bio import Entrez, SeqIO
from pprint import pprint
from datetime import datetime
from collections import Counter
from urllib.error import HTTPError

import sys


# TODO
# clean code
# double check!


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', "--taxid_input", required=True, help="File containing list with taxonomy IDs.")
    parser.add_argument('-e', "--email", required=True, help="Email address to use for NCBI access.")
    parser.add_argument('-k', "--key", required=True, help="NCBI key to allow faster access.")
    return parser.parse_args(args)


# get assembly length ("total_length") from entrez
def get_assembly_length(assembly_id):
    success = False
    for attempt in range(3):
        try:
            with Entrez.efetch(db="assembly", rettype="docsum", id=assembly_id) as entrez_handle:
                assembly_stats = Entrez.read(entrez_handle, validate=False)
            time.sleep(1)   # avoid getting blocked by ncbi
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
        sys.exit("Entrez download failed!")

    root = ET.fromstring("<root>" + str(assembly_stats['DocumentSummarySet']['DocumentSummary'][0]['Meta']) + "</root>")
    return root.find("./Stats/Stat[@category='total_length'][@sequence_tag='all']").text


def main(args=None):
    args = parse_args(args)

    # setup entrez email
    Entrez.api_key = args.key
    Entrez.email = args.email

    # output filenames
    tax_assembly_id_output = "taxa_assembly.tsv"
    fasta_output = "proteins.fasta"
    prot_assembly_id_output = "proteinid_assemblyids.txt"

    with open(args.taxid_input, 'rt') as in_file:
        taxids = [ row[0] for row in csv.reader(in_file) ]

    print("Processing the following taxonmy IDs:")
    print(taxids)

    ####################################################################################################
    # 1) for each taxId -> get all assembly IDs
    print("# taxa: ", len(taxids))
    print("for each taxon retrieve assembly IDs ...")

    success = False
    for attempt in range(3):
        try:
            with Entrez.elink(dbfrom="taxonomy", db="assembly", LinkName="taxonomy_assembly", id=taxids) as entrez_handle:
                assembly_results = Entrez.read(entrez_handle)
            time.sleep(1)   # avoid getting blocked by ncbi
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
        sys.exit("Entrez download failed!")
    # print("assembly_results: ")
    # print(assembly_results)

    # 2) for each taxon -> select one assembly (largest for now)
    print("get assembly lengths and select largest assembly for each taxon ...")
    dict_taxid_assemblyids = {}
    for tax_record in assembly_results:
        taxid = tax_record["IdList"][0]
        if len(tax_record["LinkSetDb"]) > 0:
            # get all assembly ids
            ids = [ assembly_record["Id"] for assembly_record in tax_record["LinkSetDb"][0]["Link"] ]
            # get corresponding lengths
            lengths = [get_assembly_length(id) for id in ids]
            # get id for largest assembly
            selected_assembly_id = ids[lengths.index(max(lengths))]
            dict_taxid_assemblyids[taxid] = selected_assembly_id

    # print("dict")
    # print(dict_taxid_assemblyids)

    # write taxId - assemblyId out (-> results!)
    out_file = open(tax_assembly_id_output,'w')
    for taxon in dict_taxid_assemblyids.keys():
            print(taxon, dict_taxid_assemblyids[taxon], file = out_file, flush=True)
    out_file.close()

    # 3) (selected) assembly -> nucleotide sequences
    # (maybe split here)
    assembly_ids = dict_taxid_assemblyids.values()
    print("# selected assemblies: ", len(assembly_ids))
    print("for each assembly get nucloetide sequence IDs...")

    success = False
    for attempt in range(3):
        try:
            with Entrez.elink(dbfrom="assembly", db="nuccore", LinkName="assembly_nuccore_refseq", id=assembly_ids) as entrez_handle:
                nucleotide_results = Entrez.read(entrez_handle)
            time.sleep(1)   # avoid getting blocked by ncbi
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
            sys.exit("Entrez download failed!")
    # print("nucleotide_results: ")
    # print(nucleotide_results)

    ### for each assembly get list of sequence ids
    # seq_ids = [ record["Id"] for assembly_record in nucleotide_results for record in assembly_record["LinkSetDb"][0]["Link"] ]

    dict_seqid_assemblyids = {}
    for assembly_record in nucleotide_results:
        assemblyid = assembly_record["IdList"][0]
        for record in assembly_record["LinkSetDb"][0]["Link"]:
            if record["Id"] not in dict_seqid_assemblyids:
                dict_seqid_assemblyids[record["Id"]] = [assemblyid]
            else:
                dict_seqid_assemblyids[record["Id"]].append(assemblyid)

    print("# nucleotide sequences (unique): ", len(dict_seqid_assemblyids.keys()))
    # -> # contigs

    # 4) nucelotide sequences -> proteins
    print("for each nucleotide sequence get proteins ...")

    success = False
    for attempt in range(3):
        try:
            with Entrez.elink(dbfrom="nuccore", db="protein", LinkName="nuccore_protein", id=list(dict_seqid_assemblyids.keys())) as entrez_handle:
                protein_results = Entrez.read(entrez_handle)
            time.sleep(1)   # avoid getting blocked by ncbi
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
            sys.exit("Entrez download failed!")
    # print("protein_results: ")
    # print(protein_results)

    ### for each nucleotide sequence get list of protein ids
    # protein_ids = [ record["Id"] for nuc_record in protein_results for record in nuc_record["LinkSetDb"][0]["Link"] ]

    dict_proteinid_assemblyids = {}
    for nucleotide_record in protein_results:
        seqid = nucleotide_record["IdList"][0]
        assemblyids = dict_seqid_assemblyids[seqid]
        # print("assemblyids")
        # print(assemblyids)
        if len(nucleotide_record["LinkSetDb"]) > 0:
            for protein_record in nucleotide_record["LinkSetDb"][0]["Link"]:
                if protein_record["Id"] not in dict_proteinid_assemblyids:
                    dict_proteinid_assemblyids[protein_record["Id"]] = assemblyids
                else:
                    for i in assemblyids:
                        if i not in dict_proteinid_assemblyids[protein_record["Id"]]:
                            dict_proteinid_assemblyids[protein_record["Id"]].append(i)

    # NOTE:
    # some proteins, such as 487413233, occur within multiple sequences of the assembly!
    # -> assembly only listed once!

    proteinIds = list(dict_proteinid_assemblyids.keys())

    print("# proteins (unique): ", len(proteinIds))
    # -> # proteins with refseq source (<= # IPG proteins)

    # protein_id, assembly_ids
    out_file = open(prot_assembly_id_output,'w')
    for proteinid in proteinIds:
        print(proteinid, dict_proteinid_assemblyids[proteinid], file = out_file, flush=True)
    out_file.close()

    # 5) download protein FASTAs
    # print("    download proteins ...")
    # (or if mem problem: assembly-wise)
    success = False
    for attempt in range(3):
        try:
            with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=proteinIds) as entrez_handle:
                with open(fasta_output, "w") as out_handle:
                    out_handle.write(entrez_handle.read().replace('\n\n', '\n'))
                    # ... and get rid of additional blank lines between records (not sure why they are there)
            time.sleep(1)   # avoid getting blocked by ncbi
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
            sys.exit("Entrez download failed!")

    print("Done!")


if __name__ == "__main__":
    sys.exit(main())
