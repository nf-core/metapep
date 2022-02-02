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
#import tqdm
import io
import argparse
import time

from Bio import Entrez, SeqIO
from pprint import pprint
from datetime import datetime
from collections import Counter, defaultdict
from urllib.error import HTTPError
import pandas as pd

import sys


# TODO
# clean code
# double check!


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', "--taxon", required=True, metavar='TAXON ID', help="Taxon ID to download proteins for.")
    parser.add_argument('-e', "--email", required=True, help="Email address to use for NCBI access.")
    parser.add_argument('-k', "--key", required=True, help="NCBI key to allow faster access.")
    parser.add_argument('-p', "--proteins", required=True, metavar='FILE', help="Output file (compressed) containing: protein_tmp_id, protein_sequence.")
    parser.add_argument('-f', "--fasta", required=True, metavar='FILE', help="Output fasta file (compressed).")
    parser.add_argument('-ta', "--taxa_assemblies", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: taxon_id, assembly_id.")
    parser.add_argument('-ep', "--entities_proteins", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing: protein_tmp_id, entity_name (taxon_id).")
    return parser.parse_args(args)


# get assembly length ("total_length") from entrez
def get_assembly_length(assemblyId):
    success = False
    for attempt in range(3):
        try:
            with Entrez.efetch(db="assembly", rettype="docsum", id=assemblyId) as entrez_handle:
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
        sys.exit("Entrez efetch download failed!")

    root = ET.fromstring("<root>" + str(assembly_stats['DocumentSummarySet']['DocumentSummary'][0]['Meta']) + "</root>")
    return root.find("./Stats/Stat[@category='total_length'][@sequence_tag='all']").text

def get_elink(dbfrom, db, LinkName, id):
    success = False
    for attempt in range(3):
        try:
            with Entrez.elink(dbfrom=dbfrom, db=db, LinkName=LinkName, id=id) as entrez_handle:
                results = Entrez.read(entrez_handle)
            time.sleep(1)   # avoid getting blocked by ncbi
            success = True
            return results
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(10)
            else:
                raise
    if not success:
        sys.exit("Entrez elink download failed!")

def main(args=None):
    args = parse_args(args)

    # setup entrez email
    Entrez.api_key = args.key
    Entrez.email = args.email

    taxon = args.taxon

    ####################################################################################################
    # # 1) for each taxId -> get all assembly IDs
    # print("# taxa: ", len(taxIds))
    print("retrieve assembly IDs ...")

    assembly_result = get_elink(dbfrom="taxonomy", db="assembly", LinkName="taxonomy_assembly", id=taxon)[0]

    # # 2) select one assembly (largest for now)
    print("get assembly lengths and select largest assembly ...")

    if len(assembly_result["LinkSetDb"]) > 0:
        # get all assembly ids
        ids = [ assembly_record["Id"] for assembly_record in assembly_result["LinkSetDb"][0]["Link"] ]
        # get corresponding lengths
        lengths = [get_assembly_length(id) for id in ids]
        # get id for largest assembly
        selected_assemblyId = ids[lengths.index(max(lengths))]

    # write taxId - assemblyId out
    print("taxon_id", "assembly_id", sep='\t', file=args.taxa_assemblies, flush=True)
    print(taxon, selected_assemblyId, sep='\t', file=args.taxa_assemblies, flush=True)

    # 3) (selected) assembly -> nucleotide sequences
    # (maybe split here)
    print("get nucloetide sequence IDs for assembly...")

    nucleotide_result = get_elink(dbfrom="assembly", db="nuccore", LinkName="assembly_nuccore_refseq", id=selected_assemblyId)[0]

    seqIds = set([record["Id"] for record in nucleotide_result["LinkSetDb"][0]["Link"]])

    print("# nucleotide sequences (unique): ", len(seqIds))
    # -> # contigs

    # 4) nucelotide sequences -> proteins
    print("for each nucleotide sequence get proteins ...")

    protein_results = get_elink(dbfrom="nuccore", db="protein", LinkName="nuccore_protein", id=seqIds)

    # NOTE:
    # some proteins, such as 487413233, occur within multiple sequences of the assembly!
    # -> only listed once!

    proteinIds = set([record["Id"] for protein_result in protein_results if len(protein_result["LinkSetDb"]) > 0 for record in protein_result["LinkSetDb"][0]["Link"] ])

    print("# proteins (unique): ", len(proteinIds))
    # -> # proteins with refseq source (<= # IPG proteins)

    # 5) download protein FASTAs, convert to TSV
    print("    download proteins ...")
    # (or if mem problem: assembly-wise)
    # TODO check if max. number of item that can be returned by efetch (retmax)!? compare numbers!

    # download actual fasta records and write out tsv, fasta and entities_proteins.entrez.tsv
    success = False
    for attempt in range(3):
        try:
            with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=proteinIds) as entrez_handle:
                records = list(SeqIO.parse(entrez_handle, "fasta"))
                with gzip.open(args.fasta, 'wt') as out_handle:
                    SeqIO.write(records, out_handle, 'fasta')
                with gzip.open(args.proteins, 'wt') as out_handle:
                    print("protein_tmp_id", "protein_sequence", sep='\t', file=out_handle)
                    print("protein_tmp_id", "entity_name", sep='\t', file=args.entities_proteins)
                    for record in records:
                        print(record.id, record.seq, sep='\t', file=out_handle, flush=True)
                        print(record.id, taxon, sep='\t', file=args.entities_proteins, flush=True)
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
            sys.exit("Entrez efetch download failed!")

    # NOTE for proteins the NCBI accession version ids are written out to enable a mapping to the fasta/tsv output
    # in contrast, for assemblies UIDs are written out

    print("Done!")


if __name__ == "__main__":
    sys.exit(main())
