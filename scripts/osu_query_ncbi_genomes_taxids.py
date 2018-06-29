#!/usr/bin/env python
from __future__ import print_function

__author__ = "Jimmy Saw"
"""
This script is designed to be run just once for each update of assembly_summary_genbank.txt file you can get from:
ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
The idea is to help you to get lineages first, so you can pick which organism(s) to download for your work.
Another script can then take the subset to download so that it wouldn't take too long to sift through the original table.

Usage example:
osu_query_ncbi_genomes_taxids.py -a assembly_summary_genbank.txt -e you@email.com -o test.out (this might take several hours to finish)
head -2 assembly_summary_genbank.txt > subset.txt
for i in `grep Cyanobacteria test.out | awk '{print $1}'`;do
    grep $i assembly_summary_genbank.txt >> subset.txt
done

Then run:
osu_download_ncbi_genomes.py -a subset.txt -d out_folder -e you@email.com -p Chloroflexi -n [10 | all]
"""

import argparse
import pandas as pd
from Bio import Entrez
import progressbar
from progressbar import ProgressBar

def query_taxid(tid, em):
    """
    Sends Entrez queries for tax ids and return lineage info
    :param tid:
    :param em:
    :return: tuple of (lineage, scientific name)
    """
    Entrez.email = em
    h = Entrez.efetch(db="taxonomy", id=tid)
    rec = Entrez.read(h)
    lineage = rec[0]['Lineage']
    sci_name = rec[0]['ScientificName']
    return (lineage, sci_name)

def parse_table(table, email, outfile):
    """
    Parses the table of genomes to get taxid to lineage info
    :param table:
    :param email:
    :param outfile:
    :return:
    """
    df = pd.read_csv(table, sep="\t", header=1)
    tuples = []
    for a, b in zip(df['# assembly_accession'], df['taxid']):
        tuples.append((a, b))
    pbar = ProgressBar(widgets=['[', progressbar.Timer(), ']', progressbar.Bar(), ' (', progressbar.ETA(), ') '])
    string = ""
    for i, j in enumerate(pbar(tuples)):
        taxid = j[1]
        lineage = query_taxid(taxid, email)
        string += '{0}\t{1}\n'.format(j[0], lineage[0])
    with open(outfile, "w") as outf:
        outf.write(string)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script queries NCBI genomes taxids to print a list of taxonomic lineages")
    parser.add_argument("-a", "--assembly", required=True, help="Assembly stats file")
    parser.add_argument("-e", "--email", required=True, help="Your email address")
    parser.add_argument("-o", "--outfile", required=True, help="Outfile to save to")
    args = parser.parse_args()
    parse_table(args.assembly, args.email, args.outfile)
