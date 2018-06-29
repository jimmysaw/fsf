#!/usr/bin/env python
from __future__ import print_function

__author__ = "Jimmy Saw"
"""
This script can automatically download genomes from NCBI Genbank. You will need to download "assembly_summary_genbank.txt" file
from here: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
Downloading everything in the table would take many hours, so be careful with that.
You can create a subset of genomes of interest using another script "osu_query_ncbi_genomes_taxids.py" to select just
those you are interested.

Usage examples:

## download all genomes available from NCBI
osu_download_ncbi_genomes.py -a assembly_summary_genbank.txt -d test3 -t all -e you@email.com -s download_summary.txt

## download all genomes within a subset table (selected to check for desired taxonomic affiliation; see osu_query_ncbi_genomes_taxids.py script)
osu_download_ncbi_genomes.py -a subset.txt -d test3 -t all -e you@mail.com -s download_summary.txt

## download first 5 Gammaproteobacteria from a subset table
osu_download_ncbi_genomes.py -a subset.txt -d test3 -t Gammaproteobacteria -r class -n 5 -e you@email.com -s download_summary.txt

## download all Gammaproteobacteria from a subset table
osu_download_ncbi_genomes.py -a subset.txt -d test3 -t Gammaproteobacteria -r class -n all -e you@email.com -s download_summary.txt

## download all Eukaryota from a subset table
osu_download_ncbi_genomes.py -a subset.txt -d test3 -t Eukaryota -r superkingdom -n all -e you@email.com -s download_summary.txt

## download 5 Crenarchaeota from a subset table
osu_download_ncbi_genomes.py -a subset.txt -d test3 -t Crenarchaeota -r phylum -n 5 -e you@email.com -s download_summary.txt
"""

import os
import argparse
import pandas as pd
from Bio import Entrez
import urllib

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

def query_rank(tid, em, rank, tax):
    """
    Queries NCBI taxonomy database for a given lineage info that match required taxnomic rank
    :param tid:
    :param em:
    :param rank:
    :param tax:
    :return: scientific name of the organism
    """
    Entrez.email = em
    h = Entrez.efetch(db="taxonomy", id=tid)
    rec = Entrez.read(h)
    for i in rec[0]['LineageEx']:
        if i['Rank'] == rank:
            if i['ScientificName'] == tax:
                print(i['ScientificName'])
                return i['ScientificName']

def download_genome(k, v, dd):
    """
    Downloads genomes into the specified folder
    (Currently only downloads gbff and faa files. Can be modified to download other associated files)
    (Eg: just add: [if f.endswith("_genomic.fna.gz"):] if an .fna file is required)
    :param k:
    :param v:
    :param dd:
    :return: genomes (genbank file and proteome file if available) into a specific folder
    """
    ftp_path = v[18]
    md5_path = '/'.join([ftp_path, "md5checksums.txt"])
    if not os.path.exists(dd):
        os.makedirs(dd)
    dir = os.path.join(dd, k)
    if not os.path.exists(dir):
        os.makedirs(dir)
    outfile = os.path.join(dir, "md5checksums.txt")
    if not os.path.isfile(outfile):
        urllib.urlretrieve(md5_path, outfile)
        urllib.urlcleanup()
    with open(outfile) as fh:
        lines = fh.readlines()
        for line in lines:
            c = line.split()
            f = c[1].strip().split("/")[1]
            if f.endswith("protein.faa.gz"):
                prot_out = os.path.join(dir, f)
                prot_path = '/'.join([ftp_path, f])
                if not os.path.isfile(prot_out):
                    urllib.urlretrieve(prot_path, prot_out)
                    urllib.urlcleanup()
            if f.endswith("_genomic.gbff.gz"):
                genome_out = os.path.join(dir, f)
                genome_path = '/'.join([ftp_path, f])
                if not os.path.isfile(genome_out):
                    urllib.urlretrieve(genome_path, genome_out)
                    urllib.urlcleanup()

def sample_by_rank(x, download_dir, email, tax, num, rank, summary):
    """
    Sample by rank
    :param x:
    :param download_dir:
    :param email:
    :param tax:
    :param num:
    :param rank:
    :return:
    """
    df = pd.read_csv(x, sep="\t", header=1)
    tdict = df.set_index('# assembly_accession').T.to_dict('list')
    downloaded = '{0}\t{1}\t{2}\t{3}\t{4}\n'.format("accession", "taxid", "rank", "taxonomy", "name")

    if not num == "all":
        number = 0
        target = int(num)
        print('Attempting to download {0} {1} (Rank: {2}) genomes from {3}\n'.format(target, tax, rank, x))
        for k, v in tdict.iteritems():
            if number < target:
                taxid = v[4]
                t_rank = query_rank(taxid, email, rank, tax)
                if t_rank == tax:
                    number += 1
                    print("Found", number, "of", target, "target genomes to download.")
                    lineage = query_taxid(taxid, email)
                    scientific_name = lineage[1]
                    print("Downloading", k, scientific_name, "to", download_dir, "\n")
                    download_genome(k, v, download_dir)
                    downloaded += '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(k, taxid, rank, tax, scientific_name)
        print("Downloaded first", number, "genomes found out of", target, "genomes required.")
    else:
        print('Downloading all genomes within {0} (Rank: {1}) from {2}\n'.format(tax, rank, x))
        for k, v in tdict.iteritems():
            taxid = v[4]
            t_rank = query_rank(taxid, email, rank, tax)
            if t_rank == tax:
                lineage = query_taxid(taxid, email)
                scientific_name = lineage[1]
                print("Downloading", k, scientific_name, "to", download_dir, "\n")
                download_genome(k, v, download_dir)
                downloaded += '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(k, taxid, rank, tax, scientific_name)
    with open(summary, "w") as sum_out:
        sum_out.write(downloaded)

def download_all_genomes(x, download_dir, email, summary):
    """
    Downloads all the genomes from NCBI (use with care!)
    :param x:
    :param download_dir:
    :param email:
    :return:
    """
    df = pd.read_csv(x, sep="\t", header=1)
    tdict = df.set_index('# assembly_accession').T.to_dict('list')
    downloaded = '{0}\t{1}\t{2}\n'.format("accession", "taxid", "name")

    print("Downloading all genomes within", x, "\n")
    for k, v in tdict.iteritems():
        taxid = v[4]
        lineage = query_taxid(taxid, email)
        scientific_name = lineage[1]
        print("Downloading", k, scientific_name, "to", download_dir, "\n")
        download_genome(k, v, download_dir)
        downloaded += '{0}\t{1}\t{2}\n'.format(k, taxid, scientific_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script downloads genomes from ncbi ftp site")
    parser.add_argument("-a", "--assembly", required=True, help="Assembly stats file")
    parser.add_argument("-d", "--directory", required=True, help="Folder name to save files to")
    parser.add_argument("-t", "--tax", required=True, help="Taxonomic affiliation within a rank to download (choose: all or specific phylum, class, etc)")
    parser.add_argument("-r", "--rank", help="Rank to check (eg:  superkingdom, kingdom, phylum, class, order")
    parser.add_argument("-n", "--number", help="Number to download")
    parser.add_argument("-e", "--email", required=True, help="Your email address")
    parser.add_argument("-s", "--summary", required=True, help="Download summary")
    args = parser.parse_args()
    if args.tax == 'all':
        download_all_genomes(args.assembly, args.directory, args.email, args.summary)
    else:
        sample_by_rank(args.assembly, args.directory, args.email, args.tax, args.number, args.rank, args.summary)


