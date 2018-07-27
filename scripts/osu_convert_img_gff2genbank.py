#!/usr/bin/env python

__author__ = "Jimmy Saw"

"""
Usage example:

osu_convert_img_gff2genbank.py \
    -g test.a.gff \
    -f test.a.fna \
    -p test.a.gene_product.txt \
    -l locus_id \
    -o test.gbk
"""

import re
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Data import CodonTable

def determine_codon_numbers(a, b):
    """
    Determines if the coding sequence is in multiples of 3
    :param a:
    :param b:
    :return:
    """
    length = (a - b) / 3.0
    if length.is_integer():
        return True
    else:
        return False

def determine_product_name(x, y, z):
    """
    Determine product names from gff file first. If it's not in gff file, check the product name from another file.
    :param x:
    :param y:
    :param z:
    :return: product name
    """
    product = "unknown"

    cols = x.split(";")
    for col in cols:
        if col.startswith("product="):
            product = re.sub("product=", "", col)

    if product == "unknown":
        df = pd.read_csv(y, sep="\t", header=None, names=['gene_id', 'product', 'ko'])
        for a, b in zip(df['gene_id'], df['product']):
            if a == z:
                product = b

    return product

def get_locus_tag(str):
    """
    Get locus_tag from gff file
    :param str:
    :return: locus_tag
    """
    locus_tag = "unknown"
    lst = str.split(";")
    for i in lst:
        if i.startswith("locus_tag="):
            locus_tag = re.sub("locus_tag=", "", i)
    return locus_tag

def create_genbank(gff, dna, locus, out, prod_names):
    """
    Parses gff file and DNA seq to add features to Genbank record
    :param gff:
    :param dna:
    :param locus:
    :param out:
    :return: genbank file
    """
    seqs = [i for i in SeqIO.parse(dna, "fasta")]
    gff_lines = [i.strip() for i in open(gff, "rU").readlines()]
    seq_dict = {}

    for srec in seqs:
        sequence_object = Seq(str(srec.seq), IUPAC.unambiguous_dna)
        record = SeqRecord(sequence_object, id=srec.id, name=locus, description=srec.description)
        seq_dict[srec.id] = record

    for line in gff_lines:
        col = line.split("\t")
        if col[0] in seq_dict:
            strand = 1
            if col[6].startswith("-"):
                strand = -1
            start = int(col[3])-1
            end = int(col[4])
            seq = Seq(str(seq_dict[col[0]].seq[start:end]))
            if col[2] == "CDS":
                if strand == 1:
                    trans = "Incomplete_codon_for_translation"
                    if determine_codon_numbers(start, end):
                        trans = Seq.translate(seq).strip("*")
                    f = SeqFeature(FeatureLocation(start, end), type="CDS", strand=strand)
                    f.qualifiers['translation'] = trans
                    locus_tag = get_locus_tag(col[-1])
                    f.qualifiers['locus_tag'] = locus_tag
                    product = determine_product_name(col[-1], prod_names, locus_tag)
                    f.qualifiers['product'] = product
                    seq_dict[col[0]].features.append(f)
                else:
                    trans = "Incomplete_codon_for_translation"
                    if determine_codon_numbers(start, end):
                        trans = Seq.translate(seq.reverse_complement()).strip("*")
                    f = SeqFeature(FeatureLocation(start, end), type="CDS", strand=strand)
                    f.qualifiers['translation'] = trans
                    locus_tag = get_locus_tag(col[-1])
                    f.qualifiers['locus_tag'] = locus_tag
                    product = determine_product_name(col[-1], prod_names, locus_tag)
                    f.qualifiers['product'] = product
                    seq_dict[col[0]].features.append(f)
            else:
                f = SeqFeature(FeatureLocation(start, end), type=col[2], strand=strand)
                locus_tag = get_locus_tag(col[-1])
                f.qualifiers['locus_tag'] = locus_tag
                product = determine_product_name(col[-1], prod_names, locus_tag)
                f.qualifiers['product'] = product
                seq_dict[col[0]].features.append(f)

    with open(out, "w") as outfile:
        for i, j in seq_dict.iteritems():
            SeqIO.write(j, outfile, 'genbank')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script converts gff to genbank.")
    parser.add_argument("-g", "--gff", required=True, help="GFF file.")
    parser.add_argument("-f", "--fasta", required=True, help="Sequence fasta file.")
    parser.add_argument("-p", "--products", help="Product names (if provided)")
    parser.add_argument("-l", "--locus", required=True, help="Locus id. Needs to be short. Eg: HIMB058")
    parser.add_argument("-o", "--out", required=True, help="Out file name")
    args = parser.parse_args()
    create_genbank(args.gff, args.fasta, args.locus, args.out, args.products)