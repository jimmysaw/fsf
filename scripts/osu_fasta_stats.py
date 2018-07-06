#!/usr/bin/env python

__author__ = "Jimmy Saw"

"""
Usage example:
ous_fasta_stats.py -f multi.fasta

"""

import argparse
import locale
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC

def calculate_N50(ctgsizes):
    """
    Calculates N50 statistics
    :param ctgsizes:
    :return: length of N50
    """
    totalsum = 0
    slist = []

    for s in ctgsizes:
        slist.append(s)
        totalsum += s

    halfsum = totalsum / 2.0
    sortedlist = sorted(slist, reverse=True)
    cumulative = 0
    tmp = []
    i = 0

    while i < len(sortedlist):
        cumulative += sortedlist[i]
        if cumulative >= halfsum:
            tmp.append(i)
        i += 1

    n50len = sortedlist[tmp[0]]

    return n50len

def runstats(fasta_file):
    """
    Parses and multi-fasta file and calculates stats
    :param fasta_file:
    :return: prints stats
    """
    loc = locale.setlocale(locale.LC_ALL, '')
    locale.setlocale(locale.LC_ALL, loc)

    seqs = [seq for seq in SeqIO.parse(fasta_file, "fasta")]
    lengths_list = [len(i.seq) for i in seqs]

    total_size = locale.format("%d", np.sum(lengths_list), grouping=True)
    ctg_n50 = locale.format("%d", calculate_N50(lengths_list), grouping=True)
    minsize = locale.format("%d", min(lengths_list), grouping=True)
    maxsize = locale.format("%d", max(lengths_list), grouping=True)
    gc_percent = np.average([GC(i.seq) for i in seqs])

    print '{0:25}\t{1}'.format("Total number of contigs:", locale.format("%d", len(seqs), grouping=True))
    print '{0:25}\t{1}'.format("Total size of all contigs:", total_size)
    print '{0:25}\t{1}'.format("N50 of all contigs:", ctg_n50)
    print '{0:25}\t{1}'.format("Largest contig:", maxsize)
    print '{0:25}\t{1}'.format("Smallest contig:", minsize)
    print '{0:25}\t{1:.2f}'.format("G+C % of contigs:", gc_percent)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script calculates basic statistics of a given fasta/multi-fasta"
                                                 "file")
    parser.add_argument("-f", "--fasta", required=True, help="Fasta file to check stats")
    args = parser.parse_args()
    runstats(args.fasta)
