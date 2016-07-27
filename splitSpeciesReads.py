"""
splitSpeciesReads.py

Designed to separate reads from the two species in hybrid data.
Adapted from GetGeneASEbyReads.py by Peter Combs.

Last edit: 07.26.2016

"""

from __future__ import print_function
from pysam import Samfile
from argparse import ArgumentParser, FileType
from collections import defaultdict, Counter
from os import path

def get_phase(read, snps):
    phase = None
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos + 1 in snps:
            if phase == None:
                try:
                    # 1 if alternate, -1 if reference
                    phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0 # This SNP isn't in the dataset
            else:
                try:
                    new_phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0
                if new_phase != phase:
                    return 0 # read seems misphased
    return phase

def get_snps(snpfile):
    snps = defaultdict(dict)
    if path.exists(path.join(path.dirname(snpfile), 'true_hets.tsv')):
        print("using true hets")
        true_hets = {tuple(line.strip().split()):True
                     for line in open(path.join(path.dirname(snpfile), 'true_hets.tsv'))
                    }
    else:
        true_hets = defaultdict(lambda x: True)
    if snpfile.endswith('.bed'):
        for line in open(snpfile):
            chrom, _, start, refalt = line.strip().split()
            if true_hets.get((chrom, start), True):
                snps[chrom][int(start)] = refalt.split('|')

    return snps

def split_reads(reads, ref_reads, alt_reads):
    ref_bam = Samfile(ref_reads, 'wb', template=reads)
    alt_bam = Samfile(alt_reads, 'wb', template=reads)

    prev_phase = None
    prev_read = None
    prev_qname = None

    for read in reads.fetch(until_eof=True):
        chrom = read.reference_name
        snps_on_chrom = snp_dict[chrom]
        phase = get_phase(read, snps_on_chrom)

        read_qname = read.qname
        if read_qname == prev_qname:
            if phase == prev_phase:
                if phase == 1:
                    ref_bam.write(read)
                    ref_bam.write(prev_read)
                elif phase == -1:
                    alt_bam.write(read)
                    alt_bam.write(prev_read)
        
        prev_read = read
        prev_phase = phase
        prev_qname = read_qname


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('snp_file')
    parser.add_argument('reads', type=Samfile)
    parser.add_argument('ref_reads')
    parser.add_argument('alt_reads')

    args = parser.parse_args()
    print(args)
    return args

if __name__ == "__main__":
    args = parse_args()
    snp_dict = get_snps(args.snp_file)
    split_reads(args.reads, args.ref_reads, args.alt_reads)
