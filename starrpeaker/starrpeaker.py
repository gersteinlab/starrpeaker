#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2019, Gerstein Lab"
__credits__ = ["Donghoon Lee", "Mark Gerstein"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

import argparse
import core

### LOAD ARGs ###

parser = argparse.ArgumentParser(description='STARRPeaker')

### required args
parser.add_argument('--prefix', help='Output File Prefix', required=True)
parser.add_argument('--chromsize', help='Chrom Sizes', required=True)
parser.add_argument('--blacklist', help='Blacklist Region in BED format', required=True)
parser.add_argument('-i', '--input', help='Input BAM File', required=True)
parser.add_argument('-o', '--output', help='STARR-seq BAM File', required=True)

### optional args
parser.add_argument('--length', help='Bin Length', required=False, type=int, default=500)
parser.add_argument('--step', help='Step Size', required=False, type=int, default=100)
parser.add_argument('--cov', help='Covariate BigWig Files', nargs='+', required=False)
parser.add_argument('--min', help='Minimum Template Size', required=False, type=int, default=200)
parser.add_argument('--max', help='Maximum Template Size', required=False, type=int, default=1000)
parser.add_argument('--readstart', help='Use Read Start Position instead of Fragment Center', required=False, action='store_true')
parser.add_argument('--strand', help='Use all/fwd/rev Stranded Fragments', required=False, type=str, default="all")
parser.add_argument('--threshold', help='Adjusted P-value Threshold', required=False, type=float, default=0.05)
parser.add_argument('--mode', help='Mode', required=False, type=int, default=1)
parser.add_argument('--mincov', help='Minimum Coverage', required=False, type=int, default=10)
parser.add_argument('--eq', help='Extreme Quantile to Remove', required=False, type=float, default=1e-5)
parser.add_argument('--se', help='Use Single-End instead of Paired-end Sequencing', required=False, action='store_true')
parser.add_argument('--minfc', help='Minumum Fold Change', required=False, type=float, default=1.5)

args = parser.parse_args()


def main():
    ### make genomic bin with specified window length and step size
    core.make_bin(prefix=args.prefix,
                  chromSize=args.chromsize,
                  binLength=args.length,
                  stepSize=args.step,
                  blackList=args.blacklist)

    ### process covariates
    if args.cov:
        core.proc_cov(prefix=args.prefix,
                      bedFile=args.prefix + ".bin.bed",
                      bwFiles=args.cov)

    ### process input, output bam files
    core.proc_bam(prefix=args.prefix,
                  chromSize=args.chromsize,
                  bedFile=args.prefix + ".bin.bed",
                  bamFiles=[args.input, args.output],
                  minSize=args.min,
                  maxSize=args.max,
                  readStart=args.readstart,
                  strand=args.strand,
                  singleEnd=args.se)

    ### call peaks
    if args.cov:
        core.call_peak(prefix=args.prefix,
                       chromSize=args.chromsize,
                       bedFile=args.prefix + ".bin.bed",
                       covFile=args.prefix + ".cov.tsv",
                       bctFile=args.prefix + ".bam.bct",
                       bwFile=args.prefix + ".output.bw",
                       threshold=args.threshold,
                       mode=args.mode,
                       minCoverage=args.mincov,
                       extQuantile=args.eq,
                       minFC=args.minfc)
    else:
        core.call_peak(prefix=args.prefix,
                       chromSize=args.chromsize,
                       bedFile=args.prefix + ".bin.bed",
                       covFile=None,
                       bctFile=args.prefix + ".bam.bct",
                       bwFile=args.prefix + ".output.bw",
                       threshold=args.threshold,
                       mode=args.mode,
                       minCoverage=args.mincov,
                       extQuantile=args.eq,
                       minFC=args.minfc)


if __name__ == "__main__": main()
