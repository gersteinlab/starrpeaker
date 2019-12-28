#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2019, Gerstein Lab"
__credits__ = ["Donghoon Lee"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

import argparse
import core

### LOAD ARGs ###

parser = argparse.ArgumentParser(description='STARR-Peaker')

### required args
parser.add_argument('--prefix', help='Output File Prefix', required=True)
parser.add_argument('--chromsize', help='Chrom Sizes', required=True)
parser.add_argument('--blacklist', help='Blacklist Region in BED format', required=True)
parser.add_argument('--cov', help='Covariate BigWig Files', nargs='+', required=True)
parser.add_argument('-i', '--input', help='Input BAM File', required=True)
parser.add_argument('-o', '--output', help='STARR-seq BAM File', required=True)

### optional args
parser.add_argument('--length', help='Bin Length', required=False, default=500)
parser.add_argument('--step', help='Step Size', required=False, default=100)
parser.add_argument('--min', help='Minimum Template Size', required=False, default=100)
parser.add_argument('--max', help='Maximum Template Size', required=False, default=1000)
parser.add_argument('--threshold', help='Adjusted P-value Threshold', required=False, default=0.05)
parser.add_argument('--minquantile', help='Minimum Input Quantile', required=False, default=0.2)
parser.add_argument('--mode', help='Mode', required=False, default=1)

args = parser.parse_args()


def main():
    ### make genomic bin with specified window length and step size
    core.make_bin(chromSize=args.chromsize,
                  binLength=args.length,
                  stepSize=args.step,
                  blackList=args.blacklist,
                  fileOut=args.prefix + ".bin.bed")

    ### process covariates
    core.proc_cov(bwFiles=args.cov,
                  bedFile=args.prefix + ".bin.bed",
                  fileOut=args.prefix + ".cov.tsv")

    ### process input, output bam files
    core.proc_bam(bamFiles=[args.input, args.output],
                  bedFile=args.prefix + ".bin.bed",
                  chromSize=args.chromsize,
                  fileOut=args.prefix + ".bam.bct",
                  minSize=args.min,
                  maxSize=args.max)

    ### call peaks
    core.call_peak(prefix=args.prefix,
                   bedFile=args.prefix + ".bin.bed",
                   bctFile=args.prefix + ".bam.bct",
                   covFile=args.prefix + ".cov.tsv",
                   bwFile=args.prefix + ".bam.bct.1.bw",
                   chromSize=args.chromsize,
                   threshold=args.threshold,
                   mode=args.mode)

if __name__ == "__main__": main()
