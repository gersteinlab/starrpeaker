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
parser.add_argument('-p', '--prefix', help='Output File Prefix', required=True)
parser.add_argument('--chromsize', help='Chrom Sizes', required=True)
parser.add_argument('-l', '--length', help='Bin Length', required=False, default=500)
parser.add_argument('-s', '--step', help='Step Size', required=False, default=100)
parser.add_argument('-x', '--blacklist', help='Blacklist Region in BED format', required=True)
parser.add_argument('--cov', help='Covariate BigWig Files', nargs='+', required=True)
parser.add_argument('-i', '--input', help='Input BAM File', nargs='+', required=True)
parser.add_argument('-o', '--output', help='STARR-seq BAM File', nargs='+', required=True)

### optional args
parser.add_argument('--min', help='Minimum Template Size', required=False, default=100)
parser.add_argument('--max', help='Maximum Template Size', required=False, default=1000)
parser.add_argument('-t', '--threshold', help='Adjusted P-value Threshold', required=False, default=0.05)

args = parser.parse_args()


def main():
    ### make bin with specified window length and step size
    core.make_bin(chromSize=args.chromsize,
                  binLength=args.length,
                  stepSize=args.step,
                  blackList=args.blacklist,
                  fileOut=args.prefix + ".bin.bed")

    ### process covariates
    core.proc_cov(bwFiles=args.cov,
                  bedFile=args.prefix + ".bin.bed",
                  fileOut=args.prefix + ".cov.tsv")

    ### process input, output bam
    core.proc_bam(bamFiles=[args.input, args.output],
                  bedFile=args.prefix + ".bin.bed",
                  chromSize=args.chromsize,
                  fileOut=args.prefix + ".bam.bct",
                  minSize=args.min,
                  maxSize=args.max,
                  normalize=True)

    ### call peaks
    core.call_peak(bctFile=args.prefix + ".bam.bct",
                   covFile=args.prefix + ".cov.tsv",
                   bedFile=args.prefix + ".bin.bed",
                   fileOut=args.prefix + ".peak.bed",
                   threshold=args.threshold)

    ### make signal tracks
    core.make_bigwig(chromsize=args.chromsize,
                     bedFile=args.prefix + ".bin.bed",
                     bctFile=args.prefix + ".bam.bct",
                     prefix=args.prefix)


if __name__ == "__main__": main()
