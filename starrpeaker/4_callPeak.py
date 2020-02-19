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

parser = argparse.ArgumentParser(description='Call Peaks')

### required args
parser.add_argument('--prefix', help='Output File Prefix', required=True)
parser.add_argument('--chromsize', help='Chrom Sizes', required=True)
parser.add_argument('--bed', help='BED File', required=False)
parser.add_argument('--bct', help='BCT File', required=False)
parser.add_argument('--cov', help='Covariate File', required=False)
parser.add_argument('--bw', help='BigWig File of STARR-seq Coverage', required=False)

### optional args
parser.add_argument('--threshold', help='Adjusted P-value Threshold', required=False, default=0.05)
parser.add_argument('--mode', help='Mode', required=False, default=1)

args = parser.parse_args()

if args.bed is None:
    args.bed = args.prefix + ".bin.bed"

if args.bct is None:
    args.bct = args.prefix + ".bam.bct"

if args.cov is None:
    args.cov = args.prefix + ".cov.tsv"

if args.bw is None:
    args.bw = args.prefix + ".bam.bct.1.all.bw"

if __name__ == "__main__": core.call_peak(prefix=args.prefix,
                                          bedFile=args.bed,
                                          bctFile=args.bct,
                                          covFile=args.cov,
                                          bwFile=args.bw,
                                          chromSize=args.chromsize,
                                          threshold=args.threshold,
                                          mode=args.mode)
