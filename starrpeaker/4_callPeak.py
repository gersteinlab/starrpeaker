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

### optional args
parser.add_argument('--threshold', help='Adjusted P-value Threshold', required=False, default=0.05)
parser.add_argument('--minquantile', help='Minimum Input Quantile', required=False, default=0.2)
parser.add_argument('--mode', help='Mode', required=False, default=1)

args = parser.parse_args()

if __name__ == "__main__": core.call_peak(prefix=args.prefix,
                                          bedFile=args.prefix + ".bin.bed",
                                          bctFile=args.prefix + ".bam.bct",
                                          covFile=args.prefix + ".cov.tsv",
                                          bwFile=args.prefix + ".bam.bct.1.bw",
                                          chromSize=args.chromsize,
                                          threshold=args.threshold,
                                          minInputQuantile=args.minquantile,
                                          mode=args.mode)
