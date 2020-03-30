#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2019, Gerstein Lab"
__credits__ = ["Donghoon Lee","Mark Gerstein"]
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
parser.add_argument('--bed', help='BED File', required=True)
parser.add_argument('--bct', help='BCT File', required=True)
parser.add_argument('--bw', help='BigWig File of STARR-seq Coverage', required=True)

### optional args
parser.add_argument('--cov', help='Covariate File', required=False)
parser.add_argument('--threshold', help='Adjusted P-value Threshold', required=False, type=float, default=0.05)
parser.add_argument('--mode', help='Mode', required=False, type=int, default=1)
parser.add_argument('--mincov', help='Minimum Coverage', required=False, type=int, default=10)
parser.add_argument('--eq', help='Extreme Quantile to Remove', required=False, type=float, default=1e-5)

args = parser.parse_args()

if __name__ == "__main__": core.call_peak(prefix=args.prefix,
                                          chromSize=args.chromsize,
                                          bedFile=args.bed,
                                          covFile=args.cov,
                                          bctFile=args.bct,
                                          bwFile=args.bw,
                                          threshold=args.threshold,
                                          mode=args.mode,
                                          minCoverage=args.mincov,
                                          extQuantile=args.eq)
