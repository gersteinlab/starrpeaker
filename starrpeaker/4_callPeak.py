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
parser.add_argument('-p', '--prefix', help='Output File Prefix', required=True)
parser.add_argument('--bed', help='Bin BED File', required=True)
parser.add_argument('--bct', help='Bincount BCT File', required=True)
parser.add_argument('--cov', help='Covariate File', required=True)
parser.add_argument('--bw', help='STARR-seq Output Coverage BigWig File', required=True)
parser.add_argument('-c', '--chromsize', help='Chrom Sizes', required=True)

### optional args
parser.add_argument('-t', '--threshold', help='Adjusted P-value Threshold', required=False, default=0.05)
parser.add_argument('-q', '--minquantile', help='Minimum Input Quantile', required=False, default=0.2)
parser.add_argument('-m', '--mode', help='Mode', required=False, default=1)

args = parser.parse_args()

if __name__ == "__main__": core.call_peak(prefix=args.prefix,
                                          bedFile=args.bed,
                                          bctFile=args.bct,
                                          covFile=args.cov,
                                          bwFile=args.bw,
                                          chromSize=args.chromsize,
                                          threshold=args.threshold,
                                          minInputQuantile=args.minquantile,
                                          mode=args.mode)
