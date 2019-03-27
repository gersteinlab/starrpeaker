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

parser = argparse.ArgumentParser(description='Make BigWigs')

### required args
parser.add_argument('-p', '--prefix', help='Output BigWig Prefix', required=True)
parser.add_argument('--bed', help='Bin BED File', required=True)
parser.add_argument('--bct', help='Count BCT File', required=True)
parser.add_argument('-c', '--chromsize', help='Chrom Sizes', required=True)

### optional args
parser.add_argument('--pval', help='P-value BedGraph File', required=False, default="")

args = parser.parse_args()

if __name__ == "__main__": core.make_bigwig(prefix=args.prefix,
                                            bedFile=args.bed,
                                            bctFile=args.bct,
                                            chromsize=args.chromsize,
                                            bedGraphFile=args.pval)
