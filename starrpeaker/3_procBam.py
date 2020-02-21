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

parser = argparse.ArgumentParser(description='Process BAM File(s)')

### required args
parser.add_argument('--prefix', help='Output File Prefix', required=True)
parser.add_argument('--chromsize', help='Chrom Sizes', required=True)
parser.add_argument('-i', '--input', help='STARR-seq Input BAM File', required=True)
parser.add_argument('-o', '--output', help='STARR-seq Output BAM File', required=True)

### optional args
parser.add_argument('--min', help='Minimum Insert Size', required=False, type=int, default=200)
parser.add_argument('--max', help='Maximum Insert Size', required=False, type=int, default=1000)

args = parser.parse_args()

if __name__ == "__main__": core.proc_bam(bamFiles=[args.input, args.output],
                                         bedFile=args.prefix + ".bin.bed",
                                         chromSize=args.chromsize,
                                         fileOut=args.prefix + ".bam.bct",
                                         minSize=args.min,
                                         maxSize=args.max)
