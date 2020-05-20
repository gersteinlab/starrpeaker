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

parser = argparse.ArgumentParser(description='Process BAM File(s)')

### required args
parser.add_argument('--prefix', help='Output File Prefix', required=True)
parser.add_argument('--chromsize', help='Chrom Sizes', required=True)
parser.add_argument('--bed', help='BED File', required=True)
parser.add_argument('-i', '--input', help='STARR-seq DNA Input BAM File', required=True)
parser.add_argument('-o', '--output', help='STARR-seq RNA Output BAM File', required=True)

### optional args
parser.add_argument('--min', help='Minimum Insert Size', required=False, type=int, default=200)
parser.add_argument('--max', help='Maximum Insert Size', required=False, type=int, default=1000)
parser.add_argument('--readstart', help='Use Read Start Position instead of Fragment Center', required=False, action='store_true')
parser.add_argument('--strand', help='Use all/fwd/rev Stranded Fragments', required=False, type=str, default="all")
parser.add_argument('--se', help='Use Single-End instead of Paired-end Sequencing', required=False, action='store_true')

args = parser.parse_args()

if __name__ == "__main__": core.proc_bam(prefix=args.prefix,
                                         chromSize=args.chromsize,
                                         bedFile=args.bed,
                                         bamFiles=[args.input, args.output],
                                         minSize=args.min,
                                         maxSize=args.max,
                                         readStart=args.readstart,
                                         strand=args.strand,
                                         singleEnd=args.se)
