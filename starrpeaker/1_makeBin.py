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

parser = argparse.ArgumentParser(description='Make Bin in BED Format')

### required args
parser.add_argument('--prefix', help='Output File Prefix', required=True)
parser.add_argument('--chromsize', help='Chrom Sizes', required=True)
parser.add_argument('--blacklist', help='Blacklist Region in BED format', required=True)

### optional args
parser.add_argument('--length', help='Bin Length', required=False, type=int, default=500)
parser.add_argument('--step', help='Step Size', required=False, type=int, default=100)

args = parser.parse_args()

if __name__ == "__main__": core.make_bin(prefix=args.prefix,
                                         chromSize=args.chromsize,
                                         binLength=args.length,
                                         stepSize=args.step,
                                         blackList=args.blacklist)
