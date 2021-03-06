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

parser = argparse.ArgumentParser(description='Process BigWig File(s)')

### required args

parser.add_argument('--prefix', help='Output File Prefix', required=True)
parser.add_argument('--cov', help='BigWig Files as Covariates', nargs='+', required=True)

args = parser.parse_args()

if __name__ == "__main__": core.proc_cov(prefix=args.prefix,
                                         bedFile=args.prefix + ".bin.bed",
                                         bwFiles=args.cov)
