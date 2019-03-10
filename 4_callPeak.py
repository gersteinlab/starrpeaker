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
parser.add_argument('-i', '--input', help='Input File', required=True)
parser.add_argument('-o', '--output', help='Output File', required=True)
parser.add_argument('-c', '--cov', help='Covariate File', required=True)
parser.add_argument('-b', '--bed', help='Bin BED File', required=True)
parser.add_argument('-p', '--prefix', help='Output File Prefix', required=True)

### optional args
parser.add_argument('-t', '--threshold', help='Adjusted P-value Threshold', required=False, default=0.05)
parser.add_argument('--totalinput', help='Total Input Count', required=False, default=0)
parser.add_argument('--totaloutput', help='Total Output Count', required=False, default=0)

args = parser.parse_args()

if __name__ == "__main__": core.call_peak(inputFile=args.input,
                                          outputFile=args.output,
                                          covFile=args.cov,
                                          bedFile=args.bed,
                                          fileOut=args.prefix + ".peak.bed",
                                          threshold=args.threshold,
                                          totalInput=args.totalinput,
                                          totalOutput=args.totaloutput)
