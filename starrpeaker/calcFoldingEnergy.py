#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2019, Gerstein Lab"
__credits__ = ["Donghoon Lee", "Mark Gerstein"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

import core
import argparse

### LOAD ARGs ###

parser = argparse.ArgumentParser(description='Calculate Folding Energy')

### required args
parser.add_argument('--bed', help='BED File', required=True)
parser.add_argument('--tsv', help='STARR-seq BAM File', required=True)
parser.add_argument('--linearfold', help='Path to LinearFold binary', required=True)
parser.add_argument('--genome', help='Genome FASTA File', required=True)

args = parser.parse_args()

def main():
    core.proc_fenergy(bedFile=args.bed,
                      fileOut=args.out,
                      linearfold=args.linearfold,
                      genome=args.fa)

if __name__ == "__main__": main()
