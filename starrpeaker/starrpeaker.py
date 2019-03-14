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

parser = argparse.ArgumentParser(description='STARR-Peaker')

### required args
parser.add_argument('-p', '--prefix', help='Output File Prefix', required=True)
parser.add_argument('--chromsize', help='Chrom Sizes', required=True)
parser.add_argument('-l', '--length', help='Bin Length', required=False, default=500)
parser.add_argument('-s', '--step', help='Step Size', required=False, default=100)
parser.add_argument('-x', '--blacklist', help='Blacklist Region in BED format', required=True)
parser.add_argument('--cov', help='Covariate BigWig Files', nargs='+', required=True)
parser.add_argument('-i', '--input', help='Input BAM File', nargs='+', required=True)
parser.add_argument('-o', '--output', help='STARR-seq BAM File', nargs='+', required=True)

### optional args
parser.add_argument('--min', help='Minimum Template Size', required=False, default=100)
parser.add_argument('--max', help='Maximum Template Size', required=False, default=1000)
parser.add_argument('-t', '--threshold', help='Adjusted P-value Threshold', required=False, default=0.05)

args = parser.parse_args()


def main():
    ### make bin with specified window length and step size
    core.make_bin(args.chromsize,
                  args.length,
                  args.step,
                  args.blacklist,
                  args.prefix + ".bin.bed")

    ### process covariates
    core.proc_cov(args.cov,
                  args.prefix + ".bin.bed",
                  args.prefix + ".cov.tsv")

    ### process output bam
    core.proc_bam(args.output,
                  args.prefix + ".bin.bed",
                  args.chromsize,
                  args.prefix + ".output.tsv",
                  args.min,
                  args.max)

    ### process input bam
    core.proc_bam(args.input,
                  args.prefix + ".bin.bed",
                  args.chromsize,
                  args.prefix + ".input.tsv",
                  args.min,
                  args.max)

    totalInput = core.count_total_proper_templates(args.input[0], args.min, args.max)
    totalOutput = core.count_total_proper_templates(args.output[0], args.min, args.max)

    ### call peaks
    core.call_peak(args.prefix + ".input.tsv",
                   args.prefix + ".output.tsv",
                   args.prefix + ".cov.tsv",
                   args.prefix + ".bin.bed",
                   args.prefix + ".peak.bed",
                   args.threshold,
                   totalInput,
                   totalOutput)


if __name__ == "__main__": main()
