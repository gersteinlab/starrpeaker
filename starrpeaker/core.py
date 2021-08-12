#!/usr/bin/python
from __future__ import division

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2019, Gerstein Lab"
__credits__ = ["Donghoon Lee", "Mark Gerstein"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

import numpy as np
import pandas as pd
import pybedtools
import pyBigWig
import pysam
from scipy.stats import nbinom
from scipy.special import digamma, polygamma
import statsmodels.formula.api as smf
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
import os, uuid, datetime
from itertools import compress
from subprocess import check_output, call
from multiprocessing import Pool, cpu_count
from sklearn import preprocessing


# from functools import reduce


def timestamp():
    return str(datetime.datetime.now()).split('.')[0]


def safe_remove(file):
    if os.path.exists(file):
        os.remove(file)


def safe_bedsort(input, output):
    call("sort -k1,1 -k2,2n " + input + " > " + output, shell=True)


def get_uid():
    '''

    Returns:
        8 digit unique identifier in string

    '''
    return str(uuid.uuid4())[:8]


def make_bin(prefix, chromSize, binLength, stepSize, blackList):
    '''
    makes genomic bins in BED format

    Args:
        prefix: prefix for the output
        chromSize: chrom size file
        binLength: genomic bin width
        stepSize: sliding window step
        blackList: blacklist in BED format

    Returns:
        None

    '''
    ### make sliding window
    print("[%s] Making bins" % (timestamp()))
    bin = pybedtools.BedTool().window_maker(g=chromSize, w=binLength, s=stepSize)

    ### filter blacklist region
    print("[%s] Filtering blacklist region" % (timestamp()))
    blk = pybedtools.BedTool(blackList).sort()
    out = bin.intersect(blk, v=True, sorted=True)

    ### write to file
    with open(prefix + ".bin.bed", 'w') as file:
        file.write(str(out))
    del bin, blk, out
    print("[%s] Done" % (timestamp()))


def proc_cov(prefix, bedFile, bwFiles):
    '''

    processes covariate bigWig files

    Args:
        prefix: prefix for the output file
        bedFile: genomic bin in BED format
        bwFiles: covariates in bigWig format

    Returns:
        None

    '''
    ### average bigwig over bin bed
    print("[%s] Averaging features per bin" % (timestamp()))
    mat = np.zeros(shape=(sum(1 for l in open(bedFile)), len(bwFiles)), dtype=float)
    for j, bw in enumerate(bwFiles):
        print("[%s] Processing %s" % (timestamp(), bw))
        b = pyBigWig.open(bw)
        with open(bedFile, "r") as bed:
            for i, bin in enumerate(bed.readlines()):
                chr, start, end = bin.strip().split("\t")
                val = b.stats(chr, int(start), int(end), type="mean")
                if isinstance(val[0], float):
                    mat[i][j] = val[0]
        b.close()
    np.savetxt(prefix + ".cov.tsv", mat, fmt='%.2f', delimiter="\t")
    del mat
    print("[%s] Done" % (timestamp()))


def list_chr(chromSize):
    with open(chromSize, "r") as cs:
        return [c.split("\t")[0] for c in cs]


def count_total_mapped_reads(bam):
    idxstats_by_line = [l.split("\t") for l in pysam.idxstats(bam).split("\n")]
    idxstats_by_line_clean = filter(lambda x: len(x) == 4, idxstats_by_line)
    return reduce(lambda x, y: x + y, [int(count_by_chr[2]) for count_by_chr in idxstats_by_line_clean])


def count_total_proper_templates(bam, minSize, maxSize):
    if not os.path.exists(bam + ".bai"):
        print("[%s] (Warning) Index not found: %s" % (timestamp(), bam))
        print("[%s] Indexing %s" % (timestamp(), bam))
        pysam.index(bam)
    b = pysam.AlignmentFile(bam, "rb")

    proper_pair_count = 0
    chimeric_count = 0
    template_count = 0
    proper_template_count = 0

    for read in b.fetch():
        ### read is in proper pair
        ### read is NOT chimeric read (i.e., no SA tag)
        ### read is mapped to forward strand, mate is mapped to reverse strand
        if read.is_proper_pair:
            proper_pair_count += 1
            if read.has_tag("SA"):
                chimeric_count += 1
            else:
                if not read.is_reverse and read.mate_is_reverse:
                    template_count += 1
                    if read.template_length >= int(minSize) and read.template_length <= int(maxSize):
                        proper_template_count += 1
    b.close()
    return proper_template_count


def bam_proc_worker_se(args):
    bam, chr, fid, minSize, maxSize, strand, readStart = args
    print("[%s] Processing %s" % (timestamp(), chr))

    b = pysam.AlignmentFile(bam, "rb")

    read_count_all = 0
    read_count_fwd = 0
    read_count_rev = 0
    read_count_used_all = 0
    read_count_used_fwd = 0
    read_count_used_rev = 0

    with open("tmp" + fid + chr + ".bed", "w") as s, open("tmp" + fid + chr + ".bpCount.bed", "w") as bpCount:
        for read in b.fetch(reference=chr):

            rid = read.query_name

            ### read is NOT duplicated
            if not read.is_duplicate:

                read_count_all += 1

                # r_c = b.get_reference_name(read.reference_id)
                r_s = read.reference_start
                r_e = read.reference_end

                if not read.is_reverse:  ### read FWD
                    read_count_fwd += 1

                    if strand.lower() != "rev":
                        read_count_used_all += 1
                        read_count_used_fwd += 1

                        s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r_s, r_e, rid, "+"))
                        bpCount.write(chr + '\t' + str(r_s) + '\t' + str(int(r_s) + 1) + '\n')

                elif read.is_reverse:  ### read REV
                    read_count_rev += 1

                    if strand.lower() != "fwd":
                        read_count_used_all += 1
                        read_count_used_rev += 1

                        s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r_s, r_e, rid, "-"))
                        bpCount.write(chr + '\t' + str(r_s) + '\t' + str(int(r_s) + 1) + '\n')

                else:
                    print("[%s] (Warning) read pair mapped to different chromosomes: %s" % (timestamp(), rid))

    b.close()

    safe_bedsort("tmp" + fid + chr + ".bed", "tmp" + fid + chr + ".sorted.bed")
    safe_remove("tmp" + fid + chr + ".bed")

    safe_bedsort("tmp" + fid + chr + ".bpCount.bed", "tmp" + fid + chr + ".bpCount.sorted.bed")
    safe_remove("tmp" + fid + chr + ".bpCount.bed")

    print("[%s] Finished processing %s" % (timestamp(), chr))
    return (read_count_all, read_count_fwd, read_count_rev, read_count_used_all, read_count_used_fwd, read_count_used_rev)

def bam_proc_worker(args):
    bam, chr, fid, minSize, maxSize, strand, readStart = args
    print("[%s] Processing %s" % (timestamp(), chr))

    b = pysam.AlignmentFile(bam, "rb")

    template_count_all = 0
    template_count_fwd = 0
    template_count_rev = 0
    template_count_used_all = 0
    template_count_used_fwd = 0
    template_count_used_rev = 0

    with open("tmp" + fid + chr + ".bed", "w") as s, open("tmp" + fid + chr + ".bpCount.bed", "w") as bpCount:
        r_cache = {}
        for read in b.fetch(reference=chr):

            rid = read.query_name

            ### read is NOT duplicated
            ### read IS first in pair
            ### read IS "properly paired" (read is mapped to forward strand, mate is mapped to reverse strand, and vice versa)
            ### read is NOT chimeric read (i.e., no SA tag)
            ### IF read has SA tag (SA: secondary alignment), read is ambiguous, and thus discard
            if not read.is_duplicate and read.is_proper_pair and not read.has_tag("SA"):

                if rid in r_cache:
                    template_count_all += 1

                    if read.is_read2:
                        read1 = r_cache[rid]
                        read2 = read
                    else:
                        read1 = read
                        read2 = r_cache[rid]

                    del r_cache[rid]

                    r1_c = b.get_reference_name(read1.reference_id)
                    r1_s = read1.reference_start
                    r1_e = read1.reference_end
                    r1_l = read1.template_length

                    r2_c = b.get_reference_name(read2.reference_id)
                    r2_s = read2.reference_start
                    r2_e = read2.reference_end
                    r2_l = read2.template_length

                    if r1_c == r2_c:

                        if not read1.is_reverse and read2.is_reverse:  ### read1 FWD read2 REV
                            template_count_fwd += 1

                            if int(minSize) <= r1_l and int(maxSize) >= r1_l and strand.lower() != "rev":
                                if (r2_e - r1_s) != r1_l:
                                    print("[%s] (Warning) reads are not properly paired: %s" % (timestamp(), rid))
                                else:
                                    template_count_used_all += 1
                                    template_count_used_fwd += 1

                                    if readStart:
                                        if strand.lower() == "fwd":
                                            s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r1_s, r1_e, rid, "+"))
                                        else:
                                            s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r1_s, r1_e, rid, "+"))
                                            s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r2_s, r2_e, rid, "-"))
                                        bpCount.write(chr + '\t' + str(r1_s) + '\t' + str(int(r1_s) + 1) + '\n')
                                        bpCount.write(chr + '\t' + str(int(r2_e) - 1) + '\t' + str(r2_e) + '\n')
                                    else:
                                        s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r1_s, r2_e, rid, "+"))
                                        bpCount.write(chr + '\t' + str(int(r1_s) + int((int(r2_e) - int(r1_s)) / 2)) + '\t' + str(int(r1_s) + int((int(r2_e) - int(r1_s)) / 2) + 1) + '\n')

                        elif read1.is_reverse and not read2.is_reverse:  ### read1 REV read2 FWD
                            template_count_rev += 1

                            if int(minSize) <= r2_l and int(maxSize) >= r2_l and strand.lower() != "fwd":
                                if (r1_e - r2_s) != r2_l:
                                    print("[%s] (Warning) reads are not properly paired: %s" % (timestamp(), rid))
                                else:
                                    template_count_used_all += 1
                                    template_count_used_rev += 1

                                    if readStart:
                                        if strand.lower() == "rev":
                                            s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r2_s, r2_e, rid, "-"))
                                        else:
                                            s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r1_s, r1_e, rid, "+"))
                                            s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r2_s, r2_e, rid, "-"))
                                        bpCount.write(chr + '\t' + str(r2_s) + '\t' + str(int(r2_s) + 1) + '\n')
                                        bpCount.write(chr + '\t' + str(int(r1_e) - 1) + '\t' + str(r1_e) + '\n')
                                    else:
                                        s.write("%s\t%i\t%i\t%s\t.\t%s\n" % (chr, r2_s, r1_e, rid, "-"))
                                        bpCount.write(chr + '\t' + str(int(r2_s) + int((int(r1_e) - int(r2_s)) / 2)) + '\t' + str(int(r2_s) + int((int(r1_e) - int(r2_s)) / 2) + 1) + '\n')

                    else:
                        print("[%s] (Warning) read pair mapped to different chromosomes: %s" % (timestamp(), rid))

                else:
                    r_cache[rid] = read

    b.close()

    safe_bedsort("tmp" + fid + chr + ".bed", "tmp" + fid + chr + ".sorted.bed")
    safe_remove("tmp" + fid + chr + ".bed")

    safe_bedsort("tmp" + fid + chr + ".bpCount.bed", "tmp" + fid + chr + ".bpCount.sorted.bed")
    safe_remove("tmp" + fid + chr + ".bpCount.bed")

    del r_cache
    print("[%s] Finished processing %s" % (timestamp(), chr))
    return (template_count_all, template_count_fwd, template_count_rev, template_count_used_all, template_count_used_fwd, template_count_used_rev)


def proc_bam(prefix, chromSize, bedFile, bamFiles, minSize=200, maxSize=1000, readStart=False, strand="all", singleEnd=False):
    '''

    processes alignments in BAM format

    Args:
        prefix: prefix for the output
        chromSize: chrom size file
        bedFile: bin BED file
        bamFiles: list of BAM files eg. [input.bam output.bam]
        minSize: minimum size of fragment insert to consider (default 200)
        maxSize: maximum size of fragment insert to consider (default 1000)
        readStart: count at start positions of reads
        strand: use all/fwd/rev stranded fragments

    Returns:
        writes bin count output file

    '''

    print("[%s] Counting feature depth (template or read) per genomic bin %s" % (timestamp(), bedFile))
    if readStart:
        print("[%s] Using read start position" % (timestamp()))
    else:
        print("[%s] Using fragment center" % (timestamp()))

    ### initialize numpy array
    tct = np.zeros(shape=(len(bamFiles)), dtype=int)
    mat = np.zeros(shape=(sum(1 for l in open(bedFile)), len(bamFiles)), dtype=int)

    ### random unique ID
    uid = get_uid()

    ### load bin bed file
    a = pybedtools.BedTool(bedFile)

    bamidx = {0: "input", 1: "output"}

    for j, bam in enumerate(bamFiles):
        print("[%s] Processing %s alignment file %s" % (timestamp(), bamidx[j], bam))

        if not os.path.exists(bam + ".bai"):
            print("[%s] (Warning) Index not found: %s" % (timestamp(), bam))
            print("[%s] Indexing %s" % (timestamp(), bam))
            pysam.index(bam)

        b = pysam.AlignmentFile(bam, "rb")
        total_mapped_reads = b.mapped
        b.close()

        print("[%s] Parallel processing using %i cores" % (timestamp(), cpu_count()))
        args_list = [(bam, c, uid + bamidx[j], minSize, maxSize, strand, readStart) for c in list_chr(chromSize)]
        pool = Pool(processes=cpu_count())
        if singleEnd:
            print("[%s] Using single-end mode (SE mode uses read start position and min-max options are ignored)" % (timestamp()))
            list_counts = pool.map(bam_proc_worker_se, args_list)
        else:
            list_counts = pool.map(bam_proc_worker, args_list)

        template_count_all, template_count_fwd, template_count_rev, template_count_used_all, template_count_used_fwd, template_count_used_rev = np.sum(list_counts, axis=0)

        tct[j] += template_count_used_all

        print("[%s] %s total mapped reads" % (timestamp(), '{:,}'.format(total_mapped_reads)))

        if strand.lower() == "fwd":
            print("[%s] Counting fragments mapped to forward (+) strand only" % (timestamp()))
            print("[%s] %s templates extracted (+)" % (timestamp(), '{:,}'.format(template_count_fwd)))
            print("[%s] %s templates used for count (+)" % (timestamp(), '{:,}'.format(template_count_used_fwd)))

        elif strand.lower() == "rev":
            print("[%s] Counting fragments mapped to reverse (-) strand only" % (timestamp()))
            print("[%s] %s templates extracted (-)" % (timestamp(), '{:,}'.format(template_count_rev)))
            print("[%s] %s templates used for count (-)" % (timestamp(), '{:,}'.format(template_count_used_rev)))

        else:
            print("[%s] Counting fragments mapped to both forward (+) and reverse (-) strand" % (timestamp()))
            print("[%s] %s templates extracted" % (timestamp(), '{:,}'.format(template_count_all)))
            print("[%s] %s templates used for count" % (timestamp(), '{:,}'.format(template_count_used_all)))

        print("[%s] Merging fragment bed files" % (timestamp()))
        with open(prefix + "." + bamidx[j] + ".frag.bed", "w") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + bamidx[j] + chr + ".sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + bamidx[j] + chr + ".sorted.bed")

        ### save genome coverage
        print("[%s] Making genome coverage bedGraph" % (timestamp()))
        mergedBed = pybedtools.BedTool(prefix + "." + bamidx[j] + ".frag.bed")
        mergedBed.genome_coverage(bg=True, g=chromSize).saveas(prefix + "." + bamidx[j] + ".bdg")

        ### convert bedGraph to bigWig
        print("[%s] Converting bedGraph to bigWig" % (timestamp()))
        bdg2bw(bdgFile=prefix + "." + bamidx[j] + ".bdg", bwFile=prefix + "." + bamidx[j] + ".bw", chromSize=chromSize)
        safe_remove(prefix + "." + bamidx[j] + ".bdg")

        ### merge sorted bpCount files
        print("[%s] Merging count bed files" % (timestamp()))
        with open(prefix + "." + bamidx[j] + ".bpCount.bed", "w") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + bamidx[j] + chr + ".bpCount.sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + bamidx[j] + chr + ".bpCount.sorted.bed")

        ### count fragments per bin
        print("[%s] Counting depth per bin" % (timestamp()))
        bpCountBed = pybedtools.BedTool(prefix + "." + bamidx[j] + ".bpCount.bed")
        readDepth = a.coverage(bpCountBed, sorted=True, counts=True)

        ### extract 4th column, which is read counts, and assign as numpy array
        mat[:, j] = np.array([int(l.split("\t")[3]) for l in str(readDepth).rstrip("\n").split("\n")])

        ### clean up
        safe_remove(prefix + "." + bamidx[j] + ".bpCount.bed")
        del mergedBed, bpCountBed, readDepth

    ### normalize input count, normalized input count is added to additional column
    print("[%s] Normalizing factor for input: %f" % (timestamp(), (tct[1] / tct[0])))
    normalized_input = mat[:, 0] * (tct[1] / tct[0])
    np.savetxt(prefix + ".bam.bct", np.concatenate((mat, normalized_input.reshape(-1, 1)), axis=1), fmt='%i %i %.5f', delimiter="\t")

    del a, mat, tct, normalized_input
    print("[%s] Done" % (timestamp()))


def trigamma(x):
    return polygamma(1, x)


def score(th, mu, y, w):
    return sum(w * (digamma(th + y) - digamma(th) + np.log(th) + 1 - np.log(th + mu) - (y + th) / (mu + th)))


def info(th, mu, y, w):
    return sum(w * (-trigamma(th + y) + trigamma(th) - 1 / th + 2 / (mu + th) - (y + th) / (mu + th) ** 2))


def theta(y, mu, verbose=False):
    ### MLE for theta and std. error

    ### stop iteration if delta smaller than eps
    eps = np.finfo(np.float).eps ** 0.25

    ### max iter
    limit = 20

    ### init
    weights = np.repeat(1, len(y))
    n = sum(weights)
    t0 = n / sum(weights * (y / mu - 1) ** 2)
    it = 0
    de = 1

    if (verbose): print("theta: iter %d theta = %f" % (it, t0))

    while (it < limit and abs(de) > eps):
        it += 1
        t0 = abs(t0)
        de = score(t0, mu, y, weights) / info(t0, mu, y, weights)
        t0 = t0 + de
        if (verbose): print("theta: iter %d theta = %f" % (it, t0))

    ### warning
    if (t0 < 0):
        t0 = 0
        print("warning: estimate truncated at zero")
    if (it == limit):
        print("warning: iteration limit reached")

    ### standard error
    se = np.sqrt(1 / info(t0, mu, y, weights))

    return t0, se


def call_peak(prefix, bedFile, bctFile, chromSize, bwFile, covFile=None, threshold=0.05, mode=1, minCoverage=10, extQuantile=1e-5):
    '''

    calls peak

    Args:
        (Required)
        prefix: prefix for the output
        bedFile: bin in BED format
        bctFile: fragment insert coverage in BCT format (made from procBam.py)
        chromSize: chromosome sizes
        bwFile: fragment coverage in BigWig format

        (Optional)
        covFile: covariates in TSV format
        threshold: threshold to call peak
        mode: 1 - using input as covariate 2 - using input as offset
        minCoverage: minimum coverage required for peak
        extQuantile: for removing genomic bins having extreme quantile (sequencing artifact)

    Returns:
        peak files (original peaks and center-merged final peaks)
        peak format: chr, start, end, foldChg, inputCt, outputCt, pval, qval

    '''
    print("[%s] Calling peaks" % (timestamp()))

    ### load data
    print("[%s] Loading fragment coverages" % (timestamp()))
    bct = np.loadtxt(bctFile, ndmin=2)  # BCT 0=input, 1=output, 2=normalized input

    if covFile:
        print("[%s] Loading covariates" % (timestamp()))
        cov = np.loadtxt(covFile, ndmin=2)  # COV 3:=cov

        ### scale covariates to have mean 0 and sd 1
        cov_scaled = preprocessing.scale(cov, axis=0)

        ### merge data
        mat = np.concatenate((bct[:, [1, 0, 2]], cov_scaled), axis=1)  # MAT 0=output, 1=input, 2=normalized input, 3:=cov

        del cov, cov_scaled
    else:
        print("[%s] Running without covariates" % (timestamp()))
        ### merge data
        mat = bct[:, [1, 0, 2]]  # MAT 0=output, 1=input, 2=normalized input

    del bct

    bct_o = mat[:, 0]
    bct_i = mat[:, 1]
    bct_n = mat[:, 2]

    ### non sliding bins
    nonSliding = np.zeros(mat.shape[0], dtype=bool)  ### initialize with False
    with open(bedFile, "r") as bed:
        lastchr, lastbin = "", 0
        for i, bin in enumerate(bed.readlines()):
            if bin.split("\t")[0] != lastchr:
                lastchr = bin.split("\t")[0]
                lastbin = int(bin.split("\t")[2])
                nonSliding[i] = True
            elif int(bin.split("\t")[1]) >= lastbin:
                lastbin = int(bin.split("\t")[2])
                nonSliding[i] = True

    ### remove bins with input count of zero (i.e., untested region) OR extreme values (top 0.001%, i.e., sequencing artifacts)

    minInput = minCoverage
    maxInput = np.quantile(bct_i[bct_i > 0], (1 - extQuantile))
    # print(minInput, maxInput)

    minOutput = minCoverage
    maxOutput = np.quantile(bct_o[bct_o > 0], (1 - extQuantile))
    # print(minOutput, maxOutput)

    ### calculate fold change
    fc = np.zeros(mat.shape[0])
    fc[bct_n > 0] = bct_o[bct_n > 0] / (bct_n[bct_n > 0])  ### fc = output / normalized_input

    ### train / test genomic bin ###
    trainingBin = (bct_i > minInput) & (bct_i < maxInput) & (bct_o > minOutput) & (bct_o < maxOutput) & nonSliding
    testingBin = (bct_i > minInput) & (bct_i < maxInput) & (fc > 1.5)  ### bins w/ FC > 1.5 are tested for statistical significance

    ### filtering bins
    print("[%s] Total genomic bins: %s" % (timestamp(), '{:,}'.format(mat.shape[0])))
    print("[%s] Total non-overlapping genomic bins: %s" % (timestamp(), '{:,}'.format(sum(nonSliding))))

    print("[%s] Removing bins with input counts larger than %s for training and testing" % (timestamp(), '{:,}'.format(maxInput)))
    print("[%s] Removing bins with output counts larger than %s for training" % (timestamp(), '{:,}'.format(maxOutput)))

    print("[%s] Total genomic bins used for training: %s" % (timestamp(), '{:,}'.format(mat[trainingBin, :].shape[0])))
    print("[%s] Total genomic bins used for testing: %s" % (timestamp(), '{:,}'.format(mat[testingBin, :].shape[0])))

    ### mode 2 uses "input" as offset variable
    if int(mode) == 2:
        print("[%s] Running Mode 2" % (timestamp()))

        ### remove input
        mat_model = np.delete(mat, 1, 1)

        ### formula
        x = ["x" + str(i) for i in range(1, mat_model.shape[1] - 1)]
        df = pd.DataFrame(mat_model[trainingBin, :], columns=["y", "exposure"] + x)
        formula = "y~" + "+".join(df.columns.difference(["y", "exposure"]))
        print("[%s] Fit using formula: %s" % (timestamp(), formula))

        ### Initial parameter estimation using Poisson regression
        # print("[%s] Initial estimate" % (timestamp()))
        model0 = smf.glm(formula, data=df, family=sm.families.Poisson(), offset=np.log(df["exposure"])).fit()
        # print model0.summary()

        ### Estimate theta
        th0, _ = theta(mat_model[trainingBin, :][:, 0], model0.mu)
        print("[%s] Initial estimate of theta is %f" % (timestamp(), th0))

        ### re-estimate beta with theta
        # print("[%s] Re-estimate of beta" % (timestamp()))
        model = smf.glm(formula, data=df, family=sm.families.NegativeBinomial(alpha=1 / th0), offset=np.log(df["exposure"])).fit(start_params=model0.params)
        # print model.summary()

        ### Re-estimate theta
        th, _ = theta(mat_model[trainingBin, :][:, 0], model.mu)
        print("[%s] Re-estimate of theta is %f" % (timestamp(), th))

        ### predict
        print("[%s] Predicting expected counts" % (timestamp()))

        df = pd.DataFrame(mat_model[testingBin, :], columns=["y", "exposure"] + x)
        y_hat = model.predict(df, offset=np.log(df["exposure"]))

    ### mode 1 uses "input" as covariate (default):
    else:
        print("[%s] Running Mode 1" % (timestamp()))

        ### remove normalized input
        mat_model = np.delete(mat, 2, 1)

        ### formula
        x = ["x" + str(i) for i in range(1, mat_model.shape[1])]
        df = pd.DataFrame(mat_model[trainingBin, :], columns=["y"] + x)
        formula = "y~" + "+".join(df.columns.difference(["y"]))
        print("[%s] Fit using formula: %s" % (timestamp(), formula))

        ### Initial parameter estimation using Poisson regression
        # print("[%s] Initial estimate" % (timestamp()))
        model0 = smf.glm(formula, data=df, family=sm.families.Poisson()).fit()
        # print model0.summary()

        ### Estimate theta
        th0, _ = theta(mat_model[trainingBin, :][:, 0], model0.mu)
        print("[%s] Initial estimate of theta is %f" % (timestamp(), th0))

        ### re-estimate beta with theta
        # print("[%s] Re-estimate of beta" % (timestamp()))
        model = smf.glm(formula, data=df, family=sm.families.NegativeBinomial(alpha=1 / th0)).fit(start_params=model0.params)
        # print model.summary()

        ### Re-estimate theta
        th, _ = theta(mat_model[trainingBin, :][:, 0], model.mu)
        print("[%s] Re-estimate of theta is %f" % (timestamp(), th))

        ### predict
        print("[%s] Predicting expected counts" % (timestamp()))
        df = pd.DataFrame(mat_model[testingBin, :], columns=["y"] + x)
        y_hat = model.predict(df)

    ### calculate P-value
    print("[%s] Calculating P-value" % (timestamp()))
    theta_hat = np.repeat(th, len(y_hat))
    prob = th / (th + y_hat)  ### prob=theta/(theta+mu)
    pval = 1 - nbinom.cdf(k=mat_model[testingBin, 0] - 1, n=theta_hat, p=prob)

    ### multiple testing correction
    print("[%s] Multiple testing correction" % (timestamp()))
    _, qval, _, _ = multi.multipletests(pval, method="fdr_bh")

    nlog10pval = -np.log10(pval)
    nlog10qval = -np.log10(qval)

    ### output initial peaks
    print("[%s] Output initial peaks" % (timestamp()))
    with open(prefix + ".peak.bed", "w") as out:
        with open(bedFile, "r") as bed:
            for i, bin in enumerate(list(compress(bed.readlines(), testingBin))):
                if qval[i] <= float(threshold):
                    ### chr, start, end, foldChg, inputCov, outputCov, -log10(pval), -log10(qval)
                    out.write("%s\t%.3f\t%i\t%i\t%.3f\t%.3f\n" % (bin.strip(), fc[testingBin][i], mat[testingBin, 1][i], mat[testingBin, 0][i], nlog10pval[i], nlog10qval[i]))

    ### generate various bdg tracks
    print("[%s] Generating bedGraph files" % (timestamp()))
    with open(prefix + ".pval.bdg", "w") as fp, open(prefix + ".qval.bdg", "w") as fq, open(prefix + ".fc.bdg", "w") as ff:
        with open(bedFile, "r") as bed:
            for i, bin in enumerate(bed.readlines()):
                ff.write("%s\t%.3f\n" % (bin.strip(), fc[i]))

        with open(bedFile, "r") as bed:
            for i, bin in enumerate(list(compress(bed.readlines(), testingBin))):
                fp.write("%s\t%.3f\n" % (bin.strip(), nlog10pval[i]))
                fq.write("%s\t%.3f\n" % (bin.strip(), nlog10qval[i]))

    del mat, mat_model, nlog10pval, nlog10qval, pval, qval

    ### make bigWig track
    print("[%s] Making BigWig tracks" % (timestamp()))

    with open(bedFile) as bed:
        c1, s1, e1 = bed.readline().strip().split('\t')
        c2, s2, e2 = bed.readline().strip().split('\t')

    w = int(e1) - int(s1)
    s = int(s2) - int(s1)

    bdg2bw(bdgFile=prefix + ".fc.bdg", bwFile=prefix + ".fc.bw", chromSize=chromSize, window=w, step=s)
    safe_remove(prefix + ".fc.bdg")

    bdg2bw(bdgFile=prefix + ".pval.bdg", bwFile=prefix + ".pval.bw", chromSize=chromSize, window=w, step=s)
    safe_remove(prefix + ".pval.bdg")

    bdg2bw(bdgFile=prefix + ".qval.bdg", bwFile=prefix + ".qval.bw", chromSize=chromSize, window=w, step=s)
    safe_remove(prefix + ".qval.bdg")

    ### center peak
    print("[%s] Center peaks" % (timestamp()))
    center_peak(bwFile=bwFile, peakFile=prefix + ".peak.bed", centeredPeakFile=prefix + ".peak.centered.bed")

    ### merge peak
    print("[%s] Merge peaks" % (timestamp()))
    # pybedtools.BedTool(prefix + ".peak.centered.bed").merge(c=[4, 6, 7, 8], o=["max", "max", "max", "max"]).saveas(prefix + ".peak.centered.merged.bed")
    pybedtools.BedTool(prefix + ".peak.centered.bed").merge(c=[4, 5, 6, 7, 8], o=["max", "max", "max", "max", "max"]).saveas(prefix + ".peak.centered.merged.bed")

    ### finalize peak
    print("[%s] Finalize peaks" % (timestamp()))
    peak = pd.read_csv(prefix + ".peak.centered.merged.bed", sep='\t', header=None)
    # peak.columns = ['chr', 'start', 'end', 'fc', 'cov', 'nlog10pval', 'nlog10qval']
    peak.columns = ['chr', 'start', 'end', 'fc', 'icov', 'ocov', 'nlog10pval', 'nlog10qval']
    peak['log2fc'] = np.log2(peak['fc'])
    peak['idx'] = peak.index
    peak['strand'] = "."
    peak['score'] = [min(int(round(x)), 1000) for x in peak['fc'] * 100]
    peak = peak.sort_values(by=['fc'], ascending=False)
    peak['name'] = ['peak_' + str(x) for x in peak.reset_index().index + 1]
    peak = peak.sort_values(by=['idx'])
    # final = peak[['chr', 'start', 'end', 'name', 'score', 'strand', 'fc', 'cov', 'nlog10pval', 'nlog10qval']]
    final = peak[['chr', 'start', 'end', 'name', 'score', 'strand', 'log2fc', 'icov', 'ocov', 'nlog10pval', 'nlog10qval']]
    final.to_csv(prefix + '.peak.final.bed', sep='\t', float_format='%.3f', index=False, header=False)

    ### remove intermediate peak files
    safe_remove(prefix + ".peak.centered.bed")
    safe_remove(prefix + ".peak.centered.merged.bed")

    print("[%s] Done" % (timestamp()))


def bdg2bw(bdgFile, bwFile, chromSize, window=None, step=None):
    print("[%s] Writing %s" % (timestamp(), bwFile))
    with open(chromSize) as f:
        cs = [line.strip().split('\t') for line in f.readlines()]

    bw = pyBigWig.open(bwFile, "w")
    bw.addHeader([(str(x[0]), int(x[1])) for x in cs])

    if window and step:
        print("[%s] Making bigWig using fixed interval size of %i" % (timestamp(), step))
        with open(bdgFile, "r") as bdg:
            for line in bdg:
                if len(line.strip().split("\t")) == 4:
                    chr, start, end, val = line.strip().split("\t")
                    if (int(end) - int(start)) == window:
                        bw.addEntries(chroms=[chr], starts=[(int(end) - int(window / 2) - int(step / 2))], ends=[(int(end) - int(window / 2) + int(step / 2))], values=[float(val)])
                else:
                    print("[%s] Warning: skipping bedGraph entry: %s" % (timestamp(), line.strip()))

    else:
        with open(bdgFile, "r") as bdg:
            for line in bdg:
                if len(line.strip().split("\t")) == 4:
                    chr, start, end, val = line.strip().split("\t")
                    bw.addEntries(chroms=[chr], starts=[int(start)], ends=[int(end)], values=[float(val)])
                else:
                    print("[%s] Warning: skipping bedGraph entry: %s" % (timestamp(), line.strip()))

    bw.close()


def center_peak(bwFile, peakFile, centeredPeakFile):
    '''

    centers peak file

    Args:
        bwFile: bigWig file
        peakFile: peak file in BED format
        centeredPeakFile: centered peak file name

    Returns:
        None

    '''
    bw = pyBigWig.open(bwFile)

    with open(centeredPeakFile, "w") as out:
        for p in pybedtools.BedTool(peakFile):
            chr = p[0]
            start = int(p[1])
            end = int(p[2])
            radius = int((end - start) / 2)
            other = '\t'.join(p[3:])

            # find output coverage for peak
            interval = bw.intervals(chr, start, end)

            if len(interval) == 0:
                print("[%s] Warning! No Intersect Found for peak: %s" % (timestamp(), '\t'.join(p)))
                out.write('\t'.join([chr, str(start), str(end), other]) + '\n')
            # elif (end - start) < int(windowSize):
            #     print("[%s] Warning! Smaller Peak Found: %s" % (timestamp(), '\t'.join(p)))
            else:
                cov = np.array(interval, dtype=int)
                summit = cov[np.nonzero(cov[:, 2] == np.max(cov[:, 2]))]
                center = int((summit[0, 0] + summit[-1, 1]) / 2)
                out.write('\t'.join([chr, str(center - radius), str(center + radius), other]) + '\n')

    bw.close()


def proc_fenergy(bedFile, fileOut, linearfold, genome):
    '''

    calculates folding free energy in parallel

    Args:
        bedFile: genomic bins in BED format
        fileOut: output file
        linearfold: path to linearfold
        genome: reference genome in fasta format

    Returns:
        None

    '''
    print("[%s] Calculate Folding Free Energy per bin" % (timestamp()))

    ### random unique ID
    uid = get_uid()

    # chr_seq=split_bed(bedFile,uid)
    nrow = 10000
    part_list = []

    ### split BED file
    with open(bedFile, 'r') as bin:
        part, part_num = None, 0
        for idx, line in enumerate(bin.readlines()):
            if idx % nrow == 0:
                if part:
                    part.close()
                part_num += 1
                part = open("tmp" + uid + "_" + str(part_num) + ".bed", 'w')
                part_list.append(part_num)
            part.write(line)
        if part:
            part.close()

    ### make arg for parallel computing
    prefixes = []
    for part_num in part_list:
        prefixes.append(["tmp" + uid + "_" + str(part_num), linearfold, genome])

    ### parallel compute linearfold
    print("[%s] Parallel computing using %i cores" % (timestamp(), cpu_count()))
    pool = Pool(processes=cpu_count())
    pool.map(run_linearfold, prefixes)

    ### merge bed
    print("[%s] Merging files" % (timestamp()))
    safe_remove(fileOut)
    with open(fileOut, "w") as merged:
        for part_num in part_list:

            ### merge tmp bed files
            with open("tmp" + uid + "_" + str(part_num) + ".out", "r") as t:
                if t.read(1).strip():
                    t.seek(0)
                    merged.write(t.read())

            ### delete tmp bed files
            safe_remove("tmp" + uid + "_" + str(part_num) + ".out")


def run_linearfold(arg):
    '''

    Args:
        list of arguments
        [fasta, path2linearfold, genome]

    Returns:
        NA

    '''
    pybedtools.BedTool(arg[0] + ".bed").sequence(fi=arg[2]).save_seqs(arg[0] + ".fa")
    safe_remove(arg[0] + ".bed")
    with open(arg[0] + ".fa", 'r') as fa:
        out = check_output([arg[1], "-V"], stdin=fa)
    with open(arg[0] + ".out", 'w') as fo:
        for line in out.split("\n"):
            line_split = line.strip().split(" ")
            if len(line_split) > 1:
                fo.write(line_split[1].replace("(", "").replace(")", "") + "\n")
    safe_remove(arg[0] + ".fa")


def split_bed(bedFile, uid):
    '''
    split BED file by chromosomes

    Args:
        bedFile: input BED file
        uid: random unique idenfifier

    Returns:
        list of chromosomes

    '''
    chr = ""
    chr_list = []
    part = None
    with open(bedFile, 'r') as bin:
        for line in bin.readlines():
            line_split = line.strip().split('\t')

            if line_split[0] == chr:
                part.write(line)
            else:
                if part:
                    part.close()
                chr = line_split[0]
                chr_list.append(chr)
                part = open("tmp" + uid + chr + ".bed", 'w')
                part.write(line)
        if part:
            part.close()
    return chr_list


def proc_bam_legacy(bamFiles, bedFile, chromSize, fileOut, minSize, maxSize, readStart=False):
    '''
    process bam file legacy code, retained for backward compatibility
    use 'proc_bam' instead

    Args:
        bamFiles: list of BAM files eg. [input.bam output.bam]
        bedFile: bin BED file
        chromSize: chrom size file
        fileOut: output file
        minSize: minimum size of fragment insert to consider
        maxSize: maximum size of fragment insert to consider
        readStart: count at start positions of reads

    Returns:
        writes bin count output file

    '''
    print("[%s] Counting template depth per bin %s" % (timestamp(), bedFile))

    ### initialize numpy array
    tct = np.zeros(shape=(len(bamFiles)), dtype=int)
    mat = np.zeros(shape=(sum(1 for l in open(bedFile)), len(bamFiles)), dtype=int)

    ### random unique ID
    uid = get_uid()

    ### load bin bed file
    a = pybedtools.BedTool(bedFile)

    bamidx = {0: "input", 1: "output"}

    for j, bam in enumerate(bamFiles):
        print("[%s] Processing %s" % (timestamp(), bam))

        if not os.path.exists(bam + ".bai"):
            print("[%s] (Warning) Index not found: %s" % (timestamp(), bam))
            print("[%s] Indexing %s" % (timestamp(), bam))
            pysam.index(bam)

        b = pysam.AlignmentFile(bam, "rb")

        template_count_all = 0
        template_count_fwd = 0
        template_count_rev = 0
        template_count_used_all = 0
        template_count_used_fwd = 0
        template_count_used_rev = 0

        ###

        for chr in list_chr(chromSize):

            print("[%s] Processing %s" % (timestamp(), chr))

            with open("tmp" + uid + bamidx[j] + chr + "all.bed", "w") as s, open("tmp" + uid + bamidx[j] + chr + "fwd.bed", "w") as s_fwd, open("tmp" + uid + bamidx[j] + chr + "rev.bed",
                                                                                                                                                "w") as s_rev:
                r_cache = {}
                for read in b.fetch(reference=chr):

                    rid = read.query_name
                    # r_c = b.get_reference_name(read.reference_id)

                    ### read is NOT duplicated
                    ### read IS first in pair
                    ### read IS "properly paired" (read is mapped to forward strand, mate is mapped to reverse strand, and vice versa)
                    ### read is NOT chimeric read (i.e., no SA tag)
                    ### IF read has SA tag (SA: secondary alignment), read is ambiguous, and thus discard
                    if not read.is_duplicate and read.is_proper_pair and not read.has_tag("SA"):

                        if rid in r_cache:
                            template_count_all += 1

                            if readStart:
                                if read.is_read2:
                                    read1 = r_cache[rid]
                                    read2 = read
                                else:
                                    read1 = read
                                    read2 = r_cache[rid]

                                r1_s = read1.reference_start
                                r1_e = read1.reference_end
                                r1_l = read1.template_length

                                r2_s = read2.reference_start
                                r2_e = read2.reference_end
                                r2_l = read2.template_length

                                if read2.is_reverse:  ### read1 FWD read2 REV
                                    template_count_fwd += 1

                                    if int(minSize) <= r1_l and int(maxSize) >= r1_l:
                                        template_count_used_all += 1
                                        template_count_used_fwd += 1
                                        s.write("%s\t%i\t%i\n" % (chr, r1_s, (r1_s + 1)))
                                        s.write("%s\t%i\t%i\n" % (chr, (r2_e - 1), r2_e))
                                        s_fwd.write("%s\t%i\t%i\n" % (chr, r1_s, (r1_s + 1)))
                                        s_rev.write("%s\t%i\t%i\n" % (chr, (r2_e - 1), r2_e))

                                elif read1.is_reverse:  ### read1 REV read2 FWD
                                    template_count_rev += 1

                                    if int(minSize) <= r2_l and int(maxSize) >= r2_l:
                                        template_count_used_all += 1
                                        template_count_used_rev += 1

                                        s.write("%s\t%i\t%i\n" % (chr, r2_s, (r2_s + 1)))
                                        s.write("%s\t%i\t%i\n" % (chr, (r1_e - 1), r1_e))

                                        s_fwd.write("%s\t%i\t%i\n" % (chr, r2_s, (r2_s + 1)))
                                        s_rev.write("%s\t%i\t%i\n" % (chr, (r1_e - 1), r1_e))

                            else:
                                if read.is_read2:
                                    read = r_cache[rid]

                                del r_cache[rid]

                                r_s = read.reference_start
                                r_e = read.reference_end
                                r_l = read.template_length

                                if not read.is_reverse and read.mate_is_reverse:
                                    template_count_fwd += 1

                                    if int(minSize) <= r_l and int(maxSize) >= r_l:
                                        template_count_used_all += 1
                                        template_count_used_fwd += 1
                                        s.write("%s\t%i\t%i\n" % (chr, (r_s + int(r_l / 2)), (r_s + int(r_l / 2) + 1)))
                                        s_fwd.write(
                                            "%s\t%i\t%i\n" % (chr, (r_s + int(r_l / 2)), (r_s + int(r_l / 2) + 1)))

                                elif read.is_reverse and not read.mate_is_reverse:
                                    template_count_rev += 1

                                    if int(minSize) <= -(r_l) and int(maxSize) >= -(r_l):
                                        template_count_used_all += 1
                                        template_count_used_rev += 1
                                        s.write("%s\t%i\t%i\n" % (chr, (r_e + int(r_l / 2) - 1), (r_e + int(r_l / 2))))
                                        s_rev.write(
                                            "%s\t%i\t%i\n" % (chr, (r_e + int(r_l / 2) - 1), (r_e + int(r_l / 2))))
                        else:
                            r_cache[rid] = read

            print("[%s] Sorting %s" % (timestamp(), chr))
            pybedtools.BedTool("tmp" + uid + bamidx[j] + chr + "all.bed").sort().saveas("tmp" + uid + bamidx[j] + chr + "all.sorted.bed")
            pybedtools.BedTool("tmp" + uid + bamidx[j] + chr + "fwd.bed").sort().saveas("tmp" + uid + bamidx[j] + chr + "fwd.sorted.bed")
            pybedtools.BedTool("tmp" + uid + bamidx[j] + chr + "rev.bed").sort().saveas("tmp" + uid + bamidx[j] + chr + "rev.sorted.bed")

            ### delete bed
            safe_remove("tmp" + uid + bamidx[j] + chr + "all.bed")
            safe_remove("tmp" + uid + bamidx[j] + chr + "fwd.bed")
            safe_remove("tmp" + uid + bamidx[j] + chr + "rev.bed")

        ###

        tct[j] += template_count_used_all

        print("[%s] %s mapped reads" % (timestamp(), '{:,}'.format(b.mapped)))
        print("[%s] %s templates extracted" % (timestamp(), '{:,}'.format(template_count_all)))
        print("[%s] %s templates extracted (+)" % (timestamp(), '{:,}'.format(template_count_fwd)))
        print("[%s] %s templates extracted (-)" % (timestamp(), '{:,}'.format(template_count_rev)))

        print("[%s] %s templates used for count" % (timestamp(), '{:,}'.format(template_count_used_all)))
        print("[%s] %s templates used for count (+)" % (timestamp(), '{:,}'.format(template_count_used_fwd)))
        print("[%s] %s templates used for count (-)" % (timestamp(), '{:,}'.format(template_count_used_rev)))

        b.close()

        ############################## merge ALL bed ##############################
        print("[%s] Merging ALL BED files" % (timestamp()))
        with open("tmp" + uid + bamidx[j] + "all.sorted.merged.bed", "a") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + bamidx[j] + chr + "all.sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + bamidx[j] + chr + "all.sorted.bed")

        print("[%s] Counting depth per bin" % (timestamp()))
        mergedBed = pybedtools.BedTool("tmp" + uid + bamidx[j] + "all.sorted.merged.bed")
        readDepth = a.coverage(mergedBed, sorted=True, counts=True)

        ### extract 4th column, which is read counts, and assign as numpy array
        mat[:, j] = np.array([int(l.split("\t")[3]) for l in str(readDepth).rstrip("\n").split("\n")])

        ### save genome coverage
        print("[%s] Making genome coverage bedGraph" % (timestamp()))
        binSize = int(a[0][2]) - int(a[0][1])
        mergedBed.slop(g=chromSize, b=int(binSize / 2)).genome_coverage(bg=True, g=chromSize).saveas(fileOut + "." + bamidx[j] + ".all.bdg")

        ### delete tmp merged bed files
        safe_remove("tmp" + uid + bamidx[j] + "all.sorted.merged.bed")
        del merged, mergedBed, readDepth

        ### convert bedGraph to bigWig
        print("[%s] Converting bedGraph to bigWig" % (timestamp()))
        bdg2bw(bdgFile=fileOut + "." + bamidx[j] + ".all.bdg", bwFile=fileOut + "." + bamidx[j] + ".bw", chromSize=chromSize)
        safe_remove(fileOut + "." + bamidx[j] + ".all.bdg")

        ############################## merge FWD bed ##############################
        print("[%s] Merging FWD BED files" % (timestamp()))
        with open("tmp" + uid + bamidx[j] + "fwd.sorted.merged.bed", "a") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + bamidx[j] + chr + "fwd.sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + bamidx[j] + chr + "fwd.sorted.bed")

        print("[%s] Counting depth per bin" % (timestamp()))
        mergedBed = pybedtools.BedTool("tmp" + uid + bamidx[j] + "fwd.sorted.merged.bed")

        ### save genome coverage
        print("[%s] Making genome coverage bedGraph" % (timestamp()))
        binSize = int(a[0][2]) - int(a[0][1])
        mergedBed.slop(g=chromSize, b=int(binSize / 2)).genome_coverage(bg=True, g=chromSize).saveas(fileOut + "." + bamidx[j] + ".fwd.bdg")

        ### delete tmp merged bed files
        safe_remove("tmp" + uid + bamidx[j] + "fwd.sorted.merged.bed")
        del merged, mergedBed

        ### convert bedGraph to bigWig
        print("[%s] Converting bedGraph to bigWig" % (timestamp()))
        bdg2bw(bdgFile=fileOut + "." + bamidx[j] + ".fwd.bdg", bwFile=fileOut + "." + bamidx[j] + ".fwd.bw", chromSize=chromSize)
        safe_remove(fileOut + "." + bamidx[j] + ".fwd.bdg")

        ############################## merge REV bed ##############################
        print("[%s] Merging REV BED files" % (timestamp()))
        with open("tmp" + uid + bamidx[j] + "rev.sorted.merged.bed", "a") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + bamidx[j] + chr + "rev.sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + bamidx[j] + chr + "rev.sorted.bed")

        print("[%s] Counting depth per bin" % (timestamp()))
        mergedBed = pybedtools.BedTool("tmp" + uid + bamidx[j] + "rev.sorted.merged.bed")

        ### save genome coverage
        print("[%s] Making genome coverage bedGraph" % (timestamp()))
        binSize = int(a[0][2]) - int(a[0][1])
        mergedBed.slop(g=chromSize, b=int(binSize / 2)).genome_coverage(bg=True, g=chromSize).saveas(fileOut + "." + bamidx[j] + ".rev.bdg")

        ### delete tmp merged bed files
        safe_remove("tmp" + uid + bamidx[j] + "rev.sorted.merged.bed")
        del merged, mergedBed

        ### convert bedGraph to bigWig
        print("[%s] Converting bedGraph to bigWig" % (timestamp()))
        bdg2bw(bdgFile=fileOut + "." + bamidx[j] + ".rev.bdg", bwFile=fileOut + "." + bamidx[j] + ".rev.bw", chromSize=chromSize)
        safe_remove(fileOut + "." + bamidx[j] + ".rev.bdg")

    ### normalize input count, normalized input count is added to additional column
    print("[%s] Normalizing factor for input: %f" % (timestamp(), (tct[1] / tct[0])))
    normalized_input = mat[:, 0] * (tct[1] / tct[0])
    np.savetxt(fileOut, np.concatenate((mat, normalized_input.reshape(-1, 1)), axis=1), fmt='%i %i %.5f', delimiter="\t")

    del a, mat, tct, normalized_input
    print("[%s] Done" % (timestamp()))
