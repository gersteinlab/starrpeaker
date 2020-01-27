#!/usr/bin/python
from __future__ import division

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2019, Gerstein Lab"
__credits__ = ["Donghoon Lee"]
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
import subprocess
import multiprocessing
from sklearn import preprocessing


# from functools import reduce


def timestamp():
    return str(datetime.datetime.now()).split('.')[0]


def safe_remove(file):
    if os.path.exists(file):
        os.remove(file)


def get_uid():
    return str(uuid.uuid4())[:8]


def make_bin(chromSize, binLength, stepSize, blackList, fileOut):
    ### make sliding window
    print("[%s] Making bins" % (timestamp()))
    bin = pybedtools.BedTool().window_maker(g=chromSize, w=binLength, s=stepSize)

    ### filter blacklist region
    print("[%s] Filtering blacklist region" % (timestamp()))
    blk = pybedtools.BedTool(blackList).sort()
    out = bin.intersect(blk, v=True, sorted=True)

    ### write to file
    with open(fileOut, 'w') as file:
        file.write(str(out))
    del bin, blk, out
    print("[%s] Done" % (timestamp()))


def proc_cov(bwFiles, bedFile, fileOut):
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
    np.savetxt(fileOut, mat, fmt='%.2f', delimiter="\t")
    print("[%s] Done" % (timestamp()))
    del mat


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


def proc_bam(bamFiles, bedFile, chromSize, fileOut, minSize, maxSize):
    '''

    Args:
        bamFiles: list of BAM files eg. [input.bam output.bam]
        bedFile: bin BED file
        chromSize: chrom size file
        fileOut: output file
        minSize: minimum size of fragment insert to consider
        maxSize: maximum size of fragment insert to consider

    Returns:
        writes bin count output file

    '''
    print("[%s] Counting template depth per bin %s" % (timestamp(), bedFile))

    ### initialize numpy array
    tct = np.zeros(shape=(len(bamFiles)), dtype=int)

    ### initialize numpy array
    mat = np.zeros(shape=(sum(1 for l in open(bedFile)), len(bamFiles)), dtype=int)

    ### random unique ID
    uid = get_uid()

    ### load bin bed file
    a = pybedtools.BedTool(bedFile)

    for j, bam in enumerate(bamFiles):
        print("[%s] Processing %s" % (timestamp(), bam))

        if not os.path.exists(bam + ".bai"):
            print("[%s] (Warning) Index not found: %s" % (timestamp(), bam))
            print("[%s] Indexing %s" % (timestamp(), bam))
            pysam.index(bam)

        b = pysam.AlignmentFile(bam, "rb")

        template_count = 0
        template_count_fwd = 0
        template_count_rev = 0
        template_count_used = 0
        template_count_used_fwd = 0
        template_count_used_rev = 0

        ###

        for chr in list_chr(chromSize):

            print("[%s] Processing %s" % (timestamp(), chr))

            with open("tmp" + uid + str(j) + chr + "all.bed", "w") as s, open("tmp" + uid + str(j) + chr + "fwd.bed",
                                                                              "w") as s_fwd, open(
                "tmp" + uid + str(j) + chr + "rev.bed", "w") as s_rev:
                for read in b.fetch(reference=chr):

                    ### read IS first in pair
                    ### read IS "properly paired"
                    ### read is NOT chimeric read (i.e., no SA tag)
                    ### read is mapped to forward strand, mate is mapped to reverse strand
                    if read.is_read1 and read.is_proper_pair and not read.has_tag(
                            "SA"):  ### if read has SA tag (SA: secondary alignment), read is ambiguous and thus discard

                        if not read.is_reverse and read.mate_is_reverse:
                            template_count += 1
                            template_count_fwd += 1
                            if read.template_length >= int(minSize) and read.template_length <= int(maxSize):
                                template_count_used += 1
                                s.write("%s\t%i\t%i\n" % ((b.get_reference_name(read.reference_id)),
                                                          (read.reference_start + int(read.template_length / 2)),
                                                          (read.reference_start + int(read.template_length / 2) + 1)))
                                template_count_used_fwd += 1
                                s_fwd.write("%s\t%i\t%i\n" % ((b.get_reference_name(read.reference_id)),
                                                              (read.reference_start + int(read.template_length / 2)), (
                                                                      read.reference_start + int(
                                                                  read.template_length / 2) + 1)))

                        elif read.is_reverse and not read.mate_is_reverse:
                            template_count += 1
                            template_count_rev += 1
                            if -(read.template_length) >= int(minSize) and -(read.template_length) <= int(maxSize):
                                template_count_used += 1
                                s.write("%s\t%i\t%i\n" % ((b.get_reference_name(read.reference_id)),
                                                          (read.reference_end + int(read.template_length / 2) - 1),
                                                          (read.reference_end + int(read.template_length / 2))))
                                template_count_used_rev += 1
                                s_rev.write("%s\t%i\t%i\n" % ((b.get_reference_name(read.reference_id)),
                                                              (read.reference_end + int(read.template_length / 2) - 1),
                                                              (read.reference_end + int(read.template_length / 2))))

            print("[%s] Sorting %s" % (timestamp(), chr))
            pybedtools.BedTool("tmp" + uid + str(j) + chr + "all.bed").sort().saveas(
                "tmp" + uid + str(j) + chr + "all.sorted.bed")
            pybedtools.BedTool("tmp" + uid + str(j) + chr + "fwd.bed").sort().saveas(
                "tmp" + uid + str(j) + chr + "fwd.sorted.bed")
            pybedtools.BedTool("tmp" + uid + str(j) + chr + "rev.bed").sort().saveas(
                "tmp" + uid + str(j) + chr + "rev.sorted.bed")

            ### delete bed
            safe_remove("tmp" + uid + str(j) + chr + "all.bed")
            safe_remove("tmp" + uid + str(j) + chr + "fwd.bed")
            safe_remove("tmp" + uid + str(j) + chr + "rev.bed")

        ###

        tct[j] += template_count_used

        print("[%s] Total mapped reads: %i" % (timestamp(), b.mapped))
        # print("[%s] %i reads in proper pairs" % (timestamp(), proper_pair_count))
        # print("[%s] %i chimeric reads removed" % (timestamp(), chimeric_count))
        print("[%s] %i templates extracted" % (timestamp(), template_count))
        print("[%s] %i templates extracted (+)" % (timestamp(), template_count_fwd))
        print("[%s] %i templates extracted (-)" % (timestamp(), template_count_rev))

        print("[%s] %i templates used for count" % (timestamp(), template_count_used))
        print("[%s] %i templates used for count (+)" % (timestamp(), template_count_used_fwd))
        print("[%s] %i templates used for count (-)" % (timestamp(), template_count_used_rev))

        b.close()

        ### merge ALL bed
        print("[%s] Merging ALL BED files" % (timestamp()))
        with open("tmp" + uid + str(j) + "all.sorted.merged.bed", "a") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + str(j) + chr + "all.sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + str(j) + chr + "all.sorted.bed")

        print("[%s] Counting depth per bin" % (timestamp()))
        mergedBed = pybedtools.BedTool("tmp" + uid + str(j) + "all.sorted.merged.bed")
        readDepth = a.coverage(mergedBed, sorted=True, counts=True)

        ### extract 4th column, which is read counts, and assign as numpy array
        mat[:, j] = np.array([int(l.split("\t")[3]) for l in str(readDepth).rstrip("\n").split("\n")])

        ### save genome coverage
        print("[%s] Making genome coverage bedGraph" % (timestamp()))
        binSize = int(a[0][2]) - int(a[0][1])
        mergedBed.slop(g=chromSize, b=int(binSize / 2)).genome_coverage(bg=True, g=chromSize).saveas(
            fileOut + "." + str(j) + ".all.bdg")

        ### delete tmp merged bed files
        safe_remove("tmp" + uid + str(j) + "all.sorted.merged.bed")
        del merged, mergedBed, readDepth

        ### convert bedGraph to bigWig
        print("[%s] Converting bedGraph to bigWig" % (timestamp()))
        bdg2bw(bdgFile=fileOut + "." + str(j) + ".all.bdg", bwFile=fileOut + "." + str(j) + ".all.bw",
               chromSize=chromSize)
        safe_remove(fileOut + "." + str(j) + ".all.bdg")

        ### merge FWD bed
        print("[%s] Merging FWD BED files" % (timestamp()))
        with open("tmp" + uid + str(j) + "fwd.sorted.merged.bed", "a") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + str(j) + chr + "fwd.sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + str(j) + chr + "fwd.sorted.bed")

        print("[%s] Counting depth per bin" % (timestamp()))
        mergedBed = pybedtools.BedTool("tmp" + uid + str(j) + "fwd.sorted.merged.bed")

        ### save genome coverage
        print("[%s] Making genome coverage bedGraph" % (timestamp()))
        binSize = int(a[0][2]) - int(a[0][1])
        mergedBed.slop(g=chromSize, b=int(binSize / 2)).genome_coverage(bg=True, g=chromSize).saveas(
            fileOut + "." + str(j) + ".fwd.bdg")

        ### delete tmp merged bed files
        safe_remove("tmp" + uid + str(j) + "fwd.sorted.merged.bed")
        del merged, mergedBed

        ### convert bedGraph to bigWig
        print("[%s] Converting bedGraph to bigWig" % (timestamp()))
        bdg2bw(bdgFile=fileOut + "." + str(j) + ".fwd.bdg", bwFile=fileOut + "." + str(j) + ".fwd.bw",
               chromSize=chromSize)
        safe_remove(fileOut + "." + str(j) + ".fwd.bdg")

        ### merge REV bed
        print("[%s] Merging REV BED files" % (timestamp()))
        with open("tmp" + uid + str(j) + "rev.sorted.merged.bed", "a") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + str(j) + chr + "rev.sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + str(j) + chr + "rev.sorted.bed")

        print("[%s] Counting depth per bin" % (timestamp()))
        mergedBed = pybedtools.BedTool("tmp" + uid + str(j) + "rev.sorted.merged.bed")

        ### save genome coverage
        print("[%s] Making genome coverage bedGraph" % (timestamp()))
        binSize = int(a[0][2]) - int(a[0][1])
        mergedBed.slop(g=chromSize, b=int(binSize / 2)).genome_coverage(bg=True, g=chromSize).saveas(
            fileOut + "." + str(j) + ".rev.bdg")

        ### delete tmp merged bed files
        safe_remove("tmp" + uid + str(j) + "rev.sorted.merged.bed")
        del merged, mergedBed

        ### convert bedGraph to bigWig
        print("[%s] Converting bedGraph to bigWig" % (timestamp()))
        bdg2bw(bdgFile=fileOut + "." + str(j) + ".rev.bdg", bwFile=fileOut + "." + str(j) + ".rev.bw",
               chromSize=chromSize)
        safe_remove(fileOut + "." + str(j) + ".rev.bdg")

    ### normalize input count, normalized input count is added to additional column
    normalized_input = mat[:, 0] * (tct[1] / tct[0])
    np.savetxt(fileOut, np.concatenate((mat, normalized_input.reshape(-1, 1)), axis=1), fmt='%.5f', delimiter="\t")

    del a, mat, tct, normalized_input
    print("[%s] Done" % (timestamp()))


def proc_bam_readstart(bamFiles, bedFile, chromSize, fileOut):
    print("[%s] Counting template depth per bin %s" % (timestamp(), bedFile))

    ### initialize numpy array
    tct = np.zeros(shape=(len(bamFiles)), dtype=int)

    ### initialize numpy array
    mat = np.zeros(shape=(sum(1 for l in open(bedFile)), len(bamFiles)), dtype=int)

    ### random unique ID
    uid = get_uid()

    ### load bin bed file
    a = pybedtools.BedTool(bedFile)

    for j, bam in enumerate(bamFiles):
        print("[%s] Processing %s" % (timestamp(), bam))

        if not os.path.exists(bam + ".bai"):
            print("[%s] (Warning) Index not found: %s" % (timestamp(), bam))
            print("[%s] Indexing %s" % (timestamp(), bam))
            pysam.index(bam)

        b = pysam.AlignmentFile(bam, "rb")

        proper_pair_count = 0
        chimeric_count = 0
        template_count = 0
        proper_template_count = 0

        for chr in list_chr(chromSize):

            print("[%s] Processing %s" % (timestamp(), chr))

            with open("tmp" + uid + str(j) + chr + ".bed", "w") as s:
                for read in b.fetch(reference=chr):
                    s.write("%s\t%i\t%i\n" % ((b.get_reference_name(read.reference_id)), (read.reference_start),
                                              (read.reference_start + 1)))  ## start position of read

            print("[%s] Sorting %s" % (timestamp(), chr))
            pybedtools.BedTool("tmp" + uid + str(j) + chr + ".bed").sort().saveas(
                "tmp" + uid + str(j) + chr + "sorted.bed")

            ### delete bed
            safe_remove("tmp" + uid + str(j) + chr + ".bed")

        tct[j] += proper_template_count

        print("[%s] Total mapped reads: %i" % (timestamp(), b.mapped))
        print("[%s] %i reads in proper pairs" % (timestamp(), proper_pair_count))
        print("[%s] %i chimeric reads removed" % (timestamp(), chimeric_count))
        print("[%s] %i templates extracted" % (timestamp(), template_count))
        print("[%s] %i templates used for count" % (timestamp(), proper_template_count))

        b.close()

        ### merge bed
        print("[%s] Merging BED files" % (timestamp()))
        with open("tmp" + uid + str(j) + ".merged.bed", "a") as merged:
            for chr in list_chr(chromSize):

                ### merge tmp bed files
                with open("tmp" + uid + str(j) + chr + "sorted.bed", "r") as t:
                    if t.read(1).strip():
                        t.seek(0)
                        merged.write(t.read())

                ### delete tmp bed files
                safe_remove("tmp" + uid + str(j) + chr + "sorted.bed")

        print("[%s] Counting depth per bin" % (timestamp()))
        mergedBed = pybedtools.BedTool("tmp" + uid + str(j) + ".merged.bed")
        readDepth = a.coverage(mergedBed, sorted=True, counts=True)

        ### extract 4th column, which is read counts, and assign as numpy array
        mat[:, j] = np.array([int(l.split("\t")[3]) for l in str(readDepth).rstrip("\n").split("\n")])

        ### save genome coverage
        print("[%s] Making genome coverage bedGraph" % (timestamp()))
        binSize = int(a[0][2]) - int(a[0][1])
        mergedBed.slop(g=chromSize, b=int(binSize / 2)).genome_coverage(bg=True, g=chromSize).saveas(
            fileOut + "." + str(j) + ".bdg")
        bdg2bw(fileOut + "." + str(j) + ".bdg", fileOut + "." + str(j) + ".bw", chromSize)

        ### delete tmp merged bed files
        safe_remove("tmp" + uid + str(j) + ".merged.bed")
        safe_remove(fileOut + "." + str(j) + ".bdg")
        del merged, mergedBed, readDepth

    ### normalize input count
    normalized_input = mat[:, 0] * (tct[1] / tct[0])
    np.savetxt(fileOut, np.concatenate((mat, normalized_input.reshape(-1, 1)), axis=1), fmt='%.5f', delimiter="\t")

    print("[%s] Done" % (timestamp()))
    del a, mat, tct


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


def call_peak(prefix, bedFile, bctFile, covFile, bwFile, chromSize, threshold, mode):
    '''

    Args:
        prefix: prefix for the output
        bedFile: bin in BED format
        bctFile: fragment insert coverage in BCT format (made from procBam.py)
        covFile: covariates in TSV format
        bwFile: fragment coverage in BigWig format
        chromSize: chromosome sizes
        threshold: threshold to call peak
        mode: 1 - using input as covariate 2 - using input as offset

    Returns:
        peak files (original peaks, merged peaks, and centered final peaks)

    '''
    print("[%s] Calling peaks" % (timestamp()))

    ### load data
    print("[%s] Loading fragment coverages, and covariates" % (timestamp()))
    bct = np.loadtxt(bctFile, ndmin=2)  # 0=input, 1=output, 2=normalized input
    cov = np.loadtxt(covFile, ndmin=2)  # 3:=cov

    ### scale covariates to have mean 0 and sd 1
    cov_scaled = preprocessing.scale(cov, axis=0)

    ### merge data
    mat = np.concatenate((bct[:, [1, 0, 2]], cov_scaled), axis=1)  # 0=output, 1=input, 2=normalized input, 3:=cov
    del bct, cov, cov_scaled

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

    ### remove bins with input (and output) count of zero (i.e., untested region) OR extreme values (top 1%, i.e., sequencing artifacts)
    minInput = 0
    maxInput = np.quantile(mat[:, 1], 0.99)
    nonZeroInput = (mat[:, 1] > minInput) & (mat[:, 1] < maxInput)
    nonZeroOutput = (mat[:, 0] > 0)

    ### calculate fold change
    fc = np.zeros(mat.shape[0])
    fc[mat[:, 1] > 0] = mat[mat[:, 1] > 0, 0] / (mat[mat[:, 1] > 0, 2])

    ### filtering bins
    print("[%s] Before filtering: %s" % (timestamp(), mat.shape[0]))

    print("[%s] Removing %i bins with insufficient input coverage" % (timestamp(), sum(np.invert(nonZeroInput))))
    print("[%s] Bins with sufficient input coverage: %s" % (timestamp(), mat[nonZeroInput, :].shape[0]))

    print("[%s] Removing %i sliding bins" % (timestamp(), sum(np.invert(nonSliding))))
    print("[%s] Bins with non-sliding window: %s" % (timestamp(), mat[nonSliding, :].shape[0]))

    print("[%s] After filtering: %s" % (timestamp(), mat[nonZeroInput & nonSliding, :].shape[0]))

    ### mode 2 uses "input" as offset variable
    if int(mode) == 2:
        print("[%s] Running Mode 2" % (timestamp()))
        ### remove input
        mat = np.delete(mat, 1, 1)

        ### formula
        x = ["x" + str(i) for i in range(1, mat.shape[1] - 1)]
        df = pd.DataFrame(mat[nonZeroInput & nonSliding, :], columns=["y", "exposure"] + x)
        formula = "y~" + "+".join(df.columns.difference(["y", "exposure"]))
        print("[%s] Fit using formula: %s" % (timestamp(), formula))

        ### Initial parameter estimation using Poisson regression
        # print("[%s] Initial estimate" % (timestamp()))
        model0 = smf.glm(formula, data=df, family=sm.families.Poisson(), offset=np.log(df["exposure"])).fit()
        # print model0.summary()

        ### Estimate theta
        th0, _ = theta(mat[nonZeroInput & nonSliding, :][:, 0], model0.mu)
        print("[%s] Initial estimate of theta is %f" % (timestamp(), th0))

        ### re-estimate beta with theta
        # print("[%s] Re-estimate of beta" % (timestamp()))
        model = smf.glm(formula, data=df, family=sm.families.NegativeBinomial(alpha=1 / th0),
                        offset=np.log(df["exposure"])).fit(start_params=model0.params)
        # print model.summary()

        ### Re-estimate theta
        th, _ = theta(mat[nonZeroInput & nonSliding, :][:, 0], model.mu)
        print("[%s] Re-estimate of theta is %f" % (timestamp(), th))

        ### predict
        print("[%s] Predicting expected counts for bins above a minimum threshold: %s" % (
            timestamp(), mat[nonZeroInput & nonZeroOutput, :].shape[0]))
        df = pd.DataFrame(mat[nonZeroInput & nonZeroOutput, :], columns=["y", "exposure"] + x)
        y_hat = model.predict(df, offset=np.log(df["exposure"]))

    ### mode 1 uses "input" as covariate (default):
    else:
        print("[%s] Running Mode 1" % (timestamp()))
        ### remove normalized input
        mat = np.delete(mat, 2, 1)

        ### formula
        x = ["x" + str(i) for i in range(1, mat.shape[1])]
        df = pd.DataFrame(mat[nonZeroInput & nonSliding, :], columns=["y"] + x)
        formula = "y~" + "+".join(df.columns.difference(["y"]))
        print("[%s] Fit using formula: %s" % (timestamp(), formula))

        ### Initial parameter estimation using Poisson regression
        # print("[%s] Initial estimate" % (timestamp()))
        model0 = smf.glm(formula, data=df, family=sm.families.Poisson()).fit()
        # print model0.summary()

        ### Estimate theta
        th0, _ = theta(mat[nonZeroInput & nonSliding, :][:, 0], model0.mu)
        print("[%s] Initial estimate of theta is %f" % (timestamp(), th0))

        ### re-estimate beta with theta
        # print("[%s] Re-estimate of beta" % (timestamp()))
        model = smf.glm(formula, data=df, family=sm.families.NegativeBinomial(alpha=1 / th0)).fit(
            start_params=model0.params)
        # print model.summary()

        ### Re-estimate theta
        th, _ = theta(mat[nonZeroInput & nonSliding, :][:, 0], model.mu)
        print("[%s] Re-estimate of theta is %f" % (timestamp(), th))

        ### predict
        print("[%s] Predicting expected counts for bins above a minimum threshold: %s" % (
            timestamp(), mat[nonZeroInput & nonZeroOutput, :].shape[0]))
        df = pd.DataFrame(mat[nonZeroInput & nonZeroOutput, :], columns=["y"] + x)
        y_hat = model.predict(df)

    ### calculate P-value
    print("[%s] Calculating P-value" % (timestamp()))
    theta_hat = np.repeat(th, len(y_hat))
    prob = th / (th + y_hat)  ### prob=theta/(theta+mu)
    pval = 1 - nbinom.cdf(mat[nonZeroInput & nonZeroOutput, 0] - 1, n=theta_hat, p=prob)
    del mat

    ### multiple testing correction
    print("[%s] Multiple testing correction" % (timestamp()))
    _, pval_adj, _, _ = multi.multipletests(pval, method="fdr_bh")

    p_score = -np.log10(pval)
    q_score = -np.log10(pval_adj)

    ### output peak
    print("[%s] Output initial peaks" % (timestamp()))
    with open(prefix + ".peak.bed", "w") as out:
        with open(bedFile, "r") as bed:
            for i, bin in enumerate(list(compress(bed.readlines(), nonZeroInput & nonZeroOutput))):
                if pval[i] <= float(threshold):
                    out.write("%s\t%.3f\t%.3f\t%.3f\t%.5e\t%.5e\n" % (
                        bin.strip(), fc[nonZeroInput & nonZeroOutput][i], p_score[i], q_score[i], pval[i], pval_adj[i]))

    ### output p-val track
    print("[%s] Generating P-value bedGraph" % (timestamp()))
    with open(prefix + ".pval.bdg", "w") as out:
        with open(bedFile, "r") as bed:
            for i, bin in enumerate(list(compress(bed.readlines(), nonZeroInput & nonZeroOutput))):
                out.write("%s\t%.3f\n" % (bin.strip(), abs(p_score[i])))
    del p_score, q_score, pval, pval_adj

    ### make bigWig track
    print("[%s] Making BigWig tracks" % (timestamp()))
    make_bigwig(prefix=prefix, bedFile=bedFile, bctFile=bctFile, chromSize=chromSize, bedGraphFile=prefix + ".pval.bdg")
    safe_remove(prefix + ".pval.bdg")

    ### merge peak
    print("[%s] Merge peaks" % (timestamp()))
    pybedtools.BedTool(prefix + ".peak.bed").merge(c=[4, 5, 6, 7, 8], o=["max", "max", "max", "min", "min"]).saveas(
        prefix + ".peak.merged.bed")

    ### center merged peak
    print("[%s] Finalizing peaks" % (timestamp()))
    center_peak(bwFile=bwFile,
                peakFile=prefix + ".peak.merged.bed",
                centeredPeakFile=prefix + ".peak.final.bed")

    print("[%s] Done" % (timestamp()))


def make_bigwig(prefix, bedFile, bctFile, chromSize, bedGraphFile=""):
    with open(chromSize) as f: cs = [line.strip().split('\t') for line in f.readlines()]

    bin = np.genfromtxt(bedFile, dtype=str)
    bct = np.loadtxt(bctFile, dtype=np.float32, ndmin=2)

    starts = np.array(bin[:, 1], dtype=np.int32)
    ends = np.array(bin[:, 2], dtype=np.int32)

    l = ends[0] - starts[0]
    s = starts[1] - starts[0]
    print("[%s] Using fixed interval of %i" % (timestamp(), s))

    nonoverlapping = ends - starts == l

    chroms = (np.array(bin[:, 0]))[nonoverlapping]
    starts = (np.array(bin[:, 1], dtype=np.int32) + int(l / 2) - int(s / 2))[nonoverlapping]
    ends = (np.array(bin[:, 2], dtype=np.int32) - int(l / 2) + int(s / 2))[nonoverlapping]
    val_input = np.array(bct[:, 0], dtype=np.float32)[nonoverlapping]
    val_output = np.array(bct[:, 1], dtype=np.float32)[nonoverlapping]
    val_normalized_input = np.array(bct[:, 2], dtype=np.float32)[nonoverlapping]

    chroms_fc = chroms[np.nonzero(val_normalized_input)]
    starts_fc = starts[np.nonzero(val_normalized_input)]
    ends_fc = ends[np.nonzero(val_normalized_input)]
    val_fc = val_output[np.nonzero(val_normalized_input)] / val_normalized_input[np.nonzero(val_normalized_input)]

    ### input signal

    bw0 = pyBigWig.open(prefix + ".input.bw", "w")
    bw0.addHeader([(str(x[0]), int(x[1])) for x in cs])
    bw0.addEntries(chroms=chroms, starts=starts, ends=ends, values=val_input)
    bw0.close()

    ### output signal

    bw1 = pyBigWig.open(prefix + ".output.bw", "w")
    bw1.addHeader([(str(x[0]), int(x[1])) for x in cs])
    bw1.addEntries(chroms=chroms, starts=starts, ends=ends, values=val_output)
    bw1.close()

    ### normalized input signal

    bw2 = pyBigWig.open(prefix + ".normalized_input.bw", "w")
    bw2.addHeader([(str(x[0]), int(x[1])) for x in cs])
    bw2.addEntries(chroms=chroms, starts=starts, ends=ends, values=val_normalized_input)
    bw2.close()

    ### fold change

    bw3 = pyBigWig.open(prefix + ".fc.bw", "w")
    bw3.addHeader([(str(x[0]), int(x[1])) for x in cs])
    bw3.addEntries(chroms=chroms_fc, starts=starts_fc, ends=ends_fc, values=val_fc)
    bw3.close()

    # with open(prefix+".bedGraph","w") as b:
    #     b.write("track type=bedGraph\n")
    #     for x in zip(chroms,starts,ends,val_input):
    #         b.write('\t'.join(map(str,x))+'\n')

    if bedGraphFile != "":
        print("[%s] Making P-value BigWig track" % (timestamp()))

        bin = np.genfromtxt(bedGraphFile, dtype=str)
        starts = np.array(bin[:, 1], dtype=np.int32)
        ends = np.array(bin[:, 2], dtype=np.int32)

        nonoverlapping = ends - starts == l

        chroms = (np.array(bin[:, 0]))[nonoverlapping]
        starts = (np.array(bin[:, 1], dtype=np.int32) + int(l / 2) - int(s / 2))[nonoverlapping]
        ends = (np.array(bin[:, 2], dtype=np.int32) - int(l / 2) + int(s / 2))[nonoverlapping]
        val_pval = np.array(bin[:, 3], dtype=np.float32)[nonoverlapping]

        ### pval signal

        bw4 = pyBigWig.open(prefix + ".pval.bw", "w")
        bw4.addHeader([(str(x[0]), int(x[1])) for x in cs])
        bw4.addEntries(chroms=chroms, starts=starts, ends=ends, values=val_pval)
        bw4.close()


def bdg2bw(bdgFile, bwFile, chromSize):
    with open(chromSize) as f:
        cs = [line.strip().split('\t') for line in f.readlines()]

    bw = pyBigWig.open(bwFile, "w")
    bw.addHeader([(str(x[0]), int(x[1])) for x in cs])

    with open(bdgFile, "r") as bdg:
        for line in bdg:
            if len(line.strip().split("\t")) == 4:
                chr, start, end, val = line.strip().split("\t")
                bw.addEntries(chroms=[chr], starts=[int(start)], ends=[int(end)], values=[float(val)])
            else:
                print("[%s] Warning: skipping bedGraph entry: %s" % (timestamp(), line.strip()))

    bw.close()


def center_peak(bwFile, peakFile, centeredPeakFile):
    bw = pyBigWig.open(bwFile)
    peak = pybedtools.BedTool(peakFile)

    with open(centeredPeakFile, "w") as out:
        for p in peak:
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
                out.write('\t'.join(
                    [chr, str(center - radius), str(center + radius), other]) + '\n')


def proc_fenergy(bedFile, fileOut, linearfold, genome):
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
    print("[%s] Parallel computing using %i cores" % (timestamp(), multiprocessing.cpu_count()))
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(run_linearfold, prefixes)
    pool.terminate()

    ### merge bed
    print("[%s] Merging files" % (timestamp()))
    safe_remove(fileOut)
    with open(fileOut, "a") as merged:
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
        out = subprocess.check_output([arg[1], "-V"], stdin=fa)
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
