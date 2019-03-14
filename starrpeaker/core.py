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
                    if read.template_length >= minSize and read.template_length <= maxSize:
                        proper_template_count += 1
    b.close()
    return proper_template_count


def proc_bam(bamFiles, bedFile, chromSize, fileOut, minSize, maxSize):
    print("[%s] Counting template depth per bin %s" % (timestamp(), bedFile))

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
                                if read.template_length >= minSize and read.template_length <= maxSize:
                                    proper_template_count += 1
                                    s.write("%s\t%i\t%i\n" % ((b.get_reference_name(read.reference_id)),
                                                              (read.reference_start + int(read.template_length / 2)), (
                                                                  read.reference_start + int(
                                                                      read.template_length / 2) + 1)))

            print("[%s] Sorting %s" % (timestamp(), chr))
            pybedtools.BedTool("tmp" + uid + str(j) + chr + ".bed").sort().saveas(
                "tmp" + uid + str(j) + chr + "sorted.bed")

            ### delete bed
            safe_remove("tmp" + uid + str(j) + chr + ".bed")

        print("[%s] Total mapped reads: %i" % (timestamp(), b.mapped))
        print("[%s] %i proper pairs" % (timestamp(), proper_pair_count))
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
        readDepth = a.coverage(pybedtools.BedTool("tmp" + uid + str(j) + ".merged.bed"), sorted=True, counts=True)

        ### extract 4th column, which is read counts, and assign as numpy array
        mat[:, j] = np.array([int(l.split("\t")[3]) for l in str(readDepth).rstrip("\n").split("\n")])

        ### delete tmp merged bed files
        safe_remove("tmp" + uid + str(j) + ".merged.bed")

    np.savetxt(fileOut, mat, fmt='%i', delimiter="\t")
    print("[%s] Done" % (timestamp()))
    del a, mat


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


def call_peak(inputFile, outputFile, covFile, bedFile, fileOut, threshold, totalInput, totalOutput):
    print("[%s] Calling peaks" % (timestamp()))

    ### load data
    print("[%s] Loading response, exposure, and covariates" % (timestamp()))
    input = np.loadtxt(inputFile).reshape((-1, 1))
    output = np.loadtxt(outputFile).reshape((-1, 1))
    cov = np.loadtxt(covFile)

    ### normalize input count
    if totalOutput != 0 and totalInput != 0:
        normalized_input = input * (totalOutput / totalInput)
    else:
        normalized_input = input * (sum(output) / sum(input))

    ### merge data
    mat = np.concatenate((output, normalized_input, cov), axis=1)
    del input, normalized_input, output, cov

    ### remove bins with zero input count (i.e., untested region)
    nonZeroInput = mat[:, 1] != 0

    ### non sliding bins
    nonSliding = np.zeros(mat.shape[0], dtype=bool) ### initialize with False
    with open(bedFile, "r") as bed:
        lastchr,lastbin="",0
        for i, bin in enumerate(bed.readlines()):
            if bin.split("\t")[0] != lastchr:
                lastchr=bin.split("\t")[0]
                lastbin=int(bin.split("\t")[2])
                nonSliding[i]=True
            elif int(bin.split("\t")[1]) >= lastbin:
                lastbin=int(bin.split("\t")[2])
                nonSliding[i]=True

    ### filter inputs with zero
    print("[%s] Removing %i bins with zero input" % (timestamp(), sum(np.invert(nonZeroInput))))
    print("[%s] Before filtering: %s" % (timestamp(), mat.shape))
    print("[%s] After filtering: %s" % (timestamp(), mat[nonZeroInput, :].shape))
    print("[%s] Non sliding bin: %s" % (timestamp(), mat[nonSliding, :].shape))
    print("[%s] Non zero input & Non sliding bin: %s" % (timestamp(), mat[nonZeroInput & nonSliding, :].shape))

    ### formula
    x = ["x" + str(i) for i in range(1, mat.shape[1] - 1)]
    df = pd.DataFrame(mat[nonZeroInput & nonSliding, :], columns=["y", "exposure"] + x)
    formula = "y~" + "+".join(df.columns.difference(["y", "exposure"]))
    print("[%s] Fit using formula: %s" % (timestamp(), formula))

    ### Initial parameter estimation using Poisson regression
    print("[%s] Initial estimate" % (timestamp()))
    model0 = smf.glm(formula, data=df, family=sm.families.Poisson(), offset=np.log(df["exposure"])).fit()
    # print model0.summary()

    ### Estimate theta
    th0, _ = theta(mat[nonZeroInput & nonSliding,:][:,0], model0.mu)
    print("[%s] Initial estimate of theta is %f" % (timestamp(), th0))

    ### re-estimate beta with theta
    print("[%s] Re-estimate of beta" % (timestamp()))
    model = smf.glm(formula, data=df, family=sm.families.NegativeBinomial(alpha=1 / th0),
                    offset=np.log(df["exposure"])).fit(start_params=model0.params)
    # print model.summary()

    ### Re-estimate theta
    th, _ = theta(mat[nonZeroInput & nonSliding,:][:,0], model.mu)
    print("[%s] Re-estimate of theta is %f" % (timestamp(), th))

    ### predict
    df = pd.DataFrame(mat[nonZeroInput, :], columns=["y", "exposure"] + x)
    y_hat = model.predict(df, offset=np.log(df["exposure"]))

    ### calculate P-value
    print("[%s] Calculating P-value" % (timestamp()))
    theta_hat = np.repeat(th, len(y_hat))
    prob = th / (th + y_hat)  ### prob=theta/(theta+mu)
    pval = 1 - nbinom.cdf(mat[nonZeroInput, 0] - 1, n=theta_hat, p=prob)
    del mat

    ### multiple testing correction
    _, pval_adj, _, _ = multi.multipletests(pval, method="fdr_bh")

    p_score = -np.log10(pval)
    q_score = -np.log10(pval_adj)

    ### output
    with open(fileOut, "w") as out:
        with open(bedFile, "r") as bed:
            for i, bin in enumerate(list(compress(bed.readlines(), nonZeroInput))):
                if pval_adj[i] <= float(threshold):
                    out.write("%s\t%.3f\t%.3f\t%.5e\t%.5e\n" % (
                        bin.strip(), p_score[i], q_score[i], pval[i], pval_adj[i]))
    pybedtools.BedTool(fileOut).merge(c=[4,5,6,7],o=["max","max","min","min"]).saveas(fileOut.replace(".bed","")+".merged.bed")

    print("[%s] Done" % (timestamp()))

