# starrpeaker
Peak caller for STARR-seq data

## Dependencies

* Python 2.7 (tested on Python 2.7.10)
* pysam
* pybedtools
* pyBigWig
* numpy
* scipy
* pandas
* statsmodels

## installation

```
git clone https://github.com/hoondy/starrpeaker.git
cd starrpeaker
sudo python setup.py install --record files.txt
starrpeaker -h
```

## Input

* Input alignment (BAM) file
* Output alignment (BAM) file
* Covariates (BigWig) file(s)
* Chrom Size file (i.e., https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes) 
* Blacklist (BED) file (i.e., https://www.encodeproject.org/files/ENCFF419RSJ/)

## Covariates

The peak calling algorithm models STARR-seq read counts across the genome using multiple covariates to correct for potential sequencing bias. It is recommended to include potential confounding factors into the model. These includes but not limited to mappability, GC-content, cross-species conservation, and so on.

## Usage

```
usage: starrpeaker [-h] -p PREFIX --chromsize CHROMSIZE [-l LENGTH] [-s STEP]
                   -x BLACKLIST --cov COV [COV ...] -i INPUT [INPUT ...] -o
                   OUTPUT [OUTPUT ...] [--min MIN] [--max MAX] [-t THRESHOLD]

STARR-Peaker

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Output File Prefix
  --chromsize CHROMSIZE
                        Chrom Sizes
  -l LENGTH, --length LENGTH
                        Bin Length
  -s STEP, --step STEP  Step Size
  -x BLACKLIST, --blacklist BLACKLIST
                        Blacklist Region in BED format
  --cov COV [COV ...]   Covariate BigWig Files
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Input BAM File
  -o OUTPUT [OUTPUT ...], --output OUTPUT [OUTPUT ...]
                        STARR-seq BAM File
  --min MIN             Minimum Template Size
  --max MAX             Maximum Template Size
  -t THRESHOLD, --threshold THRESHOLD
                        Adjusted P-value Threshold
```

* '-p' or '--prefix': Output File Prefix
* '--chromsize': Chrom Sizes
* '-l' or '--length': Bin Length, default=500
* '-s' or '--step': Step Size, default=100
* '-x' or '--blacklist': Blacklist Region in BED format
* '--cov': Covariate BigWig Files

* '-i' or '--input': Input BAM File
* '-o' or '--output': STARR-seq BAM File

optional args
* '--min': Minimum Template Size, default=100
* '--max': Maximum Template Size, default=1000
* '-t', '--threshold': Adjusted P-value Threshold, default=0.05

```
starrpeaker --prefix <prefix for output files> --chromsize <hg38.chrom.sizes> --blacklist <blacklist_GRCh38.bed> --cov <covariate 1: gc content> <covariate 2: mappability> <covariate 3: conservation> --input <input.bam> --output <output.bam>
```

## Output

* <prefix>.input.tsv: input count
* <prefix>.output.tsv: output count
* <prefix>.cov.tsv: covariate matrix
* <prefix>.bin.bed: bin
* <prefix>.peak.bed: peaks per bin
* <prefix>.peak.merged.bed: merged peaks
