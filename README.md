# STARRPeaker
Uniform processing pipeline and peak caller for STARR-seq data

## Changelog
### v1.2
- Added a new functionality to restrict analysis to a supplied capture panel BED file (--capture)

### v1.1
- Updated the final peak call (BED6+) ENCODE specification (https://www.encodeproject.org/documents/9f5d2b5a-bd29-4983-9c01-fab4ab8b5ea2/)
- In specific, fold change is changed to log2 fold change
- In specific, input coverage is reported along with output coverage

### v1.0
- Updated documentation for ENCODE release (https://bit.ly/whg-starr-seq)

### v1.0-rc
- Release candidate with early version of documentation

## Dependencies (version tested)
* Python 2.7 (v2.7.15)
* pysam (v0.15.3)
* pybedtools (v0.8.1)
* pyBigWig (v0.3.13)
* numpy (v1.15.4)
* scipy (v1.2.0)
* pandas (v0.24.1)
* statsmodels (v0.10.1, use v0.10.2 or earlier, new function statsmodels/tools/validation/validation.py introduced in v0.11.0 may introduce error in Python 2)
* scikit-learn (v0.20.3)

## Installation
Preferably, create a conda environment with Python 2.7
```
conda create -n starrpeaker python=2.7 pybedtools
conda activate starrpeaker
pip install git+https://github.com/gersteinlab/starrpeaker
starrpeaker -h
```

## Preprocessing
*Few notes on how alignment (BAM) files were prepared*

For each biological replicates in FASTQ format

1. Aligned paired-end reads using BWA mem (v0.7.17)
2. Removed duplicates using picard (v2.9.0)
3. Filtered alignments using SAMtools (v1.9) with the following arguments
```
samtools view -F 3852 -f 2 -q 40

# -F: exclude FLAG 3852
#    4 read unmapped (0x4)
#    8 mate unmapped (0x8)
#  256 not primary alignment (0x100)
#  512 read fails platform/vendor quality checks (0x200)
# 1024 read is PCR or optical duplicate (0x400)
# 2048 supplementary alignment (0x800)

# -f: require FLAG 2
#    2 properly aligned

# -q: exclude MAPQ less than 40
```
4. Merged biological replicates using SAMtools

## Inputs
* Input alignment (BAM) file (STARR-seq input)
* Output alignment (BAM) file (STARR-seq output)
* Covariates (BigWig) file(s)
* Chrom Size file (i.e., https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes) 
* Blacklist (BED) file (i.e., https://www.encodeproject.org/files/ENCFF419RSJ/)

## Covariates
The peak calling algorithm models STARR-seq fragment coverage across the genome using multiple covariates to correct for potential sequencing bias. It is recommended to include potential confounding variables into the model. These include but not limited to GC-content, mappability tracks, and so on.

The following covariates have been precomputed for use with STARRPeaker:
* [Download link to covariates](http://gofile.me/4kBY9/EKNVQrVHM)

Sources:
* GRCh38 genome: https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/
* GRCh38 GC content: https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw
* GRCh38 mappability: (computed using gem-library)
* GRCh38 RNA folding energy: (computed using linearfold; see MS for details)
* hg19 genome: https://www.encodeproject.org/files/male.hg19/
* hg19 GC content: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/gc5Base/hg19.gc5Base.txt.gz 
* hg19 mappability: (computed using gem-library)
* hg19 RNA folding energy: (computed using linearfold; see MS for details)

## Usage
```
usage: starrpeaker.py [-h] --prefix PREFIX --chromsize CHROMSIZE --blacklist
                      BLACKLIST -i INPUT -o OUTPUT [--length LENGTH]
                      [--step STEP] [--cov COV [COV ...]] [--min MIN]
                      [--max MAX] [--readstart] [--strand STRAND]
                      [--threshold THRESHOLD] [--mode MODE] [--mincov MINCOV]
                      [--eq EQ] [--se] [--minfc MINFC] [--capture CAPTURE]
                      [--slop SLOP]

STARRPeaker

optional arguments:
  -h, --help            show this help message and exit
  --prefix PREFIX       Output File Prefix
  --chromsize CHROMSIZE
                        Chrom Sizes
  --blacklist BLACKLIST
                        Blacklist Region in BED format
  -i INPUT, --input INPUT
                        Input BAM File
  -o OUTPUT, --output OUTPUT
                        STARR-seq BAM File
  --length LENGTH       Bin Length
  --step STEP           Step Size
  --cov COV [COV ...]   Covariate BigWig Files
  --min MIN             Minimum Template Size
  --max MAX             Maximum Template Size
  --readstart           Use Read Start Position instead of Fragment Center
  --strand STRAND       Use all/fwd/rev Stranded Fragments
  --threshold THRESHOLD
                        Adjusted P-value Threshold
  --mode MODE           Mode
  --mincov MINCOV       Minimum Coverage
  --eq EQ               Extreme Quantile to Remove
  --se                  Use Single-End instead of Paired-end Sequencing
  --minfc MINFC         Minumum Fold Change
  --capture CAPTURE     Capture Region in BED format
  --slop SLOP           Extend Capture Region in each direction
```

## Example
```
starrpeaker --prefix <prefix for output files> --chromsize <hg38.chrom.sizes> --blacklist <blacklist_GRCh38.bed> --cov <covariate 1: gc content> <covariate 2: mappability> <covariate 3: conservation> --input <input.bam> --output <output.bam> --threshold 0.05
```

## Outputs Files
* *prefix*.bin.bed: Genomic bin BED file
* *prefix*.bam.bct: Alignment counts in BST format (1st col: input, 2nd col: output, 3rd col: normalized input)
* *prefix*.cov.tsv: Covariate matrix in TSV format
* *prefix*.input.bw: Input fragment coverage in bigWig format
* *prefix*.output.bw: Output fragment coverage in bigWig format
* *prefix*.peak.bed: Initial peak calls (before centering and merging)
* *prefix*.peak.final.bed: Final peak calls
* *prefix*.peak.pval.bw: P-value track in bigWig format (-log10)
* *prefix*.peak.qval.bw: Q-value track in bigWig format (-log10)

## Final Peak Call Format (v1.1 and above; BED6+5)
* Column 1: Chromosome
* Column 2: Start position
* Column 3: End position
* Column 4: Name (peak rank based on score, 1 being the highest rank)
* Column 5: Score (integer value of "100 * fold change", maxed at 1000 per BED format specification)
* Column 6: Strand
* Column 7: Log2 Fold change (normalized output/input ratio, in log2 space)
* Column 8: Input fragment coverage (total fragments across/within replicate(s))
* Column 9: Output fragment coverage (total fragments across/within replicate(s))
* Column 10: -log10 of P-value
* Column 11: -log10 of Q-value (Benjamini-Hochberg False Discovery Rate, FDR)

*ENCODE MPRA/STARR-seq BED6+5 common file format: https://www.encodeproject.org/documents/9f5d2b5a-bd29-4983-9c01-fab4ab8b5ea2/*

## Final Peak Call Format (up to v1.0; BED6+4)
* Column 1: Chromosome
* Column 2: Start position
* Column 3: End position
* Column 4: Name (peak rank based on score, 1 being the highest rank)
* Column 5: Score (integer value of "100 * fold change", maxed at 1000 per BED format specification)
* Column 6: Strand
* Column 7: Fold change (output/normalized-input)
* Column 8: Output fragment coverage
* Column 9: -log10 of P-value
* Column 10: -log10 of Q-value (Benjamini-Hochberg False Discovery Rate, FDR)

*BED format specification: https://genome.ucsc.edu/FAQ/FAQformat.html#format1*
