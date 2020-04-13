# STARRPeaker
Uniform processing pipeline and peak caller for STARR-seq data

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

## installation
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
2. Filtered alignments using SAMtools (v1.5) with the following arguments
filter: -F 1804 exclude FLAG 1804: unmapped, next segment unmapped, secondary alignments, not passing platform q, PCR or optical duplicates; -f 2 require FLAG 2: properly aligned; -q 30 exclude MAPQ < 30; -u uncompressed output; 
3. Removed duplicates using picard (v2.9.0)
4. Merged biological replicates using SAMtools

## Inputs
* Input alignment (BAM) file (STARR-seq input)
* Output alignment (BAM) file (STARR-seq output)
* Covariates (BigWig) file(s)
* Chrom Size file (i.e., https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes) 
* Blacklist (BED) file (i.e., https://www.encodeproject.org/files/ENCFF419RSJ/)

## Covariates
The peak calling algorithm models STARR-seq fragment coverage across the genome using multiple covariates to correct for potential sequencing bias. It is recommended to include potential confounding variables into the model. These includes but not limited to GC-content, mappability tracks, and so on.

## Usage
```
usage: starrpeaker.py [-h] --prefix PREFIX --chromsize CHROMSIZE --blacklist
                      BLACKLIST -i INPUT -o OUTPUT [--length LENGTH]
                      [--step STEP] [--cov COV [COV ...]] [--min MIN]
                      [--max MAX] [--readstart] [--strand STRAND]
                      [--threshold THRESHOLD] [--mode MODE] [--mincov MINCOV]
                      [--eq EQ]

STARRPeaker

required arguments:
  --prefix PREFIX             Output File Prefix
  --chromsize CHROMSIZE       Chrom Sizes
  --blacklist BLACKLIST       Blacklist Region in BED format
  -i INPUT, --input INPUT     Input BAM File
  -o OUTPUT, --output OUTPUT  STARR-seq BAM File

optional arguments:
  -h, --help                  show this help message and exit
  --length LENGTH             Bin Length (default: 500)
  --step STEP                 Step Size (default: 100)
  --cov COV [COV ...]         Covariate BigWig Files
  --min MIN                   Minimum Template Size (default: 200)
  --max MAX                   Maximum Template Size (default: 1000)
  --readstart                 Use Read Start Position instead of Fragment Center
  --strand STRAND             Use all/fwd/rev Stranded Fragments (default: all)
  --threshold THRESHOLD       Adjusted P-value Threshold (default: 0.05)
  --mode MODE                 Mode [1 - using input as covariate (default), 2 - using input as offset]
  --mincov MINCOV             Minimum Coverage (default: 10)
  --eq EQ                     Extreme Quantile to Remove (default: 1e-5)
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

## Final Peak Call Format (BED6+4)
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
