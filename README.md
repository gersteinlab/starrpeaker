# STARRPeaker
Uniform processing pipeline and peak caller for STARR-seq data

## Dependencies
* Python 2.7 (tested on Python 2.7.10-2.7.15)
* pysam
* pybedtools
* pyBigWig
* numpy
* scipy
* pandas
* statsmodels
* sklearn (tested on 0.20.3)

## installation
Preferably, create a conda environment with Python 2.7
```
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
* *prefix*.bin.bed: bin BED file
* *prefix*.bam.bct: bin count in BST format (1st col: input, 2nd col: output, 3rd col: normalized input)
* *prefix*.cov.tsv: covariate matrix in TSV format
* *prefix*.input.bw: input fragment coverage in bigWig format
* *prefix*.output.bw: output fragment coverage in bigWig format
* *prefix*.peak.bed: initial peak calls (before centering and merging)
* *prefix*.peak.final.bed: final peak calls
* *prefix*.peak.pval.bw: p-value track in bigWig format (-log10)
* *prefix*.peak.qval.bw: q-value track in bigWig format (-log10)

## Peak File Format
* Column 1: chromosome
* Column 2: start
* Column 3: end
* Column 4: fold change (output/normalized-input)
* Column 5: input fragment count
* Column 6: output fragment count
* Column 7: p-value
* Column 8: q-value (Benjamini-Hochberg False Discovery Rate, FDR)
