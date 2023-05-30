
## Contents

- [Contents](#contents)
  - [SRA Toolkit 3.0.5](#sra-toolkit-305)
- [Fastq info](#fastq-info)
- [MultiQC](#multiqc)
  - [Trim Galore  0.6.5](#trim-galore--065)
- [Trimmomatic](#trimmomatic)
  - [FASTX-Toolkit](#fastx-toolkit)
- [RSeQC](#rseqc)
- [Rsubread/Subread](#rsubreadsubread)
  - [Featurecounts](#featurecounts)



### SRA Toolkit 3.0.5  
  \<terminal\>
[Github source](https://github.com/ncbi/sra-tools/wiki)  

## Fastq info
## MultiQC

### Trim Galore  0.6.5  
\<terminal, Bioconda\> 
[Github source](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)
[babraham bioinformat](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

Trim Galore is a wrapper around **Cutadapt** and **FastQC** to consistently apply adapter and quality trimming to FastQ files, with extra functionality for RRBS data.

## Trimmomatic

Flexible read trimming tool for Illumina NGS data

### FASTX-Toolkit   
\<terminal\>  
[Github source](https://github.com/agordon/fastx_toolkit)  

The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.
**Available Tools** 
+ **FASTQ-to-FASTA converter**:Convert FASTQ files to FASTA files.
+ **FASTQ Information**: Chart Quality Statistics and Nucleotide Distribution
+ **FASTQ/A Collapser**: Collapsing identical sequences in a FASTQ/A file into a single sequence (while maintaining reads counts)
+ **FASTQ/A Trimmer**: Shortening reads in a FASTQ or FASTQ files (removing barcodes or noise).
+ **FASTQ/A Renamer**: Renames the sequence identifiers in FASTQ/A file.
+ **FASTQ/A Clipper**: Removing sequencing adapters / linkers
+ **FASTQ/A Reverse-Complement**: Producing the Reverse-complement of each sequence in a FASTQ/FASTA file.
+ **FASTQ/A Barcode splitter**: Splitting a FASTQ/FASTA files containning multiple samples
+ **FASTA Formatte**r: changes the width of sequences line in a FASTA file
+ **FASTA Nucleotide Changer**: Convets FASTA sequences from/to RNA/DNA
+ **FASTQ Quality Filter**: Filters sequences based on quality
+ **FASTQ Quality Trimmer**: Trims (cuts) sequences based on quality
+ **FASTQ Masker**: Masks nucleotides with 'N' (or other character) based on quality


## RSeQC
[Manual](https://rseqc.sourceforge.net/)  
[Github](https://github.com/MonashBioinformaticsPlatform/RSeQC)  
[Author's ](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Identify_Strandedness_Information.html#gsc.tab=0)  
[Galaxy intro](https://github.com/galaxyproject/tools-iuc/tree/main/tools/rseqc)   
[RSeQC genome](https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/)  
An RNA-seq Quality Control Package.  
Usage Information
- **bam2fq.py**: Convert alignments in BAM or SAM format into fastq format.
- **bam2wig.py**: Convert BAM file into wig/bigWig format.
- **bam_stat.py**: Summarizing mapping statistics of a BAM or SAM file.
- **clipping_profile.py**: Calculate the distributions of clipped nucleotides across reads
- **deletion_profile.py**: Calculate the distributions of deletions across reads
- **divide_bam.py**: Equally divide BAM file (m alignments) into n parts. Each part contains roughly m/n alignments that are randomly sampled from total alignments.
- **FPKM_count.py**: Calculate raw read count, FPM (fragment per million) and FPKM (fragment per million mapped reads per kilobase exon) for each gene in BED file. Note: SAM file is not supported.
- **FPKM-UQ.py**: Calculate count, FPKM, and FPKM-UQ values defined by TCGA
- **geneBody_coverage.py**: Calculate the RNA-seq reads coverage over gene body.
- **geneBody_coverage2.py**: Calculate the RNA-seq reads coverage over gene body. This module uses bigwig file as input.
- **infer_experiment.py**:  “guess” how RNA-seq sequencing were configured, particulary how reads were stranded for strand-specific RNA-seq data, through comparing the “strandness of reads” with the “standness of transcripts”.
- **inner_distance.py**: Calculate inner distance between read pairs.
- **insertion_profile.py**: Calculate the distributions of inserted nucleotides across reads. Note that to use this funciton, CIGAR strings within SAM/BAM file should have ‘I’ operation
- **junction_annotation.py**: For a given alignment file (-i) in BAM or SAM format and a reference gene model (-r) in BED format, this program will compare detected splice junctions to reference gene model. splicing annotation is performed in two levels: splice event level and splice junction level.
- **junction_saturation.py**: It’s very important to check if current sequencing depth is deep enough to perform alternative splicing analyses. For a well annotated organism, the number of expressed genes in particular tissue is almost fixed so the number of splice junctions is also fixed. The fixed splice junctions can be predetermined from reference gene model. All (annotated) splice junctions should be rediscovered from a saturated RNA-seq data, otherwise, downstream alternative splicing analysis is problematic because low abundance splice junctions are missing. This module checks for saturation by resampling 5%, 10%, 15%, …, 95% of total alignments from BAM or SAM file, and then detects splice junctions from each subset and compares them to reference gene model.
- **mismatch_profile.py**: Calculate the distribution of mismatches across reads.
- **normalize_bigwig.py**: Visualizing is the most straightforward and effective way to QC your RNA-seq data.
- **overlay_bigwig.py**
- **read_distribution.py**
- **read_duplication.py**
- **read_GC.py**
- **read_hexamer.py**
- **read_NVC.py**
- **read_quality.py**
- **sc_bamStat.py** -h: reports single cell RNA-seq (scRNA-seq) reads mapping statistics. It needs the BAM file generated by the Cell Ranger workflow.
- **sc_seqLogo.py -h**:  It is useful to visualize the nucleotide compositions of “sample barcode”, “cell barcode” and UMI (Unique molecular identifier).
- **sc_seqQual.py**:  
- **sc_editMatrix.py**: This program generates heatmap from a FASTQ file to visualize the sequencing quality.
- **RNA_fragment_size.py**: Calculate fragment size for each gene/transcript. For each transcript, it will report : 1) Number of fragment that was used to estimate mean, median, std (see below). 2) mean of fragment size 3) median of fragment size 4) stdev of fragment size
- **RPKM_count.py**
- **RPKM_saturation.py**: It is particular useful if the input gene list is ribosomal RNA, in this situation, user can estimate how many reads are originated from ribosomal RNA. Download rRNA
- **spilt_bam.py**: spilt_bam.py
- **split_paired_bam.py**: Split bam file (pair-end) into 2 single-end bam file
- **tin.py**: This program is designed to evaluate RNA integrity at transcript level. TIN (transcript integrity number) is named in analogous to RIN (RNA integrity number). RIN (RNA integrity number) is the most widely used metric to evaluate RNA integrity at sample (or transcriptome) level. It is a very useful preventive measure to ensure good RNA quality and robust, reproducible RNA sequencing. However, it has several weaknesses:

## Rsubread/Subread
[Manual](https://subread.sourceforge.net/SubreadUsersGuide.pdf)
### Featurecounts
 a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It can be used to count both RNA-seq and genomic DNA-seq reads. It is available in the SourceForge Subread package or the Bioconductor Rsubread package.

 