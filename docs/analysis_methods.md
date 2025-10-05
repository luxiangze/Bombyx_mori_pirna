## Methods

### Small RNA Sequencing Analysis

#### Data Preprocessing and Quality Control

Raw small RNA sequencing reads were subjected to a rigorous quality control workflow. Initial quality assessment was performed using FastQC (v0.11.9). Subsequently, **fastp** (v0.24.0) was employed to trim the 3' adapter sequence (5'-AACTGTAGGCACCATCAAT-3') and filter the reads. Reads were discarded if they had a Phred quality score below 20, contained any ambiguous 'N' bases, or fell outside the 18–37 nucleotide length range. The quality of the final clean dataset was confirmed using FastQC and MultiQC (v1.14).

#### piRNA Identification and Characterization

To identify piRNAs and characterize their properties, clean reads were first converted from FASTQ to FASTA format using **seqkit** (v2.10.0). Identical reads were then collapsed into unique sequences with their corresponding read counts using **fastx_collapser** (FASTX-Toolkit v0.0.14). These unique sequences were aligned to the *Bombyx mori* reference piRNA database (piRBase) using **seqmap** (v1.0.13) under highly stringent conditions, allowing only perfect, forward-strand matches (zero mismatches).

Based on these alignments, we first performed a global characterization of the piRNA populations. The overall length distribution for each sample was plotted to assess general size profiles. Furthermore, to investigate potential post-transcriptional modifications, we analyzed nucleotide variations at the 5' and 3' termini of the mapped piRNAs relative to their reference sequences. The distribution of these terminal variations was compared between treatment and control groups to identify any systemic shifts in piRNA processing. The alignments also generated a count matrix of all confirmed piRNAs, which served as the basis for subsequent differential expression analysis.

#### Differential Expression Analysis of piRNAs

Differential expression analysis between experimental and control groups was performed using a custom R script leveraging the **edgeR** package (v3.34.1). A pairwise comparison approach was employed for each treatment versus control. piRNAs were considered significantly differentially expressed if they met the criteria of a False Discovery Rate (FDR) less than 0.05 and a total read count of at least 100 across the compared samples. The length distributions of both upregulated and downregulated piRNAs were plotted to observe any size preferences.

#### Characterization of Differentially Expressed piRNAs

To further characterize the differentially expressed piRNAs (DE-piRNAs), we analyzed their sequence motifs and genomic origins. Sequence logos for upregulated and downregulated piRNA populations were generated using a custom Python script to identify conserved nucleotide patterns. For genomic origin analysis, DE-piRNA sequences were aligned to the silkworm reference genome (NCBI: ASM3026992v2) using Bowtie. The resulting alignments were converted to BED format, and **intersectBed** (bedtools v2.30.0) was used to determine the overlap between DE-piRNA loci and annotated genomic features (e.g., CDS, UTR, intron, exon). The distribution of these features was then quantified and plotted to reveal the primary sources of the DE-piRNAs.

### Transposon Expression Analysis

#### Transposon Library Construction

To facilitate transposon analysis, an up-to-date transposable element (TE) library was constructed. TEs were identified *de novo* from the latest silkworm reference genome assembly (NCBI: ASM3026992v2) using HiTE (v3.3.3).

#### Quantification of Transposon Expression

Using this custom TE library, the expression levels of transposons were quantified from a separate paired-end RNA-seq dataset. Raw reads were processed for quality control using **fastp**, where adapters were trimmed and low-quality reads were filtered based on multiple criteria, including mean quality score, 'N' base content, and complexity. Clean reads shorter than 50 bp were discarded. The expression analysis was then performed by aligning the clean reads to the reference genome using STAR (v2.7.11b) and quantifying TE expression with the TEcount algorithm from the TEtranscripts software package (v2.2.3).
