---
title: "Intro to RNA-seq"
author: "Mary Piper, Meeta Mistry, Radhika Khetani"
date: "Monday, February 25, 2019"
---

Approximate time: 90 minutes

## Learning Objectives:

* Describe the process of RNA-seq library preparation
* Describe the Illumina sequencing method

## Introduction to RNA-seq

RNA-seq is an exciting experimental technique that is utilized to explore and/or quantify gene expression within or between conditions. 


As we know, genes provide instructions to make proteins, which perform some function within the cell. Although all cells contain the same DNA sequence, muscle cells are different from nerve cells and other types of cells because of the different genes that are turned on in these cells and the different RNAs and proteins produced. 

 <img src="../img/gene_expression_cells.png" width="600">

Different biological processes, as well as mutations, can affect which genes are turned on and which are turned off, in addition to, *how much* specific genes are turned on/off.

To make proteins, the DNA is transcribed into messenger RNA, or mRNA, which is translated by the ribosome into protein. However, some genes encode RNA that does not get translated into protein; these RNAs are called non-coding RNAs, or ncRNAs. Often these RNAs have a function in and of themselves and include rRNAs, tRNAs, and siRNAs, among others. All RNAs transcribed from genes are called transcripts.

 <img src="../img/Gene_products.png" width="400">

To be translated into proteins, the RNA must undergo processing to generate the mRNA. In the figure below, the top strand in the image represents a gene in the DNA, comprised of the untranslated regions (UTRs) and the open read frame. Genes are transcribed into pre-mRNA, which still contains the intronic sequences.  After post-transciptional processing, the introns are spliced out and a polyA tail and 5' cap are added to yield mature mRNA transcripts, which can be translated into proteins.

 <img src="../img/Gene_structure.png" width="600">

**While mRNA transcripts have a polyA tail, many of the non-coding RNA transcripts do not as the post-transcriptional processing is different for these transcripts.**

### Transcriptomics

The transcriptome is defined as a collection of all the transcript readouts present in a cell. RNA-seq data can be used to explore and/or quantify the transcriptome of an organism, which can be utilized for the following types of experiments:

- **Differential Gene Expression**: *quantitative* evaluation and comparison of transcript levels
- **Transcriptome assembly**: building the profile of transcribed regions of the genome, a *qualitative* evaluation. 
- Can be used to **help build better gene models**, and verify them using the assembly
- **Metatranscriptomics** or community transcriptome analysis


## Illumina library preparation

When starting an RNA-seq experiment, for every sample the RNA needs to be isolated and turned into a cDNA library for sequencing. Generally, ribosomal RNA represents the majority of the RNAs present in a cell, while messenger RNAs represent a small percentage of total RNA, ~2% in humans.

Therefore, if we want to study the protein-coding genes, we need to enrich for mRNA or deplete the rRNA. **For differential gene expression analysis, it is best to enrich for Poly(A)+, unless you are aiming to obtain information about long non-coding RNAs, then do a ribosomal RNA depletion.**

The workflow for library preparation is detailed in the step-by-step images below.

**Briefly, the RNA is isolated from the sample and contaminating DNA is removed with DNase.**

 <img src="../img/libraryprep_step1-2.png" width="400">
 
**The RNA sample then undergoes either selection of the mRNA (polyA selection) or depletion of the rRNA. The resulting RNA is fragmented.** 

> *The size of the target fragments in the final library is a key parameter for library construction. DNA fragmentation is typically done by physical methods (i.e., acoustic shearing and sonication) or enzymatic methods (i.e., non-specific endonuclease cocktails and transposase tagmentation reactions.*
 
 <img src="../img/libraryprep_step3.png" width="400">

**The RNA is then reverse transcribed into double-stranded cDNA and sequence adapters are then added to the ends of the fragments.**

> The cDNA libraries can be generated in a way to retain information about which strand of DNA the RNA was transcribed from. Libraries that retain this information are called stranded libraries, which are now standard with Illumina’s TruSeq stranded RNA-Seq kits. Stranded libraries should not be any more expensive than unstranded, so there is not really any reason not to acquire this additional information. 
> 
> There are 3 types of cDNA libraries available:
> 
> * Forward (secondstrand) – reads resemble the gene sequence or the secondstrand cDNA sequence
> * Reverse (firststrand) – reads resemble the complement of the gene sequence or firststrand cDNA sequence (TruSeq)
> * Unstranded

 <img src="../img/libraryprep_step4-5.png" width="400">
 
**Finally, the fragments are PCR amplified if needed, and the fragments are size selected (usually ~300-500bp) to finish the library.**

 <img src="../img/libraryprep_step6.png" width="400">
 
 *Image credit: [Martin J.A. and Wang Z., Nat. Rev. Genet. (2011) 12:671–682](https://www.nature.com/articles/nrg3068)*


> **NOTE:** This workflow is specific to Illumina sequencing, which is currently the most utilized sequencing method. But there are other long-read methods worth noting, such as:
>
> - Pacific Biosciences: http://www.pacb.com/ 
> - Oxford Nanopore (MinION): https://nanoporetech.com/ 
> - 10X Genomics: https://www.10xgenomics.com/ 
> 
> Advantages and disadvantages of these technologies can be explored in the table below:
>
> <img src="../img/long_read_tech.png" width="700">

## Illumina Sequencing

After preparation of the libraries, sequencing can be performed to generate the nucleotide sequences of the ends of the fragments, which are called **reads**. You will have the choice of sequencing a single end of the cDNA fragments (single-end reads) or both ends of the fragments (paired-end reads).

<img src="../img/paired_end_reads.png" width="500">

- SE - Single end dataset => Only Read1
- PE - Paired-end dataset => Read1 + Read2
	- can be 2 separate FastQ files or just one with interleaved pairs

Generally single-end sequencing is sufficient unless it is expected that the reads will match multiple locations on the genome (e.g. organisms with many paralogous genes), assemblies are being performed, or for splice isoform differentiation. Be aware that paired-end reads are generally 2x more expensive.

There are a variety of Illumina platforms to choose from to sequence the cDNA libraries.

<img src="../img/illumina_platforms.png" width="800">

 *Image credit: Adapted from [Illumina](www.illumina.com)*

Differences in platform can alter the length of reads generated as well as the total number of reads sequenced per run and the amount of time required to sequence the libraries. The different platforms each use a different flow cell, which is a glass surface coated with an arrangement of paired oligos that are complementary to the adapters added to your template molecules. The flow cell is where the sequencing reactions take place.

<img src="../img/flow_cells.png" width="400">

 *Image credit: Adapted from [Illumina](www.illumina.com)*
 

Let's explore how Illumina sequencing is performed:

[<img src="../img/illumina_sequencing.png" width="400">](https://www.dropbox.com/s/f4t94tcw06f9stg/Illumina%20Sequencing%20by%20Synthesis-14840.mp4?dl=0)

- Number of clusters ~= Number of reads
- Number of sequencing cycles = Length of reads

The number of cycles (length of the reads) will depend on sequencing platform used as well as your preferences.

Charges for sequencing are usually per lane of the flow cell, and usually you don’t need one lane per sample. Multiplexing allows you to sequence multiple samples per lane with addition of indices (within the Illumina adapter) or special barcodes (outside the Illumina adapter).

<img src="../img/demultiplexing.png" width="800">
