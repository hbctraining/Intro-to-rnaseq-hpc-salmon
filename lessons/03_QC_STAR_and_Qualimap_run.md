---
title: "QC with STAR and Qualimap"
author: "Meeta Mistry, Mary Piper"
date: Monday September 10, 2018
---

Approximate time: 50 minutes

## Learning Objectives:

* Running an alignment tool to generate BAM files
* Brief explanation about SAM file
* Running Qualimap to compute metrics on alignment files


## Alignment of Raw Counts

<img src="../img/workflow_align_qualimap.png">

Now that we have explored the quality of our raw reads, we can align the raw reads to the genome to explore other quality metrics, such as DNA or rRNA contamination, 5'-3' biases, and coverage biases. We perform read alignment or mapping to determine where in the genome the reads originated from. To explore these QC metrics, we need to use a traditional, splice-aware alignment tool that outputs an alignment file, BAM or SAM, with information on the genome coordinates for where the entire read mapped.

The alignment process consists of choosing an appropriate reference genome to map our reads against and performing the read alignment using one of several splice-aware alignment tools such as [STAR](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635) or [HISAT2](http://ccb.jhu.edu/software/hisat2/index.shtml). The choice of aligner is often a personal preference and also dependent on the computational resources that are available to you. 

>**NOTE:** The latest human genome build, GRCh38, contains information about alternative alleles for various locations on the genome. If using this genome then it is advisable to use the HISAT2 aligner as it is able to utilize this information during the alignment. STAR, however, will need to use a genome that does not have the alternate alleles present, as it does not have the functionality to deal with the alternate alleles.

## STAR Aligner

To determine where on the human genome our reads originated from, we will align our reads to the reference genome using [STAR](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments. 

### STAR Alignment Strategy

STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. The algorithm achieves this highly efficient mapping by performing a two-step process:

1. Seed searching
2. Clustering, stitching, and scoring

#### Seed searching

For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs):

![STAR_step1](../img/alignment_STAR_step1.png)
	
The different parts of the read that are mapped separately are called 'seeds'. So the first MMP that is mapped to the genome is called *seed1*.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be *seed2*. 

![STAR_step2](../img/alignment_STAR_step2.png)

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping.

**If STAR does not find an exact matching sequence** for each part of the read due to mismatches or indels, the previous MMPs will be extended.

![STAR_step3](../img/alignment_STAR_step3.png)

**If extension does not give a good alignment**, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

![STAR_step4](../img/alignment_STAR_step4.png)


#### Clustering, stitching, and scoring

The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of 'anchor' seeds, or seeds that are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.). 

![STAR_step5](../img/alignment_STAR_step5.png)

## Running STAR

### Set-up

To get started with this lesson, start an interactive session with 6 cores:

```bash
$ srun --pty -p short -t 0-12:00 -n 6 --mem 8G --reservation=HBC /bin/bash	
```

You should have a directory tree setup similar to that shown below. It is best practice to have all files you intend on using for your workflow present within the same directory. In our case, we have our original FASTQ files generated in the previous section. 

```bash
rnaseq
	├── logs
	├── meta
	├── raw_data
	│   ├── Irrel_kd_1.subset.fq
	│   ├── Irrel_kd_2.subset.fq
	│   ├── Irrel_kd_3.subset.fq
	│   ├── Mov10_oe_1.subset.fq
	│   ├── Mov10_oe_2.subset.fq
	│   └── Mov10_oe_3.subset.fq
	├── results
	└── scripts
```

To use the STAR aligner, load the module: 

```bash
$ module load gcc/6.2.0 star/2.5.4a
```

Aligning reads using STAR is a two step process:   

1. Create a genome index 
2. Map reads to the genome

> A quick note on shared databases for human and other commonly used model organisms. The O2 cluster has a designated directory at `/n/groups/shared_databases/` in which there are files that can be accessed by any user. These files contain, but are not limited to, genome indices for various tools, reference sequences, tool specific data, and data from public databases, such as NCBI and PDB. So when using a tool that requires a reference of sorts, it is worth taking a quick look here because chances are it's already been taken care of for you. 
>
>```bash
> $ ls -l /n/groups/shared_databases/igenome/
>```

### Creating a genome index

For this workshop we are using reads that originate from a small subsection of chromosome 1 (~300,000 reads) and so we are using only chr1 as the reference genome from hg38 without the alternative alleles. 

The basic options to **generate genome indices** using STAR are as follows:

* `--runThreadN`: number of threads
* `--runMode`: genomeGenerate mode
* `--genomeDir`: /path/to/store/genome_indices
* `--genomeFastaFiles`: /path/to/FASTA_file 
* `--sjdbGTFfile`: /path/to/GTF_file
* `--sjdbOverhang`: readlength -1

> *NOTE:* In case of reads of varying length, the ideal value for `--sjdbOverhang` is max(ReadLength)-1. In most cases, the default value of 100 will work similarly to the ideal value.

The final command to create an index can be found in the job submission script we have linked [here](../scripts/star_genome_index.sh). We have generated the genome indices for you, so that we don't get held up waiting on the generation of the indices. The index can be found in the `/n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/ensembl38_STAR_index/` directory. 

### Aligning reads

After you have the genome indices generated, you can perform the read alignment. 

Create an output directory for our alignment files:

```bash
$ cd ~/rnaseq/raw_data

$ mkdir ../results/STAR
```

### Running STAR interactively

For now, we're going to work on just one sample to set up our workflow. To start we will use the first replicate in the Mov10 over-expression group, `Mov10_oe_1.subset.fq`. Details on STAR and its functionality can be found in the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf); we encourage you to peruse through to get familiar with all available options.

The basic options for aligning reads to the genome using STAR are:

* `--runThreadN`: number of threads / cores
* `--readFilesIn`: /path/to/FASTQ_file
* `--genomeDir`: /path/to/genome_indices_directory
* `--outFileNamePrefix`: prefix for all output files

Listed below are additional parameters that we will use in our command:

* `--outSAMtype`: output filetype (SAM default)
* `--outSAMunmapped`: what to do with unmapped reads

> **NOTE:** Default filtering is applied in which the maximum number of multiple alignments allowed for a read is set to 10. If a read exceeds this number there is no alignment output. To change the default you can use `--outFilterMultimapNmax`, but for this lesson we will leave it as default. Also, note that "**STAR’s default parameters are optimized for mammalian genomes.** Other species may require significant modifications of some alignment parameters; in particular, the maximum and minimum intron sizes have to be reduced for organisms with smaller introns" [[1](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635.full.pdf+html)].

The full command is provided below for you to copy paste into your terminal. If you want to manually enter the command, it is advisable to first type out the full command in a text editor (i.e. [Sublime Text](http://www.sublimetext.com/) or [Notepad++](https://notepad-plus-plus.org/)) on your local machine and then copy paste into the terminal. This will make it easier to catch typos and make appropriate changes. 

```bash

$ STAR --genomeDir /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/ensembl38_STAR_index/ \
--runThreadN 6 \
--readFilesIn Mov10_oe_1.subset.fq \
--outFileNamePrefix ../results/STAR/Mov10_oe_1_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

```

## Assessing alignment quality

After running our single FASTQ file through the STAR aligner, you should have a number of output files in the `~/rnaseq/results/STAR` directory. Let's take a quick look at some of the files that were generated and explore their content. 

	$ cd ~/rnaseq/results/STAR
	
	$ ls -lh

What you should see, is that for each FASTQ file you have **5 output files** and a single tmp directory. Briefly, these files are described below:

* `Log.final.out` - a summary of mapping statistics for the sample
* `Aligned.sortedByCoord.out.bam` - the aligned reads, sorted by coordinate, in BAM format
* `Log.out` - a running log from STAR, with information about the run 
* `Log.progress.out` -  job progress with the number of processed reads, % of mapped reads etc., updated every ~1 minute
* `SJ.out.tab` - high confidence collapsed splice junctions in tab-delimited format. Only junctions supported by uniquely mapping reads are reported

### Alignment file format: SAM/BAM

The output we requested from the STAR aligner (using the appropriate parameters) is a BAM file. By default STAR will return a file in SAM format. BAM is a binary, compressed version of the SAM file, also known as **Sequence Alignment Map format**. The SAM file is a tab-delimited text file that contains all information from the FASTQ file, with additional alignment information for each read. We will explore these files in more detail in later sessions, but the paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

![SAM1](../img/sam_bam.png)

### Mapping statistics

To determine whether we have any contamination or biases in our data, we want to know is how well did our reads align to the reference. Rather than looking at each read alignment, it can be more useful to evaluate statistics that give a general overview for the sample. The `Log.final.out` file output from STAR contains mapping statistics. Let's use the `less` command to scroll through it: 

	$ less Mov10_oe_1_Log.final.out
	
The log file provides information on reads that 1) mapped uniquely, 2) reads that mapped to mutliple locations and 3) reads that are unmapped. Additionally, we get details on splicing, insertion and deletion. From this file the most informative statistics include the **mapping rate and the number of multimappers**. This information will be compiled by MultiQC downstream in report form.

In addition to the aligner-specific summary, we can also obtain quality metrics using tools like [Qualimap](http://qualimap.bioinfo.cipf.es/doc_html/intro.html#what-is-qualimap) or [RNASeQC](http://archive.broadinstitute.org/cancer/cga/rna-seqc). 

### Qualimap

The Qualimap tool is written in Java and R and explores the features of mapped reads and their genomic properties. Qualimap **provides an overall view of the data that helps to detect biases in the sequencing and/or mapping of the data**. The input can be one or more BAM files and the output consists of HTML or PDF reports with useful figures and tab delimited files of metrics data.

To run Qualimap, change directories to the `rnaseq` folder and make a `qualimap` folder inside the `results` directory:

```bash
$ cd ~/rnaseq

$ mkdir -p results/qualimap
```

By default, Qualimap will try to open a GUI to run Qualimap, so we need to run the `unset DISPLAY` command:

```bash
$ unset DISPLAY
```

We also need to add the location of the Qualimap tool to our PATH variable:

```bash
$ export PATH=/n/app/bcbio/dev/anaconda/bin:$PATH
```

Now we can run Qualimap on our aligned data. There are different tools or modules available through Qualimap, and the [documentation](http://qualimap.bioinfo.cipf.es/doc_html/command_line.html) details the tools and options available. We are interested in the `rnaseq` tool. To see the arguments available for this tool we can search the help:

```bash
$ qualimap rnaseq --help
```

 We will be running Qualimap with the following specifications:

- `-outdir`: output directory for html report
- `-a`: Counting algorithm - uniquely-mapped-reads(default) or proportional
- `-bam`: path/to/bam/file(s)
- `-p`: Sequencing library protocol - strand-specific-forward, strand-specific-reverse or non-strand-specific (default)
- `-gtf`: path/to/gtf/file - **needs to match the genome build and GTF used in alignment**
-  `--java-mem-size=`: set Java memory

```bash
$ qualimap rnaseq \
-outdir results/qualimap/Mov10_oe_1 \
-a proportional \
-bam results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam \
-p strand-specific-reverse \
-gtf /groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/Homo_sapiens.GRCh38.92.chr1.gtf\
--java-mem-size=8G
```

The Qualimap report should be present in our `results/qualimap` directory. To view this report we would need to transfer it over to our local computers. However, this report was generated on a subset of data on chromosome 1; it would be better to visualize the report of the full dataset for `Mov_oe_1`, which is available [here]().


---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
