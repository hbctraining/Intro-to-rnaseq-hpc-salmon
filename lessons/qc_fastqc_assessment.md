---
title: "Quality control: Assessing FASTQC results"
author: "Mary Piper, Radhika Khetani"
date: Wednesday, September 5, 2018
duration: 45 minutes
---

## Learning Objectives:

* Evaluate the quality of your NGS data using FastQC

## Assessing quality metrics	

Now that we have run and downloaded the FASTQC report, we can take a look at the metrics and assess the quality of our sequencing data.

FastQC has a really well documented [manual page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) with [detailed explanations](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) about all the plots in the report. 

A summary of all of the modules is given on the left-hand side of the report. Don't take the yellow "WARNING"s and red "FAIL"s too seriously, these are more indicative of a module that you should take a look at to make sure there is nothing worrisome about the data. 

The first module gives the basic statistics for the sample. Generally it is good to keep track of the total number of reads for each sample and make sure the read length and %GC content is expected.

<img src>

One of the most important analysis modules is the **"Per base sequence quality"** plot. This plot provides the distribution of quality scores at each position in the read across all reads. This plot can alert us to whether there were any problems occuring during sequencing, and whether we might need to contact the sequencing facility.

![FastQC_seq_qual](../img/FastQC_seq_qual.png)

The y-axis gives the quality scores, while the x-axis represents the position in the read. The color coding of the plot is denotes what are considered high, medium and low quality scores. 

For example, the box plot at nucleotide 1 shows the distribution of quality scores for the first nucleotide of all reads in the `Mov10_oe_1` sample. The yellow box represents the 25th and 75th percentiles, with the red line as the median. The whiskers are the 10th and 90th percentiles. The blue line represents the average quality score for the nucleotide. The quality scores for the first nucleotide are quite high, with nearly all reads having scores above 28.

The quality scores appear to drop going from the beginning toward the end of the reads. This is not unexpected, and there are known causes for this drop in quality. To better interpret this plot it is helpful to understand the different sequencing error profiles, including expected and unexpected quality issues.

### Quality error profiles



The **"Overrepresented sequences"** table displays the sequences (at least 20 bp) that occur in more than 0.1% of the total number of sequences. This table aids in identifying contamination, such as vector or adapter sequences. 

![FastQC_contam](../img/FastQC_contam.png)


***FastQC is just an indicator of what's going on with your data, don't take the "PASS"es and "FAIL"s too seriously.***
We recommend looking at [this post](http://bioinfo-core.org/index.php/9th_Discussion-28_October_2010) for more information on what bad plots look like and what they mean for your data.

We will go over the remaining plots in class. Remember, our report only represents a subset of reads (chromosome 1) for `Mov10_oe_1.subset.fq`, which can skew the QC results. We encourage you to look at the [full set of reads](../fastqc/Mov10oe_1-fastqc_report.html) and note how the QC results differ when using the entire dataset.
   
> **_NOTE:_** 
>The other output of FastQC is a .zip file. These .zip files need to be unpacked with the `unzip` program. If we try to `unzip` them all at once:
>
>```bash
>$ cd ~/unix_lesson/rnaseq/results/fastqc/    
>$ unzip *.zip
>```
>
>Did it work? 
>
>No, because `unzip` expects to get only one zip file. Welcome to the real world.
>We *could* do each file, one by one, but what if we have 500 files? There is a smarter way.
>We can save time by using a simple shell `for loop` to iterate through the list of files in *.zip.
>
>After you type the first line, you will get a special '>' prompt to type next lines.  
>You start with 'do', then enter your commands, then end with 'done' to execute the loop.
>
>Note that in the first line, we create a variable named `zip`.  After that, we call that variable with the syntax `$zip`. `$zip` is assigned the value of each item (file) in the list *.zip, once for each iteration of the loop.
>
>This loop is basically a simple program. When it runs
>
>```bash
>$ for zip in *.zip
> do
> unzip $zip
> done
>```
>This will run unzip once for each file (whose name is stored in the $zip variable). The contents of each file will be unpacked into a separate directory by the unzip program.
>
>The 'for loop' is interpreted as a multipart command.  If you press the up arrow on your keyboard to recall the command, it will be shown like so:
>
>```bash
>for zip in *.zip; do unzip $zip; done
>```
>
>When you check your history later, it will help you remember what you did!
>
>What information is contained in the unzipped folder?
>
>```bash
>$ ls -lh Mov10_oe_1.subset_fastqc
>$ head Mov10_oe_1.subset_fastqc/summary.txt
>```
>
>To save a record, let's `cat` all `fastqc summary.txt` files into one `full_report.txt` and move this to `~/unix_lesson/rnaseq/docs`. 
>You can use wildcards in paths as well as file names.  Do you remember how we said `cat` is really meant for concatenating text files?
>    
>```bash
>$ cat */summary.txt > ~/unix_lesson/rnaseq/logs/fastqc_summaries.txt
>```

## Quality Control (*Optional*) - Trimming 

We want to make sure that as many reads as possible map or align accurately to the genome. To ensure accuracy, only a small number of mismatches between the read sequence and the genome sequence are allowed, and any read with more than a few mismatches will be marked as being unaligned. 

Therefore, to make sure that all the reads in the dataset have a chance to map/align to the genome, unwanted information can be trimmed off from every read, one read at a time. The types of unwanted information can include one or more of the following:
- leftover adapter sequences
- known contaminants (strings of As/Ts, other sequences)
- poor quality bases at read ends

**We will not be performing this step** because:
* our data does not have an appreciable amount of leftover adapter sequences or other contaminating sequences based on FastQC.
* the alignment tool we have picked (STAR) is able to account for low-quality bases at the ends of reads when matching them to the genome.

If you need to perform trimming on your fastq data to remove unwanted sequences/bases, the recommended tool is [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html). 

Example of cutadapt usage:

```bash
$ cutadapt --adapter=AGATCGGAAGAG --minimum-length=25  -o myfile_trimmed.fastq.gz myfile.fastq.gz 
```

After trimming, cutadapt can remove any reads that are too short to ensure that you do not get spurious mapping of very short sequences to multiple locations on the genome. In addition to adapter trimming, cutadapt can trim off any low-quality bases too, but **please note that quality-based trimming is not considered best practice, since majority of the newer, recommended alignment tools can account for this.**

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
