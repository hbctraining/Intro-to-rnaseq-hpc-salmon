---
title: "The Shell"
author: "Mary Piper, Radhika Khetani, Meeta Mistry"
date: "February 2019"
---

## Learning Objectives
- Review shell commands and concepts

## Setting up

### Connecting to a *login* node on O2

Type in the following command with your username to login:

```bash
ssh username@o2.hms.harvard.edu
```

You will receive a prompt for your password, and you should type in your associated password; note that the cursor will *not move* as you type in your password.

A warning might pop up the first time you try to connect to a remote machine, type "Yes" or "Y". 

Once logged in, you should see the O2 icon, some news, and the command prompt: 

```bash
[rc_training10@login01 ~]$ 
```

> A "node" on a cluster is essentially a computer in the cluster of computers.
>
> A login node is only to enable users to log in, it is not meant to be used for any actual work/computing.

### Connecting to a *compute* node on O2

There are multiple ways to connect with, and do work on, a compute node; a compute node is where all work should be performed. To connect to a compute node, users have to interact with a job scheduler like *slurm* using commands like `srun` or `sbatch` and by specifying what resources they require.

1. The `srun` command with a few mandatory parameters will create an "interactive session" on O2. This is essentially a way for us to do work on the compute node directly from the terminal. If the connectivity to the cluster is lost in the middle of a command being run that work will be lost in an interactive session.

2. The `sbatch` command with a few mandatory parameters along with a specialized shell script will result in the script being run on a compute node. This "job" will not be accessible directly from the Terminal and will run in the background. Users do not need to remain connected to the cluster when such a job "batch job" is running.

You will get practice to run and check on batch jobs, for now we are going to start an interactive session on O2 using `srun`. 

```bash
$ srun --pty -p short -t 0-8:00 --mem 1G --reservation=HBC /bin/bash
```

In the above command the parameters we are using are requesting specific resources:
* `--pty` - Start an interactive session
* `-p short` - on the "partition" called "short" (a partition is a group of computers dedicated to certain types of jobs, long, short, high-memory, etc.)
* `-t 0-8:00` - time needed for this work: 0 days, 8 hours, 0 minutes.
* `--mem 1G` - memory needed - 1 gigabyte
* `--reservation=HBC` - *this is only for this workshop, make sure you don't use it in the future with your own accounts*
* `/bin/bash` - You want to interact with the compute node using the *bash* shell

Make sure that your command prompt is now preceded by a character string that contains the word "compute":

```bash
[rc_training10@compute-a-16-163 ~]$
```

### Copying example data folder

Your accounts were erased after last time, so we are starting fresh this time, let's copy over the data folder we worked with in the shell workshop to our new home directories:

```bash
$ cp -r /n/groups/hbctraining/unix_lesson/ .
```

****

**Exercise**

1. In the above command, what does the `.` at the end mean? Is it essential?
2. Why did we have to run the `cp` command with `-r`?
3. Is the path to the `unix_lesson/` directory a "full" path or a "relative" path?

****

## Reviewing shell commands

We are going to start this review with more exercises, this time hands on! Remember, there are likely multiple ways to do the same thing and we will try to go over some of them.

****

**Exercises**

**Shell basics**
1. Change directory into the `unix_lesson/` directory using a relative path.
2. Use the `tree` command to get a layout of `unix_lesson/`.
3. Take a quick look at the `Mov10_oe_1.subset.fq` file using `less` from `unix_lesson/` without changing directories.
4. Move up to your home directory (parent of `unix_lesson/`).
5. With a single command change directories to the `raw_fastq/` folder.
6. With a shortest possible command change directories back to the home directory.
7. What does the `~` in the command prompt mean?
8. What is the full path to your home directory?
9. List (long listing format) the contents of `/n/groups/hbctraining/intro_rnaseq_hpc/full_dataset/` using tab completion.
10. Modify the above command using the `*` wildcard to only list those files that have "oe" in their names.
11. How many and which commands have you run today?

**Searching and redirection**
12. How many lines are in `~/unix_lesson/reference_data/chr1-hg19_genes.gtf`?
 * How many of those lines have the string "MOV10" in them?
 * How many of the lines with the string "MOV10" in them have the word "exon" in them?
13. Create a new directory called `shell_review/` within the `unix_lesson/` directory.
14. Grab the lines in `~/unix_lesson/reference_data/chr1-hg19_genes.gtf` with the string "MOV10" in them and save it in the `shell_review/` directory with a new name "Mov10_hg19.gtf".
15. Use `vim` to open the newly created file `~/unix_lesson/shell_review/Mov10_hg19.gtf` and add a comment at the top specifying how this file was created and the source of the content. Save the modified file and come back to the command prompt.

**Loops and shell scripts**
16. Use the `for` loop to iterate over each FASTQ file in `~/unix_lesson/raw_fastq/` and do the following:
 * Generate a prefix to use for naming our output files
 * Print the name of the current file
 * Dump out the first 40 lines into a new file that will be saved in `~/unix_lesson/shell_review/`
17. Place the above `for` loop into a shell script using `vim` and run it.

**Permissions**
18. List `/n/groups/hbctraining/intro_rnaseq_hpc/` directory in long listing format
 * How many owners have files in this folder?
 * How many groups?
 * Are there any executable *files* in this folder?
 * What kind of access do you have to the `full_dataset/` directory?
 * What could user `mm573` do to take away your ability to look inside the `full_dataset/` directory?

**Environment variables**
19.

****

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
