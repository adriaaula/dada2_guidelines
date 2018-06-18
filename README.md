# DADA2 pipeline, adapted for our biocluster (and any cluster, in fact)

This repository contains an organized `DADA2` pipeline to be cloned and work
directly with the package. 

Users need to have some `bash` and `R` knowledge along with understanding of the `DADA2` algorithm. 

Before starting, we highly recommend that you read at least one of the following tutorials:

- [DADA2 pipeline official tutorial](https://benjjneb.github.io/dada2/tutorial.html) by *benjjneb*.
- [DADA2 example workflow](https://astrobiomike.github.io/amplicon/dada2_workflow_ex) by *astrobiomike*.

We have divided the pipeline in 4 big steps (plus an initial preprocessing step):

1. Qscore profile plots.
1. `DADA2`.
1. Merge runs and add taxonomy.
1. Check if output seqs contain duplicated reads (100% clustering), or do a lower identity clustering.

## Initial setup

To download the pipeline to your computer/cluster, open the terminal and go to your desired directory with `cd`. Then, type the following (assuming that you have `git` installed):

```sh
    git clone https://github.com/adriaaula/dada2_guidelines.git

```

(You can also download the repository from the [Github server](https://github.com/adriaaula/dada2_guidelines)).

The directory `dada_guidelines/` will be copied to your computer. It contains the following files and subdirectories:

- README.md: this readme file.
- `scripts/`: where all scripts are located.
- `analysis/`: where output files will be written. Contains the `logs/` subdirectory, where log files will be located.
- `data/`: where data files are located. It contains a vanilla dataset.

For everything to work properly, all scripts have to be submitted from the root directory of your project (`dada_guidelines/` in this case). That is, the jobs have to be run (sended to the cluster) from the root directory, not from `scripts/`. 

## Preprocessing step

As explained in *benjjneb*'s tutorial, `DADA2` needs that your sequencing data meets the following criteria: 
- Samples have been demultiplexed (one fastq file per sample).
- Non-biological nucleotides have been removed (e.g. primers, adapters).
- If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order.

In this preprocessing step, we provide the script `XXXX.sh`, which trims primers using [`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html). The script is an array of jobs that outputs the trimmed R1 and R2 files plus a log file (in `analysis/logs/cutadapt/`) for each sample.  

## 0 - Qscore plots

The first step in `DADA2` is to check the quality of your sequencing data. To do so, we provide the script `0_run-qscore.sh`, which calls the R script `0_qscore.R`. Basically, it creates 2 pdf files (forward and reverse) with the q score profile of your first 9 samples (or all your samples if your dataset is smaller). The files generated look like this:

![](https://github.com/adriaaula/dada2_guidelines/.example_graphs/forward.pdf )

![](https://github.com/adriaaula/dada2_guidelines/.example_graphs/reverse.pdf )

Inspecting the quality of your samples will help you decide where to trim them in the following step. It is advisable to always trim (10 bp if your reads are good quality), as you remove the most error-prone regions of your sequences. Given the plots above, we would trim around  **230** for the forward read and **210** for the reverse one.

As a rule of thumb, more trimming will lead to more ASVs and vice versa. 

## 1 - dada2

Having decided where to trim (don't worry, you'll most probably get it wrong the first time), we are able to jump to the next step, which is `DADA2` itself. Here, the script `1_run-dada2.sh` is provided, which calls the script `1_dada2-error-output.R`. A value of max expected error (`maxEE`) for forward and reverse reads has to be provided along the trimming regions (`truncLen`). It is important to note that all reads not reaching the specified lengths will be discarded (so, if you are losing a lot of sequences, maybe you set too high the trimming value).

The main steps of this script are the following:
- Filter (`maxEE`) and trim (`truncLen`). Creates the `data/filtered` directory where it dumps the processed files.
- Learn errors. It creates a forward and reverse pdf files with the plotted error model of your samples.
- Dereplication.
- Sample inference (`dada`).
- Merge paired reads.
- Build sequence table. Creates an ASV table (ASVs as rows, samples as columns)

The script also writes a file to track reads through the pipeline. If you see a big drop in reads in some step, maybe something went wrong.

## 2 - Merge runs & add taxonomy 

In some cases, our dataset is splitted into multiple sequencing runs. Each of them should be processed separatedly with the `dada2` algorithm, and here we will join the outputs (`seqtab` files) into a merged version with the abundance tables. 
See the [Big data](https://benjjneb.github.io/dada2/bigdata.html) for a detailed explanation. 

Additionally, the taxonomy of the various ASVs will be established with the `assignTaxonomy` and `assignSpecies` functions. 
See the [taxonomy tutorial](https://benjjneb.github.io/dada2/assign.html) for further details!

## 3 - Look for duplicates/Clustering

We detected that, in some cases, merging tables from different runs gives reads that are identical but differ in length by a few base pairs. We provide here the script `3_run_cluster_sequences.sh`, which does a 100% identity clustering of the fasta file from the previous step. Additionaly, it creates a table with ASVs to OTUs correspondence.

This script can also be used to do clustering with lower identities (97%, for example). All you have to do is to change this value in the `USEARCH` command inside the script. 
