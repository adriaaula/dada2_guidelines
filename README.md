# DADA2 pipeline, adapted for our marbits cluster (and any cluster, in fact)

This repository contains an organized `DADA2` pipeline that you can clone and work directly on it.

It provides a directories structure so you can use it as the backbone of your project. You can find more info on structuring projects [here](http://www.riffomonas.org/reproducible_research/organization/#10).

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
    git clone --depth 1 https://github.com/adriaaula/dada2_guidelines.git

```

(You can also download the repository from the [Github server](https://github.com/adriaaula/dada2_guidelines)).

The directory `dada_guidelines/` will be copied to your computer. It contains the following files and subdirectories:

- README.md: this readme file.
- `scripts/`: where all scripts are located.
- `data/`:
    - `logs/`: all log files will be dumped here.
    - `raw/`: it contains a vanilla dataset, initially. It will contain your data. We recommend that instead of copying it, you use symlinks pointing to its original directory. You can easily create a symlink like this: `ln -s /my/original/data/*fastq.gz data/raw/`.
    - `dada2/`: where output files from `DADA2` are written.

For everything to work properly, all scripts have to be submitted from the root directory of your project (`dada_guidelines/` in this case). That is, the jobs have to be run (sent to the cluster) from the root directory, not from `scripts/`. 

## Preprocessing step

As explained in *benjjneb*'s tutorial, `DADA2` needs that your sequencing data meets the following criteria: 
- Samples have been demultiplexed (one fastq file per sample).
- Non-biological nucleotides have been removed (e.g. primers, adapters).
- If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order.

In this preprocessing step, we provide the script `cutadapt.sh`, which trims primers using [`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html). The script outputs the trimmed R1 and R2 files plus a log file (in `analysis/logs/cutadapt/`) for each sample.  

## 0 - Qscore plots

The first step in `DADA2` is to check the quality of your sequencing data. To do so, we provide the script `0_run-qscore.sh`, which calls the R script `0_qscore.R`. Basically, it creates 2 pdf files (forward and reverse) with the q score profile of your first 9 samples (or all your samples if your dataset is smaller). The files generated look like this:


![](https://github.com/adriaaula/dada2_guidelines/blob/master/.examples_output/forward.png)

![](https://github.com/adriaaula/dada2_guidelines/blob/master/.examples_output/reverse.png)


Inspecting the quality of your samples will help you to decide where to trim them in the following step. It is advisable to always trim (10 bp if your reads are good quality), as you remove the most error-prone regions of your sequences. Given the plots above, we would trim around  **230** for the forward read and **210** for the reverse one.

As a rule of thumb, more trimming will lead to more ASVs and vice versa. 

## 1 - dada2

Having decided where to trim (don't worry, you'll most probably get it wrong the first time), we are able to jump to the next step, which is `DADA2` itself. Here, the script `1_run-dada2.sh` is provided, which calls the script `1_dada2-error-output.R`. A value of max expected error (`maxEE`) for forward and reverse reads has to be provided along the trimming regions (`truncLen`). It is important to note that all reads not reaching the specified lengths will be discarded (so, if you are losing a lot of sequences, maybe you set too high the trimming value).

The main steps of this script are the following:
- Filter (`maxEE`) and trim (`truncLen`). Creates the `data/filtered` directory where it dumps the processed files.
- Learn errors. It creates a forward and reverse pdf files with the plotted error model of your samples. It is vital to check the results, since a bad error model will bring problems afterwards! An example:

![](https://benjjneb.github.io/dada2/tutorial_files/figure-html/plot-errors-1.png) 

- Dereplication.
- Sample inference (`dada`).
- Merge paired reads.
- Build sequence table. Creates an ASV table (ASVs as rows, samples as columns)

The script also writes a file to track reads through the pipeline. If you see a big drop in reads in some step, maybe something went wrong.

## 2 - Merge runs & add taxonomy 

In some cases, our dataset is splitted into multiple sequencing runs. Each of them should be processed separatedly with the `dada2` algorithm, and here we will join the outputs (`seqtab` files) into a merged version with the abundance tables. 
See the [Big data](https://benjjneb.github.io/dada2/bigdata.html) for a detailed explanation. 

Additionally, the taxonomy of ASVs will be established with the `assignTaxonomy` and `assignSpecies` functions. 
See the [taxonomy tutorial](https://benjjneb.github.io/dada2/assign.html) for further details!

## 3 - Look for duplicates/Clustering

We detected that, in some cases, merging tables from different runs gives reads that are identical but differ in length by a few base pairs. We provide here the script `3_run_cluster_sequences.sh`, which does a 100% identity clustering of the fasta file from the previous step. Additionaly, it creates a table with ASVs to OTUs correspondence.

This script can also be used to do clustering with lower identities (97%, for example). All you have to do is to change this value in the `USEARCH` command inside the script. 

# Some general rules 

Some tips about the parameter selection!

- Cut your primers. `cutadapt` does the job really easy! 

- If around 25-30 % of the reads are lost in the process of ASV generation, possibly some of the parameters have to be changed. 

	* Are you sure that the primers from the FASTQ are removed?

 * What `maxee` did you specify? If this is making many reads to be lost, you can specify a bigger maxee, and even different values for the F and R reads (for example `c(2,4)`).
The algorithm will take into account the errors in the modelling phase, so this will not make your ASVs erroneus. 

 * Does the pair of reads overlapp? By how many bases? It should be >= 20 nt. 

> If you follow the tutorial, at the end of the procedure a **track analysis** is generated specifing how many reads are lost along the whole procedure. It is the best way to know where it failed. 

- In the trimming procedure, the `truncLen` cuts all the reads to an specific length and *removes* all reads being smaller.  It is important then to know the average read length, since if you go too low with the trimming you will lose too much reads. 

 * In the pipeline of DADA2 there is a quality profile, you should be aware of it in deciding where to cut.
 
 * For each run the trimming point is different, so if you are working on multiple runs each of them have to be processed separatedly and then joined together with `mergeSequenceTables`. 

 * You should have an analysis of the FASTQs. The av. length, the avg quality for each sample, and so on. Many of the problems with recovering most of the reads
stem from having a low quality sample, or the reads not being properly amplified. [`seqkit`](https://github.com/shenwei356/seqkit) is a good tool for this kind of information. 

- The taxonomy assignation is done at the Species level only if only a 100%, exact matching. This can make that some Bacteria/Eukarya present differences
in the identification at that level when comparing with OTU results. See a link explaining this in more detail [here](https://benjjneb.github.io/dada2/assign.html#species-assignment).


**Adrià & Aleix**

*Backbone team*

:skull:
