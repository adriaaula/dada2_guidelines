# DADA2 pre-processing. Q score distribution

## Load packages

library(stringr)
library(tidyverse)
library(dada2); packageVersion("dada2")

## Load Rscript args

args <- commandArgs(trailingOnly = TRUE)

## Variables to change for EACH RUN

path<- args[1] # where the fastq files are located
output <- args[2] # where we will dump the results//diagnostic
run.name <- args[3] # Run name. A directory for the values will be created

## Taking filenames

fnFs <- sort(list.files(path, pattern = 'R1.fastq'))
fnRs <- sort(list.files(path, pattern= 'R2.fastq'))

sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) # BEWARE with the field separator in strsplit() command, can change between runs

fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

### Since a qprofile dir can be not present we check and create it

dir.create(file.path(output, "0_qprofiles"), showWarnings = FALSE)
dir.create(file.path(output, "0_qprofiles",run.name), showWarnings = FALSE)

out.diag <- file.path(output, "0_qprofiles", run.name)

## Qscore distribution

# Forward
ggsave(plot = plotQualityProfile(fnFs[1:9]) +
         ggtitle("Forward reads"), 
       path= out.diag,
       device="pdf",
       filename = "forward.pdf",
       width = 410,
       height = 300,
       units = 'mm')

# Same for reverse reads! 

ggsave(plot=plotQualityProfile(fnRs[1:9]) +
         ggtitle("Reverse reads"),
       path= out.diag,
       device="pdf",
       filename = "reverse.pdf",
       width = 410,
       height = 300,
       units = 'mm')

# The scripts generates an artifactural pdf in the main directory. 
# We are going to remove it. It is a bad solution, but a solution after all!
if (file.exists('Rplots.pdf')) file.remove('Rplots.pdf')

# END