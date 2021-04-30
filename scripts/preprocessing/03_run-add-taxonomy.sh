#!/bin/sh

#SBATCH --account=<your-account>
#SBATCH --job-name=add_tax
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=data/logs/3_add_taxonomy_%J.out
#SBATCH --error=data/logs/3_add_taxonomy_%J.err

### DADA2 pipeline ######
## ~~ Trimming, error generation and DADA2 run ~~ ##

#Trimming of the reads (and filtering with maxEE),
#generation of an error model whit the probabilities
#transitions between bases (A -> T) and dada2 algorithm.
#BEFORE continue the pipeline, PLEASE, check the error plot
#(a .pdf generated in the output dir).

## ARGS ##

#[1] /seqtab/ 
#   seqtab with chimeras removed.

#[2] /output dir/ 
#   Directory output, usually data (a subdirectory will be created).

#[3] /name/ 
#   A common identifier to be sure that the output is the
#   correct one. You will thank us that :^)

#[4] /Taxonomy db/ 
#   Here you have to put your database[s] for classification.
#   if using dada's assignTaxonomy() and addSpecies() functions, put your 2 databases separated by a comma
#   if using only assignTaxonomy(), put your database   
#   if using DECIPHER, the db has to be downloaded from http://www2.decipher.codes/Downloads.html and put here

#[5] /Taxonomy classification method/ 
#   write 'decipher' if you want to use DECIPHER's IdTaxa() to classify
#   write 'dada' if you want to use the classifier included in dada2

#[6] /Confidence level of classification (DECIPHER) or minBoot (dada2)/
#   DECIPHER (default 60):       
#       Numeric specifying the confidence at which to truncate
#       the output taxonomic classifications. Lower values of threshold
#       will classify deeper into the taxonomic tree at the expense of accuracy,
#       and vice-versa for higher values of threshold.
#
#   DADA2 (default 50):
#       The minimum bootstrap confidence for assigning a taxonomic level.

module load gcc
module load R

# Example with DECIPHER

Rscript scripts/preprocessing/03_add-taxonomy.R \
    data/dada2/02_nochimera_mergeruns/blanes_project/blanes_project_seqtab_final.rds \
    data/dada2/ \
    blanes_project \
    data/assign_tax/SILVA_SSU_r132_March2018.RData \
    decipher \
    60

# Example with dada2 classifier (commented to avoid running it)

# Rscript scripts/preprocessing/03_add-taxonomy.R \
#     data/dada2/02_nochimera_mergeruns/blanes_project/blanes_project_seqtab_final.rds \
#     data/dada2/ \
#     blanes_project \
#     data/assign_tax/your_database_training.fasta,data/assing_tax/your_database_species.fasta \
#     dada \
#     50
