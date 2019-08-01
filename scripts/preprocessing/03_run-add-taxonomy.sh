#!/bin/sh

#SBATCH --account=<your-account>
#SBATCH --job-name=tax_decipher
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

#[1] /seqtab/ seqtab with chimeras removed.

#[2] /output dir/ Directory output, usually data (a subdirectory will be created).

#[3] /name/ A common identifier to be sure that the output is the
#            correct one. You will thank us that :^)

#[4] /Taxonomy db/ The db has to be downloaded from http://www2.decipher.codes/Downloads.html
#                             and saved in the biocluster.

#[5] /Confidence level of classification/ Numeric specifying the confidence at which to truncate
#                                         the output taxonomic classifications. Lower values of threshold
#                                         will classify deeper into the taxonomic tree at the expense of accuracy,
#                                         and vice-versa for higher values of threshold.

module load gcc
module load R

# remember, this is an example

Rscript scripts/preprocessing/03_add-taxonomy.R \
                    data/dada2/02_nochimera_mergeruns/blanes_project/blanes_project_seqtab_final.rds \
                    data/dada2/ \
                    blanes_project \
                    data/assign_tax/SILVA_SSU_r132_March2018.RData \
                    60
