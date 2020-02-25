#!/bin/sh

#SBATCH --account=<your-account>
#SBATCH --job-name=clustering
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=data/logs/4_clustering_%J.out
#SBATCH --error=data/logs/4_clustering_%J.err

### DADA2 pipeline ######
###### clustering ######

## ARGS ##

#[1] /seqtab/ seqtab you want to cluster.

#[2] /clustering identity/ clustering identity in a 0-100 scale.

#[3] /output directory/ directory where output files should be written.

#[4] /name/ A common identifier to be sure that the output is the
#            correct one. You will thank us that :^)

module load R

# remember, this is an example

Rscript scripts/preprocessing/04_ASV-clustering.R \
    data/dada2/02_nochimera_mergeruns/blanes_project/blanes_project_seqtab_final.rds \
    97 \
    data/dada2/ \
    blanes_project
