#!/bin/sh

#SBATCH --account=<your-account>
#SBATCH --job-name=merge_nochimera
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=data/logs/2_merge_nochimera_%J.out
#SBATCH --error=data/logs/2_merge_nochimera_%J.err

### DADA2 pipeline ######
## ~~ Trimming, error generation and DADA2 run ~~ ##

#Trimming of the reads (and filtering with maxEE),
#generation of an error model whit the probabilities
#transitions between bases (A -> T) and dada2 algorithm.
#BEFORE continue the pipeline, PLEASE, check the error plot
#(a .pdf generated in the output dir).

## ARGS ##

#[1] /seqtables/ each of them separated by a comma.
#                The path should be the complete one. With many seqtables,
#                first they are merged and then bimeras are tested.

#[2] /output dir/ Directory output, usually data (a subdirectory will be created).

#[3] /name/ A common identifier to be sure that the output is the
#            correct one. You will thank us that :^)

#[4] /trim length/ After all the processing, some unespecific amplified reads still are present in the samples.
#                   Time to cut them down. Specify with a range which read you want to keep (Example: 400,450)

#[5] /chimera removal method/ One of 'consensus' (default), 'pooled' or 'per-sample'
#                             If you used pooling in dada inference step you should use 'pooled' method

module load gcc
module load R

# remember, this is an example: you should change [1,3,4] at least

Rscript scripts/preprocessing/02_chimerarem_merge.R \
                    data/dada2/01_errors-output/blanes_project/blanes_project_seqtab.rds \
                    data/dada2/ \
                    blanes_project \
                    400,450 \
                    consensus

# If you have multple seqtabs, it should be written like this:

#data/1_errors-output/blanes_project/blanes4020,data/1_errors-output/blanes_project/blanes4300 ...
