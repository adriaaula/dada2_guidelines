#!/bin/sh

#SBATCH --account=<your-account>
#SBATCH --job-name=chimera
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=data/logs/chimera_%J.out
#SBATCH --error=data/logs/chimera_%J.err

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

#[4] /Taxonomy db (general)/ The db has to be downloaded from https://benjjneb.github.io/dada2/training.html
#                             and saved in the biocluster. It can be in .gz, no problem! 

#[5] /Taxonomy db (species)/ The db has to be downloaded from https://benjjneb.github.io/dada2/training.html
#                             and saved in the biocluster. It can be in .gz, no problem!

#[6] /trim length/ After all the processing, some unespecific amplified reads still are present in the samples. 
#                   Time to cut them down. Specify with a range which read you want to keep (Example: 400,450)


module load gcc
module load R

# remember, this is an example: you should change [1,3,6] at least

Rscript scripts/preprocessing/02_chimerarem_merge_tax.R \
                    data/dada2/01_errors-output/blanes_project/blanes_projectseqtab.rds \
                    data/dada2/ \
                    blanes_project \
                    data/assign_tax/rdp_train_set_16.fa.gz \
                    data/assign_tax/rdp_species_assignment_16.fa.gz \
                    400,450 

# If you have multple seqtabs, it should be written like this:

#data/1_errors-output/blanes_project/blanes4020,data/1_errors-output/blanes_project/blanes4300 ...
