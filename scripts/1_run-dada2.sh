#!/bin/sh

#SBATCH --account=<your-account>
#SBATCH --job-name=dada2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<num-threads>
#SBATCH --output=analysis/logs/dada2_%J.out
#SBATCH --error=analysis/logs/dada2_%J.err

### DADA2 pipeline ######
## ~~ Trimming, error generation and DADA2 run ~~ ##

#Trimming of the reads (and filtering with maxEE), 
#generation of an error model whit the probabilities 
#transitions between bases (A -> T) and dada2 algorithm. 
#BEFORE continue the pipeline, PLEASE, check the error plot
#(a .pdf generated in the output dir). 

## ARGS ##

#[1] /input dataset/, with the fastq.gz (the name of the samples 
#    has to be present on the left and separated by an underscore.
#    (multiple underscores are OK as long as the name is on the left).
#    Example: sample120412-g20-mallorca_R1.fastq.gz 

#[2] /output dir/ Usually it should be in analysis/.
#                 Best option, tidy option

#[3] /name/ A common identifier to be sure that the output is the 
#            correct one. You will thank us that :^)

#[4] /trim position/ The positions to trim, both at the left and the right. The numbers have 
#                     to be separated by a comma. Example: 240 (forward), 210 (rev). 
#                     WATCH OUT. If you trim too close to your read length, you will lose many reads. 
#                     Coded right now for the paired end case. If this is not your case,
#                     say it and we will see what can we do. 

#[5] /max expected error/ The tolerated threshold of maxEE. Too stringent, the reads don't pass. Too
#                         relaxed, the modeled error can be incorrect. Usually, a maxEE of 1,2 if the 
#                         read quality is really good, or a 2,4. 

module load gcc/4.9.0
module load R/R-3.5.0

Rscript scripts/1_dada2-error-output.R \
        data \
        analysis/ \
        blanes_project \
        230,220 \
        2,6

# In case you have multiple runs, simply run again the script with the values of interest. Example:

#Rscript scripts/1_dada2-error-output.R \
#        data/4300 \
#        analysis/ \
#        blanes_run4300 \
#        210,200 \
#        1,4

# In the 3_merge you can aggregate them!

