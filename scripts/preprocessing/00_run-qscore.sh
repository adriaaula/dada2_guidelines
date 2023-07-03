#!/bin/sh

#SBATCH --account=<your-account>
#SBATCH --job-name=qprofile
#SBATCH --ntasks=1
#SBATCH --output=data/logs/0_qprofile_%J.out
#SBATCH --error=data/logs/0_qprofile_%J.err

# ~ here you will have to include the chunk necessary for your cluster!~


### DADA2 pipeline ######
## ~~ Q profile ~~ ##

# The script generates a .pdf with the Qscore profiles
# from your first 9 samples, in order to be able to define
# the trimming length of the reads (both forward and reverse).

## ARGS ##

#[1] /input dataset/, with the fastq.gz (the name of the samples
#    has to be present on the left and separated by an underscore.
#    (multiple underscores are OK as long as the name is on the left).
#    Examples: sample120412-g20-mallorca_R1.fastq.gz
#              bl120513_primerx_R1.fastq.gz

#[2] /output dir/ A directory in which all the output will be stored.
#                 If you have copied the github project, you should have
#                 an data dir. in which it should be copied

#[3] /name/ A common identifier to keep track of the output
#           You will thank us that :^)

# If your cluster works with modules, first you
# should activate them.
#module load gcc/4.9.0
#module load Rstats/R-3.4.1

module load R

Rscript scripts/preprocessing/00_qscore.R \
        data/trimmed \
        data/dada2 \
        blanes_project

# IMPORTANT POINT !!!!
# If you want to save the ouptput//errors in a logfile, simply add
# > logout.txt 2> logout_err.txt to have all the info

# As it is written now, you can also work locally. 
# Executing this script locally will make a call to scripts/preprocessing/00_qscore.R
# and output directly the messages in the console
# The same applies to all the other scripts!
