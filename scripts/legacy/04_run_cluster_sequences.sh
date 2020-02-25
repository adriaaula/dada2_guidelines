#!/bin/sh

#SBATCH --account=<your-account>
#SBATCH --job-name=cluster
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=data/logs/cluster_%J.out
#SBATCH --error=data/logs/cluster_%J.err
 
clear
    
date
    
#### DADA2 pipeline ######
## ~~ Clustering at desired threshold~~ ##

#For some analysis, it can be interesting to cluster all 
#the sequences at a desired threshold. This script takes a fasta file
#and generates two documents: A centroid file with all the fastas of the threshold clusters, 
#and a correspondence .tsv to be able to know the ASV/OTU relationships. 

# Usearch module (where does your happy usearch live?)
module load usearch

#Parameters!
id_job="blanes_project"
input=data/dada2/02_merge-taxonomy/asv_Blanes_16S.fasta
output=data/dada2/03_clustering

mkdir ${output}

# This is an example, any threshold can be applied
usearch -cluster_smallmem ${input} \
         -id 0.97 \
         -centroids ${output}/${id_job}_otus_97.fa \
         -sortedby size \
         -uc ${output}/${id_job}_uc_clusters.uc


module load python/3.6.5 # script below is written for Python 3

# A small script to create the tsv!
python3 scripts/preprocessing/03_asv2otu.py ${input} \
                     ${output}/${id_job}_uc_clusters.uc > $output/${id_job}_correspondence.tsv 
