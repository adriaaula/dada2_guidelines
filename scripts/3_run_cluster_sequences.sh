#!/bin/bash
    
    clear
    
    date
    

#### DADA2 pipeline ######
## ~~ Clustering at desired threshold~~ ##

#For some analysis, it can be interesting to cluster all 
#the sequences at a desired threshold. This script takes a fasta file
#and generates two documents: A centroid file with all the fastas of the threshold clusters, 
#and a correspondence .tsv to be able to know the ASV/OTU relationships. 

# Usearch path (where does your happy usearch live?)
usearch=/share/apps/usearch/USEARCH_V9.2.64_64bits/usearch9.2.64_i86linux64

#Parameters!
id_job="blanes_project"
input=analysis/2_merge-taxonomy/asv_Blanes_16S.fasta
output=analysis/3_clustering

mkdir output

# This is an example, any threshold can be applied
$usearch -cluster_smallmem ${input} \
         -id 0.97 \
         -centroids ${output}/otus_97.fa \
         -sortedby size \
         -uc ${output}/${id_job}_uc_clusters.uc


module load python/3.6.0

# A small script to create the tsv!
python3 3_asv2otu.py ${input} \
                     ${output}/${id_job}_uc_clusters.uc > $output/${id_job}_correspondence.tsv 
