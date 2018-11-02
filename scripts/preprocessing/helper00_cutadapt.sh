#!/bin/bash

DATA_DIR=data # where your untrimmed files are located
PRIMER_F="TTGTACACACCGCCC" # 1389F, change it for your forward primer (5'-3')
PRIMER_R="CCTTCYGCAGGTTCACCTAC" # 1510R, change it for your reverse primer (5'-3')
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))  # primer match is >= 2/3 of primer length, taken from Fred Mahé's swarm pipeline
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))

# create trimmed and logs directory

mkdir data/trimmed
OUT_DIR=data/trimmed

mkdir -p data/logs/cutadapt
LOG_DIR=data/logs/cutadapt

# load Python for cutadapt

module load python/2.7.13

# cutadapt

### The first line in the following loop takes sample names.
### In this case, assumes that fields in filename are separated by '-' (in the option -F '-')
### and sample name is in the second field (stated by '{print $2}').
### Before runnning the script adapt these 2 options to your samples names.

for SAMPLE in $(ls ${DATA_DIR}/*fastq* | awk -F '-' '{print $2}' | sort -u); do
  cutadapt \
    --discard-untrimmed \
    --pair-filter=any \
    --minimum-length=${MIN_LENGTH} \
    -g ${PRIMER_F} \
    -G ${PRIMER_R} \
    -o ${OUT_DIR}/${SAMPLE}_trimmed_R1.fastq \
    -p ${OUT_DIR}/${SAMPLE}_trimmed_R2.fastq \
    -O ${MIN_R} \
    ${DATA_DIR}/*${SAMPLE}*_R1.fastq* \
    ${DATA_DIR}/*${SAMPLE}*_R2.fastq* \
    > ${LOG_DIR}/${SAMPLE}.log
done
