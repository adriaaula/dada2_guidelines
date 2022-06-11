#!/bin/bash
#SBATCH --account=<your-account>
#SBATCH --job-name=cutadapt
#SBATCH --ntasks=1
#SBATCH --output=data/logs/cutadapt_%J.out
#SBATCH --error=data/logs/cutadapt_%J.err

DATA_DIR=data/raw # where your untrimmed files are located
PRIMER_F="TTGTACACACCGCCC" # 1389F, change it for your forward primer (5'-3')
PRIMER_R="CCTTCYGCAGGTTCACCTAC" # 1510R, change it for your reverse primer (5'-3')
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))  # primer match is >= 2/3 of primer length, taken from Fred Mahé's swarm pipeline
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))

# create trimmed and logs directory

mkdir -p data/trimmed
OUT_DIR=data/trimmed

mkdir -p data/logs/cutadapt
LOG_DIR=data/logs/cutadapt

# load Python for cutadapt

module load cutadapt

# cutadapt

### The first line in the following loop takes sample names.
### In this case, assumes that pair identifiers are '_R1' and '_R2'.
### Change these if this is not your case.
### Alternatively you could create a file with all your samples names and iterate over it (i.e. for SAMPLE in $(cat samples.txt); do...)

for SAMPLE in $(ls ${DATA_DIR}/*fastq* | awk -F"/" '{print $NF}' |  sed 's/_R[1,2].*$//g' | sort -u)
do
  cutadapt \
    --discard-untrimmed \
    --pair-filter=any \
    --minimum-length=${MIN_LENGTH} \
    -g ${PRIMER_F} \
    -G ${PRIMER_R} \
    -o ${OUT_DIR}/${SAMPLE}_trimmed_R1.fastq.gz \
    -p ${OUT_DIR}/${SAMPLE}_trimmed_R2.fastq.gz \
    -O ${MIN_R} \
    ${DATA_DIR}/${SAMPLE}_R1*.fastq* \
    ${DATA_DIR}/${SAMPLE}_R2*.fastq* \
    > ${LOG_DIR}/${SAMPLE}.log
done

## Stats

module load seqkit 
### Additionally, we can extract the stats of the files
# If you don't need to cut the adapter, you can simply copy this last command. 

seqkit stats --all --tabular ${DATA_DIR}/*fastq* > ${DATA_DIR}/seqkit_stats_untrimmed.tsv
seqkit stats --all --tabular ${OUT_DIR}/*fastq* > ${OUT_DIR}/seqkit_stats_trimmed.tsv
