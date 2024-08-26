#!/bin/bash

#PBS -l vmem=200G,walltime=48:00:00
#PBS -p 1000

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/ref/

module load hisat/2.1.0
module load python/2.7.9

#### Input Reference files
GENOME_DATA=$WORKING_DIR/gencode.v19/GRCh37.p13.genome.fa
GENE_MAPPING=$WORKING_DIR/gencode.v19/gencode.v19.annotation.gtf


#### Output Reference files
OUTPUT_DIR=$WORKING_DIR/hisat2_reference/
mkdir -p $OUTPUT_DIR
if [ "$(ls -A $OUTPUT_DIR)" ]; then rm $OUTPUT_DIR/*; fi

GENOME_DATA_COPY=$OUTPUT_DIR/GRCh37.p13.genome.1.fa
HISAT2_INDEX=$OUTPUT_DIR/GRCh37.p13.genome



#### correct chromosome name in the genome.fa, such that it match with the chromosome name annotation.gtf file
cut -f 1 -d ' ' $GENOME_DATA > $GENOME_DATA_COPY

#### Create splice site input file and exon input file
SS=$OUTPUT_DIR/splicesites.txt
EXON=$OUTPUT_DIR/exon.txt
python /hpf/tools/centos6/hisat/2.1.0/hisat2_extract_splice_sites.py $GENE_MAPPING > $SS
python /hpf/tools/centos6/hisat/2.1.0/hisat2_extract_exons.py $GENE_MAPPING > $EXON


#### Build hisat2 specific index genome
hisat2-build -p 1 --ss $SS --exon $EXON --seed 1 $GENOME_DATA_COPY $HISAT2_INDEX

