#!/bin/bash

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/cnv_preprocess/

JOB=10_cnv_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R

# OE
OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
#if [ "$(ls -A $OE)" ]; then rm $OE/*; fi

#CANCERS=(BLCA COAD HNSC LAML LUAD PAAD READ STAD THYM BRCA ESCA KIRC LGG LUSC PCPG SARC TGCT UCEC CESC GBM KIRP LIHC OV PRAD SKCM THCA)
CANCERS=(READ-COAD BRCA OV)
GENE_GTF=$HOME_DIR/ref/gencode.v19/gencode.v19.annotation.gtf

for i in {0..0}; do

    # Input
    CANCER=${CANCERS[i]}
    CNV_DATA=$HOME_DIR/cnv_preprocess_input/$CANCER.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt

    # create temp folder for each cancer
    TEMP_FOLDER=$WORKING_DIR/temp/$CANCER/
    mkdir -p $TEMP_FOLDER
    if [ "$(ls -A $TEMP_FOLDER)" ]; then rm $TEMP_FOLDER/*; fi

    for j in {0..1}; do

        EXON_SPECIFIC=$j

        if [[ $EXON_SPECIFIC -eq 1 ]]; then
            CNV_MATRIX=$HOME_DIR/cnv_preprocess_output/${CANCER}_cnv_matrix_exon.rds
        else
            CNV_MATRIX=$HOME_DIR/cnv_preprocess_output/${CANCER}_cnv_matrix_gene.rds
        fi

        SH_JOB="module load R/3.3.2;
                Rscript --vanilla $R_SCRIPT
                --CNV_DATA $CNV_DATA
                --GENE_GTF $GENE_GTF
                --TEMP_FOLDER $TEMP_FOLDER
                --CNV_MATRIX $CNV_MATRIX
                --EXON_SPECIFIC $EXON_SPECIFIC"

        echo $SH_JOB | qsub -N $CANCER-$JOB-$EXON_SPECIFIC \
                            -l vmem=80G,walltime=23:59:59 \
                            -e $OE -o $OE
    done
done
