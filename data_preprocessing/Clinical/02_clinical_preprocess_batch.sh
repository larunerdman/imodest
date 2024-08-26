#!/bin/bash

#PBS -l vmem=2G,walltime=23:59:59,gres=localhd:0

umask 0007

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/clinical_preprocess/

#### Script Description
JOB=02_clinical_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R

OE=$WORKING_DIR/OE_$JOB/
if [ -d "$OE" ]; then rm -r $OE; fi
mkdir -p $OE


# Input data
CANCERS=(READ-COAD LIHC BLCA KIRC PAAD SARC THCA BRCA ESCA KIRP LUAD PCPG SKCM CESC GBM LUSC PRAD STAD UCEC COAD HNSC LGG OV READ TGCT)

for i in {0..0}; do

    CANCER=${CANCERS[i]}

    OUTPUT=$HOME_DIR/temp/${CANCER}_sample_ids.rds

    RNA_EXPRESSION_DATA=$HOME_DIR/RNA_preprocess_output/${CANCER}_pc_gene_expression.rds
    CLINICAL_DATA_PICKED=$HOME_DIR/clinical_preprocess_input/${CANCER}.clin.picked.txt
    MIRNA_EXPRESSION_DATA=$HOME_DIR/miRNA_preprocess_output/${CANCER}_miRNA_isoform_expression.rds
    LNCRNA_EXPRESSION_DATA=$HOME_DIR/lncRNA_preprocess_output/${CANCER}_lncRNA_gene_expression.rds
    TF_EXPRESSION_DATA=$HOME_DIR/TF_preprocess_output/${CANCER}_TF_gene_expression.rds
    CNV_MATRIX=$HOME_DIR/cnv_preprocess_output/${CANCER}_cnv_matrix_gene.rds
    METHYLATION_MATRIX=$HOME_DIR/methylation_preprocess_output/${CANCER}_methylation_matrix.rds
    SNP_MATRIX_FOLDER=$HOME_DIR/SNP_preprocess_output/${CANCER}/

    R_JOB=" module load R/3.3.2;

            Rscript --vanilla $R_SCRIPT
                    --CANCER $CANCER
                    --OUTPUT $OUTPUT
                    --RNA_EXPRESSION_DATA $RNA_EXPRESSION_DATA
                    --CLINICAL_DATA_PICKED $CLINICAL_DATA_PICKED
                    --MIRNA_EXPRESSION_DATA $MIRNA_EXPRESSION_DATA
                    --LNCRNA_EXPRESSION_DATA $LNCRNA_EXPRESSION_DATA
                    --TF_EXPRESSION_DATA $TF_EXPRESSION_DATA
                    --METHYLATION_MATRIX $METHYLATION_MATRIX
                    --CNV_MATRIX $CNV_MATRIX
                    --SNP_MATRIX_FOLDER $SNP_MATRIX_FOLDER"

    echo $R_JOB

    echo $R_JOB | qsub -N $CANCER-$JOB \
                       -o $OE -e $OE \
                       -l vmem=10G,walltime=1:00:00 \
                       -p 512
done
