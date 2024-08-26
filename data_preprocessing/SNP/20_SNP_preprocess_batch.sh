#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

umask 0007

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/projects/MoR/
WORKING_DIR=$HOME_DIR/scripts/SNP_preprocess/

JOB=20_SNP_preprocess
R_SCRIPT=$WORKING_DIR/${JOB}.R

# Some constants
NUM_NUCLEOTIDES=1000000

# Reference
GENE_GTF=$HOME_DIR/ref/gencode.v19/gencode.v19.annotation.gtf

# separate model
#DATASETS=(READ-COAD LIHC BLCA KIRC PAAD SARC THCA BRCA ESCA KIRP LUAD PCPG SKCM CESC GBM LUSC PRAD STAD UCEC COAD HNSC LGG OV READ TGCT)
# second run Sep 10
DATASETS=(ESCA LUSC KIRC)
# joint model
# DATASETS=(AVG)

for i in {0..2}; do # todo: chagne when need

        DATASET=${DATASETS[i]}

        OE=$WORKING_DIR/OE_$JOB/$DATASET/
        if [[ -d "$OE" ]]; then rm -r $OE; fi
        mkdir -p $OE

        TEMP_DIR=$WORKING_DIR/temp/$DATASET/
        if [[ -d "$TEMP_DIR" ]]; then rm -r $TEMP_DIR; fi
        mkdir -p $TEMP_DIR

        # Input folder (imputed SNP)
        # SNP_INPUT_DIR=$HOME_DIR/SNP_preprocess_input/$DATASET/genotype-calls-july-2016/imputation/
        SNP_INPUT_DIR=$HOME_DIR/data/individual_cohort/SNP_preprocess_input/$DATASET/genotype-calls-july-2016/imputation/
        # SNP_INPUT_DIR=$HOME_DIR/rawdata/gdc/merged-cancers/imputation/

        # Output data
        BIM_FILE=$TEMP_DIR/$DATASET-postqc2-postqc.bim

        # Copy and combine over to temp folder, need for later step
        > $BIM_FILE

        for j in {1..22}; do
            # cat $SNP_INPUT_DIR/$DATASET-postqc2-postqc-chr$j.bim >> $BIM_FILE
            cat $SNP_INPUT_DIR/$DATASET-cancer-postqc-postqc2-postqc-somup-chr$j.bim >> $BIM_FILE
            #cat $SNP_INPUT_DIR/all-TCGA-mj-sub-postqc2-postqc-chr$j.bim >> $BIM_FILE
        done

        # Generate SNP target matrix for each chromosome
        for j in {1..22}; do

            c=$j

            for k in {0..0}; do # EXON specific 1; gene specific 0
                EXON_SPECIFIC=$k

                R_JOB=" module load R/3.3.2;

                        Rscript --vanilla $R_SCRIPT
                                --NUM_NUCLEOTIDES $NUM_NUCLEOTIDES
                                --c $c
                                --EXON_SPECIFIC $EXON_SPECIFIC
                                --GENE_GTF $GENE_GTF
                                --BIM_FILE $BIM_FILE
                                --TEMP_DIR $TEMP_DIR"

                echo $R_JOB  qsub -N ${JOB}-${EXON_SPECIFIC} \
                                   -l mem=48G,vmem=48G,walltime=23:59:59 \
                                   -o $OE -e $OE

                echo $R_JOB | qsub -N ${JOB}-${EXON_SPECIFIC} \
                                   -l mem=48G,vmem=48G,walltime=23:59:59 \
                                   -o $OE -e $OE
            done
        done
done


