#!/bin/bash

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/methylation_preprocess/

JOB=36_methylation_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R


# Input data
CANCERS=(READ-COAD LIHC BLCA COAD BRCA KIRC LGG HNSC LUAD PRAD SKCM UCEC CESC THCA KIRP LUSC PAAD TGCT PCPG GBM THYM SARC READ LAML STAD ESCA OV STAD ESCA)
PEER_HIDDENS_SAMPLES_RATIOS=(0 0.1 0.2 0.25 0.3 0.4 0.5)



for i in {0..28}; do


    CANCER=${CANCERS[i]}

    OE=$WORKING_DIR/OE_$JOB/$CANCER/
    # if [[ -d "$OE" ]]; then rm -r "$OE"; fi
    mkdir -p "$OE"

    #OUTPUT_DIR=$HOME_DIR/methylation_preprocess_output/peer/$CANCER/
    OUTPUT_DIR=/hpf/largeprojects/agoldenb/mingjie_new/methylation_preprocess_output/peer/$CANCER/
    TEMP_DIR=$OUTPUT_DIR/correct_data_analysis/


    NUM_CORES=1

    for j in {0..0}; do

        PEER_HIDDENS_SAMPLES_RATIO=${PEER_HIDDENS_SAMPLES_RATIOS[j]}
        METHYLATION_MATRIX_CORRECTED_1=$TEMP_DIR/${CANCER}_methylation_matrix_corrected_${j}_1.rds
        METHYLATION_MATRIX_CORRECTED_2=$TEMP_DIR/${CANCER}_methylation_matrix_corrected_${j}_2.rds
        METHYLATION_MATRIX_CORRECTED=$OUTPUT_DIR/${CANCER}_methylation_matrix_corrected_${j}.rds

        if [[ -f $METHYLATION_MATRIX_CORRECTED_1 ]]; then
            R_JOB=" module load R/3.3.2;

                    Rscript --vanilla $R_SCRIPT
                            --CANCER $CANCER
                            --METHYLATION_MATRIX_CORRECTED_1 $METHYLATION_MATRIX_CORRECTED_1
                            --METHYLATION_MATRIX_CORRECTED_2 $METHYLATION_MATRIX_CORRECTED_2
                            --METHYLATION_MATRIX_CORRECTED $METHYLATION_MATRIX_CORRECTED"

            echo $R_JOB   qsub -N $CANCER-$JOB \
                               -l mem=48G,vmem=48G,walltime=1:59:59,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                               -p 100 \
                               -o $OE -e $OE

            echo $R_JOB | qsub -N $CANCER-$JOB \
                               -l mem=48G,vmem=48G,walltime=1:59:59,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                               -p 100 \
                               -o $OE -e $OE
        fi

    done

done
