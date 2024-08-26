#/bin/bash

## 	THIS SCRIPT CONVERTS .CEL FILE DATA INTO BIRDSEED CALL FILES
## 			SCRIPT BASED ON APT-PROBESET GENOTYPE MANUAL: http://www.affymetrix.com/estore/support/developer/powertools/changelog/apt-probeset-genotype.html.affx;jsessionid=19CBBF005CB30218DAF14CBECCE0D887

WORKING_DIR=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess
JOB=converting-cel-files
SH_SCRIPT=$JOB.sh

cd $WORKING_DIR

OE=$WORKING_DIR/OE-$JOB/
mkdir -p $OE
rm $OE/*

CANCERS=(LIHC BLCA COAD HNSC LAML LUAD PAAD READ STAD THYM BRCA ESCA KIRC LGG LUSC PCPG SARC TGCT UCEC CESC GBM KIRP OV PRAD SKCM THCA)

for i in {1..25}; do

INPUT_DIR=/hpf/largeprojects/agoldenb/mingjie/rawdata/gdc/CEL/${CANCERS[i]}/
OUTPUT_DIR=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/${CANCERS[i]}/genotype-calls-july-2016

mkdir -p $OUTPUT_DIR

SH_JOB="module load apt; \
apt-probeset-genotype -o $OUTPUT_DIR/genotype-calls-july-2016 \
-c $WORKING_DIR/GenomeWideSNP_6.cdf \
--set-gender-method cn-probe-chrXY-ratio \
--chrX-probes $WORKING_DIR/GenomeWideSNP_6chrX.chrXprobes \
--chrY-probes $WORKING_DIR/GenomeWideSNP_6chrYprobes.chrYprobes \
--special-snps $WORKING_DIR/GenomeWideSNP_6SpecialSNPlist.specialSNPs \
--read-models-birdseed $WORKING_DIR/GenomeWideSNP_6.v2.6.birdseed.models \
-a birdseed-v2 $INPUT_DIR/*/*.CEL"

echo $SH_JOB | qsub -N ${CANCERS[i]}-$JOB \
-l walltime=48:00:00,vmem=30G \
-o $OE -e $OE
done 