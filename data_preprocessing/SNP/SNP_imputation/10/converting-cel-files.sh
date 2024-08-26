#!/bin/bash

#PBS -l walltime=23:00:00,vmem=10g
#PBS -o /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/OE-converting-cel-files
#PBS -e /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/OE-converting-cel-files

##
## 	THIS SCRIPT CONVERTS .CEL FILE DATA INTO BIRDSEED CALL FILES
## 			SCRIPT BASED ON APT-PROBESET GENOTYPE MANUAL: http://www.affymetrix.com/estore/support/developer/powertools/changelog/apt-probeset-genotype.html.affx;jsessionid=19CBBF005CB30218DAF14CBECCE0D887
## 


echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR

cd $DATA_DIR
mkdir -p ../genotype-calls-july-2016

module load apt


apt-probeset-genotype -o ../genotype-calls-july-2016 \
					  -c $WORKING_DIR/GenomeWideSNP_6.cdf \
					  --set-gender-method cn-probe-chrXY-ratio \
					  --chrX-probes $WORKING_DIR/GenomeWideSNP_6chrX.chrXprobes \
					  --chrY-probes $WORKING_DIR/GenomeWideSNP_6chrYprobes.chrYprobes \
					  --special-snps $WORKING_DIR/GenomeWideSNP_6SpecialSNPlist.specialSNPs \
					  --read-models-birdseed $WORKING_DIR/GenomeWideSNP_6.v2.6.birdseed.models \
					  -a birdseed-v2 *.CEL
