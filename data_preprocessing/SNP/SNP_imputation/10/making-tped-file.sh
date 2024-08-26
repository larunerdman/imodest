#/bin/bash

#PBS -l vmem=25G,walltime=60:00:00
#PBS -o /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/OE-making-tped-file/
#PBS -e /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/OE-making-tped-file/

echo WORKING_DIR: $WORKING_DIR
echo CANCER: $CANCER
echo DATA_DIR: $DATA_DIR

cd $WORKING_DIR

module load python/2.7.9

python $WORKING_DIR/making-tped-file6.py \
		-c $DATA_DIR/birdseed-v2.calls.txt \
		-a $WORKING_DIR/GenomeWideSNP_6.na35.annot.csv \
		-o $DATA_DIR/$CANCER.tped
