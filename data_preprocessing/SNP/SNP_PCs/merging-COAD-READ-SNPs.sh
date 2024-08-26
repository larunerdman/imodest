#!/bin/bash

#PBS -l walltime=23:00:00,mem=45g
#PBS -o /hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/COAD-READ/OE-merge
#PBS -e /hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/COAD-READ/OE-merge

working_dir=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/COAD-READ/

setwd $working_dir

module load gtool

coad_data_dir=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/COAD/genotype-calls-july-2016/imputation/
read_data_dir=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/READ/genotype-calls-july-2016/imputation/

gtool -M --g $coad_data_dir/COAD-postqc2-postqc-chr$chr-cut.gen $read_data_dir/READ-postqc2-postqc-chr$chr-cut.gen --s $coad_data_dir/COAD-postqc2-revised-chr$chr.sample $read_data_dir/READ-postqc2-revised-chr$chr.sample --og $working_dir/READ-COAD-postqc2-chr$chr.gen --os $working_dir/READ-COAD-postqc2-chr$chr.sample
