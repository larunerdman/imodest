#!/bin/bash

#PBS -l vmem=40G,walltime=20:00:00

module load impute2/2.3.1

(>&2 echo WORKING_DIR: $WORKING_DIR
>&2 echo LINE: $LINE
>&2 echo chr: $chr
>&2 echo INPUT: $INPUT
>&2 echo DATA_DIR: $DATA_DIR)


cd $DATA_DIR

impute2 -use_prephased_g -known_haps_g $INPUT-chr$chr.haps.gz \
		-sample_known_haps_g $INPUT-chr$chr.sample \
		-h /hpf/projects/arnold/Reference-files-phase3/1000GP_Phase3_chr${chr}.hap.gz \
		-m /hpf/projects/arnold/Reference-files-phase3/genetic_map_chr${chr}_combined_b37.txt \
		-l /hpf/projects/arnold/Reference-files-phase3/1000GP_Phase3_chr${chr}.legend.gz -allow_large_regions \
		-int $LINE -strand_g $WORKING_DIR/strand-files/chr${chr}.strand \
		-o $INPUT-chr$chr-int$LINE 


 
