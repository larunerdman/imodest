#!/bin/bash

# Intro to what to pick
# http://docs.cancergenomicscloud.org/docs/tcga-metadata#section-sample
WORKING_DIR=/hpf/largeprojects/agoldenb/mingjie/clinical_preprocess_input/

cd $WORKING_DIR

CANCERS=(READ-COAD BLCA COAD HNSC LAML LUAD PAAD READ STAD THYM BRCA ESCA KIRC LGG LUSC PCPG SARC TGCT UCEC CESC GBM KIRP LIHC OV PRAD SKCM THCA)


for i in {0..0}; do
	> ${CANCERS[i]}.clin.picked.txt

	cat ${CANCERS[i]}.clin.merged.txt | grep "sample.bcr_sample_barcode" >> ${CANCERS[i]}.clin.picked.txt

	cat ${CANCERS[i]}.clin.merged.txt | grep "sample-2.bcr_sample_barcode" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "sample-3.bcr_sample_barcode" >> ${CANCERS[i]}.clin.picked.txt


	cat ${CANCERS[i]}.clin.merged.txt | grep "admin.batch_number" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "patient.icd_10" >> ${CANCERS[i]}.clin.picked.txt
	# cat ${CANCERS[i]}.clin.merged.txt | grep "patient.icd_o_3_site" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "patient.icd_o_3_histology" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "patient.gender" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "age_at_initial_pathologic_diagnosis" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "patient.histological_type" >> ${CANCERS[i]}.clin.picked.txt

	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample.tumor_nuclei_percent" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample.tumor_locations.tumor_location.site_of_disease\t" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample.tumor_locations.tumor_location.site_of_disease_description\t" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample-2.tumor_nuclei_percent" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample-2.tumor_locations.tumor_location.site_of_disease\t" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample-2.tumor_locations.tumor_location.site_of_disease_description\t" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample-3.tumor_nuclei_percent" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample-3.tumor_locations.tumor_location.site_of_disease\t" >> ${CANCERS[i]}.clin.picked.txt
	cat ${CANCERS[i]}.clin.merged.txt | grep "tumor_sample-3.tumor_locations.tumor_location.site_of_disease_description\t" >> ${CANCERS[i]}.clin.picked.txt



done
