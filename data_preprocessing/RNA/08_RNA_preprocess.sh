#!/bin/bash

umask 0007

#### Reference files
HISAT2_INDEX=$HOME_DIR/ref/hisat2_reference/grch37_snp_tran/genome_snp_tran
GENE_MAPPING=$HOME_DIR/ref/Homo_sapiens.GRCh37.75/Homo_sapiens.GRCh37.75.gtf

#### Module specification
module load trim_galore/0.4.4
module load hisat/2.1.0
module load subread/1.5.3

#### HPF compute node scratch space
TEMP_DIR=/localhd/scratch/$PBS_JOBID/
mkdir -p $TEMP_DIR


#### uncompress fastq data
fastq=$(ls $fastq_dir/* | grep ".tar.gz$")
tar -xzvf $fastq -C $TEMP_DIR


#### obtain patient sample uuid from file name
sample_uuid=$(echo $fastq | cut -f 2 -d ".")
file_size_mb=$(du -m $fastq | cut -f1)

#### Output specification
OUTPUT=$OUTPUT_DIR/${sample_uuid}_${file_size_mb}_exon_count.txt

if [[ ! $QRESUB || ! -f $OUTPUT ]]; then

  #### Read1 and Read2 for paired-end FASTQ
  MATE1=$(ls $TEMP_DIR/*_1.fastq)
  MATE2=$(ls $TEMP_DIR/*_2.fastq)


  # Trim Galore: wrapper of Cutadapt and FastQC
  #   https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#version-043

  trim_galore --quality 17 --phred33 --stringency 1 --gzip --paired --length 17 -o $TEMP_DIR $MATE1 $MATE2
  #trim_galore --quality 18 --phred33 --stringency 1 --fastqc --gzip --paired --length 18 -o $TEMP_DIR $MATE1 $MATE2

  MATE1_VALIDATED=$(ls $TEMP_DIR/*_1_val_1.fq.gz)
  MATE2_VALIDATED=$(ls $TEMP_DIR/*_2_val_2.fq.gz)


  if [[ -n "$MATE1_VALIDATED" && -n "$MATE2_VALIDATED" ]]; then

      #### Generate aligned reads to genome using hisat2
      #   https://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer
      ALIGNED_READ=$TEMP_DIR/${sample_uuid}_${file_size_mb}_accepted_hits.sam
      hisat2 --seed 1 -q --phred33 -p 8 -x $HISAT2_INDEX -1 $MATE1_VALIDATED -2 $MATE2_VALIDATED -S $ALIGNED_READ


      #### Generate exon-count using featureCounts
      #   http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf
      featureCounts -a $GENE_MAPPING -F GTF \
                    -t exon -g gene_id -f \
                    -Q 20 \
                    -s 0 \
                    -T 8 \
                    -p -B -C \
                    --tmpDir $TEMP_DIR -o $OUTPUT \
                    $ALIGNED_READ

  fi
fi


#### Clean up
rm -Rfv $TEMP_DIR
