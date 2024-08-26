#!/bin/bash

#### Expected runtime < 2h, expected vmem usage < 5G, expected localhd usage > 30G
umask 0007

#### HPF compute node scratch space
TEMP_DIR=/localhd/scratch/$PBS_JOBID/
mkdir -p $TEMP_DIR
mkdir $TEMP_DIR/result1/
mkdir $TEMP_DIR/result2/

#### Reference data
# Create Salmon quasi mapping index for gencode.v19.protein-coding_transcript and gencode.v19.lncRNA_transcript
# Salmon-0.8.2_linux_x86_64/bin/salmon index -t ../ref/gencode.v19/gencode.v19.pc_transcripts.fa -i ../ref/salmon_reference/gencode.v19.pc_transcripts_index --gencode -k 17
# Salmon-0.8.2_linux_x86_64/bin/salmon index -t ../ref/gencode.v19/gencode.v19.lncRNA_transcripts.fa -i ../ref/salmon_reference/gencode.v19.lncRNA_transcripts_index --gencode -k 17
INDEX_DIR1=$HOME_DIR/ref/salmon_reference/gencode.v19.pc_transcripts_index/
INDEX_DIR2=$HOME_DIR/ref/salmon_reference/gencode.v19.lncRNA_transcripts_index/

# Gene to transcript mapping
GTF_FILE=$HOME_DIR/ref/gencode.v19/gencode.v19.annotation.gtf

# unaligned reads for a patient sample
fastq=$(ls $fastq_dir/* | grep "rnaseq_fastq.tar$")
fastq_fname=$(ls $fastq_dir/ | grep "rnaseq_fastq.tar$")
# obtain patient sample uuid
file_size_mb=$(du -m $fastq | cut -f1)
sample_uuid=$(echo $fastq_fname | cut -f 1-4 -d "-")
OUTPUT=${sample_uuid}_${file_size_mb}

#### Output specification
OUTPUT1=$OUTPUT_DIR/${OUTPUT}_pc_gene_expression.txt
OUTPUT2=$OUTPUT_DIR/${OUTPUT}_lncRNA_gene_expression.txt

#### Main
if [[ ! $QRESUB || ! -f $OUTPUT1 || ! -f $OUTPUT2 ]]; then

    tar xvC $TEMP_DIR -f $fastq
    MATE1=$(ls $TEMP_DIR/*_1.fastq.gz)
    MATE2=$(ls $TEMP_DIR/*_2.fastq.gz)

    # gzip -d $MATE1
    # gzip -d $MATE2

    # # Read1 and Read2 for paired-end fastq
    # MATE1=$(ls $TEMP_DIR/*_1.fastq)
    # MATE2=$(ls $TEMP_DIR/*_2.fastq)

    # Trim Galore: wrapper of Cutadapt and FastQC
    #   https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#version-043
    module load trim_galore/0.4.4

    trim_galore --quality 17 --phred33 --stringency 1 --gzip --paired --length 17 -o $TEMP_DIR $MATE1 $MATE2
    #trim_galore --quality 18 --phred33 --stringency 1 --fastqc --gzip --paired --length 18 -o $TEMP_DIR $MATE1 $MATE2

    MATE1_VALIDATED=$(ls $TEMP_DIR/*_1_val_1.fq.gz)
    MATE2_VALIDATED=$(ls $TEMP_DIR/*_2_val_2.fq.gz)

    if [[ -n "$MATE1_VALIDATED" && -n "$MATE2_VALIDATED" ]]; then
        # Quantify protein-coding transcripts
        $WORKING_DIR/Salmon-0.8.2_linux_x86_64/bin/salmon quant -i $INDEX_DIR1 -l A \
                    -1 $MATE1_VALIDATED -2 $MATE2_VALIDATED -o $TEMP_DIR/result1/  \
                    --seqBias  --gcBias -p $NUM_CORES \
                    -g $GTF_FILE --useVBOpt

        # Quantify lncRNA transcripts
        $WORKING_DIR/Salmon-0.8.2_linux_x86_64/bin/salmon quant -i $INDEX_DIR2 -l A \
                    -1 $MATE1_VALIDATED -2 $MATE2_VALIDATED -o $TEMP_DIR/result2/ \
                    --seqBias  --gcBias -p $NUM_CORES \
                    -g $GTF_FILE --useVBOpt

        # Write to output
        cp $TEMP_DIR/result1/quant.genes.sf $OUTPUT1
        cp $TEMP_DIR/result2/quant.genes.sf $OUTPUT2
    fi
fi

# Clean up compute node scratch space
rm -Rfv $TEMP_DIR
