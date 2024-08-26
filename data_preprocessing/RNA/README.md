# RNAseq data preprocessing
## 05 - Download RNAseq FASTQ data from gdc
Protocol:
1. Build Salmon quasi mapping index based on gencode.v19.protein-coding_transcript and gencode.v19.lncRNA_transcript
[https://www.gencodegenes.org/releases/19.html]

2. Download RNASeq HiSeq paired-end sequencing raw reads (FASTQ) from GDC
[https://portal.gdc.cancer.gov/]

3. Use trim_galore/0.4.4 do quality trimming (Phred quality score threshold of 17) and adapter trimming. Short paired end reads are then removed (threshold 17 bps)
[https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md]

4. Use Salmon/0.82 to quantify gene expression from processed paired end RNAseq data, together sequence-specific bias correction and fragment GC bias correction are performed
[https://combine-lab.github.io/salmon/getting_started/]

## 07 - Use Salmon and do gene expression quantification from RNAseq
### prerequisite: 05
## 08 - Do exon count from RNAseq
### prerequisite: 05
## 20 - Collect gene expression and do data correction


