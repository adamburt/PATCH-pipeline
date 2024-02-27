#!/bin/bash -l


bam2fq=bamtofastq
cd /scratch

#Fastq quality control & trimming
echo Running fastqc
mkdir /scratch/processed_data/fastqc
fastqc  /scratch/processed_data/bam_processed/*_F1.fq.gz --outdir /scratch/processed_data/fastqc
fastqc  /scratch/processed_data/bam_processed/*_F2.fq.gz --outdir /scratch/processed_data/fastqc


for fastq in /scratch/processed_data/fastqc/*_F1.fq.gz 
do
        base=${fastq%%_F1.fq.gz}
        java -jar /scratch/prj/cb_microbiome/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        ${base}_F1.fq.gz \
        ${base}_F2.fq.gz \
        ${base}_trimmed_F1.fq.gz \
        ${base}_F1_UP.fq.gz \
        ${base}_trimmed_F2.fq.gz \
        ${base}_F2_UP.fq.gz \
        LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28 
done

