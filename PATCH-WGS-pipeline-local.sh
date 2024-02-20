#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/%j.out
#SBATCH --job-name=blastn-redo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20000
#SBATCH --output=blastn_2021.out.%j

#Convert bam to fastq
#export bam2fq=/scratch/prj/cb_microbiome/tools/biobambam2/2.0.87-release-20180301132713/x86_64-etch-linux-gnu/bin/bamtofastq  > /dev/null 2>&1

cd /scratch
bam=$(echo /scratch/prj/cb_microbiome/bam_files/*.bam)
echo Running bam2fq
$bam2fq \
    collate=1 \
    exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
    filename=${bam} \
    inputformat=bam \
    F=${bam%%.bam}_F1.fq.gz \
    F2=${bam%%.bam}_F2.fq.gz \
    S=${bam%%.bam}_s.fq.gz \
    0=${bam%%.bam}_0.fq.gz \
    02=${bam%%.bam}_02.fq.gz \
    tryoq=1 \
    gz=1 \
    exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
    level=5


#Fastq quality control & trimming 
echo Running fastqc
mkdir fastqc
fastqc  *_F1.fq.gz --outdir fastqc
fastqc  *_F2.fq.gz --outdir fastqc


for fastq in *_F1.fq.gz 
do
        base=${fastq%%_F1.fq.gz}
        trimmomatic PE \
        ${base}_F1.fq.gz \
        ${base}_F2.fq.gz \
        ${base}_trimmed_F1.fq.gz \
        ${base}_F1_UP.fq.gz \
        ${base}_trimmed_F2.fq.gz \
        ${base}_F2_UP.fq.gz \
        LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28 
done

fq1=*_trimmed_F1.fq.gz
fq2=*_trimmed_F2.fq.gz

#Align to combined reference genome

bwa mem /scratch/prj/cb_microbiome/databases/kraken/databases/refseq/standard/library/human/library.fna ${fq1} ${fq2} -t 10 | samtools sort -o ${sample}_host_aligned.bam - 

#Extract "unmapped"/ non-human reads using samtools flags
samtools view -F 4 ${sample}_host_aligned.bam > ${sample}_host_aligned.bam_unmapped.bam
echo ${sample} Unmapped extracted

#Sort the bam file using samtools
samtools sort ${sample}_unmapped.bam > ${sample}_unmapped_sorted.bam

#Convert bam file to fastq format using bedtools
bedtools bamtofastq -i ${sample}_unmapped_sorted.bam -fq ${sample}_unmapped_F1.fq -fq2 ${sample}_unmapped_F2.fq

echo ${sample} Unmapped bam to fastq

#Denovo assemlbly of unmapped reads using SPAdes
spades -1 ${sample}_unmapped_F1.fq -2 ${sample}_unmapped_F2.fq --only-assembler -o spades

#Run kraken for pathogen classification 

mkdir kraken

kraken2 --use-names  --db /scratch/prj/cb_microbiome/databases/kraken/databases/refseq/standard spades/transcripts.fasta \
--report kraken/kraken_report.txt --classified-out kraken/kraken_classifications.txt >> kraken/output_kraken.txt

grep -e 'pathogen_of_interest' kraken/output_kraken.txt | awk '{print $2}' > kraken/pathogen_nodes.txt

seqtk subseq spades/transcripts.fasta kraken/pathogen_nodes.txt > kraken/pathogen_sequences.fasta

blastn -query kraken_2021/fuso_seq_21.fasta -db pathogens_of_interest_ref_database.fa \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
-max_target_seqs 1 -max_hsps 1 -out kraken/pathogen_gene_annotation.blastn


#Using 
mkdir centrifuge
centrifuge -x /scratch/users/k1802884/azure/radhika/centrifuge/p_compressed+h+v \
-f spades/transcripts.fasta \
--report-file centrifuge/centrifuge_report.txt -S centrifuge/centrifuge_output.txt

grep -e 'pathogen_of_interest' centrifuge/centrifuge_report.txt | awk '{print $3}' | sort | uniq > centrifuge/tax_id_list.txt
awk -F' ' 'NR==FNR{c[$1]++;next};c[$3]' tax_id_list.txt entrifuge/centrifuge_output.txt | awk '{print $1}' | sort -u | uniq > centrifuge/pathogen_nodes.txt
seqtk subseq spades/transcripts.fasta centrifuge/pathogen_nodes.txt > centrifuge/pathogen_sequences.fasta

blastn -query centrifuge/pathogen_sequences.fasta -db /pathogens_of_interest_ref_database.fa \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
	-max_target_seqs 1 -max_hsps 1 -out centrifuge/pathogen_gene_annotation.blastn

#Run BLASTn for pathogen classification 

mkdir blastn

blastn -db /scratch/users/k1802884/azure/radhika/blast_nt/nt -query spades/transcripts.fasta -num_threads 16 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
-max_target_seqs 1 -max_hsps 1 -out blastn/blastn_output.txt

grep -e 'pathogen_of_interest' blastn/blastn_output.txt | awk '{print $1}' | sort -u | uniq > blastn/pathogen_nodes.txt
seqtk subseq spades/transcripts.fasta blastn/pathogen_nodes.txt > blastn/pathogen_reads.fasta

blastn -query blastn/pathogen_reads.fasta -db pathogens_of_interest_ref_database.fa \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
-max_target_seqs 1 -max_hsps 1 -out blastn/pathogen_gene_annotation.blastn

