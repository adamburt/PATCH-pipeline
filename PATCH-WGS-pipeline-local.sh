#!/bin/bash -l
if [ $# -eq 0 ]
  then
    echo "No UUID supplied, please supply the UUID:  PATCH-WGS-piepline-local.sh <UUID>"
    exit -1
fi
uuid=$1
bam_files_location="/scratch/bam_files"
output_fastq_files="/scratch/processed_data/fastqc"
trimmomatic_jar_file="/scratch/files/trimmomatic/trimmomatic-0.39.jar"
human_genome_lib="/scratch/files/human_genome/library.fna"
kraken_standard_db="/scratch/files/human_genome/library.fna"
bam_aligned_folder="/scratch/processed_data/bam_aligned"
spades_output_folder="/scratch/processed_data/spades"
kraken_directory="/scratch/processed_data/kraken"

# Convert BAM to FastQ
bam2fq=bamtofastq
cd ${bam_files_location}
bam="${bam_files_location}/${uuid}.bam"
echo "Running bam2fq on file ${uuid}.bam"
$bam2fq \
        collate=1 \
        exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
        filename=${bam} \
        outputdir=${output_fastq_files}/ \
        inputformat=bam \
        F="${output_fastq_files}/${bam%%.bam}_F1.fq.gz" \
        F2="${output_fastq_files}/${bam%%.bam}_F2.fq.gz" \
        S="${output_fastq_files}/${bam%%.bam}_s.fq.gz" \
        tryoq=1 \
        gz=1 \
        level=5

#Fastq quality control & trimming
echo Running fastqc
for file in "${output_fastq_files}/*_F*.fq.gz"
do
    echo Running fastqc on ${file}
    fastqc ${file} -o ${output_fastq_files}
done

# Run Trimmomatic
echo Running Trimmomatic
for fastq in ${output_fastq_files}/*_F1.fq.gz 
do
        base=${fastq%%_F1.fq.gz}
        java -jar ${trimmomatic_jar_file} PE \
        ${base}_F1.fq.gz \
        ${base}_F2.fq.gz \
        ${base}_trimmed_F1.fq.gz \
        ${base}_F1_UP.fq.gz \
        ${base}_trimmed_F2.fq.gz \
        ${base}_F2_UP.fq.gz \
        LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28 
done

#Align to combined reference genome
echo Aligning to combined reference genome
fq1=${output_fastq_files}/*_trimmed_F1.fq.gz
fq2=${output_fastq_files}/*_trimmed_F2.fq.gz
bwa mem ${human_genome_lib} ${fq1} ${fq2} -t 12 -o "${output_fastq_files}/sam_file_temp.sam"


# Sort aligned data
echo Sorting aligned data
samtools sort -o "${output_fastq_files}/${uuid}_host_aligned.bam" "${output_fastq_files}/sam_file_temp.sam"

#Extract "unmapped"/ non-human reads using samtools flags
echo Extracting unmapped data
samtools view -F 4 ${bam_aligned_folder}/host_aligned.bam > ${bam_aligned_folder}/unmapped.bam

#Sort the bam file using samtools
samtools sort ${bam_aligned_folder}/unmapped.bam > ${bam_aligned_folder}/unmapped_sorted.bam

#Convert bam file to fastq format using bedtools
bedtools bamtofastq -i ${bam_aligned_folder}/unmapped_sorted.bam -fq ${bam_aligned_folder}/unmapped_F1.fq -fq2 ${bam_aligned_folder}/unmapped_F2.fq

echo Unmapped bam to fastq

#Denovo assemlbly of unmapped reads using SPAdes
echo Running spades
spades -1 ${bam_aligned_folder}/unmapped_F1.fq -2 ${bam_aligned_folder}/unmapped_F2.fq --only-assembler -o ${spades_output_folder}

#Run kraken for pathogen classification 
echo Running Kraken

kraken2 --use-names  --db ${kraken_standard_db} \
--report ${kraken_directory}/kraken_report.txt --classified-out ${kraken_directory}/kraken_classifications.txt >> ${kraken_directory}/output_kraken.txt

grep -e 'pathogen_of_interest' ${kraken_directory}/output_kraken.txt | awk '{print $2}' > ${kraken_directory}/pathogen_nodes.txt

seqtk subseq ${spades_output_folder}/transcripts.fasta ${kraken_directory}/pathogen_nodes.txt > ${kraken_directory}/pathogen_sequences.fasta


'''
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
'''