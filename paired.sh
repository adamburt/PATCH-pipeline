#!/bin/bash -l
send_mail="/home/create/send_mail.sh"
token=$(cat /home/create/token.txt)
bam_files_location="/scratch/bam_files"
output_fastq_files="/scratch/processed_data/fastqc"
trimmomatic_jar_file="/scratch/prj/cb_microbiome/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
kraken_standard_db="/scratch/prj/cb_microbiome/databases/kraken/databases/refseq/standard/library/human/library.fna"
bam_aligned_folder="/scratch/processed_data/bam_aligned"
spades_output_folder="/scratch/processed_data/spades"
kraken2_db_dir="/scratch/databases/kraken2_db"
kraken_directory="/scratch/processed_data/kraken"
centrifuge_dir="/scratch/processsed_data/centrifuge"
centrifuge_db_dir="/scratch/databases/centrifuge"
#splice_sites_file="/scratch/databases/GRCh38/splice_sites/splice_sites.tsv"
splice_sites_file="/scratch/databases/hisat/splicesites.tsv"
#grch38_index_folder="/scratch/databases/GRCh38/index"
grch38_index_folder="/scratch/databases/hisat/index"
blastn_dir="/scratch/processed_data/blastn"
uuid=""
only_download=false
show_help=false
reset=false

# Parse command line arguments
while getopts ":u:dhrp:" opt; do
  case $opt in
    u)
      uuid="$OPTARG"
      ;;
    d)
      only_download=true
      ;;
    p)
      pathogen_of_interest="$OPTARG"
      ;;
    h)
      show_help=true
      ;;
    r)
      reset=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if -h argument is present
if $show_help; then
  echo "Arguments:
    -u <UUID>                           required
    -d download the file then exit      optional
    -r reset the progress for UUID      optional
"
  exit 1
fi

# Check if -u argument is present
if [ -z "$uuid" ]; then
  echo "Error: -u argument is required. Please specify a UUID" >&2
  exit 1
fi


if $reset; then
    echo "Resetting progress for UUID ${uuid}"
    rm ${bam_aligned_folder}/${uuuid}* > /dev/null 2>&1
    rm ${output_fastq_files}/${uuid}* > /dev/null 2>&1 
    rm ${spades_output_folder}/${uuid}* > /dev/null 2>&1
    rm ${kraken_directory}/${uuid}* > /dev/null 2>&1
    rm ${blastn_dir}/${uuid}* > /dev/null 2>&1
    rm -rvf ${blastn_dir}/${uuid} > /dev/null 2>&1
    rm -rvf ${spades_output_folder}/${uuid} > /dev/null 2>&1

    exit 0
fi

if [ ! -f "${bam_files_location}/${uuid}.bam" ]; then
    echo "The file ${uuid}.bam does not exist."
    echo "
#############################
### DOWNLOADING BAM FILES ###
#############################
"

    url='https://api.gdc.cancer.gov/data/'$uuid
    echo 'Downloading file with UUID '$uuid':'
    curl -J -H "X-Auth-Token: $token" -o "${bam_files_location}/${uuid}.bam" $url
    echo 'Download complete'

    chmod 775 "${bam_files_location}/${uuid}.bam"

    $send_mail "The file with UUID ${uuid} was downloaded successfully."

    if $only_download; then
        exit 0
    fi
else
    if $only_download; then
        echo "File ${uuid} is already downloaded"
        exit 0
    fi
fi



echo "
##############################
### COVERTING BAM -> FASTQ ###
##############################
"

if ! [ -f ${output_fastq_files}/${uuid}_F1.fq ]; then
  if ! [ -f ${output_fastq_files}/${uuid}_F2.fq ]; then
    echo "Both ${output_fastq_files}/${uuid}_F1.fq} and ${output_fastq_files}/${uuid}_F2.fq already exist"
  fi
fi


# Convert BAM to FastQ
echo "Running bam2fq on file ${uuid}.bam"
bedtools bamtofastq -i ${bam_files_location}/${uuid}.bam -fq ${output_fastq_files}/${uuid}_F1.fq -fq2 ${output_fastq_files}/${uuid}_F2.fq > /dev/null 2>&1

$send_mail "UUID ${uuid} BAM to FASTQ conversion completed"

#Check quality of fastq file using fastqc
echo "Checking quality of ${uuid}.fq"
fastqc ${output_fastq_files}/${uuid}_F1.fq --outdir ${output_fastq_files}
fastqc ${output_fastq_files}/${uuid}_F2.fq --outdir ${output_fastq_files}

cd ${output_fastq_files}


#Fastq quality control & trimming
echo "
#################################
### Running fastqc - trimming ###
#################################
"

java -jar ${trimmomatic_jar_file} PE \
${output_fastq_files}/${uuid}_F1.fq \
${output_fastq_files}/${uuid}_F2.fq \
${output_fastq_files}/${uuid}_trimmed_F1.fq \
${output_fastq_files}/${uuid}_trimmed_F1_UP.fq \
${output_fastq_files}/${uuid}_trimmed_F2.fq \
${output_fastq_files}/${uuid}_trimmed_F2_UP.fq \
LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28

fastqc ${output_fastq_files}/${uuid}_trimmed_F1.fq --outdir ${output_fastq_files}
fastqc ${output_fastq_files}/${uuid}_trimmed_F2.fq --outdir ${output_fastq_files}

$send_mail "Trimming has been completed for ${uuid}"

echo "
#################################
### Hisat alignment beginning ###
#################################
"
#Align to combined reference genome
hisat2 -x ${grch38_index_folder}/genome --known-splicesite-infile ${splice_sites_file} -1 ${output_fastq_files}/${uuid}_trimmed_F1.fq -2 ${output_fastq_files}/${uuid}_trimmed_F2.fq -S ${output_fastq_files}/${uuid}_trimmed_sam
samtools view -bS -o ${output_fastq_files}/${uuid}_host_aligned.bam ${output_fastq_files}/${uuid}_trimmed_sam
$send_mail "Hisat alignment has been completed for ${uuid}"

echo "
##############################################################
### Extract "unmapped" / non-human reads, sort and convert ###
##############################################################
"
#Extract unmapped/ non-human reads using samtools flags
samtools view -F 4 ${output_fastq_files}/${uuid}_host_aligned.bam -o ${output_fastq_files}/${uuid}_unmapped.bam
samtools sort ${output_fastq_files}/${uuid}_unmapped.bam -o ${output_fastq_files}/${uuid}_unmapped_sorted.bam
bedtools bamtofastq -i ${output_fastq_files}/${uuid}_unmapped_sorted.bam -fq ${output_fastq_files}/${uuid}_unmapped_F1.fq -fq2 ${output_fastq_files}/${uuid}_unmapped_F2.fq > /dev/null 2>&1
$send_mail "Unmapped extraction, sorting and conversion completed for ${uuid}"


echo "
#########################
### Spades processing ###
#########################
"
#Denovo assembly of unmapped reads using SPAdes
mkdir ${spades_output_folder}/${uuid}
rnaspades -1 ${output_fastq_files}/${uuid}_unmapped_F1.fq -2 ${output_fastq_files}/${uuid}_unmapped_F2.fq --only-assembler -o ${spades_output_folder}/${uuid}

echo "
#####################################
### Running Kraken classification ###
#####################################
"
kraken2 --use-names --db ${kraken2_db_dir} --report ${kraken_directory}/${uuid}_kraken_report.txt --classified-out ${kraken_directory}/${uuid}_kraken_classifications.txt --output ${kraken_directory}/${uuid}_output_kraken.txt ${spades_output_folder}/${uuid}/transcripts.fasta

echo "
########################################
### Extracting pathogens of interest ###
########################################
"

# Check if -p argument is present
if [[ ${#pathogen_of_interest} -gt 0 ]] ; then
  echo "Pathogen = '${pathogen_of_interest}'"
  grep -e '${pathogen_of_interest}' ${kraken_directory}/${uuid}_output_kraken.txt | awk '{print $2}' > ${kraken_directory}/${uuid}_pathogen_nodes.txt
else
  echo "Pathogen = ALL"
   
fi

seqtk subseq ${spades_output_folder}/${uuid}/transcripts.fasta ${kraken_directory}/${uuid}_pathogen_nodes.txt > ${kraken_directory}/${uuid}_pathogen_sequences.fasta

echo "
################################
### Building BLASTN database ###
################################
"
mkdir ${blastn_dir}/${uuid}
makeblastdb -in ${kraken_directory}/${uuid}_pathogen_sequences.fasta -dbtype nucl -out ${blastn_dir}/${uuid}/blastndb


blastn -query ${kraken_directory}/${uuid}_pathogen_sequences.fasta -db ${blastn_dir}/${uuid}/blastndb -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out ${kraken_directory}/${uuid}_pathogen_gene_annotation.blastn

echo "
########################################
### Running Centrifuge classification ##
########################################
"
centrifuge -x ${centrifuge_db_dir}/p_compressed+h+v -f ${spades_output_folder}/${uuid}/transcripts.fasta --report-file ${centrifuge_dir}/${uuid}_centrifuge_report.txt -S ${centrifuge_dir}/${uuid}_centrifuge_output.txt
if [[ ${#pathogen_of_interest} -gt 0 ]] ; then
  grep -e '${pathogen_of_interest}' ${centrifuge_dir}/${uuid}_centrifuge_report.txt | awk '{print $3}' | sort | uniq > ${centrifuge_dir}/${uuid}_tax_id_list.txt
else
  cat ${centrifuge_dir}/${uuid}_centrifuge_report.txt | awk '{print $3}' | sort | uniq > ${centrifuge_dir}/${uuid}_tax_id_list.txt
fi

awk -F' ' 'NR==FNR{c[$1]++;next};c[$3]' ${centrifuge_dir}/${uuid}_tax_id_list.txt ${centrifuge_dir}/${uuid}_centrifuge_output.txt | awk '{print $1}' | sort -u | uniq > ${centrifuge_dir}/${uuid}_pathogen_nodes.txt
seqtk subseq ${spades_output_folder}/${uuid}/transcripts.fasta ${centrifuge_dir}/${uuid}_pathogen_nodes.txt > ${centrifuge_dir}/${uuid}_pathogen_sequences.fasta
blastn -query ${centrifuge_dir}/${uuid}_pathogen_sequences.fasta -db ${blastn_dir}/${uuid}_blastndb -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out ${centrifuge_dir}/${uuid}_pathogen_gene_annotation.blastn



echo "
####################################
### Running BLASTN classification ##
####################################
"

blastn --db ${kraken2_db_dir} -query ${spades_output_folder}/${uuid}_transcripts.fasta -num_threads 6 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out ${blastn_dir}/${uuid}_blastn_output.txt
if [[ ${#pathogen_of_interest} -gt 0 ]] ; then
  grep -e '${pathogen_of_interest}' ${blastn_dir}/${uuid}_blastn_output.txt | awk '{print $1}' | sort -u | uniq > ${blastn_dir}/${uuid}_pathogen_nodes.txt
else
  cat ${blastn_dir}/${uuid}_blastn_output.txt | awk '{print $1}' | sort -u | uniq > ${blastn_dir}/${uuid}_pathogen_nodes.txt
fi
seqtk subseq ${spades_output_folder}/${uuid}_transcripts.fasta ${blastn_dir}/${uuid}_pathogen_nodes.txt > ${blastn_dir}/${uuid}_pathogen_reads.fasta

blastn -query ${blastn_output}/${uuid}_pathogen_reads.fasta --db ${kraken2_db_dir} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out ${blastn_output}/${uuid}_pathogen_gene_annotation.blastn
