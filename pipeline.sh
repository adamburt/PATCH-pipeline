#!/bin/bash -l
send_mail="/home/create/send_mail.sh"
token=$(cat /home/create/token.txt)
bam_files_location="/scratch/bam_files"
output_fastq_files="/scratch/processed_data/fastqc"
trimmomatic_jar_file="/scratch/prj/cb_microbiome/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
human_genome_lib="/scratch/prj/cb_microbiome/databases/kraken/databases/refseq/standard/library/human/library.fna"
kraken_db_out="/scratch/prj/cb_microbiome/databases/kraken/databases/refseq/standard"
kraken_standard_db="/scratch/prj/cb_microbiome/databases/kraken/databases/refseq/standard/library/human/library.fna"
bam_aligned_folder="/scratch/processed_data/bam_aligned"
spades_output_folder="/scratch/processed_data/spades"
kraken_directory="/scratch/processed_data/kraken"
uuid=""
only_download=false
show_help=false
reset=false

# Parse command line arguments
while getopts ":u:dhr" opt; do
  case $opt in
    u)
      uuid="$OPTARG"
      ;;
    d)
      only_download=true
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
    rm ${bam_aligned_folder}/${uuuid}*.* > /dev/null 2>&1
    rm ${output_fastq_files}/${uuid}*.* > /dev/null 2>&1 
    rm ${output_fastq_files}/${uuid}*.* > /dev/null 2>&1
    

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

# Convert BAM to FastQ
bam2fq=bamtofastq
cd ${bam_files_location}
bam="${bam_files_location}/${uuid}.bam"
echo "Running bam2fq on file ${bam}"
$bam2fq \
        collate=1 \
        exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
        filename=${bam} \
        outputdir=${output_fastq_files}/ \
        inputformat=bam \
        F="${output_fastq_files}/${uuid}.bam_F1.fq.gz" \
        F2="${output_fastq_files}/${uuid}.bam_F2.fq.gz" \
        S="${output_fastq_files}/${uuid}.bam_s.fq.gz" \
        tryoq=1 \
        gz=1 \
        level=5

$send_mail "UUID ${uuid} BAM to FASTQ conversion completed"

#Fastq quality control & trimming
echo "
#################################
### Running fastqc - trimming ###
#################################
"
for file in "${output_fastq_files}/${uuid}.bam_F*.fq.gz"
do
    echo Running fastqc on ${file}
    fastqc ${file} -o ${output_fastq_files}
done

java -jar ${trimmomatic_jar_file} PE \
${output_fastq_files}/${uuid}.bam_F1.fq.gz \
${output_fastq_files}/${uuid}.bam_F2.fq.gz \
${output_fastq_files}/${uuid}.bam_trimmed_F1.fq.gz \
${output_fastq_files}/${uuid}.bam_F1_UP.fq.gz \
${output_fastq_files}/${uuid}.bam_trimmed_F2.fq.gz \
${output_fastq_files}/${uuid}.bam_F2_UP.fq.gz \
LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28 

$send_mail "Trimming has been completed for ${uuid}"

echo "
#############################################
### Aligning to combined reference genome ###
#############################################
"
#Align to combined reference genome
fq1=${output_fastq_files}/${uuid}.bam_trimmed_F1.fq.gz
fq2=${output_fastq_files}/${uuid}.bam_trimmed_F2.fq.gz
cd ${bam_aligned_folder}
bwa mem -t 12 -o "${uuid}_mem_temp_unsorted.bam" ${human_genome_lib} ${fn1} ${fn2} 
samtools sort -o "${uuid}_host_aligned.bam" "${bam_aligned_folder}/${uuid}_mem_temp_unsorted.bam"

$send_mail "Alignment has completed for ${uuid}"

echo "
################################
### Extracting unmapped data ###
################################
"
samtools view -F 4 -o "${bam_aligned_folder}/${uuid}_host_aligned_unmapped.bam" "${bam_aligned_folder}/${uuid}_host_aligned.bam"

samtools sort -o "${bam_aligned_folder}/${uuid}_host_aligned_unmapped_sorted.bam" "${bam_aligned_folder}/${uuid}_host_aligned_unmapped.bam"

$send_mail "Unmapped data extraction has completed for ${uuid}"


echo "
##############################
### Unmapping bam to fastq ###
##############################
"
bedtools bamtofastq -i "${bam_aligned_folder}/${uuid}_host_aligned_unmapped_sorted.bam" -fq "${bam_aligned_folder}/${uuid}_host_aligned_unmapped_F1.fq" -fq2 "${bam_aligned_folder}/${uuid}_host_aligned_unmapped_F2.fq"
$send_mail "Unmapping bam to fastq completed for ${uuid}"


echo "
###########################################
### Starting assembly of unmapped reads ###
###########################################
"
spades -1 "${bam_aligned_folder}/${uuid}_host_aligned_unmapped_F1.fq" -2 "${bam_aligned_folder}/${uuid}_host_aligned_unmapped_F2.fq" --only-assembler -o ${spades_output_folder} > /dev/null 2>&1


echo "
###########################################
### Starting classification with kraken ###
###########################################
"
kraken2 --use-names --db ${kraken_db_out} "${spades_output_folder}/contigs.fasta" --report "${kraken_directory}/${uuid}_kraken_report.txt" --classified-out "${kraken_directory}/${uuid}_kraken_classifications.txt" >> "${kraken_directory}/${uuid}_output_kraken.txt"