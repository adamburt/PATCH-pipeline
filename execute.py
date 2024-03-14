from __future__ import annotations
import argparse
import os
import subprocess
import signal
import coloredlogs
import logging
import sys
from tqdm import tqdm
import requests
from beaupy import confirm
from beaupy.spinners import Spinner, DOTS, ARC
import zipfile


class Env():
    def __init__(self, mode: str, uuid: str, base_dir: str, bam_files_dir: str, trimmomatic_jar: str, hisat_db_dir: str, splice_sites_file: str, kraken_db_dir: str, pathogen: str, centrifuge_db_dir: str, token_file: str, start: int = 1, reset: bool = False):
        self.mode: str = mode
        self.uuid: str = uuid
        self.base_dir: str = os.path.join(base_dir, uuid)
        self.start: int = int(start)
        self.bam_files_dir: str = bam_files_dir
        self.fastq_dir: str = os.path.join(self.base_dir, "fastq_output")
        self.bam_aligned_dir: str = os.path.join(self.base_dir, "bam_aligned")
        self.trimmomatic_jar: str = trimmomatic_jar
        self.spades_output_dir: str = os.path.join(
            self.base_dir, "spades_output")
        self.kraken_output_dir: str = os.path.join(
            self.base_dir, "kraken2_output")
        self.blastn_output_dir: str = os.path.join(
            self.base_dir, "blastn_output")
        self.centrifuge_db_dir: str = centrifuge_db_dir
        self.centrifuge_output_dir: str = os.path.join(
            self.base_dir, "centrifuge_output")
        self.hisat_db_dir: str = hisat_db_dir
        self.splice_sites_file: str = splice_sites_file
        self.kraken_db_dir: str = kraken_db_dir
        self.pathogen: str = pathogen
        self.token_file: str = token_file
        self.token: str = ""
        self.reset: bool = reset
        self.__bam_exists__: bool = False

    def validate(self) -> bool:
        success = True

        # Validate we have all options specified
        check_items = [
            {self.bam_files_dir:
                "You must provide a directory that will contain the bam files (--bam-files-dir)"},
            {self.centrifuge_db_dir:
                "You must provide a directory that contains the centrifuge database (--centrifuge-db-dir)"},
            {self.hisat_db_dir:
                "You must provide a directory that contains the Hisat database files (--hisat-db-dir)"},
            {self.kraken_db_dir:
                "You must provide a directory that contains the kraken2 database files (--kraken-db-dir)"},
            {self.splice_sites_file:
                "You must provide the file that contains splice sites information (--splice-sites-file)"},
            {self.token_file:
                "You must provide the file that contains the token for use with downloading BAM files (--token)"}
        ]
        for item in check_items:
            for c, v in item.items():
                if not c:
                    success = False
                    log.error(v)
        if not success:
            return success

        for p in [self.base_dir, self.bam_files_dir, self.fastq_dir, self.bam_aligned_dir, self.spades_output_dir, self.kraken_output_dir, self.blastn_output_dir, self.centrifuge_output_dir]:
            if not os.path.exists(p):
                try:
                    os.makedirs(p, exist_ok=True)
                except Exception as err:
                    success = False
                    log.error(
                        f"The directory ({p}) does not exist")

        if not os.path.exists(self.trimmomatic_jar):
            if confirm(f'Could not find the jar file {self.trimmomatic_jar}, would you like to download it?'):
                os.makedirs(os.path.dirname(
                    self.trimmomatic_jar), exist_ok=True)
                url = "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"
                filepath = f"{self.trimmomatic_jar.rstrip('.jar')}.zip"
                download_file(url, filepath, "Trimmomatic", keep=False)
                log.info(f"Extracting {filepath}")
                with zipfile.ZipFile(filepath, 'r') as zip_ref:
                    file_list = zip_ref.namelist()

                    common_prefix = os.path.commonprefix(file_list)

                    for file in file_list:
                        if file == common_prefix:
                            continue
                        extracted_path = os.path.join(os.path.dirname(
                            self.trimmomatic_jar), file[len(common_prefix):])

                        # Ensure parent directories exist
                        os.makedirs(os.path.dirname(
                            extracted_path), exist_ok=True)

                        # Extract the file
                        try:
                            with zip_ref.open(file) as source, open(extracted_path, 'wb') as dest:
                                dest.write(source.read())
                        except IsADirectoryError:
                            os.makedirs(extracted_path, exist_ok=True)
                os.remove(filepath)
            else:
                success = False
                log.error(
                    f"The jar file ({self.trimmomatic_jar}) does not exist")

        if not os.path.exists(self.kraken_db_dir):
            success = False
            log.error(
                f"The directory ({self.kraken_output_dir}) does not exist")

        if not os.path.exists(self.hisat_db_dir):
            if confirm(f'Could not find the Hisat folder {self.hisat_db_dir}, would you like to create it and download the DB?'):
                os.makedirs(self.hisat_db_dir, exist_ok=True)
                pass
            else:
                success = False
                log.error(
                    f"The directory ({self.hisat_db_dir}) does not exist")

        if not os.path.exists(self.splice_sites_file):
            success = False
            log.warning(
                f"The file ({self.splice_sites_file}) does not exist")

        if not os.path.exists(self.token_file):
            log.warning(
                f"The token file ({self.token_file}) does not exist, it will be required if downloading BAM files")
        else:
            with open(self.token_file, "r") as fp:
                data = fp.read()
                token = data.strip()
                self.token = token

        return success


class Logger():
    def __init__(self, name: str = __name__):
        self.log = logging.getLogger(name)
        self.log.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
        sh = logging.StreamHandler()
        sh.setFormatter(formatter)
        sh.setLevel(logging.INFO)
        self.log.addHandler(sh)
        coloredlogs.install(level='INFO', logger=self.log,
                            fmt='%(asctime)s - %(levelname)s - %(message)s')

    def set_level(self, level: logging.Logger = logging.INFO):
        self.log.setLevel(level)
        coloredlogs.install(level=level, logger=self.log,
                            fmt='%(asctime)s - %(levelname)s - %(message)s')


logger = Logger()
log = logger.log


def signal_handler(signal, frame):
    log.info("\n\nStopping due to CTRL+C")
    sys.exit(1)


def execute(cmd: list, show_output: bool = True) -> tuple:
    first_cmd = cmd[0] if isinstance(
        cmd, list) else cmd if isintance(cmd, str) else ""
    working = Spinner(DOTS, f"Executing {first_cmd}")
    working.start()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    if show_output:
        while proc.poll() == None:
            for stdout_line in iter(proc.stdout.readline, ""):
                print(stdout_line)
    proc.stdout.close()
    return_code = proc.wait()
    output_msg, error_msg = proc.communicate()
    working.stop()
    return return_code, output_msg, error_msg


def download_file(url: str, filepath: str, filename: str, params: dict = {}, headers: dict = {}, keep: bool = True):

    res = requests.get(url, headers=headers, stream=True)
    total_size = int(res.headers.get("content-length", 0))
    block_size = 1024

    with tqdm(total=total_size, unit="B", unit_scale=True, colour="green", desc=f"Downloading {filename}", leave=keep) as progress_bar:
        with open(filepath, "wb") as file:
            for data in res.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)


def single(env: Env):
    pass


def paired(env: Env):

    # Download BAM files
    if not env.start > 1:
        log.info('### Downloading bam files ###')
        if env.__bam_exists__:
            log.warning(
                f"{env.uuid}.bam already exists, overwriting. To avoid this in the future use --start 2")
        if not os.path.exists(env.token_file):
            log.error(f"The token file ({env.token_file}) does not exist")
            sys.exit(-1)
        if not env.token:
            log.error(f"No token could be extracted from {env.token_file}")
            sys.exit(-1)
        url = f"https://api.gdc.cancer.gov/data/{env.uuid}"
        headers = {
            "X-Auth-Token": env.token
        }
        filepath = os.path.join(env.bam_files_dir, f"{env.uuid}.bam")
        download_file(url, filepath, headers=headers)
        log.info('Download complete')
    else:
        log.info("Skipping download")

    # Convert BAM files
    uuid_f1 = f"{env.fastq_dir}/{env.uuid}_F1.fq"
    uuid_f2 = f"{env.fastq_dir}/{env.uuid}_F2.fq"
    if not env.start > 2:
        log.info(f'### Converting {env.uuid}.bam file to FASTA format ###')
        bam_source = os.path.join(env.bam_files_dir, f"{env.uuid}.bam")

        if os.path.exists(os.path.join(env.fastq_dir, f"{env.uuid}_F1.fq")) and os.path.exists(os.path.join(env.fastq_dir, f"{env.uuid}_F2.fq")):
            log.warning(
                f"{env.fastq_dir}/{env.uuid}_F1.fq already exists, overwriting. To avoid this in the future use --start 3")
            log.warning(
                f"{env.fastq_dir}/{env.uuid}_F2.fq already exists, overwriting. To avoid this in the future use --start 3")
        error_code, output, errors = execute([
            "bedtools",
            "bamtofastq",
            "-i",
            bam_source,
            "-fq",
            uuid_f1,
            "-fq2",
            uuid_f2
        ])
        log.info("Conversion complete")
        if errors:
            log.error(errors)

        log.info(f"Conversion complete for {env.uuid}.bam")
    else:
        log.info("Skipping BAM -> FASTA conversion")

    # Check quality
    if not env.start > 3:
        log.info(f"### Checking quality of F*.fq files ###")
        error_code, output, errors = execute([
            "fastqc",
            uuid_f1,
            "--outdir",
            env.fastq_dir
        ])
        if error_code != 0:
            log.error(errors)
        error_code, output, errors = execute([
            "fastqc",
            uuid_f2,
            "--outdir",
            env.fastq_dir
        ])
        if error_code != 0:
            log.error(errors)
        log.info("Quality check complete")
    else:
        log.info("Skipping fastqc quality check")

    #  Trimming
    trimmed_f1 = f"{env.fastq_dir}/{env.uuid}_trimmed_F1.fq"
    trimmed_f1_up = f"{env.fastq_dir}/{env.uuid}_trimmed_F1_UP.fq"
    trimmed_f2 = f"{env.fastq_dir}/{env.uuid}_trimmed_F2.fq"
    trimmed_f2_up = f"{env.fastq_dir}/{env.uuid}_trimmed_F2_UP.fq"
    if not env.start > 4:
        log.info(f"### Trimming files with Trimmomatic ###")
        error_code, output, errors = execute(
            ["java",
             "-jar",
             env.trimmomatic_jar,
             "PE",
             uuid_f1,
             uuid_f2,
             trimmed_f1,
             trimmed_f1_up,
             trimmed_f2,
             trimmed_f2_up,
             "LEADING:28",
             "TRAILING:28",
             "SLIDINGWINDOW:4:28",
             "MINLEN:28"
             ])
        if error_code != 0:
            log.error(errors)
        log.info("Trimming complete")

        log.info("### Checking quality with fastq ###")
        error_code, output_msg, errors = execute([
            "fastqc",
            trimmed_f1,
            "--outdir",
            f"{env.fastq_dir}"
        ])
        if error_code != 0:
            log.error(errors)

        error_code, output_msg, errors = execute([
            "fastqc",
            trimmed_f2,
            "--outdir",
            f"{env.fastq_dir}"
        ])
        if error_code != 0:
            log.error(errors)
        log.info("Quality check complete")
    else:
        log.info("Skipping trimming with Trimmomatic")

    #  Hisat
    trimmed_sam = f"{env.fastq_dir}/{env.uuid}_trimmed_sam"
    host_aligned = f"{env.fastq_dir}/{env.uuid}_host_aligned.bam"
    if not env.start > 5:
        log.info(f"### Running Hisat and alignment ###")
        error_code, output_msg, errors = execute([
            "hisat2",
            "-x",
            f"{env.hisat_db_dir}/genome",
            "--known-splicesite-infile",
            env.splice_sites_file,
            "-1",
            trimmed_f1,
            "-2",
            trimmed_f2,
            "-S",
            trimmed_sam
        ])

        if error_code != 0:
            log.error(errors)
        log.info(f"Hisat completed")
        log.info(f"Executing samtools")

        error_code, output_msg, errors = execute([
            "samtools",
            "view",
            "-bS",
            "-o",
            host_aligned,
            trimmed_sam
        ])
        if error_code != 0:
            log.error(errors)

        log.info('Samtools complete')
    else:
        log.info("Skipping Hisat alignment")

    # Extract unmapped reads
    unmapped_bam = f"{env.fastq_dir}/{env.uuid}_unmapped.bam"
    unmapped_sorted = f"{env.fastq_dir}/{env.uuid}_unmapped_sorted.bam"
    unmapped_f1 = f"{env.fastq_dir}/{env.uuid}_unmapped_F1.fq"
    unmapped_f2 = f"{env.fastq_dir}/{env.uuid}_unmapped_F2.fq"
    if not env.start > 6:

        log.info("### Extracting unmapped / non-human reads ###")
        error_code, output_msg, errors = execute([
            "samtools",
            "view",
            "-F",
            "4",
            "host_aligned",
            "-o",
            unmapped_bam
        ])

        if error_code != 0:
            log.error(errors)
        log.info('Extraction complete')

        log.info("Sorting results")

        error_code, output_msg, errors = execute([
            "samtools",
            "sort",
            unmapped_bam,
            "-o",
            unmapped_sorted
        ])

        if error_code != 0:
            log.error(errors)
        log.info("Sorting complete")

        log.info("Converting results to fastq")

        error_code, output_msg, errors = execute([
            "bedtools",
            "bamtofastq",
            "-i",
            unmapped_sorted,
            "-fq",
            unmapped_f1,
            "-fq2",
            unmapped_f2
        ])

        if error_code != 0:
            log.error(errors)
        log.info('Conversion complete')
    else:
        log.info("Skipping extraction")

    if not env.start > 7:
        log.info("### Running assembly of unmapped reads using SPAdes ###")

        os.mkdirs(os.path.join(
            env.spades_output_dir, env.uuid), exist_ok=True)
        error_code, output_msg, errors = execute([
            "rnaspades",
            "-1",
            unmapped_f1,
            "-2",
            unmapped_f2,
            "--only-assembler",
            "-o",
            f"{env.fastq_dir}/{env.uuid}"
        ])

        if error_code != 0:
            log.error(errors)
        log.info('Assembly complete')
    else:
        log.info("Skipping assembly with SPAdes")

    kraken_report = f"{env.kraken_output_dir}/{env.uuid}_kraken_report.txt"
    kraken_classifications = f"{env.kraken_output_dir}/{env.uuid}_kraken_classifications.txt"
    kraken_output = f"{env.kraken_output_dir}/{env.uuid}_output_kraken.txt"
    spades_transcripts = f"{env.spades_output_dir}/{env.uuid}/transcripts.fasta"
    kraken_pathogen_nodes = f"{env.kraken_output_dir}/{env.uuid}_pathogen_nodes.txt"
    kraken_pathogen_transcripts = f"{env.kraken_output_dir}/{env.uuid}_pathogen_sequences.fasta"
    if not env.start > 8:
        log.info("### Running Kraken2 classification ###")

        error_code, output_msg, errors = execute([
            "kraken2",
            "--use-names",
            "--db",
            kraken_db_dir,
            "--report",
            kraken_report,
            "--classified-out",
            kraken_classifications,
            "--output",
            kraken_output,
            spades_transcripts
        ])

        if error_code != 0:
            log.error(errors)
        log.info('Kraken2 classification complete')
    else:
        log.info("Skipping Kraken2 classification")

    # Pathogen extraction
    if not env.start > 9:
        log.info("### Starting pathogen extraction ###")

        if env.pathogen:
            error_code, output_msg, errors = execute([
                "grep",
                "-e",
                f"'{env.pathogen}'",
                kraken_output,
                "|",
                "awk",
                "'{print $2}'",
                ">",
                kraken_pathogen_nodes
            ])
        else:
            error_code, output_msg, errors = execute([
                "cat",
                kraken_output,
                "|",
                "awk",
                "'{print $2}'",
                ">",
                kraken_pathogen_nodes
            ])

        if error_code != 0:
            log.error(errors)
        log.info('Pathogen extraction complete complete')

        log.info('### Sequencing ###')

        error_code, output_msg, errors = execute([
            "seqtk",
            "subseq",
            spades_transcripts,
            kraken_pathogen_nodes,
            ">",
            kraken_pathogen_transcripts
        ])

        if error_code != 0:
            log.error(errors)
        log.info('Sequencing complete')
    else:
        log.info("Skipping pathogen extraction")

    # Blastn DB creation
    blastn_db = f"{env.blastn_output_dir}/{env.uuid}/blastndb"
    kraken_pathogen_gene = f"{env.kraken_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
    if not env.start > 11:
        log.info("### Creating Blastn DB ###")

        os.mkdirs(os.path.join(env.blastn_output_dir,
                               env.uuid), exists_ok=True)
        error_code, output_msg, errors = execute([
            "makeblastdb",
            "-in",
            kraken_pathogen_transcripts,
            "-dbtype",
            "nucl",
            "-out",
            blastn_db
        ])

        if error_code != 0:
            log.error(errors)
        log.info('DB build complete')
    else:
        log.info("Skipping Blastn DB build")

    # Blastn DB query
    if not env.start > 12:
        log.info("### Starting Blastn query ###")

        error_code, output_msg, errors = execute([
            "blastn",
            "-query",
            kraken_pathogen_transcripts,
            "-db",
            blastn_db,
            "-outfmt",
            "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'",
            "-max_target_seqs",
            "1",
            "-max_hsps",
            "1",
            "-out",
            kraken_pathogen_gene
        ])

        if error_code != 0:
            log.error(errors)
        log.info('DB query complete')
    else:
        log.info("Skipping Blastn query")

    centrifuge_report = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_report.txt",
    tax_id_list = f"{env.centrifuge_output_dir}/{env.uuid}_tax_id_list.txt"
    centrifuge_output = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_output.txt"
    centrifuge_pathogen_nodes = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_nodes.txt"
    centrifuge_seq_fasta = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_sequences.fasta"
    centrifuge_compressed = f"{centrifuge_db_dir}/p_compressed+h+v"
    if not env.start > 13:
        log.info("### Running Centrifuge classification ###")

        error_code, output_msg, errors = execute([
            "centrifuge",
            "-x",
            centrifuge_compressed,
            "-f",
            spades_transcripts,
            "--report-file",
            centrifuge_report,
            "-S",
            centrifuge_output
        ])

        if error_code != 0:
            log.error(errors)

        info.log('### Extracting pathogens ###')

        if env.pathogen:
            error_code, output_msg, errors = execute([
                "grep",
                "-e",
                f"'{env.pathogen}'",
                centrifuge_report,
                "|",
                "awk",
                "'{print $3}'",
                "|",
                "sort",
                "|",
                "uniq",
                ">",
                tax_id_list
            ])
        else:
            error_code, output_msg, errors = execute([
                "cat",
                centrifuge_report,
                "|",
                "awk",
                "'{print $3}'",
                "|",
                "sort",
                "|",
                "uniq",
                ">",
                tax_id_list
            ])
        if error_code != 0:
            log.error(errors)

        log.info('Pathogen extraction complete complete')
        log.info('### Running seqtk and blastn ###')

        error_code, output_msg, errors = execute([
            "awk",
            "-F' '",
            "'NR==FNR{c[$1]++;next};c[$3]'",
            tax_id_list,
            centrifuge_output,
            "|",
            "awk",
            "'{print $1}'",
            "|",
            "sort",
            "-u",
            "|",
            "uniq",
            ">",
            centrifuge_pathogen_nodes
        ])
        if error_code != 0:
            log.error(errors)
        error_code, output_msg, errors = execute([
            "seqtk",
            "subseq",
            spades_transcripts,
            centrifuge_pathogen_nodes,
            ">",
            centrifuge_seq_fasta
        ])
        if error_code != 0:
            log.error(errors)
        error_code, output_msg, errors = execute([
            "blastn",
            "-query",
            centrifuge_seq_fasta,
            "-db",
            f"{env.blastn_db_dir}/{env.uuid}_blastndb",
            "-outfmt",
            "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'",
            "-max_target_seqs",
            "1",
            "-max_hsps",
            "1",
            "-out",
            f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
        ])

        if error_code != 0:
            log.error(errors)
        info.log('Seqtk and blastn complete')
        log.info('Centrifuge classification complete')

    else:
        log.info("Skipping Centrifuge classification")

    if not env.start > 14:
        pass


def pipeline(env: Env):
    pass


def reset_progress(env: Env):
    if confirm("Would you like to remove the downloaded bam files?"):
        try:
            os.remove(os.path.join(env.bam_files_dir, f"{env.uuid}.bam"))
        except FileNotFoundError:
            pass
        except Exception as err:
            log.error(err)
    if confirm(f"Would you like to remove the processed files (this will delete the folder '{env.base_dir}')?"):
        error_code, output, errors = execute(
            ["rm", "-rvf", f"{env.base_dir}"], show_output=False)
    log.info(f"Complete")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reset', action="store_true",
                        help="Resets progress (except dowloading BAM file) for the UUID specified.")
    parser.add_argument('mode', metavar="MODE", choices=[
                        'single', 'paired', 'pipeline'], action='store', help="The mode to run, either single, paired or pipeline.")
    parser.add_argument('--uuid', action='store',
                        required=True, help='The UUID to process.')
    parser.add_argument('--base_dir', action='store',
                        required=True, help='The base directory under which all processed files are stored.')
    parser.add_argument('--token', metavar="FILE", action='store',
                        help='The file containing the API Key to download BAM files.')
    parser.add_argument('--start', metavar="LEVEL",
                        choices=["1", "2", "3", "4", "5", "6"], action="store", help="Start at a later stage of processing:\n\n1=Beginning, 2=Skip downloading, 3=Skip conversion to fastq (bedtools), 4=Skip quality checking (fastqc), 5=Skip trimming (Trimmomatic), 6=Skip alignment (Hisat2), 7=Skip extraction (samtools and bedtools), 8=Skip assembly (SPAdes), 9=Skiup Kraken2 classification, 10=Skip pathogen extraction, 11=Skip Blastn db creation, 12=Skip Blastn DB query, 13=Skip Blastn query DB, 14=Skip Centrifuge classification, 15=Skip Blastn classification")
    parser.add_argument('--bam-files-dir', metavar="DIRECTORY",
                        required=True, action='store', help="Specify the bam files directory.")
    parser.add_argument('--trimmomatic-jar', metavar="FILE",
                        action='store', help="Specify the Trimmomatic jar file.")
    parser.add_argument('--hisat-db-dir', metavar="FILE",
                        action='store', help="Specify the Hisat DB directory.")
    parser.add_argument('--splice-sites-file', metavar="FILE",
                        action='store', help="Specify the splice-sites tsv file.")
    parser.add_argument('--kraken-db-dir', metavar="FILE",
                        action='store', help="Specify the Kraken2 DB directory.")
    parser.add_argument('--centrifuge-db-dir', metavar="FILE",
                        action='store', help="Specify the Centrifuge DB directory.")
    parser.add_argument('--pathogen', metavar="NAME",
                        action="store", help="Specify the pathogen of interest.")
    parser.add_argument('--debug', action="store_true",
                        help="Enable debug logging")

    try:
        args = parser.parse_args()
    except Exception as err:
        log.error(err)
        sys.exit(-1)
    return args


def main():
    args = parse_arguments()
    if args.debug:
        logger.set_level(logging.DEBUG)
        log.debug("Debug logging enabled")

    env = Env(
        args.mode,
        args.uuid,
        args.base_dir,
        args.bam_files_dir,
        args.trimmomatic_jar,
        args.hisat_db_dir,
        args.splice_sites_file,
        args.kraken_db_dir,
        args.pathogen,
        args.centrifuge_db_dir,
        args.token,
        start=args.start if args.start else 1,
        reset=args.reset if args.reset else False
    )

    if env.reset:
        log.warning(f"Resetting all progress for UUID {env.uuid}")
        reset_progress(env)
        sys.exit(0)

    if not env.validate():
        sys.exit(-1)

    log.debug(f"BAM_FILES_DIR={env.bam_files_dir}")
    log.debug(f"FASTQC_DIR={env.fastq_dir}")
    log.debug(f"TRIMMOMATIC_JAR_FILE={env.trimmomatic_jar}")

    log.info(f"##### Running in {env.mode} mode #####")
    if env.mode == 'single':
        single(env)
    elif env.mode == 'paired':
        paired(env)
    elif env.mode == 'pipeline':
        piepline(env)


if __name__ in ['__main__', 'builtin', '__builtins__']:
    main()
