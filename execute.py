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
from beaupy import confirm, prompt, console, select
from beaupy.spinners import Spinner, DOTS, ARC
import zipfile
import gzip
import tarfile


class Env():
    def __init__(self, mode: str, uuid: str, base_dir: str, bam_files_dir: str, trimmomatic_jar: str, hisat_db_dir: str, splice_sites_file: str, kraken_db_dir: str, pathogen: str, centrifuge_db_dir: str, token_file: str, start: int = 1, reset: bool = False, skip_quality: bool = False):
        self.mode: str = mode
        self.uuid: str = uuid
        if base_dir:
            self.base_dir: str = os.path.join(base_dir, uuid)
        else:
            self.base_dir = ""
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
        self.skip_quality: bool = skip_quality
        self.token_file: str = token_file
        self.token: str = ""
        self.reset: bool = reset
        self.__bam_exists__: bool = False

    def validate(self) -> bool:
        success = True

        # Validate we have all options specified
        check_items = [
            {self.uuid:
                "You must provide a UUID to process (--uuid)"},
            {self.base_dir:
                "You must provide a base directory that will contain results files (--base-dir)"},
            {self.bam_files_dir:
                "You must provide a directory that will contain the bam files (--bam-files-dir)"},
            {self.trimmomatic_jar:
                "You must provide the file location to the Trimmomatic jar file (--trimmomatic-jar)"},
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
                self.create_trimmomatic_jar()
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

        if os.path.exists(self.token_file):
            with open(self.token_file, "r") as fp:
                data = fp.read()
                token = data.strip()
                self.token = token
        return success

    def create_trimmomatic_jar(self, filedir: str) -> str:
        jar_file = os.path.join(filedir, "trimmomatic.jar")
        zip_file = jar_file.replace(".jar", ".zip")

        if not os.path.exists(filedir):
            try:
                os.makedirs(filedir, exist_ok=True)
            except Exception as err:
                log.error(err)
                sys.exit(-1)

        url = "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"
        download_file(url, zip_file, os.path.basename(zip_file), keep=False)
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            file_list = zip_ref.namelist()

            common_prefix = os.path.commonprefix(file_list)

            for file in file_list:
                if file == common_prefix:
                    continue
                extracted_path = os.path.join(
                    filedir, file[len(common_prefix):])

                # Ensure parent directories exist
                os.makedirs(os.path.dirname(
                    extracted_path), exist_ok=True)

                # Extract the file
                try:
                    with zip_ref.open(file) as source, open(extracted_path, 'wb') as dest:
                        dest.write(source.read())
                except IsADirectoryError:
                    os.makedirs(extracted_path, exist_ok=True)
        try:
            os.remove(zip_file)
        except:
            pass
        return jar_file, []

    def create_centrifuge_db(self, filedir):
        gunzip_file = os.path.join(filedir, "genome.tar.gz")
        if not os.path.exists(filedir):
            os.makedirs(filedir, exist_ok=True)
        error_code, output_msg, errors = execute([
            "chmod",
            "-R",
            "775",
            filedir
        ])
        console.print(
            "Please select the Centrifuge database to download and extract:")
        choices = {
            "Refseq: bacteria, archaea, viral, human (compressed)": "https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz",
            "Refseq: bacteria, archaea, viral, human": "https://genome-idx.s3.amazonaws.com/centrifuge/p%2Bh%2Bv.tar.gz"
        }
        choice = select(list(choices.keys()))
        if choice:
            url = choices[choice]
            download_file(url, gunzip_file, os.path.basename(
                gunzip_file), keep=False)

            extr = Spinner(DOTS, "Extracting database...")
            extr.start()
            with tarfile.open(gunzip_file, 'r:gz') as tar:
                tar.extractall(path=os.path.dirname(gunzip_file))
            extr.stop()
            try:
                os.remove(gunzip_file)
            except:
                pass
        return filedir, []

    def create_kraken_db(self, filedir):
        gunzip_file = os.path.join(filedir, "kraken.tar.gz")
        if not os.path.exists(filedir):
            os.makedirs(filedir, exist_ok=True)
        error_code, output_msg, errors = execute([
            "chmod",
            "-R",
            "775",
            filedir
        ])
        console.print(
            "Please select the Kraken database to download and extract:")
        choices = {
            "Standard (55GB)": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz",
            "Standard-8 (8GB)": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz",
            "Standard-16 (16GB)": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240112.tar.gz",
            "nt Database (550GB)": "https://genome-idx.s3.amazonaws.com/kraken/k2_nt_20231129.tar.gz"
        }
        choice = select(list(choices.keys()))
        if choice:
            url = choices[choice]

            download_file(url, gunzip_file, os.path.basename(
                gunzip_file), keep=False)

            extr = Spinner(DOTS, "Extracting database...")
            extr.start()
            with tarfile.open(gunzip_file, 'r:gz') as tar:
                tar.extractall(path=os.path.dirname(gunzip_file))
            extr.stop()
            try:
                os.remove(gunzip_file)
            except:
                pass
        return filedir, []

    def create_hisat_db(self, filedir):
        gunzip_file = os.path.join(filedir, "genome.fa.gz")
        genome_file = os.path.join(filedir, "genome.fa")
        database_file = os.path.join(filedir, "genome")
        if not os.path.exists(filedir):
            os.makedirs(filedir, exist_ok=True)
        error_code, output_msg, errors = execute([
            "chmod",
            "-R",
            "775",
            filedir
        ])
        choices = {
            "Genome sequence (GRCh38.p14)": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.p14.genome.fa.gz"
        }
        console.print(
            "Please select the Hisat database to download and extract:")
        choice = select(list(choices.keys()))
        if choice:
            url = choices[choice]

        download_file(url, gunzip_file, os.path.basename(
            gunzip_file), keep=False)

        extr = Spinner(DOTS, "Extracting database...")
        extr.start()
        with gzip.open(gunzip_file, 'rb') as f_in:
            with open(gunzip_file.replace(".gz", ""), 'wb') as f_out:
                f_out.write(f_in.read())
        try:
            os.remove(gunzip_file)
        except:
            pass
        extr.stop()

        followups = [
            f"hisat2-build -p 4 {genome_file} {database_file}",
            f"rm {genome_file}"
        ]

        return database_file, followups

    def create_splicesites_file(self, filedir):
        gunzip_file = os.path.join(filedir, "gencode_annotations.gtf.gz")
        gtf_file = os.path.join(filedir, "gencode_annotations.gtf")
        final_tsv = os.path.join(filedir, "splicesites.tsv")
        if not os.path.exists(filedir):
            os.makedirs(filedir, exist_ok=True)
        error_code, output_msg, errors = execute([
            "chmod",
            "-R",
            "775",
            filedir
        ])
        console.print(
            "Please select the GRHc GTF files to use to create the splicesites file:")
        choices = {
            "EBI - Comprehensive - CHR": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz",
            "EBI - Comprehensive - ALL": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz",
            "EBI - Comprehensive - PRI": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz",
            "EBI - Basic gene annotation - CHR": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz",
            "EBI - Basic gene annotation - ALL": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.chr_patch_hapl_scaff.basic.annotation.gtf.gz",
            "EBI - Basic gene annotation - PRI": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.basic.annotation.gtf.gz",
            "NCBI Genomic": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
        }
        choice = select(list(choices.keys()))
        if choice:
            url = choices[choice]
            download_file(url, gunzip_file, os.path.basename(
                gunzip_file), keep=False)

            extr = Spinner(DOTS, "Unzipping...")
            extr.start()
            with gzip.open(gunzip_file, 'rb') as f_in:
                with open(gtf_file, 'wb') as f_out:
                    f_out.write(f_in.read())
            extr.stop()
            try:
                os.remove(gunzip_file)
            except:
                pass

            # Extract Splicesites
            splice_sites = []
            extr = Spinner(DOTS, "Extracting sites information...")
            extr.start()
            with open(gtf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        fields = line.strip().split('\t')
                        if fields[2] == 'exon':
                            chrom = fields[0]
                            start = int(fields[3])
                            end = int(fields[4])
                            strand = fields[6]

                            splice_sites.append((chrom, start, end, strand))

            splice_sites.sort(key=lambda x: (x[0], x[3], x[1]))

            with open(final_tsv, 'w') as f_out:
                for i in range(len(splice_sites) - 1):
                    if splice_sites[i][0] == splice_sites[i + 1][0] and splice_sites[i][3] == splice_sites[i + 1][3]:
                        splice_site = splice_sites[i][0], splice_sites[i][2], splice_sites[i +
                                                                                           1][1], splice_sites[i][3]
                        f_out.write('\t'.join(map(str, splice_site)) + '\n')
            extr.stop()
        return final_tsv, []


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
        cmd, list) else cmd if isinstance(cmd, str) else ""
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
        download_file(url, filepath, f"{env.uuid}.bam", headers=headers)
        log.info('Download complete')
        log.info(f'Saved BAM file to {filepath}')
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
        log.info(f'Saved fastq file to {uuid_f1}')
        log.info(f'Saved fastq file to {uuid_f2}')
        if errors:
            log.error(errors)
            sys.exit(-1)

        log.info(f"Conversion complete for {env.uuid}.bam")
    else:
        log.info("Skipping BAM -> FASTA conversion")

    # Check quality
    if not env.start > 3:
        if not env.skip_quality:
            log.info(f"### Checking quality of F*.fq files ###")
            error_code, output, errors = execute([
                "fastqc",
                uuid_f1,
                "--outdir",
                env.fastq_dir
            ])
            if error_code != 0:
                log.error(errors)
                sys.exit(-1)
            error_code, output, errors = execute([
                "fastqc",
                uuid_f2,
                "--outdir",
                env.fastq_dir
            ])
            if error_code != 0:
                log.error(errors)
                sys.exit(-1)
            log.info("Quality check complete")
        else:
            log.info("Skipping quality check due to --skip-quality")
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
            sys.exit(-1)
        log.info("Trimming complete")
        log.info(f'Saved fastq file to {trimmed_f1}')
        log.info(f'Saved fastq file to {trimmed_f1_up}')
        log.info(f'Saved fastq file to {trimmed_f2}')
        log.info(f'Saved fastq file to {trimmed_f2_up}')

        if not env.skip_quality:
            log.info("### Checking quality with fastq ###")
            error_code, output_msg, errors = execute([
                "fastqc",
                trimmed_f1,
                "--outdir",
                f"{env.fastq_dir}"
            ])
            if error_code != 0:
                log.error(errors)
                sys.exit(-1)

            error_code, output_msg, errors = execute([
                "fastqc",
                trimmed_f2,
                "--outdir",
                f"{env.fastq_dir}"
            ])
            if error_code != 0:
                log.error(errors)
                sys.exit(-1)
            log.info("Quality check complete")
        else:
            log.info("Skipping quality check due to --skip-quality")
    else:
        log.info("Skipping trimming with Trimmomatic")

    #  Hisat
    trimmed_sam = f"{env.fastq_dir}/{env.uuid}_trimmed_sam.sam"
    trimmed_sam_compresssed = f"{env.fastq_dir}/{env.uuid}_trimmed_sam.sam.bz"
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
        ], show_output=False)

        if error_code != 0:
            log.error(errors)
            sys.exit(-1)
        log.info(f"Hisat completed")
        log.info(f'Saved file to {trimmed_sam}')

        log.info(f"Compressing sam file")
        error_code, output_msg, errors = execute([
            "samtools",
            "view",
            "-h",
            "-o",
            trimmed_sam_compresssed,
            trimmed_sam
        ])
        if error_code != 0:
            log.error(errors)
            sys.exit(-1)
        log.info(f"Compressing complete")
        log.info(f'Saved file to {trimmed_sam_compresssed}')

        log.info(f"Converting to bam file")
        error_code, output_msg, errors = execute([
            "samtools",
            "view",
            "-bS",
            "-o",
            host_aligned,
            trimmed_sam_compresssed
        ])
        if error_code != 0:
            log.error(errors)
            sys.exit(-1)
        log.info('Samtools complete')
        log.info(f'Saved file to {host_aligned}')

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
            host_aligned,
            "-o",
            unmapped_bam
        ])

        if error_code != 0:
            log.error(errors)
            sys.exit(-1)
        log.info('Extraction complete')
        log.info(f'Saved file to {unmapped_bam}')

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
            sys.exit(-1)
        log.info("Sorting complete")
        log.info(f'Saved file to {unmapped_sorted}')

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
            sys.exit(-1)
        log.info('Conversion complete')
        log.info(f'Saved fastq file to {unmapped_f1}')
        log.info(f'Saved fastq file to {unmapped_f2}')
    else:
        log.info("Skipping extraction")

    if not env.start > 7:
        log.info("### Running assembly of unmapped reads using SPAdes ###")

        error_code, output_msg, errors = execute([
            "rnaspades",
            "-1",
            unmapped_f1,
            "-2",
            unmapped_f2,
            "--only-assembler",
            "-o",
            env.spades_output_dir
        ])

        if error_code != 0:
            log.error(errors)
            sys.exit(-1)
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
            env.kraken_db_dir,
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
            sys.exit(-1)
        log.info('Kraken2 classification complete')
        log.info(
            f'Saved Kraken classifications file to {kraken_classifications}')
        log.info(f'Saved Kraken output file to {kraken_output}')
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
            sys.exit(-1)
        log.info('Pathogen extraction complete complete')
        log.info(f'Saved pathogen nodes file to {kraken_pathogen_nodes}')
    else:
        log.info("Skipping pathogen extraction")

    if not env.start > 10:
        log.info('### Sequencing pathogens ###')
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
            sys.exit(-1)
        log.info('Sequencing complete')
        log.info(
            f'Saved pathogen transcripts file to {kraken_pathogen_transcripts}')
    else:
        log.info("Skipping pathogen sequencing")

    # Blastn DB creation
    blastn_db = f"{env.blastn_output_dir}/{env.uuid}/blastndb"
    kraken_pathogen_gene = f"{env.kraken_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
    if not env.start > 11:
        log.info("### Creating Blastn DB ###")

        os.makedirs(os.path.join(env.blastn_output_dir,
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
            sys.exit(-1)
        log.info('DB build complete')
        log.info(f'Saved Blastn DB to {blastn_db}')
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
            sys.exit(-1)
        log.info('DB query complete')
        log.info(f'Saved kraken pathogen gene file to {kraken_pathogen_gene}')
    else:
        log.info("Skipping Blastn query")

    centrifuge_report = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_report.txt",
    tax_id_list = f"{env.centrifuge_output_dir}/{env.uuid}_tax_id_list.txt"
    centrifuge_output = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_output.txt"
    centrifuge_pathogen_nodes = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_nodes.txt"
    centrifuge_seq_fasta = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_sequences.fasta"
    centrifuge_compressed = f"{env.centrifuge_db_dir}/p_compressed+h+v"
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
            sys.exit(-1)

        log.info('### Extracting pathogens ###')
        log.info(f'Saved Centrifuge report to {centrifuge_report}')

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
            sys.exit(-1)
        log.info('Pathogen extraction complete complete')
        log.info(f'Saved tax ID file to {tax_id_list}')

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
            sys.exit(-1)
        log.info(
            f'Saved Centrifuge pathogen nodes to {centrifuge_pathogen_nodes}')

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
            sys.exit(-1)
        log.info(
            f'Saved Centrifuge sequenced fasta file to {centrifuge_seq_fasta}')

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
            sys.exit(-1)
        log.info('Seqtk and blastn complete')
        log.info(
            f'Saved Centrifuge pathogen gene annotation to {env.centrifuge_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn')
        log.info('Centrifuge classification complete')

    else:
        log.info("Skipping Centrifuge classification")

    if not env.start > 14:
        pass


def pipeline(env: Env):
    pass


def build_environment(env: Env):
    console.clear()
    console.print("### Environment builder ###")

    required_directories = {
        "trimmomatic_jar": {
            "name": "Trimmomatic jar",
            "desc": "This is required for running Trimmomatic",
            "arg": "--trimmomatic-jar",
            "func": env.create_trimmomatic_jar,
            "complete": False,
            "new_path": "",
            "followups": []
        },
        "centrifuge_db_dir": {
            "name": "Centrifuge DB directory",
            "desc": "This is required for Centrifuge classification",
            "arg": "--centrifuge-db-dir",
            "func": env.create_centrifuge_db,
            "complete": False,
            "new_path": "",
            "followups": []
        },
        "kraken_db_dir": {
            "name": "Kraken2 DB folder",
            "desc": "This is required for Kraken classification",
            "arg": "--kraken-db-dir",
            "func": env.create_kraken_db,
            "complete": False,
            "new_path": "",
            "followups": []
        },
        "hisat_db_dir": {
            "name": "Hisat DB directory",
            "desc": "This is required for Hisat processing",
            "arg": "--hisat-db-dir",
            "func": env.create_hisat_db,
            "complete": False,
            "new_path": "",
            "followups": []
        },
        "splicesites_file": {
            "name": "Splicesites file",
            "desc": "This is required for Hisat processing",
            "arg": "--splice-sites-file",
            "func": env.create_splicesites_file,
            "complete": False,
            "new_path": "",
            "followups": []
        }
    }

    for path, details in required_directories.items():
        name = details.get('name')
        desc = details.get('desc')
        arg = details.get('arg')
        func = details.get('func')
        new_path = prompt(f"Please provide the directory for {name}:")
        if new_path:
            required_directories[path]['new_path'] = new_path

    for path, details in required_directories.items():
        name = details.get('name')
        desc = details.get('desc')
        arg = details.get('arg')
        func = details.get('func')
        new_path = details.get('new_path')
        if new_path:
            if func:
                file_path, followups = func(new_path)
                required_directories[path]['new_path'] = file_path
                required_directories[path]['complete'] = True
                required_directories[path]['followups'] = followups

    followups = [v.get('followups')
                 for k, v in required_directories.items() if v.get('followups')]
    console.clear()
    if followups:
        with open('followup.sh', 'w') as pf:
            for k, v in required_directories.items():
                if v.get('followups'):
                    [pf.write(f"{x}\n") for x in v.get('followups')]
        console.print("To complete the build please execute:\n\n./followup.sh")
    console.print(
        "\n### Complete! ###\n\nEnsure you run any follow up commands (if mentioned above) and then the following file / directories can be used:\n")
    for path, details in required_directories.items():
        name = details.get('name')
        desc = details.get('desc')
        arg = details.get('arg')
        func = details.get('func')
        new_path = details.get('new_path')
        if new_path:
            console.print(f"{arg} {new_path}")


def test_environment(env: Env):
    commands = [
        "hisat2-build",
        "bedtools",
        "fastqc",
        "hisat2",
        "samtools",
        "rnaspades",
        "kraken2",
        "seqtk",
        "makeblastdb",
        "blastn",
        "centrifuge"
    ]
    for command in commands:
        error_code, output_msg, errors = execute([
            command
        ],
            show_output=False)
        if error_code == 2:
            log.error(f"{command} - FAILED: The program was not found")
        else:
            log.info(f"{command} - PASSED")


def reset_progress(env: Env):
    if not env.start > 1:
        log.info("Removing BAM files")
        try:
            os.remove(os.path.join(env.bam_files_dir, f"{env.uuid}.bam"))
        except FileNotFoundError:
            pass
        except Exception as err:
            log.error(err)
    else:
        log.info("Skipping BAM files removal")

    if not env.start > 2:
        log.info("Removing BAM->FASTQ files")
        for f in [
            os.path.join(env.fastq_dir, f"{env.uuid}_F1.fq"),
            os.path.join(env.fastq_dir, f"{env.uuid}_F2.fq")
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping BAM->FASTQ converted files")

    if not env.start > 3:
        log.info("Removing FASTQ quality check files")
        for f in [
            os.path.join(env.fastq_dir, f"{env.uuid}_F1_fastqc.html"),
            os.path.join(env.fastq_dir, f"{env.uuid}_F1_fastqc.zip"),
            os.path.join(env.fastq_dir, f"{env.uuid}_F2_fastqc.html"),
            os.path.join(env.fastq_dir, f"{env.uuid}_F2_fastqc.zip"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping FASTQ quality check files")

    if not env.start > 4:
        log.info("Removing trimmed fastq files")
        for f in [
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_F1.fq"),
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_F1_fastqc.html"),
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_F1_fastqc.zip"),
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_F1_UP.fq"),
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_F2.fq"),
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_F2_fastqc.html"),
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_F2_fastqc.zip"),
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_F2_UP.fq"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping trimmed fastq files")

    if not env.start > 5:
        log.info("Removing Hisat alignment files")
        for f in [
            os.path.join(env.fastq_dir, f"{env.uuid}_trimmed_sam"),
            os.path.join(env.fastq_dir, f"{env.uuid}_host_aligned.bam"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping Hisat alignment files")

    if not env.start > 6:
        log.info("Removing extracted / non-human reads files")
        for f in [
            os.path.join(env.fastq_dir, f"{env.uuid}_unmapped.bam"),
            os.path.join(env.fastq_dir, f"{env.uuid}_unmapped_sorted.bam"),
            os.path.join(env.fastq_dir, f"{env.uuid}_unmapped_F1.fq"),
            os.path.join(env.fastq_dir, f"{env.uuid}_unmapped_F2.fq"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping extracted / non-humn reads files")

    if not env.start > 7:
        log.info("Removing assembly of unmapped files")
        execute(
            [
                "rm",
                "-rvf",
                env.spades_output_dir
            ]
        )
    else:
        log.info("Skipping assembly of unmapped files")

    if not env.start > 8:
        log.info("Removing Kraken classification files")
        for f in [
            os.path.join(env.kraken_output_dir,
                         f"{env.uuid}_kraken_report.txt"),
            os.path.join(env.kraken_output_dir,
                         f"{env.uuid}_kraken_classifications.txt"),
            os.path.join(env.kraken_output_dir,
                         f"{env.uuid}_kraken_output.txt"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping Kraken classification files")

    if not env.start > 9:
        log.info("Removing Kraken pathogen nodes files")
        for f in [
            os.path.join(env.kraken_output_dir,
                         f"{env.uuid}_pathogen_nodes.txt"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping Kraken pathogen nodes files")

    if not env.start > 10:
        log.info("Removing Kraken pathogen nodes files")
        for f in [
            os.path.join(env.kraken_output_dir,
                         f"{env.uuid}_pathogen_sequences.fasta"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping Kraken classification files")

    if not env.start > 11:
        log.info("Removing Blastn DB files")
        try:
            execute([
                "rm",
                "-rvf",
                os.path.join(env.blastn_output_dir, f"{env.uuid}/blastndb")
            ])
        except:
            pass
    else:
        log.info("Skipping blastn DB files")

    if not env.start > 12:
        log.info("Removing Blastn query files")
        for f in [
            os.path.join(env.kraken_output_dir,
                         f"{env.uuid}_pathogen_gene_annotation.blastn"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping Blastn query files")

    if not env.start > 13:
        log.info("Removing Centrifuge classification files")
        for f in [
            os.path.join(env.centrifuge_output_dir,
                         f"{env.uuid}_centrifuge_report.txt"),
            os.path.join(env.centrifuge_output_dir,
                         f"{env.uuid}_tax_id_list.txt"),
            os.path.join(env.centrifuge_output_dir,
                         f"{env.uuid}_centrifuge_output.txt"),
            os.path.join(env.centrifuge_output_dir,
                         f"{env.uuid}_pathogen_nodes.txt"),
            os.path.join(env.centrifuge_output_dir,
                         f"{env.uuid}_pathogen_sequences.fasta"),
        ]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
            except Exception as err:
                log.error(err)
    else:
        log.info("Skipping Centrifuge classification files")

    log.info(f"Complete")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reset', action="store_true",
                        help="Resets progress (except dowloading BAM file) for the UUID specified.")
    parser.add_argument('mode', metavar="MODE", choices=[
                        'single', 'paired', 'pipeline', 'build', 'test'], action='store', help="The mode to run, either single, paired, pipeline, build or test. Build will build the environment. Test will test the environment.")
    parser.add_argument('--uuid', action='store',
                        help='The UUID to process.')
    parser.add_argument('--base-dir', action='store',
                        help='The base directory under which all processed files are stored.')
    parser.add_argument('--token', metavar="FILE", action='store',
                        help='The file containing the API Key to download BAM files.')
    parser.add_argument('--start', metavar="LEVEL",
                        choices=["1", "2", "3", "4", "5", "6"], action="store", help="Start at a later stage of processing:\n\n1=Beginning, 2=Skip downloading, 3=Skip conversion to fastq (bedtools), 4=Skip quality checking (fastqc), 5=Skip trimming (Trimmomatic), 6=Skip alignment (Hisat2), 7=Skip extraction (samtools and bedtools), 8=Skip assembly (SPAdes), 9=Skiup Kraken2 classification, 10=Skip pathogen extraction, 11=Skip Blastn db creation, 12=Skip Blastn DB query, 13=Skip Blastn query DB, 14=Skip Centrifuge classification, 15=Skip Blastn classification")
    parser.add_argument('--bam-files-dir', metavar="DIRECTORY",
                        action='store', help="Specify the bam files directory.")
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
    parser.add_argument('--skip-quality', action="store_true",
                        help="Skips all quality checking using fastq")

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
        build = args.build

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
        reset=args.reset if args.reset else False,
        skip_quality=args.skip_quality if args.skip_quality else False
    )

    if env.mode == 'test':
        test_environment(env)
        sys.exit(0)

    if env.mode == 'build':
        build_environment(env)
        sys.exit(0)

    if env.reset:
        reset_progress(env)
        sys.exit(0)

    if not env.validate():
        sys.exit(-1)

    log.info(f"##### Running in {env.mode} mode #####")
    if env.mode == 'single':
        single(env)
    elif env.mode == 'paired':
        paired(env)
    elif env.mode == 'pipeline':
        pipeline(env)


if __name__ in ['__main__', 'builtin', '__builtins__']:
    main()
