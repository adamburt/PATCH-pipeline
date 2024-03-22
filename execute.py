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
from io import TextIOWrapper
import multiprocessing


class Env():
    def __init__(self, mode: str, uuid: str, base_dir: str, bam_files_dir: str, trimmomatic_jar: str, hisat_db_dir: str, splice_sites_file: str, kraken_db_dir: str, human_genome_lib: str, pathogen: str, centrifuge_db_dir: str, token_file: str, start: int = 1, reset: bool = False, skip_quality: bool = False):
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
        self.human_genome_lib: str = human_genome_lib
        self.pathogen: str = pathogen
        self.skip_quality: bool = skip_quality
        self.token_file: str = token_file
        self.token: str = ""
        self.reset: bool = reset
        self.__bam_exists__: bool = False
        self.cpu_count: int = 1

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
            {self.kraken_db_dir:
                "You must provide a directory that contains the kraken2 database files (--kraken-db-dir)"},
            {self.token_file:
                "You must provide the file that contains the token for use with downloading BAM files (--token)"}
        ]
        for item in check_items:
            for c, v in item.items():
                if not c:
                    success = False
                    log.error(v)
        if not self.hisat_db_dir and not self.human_genome_lib:
            success = False
            log.error("You must provide at least one of a hisat DB directory (--hisat-db-dir), or a human genome library (--human-genome). Human genome wins if both are provided.")
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

        if self.splice_sites_file and not os.path.exists(self.splice_sites_file):
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
        print(
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
        print(
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
        print(
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
        print(
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
                for line in tqdm(f):
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

    def create_human_genome(self, filedir):
        gunzip_file = os.path.join(filedir, "library.fna.gz")
        library_file = os.path.join(filedir, "library.fna")
        choices = {
            "NCBI - GRCh38 latest": "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz",
        }
        if not os.path.exists(filedir):
            os.makedirs(filedir, exist_ok=True)
        error_code, output_msg, errors = execute([
            "chmod",
            "-R",
            "775",
            filedir
        ])

        print(
            "Please select the Human genome file to use:")
        choice = select(list(choices.keys()))
        if choice:
            url = choices[choice]

            download_file(url, gunzip_file, os.path.basename(
                gunzip_file), keep=False)

            extr = Spinner(DOTS, "Extracting library...")
            extr.start()
            with gzip.open(gunzip_file, 'rb') as f_in:
                with open(gunzip_file.replace(".gz", ""), 'wb') as f_out:
                    f_out.write(f_in.read())
            try:
                os.remove(gunzip_file)
            except:
                pass
            extr.stop()
        return library_file, [f"bwa index {library_file}"]

        pass


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
    log.info(f"Executing: {' '.join(cmd)}")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    if show_output:
        while proc.poll() == None:
            for stdout_line in iter(proc.stdout.readline, ""):
                print(stdout_line)
    proc.stdout.close()
    return_code = proc.wait()
    output_msg, error_msg = proc.communicate()
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


def add_to_output(fp: TextIOWrapper, command: str, header: str = "", echo: str = "", complete_msg: list = [], supress_outpout: bool = False):
    if header:
        fp.write(f"\n### {header} ###\n")
    if echo:
        s = ""
        for x in range(0, len(echo) + 8):
            s = f"{s}#"
        fp.write('echo ""\n')
        fp.write(f'echo "{s}"\n')
        fp.write(f'echo "### {echo} ###"\n')
        fp.write(f'echo "{s}"\n')
    fp.write(command)
    if supress_outpout:
        fp.write("  > /dev/null 2>&1")
    fp.write("\n")
    if complete_msg:
        for msg in complete_msg:
            fp.write(f'echo "{msg}"\n')


def download_bam_files(env: Env):
    log.info('### Downloading bam files ###')
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


def single(env: Env):

    # Create output file
    with open("execute.sh", "w") as fp:

        # Download BAM files
        if env.__bam_exists__:
            log.warning(
                f"{env.uuid}.bam already exists, overwriting. To avoid this in the future use --start 2")
        url = f"https://api.gdc.cancer.gov/data/{env.uuid}"
        add_to_output(
            fp,
            f'curl -J -H "X-Auth-Token: {env.token}" -o "{env.bam_files_dir}/{env.uuid}.bam" {url}',
            header="Download BAM files",
            echo="Downloading BAM files",
            complete_msg=[
                f"BAM file saved to {env.bam_files_dir}/{env.uuid}.bam"
            ]
        )

        # Convert BAM files
        uuid_f1 = f"{env.fastq_dir}/{env.uuid}.fq"
        log.info(f'### Converting {env.uuid}.bam file to FASTA format ###')
        bam_source = os.path.join(env.bam_files_dir, f"{env.uuid}.bam")

        if os.path.exists(os.path.join(env.fastq_dir, f"{env.uuid}_F1.fq")) and os.path.exists(os.path.join(env.fastq_dir, f"{env.uuid}_F2.fq")):
            log.warning(
                f"{env.fastq_dir}/{env.uuid}_F1.fq already exists, overwriting. To avoid this in the future use --start 3")

        add_to_output(
            fp,
            f"bedtools bamtofastq -i {bam_source} -fq {uuid_f1}",
            header="Convert BAM to FASTQ",
            echo=f"Converting {env.uuid}.bam to FASTQ",
            complete_msg=[
                f'echo "Saved fastq file to {uuid_f1}"'
            ]
        )

        # Check quality
        if not env.skip_quality:
            add_to_output(
                fp,
                f"fastqc {uuid_f1} --outdir {env.fastq_dir}",
                header="Quality check with FASTQ",
                echo=f"Checking quality of {os.path.basename(uuid_f1)} with FASTQ"
            )
        else:
            log.info("Skipping quality check due to --skip-quality")

        # Trimming
        trimmed_f1 = f"{env.fastq_dir}/{env.uuid}_trimmed.fq"
        add_to_output(
            fp,
            f"java -jar {env.trimmomatic_jar} SE {uuid_f1} {trimmed_f1} -threads {env.cpu_count} LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28",
            header="Trimming files with Trimmomatic",
            echo="Trimming files with Trimmomatic",
            complete_msg=[
                f"Saved fastq file to {trimmed_f1}"
            ]
        )

        if not env.skip_quality:
            add_to_output(
                fp,
                f"fastqc {trimmed_f1} --outdir {env.fastq_dir}",
                header=f"Quality check {trimmed_f1} with FASTQ",
                echo=f"Checking quality of {os.path.basename(trimmed_f1)} with FASTQ"
            )
        else:
            log.info("Skipping quality check due to --skip-quality")

        # Hisat or bwa
        trimmed_sam = f"{env.fastq_dir}/{env.uuid}_trimmed_sam.sam"
        host_aligned = f"{env.fastq_dir}/{env.uuid}_host_aligned.bam"
        use_bwa = False

        if env.human_genome_lib and env.hisat_db_dir:
            print("You have specified both hisat directories and bwa directories. Which tool would you like to use for alignment (default: bwa)?")
            selection = select(["bwa", "hisat"])
            if selection:
                if selection == "bwa":
                    use_bwa = True
                else:
                    use_bwa = False
            else:
                use_bwa = True
        else:
            if env.human_genome_lib:
                use_bwa = True
            else:
                use_bwa = False

        # BWA
        if use_bwa:
            command = f"bwa mem -t {env.cpu_count} {env.human_genome_lib} {trimmed_f1} -o {trimmed_sam}"

            add_to_output(
                fp,
                command,
                header="Run BWA and alignment",
                echo="Running BWA and alignment",
                complete_msg=[
                    f"Saved file to {trimmed_sam}"
                ]
            )

            add_to_output(
                fp,
                f"samtools sort -o {host_aligned} {trimmed_sam}",
                header="Converting to bam file",
                echo=f"Converting {os.path.basename(trimmed_sam)} to bam file",
                complete_msg=[
                    f"Saved file to {host_aligned}"
                ]
            )

        # Hisat
        else:
            command = f"hisat2 -x {env.hisat_db_dir}/genome"
            if env.splice_sites_file:
                command = f"{command} --known-splicesite-infile {env.splice_sites_file}"
            command = f"{command} -U {trimmed_f1} -S {trimmed_sam}"

            add_to_output(
                fp,
                command,
                header="Running Hisat and alignment",
                echo="Running Hisat and alignment",
                complete_msg=[
                    f"Saved file to {trimmed_sam}"
                ]
            )

            add_to_output(
                fp,
                f"samtools view -bS -o {host_aligned} {trimmed_sam}",
                header="Converting to bam file",
                echo=f"Converting {os.path.basename(trimmed_sam)} to bam file",
                complete_msg=[
                    f"Saved file to {host_aligned}"
                ]
            )

        # Extract unmapped reads
        unmapped_bam = f"{env.fastq_dir}/{env.uuid}_unmapped.bam"
        unmapped_sorted = f"{env.fastq_dir}/{env.uuid}_unmapped_sorted.bam"
        unmapped_sorted_fastq = f"{env.fastq_dir}/{env.uuid}_unmapped_sorted.fq"

        add_to_output(
            fp,
            f"samtools view -F 4 {host_aligned} -o {unmapped_bam}",
            header="Extracting unmapped / non-human reads",
            echo="Extracting unmapped / non-human reads",
            complete_msg=[
                f"Saved file to {unmapped_bam}"
            ]
        )

        add_to_output(
            fp,
            f"samtools sort {unmapped_bam} -o {unmapped_sorted}",
            header="Sorting results",
            echo="Sorting results",
            complete_msg=[
                f"Saved file to {unmapped_sorted}"
            ]
        )

        add_to_output(
            fp,
            f"bedtools bamtofastq -i {unmapped_sorted} -fq {unmapped_sorted_fastq}",
            header="Converting results to fastq",
            echo="Converting results to fastq",
            complete_msg=[
                f"Saved file to {unmapped_sorted_fastq}"
            ]
        )

        # SPAdes assembly
        add_to_output(
            fp,
            f"spades -s {unmapped_sorted_fastq} --only-assembler -o {env.spades_output_dir}",
            header="SPAdes assembly",
            echo="Running assembly of unmapped reads using SPAdes",
            complete_msg=[
                "Assembly complete"
            ]
        )

        #  Kraken classification
        kraken_report = f"{env.kraken_output_dir}/{env.uuid}_kraken_report.txt"
        kraken_classifications = f"{env.kraken_output_dir}/{env.uuid}_kraken_classifications.txt"
        kraken_output = f"{env.kraken_output_dir}/{env.uuid}_output_kraken.txt"
        spades_transcripts = f"{env.spades_output_dir}/{env.uuid}/transcripts.fasta"

        add_to_output(
            fp,
            f"kraken2 --use-names --db {env.kraken_db_dir} --report {kraken_report} --classified-out {kraken_classifications} --output {kraken_output} {spades_transcripts}",
            header="Kraken classification",
            echo="Running Kraken classification",
            complete_msg=[
                f"Saved Kraken classifications file to {kraken_classifications}"
            ]
        )

        # Pathogen extraction
        kraken_pathogen_nodes = f"{env.kraken_output_dir}/{env.uuid}_pathogen_nodes.txt"

        if env.pathogen:
            command = f"grep -e {env.pathogen} {kraken_output} | awk '%REPLACE%' > {kraken_pathogen_nodes}".replace(
                "%REPLACE%", "{print $2}")
        else:
            command = f"cat {kraken_output} | awk '%REPLACE%' > {kraken_pathogen_nodes}".replace(
                "%REPLACE%", "{print $2}")
        log.info("### Starting pathogen extraction ###")

        add_to_output(
            fp,
            command,
            header="Pathogen extraction",
            echo="Starting pathogen extraction",
            complete_msg=[
                "Pathogen extraction complete complete",
                f"Saved pathogen nodes file to {kraken_pathogen_nodes}"
            ]
        )

        # Sequence pathogens
        kraken_pathogen_transcripts = f"{env.kraken_output_dir}/{env.uuid}_pathogen_sequences.fasta"
        add_to_output(
            fp,
            f"seqtk subseq {spades_transcripts} {kraken_pathogen_nodes} > {kraken_pathogen_transcripts}",
            header="Sequencing pathogens",
            echo="Sequencing pathogens",
            complete_msg=[
                "Sequencing complete",
                f"Saved pathogen transcripts file to {kraken_pathogen_transcripts}"
            ]
        )

        # Blastn DB creation
        blastn_db = f"{env.blastn_output_dir}/{env.uuid}/blastndb"
        os.makedirs(
            os.path.join(
                env.blastn_output_dir,
                env.uuid
            ),
            exist_ok=True
        )
        blastn_db = f"{env.blastn_output_dir}/{env.uuid}/blastndb"
        kraken_pathogen_gene = f"{env.kraken_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
        add_to_output(
            fp,
            f"makeblastdb -in {kraken_pathogen_transcripts} -dbtype nucl -out {blastn_db}",
            header="Create Blastn DB",
            echo="Creating Blastn DB",
            complete_msg=[
                "Blastn DB build complete",
                f"Saved Blastn DB to {blastn_db}"
            ]
        )

        # Blastn DB query
        kraken_pathogen_gene = f"{env.kraken_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
        add_to_output(
            fp,
            f"blastn -query {kraken_pathogen_transcripts} -db {blastn_db} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out {kraken_pathogen_gene}",
            header="Blastn query",
            echo="Starting Blastn query",
            complete_msg=[
                "Blastm DB query complete",
                f"Saved kraken pathogen gene file to {kraken_pathogen_gene}"
            ]
        )

        # Centrifuge classification
        centrifuge_report = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_report.txt",
        tax_id_list = f"{env.centrifuge_output_dir}/{env.uuid}_tax_id_list.txt"
        centrifuge_output = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_output.txt"
        centrifuge_pathogen_nodes = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_nodes.txt"
        centrifuge_seq_fasta = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_sequences.fasta"
        centrifuge_compressed = f"{env.centrifuge_db_dir}/p_compressed+h+v"

        add_to_output(
            fp,
            f"centrifuge -x {centrifuge_compressed} -f {spades_transcripts} --report-file {centrifuge_report} -S {centrifuge_output}",
            header="Centrifuge classification",
            echo="Running Centrifuge classification",
            complete_msg=[
                "Classifications complete",
                f"Saved Centrifuge report to {centrifuge_report}"
            ]
        )

        if env.pathogen:
            command = f"grep -e '{env.pathogen}' {centrifuge_report} | awk '%REPLACE%' | sort | uniq > {tax_id_list}".replace(
                "%REPLACE%", "{print $3}")
        else:
            command = f"cat {centrifuge_report} | awk '%REPLACE%' |  sort | uniq > {tax_id_list}".replace(
                "%REPLACE%", "{print $3}")
        add_to_output(
            fp,
            command,
            header="Extract pathogens",
            echo="Extracting pathogens",
            complete_msg=[
                f"Saved tax ID file to {tax_id_list}"
            ]
        )

        add_to_output(
            fp,
            f"awk -F' ' 'NR==FNR%REPLACE1%' {tax_id_list} {centrifuge_output} | awk '%REPLACE2%' | sort -u | uniq > {centrifuge_pathogen_nodes}".replace(
                "%REPLACE1%", "{c[$1]++;next};c[$3]").replace("%REPLACE2%", "{print $1}"),
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge pathogen nodes to {centrifuge_pathogen_nodes}"
            ]
        )

        add_to_output(
            fp,
            f"seqtk subseq {spades_transcripts} {centrifuge_pathogen_nodes} > {centrifuge_seq_fasta}",
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge sequenced fasta file to {centrifuge_seq_fasta}"
            ]
        )

        add_to_output(
            fp,
            f"blastn -query {centrifuge_seq_fasta} -db {env.blastn_output_dir}/{env.uuid}_blastndb -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out {env.centrifuge_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn",
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge pathogen gene annotation to {env.centrifuge_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
            ]
        )

    log.info("Saved output to execute.sh")


def paired(env: Env):

    # Create output file
    with open("execute.sh", "w") as fp:

        # Download BAM files
        if env.__bam_exists__:
            log.warning(
                f"{env.uuid}.bam already exists, overwriting. To avoid this in the future use --start 2")
        url = f"https://api.gdc.cancer.gov/data/{env.uuid}"
        add_to_output(
            fp,
            f'curl -J -H "X-Auth-Token: {env.token}" -o "{env.bam_files_dir}/{env.uuid}.bam" {url}',
            header="Download BAM files",
            echo="Downloading BAM files",
            complete_msg=[
                f"BAM file saved to {env.bam_files_dir}/{env.uuid}.bam"
            ]
        )

        # Convert BAM files
        uuid_f1 = f"{env.fastq_dir}/{env.uuid}_F1.fq"
        uuid_f2 = f"{env.fastq_dir}/{env.uuid}_F2.fq"
        bam_source = os.path.join(env.bam_files_dir, f"{env.uuid}.bam")
        bam_aligned_file = os.path.join(env.bam_aligned_dir, f"{env.uuid}.bam")

        if os.path.exists(os.path.join(env.fastq_dir, f"{env.uuid}_F1.fq")) and os.path.exists(os.path.join(env.fastq_dir, f"{env.uuid}_F2.fq")):
            log.warning(
                f"{env.fastq_dir}/{env.uuid}_F1.fq already exists, overwriting. To avoid this in the future use --start 3")
            log.warning(
                f"{env.fastq_dir}/{env.uuid}_F2.fq already exists, overwriting. To avoid this in the future use --start 3")
        add_to_output(
            fp,
            f"samtools sort -n {bam_source} -o {bam_aligned_file}",
            header="Sort BAM file",
            echo=f"Sorting {bam_source} file",
            complete_msg=[
                f'echo "Saved sorted BAM file to {bam_aligned_file}"'
            ]
        )
        add_to_output(
            fp,
            f"bedtools bamtofastq -i {bam_aligned_file} -fq {uuid_f1} -fq2 {uuid_f2}",
            header="Convert BAM to FASTQ",
            echo=f"Converting {env.uuid}.bam to FASTQ",
            complete_msg=[
                f'echo "Saved fastq file to {uuid_f1}"',
                f'echo "Saved fastq file to {uuid_f2}"'
            ]
        )

        # Check quality
        if not env.skip_quality:
            add_to_output(
                fp,
                f"fastqc {uuid_f1} --outdir {env.fastq_dir}",
                header="Quality check with FASTQ",
                echo=f"Checking quality of {os.path.basename(uuid_f1)} with FASTQ"
            )
            add_to_output(
                fp,
                f"fastqc {uuid_f2} --outdir {env.fastq_dir}",
                header="Quality check with FASTQ",
                echo=f"Checking quality of {os.path.basename(uuid_f2)} with FASTQ"
            )
        else:
            log.info("Skipping quality check due to --skip-quality")

        # Trimming
        trimmed_f1 = f"{env.fastq_dir}/{env.uuid}_trimmed_F1.fq"
        trimmed_f1_up = f"{env.fastq_dir}/{env.uuid}_trimmed_F1_UP.fq"
        trimmed_f2 = f"{env.fastq_dir}/{env.uuid}_trimmed_F2.fq"
        trimmed_f2_up = f"{env.fastq_dir}/{env.uuid}_trimmed_F2_UP.fq"

        add_to_output(
            fp,
            f"java -jar {env.trimmomatic_jar} PE {uuid_f1} {uuid_f2} {trimmed_f1} {trimmed_f1_up} {trimmed_f2} {trimmed_f2_up} -threads {env.cpu_count} -threads {env.cpu_count} LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28",
            header="Trimming files with Trimmomatic",
            echo="Trimming files with Trimmomatic",
            complete_msg=[
                f"Saved fastq file to {trimmed_f1}",
                f"Saved fastq file to {trimmed_f1_up}",
                f"Saved fastq file to {trimmed_f2}",
                f"Saved fastq file to {trimmed_f2_up}",
            ]
        )

        # Quality check
        if not env.skip_quality:
            add_to_output(
                fp,
                f"fastqc {trimmed_f1} --outdir {env.fastq_dir}",
                header=f"Quality check {trimmed_f1} with FASTQ",
                echo=f"Checking quality of {os.path.basename(trimmed_f1)} with FASTQ"
            )
            add_to_output(
                fp,
                f"fastqc {trimmed_f2} --outdir {env.fastq_dir}",
                header=f"Quality check {trimmed_f2} with FASTQ",
                echo=f"Checking quality of {os.path.basename(trimmed_f2)} with FASTQ"
            )
        else:
            log.info("Skipping quality check due to --skip-quality")

        # Hisat or bwa
        trimmed_sam = f"{env.fastq_dir}/{env.uuid}_trimmed_sam.sam"
        host_aligned = f"{env.fastq_dir}/{env.uuid}_host_aligned.bam"

        use_bwa = False

        if env.human_genome_lib and env.hisat_db_dir:
            print("You have specified both hisat directories and bwa directories. Which tool would you like to use for alignment (default: bwa)?")
            selection = select(["bwa", "hisat"])
            if selection:
                if selection == "bwa":
                    use_bwa = True
                else:
                    use_bwa = False
            else:
                use_bwa = True
        else:
            if env.human_genome_lib:
                use_bwa = True
            else:
                use_bwa = False

        # BWA
        if use_bwa:
            trimmed_sam = f"{env.fastq_dir}/{env.uuid}_trimmed_sam.sam"
            host_aligned = f"{env.fastq_dir}/{env.uuid}_host_aligned.bam"

            command = f"bwa mem -t {env.cpu_count} {env.human_genome_lib} {trimmed_f1} {trimmed_f2} -o {trimmed_sam}"

            add_to_output(
                fp,
                command,
                header="Run BWA and alignment",
                echo="Running BWA and alignment",
                complete_msg=[
                    f"Saved file to {trimmed_sam}"
                ]
            )

            add_to_output(
                fp,
                f"samtools sort -o {host_aligned} {trimmed_sam}",
                header="Converting to bam file",
                echo=f"Converting {os.path.basename(trimmed_sam)} to bam file",
                complete_msg=[
                    f"Saved file to {host_aligned}"
                ]
            )

        # Hisat
        else:
            command = f"hisat2 -x {env.hisat_db_dir}/genome"
            if env.splice_sites_file:
                command = f"{command} --known-splicesite-infile {env.splice_sites_file}"
            command = f"{command} -1 {trimmed_f1} -2 {trimmed_f2} -S {trimmed_sam}"
            add_to_output(
                fp,
                command,
                header="Running Hisat and alignment",
                echo="Running Hisat and alignment",
                complete_msg=[
                    f"Saved file to {trimmed_sam}"
                ]
            )
            add_to_output(
                fp,
                f"samtools view -bS -o {host_aligned} {trimmed_sam}",
                header="Converting to bam file",
                echo=f"Converting {os.path.basename(trimmed_sam)} to bam file",
                complete_msg=[
                    f"Saved file to {host_aligned}"
                ]
            )

        # Extract unmapped reads
        unmapped_bam = f"{env.fastq_dir}/{env.uuid}_unmapped.bam"
        unmapped_sorted = f"{env.fastq_dir}/{env.uuid}_unmapped_sorted.bam"
        unmapped_f1 = f"{env.fastq_dir}/{env.uuid}_unmapped_F1.fq"
        unmapped_f2 = f"{env.fastq_dir}/{env.uuid}_unmapped_F2.fq"

        add_to_output(
            fp,
            f"samtools view -F 4 {host_aligned} -o {unmapped_bam}",
            header="Extracting unmapped / non-human reads",
            echo="Extracting unmapped / non-human reads",
            complete_msg=[
                f"Saved file to {unmapped_bam}"
            ]
        )

        add_to_output(
            fp,
            f"samtools sort {unmapped_bam} -o {unmapped_sorted}",
            header="Sorting results",
            echo="Sorting results",
            complete_msg=[
                f"Saved file to {unmapped_sorted}"
            ]
        )

        add_to_output(
            fp,
            f"bedtools bamtofastq -i {unmapped_sorted} -fq {unmapped_f1} -fq2 {unmapped_f2}",
            header="Converting results to fastq",
            echo="Converting results to fastq",
            complete_msg=[
                f"Saved file to {unmapped_f1}",
                f"Saved file to {unmapped_f2}"
            ]
        )

        # SPAdes assembly
        add_to_output(
            fp,
            f"rnaspades -1 {unmapped_f1} -2 {unmapped_f2} --only-assembler -o {env.spades_output_dir}",
            header="SPAdes assembly",
            echo="Running assembly of unmapped reads using SPAdes",
            complete_msg=[
                "Assembly complete"
            ]
        )

        #  Kraken classification
        kraken_report = f"{env.kraken_output_dir}/{env.uuid}_kraken_report.txt"
        kraken_classifications = f"{env.kraken_output_dir}/{env.uuid}_kraken_classifications.txt"
        kraken_output = f"{env.kraken_output_dir}/{env.uuid}_output_kraken.txt"
        spades_transcripts = f"{env.spades_output_dir}/{env.uuid}/transcripts.fasta"
        kraken_pathogen_nodes = f"{env.kraken_output_dir}/{env.uuid}_pathogen_nodes.txt"
        kraken_pathogen_transcripts = f"{env.kraken_output_dir}/{env.uuid}_pathogen_sequences.fasta"

        add_to_output(
            fp,
            f"kraken2 --use-names --db {env.kraken_db_dir} --report {kraken_report} --classified-out {kraken_classifications} --output {kraken_output} {spades_transcripts}",
            header="Kraken classification",
            echo="Running Kraken classification",
            complete_msg=[
                f"Saved Kraken classifications file to {kraken_classifications}"
            ]
        )

        # Pathogen extraction
        if env.pathogen:
            command = f"grep -e {env.pathogen} {kraken_output} | awk '%REPLACE%' > {kraken_pathogen_nodes}".replace(
                "%REPLACE%", "{print $2}")
        else:
            command = f"cat {kraken_output} | awk '%REPLACE%' > {kraken_pathogen_nodes}".replace(
                "%REPLACE%", "{print $2}")
        add_to_output(
            fp,
            command,
            header="Pathogen extraction",
            echo="Starting pathogen extraction",
            complete_msg=[
                "Pathogen extraction complete complete",
                f"Saved pathogen nodes file to {kraken_pathogen_nodes}"
            ]
        )

        # Sequence pathogens
        add_to_output(
            fp,
            f"seqtk subseq {spades_transcripts} {kraken_pathogen_nodes} > {kraken_pathogen_transcripts}",
            header="Sequencing pathogens",
            echo="Sequencing pathogens",
            complete_msg=[
                "Sequencing complete",
                f"Saved pathogen transcripts file to {kraken_pathogen_transcripts}"
            ]
        )

        # Blastn DB creation
        os.makedirs(
            os.path.join(
                env.blastn_output_dir,
                env.uuid
            ),
            exist_ok=True
        )
        blastn_db = f"{env.blastn_output_dir}/{env.uuid}/blastndb"
        kraken_pathogen_gene = f"{env.kraken_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
        add_to_output(
            fp,
            f"makeblastdb -in {kraken_pathogen_transcripts} -dbtype nucl -out {blastn_db}",
            header="Create Blastn DB",
            echo="Creating Blastn DB",
            complete_msg=[
                "Blastn DB build complete",
                f"Saved Blastn DB to {blastn_db}"
            ]
        )

        # Blastn DB query

        add_to_output(
            fp,
            f"blastn -query {kraken_pathogen_transcripts} -db {blastn_db} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out {kraken_pathogen_gene}",
            header="Blastn query",
            echo="Starting Blastn query",
            complete_msg=[
                "Blastm DB query complete",
                f"Saved kraken pathogen gene file to {kraken_pathogen_gene}"
            ]
        )

        # Centrifuge classification
        centrifuge_report = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_report.txt"
        tax_id_list = f"{env.centrifuge_output_dir}/{env.uuid}_tax_id_list.txt"
        centrifuge_output = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_output.txt"
        centrifuge_pathogen_nodes = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_nodes.txt"
        centrifuge_seq_fasta = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_sequences.fasta"
        centrifuge_compressed = f"{env.centrifuge_db_dir}/p_compressed+h+v"

        add_to_output(
            fp,
            f"centrifuge -x {centrifuge_compressed} -f {spades_transcripts} --report-file {centrifuge_report} -S {centrifuge_output}",
            header="Centrifuge classification",
            echo="Running Centrifuge classification",
            complete_msg=[
                "Classifications complete",
                f"Saved Centrifuge report to {centrifuge_report}"
            ]
        )

        if env.pathogen:
            command = f"grep -e '{env.pathogen}' {centrifuge_report} | awk '%REPLACE%' | sort | uniq > {tax_id_list}".replace(
                "%REPLACE%", "{print $3}")
        else:
            command = f"cat {centrifuge_report} | awk '%REPLACE%' |  sort | uniq > {tax_id_list}".replace(
                "%REPLACE%", "{print $3}")
        add_to_output(
            fp,
            command,
            header="Extract pathogens",
            echo="Extracting pathogens",
            complete_msg=[
                f"Saved tax ID file to {tax_id_list}"
            ]
        )

        add_to_output(
            fp,
            f"awk -F' ' 'NR==FNR%REPLACE1%' {tax_id_list} {centrifuge_output} | awk '%REPLACE2%' | sort -u | uniq > {centrifuge_pathogen_nodes}".replace(
                "%REPLACE1%", "{c[$1]++;next};c[$3]").replace("%REPLACE2%", "{print $1}"),
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge pathogen nodes to {centrifuge_pathogen_nodes}"
            ]
        )

        add_to_output(
            fp,
            f"seqtk subseq {spades_transcripts} {centrifuge_pathogen_nodes} > {centrifuge_seq_fasta}",
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge sequenced fasta file to {centrifuge_seq_fasta}"
            ]
        )

        add_to_output(
            fp,
            f"blastn -query {centrifuge_seq_fasta} -db {env.blastn_output_dir}/{env.uuid}_blastndb -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out {env.centrifuge_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn",
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge pathogen gene annotation to {env.centrifuge_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
            ]
        )

    log.info("Saved output to execute.sh")


def pipeline(env: Env):

    return

    # Create output file
    with open("execute.sh", "w") as fp:

        # Download BAM files
        if env.__bam_exists__:
            log.warning(
                f"{env.uuid}.bam already exists, overwriting. To avoid this in the future use --start 2")
        url = f"https://api.gdc.cancer.gov/data/{env.uuid}"
        add_to_output(
            fp,
            f'curl -J -H "X-Auth-Token: {env.token}" -o "{env.bam_files_dir}/{env.uuid}.bam" {url}',
            header="Download BAM files",
            echo="Downloading BAM files",
            complete_msg=[
                f"BAM file saved to {env.bam_files_dir}/{env.uuid}.bam"
            ]
        )

        # Convert BAM files
        uuid_f1 = f"{env.fastq_dir}/{env.uuid}_F1.fq.gz"
        uuid_f2 = f"{env.fastq_dir}/{env.uuid}_F2.fq.gz"
        uuid_s = f"{env.fastq_dir}/{env.uuid}_s.fq.gz"
        uuid_0 = f"{env.fastq_dir}/{env.uuid}_0.fq.gz"
        uuid_02 = f"{env.fastq_dir}/{env.uuid}_02.fq.gz"
        bam_source = os.path.join(env.bam_files_dir, f"{env.uuid}.bam")
        bam_aligned_file = os.path.join(env.bam_aligned_dir, f"{env.uuid}.bam")

        if os.path.exists(os.path.join(env.fastq_dir, f"{env.uuid}_F1.fq")) and os.path.exists(os.path.join(env.fastq_dir, f"{env.uuid}_F2.fq")):
            log.warning(
                f"{env.fastq_dir}/{env.uuid}_F1.fq already exists, overwriting. To avoid this in the future use --start 3")
            log.warning(
                f"{env.fastq_dir}/{env.uuid}_F2.fq already exists, overwriting. To avoid this in the future use --start 3")
        add_to_output(
            fp,
            f"samtools sort -n {bam_source} -o {bam_aligned_file}",
            header="Sort BAM file",
            echo=f"Sorting {bam_source} file",
            complete_msg=[
                f'echo "Saved sorted BAM file to {bam_aligned_file}"'
            ]
        )
        add_to_output(
            fp,
            f"bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename={bam_aligned_file} F={uuid_f1} F2={uuid_f2} S={uuid_s} inputformat=bam 0={uuid_0} 02={uuid_02} tryoq=1 gz=1 level=5",
            header="Convert BAM to FASTQ",
            echo=f"Converting {env.uuid}.bam to FASTQ",
            complete_msg=[
                f'echo "Saved fastq file to {uuid_f1}"',
                f'echo "Saved fastq file to {uuid_f2}"'
            ]
        )

        # Check quality
        if not env.skip_quality:
            add_to_output(
                fp,
                f"fastqc {uuid_f1} --outdir {env.fastq_dir}",
                header="Quality check with FASTQ",
                echo=f"Checking quality of {os.path.basename(uuid_f1)} with FASTQ"
            )
            add_to_output(
                fp,
                f"fastqc {uuid_f2} --outdir {env.fastq_dir}",
                header="Quality check with FASTQ",
                echo=f"Checking quality of {os.path.basename(uuid_f2)} with FASTQ"
            )
        else:
            log.info("Skipping quality check due to --skip-quality")

        # Trimming
        trimmed_f1 = f"{env.fastq_dir}/{env.uuid}_trimmed_F1.fq"
        trimmed_f1_up = f"{env.fastq_dir}/{env.uuid}_trimmed_F1_UP.fq"
        trimmed_f2 = f"{env.fastq_dir}/{env.uuid}_trimmed_F2.fq"
        trimmed_f2_up = f"{env.fastq_dir}/{env.uuid}_trimmed_F2_UP.fq"

        add_to_output(
            fp,
            f"java -jar {env.trimmomatic_jar} PE {uuid_f1} {uuid_f2} {trimmed_f1} {trimmed_f1_up} {trimmed_f2} {trimmed_f2_up} -threads {env.cpu_count} LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:28",
            header="Trimming files with Trimmomatic",
            echo="Trimming files with Trimmomatic",
            complete_msg=[
                f"Saved fastq file to {trimmed_f1}",
                f"Saved fastq file to {trimmed_f1_up}",
                f"Saved fastq file to {trimmed_f2}",
                f"Saved fastq file to {trimmed_f2_up}",
            ]
        )

        # Quality check
        if not env.skip_quality:
            add_to_output(
                fp,
                f"fastqc {trimmed_f1} --outdir {env.fastq_dir}",
                header=f"Quality check {trimmed_f1} with FASTQ",
                echo=f"Checking quality of {os.path.basename(trimmed_f1)} with FASTQ"
            )
            add_to_output(
                fp,
                f"fastqc {trimmed_f2} --outdir {env.fastq_dir}",
                header=f"Quality check {trimmed_f2} with FASTQ",
                echo=f"Checking quality of {os.path.basename(trimmed_f2)} with FASTQ"
            )
        else:
            log.info("Skipping quality check due to --skip-quality")

        # Align to combined reference genome
        trimmed_sam = f"{env.fastq_dir}/{env.uuid}_trimmed_sam.sam"
        host_aligned = f"{env.fastq_dir}/{env.uuid}_host_aligned.bam"
        command = f"bwa mem {env.kraken_db_dir} -x {env.hisat_db_dir}/genome"
        add_to_output(
            fp,
            command,
            header="Running Hisat and alignment",
            echo="Running Hisat and alignment",
            complete_msg=[
                f"Saved file to {trimmed_sam}"
            ]
        )
        add_to_output(
            fp,
            f"samtools view -bS -o {host_aligned} {trimmed_sam}",
            header="Converting to bam file",
            echo=f"Converting {os.path.basename(trimmed_sam)} to bam file",
            complete_msg=[
                f"Saved file to {host_aligned}"
            ]
        )

        # Extract unmapped reads
        unmapped_bam = f"{env.fastq_dir}/{env.uuid}_unmapped.bam"
        unmapped_sorted = f"{env.fastq_dir}/{env.uuid}_unmapped_sorted.bam"
        unmapped_f1 = f"{env.fastq_dir}/{env.uuid}_unmapped_F1.fq"
        unmapped_f2 = f"{env.fastq_dir}/{env.uuid}_unmapped_F2.fq"

        add_to_output(
            fp,
            f"samtools view -F 4 {host_aligned} -o {unmapped_bam}",
            header="Extracting unmapped / non-human reads",
            echo="Extracting unmapped / non-human reads",
            complete_msg=[
                f"Saved file to {unmapped_bam}"
            ]
        )

        add_to_output(
            fp,
            f"samtools sort {unmapped_bam} -o {unmapped_sorted}",
            header="Sorting results",
            echo="Sorting results",
            complete_msg=[
                f"Saved file to {unmapped_sorted}"
            ]
        )

        add_to_output(
            fp,
            f"bedtools bamtofastq -i {unmapped_sorted} -fq {unmapped_f1} -fq2 {unmapped_f2}",
            header="Converting results to fastq",
            echo="Converting results to fastq",
            complete_msg=[
                f"Saved file to {unmapped_f1}",
                f"Saved file to {unmapped_f2}"
            ]
        )

        # SPAdes assembly
        add_to_output(
            fp,
            f"rnaspades -1 {unmapped_f1} -2 {unmapped_f2} --only-assembler -o {env.spades_output_dir}",
            header="SPAdes assembly",
            echo="Running assembly of unmapped reads using SPAdes",
            complete_msg=[
                "Assembly complete"
            ]
        )

        #  Kraken classification
        kraken_report = f"{env.kraken_output_dir}/{env.uuid}_kraken_report.txt"
        kraken_classifications = f"{env.kraken_output_dir}/{env.uuid}_kraken_classifications.txt"
        kraken_output = f"{env.kraken_output_dir}/{env.uuid}_output_kraken.txt"
        spades_transcripts = f"{env.spades_output_dir}/{env.uuid}/transcripts.fasta"
        kraken_pathogen_nodes = f"{env.kraken_output_dir}/{env.uuid}_pathogen_nodes.txt"
        kraken_pathogen_transcripts = f"{env.kraken_output_dir}/{env.uuid}_pathogen_sequences.fasta"

        add_to_output(
            fp,
            f"kraken2 --use-names --db {env.kraken_db_dir} --report {kraken_report} --classified-out {kraken_classifications} --output {kraken_output} {spades_transcripts}",
            header="Kraken classification",
            echo="Running Kraken classification",
            complete_msg=[
                f"Saved Kraken classifications file to {kraken_classifications}"
            ]
        )

        # Pathogen extraction
        if env.pathogen:
            command = f"grep -e {env.pathogen} {kraken_output} | awk '%REPLACE%' > {kraken_pathogen_nodes}".replace(
                "%REPLACE%", "{print $2}")
        else:
            command = f"cat {kraken_output} | awk '%REPLACE%' > {kraken_pathogen_nodes}".replace(
                "%REPLACE%", "{print $2}")
        add_to_output(
            fp,
            command,
            header="Pathogen extraction",
            echo="Starting pathogen extraction",
            complete_msg=[
                "Pathogen extraction complete complete",
                f"Saved pathogen nodes file to {kraken_pathogen_nodes}"
            ]
        )

        # Sequence pathogens
        add_to_output(
            fp,
            f"seqtk subseq {spades_transcripts} {kraken_pathogen_nodes} > {kraken_pathogen_transcripts}",
            header="Sequencing pathogens",
            echo="Sequencing pathogens",
            complete_msg=[
                "Sequencing complete",
                f"Saved pathogen transcripts file to {kraken_pathogen_transcripts}"
            ]
        )

        # Blastn DB creation
        os.makedirs(
            os.path.join(
                env.blastn_output_dir,
                env.uuid
            ),
            exist_ok=True
        )
        blastn_db = f"{env.blastn_output_dir}/{env.uuid}/blastndb"
        kraken_pathogen_gene = f"{env.kraken_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
        add_to_output(
            fp,
            f"makeblastdb -in {kraken_pathogen_transcripts} -dbtype nucl -out {blastn_db}",
            header="Create Blastn DB",
            echo="Creating Blastn DB",
            complete_msg=[
                "Blastn DB build complete",
                f"Saved Blastn DB to {blastn_db}"
            ]
        )

        # Blastn DB query

        add_to_output(
            fp,
            f"blastn -query {kraken_pathogen_transcripts} -db {blastn_db} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out {kraken_pathogen_gene}",
            header="Blastn query",
            echo="Starting Blastn query",
            complete_msg=[
                "Blastm DB query complete",
                f"Saved kraken pathogen gene file to {kraken_pathogen_gene}"
            ]
        )

        # Centrifuge classification
        centrifuge_report = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_report.txt"
        tax_id_list = f"{env.centrifuge_output_dir}/{env.uuid}_tax_id_list.txt"
        centrifuge_output = f"{env.centrifuge_output_dir}/{env.uuid}_centrifuge_output.txt"
        centrifuge_pathogen_nodes = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_nodes.txt"
        centrifuge_seq_fasta = f"{env.centrifuge_output_dir}/{env.uuid}_pathogen_sequences.fasta"
        centrifuge_compressed = f"{env.centrifuge_db_dir}/p_compressed+h+v"

        add_to_output(
            fp,
            f"centrifuge -x {centrifuge_compressed} -f {spades_transcripts} --report-file {centrifuge_report} -S {centrifuge_output}",
            header="Centrifuge classification",
            echo="Running Centrifuge classification",
            complete_msg=[
                "Classifications complete",
                f"Saved Centrifuge report to {centrifuge_report}"
            ]
        )

        if env.pathogen:
            command = f"grep -e '{env.pathogen}' {centrifuge_report} | awk '%REPLACE%' | sort | uniq > {tax_id_list}".replace(
                "%REPLACE%", "{print $3}")
        else:
            command = f"cat {centrifuge_report} | awk '%REPLACE%' |  sort | uniq > {tax_id_list}".replace(
                "%REPLACE%", "{print $3}")
        add_to_output(
            fp,
            command,
            header="Extract pathogens",
            echo="Extracting pathogens",
            complete_msg=[
                f"Saved tax ID file to {tax_id_list}"
            ]
        )

        add_to_output(
            fp,
            f"awk -F' ' 'NR==FNR%REPLACE1%' {tax_id_list} {centrifuge_output} | awk '%REPLACE2%' | sort -u | uniq > {centrifuge_pathogen_nodes}".replace(
                "%REPLACE1%", "{c[$1]++;next};c[$3]").replace("%REPLACE2%", "{print $1}"),
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge pathogen nodes to {centrifuge_pathogen_nodes}"
            ]
        )

        add_to_output(
            fp,
            f"seqtk subseq {spades_transcripts} {centrifuge_pathogen_nodes} > {centrifuge_seq_fasta}",
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge sequenced fasta file to {centrifuge_seq_fasta}"
            ]
        )

        add_to_output(
            fp,
            f"blastn -query {centrifuge_seq_fasta} -db {env.blastn_output_dir}/{env.uuid}_blastndb -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 1 -max_hsps 1 -out {env.centrifuge_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn",
            header="Run seqtk and blastn",
            echo="Running seqtk and blastn",
            complete_msg=[
                f"Saved Centrifuge pathogen gene annotation to {env.centrifuge_output_dir}/{env.uuid}_pathogen_gene_annotation.blastn"
            ]
        )

    log.info("Saved output to execute.sh")


def build_environment(env: Env):
    console.clear()
    print("### Environment builder ###")

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
        "human_genome": {
            "name": "Human genome reference library",
            "desc": "This is required for aligning to a reference genome",
            "arg": "--human-genome",
            "func": env.create_human_genome,
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
        print("To complete the build please execute:\n\n./followup.sh")
    print(
        "\n### Complete! ###\n\nEnsure you run any follow up commands (if mentioned above) and then the following file / directories can be used:\n")
    for path, details in required_directories.items():
        name = details.get('name')
        desc = details.get('desc')
        arg = details.get('arg')
        func = details.get('func')
        new_path = details.get('new_path')
        if new_path:
            print(f"{arg} {new_path}")


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
    parser.add_argument('--human-genome', metavar="FILE",
                        action='store', help="Specify the human genome library fna file.")
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
        args.human_genome,
        args.pathogen,
        args.centrifuge_db_dir,
        args.token,
        start=args.start if args.start else 1,
        reset=args.reset if args.reset else False,
        skip_quality=args.skip_quality if args.skip_quality else False
    )
    cpu_count = multiprocessing.cpu_count()
    env.cpu_count = cpu_count

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
