import re
import warnings
import pandas as pd
import gzip
from Bio import SeqIO
from pathlib import Path
from tqdm.notebook import tqdm

fastq_dir = Path('/store/sdsc/sd29/med_data_wp3/AERMC')

# Modify the pattern to match gzipped files with ".fastq.gz" extension
forward_reads_files = list(fastq_dir.glob("*R1.fastq.gz"))
reverse_reads_files = list(fastq_dir.glob("*R2.fastq.gz"))

forward_reads_files.sort()
reverse_reads_files.sort()

assert len(forward_reads_files) == len(reverse_reads_files)

def get_sample_name(filename: str) -> str:
    regex = r'\d+_NB\d+_A_L1-4_AIMI-\d+_(R[12])\.fastq\.gz'
    matches = re.findall(pattern=regex, string=filename)
    if len(matches) != 1:
        print(f"Filename '{filename}' does not match the pattern.")
        return ''  # Or any other handling you need

    sample_name = matches[0].split("_R")[0]
    return sample_name

# Your filenames, corrected
all_filenames = [
    "210920_SN1126_A_L001_AIMI-393_R1.fastq.gz",    "220201_NB501473_A_L1-4_AIMI-469_R1.fastq.gz",
    "210920_SN1126_A_L001_AIMI-393_R2.fastq.gz",    "220201_NB501473_A_L1-4_AIMI-469_R2.fastq.gz",
    "210927_SN234_A_L001_AIMI-393_R1.fastq.gz",     "220201_NB501473_A_L1-4_AIMI-470_R1.fastq.gz",
    "210927_SN234_A_L001_AIMI-393_R2.fastq.gz",     "220201_NB501473_A_L1-4_AIMI-470_R2.fastq.gz",
    "210930_SN6662_A_L001_AIMI-391_R1.fastq.gz",    "220201_NB501473_A_L1-4_AIMI-471_R1.fastq.gz",
    "210930_SN6662_A_L001_AIMI-391_R2.fastq.gz",    "220201_NB501473_A_L1-4_AIMI-471_R2.fastq.gz",
    "211001_SN6662_A_L001_AIMI-392_R1.fastq.gz",    "220201_NB501473_A_L1-4_AIMI-472_R1.fastq.gz",
    "211001_SN6662_A_L001_AIMI-392_R2.fastq.gz",    "220201_NB501473_A_L1-4_AIMI-472_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-449_R1.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-473_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-449_R2.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-473_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-452_R1.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-474_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-452_R2.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-474_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-453_R1.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-475_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-453_R2.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-475_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-455_R1.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-476_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-455_R2.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-476_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-458_R1.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-477_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-458_R2.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-477_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-468_R1.fastq.gz",  "220201_NB501473_A_L1-4_AIMI-468_R2
