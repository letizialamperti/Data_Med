import re
import warnings
import pandas as pd
import gzip
from Bio import SeqIO
from pathlib import Path
from tqdm.notebook import tqdm

fastq_dir = Path('/store/sdsc/sd29/med_data_wp3/AERMC')

# Your filenames, each on a new line
all_filenames = [
    "210920_SN1126_A_L001_AIMI-393_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-469_R1.fastq.gz",
    "210920_SN1126_A_L001_AIMI-393_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-469_R2.fastq.gz",
    "210927_SN234_A_L001_AIMI-393_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-470_R1.fastq.gz",
    "210927_SN234_A_L001_AIMI-393_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-470_R2.fastq.gz",
    "210930_SN6662_A_L001_AIMI-391_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-471_R1.fastq.gz",
    "210930_SN6662_A_L001_AIMI-391_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-471_R2.fastq.gz",
    "211001_SN6662_A_L001_AIMI-392_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-472_R1.fastq.gz",
    "211001_SN6662_A_L001_AIMI-392_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-472_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-449_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-473_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-449_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-473_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-452_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-474_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-452_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-474_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-453_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-475_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-453_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-475_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-455_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-476_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-455_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-476_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-458_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-477_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-458_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-477_R2.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-468_R1.fastq.gz",
    "220201_NB501473_A_L1-4_AIMI-468_R2.fastq.gz"
]

forward_reads_files = [fastq_dir / filename for filename in all_filenames if "_R1.fastq.gz" in filename]
reverse_reads_files = [fastq_dir / filename for filename in all_filenames if "_R2.fastq.gz" in filename]

forward_reads_files.sort()
reverse_reads_files.sort()

assert len(forward_reads_files) == len(reverse_reads_files)

num_files = len(forward_reads_files)

print("Forward files:", forward_reads_files)
print("Reverse files:", reverse_reads_files)

# Inside the loop
for i in tqdm(range(num_files)):
    file = forward_reads_files[i]
    print("Processing:", file)
    # Extract sample name from the filename directly
    sample_name = file.stem.split('_R')[0]
    ids = []
    seqs = []

    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, format='fastq'):
            id = record.id
            seq = record.seq.lower()
            
            ids.append(id)
            seqs.append(str(seq))

    df = pd.DataFrame(data=seqs, index=ids, columns=['Forward'])

    df['Reverse'] = ''
    file = reverse_reads_files[i]
    # Extract sample name from the filename directly
    assert sample_name == file.stem.split('_R')[0]

    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, format='fastq'):
            id = record.id
            seq = record.seq.lower()

            if id in df.index:
                df.at[id, 'Reverse'] = str(seq)
            else:
                warnings.warn("ID of reverse read could not be matched to any forward read.")
    
    store_dir = Path('/scratch/snx3000/llampert/MED_CSV/AERMC')
    store_dir.mkdir(parents=True, exist_ok=True)
    save_file = store_dir / f'{sample_name}.csv'
    df = df.dropna()
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
