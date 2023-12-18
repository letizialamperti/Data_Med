import re
import warnings
import pandas as pd
import gzip
from Bio import SeqIO
from pathlib import Path
from tqdm.notebook import tqdm

fastq_dir = Path('/store/sdsc/sd29/med_data_wp3/AERMC_Banyuls')

# Filenames list
all_filenames = [
    "AIMI%2FAIMI-20211012a%2Fdata%2F211015_SN1126_A_L001_AIMI-402_R1.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211019_SN234_A_L001_AIMI-407_R1.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211015_SN1126_A_L001_AIMI-402_R2.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211019_SN234_A_L001_AIMI-407_R2.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211015_SN6662_A_L001_AIMI-404_R1.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211026_SN7280_A_L001_AIMI-403_R1.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211015_SN6662_A_L001_AIMI-404_R2.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211026_SN7280_A_L001_AIMI-403_R2.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211015_SN7280_A_L001_AIMI-406_R1.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211027_SN234_A_L001_AIMI-405_R1.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211015_SN7280_A_L001_AIMI-406_R2.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211027_SN234_A_L001_AIMI-405_R2.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211019_SN1126_A_L001_AIMI-408_R1.fastq.gz",
    "AIMI%2FAIMI-20211012a%2Fdata%2F211019_SN1126_A_L001_AIMI-408_R2.fastq.gz"
    # Add other filenames here...
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
    
    store_dir = Path('/scratch/snx3000/llampert/MED_CSV/AERMC_Banyuls')
    store_dir.mkdir(parents=True, exist_ok=True)
    save_file = store_dir / f'{sample_name}.csv'
    df = df.dropna()
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
