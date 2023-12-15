import re
import warnings
import pandas as pd
import gzip
from Bio import SeqIO
from pathlib import Path
from tqdm.notebook import tqdm

fastq_dir = Path('/store/sdsc/sd29/med_data_wp3/AERMC_Banyuls')

# Filenames list
filenames = [
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

# Modify the pattern to match filenames
def get_sample_name(filename: str) -> str:
    regex = r'\d+_SN\d+_A_L001_AIMI-\d+_R\d'
    matches = re.findall(pattern=regex, string=filename)
    assert len(matches) == 1

    sample_name = matches[0]
    return sample_name

num_files = len(filenames)

print("Filenames:", filenames)

# Inside the loop
for i in tqdm(range(num_files)):
    file = Path(filenames[i])
    print("Processing:", file)
    
for i in tqdm(range(num_files)):
    file = Path(filenames[i])
    sample_name = get_sample_name(file.name)
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
    reverse_file = file.parent / (sample_name + "_R2.fastq.gz")
    assert sample_name == get_sample_name(reverse_file.name)

    with gzip.open(reverse_file, 'rt') as handle:
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

