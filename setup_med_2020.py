import re
import warnings
import pandas as pd
import gzip
from Bio import SeqIO
from pathlib import Path
from tqdm.notebook import tqdm

fastq_dir = Path('/store/sdsc/sd29/med_data_wp3/med_2020')

# Modify the pattern to match gzipped files with ".fastq.gz" extension
forward_reads_files = list(fastq_dir.glob("*R1.fastq.gz"))
reverse_reads_files = list(fastq_dir.glob("*R2.fastq.gz"))

forward_reads_files.sort()
reverse_reads_files.sort()

assert len(forward_reads_files) == len(reverse_reads_files)

def get_sample_name(filename: str) -> str:
    regex = fr'\d+_SN\d+_A_L001_AIMI-\d+_R\d\.fastq\.gz'
    matches = re.findall(pattern=regex, string=filename)
    assert len(matches) == 1

    sample_name = matches[0].split("_R")[0]
    return sample_name

num_files = len(forward_reads_files)

print("Forward files:", forward_reads_files)
print("Reverse files:", reverse_reads_files)

# Inside the loop
for i in tqdm(range(num_files)):
    file = forward_reads_files[i]
    print("Processing:", file)
    
for i in tqdm(range(num_files)):
    file = forward_reads_files[i]
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
    file = reverse_reads_files[i]
    assert sample_name == get_sample_name(file.name)

    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, format='fastq'):
            id = record.id
            seq = record.seq.lower()

            if id in df.index:
                df.at[id, 'Reverse'] = str(seq)
            else:
                warnings.warn("ID of reverse read could not be matched to any forward read.")
    
    store_dir = Path('/scratch/snx3000/llampert/MED_CSV')
    store_dir.mkdir(parents=True, exist_ok=True)
    save_file = store_dir / f'{sample_name}.csv'
    df = df.dropna()
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
