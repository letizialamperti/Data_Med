import re
import sys  
import warnings
import pandas as pd
import gzip
from Bio import SeqIO
from pathlib import Path
from tqdm.notebook import tqdm

# Check if the directory name argument is provided
if len(sys.argv) < 2:
    print("Please provide the directory name as an argument.")
    sys.exit(1)

# Get the directory name from the command-line argument
directory_name = sys.argv[1]
fastq_dir = Path(f'/store/sdsc/sd29/med_data_wp3/{directory_name}')

# Read the filenames from the text file
filename = f'/users/llampert/Data_Med/{directory_name}.txt'  # Update the path accordingly
all_filenames = []

with open(filename, 'r') as file:
    for line in file:
        # Append each line (filename) to the list
        all_filenames.append(line.strip())

# Print the list of filenames
print(all_filenames)

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

    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, format='fastq'):
            id = record.id
            seq = record.seq.lower()

            if id in df.index:
                df.at[id, 'Reverse'] = str(seq)
            else:
                warnings.warn("ID of reverse read could not be matched to any forward read.")
    

    store_dir = Path(f'/scratch/snx3000/llampert/MED_CSV/{directory_name}')
    store_dir.mkdir(parents=True, exist_ok=True)
    save_file = store_dir / f'{sample_name}.csv'
    df = df.dropna()
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
