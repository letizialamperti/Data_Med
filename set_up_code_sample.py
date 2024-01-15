import re
import sys
import warnings
import pandas as pd
import gzip
from Bio import SeqIO
from pathlib import Path
from tqdm.notebook import tqdm
import os

# Check if the directory name argument is provided
if len(sys.argv) < 2:
    print("Please provide the directory name as an argument.")
    sys.exit(1)

# Get the directory name from the command-line argument
directory_name = sys.argv[1]
fastq_dir = Path(f'/store/sdsc/sd29/med_data_wp3/{directory_name}')

# Find the Excel file case-insensitively
excel_file = None
for entry in os.listdir(fastq_dir):
    if entry.lower().startswith("corr_tags") and entry.lower().endswith(f".xlsx"):
        excel_file = fastq_dir / entry
        break

if excel_file is None or not excel_file.exists():
    print(f"Excel file not found for {directory_name}.")
    sys.exit(1)

# Load reference codes and samples from the Excel file
reference_df = pd.read_excel(excel_file)
reference_codes = reference_df['ReferenceCode'].tolist()
samples_to_exclude = ['other', 'OTHER']  # Add more if needed

# Read the filenames from the text file
filename = f'/users/llampert/Data_Med/Fieldworks_refs/{directory_name}.txt'  # Update the path accordingly
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

# Dictionary to store dataframes for each sample
sample_dataframes = {}

# Inside the loop
for i in tqdm(range(num_files)):
    file = forward_reads_files[i]
    print("Processing:", file)
    # Extract sample name from the filename directly
    sample_name = file.stem.split('_R')[0]

    # Check if the sample name should be excluded
    if any(sample_name.lower().startswith(prefix) for prefix in samples_to_exclude):
        print(f"Skipping {sample_name} as it should be excluded.")
        continue

    # Check if the sample name exists in the dictionary, if not create a new dataframe
    if sample_name not in sample_dataframes:
        sample_dataframes[sample_name] = {'ids': [], 'seqs_forward': [], 'seqs_reverse': []}

    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, format='fastq'):
            id = record.id
            seq = record.seq.lower()

            # Check if the ID contains any of the reference codes
            if any(ref_code in id for ref_code in reference_codes):
                sample_dataframes[sample_name]['ids'].append(id)
                sample_dataframes[sample_name]['seqs_forward'].append(str(seq))

    file = reverse_reads_files[i]
    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, format='fastq'):
            id = record.id
            seq = record.seq.lower()

            # Check if the ID contains any of the reference codes
            if id in sample_dataframes[sample_name]['ids']:
                sample_dataframes[sample_name]['seqs_reverse'].append(str(seq))
            else:
                warnings.warn("ID of reverse read could not be matched to any forward read.")

# Save dataframes to CSV files
for sample_name, data in sample_dataframes.items():
    df = pd.DataFrame(data={'Forward': data['seqs_forward'], 'Reverse': data['seqs_reverse']}, index=data['ids'])

    store_dir = Path(f'/scratch/snx3000/llampert/MED_SAMPLES_CSV/{directory_name}')
    store_dir.mkdir(parents=True, exist_ok=True)

    # Save the CSV file with the sample name
    save_file = store_dir / f'{sample_name}.csv'
    df = df.dropna()
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
