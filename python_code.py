import re
import sys
import warnings
import os
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

# Find the Excel file that starts with 'corr_tags'
excel_file = None
for entry in os.listdir(fastq_dir):
    if entry.lower().startswith("corr_tags") and entry.lower().endswith(".xlsx"):
        excel_file = fastq_dir / entry
        break

if excel_file is None:
    print("Excel file not found. Exiting.")
    sys.exit(1)

# Load reference codes and samples from the Excel file
reference_df = pd.read_excel(excel_file)

# Use 'Sample' and 'TAG' as the column names (case-insensitive)
sample_column_name = reference_df.columns[reference_df.columns.str.lower() == 'sample'].tolist()
tag_column_name = reference_df.columns[reference_df.columns.str.lower() == 'tag'].tolist()

if not sample_column_name or not tag_column_name:
    print("Column 'Sample' or 'TAG' not found in the Excel file. Exiting.")
    sys.exit(1)

reference_df.rename(columns={tag_column_name[0]: 'TAG', sample_column_name[0]: 'SAMPLE'}, inplace=True)

samples_to_exclude = ['other', 'OTHER', 'Other']  # Add more if needed

# Read the filenames from the text file
filename = f'/users/llampert/Data_Med/Fieldworks_refs/{directory_name}.txt'  # Update the path accordingly
all_filenames = []

with open(filename, 'r') as file:
    for line in file:
        # Append each line (filename) to the list
        all_filenames.append(line.strip())

# Dictionary to store dataframes for each sample
sample_dataframes = {}

# Inside the loop
for i in tqdm(range(len(all_filenames))):
    forward_file = fastq_dir / all_filenames[i].replace(".gz", "")
    reverse_file = fastq_dir / all_filenames[i].replace("_R1.fastq.gz", "_R2.fastq.gz")

    print("Processing:", forward_file)

    with gzip.open(forward_file, 'rt') as forward_handle, gzip.open(reverse_file, 'rt') as reverse_handle:
        for forward_record, reverse_record in zip(SeqIO.parse(forward_handle, format='fastq'),
                                                  SeqIO.parse(reverse_handle, format='fastq')):
            forward_id = forward_record.id
            forward_seq = forward_record.seq.lower()

            reverse_id = reverse_record.id
            reverse_seq = reverse_record.seq.lower()

            # Extract tag from the ID
            tag_match = re.match(r'^(\S+)', forward_id)
            tag = tag_match.group(1) if tag_match else None

            # Check if the tag is in the Excel file
            if tag is not None and tag in reference_df['TAG'].tolist():
                # Get the corresponding sample name
                sample_name = reference_df.loc[reference_df['TAG'] == tag, 'SAMPLE'].values[0]

                # Check if the sample name should be excluded
                if any(sample_name.lower().startswith(prefix) for prefix in samples_to_exclude):
                    print(f"Skipping {sample_name} as it should be excluded.")
                    continue

                # Check if the sample name exists in the dictionary, if not create a new dataframe
                if sample_name not in sample_dataframes:
                    sample_dataframes[sample_name] = {'ids': [], 'seqs_forward': [], 'seqs_reverse': []}

                sample_dataframes[sample_name]['ids'].append(forward_id)
                sample_dataframes[sample_name]['seqs_forward'].append(str(forward_seq))
                sample_dataframes[sample_name]['seqs_reverse'].append(str(reverse_seq))

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
