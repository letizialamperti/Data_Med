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

# Define forward and reverse reads files
forward_reads_files = [fastq_dir / filename for filename in all_filenames if "_R1.fastq.gz" in filename]
reverse_reads_files = [fastq_dir / filename for filename in all_filenames if "_R2.fastq.gz" in filename]

forward_reads_files.sort()
reverse_reads_files.sort()

num_files = len(forward_reads_files)  # Define num_files

# Inside the loop
for i in tqdm(range(num_files)):
    file = forward_reads_files[i]
    print("Processing:", file)

    # Extract RUN name from the filename
    run_name_match = re.search(r'RUN_(.*?)_R1', file.stem)
    if not run_name_match:
        print(f"Skipping {file}: RUN name not found in filename.")
        continue

    run_name = run_name_match.group(1)

    # Subselect metadata from Excel for the same RUN
    run_metadata = reference_df[reference_df['RUN'].str.contains(run_name, case=False, na=False)]

    # Check if any metadata is found for the current RUN
    if not run_metadata.empty:
        # Extract unique tags for the current RUN
        unique_tags = run_metadata['TAG'].unique()

        # Extract sample name from the filename directly
        sample_name = run_metadata['SAMPLE'].iloc[0]

        # Extract unique sample name (excluding repeats before underscore)
        unique_sample_name = re.sub(r'_\d+', '', sample_name)

        # Check if the sample name should be excluded
        if any(unique_sample_name.lower().startswith(prefix) for prefix in samples_to_exclude):
            print(f"Skipping {sample_name} as it should be excluded.")
            continue

        # Check if the unique sample name exists in the dictionary, if not create a new dataframe
        if unique_sample_name not in unique_sample_dataframes:
            unique_sample_dataframes[unique_sample_name] = {'ids': [], 'seqs_forward': [], 'seqs_reverse': []}

        with gzip.open(file, 'rt') as handle:
            for record in SeqIO.parse(handle, format='fastq'):
                id = record.id
                seq = record.seq.lower()

                # Check if the ID contains any of the unique tags for the current RUN
                if any(tag in id for tag in unique_tags):
                    unique_sample_dataframes[unique_sample_name]['ids'].append(id)
                    unique_sample_dataframes[unique_sample_name]['seqs_forward'].append(str(seq))

        file = reverse_reads_files[i]
        with gzip.open(file, 'rt') as handle:
            for record in SeqIO.parse(handle, format='fastq'):
                id = record.id
                seq = record.seq.lower()

                # Check if the ID contains any of the unique tags for the current RUN
                if id in unique_sample_dataframes[unique_sample_name]['ids']:
                    unique_sample_dataframes[unique_sample_name]['seqs_reverse'].append(str(seq))
                else:
                    warnings.warn("ID of reverse read could not be matched to any forward read.")

# Save dataframes to CSV files
for unique_sample_name, data in tqdm(unique_sample_dataframes.items(), desc="Saving CSVs"):
    df = pd.DataFrame(data={'Forward': data['seqs_forward'], 'Reverse': data['seqs_reverse']}, index=data['ids'])

    store_dir = Path(f'/scratch/snx3000/llampert/MED_SAMPLES_CSV/{directory_name}')
    store_dir.mkdir(parents=True, exist_ok=True)

    # Save the CSV file with the unique sample name
    save_file = store_dir / f'{unique_sample_name}.csv'
    df = df.dropna()
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
