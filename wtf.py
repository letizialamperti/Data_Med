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

# Initialize sample_dataframes dictionary
sample_dataframes = {}

# Use 'Sample' as the column name (case-insensitive)
sample_column_name = reference_df.columns[reference_df.columns.str.lower() == 'sample'].tolist()

if not sample_column_name:
    print("Column 'Sample' not found in the Excel file. Exiting.")
    sys.exit(1)

reference_codes = reference_df[sample_column_name[0]].tolist()
samples_to_exclude = ['other', 'OTHER', 'Other']  # Add more if needed

# Grouping files and corresponding Excel rows by RUN
file_groups = {}
for file in os.listdir(fastq_dir):
    if file.endswith(("_R1.fastq.gz", "_R2.fastq.gz")):
        run_match = re.search(r'RUN_(.*?)_', file)
        if run_match:
            run_name = run_match.group(1)
            if run_name not in file_groups:
                file_groups[run_name] = {'files': [], 'excel_rows': []}
            file_groups[run_name]['files'].append(file)
            file_groups[run_name]['excel_rows'].append(reference_df[reference_df['RUN'] == run_name])

# Loop over file groups
for run_name, file_group in tqdm(file_groups.items(), desc="Processing RUNs"):
    print(f"Processing RUN: {run_name}")
    
    # Combine Excel rows for the current RUN
    run_excel_df = pd.concat(file_group['excel_rows'], ignore_index=True)

    # Loop over files in the current RUN
    for file in tqdm(file_group['files'], desc="Processing files", leave=False):
        # Extract sample name from the filename
        sample_name = file.stem.split('_R')[0]

        # Check if the sample name should be excluded
        if any(sample_name.lower().startswith(prefix) for prefix in samples_to_exclude):
            print(f"Skipping {sample_name} as it should be excluded.")
            continue

        # Check if the sample name exists in the dictionary, if not create a new dataframe
        if sample_name not in sample_dataframes:
            sample_dataframes[sample_name] = {'ids': [], 'seqs_forward': [], 'seqs_reverse': []}

        # Process the corresponding Excel rows for the current sample and RUN
        for _, excel_row in run_excel_df.iterrows():
            reference_code = excel_row['TAG']
            # Process forward reads
            forward_file = fastq_dir / f"{run_name}_{reference_code}_R1.fastq.gz"
            with gzip.open(forward_file, 'rt') as handle:
                for record in SeqIO.parse(handle, format='fastq'):
                    id = record.id
                    seq = record.seq.lower()
                    # Check if the ID contains any of the reference codes
                    if any(ref_code in id for ref_code in reference_codes):
                        sample_dataframes[sample_name]['ids'].append(id)
                        sample_dataframes[sample_name]['seqs_forward'].append(str(seq))

            # Process reverse reads
            reverse_file = fastq_dir / f"{run_name}_{reference_code}_R2.fastq.gz"
            with gzip.open(reverse_file, 'rt') as handle:
                for record in SeqIO.parse(handle, format='fastq'):
                    id = record.id
                    seq = record.seq.lower()
                    # Check if the ID contains any of the reference codes
                    if id in sample_dataframes[sample_name]['ids']:
                        sample_dataframes[sample_name]['seqs_reverse'].append(str(seq))
                    else:
                        warnings.warn("ID of reverse read could not be matched to any forward read.")

# Save dataframes to CSV files
for sample_name, data in tqdm(sample_dataframes.items(), desc="Saving CSVs"):
    df = pd.DataFrame(data={'Forward': data['seqs_forward'], 'Reverse': data['seqs_reverse']}, index=data['ids'])
    
    store_dir = Path(f'/scratch/snx3000/llampert/MED_SAMPLES_CSV/{directory_name}')
    store_dir.mkdir(parents=True, exist_ok=True)

    # Save the CSV file with the sample name
    save_file = store_dir / f'{sample_name}.csv'
    df = df.dropna()
    df = df.sample(frac=1)  # randomize
    df.to_csv(save_file, index=False)
