import re
import sys
import logging
import warnings
import os
import pandas as pd
import gzip
from Bio import SeqIO
from pathlib import Path
from tqdm.notebook import tqdm

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def exit_with_error(message):
    logging.error(message)
    sys.exit(1)

def locate_excel_file(directory):
    for entry in os.listdir(directory):
        if entry.lower().startswith("corr_tags") and entry.lower().endswith(".xlsx"):
            return Path(directory) / entry
    return None

def load_metadata(excel_file):
    try:
        reference_df = pd.read_excel(excel_file)
        sample_column_name = reference_df.columns[reference_df.columns.str.lower() == 'sample'].tolist()
        tag_column_name = reference_df.columns[reference_df.columns.str.lower() == 'tag'].tolist()
        if not sample_column_name or not tag_column_name:
            exit_with_error("Column 'Sample' or 'TAG' not found in the Excel file. Exiting.")
        reference_df.rename(columns={tag_column_name[0]: 'TAG', sample_column_name[0]: 'SAMPLE'}, inplace=True)
        return reference_df
    except Exception as e:
        exit_with_error(f"Error loading metadata: {str(e)}")

def read_filename_list(filename):
    try:
        with open(filename, 'r') as file:
            return [line.strip() for line in file]
    except Exception as e:
        exit_with_error(f"Error reading filename list: {str(e)}")

def process_file(directory, filename, reference_df, samples_to_exclude, unique_sample_dataframes):
    logging.info("Processing: %s", filename)
    run_name_match = re.search(r'_(.{8})_[R12]', filename)
    if not run_name_match:
        logging.warning(f"Skipping {filename}: RUN name not found in filename.")
        return

    run_name = run_name_match.group(1)
    logging.info("Extracted RUN name: %s", run_name)

    run_metadata = reference_df[reference_df['RUN'].str.contains(run_name, case=False, na=False)]
    logging.info("Run metadata: %s", run_metadata)

    if run_metadata.empty:
        logging.warning("No metadata found for the current RUN.")
        return

    unique_tags = run_metadata['TAG'].unique()
    sample_name = run_metadata['SAMPLE'].iloc[0]
    unique_sample_name = re.sub(r'_\d+', '', sample_name)

    if any(unique_sample_name.lower().startswith(prefix) for prefix in samples_to_exclude):
        return

    if unique_sample_name not in unique_sample_dataframes:
        unique_sample_dataframes[unique_sample_name] = {'ids': [], 'seqs_forward': [], 'seqs_reverse': []}

    forward_file = directory / f"{unique_sample_name}_R1.fastq.gz"
    reverse_file = directory / f"{unique_sample_name}_R2.fastq.gz"

    try:
        with gzip.open(forward_file, 'rt') as handle:
            for record in SeqIO.parse(handle, format='fastq'):
                id = record.id
                seq = record.seq.lower()
                if any(tag in id for tag in unique_tags):
                    unique_sample_dataframes[unique_sample_name]['ids'].append(id)
                    unique_sample_dataframes[unique_sample_name]['seqs_forward'].append(str(seq))

        with gzip.open(reverse_file, 'rt') as handle:
            for record in SeqIO.parse(handle, format='fastq'):
                id = record.id
                seq = record.seq.lower()
                if id in unique_sample_dataframes[unique_sample_name]['ids']:
                    unique_sample_dataframes[unique_sample_name]['seqs_reverse'].append(str(seq))
                else:
                    warnings.warn("ID of reverse read could not be matched to any forward read.")
    except Exception as e:
        logging.warning(f"Error processing file {filename}: {str(e)}")

def save_to_csv(unique_sample_dataframes, directory):
    logging.info("Saving CSVs")
    for unique_sample_name, data in tqdm(unique_sample_dataframes.items(), desc="Saving CSVs"):
        df = pd.DataFrame(data={'Forward': data['seqs_forward'], 'Reverse': data['seqs_reverse']}, index=data['ids'])
        store_dir = Path(directory) / 'MED_SAMPLES_CSV' / unique_sample_name
        store_dir.mkdir(parents=True, exist_ok=True)
        save_file = store_dir / 'output.csv'
        df = df.dropna()
        df = df.sample(frac=1)  # randomize
        df.to_csv(save_file, index=False)

def main():
    setup_logging()

    if len(sys.argv) < 2:
        exit_with_error("Please provide the directory name as an argument.")

    directory_name = sys.argv[1]
    fastq_dir = Path(f'/store/sdsc/sd29/med_data_wp3/{directory_name}')

    excel_file = locate_excel_file(fastq_dir)
    if excel_file is None:
        exit_with_error("Excel file not found. Exiting.")

    reference_df = load_metadata(excel_file)

    filename_list_path = f'/users/llampert/Data_Med/Fieldworks_refs/{directory_name}.txt'
    all_filenames = read_filename_list(filename_list_path)
    num_files = len(all_filenames)

    unique_sample_dataframes = {}

    for i in tqdm(range(num_files)):
        process_file(fastq_dir, all_filenames[i], reference_df, ['other', 'OTHER', 'Other'], unique_sample_dataframes)

    save_to_csv(unique_sample_dataframes, fastq_dir)

if __name__ == "__main__":
    main()
