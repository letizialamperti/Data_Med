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

def read_filename_list(filename):
    with open(filename, 'r') as file:
        return [line.strip() for line in file]

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

def extract_base_sample_name(sample_name):
    return re.sub(r'_\d+$', '', sample_name)


def process_file(directory, filename, reference_df, samples_to_exclude, unique_sample_dataframes):
    logging.info("Processing: %s", filename)

    forward_file = directory / f"{filename[:-12]}_R1.fastq.gz"
    reverse_file = directory / f"{filename[:-12]}_R2.fastq.gz"

    if not forward_file.exists() or not reverse_file.exists():
        error_message = f"Error: One or both of the FASTQ files do not exist for {filename}."
        logging.error(error_message)
        raise FileNotFoundError(error_message)

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

    for idx, row in run_metadata.iterrows():
        # Get the sample name and the respective tag from each row
        sample_name = row['SAMPLE']
        tag = row['TAG']

        # Use the mapping to get the cleaned-up version
        unique_sample_name = extract_base_sample_name(sample_name)

        if any(unique_sample_name.lower().startswith(prefix) for prefix in samples_to_exclude):
            continue

        if unique_sample_name not in unique_sample_dataframes:
            unique_sample_dataframes[unique_sample_name] = {'tags': {tag: {'seqs_forward': [], 'seqs_reverse': []}}}
        elif tag not in unique_sample_dataframes[unique_sample_name]['tags']:
            unique_sample_dataframes[unique_sample_name]['tags'][tag] = {'seqs_forward': [], 'seqs_reverse': []}

        try:
            with gzip.open(forward_file, 'rt') as handle:
                for record in SeqIO.parse(handle, format='fastq'):
                    id = record.id
                    seq = record.seq.lower()
                    if any(tag in id for tag in [tag]):
                        unique_sample_dataframes[unique_sample_name]['tags'][tag]['seqs_forward'].append(str(seq))

            with gzip.open(reverse_file, 'rt') as handle:
                for record in SeqIO.parse(handle, format='fastq'):
                    id = record.id
                    seq = record.seq.lower()
                    if any(tag in id for tag in [tag]):
                        unique_sample_dataframes[unique_sample_name]['tags'][tag]['seqs_reverse'].append(str(seq))
                    else:
                        warnings.warn("ID of reverse read could not be matched to any forward read.")

        except Exception as e:
            logging.warning(f"Error processing file {filename}: {str(e)}")

    logging.info("Finished processing: %s", filename)


def save_to_csv(unique_sample_dataframes, directory):
    logging.info("Saving CSVs")
    
    for unique_sample_name, data in unique_sample_dataframes.items():
        try:
            # Extract data from the 'tags' dictionary
            tags_data = data.get('tags', {})
            
            # Create an empty DataFrame to store combined data from all tags
            combined_df = pd.DataFrame(columns=['Forward', 'Reverse'])

            # Iterate through each tag and append data to the combined DataFrame
            for tag, tag_data in tags_data.items():
                df = pd.DataFrame(data={'Forward': tag_data['seqs_forward'], 'Reverse': tag_data['seqs_reverse']})
                combined_df = combined_df.append(df, ignore_index=True)

            # Debugging statements
            logging.info(f"Combined CSV DataFrame size for {unique_sample_name}: {combined_df.shape}")

            # Save the combined DataFrame to a CSV file
            store_dir = Path('/scratch/snx3000/llampert/MED_SAMPLES_CSV') / directory
            store_dir.mkdir(parents=True, exist_ok=True)
            save_file = store_dir / f'{unique_sample_name}.csv'

            # Debugging statements
            logging.info(f"Saving combined CSV for {unique_sample_name} to {save_file}")

            combined_df = combined_df.dropna()
            combined_df = combined_df.sample(frac=1)  # randomize
            combined_df.to_csv(save_file, index=False)
            logging.info(f"Combined CSV saved for {unique_sample_name}")
        except Exception as e:
            logging.error(f"Error saving combined CSV for {unique_sample_name}: {str(e)}")


def main():
    setup_logging()

    if len(sys.argv) < 2:
        exit_with_error("Please provide the directory name as an argument.")

    directory_name = sys.argv[1]
    fastq_dir = Path(f'/store/sdsc/sd29/med_data_wp3/{directory_name}')
    store_dir = Path('/scratch/snx3000/llampert/MED_SAMPLES_CSV')

    excel_file = locate_excel_file(fastq_dir)
    if excel_file is None:
        exit_with_error("Excel file not found. Exiting.")

    reference_df = load_metadata(excel_file)
  

    filename_list_path = f'/users/llampert/Data_Med/Fieldworks_refs/{directory_name}.txt'
    all_filenames = read_filename_list(filename_list_path)

    # Define num_files after reading the filenames
    num_files = len(all_filenames)

    unique_sample_dataframes = {}
    logging.info(f"All Filenames: {all_filenames}, num_files : {num_files}")

    for i, filename in tqdm(enumerate(all_filenames), desc="Processing R1 files", total=num_files):
        if "_R1.fastq.gz" not in filename:
            # Skip files that are not R1 files
            continue

        logging.info(f"Processing R1 file ({i+1}/{num_files}): {filename}")
        process_file(fastq_dir, filename, reference_df, ['other', 'OTHER', 'Other'], unique_sample_dataframes)

    save_to_csv(unique_sample_dataframes, directory_name)

if __name__ == "__main__":
    main()
