import re
import sys
import logging
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
        reference_df = pd.read_excel(excel_file, dtype=str)
        reference_df.columns = reference_df.columns.str.upper()
        return reference_df
    except Exception as e:
        exit_with_error(f"Error loading metadata: {str(e)}")

def extract_base_sample_name(sample_name):
    return re.sub(r'_\d+$', '', sample_name)

def process_file(directory, filename, reference_df, unique_sample_dataframes):
    logging.info("Processing: %s", filename)

    # Extract the RUN name from the filename using a regular expression
    run_name_match =  filename[:8]
    if not run_name_match:
        logging.warning(f"Skipping {filename}: RUN name not found in filename.")
        return

    run_name = run_name_match.group(1)
    logging.info("Extracted RUN name: %s", run_name)

    # Filter metadata based on the extracted RUN name
    run_metadata = reference_df[reference_df['RUN'].str.contains(run_name, case=False, na=False)]
    logging.info("Run metadata: %s", run_metadata)

    if run_metadata.empty:
        logging.warning("No metadata found for the current RUN.")
        return

    # Iterate through each row in the filtered metadata for the current RUN
    for idx, row in run_metadata.iterrows():
        sample_name = row['SAMPLE']
        tag = row['TAG']
        unique_sample_name = extract_base_sample_name(sample_name)

        if unique_sample_name not in unique_sample_dataframes:
            unique_sample_dataframes[unique_sample_name] = {'tags': {tag: {'seqs_forward': [], 'seqs_reverse': []}}}
        elif tag not in unique_sample_dataframes[unique_sample_name]['tags']:
            unique_sample_dataframes[unique_sample_name]['tags'][tag] = {'seqs_forward': [], 'seqs_reverse': []}

        try:
            # Process the forward FASTQ file
            forward_file = directory / f"{filename[:-12]}_R1.fastq.gz"
            with gzip.open(forward_file, 'rt') as forward_handle:
                for forward_record in SeqIO.parse(forward_handle, format='fastq'):
                    forward_id = forward_record.id
                    forward_seq = forward_record.seq.lower()

                    # Check if the tag is in the ID
                    if tag in forward_id:
                        unique_sample_dataframes[unique_sample_name]['tags'][tag]['seqs_forward'].append(str(forward_seq))

            # Process the reverse FASTQ file
            reverse_file = directory / f"{filename[:-12]}_R2.fastq.gz"
            with gzip.open(reverse_file, 'rt') as reverse_handle:
                for reverse_record in SeqIO.parse(reverse_handle, format='fastq'):
                    reverse_id = reverse_record.id
                    reverse_seq = reverse_record.seq.lower()

                    # Check if the tag is in the ID
                    if tag in reverse_id:
                        unique_sample_dataframes[unique_sample_name]['tags'][tag]['seqs_reverse'].append(str(reverse_seq))

        except Exception as e:
            logging.warning(f"Error processing file {filename}: {str(e)}")

    logging.info("Finished processing: %s", filename)

def save_to_csv(unique_sample_dataframes, directory_name):
    logging.info("Saving CSVs")

    store_dir = Path('/scratch/snx3000/llampert/MED_SAMPLES_CSV') / directory_name
    store_dir.mkdir(parents=True, exist_ok=True)

    for unique_sample_name, data in tqdm(unique_sample_dataframes.items(), desc="Saving CSVs"):
        for tag, tag_data in data['tags'].items():
            try:
                # Create a DataFrame with forward and reverse sequences
                combined_df = pd.DataFrame(data={'Forward': tag_data['seqs_forward'], 'Reverse': tag_data['seqs_reverse']})

                # Debugging statements
                logging.info(f"Combined CSV DataFrame size for {unique_sample_name} ({tag}): {combined_df.shape}")

                # Save the combined DataFrame to a CSV file
                save_file = store_dir / f'{unique_sample_name}_{tag}.csv'

                # Debugging statements
                logging.info(f"Saving combined CSV for {unique_sample_name} ({tag}) to {save_file}")

                combined_df = combined_df.dropna()
                combined_df = combined_df.sample(frac=1)  # randomize
                combined_df.to_csv(save_file, index=False)
                logging.info(f"Combined CSV saved for {unique_sample_name} ({tag})")

            except Exception as e:
                logging.error(f"Error saving combined CSV for {unique_sample_name} ({tag}): {str(e)}")

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
    logging.info(f"All Filenames: {all_filenames}, num_files: {num_files}")

    for _, row in reference_df.iterrows():
        sample_name = row['SAMPLE']
        tag = row['TAG']
        unique_sample_name = extract_base_sample_name(sample_name)

        if unique_sample_name not in unique_sample_dataframes:
            unique_sample_dataframes[unique_sample_name] = {'tags': {tag: {'seqs_forward': [], 'seqs_reverse': []}}}
        elif tag not in unique_sample_dataframes[unique_sample_name]['tags']:
            unique_sample_dataframes[unique_sample_name]['tags'][tag] = {'seqs_forward': [], 'seqs_reverse': []}

    logging.info("Initialized unique_sample_dataframes with all samples from reference_df.")

    for i, filename in tqdm(enumerate(all_filenames), desc="Processing files", total=num_files):
        if "_R1.fastq.gz" not in filename:
            continue

        logging.info(f"Processing file ({i+1}/{num_files}): {filename}")
        process_file(fastq_dir, filename, reference_df, unique_sample_dataframes)

    save_to_csv(unique_sample_dataframes, directory_name)

if __name__ == "__main__":
    main()
