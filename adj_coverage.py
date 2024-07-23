import ribopy
from ribopy import Ribo
from functions import get_cds_range_lookup, get_asite_offset, get_psite_offset
import numpy as np
import multiprocessing
import time
import pickle
import logging
import math
import gzip

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
print("yipee")

# Function to process a single transcript and return coverage
def process_transcript(transcript, exp, min_len, max_len, alias, cds_range, offset, ribo_path):
    try:
        # Initialize a new Ribo object within the worker process
        if alias == True:
            ribo_object = Ribo(ribo_path, alias=ribopy.api.alias.apris_human_alias)
        else:
            ribo_object = Ribo(ribo_path)
        
        start, stop = cds_range[transcript]

        coverages = []
        for i in range(min_len, max_len + 1): 
            if offset[i] <= start:
                coverage = ribo_object.get_coverage(experiment=exp, range_lower=i, range_upper=i, alias=alias)\
                           [transcript][start - offset[i] : stop - offset[i]]
                coverages.append(coverage)
            else:
                coverage = ribo_object.get_coverage(experiment=exp, range_lower=i, range_upper=i, alias=alias)\
                           [transcript][: stop - offset[i]]
                coverage = np.concatenate((np.zeros(offset[i] - start), coverage))
                coverages.append(coverage)
                
        coverage = sum(coverages, np.zeros_like(coverages[0]))
        
        return transcript, coverage
    except Exception as e:
        logging.error(f"Error processing transcript {transcript}: {e}")
        return transcript, None

# Wrapper function to pass to multiprocessing.Pool.imap_unordered
def process_wrapper(args):
    return process_transcript(*args)

if __name__ == '__main__':
    try:
        # Initialize variables
        min_len = int(input('Enter minimum read length to be analyzed: '))
        max_len = int(input('Enter maximum read length to be analyzed: '))
        alias_int = int(input('Enter 1 for mouse or 2 for other: '))
        ribo_path = input('Enter ribo file path, e.g., \'/home/all.ribo\': ')
        offset_mode = int(input('Enter 1 for P-site Offset, Enter 2 for A-site Offset'))
        if alias_int == 1:
            alias = True
            ribo_object = Ribo(ribo_path, alias=ribopy.api.alias.apris_human_alias)
        else: 
            alias = False
            ribo_object = Ribo(ribo_path)
        cds_range = get_cds_range_lookup(ribo_object)

        all_coverage_dict = {}
        for exp in ribo_object.experiments:
            logging.info(f"Starting {exp}...")
            coverage_dict = {}
            if offset_mode == 1:
                offset = get_psite_offset(ribo_object, exp, min_len, max_len)
            if offset_mode == 2:
                offset = get_asite_offset(ribo_object, exp, min_len, max_len)

            # Parallelize transcript processing
            with multiprocessing.Pool() as pool:
                for transcript, coverage in pool.imap_unordered(
                        process_wrapper,
                        [(t, exp, min_len, max_len, alias, cds_range, offset, ribo_path) for t in ribo_object.transcript_names]
                    ):
                    # Accumulate the coverage for each transcript in the dictionary
                    coverage_dict[transcript] = coverage
            all_coverage_dict[exp] = coverage_dict

        output_file = 'coverage.pkl.gz'
        with gzip.open(output_file, 'wb') as f:
            pickle.dump(all_coverage_dict, f)

        logging.info(f"Saved as {output_file}.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")