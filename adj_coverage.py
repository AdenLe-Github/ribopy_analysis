import ribopy
from ribopy import Ribo
from functions import get_sequence, get_cds_range_lookup, get_psite_offset
import numpy as np
import multiprocessing
import time
import pickle

ribo_path = f'/home/reiko/ribopy_analysis/mouse/all.ribo'
reference_path = '/home/reiko/ribopy_analysis/mouse/appris_mouse_v2_selected.fa.gz'
exp = 'WT_control_A'
min_len = 25
max_len = 31
alias = True
ribo_object = Ribo(ribo_path, alias = ribopy.api.alias.apris_human_alias)


def get_adj_coverage(ribo_object, transcript, exp, min_len, max_len, alias):
    start, stop = get_cds_range_lookup(ribo_object)[transcript]
    offset = get_psite_offset(ribo_object, exp, min_len, max_len)

    if start < max(offset.values()):
        return None
    
    coverages = [
        ribo_object.get_coverage(experiment=exp, range_lower=i, range_upper=i, alias=alias)
        [transcript][start - offset[i] : stop - offset[i]]
        for i in range(min_len, max_len + 1)
    ]

    coverage = sum(coverages, np.zeros_like(coverages[0]))
    return coverage


def process_transcript(transcript):
    coverage = get_adj_coverage(ribo_object, transcript, exp, min_len, max_len, alias)
    return coverage

if __name__ == '__main__':
    start_time = time.time()

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    results = []
    for transcript in ribo_object.transcript_names:
        if alias == True:
            transcript = transcript.split("|")[4]
        result = pool.apply_async(process_transcript, args=(transcript,))
        results.append((transcript, result))  # Store both transcript and its result

    pool.close()
    pool.join()

    adj_coverage_dict = {}
    for transcript, result in results:
        coverage = result.get()  # Get the actual result
        adj_coverage_dict[transcript] = coverage

    with open(f'/home/reiko/ribopy_analysis/mouse/coverage/adj_coverage_{exp}.pkl', 'wb') as f:
        pickle.dump(adj_coverage_dict, f)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")
