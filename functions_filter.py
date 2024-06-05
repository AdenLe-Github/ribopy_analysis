import pickle
import gzip
from scipy.stats import zscore

def get_filtered_transcripts(pkl_gz_path, experiments, top_n):
    with gzip.open(pkl_gz_path, 'rb') as f:
        coverage_dict = pickle.load(f)

    all_transcript_list = []
    for exp in experiments:
        exp_transcripts = [(transcript, sum(coverage) / len(coverage)) for transcript, coverage in coverage_dict[exp].items() if coverage is not None]
        exp_transcripts.sort(key=lambda x: x[1], reverse=True)  # Sort transcripts based on density
        top_transcripts = [transcript for transcript, density in exp_transcripts[:top_n]]  # Select top 200 transcripts
        all_transcript_list.append(set(top_transcripts))  # Convert to set for efficient intersection
    
    if all_transcript_list:
        filtered_transcripts = set.intersection(*all_transcript_list)
    else:
        filtered_transcripts = set()  # Return an empty set if there are no experiments
    
    return filtered_transcripts

def get_filtered_zscores(pkl_gz_path, transcripts, exp):
    with gzip.open(pkl_gz_path, 'rb') as f:
        coverage_dict = pickle.load(f)
    zscores = {transcript: zscore(coverage) for transcript, coverage in coverage_dict[exp].items() 
               if coverage is not None and transcript in transcripts}
    return zscores