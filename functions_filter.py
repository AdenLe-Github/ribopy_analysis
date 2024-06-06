import pickle
import gzip
from scipy.stats import zscore

def get_filtered_transcripts(pkl_gz_path, experiments, top_n):
    """
    Finds top_n number of highly expressed transcripts based on the coverage density, calculated as the number of reads per nucleotide.
    If more than one experiment is given, the intersection of these filtered transcripts of the given experiments will be returned.

    Parameters:
        pkl_gz_path (str): File path to the gzipped pickle file generated using adj_coverage.py.
        experiments (str): Experiments used in analysis, specifically in codon_heatmaps.py.
        top_n (int): Number of transcripts with the highest coverage density.

    Returns:
        list: A list of transcripts with the highest coverage density. 
    """
    with gzip.open(pkl_gz_path, 'rb') as f:
        coverage_dict = pickle.load(f)

    all_transcript_list = []
    for exp in experiments:
        exp_transcripts = [(transcript, sum(coverage) / len(coverage)) for transcript, coverage in coverage_dict[exp].items() if coverage is not None]
        exp_transcripts.sort(key=lambda x: x[1], reverse=True)  # Sort transcripts based on density
        top_transcripts = [transcript for transcript, density in exp_transcripts[:top_n]]  # Select top n transcripts
        all_transcript_list.append(set(top_transcripts))  # Convert to set for efficient intersection
    
    if all_transcript_list:
        filtered_transcripts = set.intersection(*all_transcript_list)
    else:
        filtered_transcripts = set() 
    
    return filtered_transcripts

def get_filtered_zscores(pkl_gz_path, transcripts, exp):
    """
    Normalizes coverage data into z-scores for the given transcripts for the given experiment.

    Parameters:
        pkl_gz_path (str): File path to the gzipped pickle file generated using adj_coverage.py.
        Transcripts (list): List of transcripts used in analysis, specifically in codon_heatmaps.py.
        exp (str): Experiment name.

    Returns:
        dict: A dictionary mapping transcript to z-scores of coverage data. 
    """
    with gzip.open(pkl_gz_path, 'rb') as f:
        coverage_dict = pickle.load(f)
    zscores = {transcript: zscore(coverage) for transcript, coverage in coverage_dict[exp].items() 
               if coverage is not None and transcript in transcripts}
    return zscores