import ribopy
from ribopy import Ribo
from functions import get_cds_range_lookup, get_sequence
from functions_filter import get_filtered_transcripts, get_filtered_zscores
import pickle
from scipy.stats import zscore
import numpy as np
import pandas as pd
from collections import defaultdict
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import gzip

def calculate_threshold(zscores, percentile):
    """
    Calculates threshold to determine stall sites based on the given percentile. 

    Parameters:
        zscores (dict): Dictionary mapping transcript to an array of the respective coverage data normalized into z-scores. 
        percentile (float): Percentile to determine threshold. In codon_heatmaps.py, this is defaulted to the 99th percentile.

    Returns:
        float: The z-score cutoff to determine stall sites.
    """
    all_zscores = [score for z_scores in zscores.values() for score in z_scores[18:-15]]
    return np.percentile(all_zscores, percentile)

def find_common_stall_sites(replicates, pickle_path, transcripts, percentile):
    """
    Finds common stall sites across replicates.

    Parameters:
        replicates (list): List of replicates.
        pickle_path (str): File path to gzipped pickle file.
        transcripts (list): List of transcripts.
        percentile (float): Percentile to determine threshold.

    Returns:
        dict: Dictionary mapping transcript to an array of booleans, each value representing the presence of a common stall site at each nucleotide position. 
              True represents stall site.
    """
    common_stall_sites = {}
    for exp in replicates:
        zscores = get_filtered_zscores(pickle_path, transcripts, exp)        
        threshold = calculate_threshold(zscores, percentile)

        for transcript, zscore in zscores.items():
            stall_sites = zscore > threshold
            if transcript not in common_stall_sites:
                common_stall_sites[transcript] = [stall_sites]
            else:
                common_stall_sites[transcript].append(stall_sites)
    
    # Intersection of stall sites
    for transcript, stall_sites_list in common_stall_sites.items():
        common_stall_sites[transcript] = np.all(stall_sites_list, axis=0)
    
    return common_stall_sites

def collect_stall_sequences(common_stall_sites, sequence, cds_range):
    """
    Finds codon sequences at and around stall sites.

    Parameters: 
        common_stall_sites (dict): Dictionary mapping transcript to an array of booleans at each nucleotide position representing stall sites from find_common_stall_sites().
        sequence (dict): Dictionary mapping transcript to nucleotide sequence from get_sequence().
        cds_range (dict): Dictionary mapping transcript to CDS range (start, stop) from get_cds_range_lookup().

    Returns:
        list: A list of stall sequences including 5 codons downstream and 5 codons upstream of the stall site codon.
    """
    stall_sequences = []
    for transcript, stall_sites in common_stall_sites.items():
        start, stop = cds_range[transcript]
        for i in range(start + 18, stop - 15, 3):
            if stall_sites[i - start:i - start + 3].any():
                stall_sequences.append(sequence[transcript][i - 15:i + 18])
    return stall_sequences

def create_raw_heatmap(stall_sequences):
    """
    Creates a DataFrame for the heatmap with the raw counts of codons at/around the stall sites. 

    Parameters:
        stall_sequences (list): List of stall sequences from collect_stall_sequences().
    
    Returns:
        DataFrame: DataFrame with columns from -5 to 5 and rows of all the codons at/around the stall sites.
    """
    df = pd.DataFrame(columns=range(-5, 6), index=range(len(stall_sequences)))
    for i, seq in enumerate(stall_sequences):
        df.loc[i] = [seq[j:j+3] for j in range(0, len(seq), 3)]
    return df.apply(pd.Series.value_counts).fillna(0)

def normalize_heatmap(raw_heatmap, sequence, cds_range, transcripts):
    """
    Creates a DataFrame for the heatmap, normalized by the total number of codons across all given transcripts within the CDS. 

    Parameters:
        raw_heatmap (DataFrame): DataFrame from create_raw_heatmap().
        sequence (dict): Dictionary mapping transcript to nucleotide sequence from get_sequence().
        cds_range (dict): Dictionary mapping transcript to CDS range (start, stop) from get_cds_range_lookup().
        transcripts (list): List of transcripts for analysis.
    
    Returns:
        DataFrame: Normalized DataFrame of raw_heatmap using total number of codons across all given transcripts within the CDS.
    """
    codon_counts = {}
    for transcript in transcripts:
        start, stop = cds_range[transcript]
        for i in range(start + 3, stop, 3):
            codon = sequence[transcript][i:i+3]

            codon_counts[codon] = codon_counts.get(codon, 0) + 1
    
    total_codons = sum(codon_counts.values())
    norm_counts = {codon: count / total_codons for codon, count in codon_counts.items()}
    # Reindex raw_heatmap to include all codons from norm_counts
    all_codons = list(norm_counts.keys())
    norm_heatmap = raw_heatmap.reindex(all_codons).fillna(0)
    norm_heatmap = norm_heatmap.div(norm_heatmap.sum(axis=0), axis=1)
    norm_heatmap = norm_heatmap.sort_index(ascending=False)
    
    for index, row in norm_heatmap.iterrows():
        norm_heatmap.loc[index] = row - norm_counts.get(index, 0)
        
    return norm_heatmap

def get_heatmap_df(raw_heatmap, sequence, cds_range, transcripts):
    """
    Retrieve DataFrame of the raw heatmap, also containing the total number of each codon within the CDS across all given transcripots.

    Parameters:
        raw_heatmap (DataFrame): DataFrame from create_raw_heatmap().
        sequence (dict): Dictionary mapping transcript to nucleotide sequence from get_sequence().
        cds_range (dict): Dictionary mapping transcript to CDS range (start, stop) from get_cds_range_lookup().
        transcripts (list): List of transcripts for analysis.
    
    Returns:
        DataFrame: DataFrame of raw counts of codons at/around stall sites. 
                   There is an additional column containing the total number of each codon within the CDS across all given transcripots.
    """
    codon_counts = {}
    for transcript in transcripts:
        start, stop = cds_range[transcript]
        for i in range(start + 3, stop, 3):
            codon = sequence[transcript][i:i+3]
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
    
    all_codons = list(codon_counts.keys())
    df = raw_heatmap.reindex(all_codons).fillna(0)
    df['Transcriptome'] = df.index.map(codon_counts)
    df = df.sort_index(ascending=True)
    return df

def get_stall_sites_df(common_stall_sites, cds_range, sequence):
    """
    Retrieves transcript, nucleotide position, and sequence of the stall sites.

    Parameters:
        common_stall_sites (dict): Dictionary mapping transcript to an array of booleans at each nucleotide position representing stall sites from find_common_stall_sites().
        cds_range (dict): Dictionary mapping transcript to CDS range (start, stop) from get_cds_range_lookup().
        sequence (dict): Dictionary mapping transcript to nucleotide sequence from get_sequence().
    """
    df_list = []
    
    for transcript, stall_sites in common_stall_sites.items():
        start, stop = cds_range[transcript]
        
        for i in range(start + 18, stop - 15, 3):
            if any(stall_sites[i - start:i - start + 3]):
                sequence_window = sequence[transcript][i - 15:i + 18]
                codons = [sequence_window[j:j + 3] for j in range(0, len(sequence_window), 3)]
                
                df_temp = pd.DataFrame({
                    'Transcript': [transcript],
                    'StallSite': [i + 1],
                    'Codons': [codons]
                })
                df_list.append(df_temp)
    
    df = pd.concat(df_list, ignore_index=True)
    return df