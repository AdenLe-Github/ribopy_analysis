import ribopy
from ribopy import Ribo
import pickle
import gzip
import pandas as pd
from collections import defaultdict
from functions import get_cds_range_lookup, get_sequence

ribo_path = input('Enter ribo file path, e.g., \'/home/all.ribo\': ')
coverage_path = input('Enter pickle file path, e.g., \'/home/coverage.pkl.gz\': ')
reference_file_path = input('Enter reference file path: ')
alias_int = int(input('Enter 1 for mouse or 2 for other: '))

if alias_int == 1:
    alias = True
    ribo_object = Ribo(ribo_path, alias=ribopy.api.alias.apris_human_alias)
else: 
    alias = False
    ribo_object = Ribo(ribo_path)

with gzip.open(coverage_path, 'rb') as f:
    coverage_dict = pickle.load(f)

cds_range = get_cds_range_lookup(ribo_object)
sequence = get_sequence(ribo_object, reference_file_path, alias)

df_codon_occ = pd.DataFrame()
for exp in coverage_dict.keys():
    codon_occ = defaultdict(int)
    for transcript, coverage in coverage_dict[exp].items():
        if coverage is not None:
            start, stop = cds_range[transcript]
            cds_sequence = sequence[transcript][start: stop]

            for i in range(0, len(cds_sequence), 3):
                codon = cds_sequence[i: i + 3]
                count = sum(coverage[i: i + 3])  
                codon_occ[codon] += count

    sorted_codon_occ = {k: codon_occ[k] for k in sorted(codon_occ)}
    
    df_temp = pd.DataFrame(list(sorted_codon_occ.items()), columns=['Codon', exp])
    
    if df_codon_occ.empty:
        df_codon_occ = df_temp
                        
        transcriptome_codon_dist = defaultdict(int)
        for transcript, coverage in coverage_dict[exp].items():
            if coverage is not None:
                start, stop = cds_range[transcript]
                cds_sequence = sequence[transcript][start: stop]

                for i in range(0, len(cds_sequence), 3):
                    codon = cds_sequence[i: i + 3]
                    transcriptome_codon_dist[codon] += 1
        sorted_codon_dist = {j: transcriptome_codon_dist[j] for j in sorted(transcriptome_codon_dist)}
        df_codon_dist = pd.DataFrame(list(sorted_codon_dist.items()), columns=['Codon', 'Transcriptome'])
        df_codon_occ = pd.merge(df_codon_dist, df_temp, on='Codon', how='inner')    
        
    else:
        df_codon_occ = pd.merge(df_codon_occ, df_temp, on='Codon', how='outer')

output_file = 'codon_occupancy.csv'
df_codon_occ.to_csv(output_file, index=False)
print(f'Saved as {output_file}.')
