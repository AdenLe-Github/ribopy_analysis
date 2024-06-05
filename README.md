# 1. Generate pickle file

Input: 
* Minimum read length
* Maximum read length
* Organism: For mouse, input 1. For other organisms, input 2.
* Ribo file path
Output
* Gzipped pickle file: '{Experiment : {Transcript : Adjusted coverage array}}'

Run the following command:
```
python3 adj_coverage.py
```

# 2. Codon occupancy

Run the following command:
```
python3 codon_occupancy.py
```

# 3. Codon heatmap

Run the following command:
```
python3 codon_heatmap_v4.py
```
