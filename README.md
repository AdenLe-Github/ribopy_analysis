Clone repository:

```
git clone https://github.com/reikostachibana/ribopy_analysis
```

# Install libraries:

```
pip install -r requirements.txt
```

# 1. Generate pickle file



Run the following command:
```
python3 adj_coverage.py
```

Input: 
* Organism: For mouse, input 1. For other organisms, input 2.
* Minimum read length of the coverage data of read lengths to analyze. 
* Maximum read length of the read lengths
  * For example, to analyze standard footprints of 28 nucleotides, I chose a range of 26-30 nucleotides.
* Ribo file path

Output
* Gzipped pickle file containing a dictionary of the adjusted coverage data: {Experiment : {Transcript : Adjusted coverage array}}.

# 2. Codon occupancy

Run the following command:
```
python3 codon_occupancy.py
```

Input:
* Ribo file path
* Pickle file path
* Reference file path
* Organism

Output:
* CSV file containing the raw counts of codons of all footprints

# 3. Codon heatmap

Run the following command:
```
python3 codon_heatmap_v4.py
```

Input:
* Ribo file
* Pickle file
* Reference file
* Organism
* Experiments
* Number of transcripts
* Percentile to determine threshold for stall sites

Output
* Codon heatmaps
* Excel worksheet of 1) raw counts of the stall site codons and 2) the stall sites' transcripts, positions, and nucleotide sequences
