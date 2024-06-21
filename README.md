# RiboPy Analysis

This workflow analysis is based on the [RiboFlow](https://github.com/ribosomeprofiling/riboflow) and [RiboPy](https://github.com/ribosomeprofiling/ribopy) ecosystem. Ribosome profiling data is saved as a ribo file using RiboFlow, which can then be read using the RiboPy Python environment. This workflow utilizes the capabilities of this ecosystem to:
1) Create a gzipped pickle file containing coverage data for the read lengths to be analyzed, offset to the P-site for subsequent analysis.
2) Quantify the cumulative number of footprints at each codon within the codon sequence across all transcripts.
3) Identify potential stall sites during translation, as well as their motif biases.

# Installation

Prerequisites:
* A ribo file generated using [RiboFlow](https://github.com/ribosomeprofiling/riboflow)
  * An example ribo file for _C. elegans_ is found in `example_data_c_elegans` as `all.ribo`.
* A reference file.
  * [Yeast](https://github.com/ribosomeprofiling/yeast_reference)
  * [Mouse](https://github.com/ribosomeprofiling/mouse_reference)
  * Reference file for _C. elegans_ can be found in `example_data_c_elegans` as `appris_celegans_v1_selected_new.fa`.
* Python
 
Clone repository:
```
git clone https://github.com/reikostachibana/ribopy_analysis
```

Install libraries:
```
pip install -r requirements.txt
```

# 1. Generate gzipped pickle file

The script `adj_coverage.py` processes ribosome profiling data from a Ribo file. It performs P-site offsetting based on selected read lengths and outputs the adjusted coverage data into a gzipped pickle file for efficient storage and retrieval for subsequent analysis. 

Run the following command:
```
python3 adj_coverage.py
```

Input (prompted by the console): 
* Minimum read length of the coverage data of read lengths to analyze. 
* Maximum read length of the coverage data of read lengths to analyze.
  * For the _C. elegans_ example, the minimum read length was set to 29 and the maximum to 33, which was determined experimentally based on the RPF length distribution of the coding region.
* Organism: For mouse, input 1. For other organisms, input 2.
  * Example input for _C. elegans_: `2`.
* Ribo file path.
  * Example input for the _C. elegans_ ribo file path: `./example_data_c_elegans/all.ribo`.
 
Example input sequence:
```
Enter minimum read length to be analyzed: 29
Enter maximum read length to be analyzed: 33
Enter 1 for mouse or 2 for other: 2
Enter ribo file path, e.g., '/home/all.ribo': ./example_data_c_elegans/all.ribo
```

As each experiment is processed, the console should give a confirmation. Example:
```
2024-06-13 19:09:57,550 - DEBUG - Creating converter from 3 to 5
2024-06-13 19:09:57,598 - INFO - Starting CHX_rep1...
2024-06-13 20:05:22,385 - INFO - Starting CHX_rep2...
2024-06-13 20:05:16,725 - INFO - Starting CHX_rep3...
```

Output:
* Gzipped pickle file containing a dictionary of the adjusted coverage data: {Experiment : {Transcript : Adjusted coverage array}}.
  * This is automatically saved as `coverage.pkl.gz` in the working directory.
  * The gzipped pickle file will be used in subsequent analysis.

Error handling:
* Upon successful completion, you should see a message in the console: `Saved as coverage.pkl.gz`.
* Check the directory to confirm the presence of the `coverage.pkl.gz` file.
  * If using the example _C. elegans_ Ribo file, compare with the example gzipped pickle file in `example_data_c_elegans` saved as `example_coverage.pkl.gz`.
* The script's runtime may vary depending on the size of the ribo file and the number of transcripts. For large datasets, this process may take hours.
  * For the example _C. elegans_ dataset, each experiment took around an hour to process, totaling to around six hours. This will vary based on your system's hardware and configuration.
* If you encounter errors, check the console for detailed messages. Ensure all inputs are correct and that you have the necessary permissions to read the ribo file and write to the directory.

More about the script:


# 2. Codon occupancy

Run the following command:
```
python3 codon_occupancy.py
```

Input:
* Ribo file path.
  * Example input for _C. elegans_: `./example_data_c_elegans/all.ribo`.
* Pickle file path.
  * Example input for _C. elegans_: `./example_data_c_elegans/coverage.pkl.gz`. 
* Reference file path.
  * Example input for _C. elegans_: `./example_data_c_elegans/appris_celegans_v1_selected_new.fa`.
* Organism: 1 for mouse, 2 for other.
  * Example input for _C. elegans_: `2`.

Output:
* A CSV file containing the raw counts of codons at the P-site of all footprints within the coding region across all transcripts. Automatically saved as `codon_occupancy.csv` in the working directory.
  * Rows: codons
  * Columns: experiments
    * The first column represents the cumulative total of codons in the coding sequence across all transcripts.
  * Values: raw counts of codons at the P-site of all transcripts within the coding region across all transcripts.
 
Error handling:
* Upon successful completion, you should see a message in the console: `Saved as codon_occupancy.csv`.
  * If using the example _C. elegans_ dataset, there is an example CSV file `codon_occupancy.csv` in the folder `example_output_c_elegans`.

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
