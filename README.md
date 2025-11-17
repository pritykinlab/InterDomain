# InterDomain

InterDomain is a Python package for detecting metadomains in Hi-C contact matrices. It supports both intra-chromosomal and inter-chromosomal metadomain calling, allowing researchers to identify significant genomic interactions within and between chromosomes.

---

## Installation

### Requirements

- **Python**: Version 3.6 or higher
- **Dependencies**: The following Python packages are required:
  - numpy
  - pandas
  - cooler
  - scipy
  - matplotlib

### Installation Steps

**Clone the Repository**:

   ```bash
   git clone https://github.com/yourusername/InterDomain.git
   cd InterDomain
   pip install .
   ```
     
## Tutorial

This short tutorial walks through calling metadomains on a Treg dataset and plotting the top results.

1) Download the Treg `.mcool` file

- Download the Treg `.mcool` file from <LINK> and note its path on your machine. If using a multi-resolution `.mcool`, reference a specific resolution using a Cooler URI, e.g. `path/to/Treg_all.mcool::/resolutions/50000`.

2) Run intra-chromosomal metadomain calling

```bash
call_metadomains_intra path/to/Treg_all.mcool::/resolutions/50000 \
  --label Treg \
  --n_workers 8 \
  --save_intermediates \
  --filter_width 1 \
  --filter_n 15 \
  --cutoff 0
```

3) Run inter-chromosomal metadomain calling

```bash
call_metadomains_inter path/to/Treg_all.mcool::/resolutions/50000 \
  --label Treg \
  --n_workers 8 \
  --save_intermediates \
  --filter_width 3 \
  --filter_n 35
```

4) Plot the top metadomains

- Replace `--type` with `intra` or `inter` depending on which results you want to visualize.

```bash
interdomain_plot \
  --top_n 20 \
  --output_dir bedfile_output/ \
  --type inter
```

Notes:
- Adjust `--n_workers` to match your available CPU cores.
- If you used a `.cool` file instead of `.mcool`, pass its path directly (no resolution suffix).
## Usage

### Intra and Inter-chromosomal Metadomain Calling

- To perform intra-chromosomal metadomain calling, use the call_metadomains_intra command:
  - call_metadomains_intra path/to/your_file.cool 
- To perform inter-chromosomal metadomain calling, use the call_metadomains_inter command:
  - call_metadomains_inter path/to/your_file.cool 

The commands are similar for both intra and inter metadomains, with slight modifications to the hyperparameters. The following are the main arguments for both commands:

- --n_workers: Number of worker processes to use (default: 1; max: number of chromosomes).
- --label: Label for the output files (default: 'test').
- --save_intermediates: Save intermediate matrices for plotting and inspection of results (default: False).
- --useSigma: Use sigma for smoothing (intra default: False; inter default: True; recommended for sparser matrices and interchromosomal matrices).
- --sigma: Sigma value for smoothing if --useSigma is set (default: 0.75).
- --prominence: Prominence threshold for peak detection (default: 4).
- --minweight: Minimum weight for balancing (default: 0.0025).
- --pc: Pseudocount for observed/expected calculations (default: 1e-6).
- --pmin: Minimum p-value for statistical tests (default: 1e-300).
- --frac_min_valid: Minimum fraction of valid pixels for filters (default: 0).
- --filter_n: Size of the outside filter (default: 35).
- --filter_width: Size of the inside filter (intra default: 1; inter default: 3;).
- --output_dir: Directory to save output files (default: 'bedfile_output/').
- --cutoff: intra only; determines the minimum size of a metadomain (default: 2,000,000).


### Outputs

Results will be saved in the specified directory (--output_dir; default: 'bedfile_output/'). The output files will contain the metadomain coordinates and other relevant information, and are saved as .tsv files. Intra and inter-chromosomal results will be saved in separate directories (/intra/ and /inter/).

Intermediate results can be saved in --save_intermediates, which can allow for better visualization of intermediate steps in the metadomain calling process.

Intermediate files are as follows:
- Observed and expected matrices
- Balanced matrices
- Smoothed matrices
- Peak detection results
- -LogP values



### Plotting Metadomains

After running the metadomain calling, you can visualize the most significant metadomains using the interdomain_plot command. This command requires :
- --top_n: specify the number of metadomains to plot
- --output_dir: specify what the output dir was for your dataset
- --type: specify whether to plot intrachromosomal or interchromosomal metadomains
