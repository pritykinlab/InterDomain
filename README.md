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
      
## Usage

### Intra and Inter-chromosomal Metadomain Calling

- To perform intra-chromosomal metadomain calling, use the call_metadomains_intra command:
  - call_metadomains_intra path/to/your_file.cool 
- To perform inter-chromosomal metadomain calling, use the call_metadomains_inter command:
  - call_metadomains_inter path/to/your_file.cool 

The commands are very similar with slight modifications to the hyperparameters. The following are the optional arguments for both commands:

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
