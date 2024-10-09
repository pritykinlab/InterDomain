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
