G-Conformers

G-Conformers is a high-throughput 3D conformer generation tool for molecules represented by SMILES. It leverages RDKit and UFF force fields to generate and optimize conformers using a simple evolutionary approach. The tool supports multiprocessing and progress tracking, making it suitable for computational chemistry and drug discovery workflows.

Features

- Convert SMILES strings to 3D conformers
- Optimize using Universal Force Field (UFF)
- Select lowest energy conformer using a genetic-like algorithm
- Parallel processing with real-time progress bar (`tqdm`)
- Outputs:
  - `.sdf` file with 3D structures
  - `.csv` file summarizing conformer energy and ID

Requirements

- Python 3.7 or higher
- See `requirements.txt`

Installation

Using pip (not recommended for RDKit):
```bash
pip install -r requirements.txt

conda create -n gconformers python=3.9 -y
conda activate gconformers
conda install -c conda-forge rdkit tqdm -y
