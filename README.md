# TEB - Tripeptide Environment Bias

Analysis of inter-residue sidechain vector angles in helical tripeptide environments.

## What this does

1. Runs STRIDE on PDB structures to assign secondary structure
2. Extracts tripeptide windows centered on a target amino acid
3. Computes signed angles between adjacent CA→sidechain centroid vectors
4. Generates density plots grouped by neighbor residue size

## Quick start

```bash
bash run.sh /path/to/pdb/files ARG
```

## Structure

```
src/            - processing scripts
pipeline/       - snakemake workflow files  
output/         - results and plots
settings.yml   - configuration
run.sh          - orchestrator
```

## Requirements

- Python 3.8+ with numpy, biopython, pandas, matplotlib, scipy, tqdm
- STRIDE (secondary structure tool)
- Snakemake
- (optional) R with ggplot2 for alternative plotting
