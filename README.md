# Workflow README

## Overview
This workflow is designed to process and analyze metagenomic assemblies, specifically focusing on human mitochondrial DNA (mtDNA). It includes downloading assemblies, extracting data, performing BLAST searches, generating variant call files (VCFs), and conducting haplogroup classification. The final output provides summary statistics for all processed samples.

The metagenomic assemblies were carried out by Pasolli et al. (2019 - https://doi.org/10.1016/j.cell.2019.01.001) and the data are publicly available (http://segatalab.cibio.unitn.it/data/Pasolli_et_al.html).

## Prerequisites
Before running the workflow, ensure that the following dependencies are installed:
- **Python**: Required for processing input files and running scripts.
- **Snakemake**: Workflow management system.
- **BLAST+**: For creating and using BLAST databases.
- **BWA**: For sequence alignment.
- **Samtools**: For processing BAM files.
- **BCFtools**: For variant calling.
- **HaploGrep3**: For haplogroup classification.
- **SeqKit**: For generating sequence statistics.

## Input Files
1. **NCBIaccession.txt**: A tab-separated file containing `sampleID` in the format `study__sample`. This file is used to extract study and sample identifiers.
2. **download_metagenomic_assemblies.sh**: A bash script contains `wget` commands to download the metagenomic assemblies of each study analyzed by Pasolli et al. (2019).
3. **Reference mtDNA File**: The workflow will download the rCRS reference genome (`rCRS.fasta`) from Phylotree.

## Output Files
The workflow generates the following outputs:

- Summary statistics for all samples:
  - Sequence stats: `summary/all_samples_rCRS_stats.tsv`
  - Haplogroup classifications: `summary/all_samples_rCRS_hg.tsv`
  - Haplogroup quality checks: `summary/all_samples_rCRS_hc.tsv`

## Usage Instructions

1. **Set Up Input File**:
   - Ensure that the file `NCBIaccession.txt` is formatted correctly with a column named `sampleID`, containing values in the format `study__sample`.

2. **Run Snakemake**:
   Execute the workflow using Snakemake:
   ```shell
   snakemake -s assembly_check.smk
   ```
   The workflow can be used with queuing systems like SLURM. In this case use (or adapt it to your SLURM profile):
   ```shell
   snakemake -s assembly_check.smk --profile slurm
   ```

3. **Check Outputs**:
   - Processed data will be available in respective directories (`contigs/`, `results/`, and `summary/`).

## Notes
- Modify resource allocations (e.g., memory or threads) in rules as needed based on your system's capabilities. All rules are set to use one thread and 10Gb of RAM.
- Ensure that external tools like HaploGrep3 and Haplocheck are properly installed and accessible in your PATH.


## Contact Information
For questions or issues regarding this workflow, please contact Mohamed S. Sarhan (mohamed.sarhan@eurac.edu).