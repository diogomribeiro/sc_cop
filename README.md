# Single cell local gene co-expression pairs (COP) project

This repository contains scripts for data processing, analysis and figure generation data for our preprint:

[Ribeiro DM, Ziyani C, Delaneau O. **Shared regulation and functional relevance of local gene co-expression revealed by single cell analysis.** bioRxiv (2021).](https://www.biorxiv.org/content/10.1101/2021.12.14.472573v1)

[![DOI](https://zenodo.org/badge/437784059.svg)](https://zenodo.org/badge/latestdoi/437784059)

## Analysis scripts
- **CODer.py** : Script to identify local co-expressed gene pairs (COPs) given a gene expression matrix. Example usage: python3 CODer.py expression_matrix.bed output_folder 1000 --fdrCutoff 0.01
- **ShareSeqCoex.py** Script to identify co-expressed gene-peak pairs from single cell a atac-seq and gene expression in same cells. Example usage: python3 ShareSeqCoex.py gene_matrix.tsv peak_matrix.tsv gencode_v19.bed output_file.txt
- ** cuomo / share_seq / sarkar ** : Folder with all scripts used to produce figures for the paper, split by dataset

## Data availability
Data on co-expressed genes discovered here are available for consultation and download through the [LoCOP DB](http://glcoex.unil.ch) database. 

## License
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
LoCOP is available under a MIT license. For more information please see the [LICENSE](LICENSE).
