# sc_cop
Single cell local gene co-expression project

This repository contains scripts for data processing, analysis and figure generation data for our preprint:

Ribeiro DM, Ziyani C, Delaneau O. **Shared regulation and functional relevance of local gene co-expression revealed by single cell analysis.**(2021). 

## Analysis scripts
- **CODer.py** : Script to identify local co-expressed gene pairs (COPs) given a gene expression matrix. Example usage: python3 CODer.py expression_matrix.bed output_folder 1000 --fdrCutoff 0.01
- **manuscript_figures** : Folder with all scripts used to produce figures for the paper. This includes the code to perform logistic regression analysis.

## Data availability
Data on co-expressed genes discovered here are available for consultation and download through the [LoCOP DB](http://glcoex.unil.ch) database. 

## License
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
LoCOP is available under a MIT license. For more information please see the [LICENSE](LICENSE).
