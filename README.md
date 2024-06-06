# Generalized bayesian factor analysis for integrative clustering with applications to multi-omics data (Conference on Data Science and Advanced Analytics, 2018)
Eun Jeong Min, Changgee Chang & Qi Long

## Abstract 
We first develop a generalized Bayesian factor analysis (GBFA) framework, a
sparse Bayesian factor analysis which can take into account the graph information. Our GBFA
framework employs the spike and slab lasso (SSL) prior to impose sparsity on the factor loadings
and the Markov random field (MRF) prior to encourage smoothing over the adjacent factor
loadings, which establishes a unified shrinkage adaptive to the loading size and the graph
structure.

## Code Description
- DWL.R contains the supporting functions for GBFA.R.   
- GBFA.R contains the main function for the model BGB.
- gbfa_code.R contains the sample code for the model implementation and model selection. 

## Link
[Original paper](https://ieeexplore.ieee.org/document/8631499)
