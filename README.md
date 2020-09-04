# Subhajit-drug-disease
 
This repository contains the code to perform semantic similarity and proximity score analyses comparing drug and disease signatures, where each signature is encoded as a list of genes. This code is still a work in progress.

## Semantic similarity

A brief description of each script/notebook.

### functions.R

Contains the underlying functions needed for the various analyses so that they can be shared across notebooks.

### SemSim_Timing.Rmd

Establish some timing benchmarks for the different semantic similarity methods so as to help determine which metric to ultimately use.


### SemSim_Null.Rmd

Generate a null distribution of semantic similarity values against which we can then compute a p-value for the true drug-disease similarities. Note that this is very computationally intensive and hence for the inital approach I've computed permutation nulls of size 100, which is fairly small. There may be ways to improve this but the current approach will give us some results to get started with.

### Semantic_Similarity.Rmd

Compute the semantic similarities between the drug and gene lists and perform some downstream analyses.

## Proximity

A brief description of each script/notebook.

### functions.R

Contains the underlying functions needed for the various analyses so that they can be shared across notebooks.

### Generate_Networks.Rmd

Wrangle the PPI and GCN background network files into graph objects that we can then compute on.

