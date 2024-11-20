# BF591-RShiny-App

**Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls**: This dataset profiled gene expression with RNASeq in post-mortem human dorsolateral prefrontal cortex from patients who died from Huntington’s Disease and age- and sex-matched neurologically healthy controls.

I engineered an R Shiny app for comprehensive molecular analysis of Huntington's Disease based on the above study, featuring sample exploration, differential expression analysis, correlation networks, and gene set enrichment, enabling direct comparisons with healthy cohorts.

## Tab 1: Sample Information Exploration

Summary table to give an overview of each column <br>
Displays a table of data and histograms of numeric data <br>

## Tab 2: Counts Matrix Exploration

Provides a summary table of counts matrix table <br>
Scatter plot to show user input filtered genes <br>
Clustered heatmap of counts post-filter <br>
User selected PCA plot <br>

## Tab 3: Differential Expression

Plot of user defined x or y axis variables with threshold colored p-adj data points <br>
Data table of resulting values after filtering by threshold <br>

## Tab 4: Gene Set Enrichment Analysis Exploration

Normalized enrichment score (NES) plot of up and down-regulated genes <br>
Table of results with filtering options and download feature <br>
NES scatterplot with threshold colored p-adj data points <br>

## File upload for each tab:
Counts.csv <br>
DE.csv <br>
Metadata.csv

Data obtained from: [mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810)
