# HuntOmics Explorer

HuntOmics Explorer is an R Shiny application developed for the comprehensive molecular analysis of Huntington's Disease (HD). It is specifically designed to analyze and visualize data from the study titled "mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals" ([GSE64810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810)).

## Features

The application offers the following functionalities:

1. **Sample Information Exploration**: Provides an overview of sample metadata, including summary statistics and visualizations.

2. **Counts Matrix Exploration**: Allows users to explore gene expression counts, visualize scatter plots, clustered heatmaps, and perform Principal Component Analysis (PCA).

3. **Differential Expression Analysis**: Enables identification of differentially expressed genes between HD and control samples, with visualization options.

4. **Gene Set Enrichment Analysis**: Assists in identifying enriched pathways and biological processes associated with HD.

## Installation

To run HuntOmics Explorer locally, follow these steps:

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/N3ha-Rao/HuntOmics-Explorer.git
   ```

2. **Navigate to the Directory**:

   ```bash
   cd HuntOmics-Explorer
   ```

3. **Install Required Packages**:

   Ensure you have R and RStudio installed. Then, install the necessary packages.
   

5. **Run the Application**:

   ```R
   library(shiny)
   runApp("N3ha-Rao-App.R")
   ```

## Usage

Upon launching the application, you will encounter several tabs, each designed for specific analyses:

1. **Sample Information Exploration**:

   - **Input**: Upload the `Metadata.csv` file containing sample information.
   - **Features**:
     - Displays a summary table of the metadata.
     - Generates histograms for numeric variables.

2. **Counts Matrix Exploration**:

   - **Input**: Upload the `Counts.csv` file containing gene expression counts.
   - **Features**:
     - Provides a summary of the counts matrix.
     - Offers scatter plots for user-selected genes.
     - Displays clustered heatmaps of filtered counts.
     - Performs PCA and visualizes the results.

3. **Differential Expression Analysis**:

   - **Input**: Upload the `DE.csv` file containing differential expression results.
   - **Features**:
     - Displays a table of differentially expressed genes.
     - Generates volcano plots to visualize gene expression changes.

4. **Gene Set Enrichment Analysis**:

   - **Input**: Uses differential expression results from the Differential Expression Analysis tab.
   - **Features**:
     - Identifies enriched Gene Ontology (GO) terms and pathways.
     - Displays bar plots and dot plots of enriched terms.

## Data Preparation

Ensure that your input files (`Metadata.csv`, `Counts.csv`, and `DE.csv`) are formatted correctly and correspond to the data from the GSE64810 study. The application expects specific column names and data structures as outlined in the `N3ha-Rao-App.R` script.

## Acknowledgements

Special thanks to the authors of the GSE64810 study for providing the data utilized in this application.
