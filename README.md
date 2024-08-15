# Bachelor Project: Comparison of single cell RNA sequencing data integration methods with application to breast cancer data

This repository contains the code for my Bachelor thesis project, developed in R using the Seurat toolkit, as well as my Bachelor thesis. The project compares four different data integration methodologies using gene expression profiles from various breast cancer samples. 

## Project Overview

The goal of this project is to objectively evaluate the qualitative aspects of different data integration methods. I use a panel of statistical and computational metrics to analyze gene expression profiles from breast cancer samples that vary in severity and metastasis. Additionally, I have developed a function to benchmark these integration procedures and provide insightful visualizations.

## Installation

To run this project, ensure you have R installed. You can install the required packages using the following command:

```R
install.packages(c("Seurat", "readxl", "harmony", "Rcpp", "data.table", "rhdf5", "MetBrewer", "RColorBrewer", "patchwork", "lubridate", "forcats", "stringr", "dplyr", "purrr", "readr", "tibble", "tidyverse", "tidyr", "TSCAN", "destiny", "slingshot", "TrajectoryUtils", "princurve", "cluster", "pls", "monocle3", "SingleCellExperiment", "SummarizedExperiment", "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "MatrixGenerics", "matrixStats", "Biobase", "BiocGenerics", "lisi", "factoextra", "magrittr", "aricode", "kBET", "SeuratDisk", "ggplot2", "SeuratWrappers", "SeuratData", "SeuratObject"))
```

## Main Dependencies

Seurat: 4.3.0
harmony: 0.1.1
tidyverse: 2.0.0
monocle3: 1.3.1

## Data

Due to size of the datasets, the data will not be uploaded. The two datasets are available in the following links:

- [GEO GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)
- [GEO GSE140819](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140819)


## License

[MIT](https://choosealicense.com/licenses/mit/)
