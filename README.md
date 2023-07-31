# Analysis of transcriptomic phenotypes of lung macrophages in severe COVID-19 and SARS-CoV-2 stimulated primary monocytes.

This repository contains the code used to analyze the phenotype of macrophages from bronchoalveolar lavage (BAL) fluid as well as primary monocytes stimulated with, among others, SARS-CoV-2 as shown in 

[Wendisch, D., Dietrich, O., Mari, T., von Stillfried, S., SARS-CoV-2 infection triggers profibrotic macrophage responses and lung fibrosis, _Cell_ (2021), doi: https://doi.org/10.1016/j.cell.2021.11.033.](https://www.cell.com/cell/fulltext/S0092-8674(21)01383-0)

## Data accessibility
Count matrices (matrix.mtx, features.tsv, barcodes.tsv) as well as Seurat objects are stored on the DESY cloud service

> https://syncandshare.desy.de/index.php/s/c4kdBJaBjRSGLR7

Count matrices are downloaded directly via R as part of the analysis workflow, manual download is not necessary. 

### Structure 

* `*-counts.h5` - Count matrix in HDF5 format
* `*-counts` - Directory containing
   * `matrix.mtx`
   * `features.tsv`
   * `barcodes.tsv`
* `monocyte-clusters-*.tsv` - Souporcell clusters to demultiplex samples (Numbers on file match barcodes)
* `*.Rds` - Seurat object with analysis intermediates (normalized data, embeddings, annotations)

## Analysis workflow
To run the analyses please run the following steps:

1. Clone the git repository
   ```
   git clone https://github.com/saliba-lab/covid19-profibrotic-macrophages.git
   cd covid19-profibrotic-macrophages
   ```
1. Install packages via conda
   ```
   conda env create -f envs/default.yml
   conda activate covid19-profibrotic-macrophages
   ```
1. Run analyses
   ```
   Rscript R/BAL-macrophages.R
   Rscript R/Monocytes.R
   ```
   
License: [MIT](https://github.com/OliverDietrich/COVID-19_profibrotic-macrophage-responses/blob/main/LICENSE)
