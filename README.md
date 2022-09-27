# Analysis of transcriptomic phenotypes of lung macrophages in severe COVID-19 and SARS-CoV-2 stimulated primary monocytes.

This repository contains the code used to analyze the phenotype of macrophages from bronchoalveolar lavage (BAL) fluid  ([here](https://github.com/OliverDietrich/SARS-CoV-2-infection-triggers-profibrotic-macrophage-responses-and-lung-fibrosis/blob/main/R/BAL-macrophages.R)) as well as primary monocytes stimulated with, among others, SARS-CoV-2 ([here](https://github.com/OliverDietrich/SARS-CoV-2-infection-triggers-profibrotic-macrophage-responses-and-lung-fibrosis/blob/main/R/Monocytes.R)) as shown in 

[Wendisch, D., Dietrich, O., Mari, T., von Stillfried, S., SARS-CoV-2 infection triggers profibrotic macrophage responses and lung fibrosis, _Cell_ (2021), doi: https://doi.org/10.1016/j.cell.2021.11.033.](https://www.cell.com/cell/fulltext/S0092-8674(21)01383-0)

## Analysis workflow
To run the analyses please run the following steps:

   1. Clone the git repository
   ```
   git clone https://github.com/OliverDietrich/covid19.macrophages.git
   cd covid19.macrophages
   ```
   1. Install packages via conda
   ```
   conda env create -f envs/default.yml
   conda activate covid19-macrophages
   ```
   1. Run analyses
   ```
   Rscript R/BAL-macrophages.R
   Rscript R/Monocytes.R
   ```

## Data accessibility
Count matrices (matrix.mtx, features.tsv, barcodes.tsv) as well as Seurat objects are stored on the Nubes cloud service

> https://nubes.helmholtz-berlin.de/s/XrM8igTzFTFSoio

Count matrices are downloaded directly via R as part of the analysis workflow, manual download is not necessary. 

License: [MIT](https://github.com/OliverDietrich/COVID-19_profibrotic-macrophage-responses/blob/main/LICENSE)
