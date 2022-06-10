# Global varibles --------------------------------------------------------------

dest <- "data/raw/BAL/"
dir.create(dest, recursive = TRUE)

# Download files ---------------------------------------------------------------

files <- list(
  barcodes = "https://nubes.helmholtz-berlin.de/s/8DYsAwEMD3CniNj/download",
  features = "https://nubes.helmholtz-berlin.de/s/rDAk9Sa7n8rGGG6/download",
  matrix   = "https://nubes.helmholtz-berlin.de/s/dS39WttqQz7qedp/download",
  coldata  = "https://nubes.helmholtz-berlin.de/s/22QDAkENpr3EoBk/download",
  genesets = "https://nubes.helmholtz-berlin.de/s/25ZGCXAHdBffZCP/download"
)

# Extend download timeout
options(timeout = Inf)

# Download iteratively
for (i in names(files)) {
  file <- paste0(dest, i)
  download.file(files[[i]], file)
}
