# Global varibles --------------------------------------------------------------

dest <- "data/raw/monocytes/"
dir.create(dest, recursive = TRUE)

# Download files ---------------------------------------------------------------

files <- list(
  barcodes = "https://nubes.helmholtz-berlin.de/s/kr3pDEnTs8Zig59/download",
  features = "https://nubes.helmholtz-berlin.de/s/kDDPYYTjj7CMfEB/download",
  matrix   = "https://nubes.helmholtz-berlin.de/s/ecHWyJy3eXxgYYc/download",
  donor_1  = "https://nubes.helmholtz-berlin.de/s/xsYCzt6wpMDZafA/download",
  donor_2  = "https://nubes.helmholtz-berlin.de/s/FoGG7wAJnDaZAXe/download"
)

for (i in names(files)) {
  file <- paste0(dest, i)
  download.file(files[[i]], file)
}
