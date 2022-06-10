# Global varibles --------------------------------------------------------------

dest <- "data/monocytes.Rds"
raw <- "data/raw/monocytes/"

# Load data --------------------------------------------------------------------

list.files(raw)

# Matrix
ds <- as(Matrix::readMM(paste0(raw, "matrix")), "dgCMatrix")
colnames(ds) <- readLines(paste0(raw, "barcodes"))
features <- read.table(paste0(raw, "features"), sep = "\t")
names(features) <- c("id", "name", "type")
rownames(ds) <- features$name

# Donors
donors <- list(
  read.table(paste0(raw, "donor_1"), header = TRUE),
  read.table(paste0(raw, "donor_2"), header = TRUE)
)
donors[[2]]$barcode <- stringr::str_replace(
  donors[[2]]$barcode, "1", "2"
)
donors <- dplyr::bind_rows(donors)

# Split matrix by count type ---------------------------------------------------

# Hashtag oligos
index <- which(
  stringr::str_detect(features$type, "Antibody Capture") & 
    stringr::str_detect(features$id, "Hashtag")
)
hto <- ds[index, ]

# CITEseq antibodies
index <- which(
  stringr::str_detect(features$type, "Antibody Capture") & 
    !stringr::str_detect(features$id, "Hashtag")
)
adt <- ds[index, ]

# Viral mRNAs
index <- which(
  stringr::str_detect(features$type, "Gene Expression") & 
    stringr::str_detect(features$id, "SCoV2")
)
scov2 <- ds[index, ]

# Host mRNAs
index <- which(
  stringr::str_detect(features$type, "Gene Expression") & 
    !stringr::str_detect(features$id, "SCoV2")
)
ds <- ds[index, ]

# Create Seurat object ---------------------------------------------------------

ds <- Seurat::CreateSeuratObject(
  counts = ds, project = "monocytes", assay = "RNA"
)

# Add other assays
ds[["ADT"]] <- Seurat::CreateAssayObject(adt)
ds[["HTO"]] <- Seurat::CreateAssayObject(hto)
ds[["SCoV2"]] <- Seurat::CreateAssayObject(scov2)

# Add meta data
all(donors$barcode == colnames(ds))
ds$donor <- factor(as.character(c(
  "0" = "A", "1" = "B", "1/0" = "Doublet", "0/1" = "Doublet"
)[donors$assignment]))

# Store features table in object
ds@misc$features <- features

# Save data set ----------------------------------------------------------------
saveRDS(ds, dest)
