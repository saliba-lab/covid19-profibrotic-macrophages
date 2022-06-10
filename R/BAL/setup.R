# Global varibles --------------------------------------------------------------

dest <- "data/BAL.Rds"
raw <- "data/raw/BAL/"

# Load data --------------------------------------------------------------------

list.files(raw)

# Matrix
ds <- as(Matrix::readMM(paste0(raw, "matrix")), "dgCMatrix")
colnames(ds) <- readLines(paste0(raw, "barcodes"))
features <- read.table(paste0(raw, "features"), sep = "\t")
names(features) <- c("id", "name", "type")
rownames(ds) <- features$name

# Metadata
coldata <- readxl::read_excel(paste0(raw, "coldata"))

# Split matrix by count type ---------------------------------------------------

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
  counts = ds, project = "BAL", assay = "RNA"
)

# Add other assays
ds[["SCoV2"]] <- Seurat::CreateAssayObject(scov2)

# Add meta data
index <- as.numeric(stringr::str_split(colnames(ds), "-", simplify = T)[, 2])
for (i in names(coldata)) {
  ds[[i]] <- factor(coldata[[i]][index])
}

index <- order(as.numeric(stringr::str_remove(coldata$dpso, "day ")))
ds$patient <- factor(as.character(ds$patient), coldata$patient[index])
ds$dpso <- factor(as.character(ds$dpso), unique(coldata$dpso[index]))

# Store features table in object
ds@misc$features <- features

# Save data set ----------------------------------------------------------------
saveRDS(ds, dest)


