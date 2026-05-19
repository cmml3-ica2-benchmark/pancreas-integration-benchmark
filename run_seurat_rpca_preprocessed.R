suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Read the exact notebook exports used for Seurat RPCA.
# The matrix is genes x cells, and metadata rownames must match barcodes.
input_dir <- "seurat_input"
matrix_path <- file.path(input_dir, "matrix.mtx")
genes_path <- file.path(input_dir, "genes.tsv")
cells_path <- file.path(input_dir, "barcodes.tsv")
meta_path <- file.path(input_dir, "metadata.csv")

required_inputs <- c(matrix_path, genes_path, cells_path, meta_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required Seurat RPCA inputs: ", paste(missing_inputs, collapse = ", "))
}

mat <- readMM(matrix_path)
genes <- read.table(genes_path, stringsAsFactors = FALSE)[, 1]
cells <- read.table(cells_path, stringsAsFactors = FALSE)[, 1]
meta <- read.csv(meta_path, row.names = 1)

if (anyDuplicated(genes) > 0) {
  message("Duplicate gene names detected; making them unique after replacing underscores.")
}
genes <- make.unique(gsub("_", "-", genes))

if (nrow(mat) != length(genes)) {
  stop("Gene count mismatch between matrix.mtx and genes.tsv.")
}
if (ncol(mat) != length(cells)) {
  stop("Cell count mismatch between matrix.mtx and barcodes.tsv.")
}
if (anyDuplicated(cells) > 0) {
  stop("Duplicated cell barcodes found in barcodes.tsv.")
}

required_meta_cols <- c("batch", "celltype")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
if (length(missing_meta_cols) > 0) {
  stop("metadata.csv is missing required columns: ", paste(missing_meta_cols, collapse = ", "))
}
if (!all(cells %in% rownames(meta))) {
  missing_meta_rows <- setdiff(cells, rownames(meta))
  stop("metadata.csv is missing rows for exported cells, for example: ", paste(head(missing_meta_rows, 5), collapse = ", "))
}

rownames(mat) <- genes
colnames(mat) <- cells
meta <- meta[cells, , drop = FALSE]

if (!identical(rownames(meta), cells)) {
  stop("Failed to preserve notebook cell order when aligning metadata.")
}

mat <- as(mat, "CsparseMatrix")
message("Loaded matrix: ", nrow(mat), " genes x ", ncol(mat), " cells")
message("Metadata columns: ", paste(colnames(meta), collapse = ", "))

# Remove genes with no variance globally so RPCA does not waste features on flat rows.
row_means <- Matrix::rowMeans(mat)
row_means_sq <- Matrix::rowMeans(mat^2)
row_vars <- row_means_sq - row_means^2
row_vars[is.na(row_vars)] <- 0

keep_genes <- row_vars > 0
mat <- mat[keep_genes, ]
message("Genes retained after zero-variance filter: ", nrow(mat))

# Seurat still needs a non-negative counts slot even when the notebook matrix is preprocessed.
# Negative values are clipped only in the counts layer; the original matrix is kept in the data layer.
mat_counts <- mat
if (length(mat_counts@x) > 0) {
  mat_counts@x[mat_counts@x < 0] <- 0
}

obj <- CreateSeuratObject(
  counts = mat_counts,
  meta.data = meta,
  min.cells = 0,
  min.features = 0
)

# Keep the original exported matrix in the data layer so the workflow does not renormalize it.
obj <- SetAssayData(
  object = obj,
  assay = "RNA",
  layer = "data",
  new.data = mat
)

if (!identical(colnames(obj), cells)) {
  stop("CreateSeuratObject changed the exported cell order unexpectedly.")
}

obj.list <- SplitObject(obj, split.by = "batch")
message("Batches in Seurat object: ", paste(names(obj.list), collapse = ", "))

sparse_row_var <- function(m) {
  mu <- Matrix::rowMeans(m)
  mu2 <- Matrix::rowMeans(m^2)
  v <- mu2 - mu^2
  v[is.na(v)] <- 0
  return(v)
}

# Use genes that vary in every batch so RPCA anchors are defined on shared informative features.
vars_list <- lapply(obj.list, function(x) {
  d <- GetAssayData(x, assay = "RNA", layer = "data")
  sparse_row_var(d)
})

common_features <- Reduce(intersect, lapply(vars_list, function(v) names(v[v > 0])))
avg_var <- Reduce("+", lapply(vars_list, function(v) v[common_features])) / length(vars_list)
features <- names(sort(avg_var, decreasing = TRUE))
features <- features[1:min(1000, length(features))]

message("Number of integration features: ", length(features))
if (length(features) < 50) {
  stop("Too few variable features for Seurat integration.")
}

# Match the report benchmark by keeping 20 PCs for RPCA and downstream UMAP.
npcs_use <- 20

obj.list <- lapply(obj.list, function(x) {
  VariableFeatures(x) <- features
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, npcs = npcs_use, verbose = FALSE)
  x
})

anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:npcs_use
)

combined <- IntegrateData(
  anchorset = anchors,
  dims = 1:npcs_use
)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = npcs_use, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:npcs_use)

pca_out <- Embeddings(combined, reduction = "pca")
umap_out <- Embeddings(combined, reduction = "umap")
meta_out <- combined@meta.data

# Check that the exported files keep the integrated cell order unchanged for notebook import.
if (!identical(rownames(pca_out), colnames(combined))) {
  stop("PCA export rownames do not match the integrated Seurat cell order.")
}
if (!identical(rownames(umap_out), colnames(combined))) {
  stop("UMAP export rownames do not match the integrated Seurat cell order.")
}
if (!identical(rownames(meta_out), colnames(combined))) {
  stop("Metadata export rownames do not match the integrated Seurat cell order.")
}

write.csv(pca_out, "seurat_rpca_pca.csv")
write.csv(umap_out, "seurat_rpca_umap.csv")
write.csv(meta_out, "seurat_rpca_metadata.csv")

exported_files <- c(
  "seurat_rpca_pca.csv",
  "seurat_rpca_umap.csv",
  "seurat_rpca_metadata.csv"
)
missing_outputs <- exported_files[!file.exists(exported_files)]
if (length(missing_outputs) > 0) {
  stop("Failed to write expected RPCA exports: ", paste(missing_outputs, collapse = ", "))
}

message("Exported files:")
for (f in exported_files) {
  message(" - ", f)
}
