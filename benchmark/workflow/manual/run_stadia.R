# NOTE: STADIA only works for batch correction

library(stringr)
cwd <- getwd()
if (grepl("workflow/manual", cwd)) {
  cwd <- str_replace(cwd, "workflow/manual", "")
}
# Set working directory
setwd(cwd)
print(cwd)
source(".Rprofile")

## % Load libraries
library(Seurat)
# library(stringr)
library(ggplot2)
library(patchwork)
library(stadia)
options(Seurat.object.assay.version = "v3")  # must use v3 mode in Seurat v5 for compatibility

input_file_1 <- "./data/human_mimic_10xPBMC/data/mimic_results/pbmc_batch/batch1.h5ad"
input_file_2 <- "./data/human_mimic_10xPBMC/data/mimic_results/pbmc_batch/batch2.h5ad"

# load data
adata_1 <- read.h5ad(input_file_1)
adata_2 <- read.h5ad(input_file_2)
for (i in 1:2) {
  adata <- get(paste0("adata_", i))
  mat <- Matrix::t(adata$X)
  slice <- CreateSeuratObject(counts = mat, meta.data = adata$obs, assay = "Spatial")
  # get spatial info
  colnames(adata$obsm$spatial) <- c('x', 'y')
  slice$row <- adata$obsm$spatial[,'x']
  slice$col <- adata$obsm$spatial[,'y']
  assign(paste0("slice_", i), slice)
}

slices <- list(slice1 = slice_1, slice2 = slice_2)  # a trick


# run model
K <- 7
d <- 35
etas <- 0.15
set.seed(123)

system.time({
    ## set hyperparameters
    hyper <- HyperParameters(slices, d = d, eta = etas)
    ## run model
    stadia_res <- stadia(
        slices,
        hyper,
        dim = d,
        n_cluster = K,
        platform = "others",
        em.maxiter = 10,
        adj.cutoff = 0.03,
        min.features = 0,
        min.spots = 0,
    )
})

saveRDS(stadia_res, './workflow/manual/results/mimic_STADIA_res.rds')
stadia_res$c_vec
