# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: R_systerm
#     language: R
#     name: r_systerm
# ---

# %% [markdown]
# # Run BASS
#
# We should run BASS manually due to it do not produce any embeddings but directly the final clustering results.

# %% vscode={"languageId": "r"}
library(stringr)
cwd <- getwd()
cwd <- str_replace(cwd, "workflow/manual", "")
# Set working directory
setwd(cwd)
print(cwd)
source(".Rprofile")

# %% vscode={"languageId": "r"}
library(Seurat)
# library(stringr)
# library(tidyverse)
library(ggplot2)
library(patchwork)
library(BASS)

# %% [markdown]
# ## 10X dataset

# %% tags=["parameters"] vscode={"languageId": "r"}
# parameter tag
input_file <- "./input/mimic/merged.h5ad"

# %% vscode={"languageId": "r"}
adata <- read.h5ad(input_file)

# %% vscode={"languageId": "r"}
cnts <- Matrix::t(adata$X)
cnts <- list(slice1 = cnts)


xys <- adata$obsm$spatial
colnames(xys) <- c("x", "y")
rownames(xys) <- rownames(adata$obs)
xys <- list(slice1 = xys)

# %% vscode={"languageId": "r"}
start <- proc.time()

# hyper-parameters
C <- 3# number of cell types
R <- 3 # number of spatial domains

# run BASS model
set.seed(0)
# Set up BASS object
BASS <- createBASSObject(cnts, xys, C = C, R = R, beta_method = "SW")

# Data pre-processing:
# 1.Library size normalization followed with a log2 transformation
# 2.Dimension reduction with PCA after standardizing all the genes
# 3.Batch effect adjustment using the Harmony package
BASS <- BASS.preprocess(BASS, geneSelect = "hvgs")

# Run BASS algorithm
BASS <- BASS.run(BASS)

# post-process posterior samples:
# 1.Adjust for label switching with the ECR-1 algorithm
# 2.Summarize the posterior samples to obtain the cell type labels, spatial
#   domain labels, and cell type proportion matrix estimate
BASS <- BASS.postprocess(BASS)

clabels <- BASS@results$c # cell type clusters
zlabels <- BASS@results$z # spatial domain labels
pi_est <- BASS@results$pi # cell type composition matrix


run_time <- as.numeric((proc.time() - start)[3])
print(paste0("Run time: ", run_time, " seconds"))

# %% [markdown]
# ## Xenium dataset

# %% vscode={"languageId": "r"}
# parameter tag
input_file <- "./input/mimic/merged.h5ad"
time_file <- ""
