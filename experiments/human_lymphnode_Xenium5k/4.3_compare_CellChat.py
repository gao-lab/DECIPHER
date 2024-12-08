# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: R_systerm
#     language: R
#     name: r_systerm
# ---

# %% [markdown]
# # Compare with CellChat
#
# We follow the CellChat spatial tutorial at https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.html

# %% vscode={"languageId": "r"}
source("renv/activate.R")  # this is very slow in cluster

# %% vscode={"languageId": "r"}
# renv::install("pak")
# pak::pkg_install("saeyslab/nichenetr")
# pak::pkg_install("jinworks/CellChat")
# pak::pkg_install("hdf5r")

# %% vscode={"languageId": "r"}
library(CellChat)
library(Seurat)
library(patchwork)
library(tidyverse)
library(future)


# %% vscode={"languageId": "r"}
path <- './data/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/cell_feature_matrix'
# load data
xenium.obj <- Read10X(path)
print(names(xenium.obj))
xenium.obj <- xenium.obj[[1]] # only get the valid matrix
# buid seurat object
seurat.obj <- CreateSeuratObject(counts = xenium.obj, project = "Xenium_Prime_Human_Lymph_Node_Reactive_FFPE")

# %% vscode={"languageId": "r"}
seurat.obj <- CreateSeuratObject(counts = xenium.obj, project = "Xenium_Prime_Human_Lymph_Node_Reactive_FFPE")
rownames(seurat.obj)

# %% vscode={"languageId": "r"}
saveRDS(seurat.obj, file = './results/raw_seurat_obj.rds')

# %% vscode={"languageId": "r"}
seurat.obj <- readRDS('./results/raw_seurat_obj.rds')

# %% vscode={"languageId": "r"}
meta <- read_csv('./data/obs.csv') %>% as.data.frame()
rownames(meta) <- meta$cell_id

seurat.obj@meta.data <- meta
seurat.obj@meta.data %>% head(2)

# %% vscode={"languageId": "r"}
seurat.obj <- NormalizeData(seurat.obj)

# %% [markdown]
# ## CellChat

# %% vscode={"languageId": "r"}
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# %% vscode={"languageId": "r"}
spatial.locs = meta %>% select( x_centroid, y_centroid) %>% as.data.frame() %>% as.matrix()

# %% vscode={"languageId": "r"}
# following https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html#seqfishmerfishstarmap
conversion.factor = 1
spot.size = 10 # use the typical human cell size
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input = Seurat::GetAssayData(seurat.obj, layer = "data")  # normalized data matrix


cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

# %% vscode={"languageId": "r"}
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 

# set the used database in the object
cellchat@DB <- CellChatDB.use

# %% vscode={"languageId": "r"}
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

future::plan("multisession", workers = 12)
options(future.globals.maxSize= 300 * 1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
future::plan('sequential')

# %% vscode={"languageId": "r"}
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.2,
                              contact.dependent = TRUE, contact.range = 100)

# %% vscode={"languageId": "r"}
cellchat <- filterCommunication(cellchat, min.cells = 10)


# %% vscode={"languageId": "r"}
cellchat <- computeCommunProbPathway(cellchat)


# %% vscode={"languageId": "r"}
cellchat <- aggregateNet(cellchat)


# %% vscode={"languageId": "r"}
?netVisual_bubble

# %% vscode={"languageId": "r"}
# chang the figure size
options(repr.plot.width=5, repr.plot.height=3)
# exclude 'Bcell' in sources
sources <- seurat.obj$cell_type %>% unique() %>% setdiff(c('Bcell', 'low_quality'))
netVisual_bubble(cellchat, sources.use = sources,targets.use = c('Bcell'))


# %% vscode={"languageId": "r"}
# get the version of all packages
sessionInfo()

# %% vscode={"languageId": "r"}
