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
# # Compare with NicheNet

# %% vscode={"languageId": "r"}
source("renv/activate.R")  # this is very slow in cluster

# %% vscode={"languageId": "r"}
# renv::install("pak")
# pak::pkg_install("saeyslab/nichenetr")
# pak::pkg_install("jinworks/CellChat")
# pak::pkg_install("hdf5r")
pak::pkg_install("immunogenomics/presto")

# %% vscode={"languageId": "r"}
library(nichenetr)
library(Seurat)
library(patchwork)
library(tidyverse)


# %% vscode={"languageId": "r"}
seuratObj <- readRDS('./results/raw_seurat_obj.rds')

# %% vscode={"languageId": "r"}
meta <- read_csv('./data/obs.csv') %>% as.data.frame()
rownames(meta) <- meta$cell_id

seuratObj@meta.data <- meta
seuratObj@meta.data %>% head(2)

# %% vscode={"languageId": "r"}
seuratObj <- NormalizeData(seuratObj)

# %% vscode={"languageId": "r"}
lr_network <- readRDS('./data/nichenet/lr_network_human_21122021.rds')
ligand_target_matrix <- readRDS('./data/nichenet/ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS('./data/nichenet/weighted_networks_nsga2r_final.rds')

# %% [markdown]
# ## Run NicheNet

# %% vscode={"languageId": "r"}
lr_network <- lr_network %>% distinct(from, to)


# %% vscode={"languageId": "r"}
seuratObj@meta.data$cell_type %>% table()

# %% vscode={"languageId": "r"}
Idents(seuratObj) <- 'cell_type'

receiver = "Bcell"
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)


# %% vscode={"languageId": "r"}
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

# %% vscode={"languageId": "r"}
sender_celltypes <- seuratObj@meta.data$cell_type %>% unique()

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 8492
length(potential_ligands)
## [1] 483
length(potential_ligands_focused)
## [1] 127

# %% vscode={"languageId": "r"}
# here we use different spatial niche as different conditions
condition_oi <-  "5"
condition_reference <- "8"

seurat_obj_receiver <- subset(seuratObj, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "leiden",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

# %% vscode={"languageId": "r"}
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


# %% vscode={"languageId": "r"}
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

# %% vscode={"languageId": "r"}
p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(15, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

# %% vscode={"languageId": "r"}
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

# %% vscode={"languageId": "r"}
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 


vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))

# %% vscode={"languageId": "r"}
