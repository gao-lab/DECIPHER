# ---
# jupyter:
#   jupytext:
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

# %% vscode={"languageId": "r"}
library(ggplot2)
source('./utils.R')

# %% [markdown]
# Plot the clustering benchmark results

# %% vscode={"languageId": "r"}
plotSingleTaskCluster(
    csv_metrics_path = "../results/scIB_cluster_vis.csv",
    outdir = "./results/cluster",
    weight_batch = 0.5
)

# %% [markdown]
# Plot the batch effect correction benchmark results

# %% vscode={"languageId": "r"}
plotSingleTaskBatch(
    csv_metrics_path = "../results/scIB_batch_vis.csv",
    outdir = "./results/batch",
    weight_batch = 0.5
)
