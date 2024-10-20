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

# %% vscode={"languageId": "r"}
library(tidyverse)
library(ggpubr)

# %% vscode={"languageId": "r"}
df <- read_csv('./results/run_time.csv')
df

# %% vscode={"languageId": "r"}
options(repr.plot.width=6, repr.plot.height=5)
df %>% pivot_longer(cols = 2:3, names_to = 'stage', values_to = 'time') %>%
    mutate(stage = factor(stage, levels = c('run_model', 'data_preprocess'))) %>%
    ggbarplot(x = 'GPUs', y = 'time', fill = 'stage', palette = "Paired", label = T) %>% 
    ggpar(legend.title = 'Stage', legend = 'right', ylab = 'Run time (minutes)', font.x = 15, font.y = 15)

# %% vscode={"languageId": "r"}
