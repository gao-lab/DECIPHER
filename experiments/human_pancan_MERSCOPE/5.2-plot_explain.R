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
library(ggpubr)
library(tidyverse)

# %% vscode={"languageId": "r"}
df <- read_csv('./results/explain_results.csv')
df <- df %>%
  filter(!Index %in% c(0,1)) %>%
  filter(!cell_type %in% c('Hepatocyte' , 'pneumocyte', 'unknown')) %>%
  mutate(Index = Index - 1) %>%
  mutate(Sample = paste0(Tissue_1, '-', Index)) %>%
  group_by(Index) %>%
  arrange(test_r2_score, .by_group = TRUE) %>%
  mutate(rank = row_number())
head(df)

# %% vscode={"languageId": "r"}
options(repr.plot.width=12, repr.plot.height=9)
ggbarplot(df, x = 'cell_type', y = 'test_r2_score', fill = 'Tissue_1', color = 'Tissue_1', palette = 'jco') %>%
    facet(facet.by = 'Sample', scales = 'free_y') %>%
    ggpar(legend = 'right', legend.title = 'Cancer', xlab = '', ylab = "Variance explained (R2)", x.text.angle = 45) +
    geom_hline(yintercept = 0, linetype = 2)

# %% vscode={"languageId": "r"}
ggbarplot(df, x = 'cell_type', y = 'rank', fill = 'Tissue_1', color = 'Tissue_1', palette = 'Set1') %>% facet(facet.by = 'Sample')

# %% vscode={"languageId": "r"}
options(repr.plot.width=6, repr.plot.height=5)
ggboxplot(df, x = 'cell_type', y = 'test_r2_score', fill = 'cell_type', palette = 'jco') %>%
    ggpar(xlab = '', legend.title = 'Cell type', legend = 'right', ylab = "Variance explained (R2)", x.text.angle = 45 ) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_compare_means(method = "anova", label.x = 3.5, label.y =  0.09, )


# %% vscode={"languageId": "r"}
ggboxplot(df, x = 'cell_type', y = 'rank', fill = 'cell_type', palette = 'jco') %>%
    ggpar(xlab = '', legend.title = 'Cell type', legend = 'right', ylab = "Variance explained (R2)",  x.text.angle = 45) +
    stat_compare_means(method = "anova", label.x = 3.5, label.y =  6.5, )

# %% vscode={"languageId": "r"}
