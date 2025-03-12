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
# this is R code
library(ggpubr)
library(tidyverse)
library(patchwork)

# %% [markdown]
# ## 10X dataset

# %% vscode={"languageId": "r"}
data <- read_csv('../results/benchmark_results_final.csv')

data$method <- gsub('decipher', 'DECIPHER', data$method)
data$method <- gsub('scvi', 'scVI', data$method)
data$method <- gsub('scanpy', 'Scanpy', data$method)
data$method <- gsub('banksy', 'Banksy', data$method)
data$method <- gsub('stagate', 'STAGATE', data$method)
data$method <- gsub('slat', 'SLAT', data$method)
data$method <- gsub('harmony', 'Harmony', data$method)
data$method <- gsub('scniche_raw', 'scNiche', data$method)
data$method <- gsub('scniche', 'scNiche-minibatch', data$method)
data <- data %>% filter(dataset == 'human_pbmc_10x_mimic') %>% filter(method != 'scNiche')

data$method <- factor(data$method, levels = c('DECIPHER', 'scVI', 'Scanpy', 'Banksy', 'STAGATE', 'SLAT', 'scNiche-minibatch', 'Harmony'))

# # replace ':' by '-' in all column names
# data <- data %>% rename_all(~ gsub(":", "__", .))

head(data)

# %% vscode={"languageId": "r"}
p_list <- list()
for(var in colnames(data)){
    if( grepl('batch_', var) ){
        var_name <- gsub('_', ' ', var)
        var_name <- gsub('gex', 'Omics', var_name)
        var_name <- gsub('nbr', 'Spatial', var_name)
        var_name <- gsub('asw', 'ASW', var_name)
        var_name <- gsub('gc', 'GC', var_name)
        var_name <- gsub('ilisi', 'iLISI', var_name)
        var_name <- gsub('kbet', 'kBET', var_name)
        p <- ggbarplot(data, x="method", y=var, fill='method', color='method', add=c("mean_se","point"), alpha=0.5, ylim=c(0, 1),
                    error.plot="pointrange", palette='aaas', xlab=F ,ylab=var_name, size=1.5, label=F,
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            # stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=50, vjust=1, hjust=1), strip.text.x=element_text(size=18)) + rremove('legend')
        p_list[[var]] <- p
    }
}

options(repr.plot.width=15, repr.plot.height=10)
(p_list[[5]] | p_list[[6]] | p_list[[7]] | p_list[[8]]) / (p_list[[1]] | p_list[[2]] | p_list[[3]] | p_list[[4]])

# %% vscode={"languageId": "r"}
data <- data %>% filter(method != 'Harmony')


# %% vscode={"languageId": "r"}
# options(repr.plot.width=5, repr.plot.height=5)
# ggscatter(data = data, x = 'gex_nmi', y = 'nbr_nmi', color = 'method',
#             xlim = c(0, 0.8), ylim = c(0, 0.85), palette='aaas')
# ggscatter(data = data, x = 'gex_ari', y = 'nbr_ari', color = 'method',
#             xlim = c(0, 0.8), ylim = c(0, 0.85), palette='aaas')

# %% vscode={"languageId": "r"}
options(repr.plot.width=8, repr.plot.height=5)

pivot_longer(data, cols = c(0:14), names_to = "metric", values_to = 'score') %>%
    filter(metric %in% c('nbr_nmi', 'nbr_ari')) %>%
    ggbarplot("metric", 'score', fill='method', color='method', add=c("mean_se","point"), alpha=0.5, ylim=c(0, 1),
        position = position_dodge(),palette='aaas', xlab=F ,ylab='Score', size=1.5, label=F, legend.title='Method'
        ) %>%
        ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
        # stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
        theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))

pivot_longer(data, cols = c(0:14), names_to = "metric", values_to = 'score') %>%
    filter(metric %in% c('gex_nmi', 'gex_ari')) %>%
    ggbarplot("metric", 'score', fill='method', color='method', add=c("mean_se","point"), alpha=0.5, ylim=c(0, 1),
        position = position_dodge(),palette='aaas', xlab=F ,ylab='Score', size=1.5, label=F, legend.title='Method'
        ) %>%
        ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
        # stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
        theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))

# %% [markdown]
# ## MERFISH brain

# %% vscode={"languageId": "r"}
data <- read_csv('../results/benchmark_results_final.csv')
# add NA row for stagate
data[nrow(data) + 1,] <- list(0,0,0,0,0,0,0,0,0,0,0,0,0,0, 'human_brainaging_merfish', 0, 'stagate') # stagate can not run on merfish
data[nrow(data) + 1,] <- list(0,0,0,0,0,0,0,0,0,0,0,0,0,0, 'human_brainaging_merfish', 0, 'scniche_raw') # scniche can not run on merfish


data$method <- gsub('decipher', 'DECIPHER', data$method)
data$method <- gsub('scvi', 'scVI', data$method)
data$method <- gsub('scanpy', 'Scanpy', data$method)
data$method <- gsub('banksy', 'Banksy', data$method)
data$method <- gsub('stagate', 'STAGATE', data$method)
data$method <- gsub('slat', 'SLAT', data$method)
data$method <- gsub('harmony', 'Harmony', data$method)
data$method <- gsub('scniche_raw', 'scNiche', data$method)
data$method <- gsub('scniche', 'scNiche-minibatch', data$method)
data <- data %>% filter(dataset == 'human_brainaging_merfish') %>% filter(method != 'Harmony') %>% filter(method != 'scNiche')


data$method <- factor(data$method, levels = c('DECIPHER', 'scVI', 'Scanpy', 'Banksy', 'STAGATE', 'SLAT', 'scNiche-minibatch', 'Harmony'))

head(data)

# %% vscode={"languageId": "r"}
options(repr.plot.width=8, repr.plot.height=5)


pivot_longer(data, cols = c(0:14), names_to = "metric", values_to = 'score') %>%
    filter(metric %in% c('nbr_nmi', 'nbr_ari')) %>%
    ggbarplot("metric", 'score', fill='method', color='method', add=c("mean_se","point"), alpha=0.5, ylim=c(0, 1),
        position = position_dodge(),palette='aaas', xlab=F ,ylab='Score', size=1.5, label=F, legend.title='Method'
        ) %>%
        ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
        # stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
        theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))


pivot_longer(data, cols = c(0:14), names_to = "metric", values_to = 'score') %>%
    filter(metric %in% c('gex_nmi', 'gex_ari')) %>%
    ggbarplot("metric", 'score', fill='method', color='method', add=c("mean_se","point"), alpha=0.5, ylim=c(0, 1),
        position = position_dodge(),palette='aaas', xlab=F ,ylab='Score', size=1.5, label=F, legend.title='Method'
        ) %>%
        ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
        # stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
        theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))

# %% [markdown]
# ## Xenium breast tumor

# %% vscode={"languageId": "r"}
data <- read_csv('../results/benchmark_results_final.csv')
data[nrow(data) + 1,] <- list(0,0,0,0,0,0,0,0,0,0,0,0,0,0, 'human_breastcancer_xenium', 0, 'scniche_raw') # scniche can not run on merfish


data$method <- gsub('decipher', 'DECIPHER', data$method)
data$method <- gsub('scvi', 'scVI', data$method)
data$method <- gsub('scanpy', 'Scanpy', data$method)
data$method <- gsub('banksy', 'Banksy', data$method)
data$method <- gsub('stagate', 'STAGATE', data$method)
data$method <- gsub('slat', 'SLAT', data$method)
data$method <- gsub('harmony', 'Harmony', data$method)
data$method <- gsub('scniche_raw', 'scNiche', data$method)
data$method <- gsub('scniche', 'scNiche-minibatch', data$method)
data <- data %>% filter(dataset == 'human_breastcancer_xenium') %>% filter(method != 'Harmony') %>% filter(method != 'scNiche')

data$method <- factor(data$method, levels = c('DECIPHER', 'scVI', 'Scanpy', 'Banksy', 'STAGATE', 'SLAT', 'scNiche-minibatch', 'Harmony'))



tail(data)

# %% vscode={"languageId": "r"}
options(repr.plot.width=8, repr.plot.height=5)

pivot_longer(data, cols = c(0:14), names_to = "metric", values_to = 'score') %>%
    filter(metric %in% c('nbr_nmi', 'nbr_ari')) %>%
    ggbarplot("metric", 'score', fill='method', color='method', add=c("mean_se","point"), alpha=0.5, ylim=c(0, 1),
        position = position_dodge(),palette='aaas', xlab=F ,ylab='Score', size=1.5, label=F, legend.title='Method'
        ) %>%
        ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
        # stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
        theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))


pivot_longer(data, cols = c(0:14), names_to = "metric", values_to = 'score') %>%
    filter(metric %in% c('gex_nmi', 'gex_ari')) %>%
    ggbarplot("metric", 'score', fill='method', color='method', add=c("mean_se","point"), alpha=0.5, ylim=c(0, 1),
        position = position_dodge(),palette='aaas', xlab=F ,ylab='Score', size=1.5, label=F, legend.title='Method'
        ) %>%
        ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
        # stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
        theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))
