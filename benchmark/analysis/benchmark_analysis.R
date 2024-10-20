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
data <- read_csv('../results/benchmark_results.csv')
data$method <- factor(data$method, levels = c('spider', 'scvi', 'harmony', 'scanpy', 'banksy', 'stagate', 'slat'))

# # replace ':' by '-' in all column names
# data <- data %>% rename_all(~ gsub(":", "__", .))

head(data)

# %% vscode={"languageId": "r"}
p_list <- list()
for(var in colnames(data)){
    if( grepl('batch__', var) ){
        p <- ggbarplot(data, x="method", y=var, fill='method', color='method', add="mean_se", error.plot="pointrange",
                    palette='aaas', xlab=F ,ylab=var, size=1.5, lab.size=5, label=F, label.pos="out",
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))
        p_list[[var]] <- p
    }
}

options(repr.plot.width=20, repr.plot.height=8)
(p_list[[1]] | p_list[[3]] | p_list[[5]] | p_list[[7]]) / (p_list[[2]] | p_list[[4]] | p_list[[6]] | p_list[[8]])

# %% vscode={"languageId": "r"}
p_list <- list()
for(var in colnames(data)){
    if( grepl('gex__', var) & !grepl('nbr', var) & !grepl('asw', var)){
        print(var)
        p <- ggbarplot(data, x="method", y=var, fill='method', color='method', add="mean_se", error.plot="pointrange",
                    palette='aaas', xlab=F ,ylab=var, size=1.5, lab.size=5, label=F, label.pos="out",
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))
        p_list[[var]] <- p
    }
}

options(repr.plot.width=11, repr.plot.height=4)
p_list[[1]] | p_list[[2]]

# %% vscode={"languageId": "r"}
p_list <- list()
for(var in colnames(data)){
    if( grepl('nbr__', var) & !grepl('gex', var) & !grepl('asw', var)){
        print(var)
        p <- ggbarplot(data, x="method", y=var, fill='method', color='method', add="mean_se", error.plot="pointrange",
                    palette='aaas', xlab=F ,ylab=var, size=1.5, lab.size=5, label=F, label.pos="out",
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))
        p_list[[var]] <- p
    }
}

options(repr.plot.width=11, repr.plot.height=4)
p_list[[1]] | p_list[[2]]

# %% vscode={"languageId": "r"}
options(repr.plot.width=5, repr.plot.height=4)
ggbarplot(data, x="method", y='run_time', fill='method', color='method', add="mean_se", error.plot="pointrange",
                    palette='aaas', xlab=F ,ylab=var, size=1.5, lab.size=5, label=F, label.pos="out",
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            stat_summary(fun.data = function(x) data.frame(y = mean(x) + 20, label = paste("", round(mean(x), 1))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))

# %% [markdown]
# ## MERFISH brain

# %% vscode={"languageId": "r"}
data <- read_csv('../results/benchmark_results.csv')
data$method <- factor(data$method, levels = c('spider', 'scvi', 'harmony', 'scanpy', 'banksy', 'stagate', 'slat'))
data <- data %>% filter(dataset == 'human_brainaging_merfish') %>% filter(method != 'harmony')

head(data)

# %% vscode={"languageId": "r"}
p_list <- list()
for(var in colnames(data)){
    if( grepl('gex__', var) & !grepl('nbr', var) & !grepl('asw', var)){
        print(var)
        p <- ggbarplot(data, x="method", y=var, fill='method', color='method', add="mean_se", error.plot="pointrange",
                    palette='aaas', xlab=F ,ylab=var, size=1.5, lab.size=5, label=F, label.pos="out",
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))
        p_list[[var]] <- p
    }
}

options(repr.plot.width=10, repr.plot.height=4)
p_list[[1]] | p_list[[2]]

# %% vscode={"languageId": "r"}
p_list <- list()
for(var in colnames(data)){
    if( grepl('nbr__', var) & !grepl('gex', var) & !grepl('asw', var)){
        print(var)
        p <- ggbarplot(data, x="method", y=var, fill='method', color='method', add="mean_se", error.plot="pointrange",
                    palette='aaas', xlab=F ,ylab=var, size=1.5, lab.size=5, label=F, label.pos="out",
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))
        p_list[[var]] <- p
    }
}

options(repr.plot.width=10, repr.plot.height=4)
p_list[[1]] | p_list[[2]]

# %% [markdown]
# ## Xenium breast tumor

# %% vscode={"languageId": "r"}
data <- read_csv('../results/benchmark_results.csv')
data$method <- factor(data$method, levels = c('spider', 'scvi', 'harmony', 'scanpy', 'banksy', 'stagate', 'slat'))
data <- data %>% filter(dataset == 'human_breastcancer_xenium') %>% filter(method != 'harmony')


tail(data)

# %% vscode={"languageId": "r"}
p_list <- list()
for(var in colnames(data)){
    if( grepl('gex__', var) & !grepl('nbr', var) & !grepl('asw', var)){
        print(var)
        p <- ggbarplot(data, x="method", y=var, fill='method', color='method', add="mean_se", error.plot="pointrange",
                    palette='aaas', xlab=F ,ylab=var, size=1.5, lab.size=5, label=F, label.pos="out",
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))
        p_list[[var]] <- p
    }
}

options(repr.plot.width=10, repr.plot.height=4)
p_list[[1]] | p_list[[2]]

# %% vscode={"languageId": "r"}
p_list <- list()
for(var in colnames(data)){
    if( grepl('nbr__', var) & !grepl('gex', var) & !grepl('asw', var)){
        print(var)
        p <- ggbarplot(data, x="method", y=var, fill='method', color='method', add="mean_se", error.plot="pointrange",
                    palette='aaas', xlab=F ,ylab=var, size=1.5, lab.size=5, label=F, label.pos="out",
                    legend.title='Method') %>%
            ggpar(legend='right', font.legend=16, font.subtitle=16, font.xtickslab=18, font.ytickslab=16, font.y=18) +
            stat_summary(fun.data = function(x) data.frame(y = mean(x) + 0.05, label = paste("", round(mean(x), 2))), geom="text", size=5.5) +
            theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5), strip.text.x=element_text(size=18))
        p_list[[var]] <- p
    }
}

options(repr.plot.width=10, repr.plot.height=4)
p_list[[1]] | p_list[[2]]
