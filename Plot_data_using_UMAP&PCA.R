rm(list =ls())

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(FactoMineR)
library(factoextra)
library(scrime)
library(Hmisc)
library(purrr)
library(uwot)
library(usefun)

# read data
snp_tbl = readr::read_csv(file = 'data/snp_tbl.csv')
X = snp_tbl %>% select(-c(1,2,3,4,5)) # keep just the SNP data
X_fct = X %>% purrr::modify_if(is.numeric, as.factor)
# fill in NA's with most frequent level (takes some minutes)
X_fct_complete = X_fct %>%
  mutate(across(.cols = everything(), ~impute(.x)))
# complete also the numerical SNP data
X_complete = X_fct_complete %>% purrr::modify_if(is.factor, as.numeric)

# save SNP tables
saveRDS(list(x = X, x_fct = X_fct, x_com = X_complete, x_fct_com = X_fct_complete), file = 'data/snp_data_list.rds')

# read saved data
if (FALSE) {
  snp_data = readRDS(file = 'data/snp_data_list.rds')
  X = snp_data$x
  X_fct = snp_data$x_fct
  X_complete = snp_data$x_com
  X_fct_complete = snp_data$x_fct_com
}

# check the substituted values
indxs = which(is.na(X_fct), arr.ind = TRUE)
imputed_values = list()
index = 1
for(i in 1:nrow(indxs)) {
  imputed_values[[index]] = X_fct_complete[indxs[i,'row'], indxs[i,'col']]
  index = index + 1
}
# some SNPs had `1` as the most frequent value and NA's were substituted with it (16340)
imputed_values %>% unlist(use.names = F) %>% table()

# convert to binary it might come in handy
X_bin = scrime::snp2bin(mat = X_complete, domrec = TRUE)
#X_bin2 = scrime::snp2bin(mat = X_complete, domrec = FALSE) # higher dimensionality


# PCA & MCA
set.seed(42) #for reproducibility
pca_res_unit_scale = FactoMineR::PCA(X_complete, scale.unit = TRUE, graph = FALSE) # fast
set.seed(42)
pca_res = FactoMineR::PCA(X_complete, scale.unit = FALSE, graph = FALSE) # fast
set.seed(42)
pca_res_bin = FactoMineR::PCA(X_bin, scale.unit = FALSE, graph = FALSE)

# PCA & MCA figures

## 12 distinct colors for population area
pop_col = c("#89100e", "#e5191d", "#ed7175", "#25496b", "#3a7db9", "#8db7db", "#235d22", "#54ad49",
                  "#9fd39b", "#414141", "#999999", "#cacaca")

geo_col = c("#e41a1c", "#377eb8", "#4daf4a", "#999999")
snp_tbl$geography <- factor(snp_tbl$geography,levels = c("RA", "VE", "MD", "NN" ))

## PCA unit scale
factoextra::fviz_pca_ind(pca_res_unit_scale, geom.ind = "point",
  col.ind = snp_tbl$geography %>% as.factor(), addEllipses = TRUE,
  legend.title = "Geography", pointshape = 20, palette = geo_col,
  title = "PCA (Unit scale)")
ggsave('fig/pca_unit_geo.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res_unit_scale, geom.ind = "point",
  col.ind = snp_tbl$sex %>% as.factor(), addEllipses = TRUE,
  legend.title = "Sex", pointshape = 20,
  title = "PCA (Unit scale)")
ggsave('fig/pca_unit_sex.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res_unit_scale, geom.ind = "point",
  col.ind = snp_tbl$poparea %>% as.factor(), addEllipses = TRUE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA (Unit scale)")
ggsave('fig/pca_unit_pop.pdf', width = 7, height = 5)

## PCA, no data scaling
factoextra::fviz_pca_ind(pca_res, geom.ind = "point",
  col.ind = snp_tbl$geography %>% as.factor(), addEllipses = TRUE,
  legend.title = "Geography", pointshape = 20,
  title = "PCA")
ggsave('fig/pca_geo.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res, geom.ind = "point",
  col.ind = snp_tbl$sex %>% as.factor(), addEllipses = TRUE,
  legend.title = "Sex", pointshape = 20,
  title = "PCA")
ggsave('fig/pca_sex.pdf', width = 7, height = 5)

 factoextra::fviz_pca_ind(pca_res, geom.ind = "point",
  col.ind = snp_tbl$poparea %>% as.factor(), addEllipses = TRUE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA")
ggsave('fig/pca_pop.pdf', width = 7, height = 5)

## PCA, binary data
factoextra::fviz_pca_ind(pca_res_bin, geom.ind = "point",
  col.ind = snp_tbl$geography %>% as.factor(), addEllipses = TRUE,
  legend.title = "Geography", pointshape = 20, palette = geo_col,
  title = "PCA (binary data)")
ggsave('fig/pca_bin_geo.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res_bin, geom.ind = "point",
  col.ind = snp_tbl$sex %>% as.factor(), addEllipses = TRUE,
  legend.title = "Sex", pointshape = 20,
  title = "PCA (binary data)")
ggsave('fig/pca_bin_sex.pdf', width = 7, height = 5)

 factoextra::fviz_pca_ind(pca_res_bin, geom.ind = "point",
  col.ind = snp_tbl$poparea %>% as.factor(), addEllipses = TRUE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA (binary data)")
ggsave('fig/pca_bin_pop.pdf', width = 7, height = 5)


# UMAP (numeric SNPs)
neighbors = c(13,15,20, 23, 27)
metrics = c('euclidean', 'cosine', 'manhattan', 'hamming')

for (nn in neighbors) {
  for (met in metrics) {
    print(paste0('Neighbors: ', nn, ", Metric: ", met))

    set.seed(42)
    umap_num_res = uwot::umap(X = X_complete, n_neighbors = nn,
      metric = met, verbose = TRUE)

    # UMAP figure (Geography coloring)
    umap_num_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = snp_tbl$geography %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = geo_col)+
      labs(title = paste0("UMAP (Numeric SNPs, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Geography",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/umap_num_', nn, '_', met,'_geo.pdf'), width = 7, height = 5)

    # UMAP figure (Sex coloring)
    umap_num_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = snp_tbl$sex %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      labs(title = paste0("UMAP (Numeric SNPs, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Sex",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/umap_num_', nn, '_', met,'_sex.pdf'), width = 7, height = 5)

    # UMAP figure (Population coloring)
    umap_num_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = snp_tbl$poparea %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = pop_col) +
      labs(title = paste0("UMAP (Numeric SNPs, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Pop",
        label.theme = element_text(size = 10),
        override.aes = list(size = 3)))
    ggsave(paste0('fig/umap_num_', nn, '_', met,'_pop.pdf'), width = 7, height = 5)
  }
}

# UMAP (binary SNPs)
for (nn in neighbors) {
  for (met in metrics) {
    print(paste0('Neighbors: ', nn, ", Metric: ", met))

    set.seed(42)
    umap_bin_res = uwot::umap(X = X_bin, n_neighbors = nn,
      metric = met, verbose = TRUE)

    # UMAP figure (Population coloring)
    umap_bin_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = snp_tbl$poparea %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = pop_col) +
      labs(title = paste0("UMAP (Binary SNPs, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Pop",
        label.theme = element_text(size = 10),
        override.aes = list(size = 3)))
    ggsave(paste0('fig/umap_bin_', nn, '_', met,'_pop.pdf'), width = 7, height = 5)

    # UMAP figure (Geography coloring)
    umap_bin_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = snp_tbl$geography %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = geo_col)+
      labs(title = paste0("UMAP (Binary SNPs, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Geography",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/umap_bin_', nn, '_', met,'_geo.pdf'), width = 7, height = 5)

    # UMAP figure (Sex coloring)
    umap_bin_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = snp_tbl$sex %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      labs(title = paste0("UMAP (Binary SNPs, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Sex",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/umap_bin_', nn, '_', met,'_sex.pdf'), width = 7, height = 5)
  }
}

writeLines(capture.output(xfun::session_info()), 'data/session_info.txt')
