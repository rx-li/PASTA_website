# Simulation data 
Our simulation based on a list of stere-seq datasets which can be downloaded at 
<https://db.cngb.org/stomics/artista/download/>. 

# Basic simulation code for PASTA 
```
import pandas as pd
import numpy as np
import scanpy as sc
import random
import anndata
import scipy
import statistics
import os, sys
sys.path.append('./pasta')
import __init__
import _version
import optimizer
import mapper
import utils 

ss_df = anndata.read_h5ad("path_to_the_h5ad_file")

pearson_cor = []
spearman_cor = []
mse = []
random.seed(233)

for i in range(30):
    #random sample 3000 cells and 1000 genes as spatial
    cells_n = random.sample(range(ss_df.shape[0]), 3000) ## modify if different region/celltype/subject
    genes_n = random.sample(range(ss_df.shape[1]), 1000)
    sp_adata = ss_df[cells_n, genes_n]
    true_adata = ss_df[cells_n, ]

    ss_df_2 = ss_df[~ss_df.obs.index.isin(sp_adata.obs.index), :]
    cells_n2 = random.sample(range(ss_df_2.shape[0]), 6000)
    sc_adata = ss_df_2[cells_n2, ]
    # remove all-zero-valued genes
    sc.pp.filter_genes(sp_adata, min_cells=1)
    sc.pp.filter_genes(sc_adata, min_cells=1)
    sc.pp.filter_genes(true_adata, min_cells=1)

    # get dist information
    dist = sp_adata.obsm["spatial"]
    dist = pd.DataFrame.from_records(dist)

    # random sample pathway
    idx_n = random.sample(range(sp_adata.shape[1]), 50)
    pathway = pd.DataFrame(sp_adata.var.index[idx_n])

    # train the model
    mapper.pp_adatas(sc_adata, sp_adata)
    ad_map = mapper.mapping(sc_adata, sp_adata, genes, sp_coords=coords, ncell_thres=10,
        sp_celltypes=cluster["Cluster"], lambda_g2=2, num_epochs=500)
    pthw_exp = utils.project_genes(adata_map=ad_map, adata_sc=sc_adata, pthw=genes)

    true_df = true_adata.to_df()
    true = true_df.loc[:, true_df.columns.isin(pathway["new_name"].str.lower())]
    true_v = true.sum(axis=1)
    pearson_cor.append(scipy.stats.pearsonr(pred_v, true_v)[0])
    spearman_cor.append(scipy.stats.spearmanr(pred_v, true_v)[0])
    mse.append(((pred_v - true_v)**2).mean())
```

# Basic simulation code for Tangram 
```
import pandas as pd
import numpy as np
import scanpy as sc
import random
import os, sys
sys.path.append('./')  # uncomment for local import
import tangram as tg
import anndata
import scipy
import statistics

ss_df = anndata.read_h5ad("path_to_the_h5ad_file")

pearson_cor = []
spearman_cor = []
mse = []
random.seed(233)

for i in range(30):
    #random sample 3000 cells and 1000 genes as spatial
    cells_n = random.sample(range(ss_df.shape[0]), 3000) ## modify if different region/celltype/subject
    genes_n = random.sample(range(ss_df.shape[1]), 1000)
    sp_adata = ss_df[cells_n, genes_n]
    true_adata = ss_df[cells_n, ]

    ss_df_2 = ss_df[~ss_df.obs.index.isin(sp_adata.obs.index), :]
    cells_n2 = random.sample(range(ss_df_2.shape[0]), 6000)
    sc_adata = ss_df_2[cells_n2, ]
    # remove all-zero-valued genes
    sc.pp.filter_genes(sp_adata, min_cells=1)
    sc.pp.filter_genes(sc_adata, min_cells=1)
    sc.pp.filter_genes(true_adata, min_cells=1)

    # get dist information
    dist = sp_adata.obsm["spatial"]
    dist = pd.DataFrame.from_records(dist)

    # random sample pathway
    idx_n = random.sample(range(sp_adata.shape[1]), 50)
    pathway = pd.DataFrame(sp_adata.var.index[idx_n])

    tg.pp_adatas(sc_adata, sp_adata)
    ad_map = tg.map_cells_to_space(sc_adata, sp_adata, num_epochs=1000)
    ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=sc_adata)
    ge_df = ad_ge.to_df()

    pred = ge_df.loc[:, ge_df.columns.isin(pathway["new_name"].str.lower())]
    pred_v = pred.sum(axis=1)
    true_df = true_adata.to_df()
    true = true_df.loc[:, true_df.columns.isin(pathway["new_name"])]
    true_v = true.sum(axis=1)

    pearson_cor.append(scipy.stats.pearsonr(pred_v, true_v)[0])
    spearman_cor.append(scipy.stats.spearmanr(pred_v, true_v)[0])
    mse.append(((pred_v - true_v)**2).mean())
```

# Basic simulation code for Seurat (run in R)
```
library(Seurat)

data <- readRDS("path_to_the_data.rds") # stereo-seq data has RDS file which can be read in R 

pearson_r <- NULL
spearman_r <- NULL
mse <- NULL
for(i in seq(30)){
        set.seed(i)

        # sample ST and scRNA-seq data
        cell_id <- sample(seq(dim(data)[2]), 3000)
        gene_id <- sample(seq(dim(data)[1]), 1000)
        cell_id2 <- sample(seq(dim(data)[2])[-cell_id], 3000)
        spatial_d <- data@assays[["SCT"]]@counts[gene_id, cell_id]
        sc_d <- data@assays[["SCT"]]@counts[, cell_id2]

        # run seurat
        sp_seurat <- CreateSeuratObject(counts = spatial_d, project = 'sp', assay = 'RNA')
        sp_seurat <- NormalizeData(sp_seurat)
        sp_seurat <- ScaleData(sp_seurat)

        train_seu <- NormalizeData(object = train_seu)
        train_seu <- FindVariableFeatures(object = train_seu, nfeatures = 2000)
        train_seu <- ScaleData(object = train_seu)
        train_seu <- RunPCA(object = train_seu, npcs = 50, verbose = FALSE)
        train_seu <- RunUMAP(object = train_seu, dims = 1:50, nneighbors = 5)

        features <- rownames(sp_seurat)

                anchors <- FindTransferAnchors(
                                       reference = train_seu,
                                       query = sp_seurat,
                                       features = features,
                                       dims = 1:30,
                                       reduction = 'cca')
        refdata <- GetAssayData(object = train_seu)
        imputation <- TransferData(anchorset = anchors,
                                     refdata = refdata,
                                     query = sp_seurat,
                                     weight.reduction = 'pca')
        pred_d <- imputation@assays[["id"]]@data

        pathway <- c(sample(rownames(st_d), 70)) # if all observed in st
        pred <- pred_d[rownames(pred_d) %in% pathway, ]
        pred_v <- colSums(pred)

        true_d <- data@assays[["SCT"]]@counts[, cell_id]
        true_f <- true_d[rownames(true_d) %in% pathway, ]
        true_v <- colSums(true)

        pearson_r <- c(pearson_r, cor(pred_v, true_v))
        spearman_r <- c(spearman_r, cor(pred_v, true_v, method="spearman"))
        mse <- c(mse, mean((pred_v - true_v)^2))
```