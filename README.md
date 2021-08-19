<img src="https://github.com/Drizzle-Zhang/scMAGIC/raw/main/figures/Logo.png" width="200">

# scMAGIC

## Contents

[Overview](#Overview)

[Installation](#Installation)

[Tutorial](#Tutorial)

[Citation]()

[Contact]()

## Overview

scMAGIC (**S**ingle **C**ell annotation using **MA**rker **G**enes **I**dentification and two rounds of reference-based **C**lassification(RBC)) is a novel single-cell classification tool using well-annotated bulk or single-cell RNA-seq data as the reference. The main idea of scMAGIC is to first identify the most likely cell type in the reference for a single cell based on expression profile correlation, and then verify its cell identity by using the corresponding marker genes of the reference cell type. Furthermore , scMAGIC conducts a second-round of RBC by using target single cells with verified cell identity as the new reference to overcome the significant batch effects between the reference and the query data. Therefore, scMAGIC's advantage is especially significant when the cell types in the query dataset are not completely covered by the reference dataset and when there exist significant batch effects between the reference and the query datasets. Moreover, when no reference dataset is available, scMAGIC can annotate query cells with reasonably high accuracy by using an atlas dataset as the reference.<img src="https://github.com/Drizzle-Zhang/scMAGIC/raw/main/figures/workflow.png" width="900">

## Installation

#### Installing dependency package

Please install following R packages before using scMAGIC (Environment: R 4.0.0) :

```R
install.packages('parallel')        # 4.0.0
install.packages('pcaPP')           # 1.9-73
BiocManager::install('Seurat')      # 3.2.0
BiocManager::install('limma')       # 3.44.3
BiocManager::install('RUVSeq')      # 1.22.0
install.packages('mclust')          # 5.4.6
install.packages('homologene')      # 1.4.68.19.3.27
```

#### Installing scMAGIC

You can install scMAGIC from Github:

```
if (!requireNamespace("devtools", quietly = TRUE)) {
	install.packages("devtools")
}
devtools::install_github("Drizzle-Zhang/scMAGIC")
```

## Tutorial

### Ⅰ. Reference-based annotation

#### Download datasets

The reference dataset includes 1727 mouse primary visual cortex cells profiled by SMART-seq. The target dataset is from a study profiling mouse hypothalamic arcuate-median eminence complex cells by Drop-seq, here we sample randomly 2000 cells to make a simple test.

```shell
cd your_path
# download reference dataset
wget https://github.com/Drizzle-Zhang/scMAGIC_scripts/raw/main/data/Tasic.Rdata
# download target dataset
wget https://github.com/Drizzle-Zhang/scMAGIC_scripts/raw/main/data/Campbell_2k.Rdata
```

#### Library and load data

```R
library(scMAGIC)

setwd('your_path')
# load reference dataset
list.Ref <- readRDS('Tasic.Rdata')
ref.mtx <- list.Ref$mat_exp
ref.labels <-list.Ref$label[, 1]
# load target dataset
list.demo <- readRDS('Campbell_2k.Rdata')
exp_sc_mat <- list.demo$mat_exp
label_sc <-list.demo$label
```

#### Cell type classification by scMAGIC

```R
output.scMAGIC <- scMAGIC(exp_sc_mat, ref.mtx, ref.labels, num_threads = 4)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
```

#### Visualization

```R
library(Seurat)
Obj.seurat <- CreateSeuratObject(counts = exp_sc_mat)
Obj.seurat <- NormalizeData(Obj.seurat)
Obj.seurat <- FindVariableFeatures(Obj.seurat, nfeatures = 2000)
Obj.seurat <- ScaleData(Obj.seurat)
Obj.seurat <- RunPCA(Obj.seurat)
Obj.seurat <- RunUMAP(Obj.seurat, dims = 1:50)
Obj.seurat@meta.data$original.label <- label_sc
Obj.seurat@meta.data$pred.tag <- pred.scMAGIC

library(ggplot2)
DimPlot(Obj.seurat, reduction = "umap",
        label = T, repel = T, group.by = 'original.label') +
    labs(title = 'True labels') + theme(plot.title = element_text(hjust = 0.5))
```

<img src="https://github.com/Drizzle-Zhang/scMAGIC/raw/main/figures/TrueLabels.png" width="500">

```R
DimPlot(Obj.seurat, reduction = "umap",
        label = T, repel = T, group.by = 'pred.tag') +
    labs(title = 'Prediction labels') + theme(plot.title = element_text(hjust = 0.5))
```

<img src="https://github.com/Drizzle-Zhang/scMAGIC/raw/main/figures/PredLabels.png" width="500">

### Ⅱ. Annotation of large dataset

If the target dataset contains more than 5000 cells, scMAGIC will accelerate the computation by merge similar cells automatically. Here, we identify cell types of more than 20000 cells in Campbell dataset.

#### Download dataset

```shell
wget https://github.com/Drizzle-Zhang/scMAGIC_scripts/raw/main/data/Campbell.Rdata
```

#### Run scMAGIC

```R
library(scMAGIC)
# load target dataset 
list.target <- readRDS('Campbell.Rdata')
exp_sc_mat <- list.target$mat_exp
label_sc <-list.target$label
# print runtime
time1 <- Sys.time()
output.scMAGIC <- scMAGIC(exp_sc_mat, ref.mtx, ref.labels, num_threads = 10)
time2 <- Sys.time()
difftime(time2, time1, units = 'mins')
# Time difference of 3.281353 mins

# classification results
table(label_sc, pred.scMAGIC)
```

Results of annotation are showed in three UMAP plots as follows: the left is cell type labels from Campbell et al; the right is scMAGIC assignments using Tasic dataset as reference.

<img src="https://github.com/Drizzle-Zhang/scMAGIC/blob/main/figures/large_dataset.png" width="900">

### Ⅲ. Reference-free annotation

In an exploratory study, users either do not have any knowledge about the cell types included in their study, or do not have a reference expression matrix available to use. In this case, an Atlas cell expression matrix can be used as the reference matrix, and scMAGIC can be applied to make an tentative classification of the target dataset. As a example, we use MCA as the reference to annotate the mouse neocortex dataset.

#### Download dataset

```shell
wget https://github.com/Drizzle-Zhang/scMAGIC_scripts/raw/main/data/MouseNeocortex.Rdata
```

#### Run scMAGIC

```R
library(scMAGIC)
# load target dataset 
list.target <- readRDS('MouseNeocortex.Rdata')
exp_sc_mat <- list.target$mat_exp
label_sc <-list.target$label
# load MCA
data("MCA_ref")
# run scMAGIC
output.scMAGIC <- scMAGIC(exp_sc_mat, MCA_ref,
                          type_ref = 'sum-counts', use_RUVseq = F,
                          min_cell = 5, num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
# classification results
table(label_sc, pred.scMAGIC)
```

Heatmap of reference-free annotation is showed as follows.

<img src="https://github.com/Drizzle-Zhang/scMAGIC/blob/main/figures/heatmap_MCA_MouseNeocortex.png" width="500">

## Contact

Please contact: 

Yu Zhang: zhang_yu18@fudan.edu.cn

