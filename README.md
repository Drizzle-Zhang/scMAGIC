---
typora-root-url: figures
---

<img src="/Logo.png" style="zoom: 67%;" />

# scMAGIC

## Contents

[Overview](#Overview)

[Installation](#Installation)

[Tutorial](#Tutorial)

[Citation]()

[Contact]()

## Overview

scMAGIC (**S**ingle **C**ell classification based on **MA**ker **G**enes **I**dentification and expression profile **C**orrelation) is a novel single-cell classification tool using well-annotated bulk or single-cell RNA-seq data as the reference. The main idea of scMAGIC is to first identify the most likely cell type in the reference for a single cell based on expression profile correlation, and then verify its cell identity by using the corresponding marker genes of the reference cell type. Furthermore , scMAGIC conducts a second-round of reference-based classification by using target single cells with verified cell identity as the new reference to overcome the significant batch effects between the reference and the query data.

![workflow](D:\scRef\github\figures\workflow.png)

## Installation

#### Installing dependency package

Please install following R packages before using scMAGIC (Environment: R 4.0.0) :

```R
install.packages('parallel')      # 4.0.0
install.packages('pcaPP')         # 1.9-73
install.packages('Seurat')        # 3.2.0
install.packages('limma')         # 3.44.3
install.packages('RUVSeq')        # 1.22.0
install.packages('mclust')        # 5.4.6
install.packages('homologene')    # 1.4.68.19.3.27
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

<img src="D:\scRef\github\figures\TrueLabels.png" alt="TrueLabels" style="zoom: 33%;" />

```R
DimPlot(Obj.seurat, reduction = "umap",
        label = T, repel = T, group.by = 'pred.tag') +
    labs(title = 'Prediction labels') + theme(plot.title = element_text(hjust = 0.5))
```

<img src="D:\scRef\github\figures\PredLabels.png" alt="PredLabels" style="zoom:33%;" />

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

### Ⅲ. Inferring cell types of unassigned cells using Cell Atlas

```R
# extract unassigned cells
cell_id.unassigned <- row.names(output.scMAGIC[output.scMAGIC$scMAGIC.tag == 'Unassigned',])
exp.unassigned <- exp_sc_mat[, cell_id.unassigned]
label.unassigned <- list.target$label[cell_id.unassigned,]
# use scMAGIC to annotate the cells
data("MCA_ref")
output.unassigned <- scMAGIC(exp.unassigned, MCA_ref,
                             type_ref = 'sum-counts', use_RUVseq = F,
                             corr_use_HVGene1 = 2000, corr_use_HVGene2 = NULL,
                             num_threads = 8)
# classification results
table(label.unassigned, output.unassigned$scMAGIC.tag)
# combine results of reference-based annotation and atlas-based inference
output.new <- rbind(output.scMAGIC[output.scMAGIC$scMAGIC.tag != 'Unassigned',], output.unassigned)
pred.new <- pred.new[rownames(output.scMAGIC), 'scMAGIC.tag']
```

Results of Ⅱ and Ⅲ is showed in three UMAP plots as follows: the left is cell type labels from Campbell et al; the middle is scMAGIC assignments using Tasic dataset as reference; the right is combination of reference-based annotation and atlas-based inference.

![fig23](/fig23.png)

### Ⅳ. Reference-free annotation

In an exploratory study, users either do not have any knowledge about the cell types included in their study, or do not have a reference expression matrix available to use. In this case, an Atlas cell expression matrix can be used as the reference matrix, and scMAGIC can be applied to make an tentative classification of the target dataset. As a example, we use MCA as the reference to annotate the Campbell dataset.

```R
library(scMAGIC)
# load target dataset 
list.target <- readRDS('Campbell.Rdata')
exp_sc_mat <- list.target$mat_exp
label_sc <-list.target$label[,1]
# load MCA
data("MCA_ref")
# run scMAGIC
output.scMAGIC <- scMAGIC(exp_sc_mat, MCA_ref,
                          type_ref = 'sum-counts', use_RUVseq = F,
                          num_threads = 10)
pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
# classification results
table(true.tags, pred.scMAGIC)
```

Heatmap of reference-free annotation is showed as follows.

![heatmap_MCA_Campbell_scMAGIC](/heatmap_MCA_Campbell_scMAGIC.png)

## Contact

Please contact: 

Yu Zhang: zhang_yu18@fudan.edu.cn
