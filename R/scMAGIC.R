########################################
#Author: Yu Zhang
#Email: zhang_yu18@fudan.edu.cn
#######################################

get_overlap_genes <- function(exp_sc_mat, exp_ref_mat) {
    # exp_ref_mat <- as.data.frame(as.matrix(exp_ref_mat))
    # exp_sc_mat <- as.data.frame(as.matrix(exp_sc_mat))
    # get overlap genes
    exp_sc_mat <- exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat <- exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc <- rownames(exp_sc_mat)
    gene_ref <- rownames(exp_ref_mat)
    gene_over <- gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat <- exp_sc_mat[gene_over,]
    exp_ref_mat <- exp_ref_mat[gene_over,]

    out.overlap <- list()
    out.overlap$exp_sc_mat <- exp_sc_mat
    out.overlap$exp_ref_mat <- exp_ref_mat
    out.overlap$gene_over <- gene_over
    return(out.overlap)

}


.get_high_variance_genes <- function(exp_ref_mat, num.genes = 2000, type_ref = 'sum-counts') {
    if (type_ref %in% c('sum-counts', 'sc-counts')) {
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                    scale.factor = 1e6, verbose = F)
    }
    if (type_ref %in% c('fpkm', 'tpm', 'rpkm')) {
        exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref@assays$RNA@data <- exp_ref_mat
    }
    seurat.Ref <- FindVariableFeatures(
        seurat.Ref,
        selection.method = "vst",
        nfeatures = num.genes,
        verbose = F
    )
    return(VariableFeatures(seurat.Ref))

}


.one_multinomial <- function(i, exp_sc_mat, exp_ref_mat, colname_ref,
                             verbose, print_step) {
    delta <- 0.5
    Refprob <- function(exp_sc, exp_ref) {
        log_p_sc_given_ref <- dmultinom(x = exp_sc, log = T, prob = exp_ref)
        return(log_p_sc_given_ref)
    }
    #################
    exp_sc <- as.array(exp_sc_mat[, i])
    log_p_sc_given_ref_list <- numeric(length = length(colname_ref))
    j = 1
    while (j <= length(colname_ref)) {
        exp_ref <- as.array(exp_ref_mat[, j])
        #####
        exp_ref[which(exp_ref == 0)] <- delta * min(exp_ref[which(exp_ref > 0)])
        #####
        log_p_sc_given_ref <- Refprob(exp_sc, exp_ref)
        log_p_sc_given_ref_list[j] <- log_p_sc_given_ref
        j = j + 1
    }
    ################################
    if (verbose) {
        if (i %% print_step == 1) {
            print(i)
        }
    }

    return(log_p_sc_given_ref_list)

}


.get_log_p_sc_given_ref <- function(exp_sc_mat, exp_ref_mat, num_threads = 4,
                                    print_step = 100, gene_overlap = FALSE,
                                    verbose = FALSE) {
    library(parallel, verbose = F)

    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name

    #Step 1. get overlapped genes
    if (gene_overlap) {
        if (verbose) {
            print('Gene number of exp_sc_mat:')
            print(nrow(exp_sc_mat))
            print('Gene number of exp_ref_mat:')
            print(nrow(exp_ref_mat))
        }
        out.overlap <- get_overlap_genes(exp_sc_mat, exp_ref_mat)
        exp_sc_mat <- out.overlap$exp_sc_mat
        exp_ref_mat <- out.overlap$exp_ref_mat
        if (verbose) {
            print('Number of overlapped genes:')
            print(nrow(exp_sc_mat))
        }
    }

    ###############
    colname_sc <- colnames(exp_sc_mat)
    colname_ref <- colnames(exp_ref_mat)

    #Step 2. calculate prob
    cl = makeCluster(num_threads, outfile = '')
    RUN <- parLapply(
        cl = cl,
        1:length(exp_sc_mat[1,]),
        .one_multinomial,
        exp_sc_mat = exp_sc_mat,
        exp_ref_mat = exp_ref_mat,
        colname_ref = colname_ref,
        verbose = verbose,
        print_step = print_step
    )
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=num_threads)
    LOG_P_SC_GIVEN_REF = c()
    for (log_p_sc_given_ref_list in RUN) {
        LOG_P_SC_GIVEN_REF <-
            cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)
    }
    #######################################
    rownames(LOG_P_SC_GIVEN_REF) <- colname_ref
    colnames(LOG_P_SC_GIVEN_REF) <- colname_sc
    ######2019.02.16 start ######
    LOG_P_SC_GIVEN_REF[which(is.na(LOG_P_SC_GIVEN_REF))] <- min(LOG_P_SC_GIVEN_REF)
    ######2019.02.16 end ######
    return(LOG_P_SC_GIVEN_REF)

}


.one_get_corr <- function(barcode, exp_sc_mat, exp_ref_mat, colname_ref,
                          method, verbose, print_step) {
    # calculate a correlation
    exp_sc <- as.array(exp_sc_mat[, barcode])
    log_p_sc_given_ref_list <- numeric(length = length(colname_ref))
    j <- 1
    while (j <= length(colname_ref)) {
        exp_ref <- as.array(exp_ref_mat[, j])
        if (method == 'kendall') {
            log_p_sc_given_ref <- cor.fk(exp_sc, exp_ref)
        } else {
            if (method == 'cosine') {
                log_p_sc_given_ref <-
                    sum(exp_sc * exp_ref) / sqrt(sum(exp_sc ^ 2) * sum(exp_ref ^ 2))
            } else {
                log_p_sc_given_ref <- cor(exp_sc, exp_ref, method = method)
            }
        }
        log_p_sc_given_ref_list[j] <- log_p_sc_given_ref
        j <- j + 1
    }
    ################################
    if (verbose) {
        if (i %% print_step == 1) {
            print(i)
        }
    }

    gc()
    return(log_p_sc_given_ref_list)

}


.get_cor <- function(exp_sc_mat, exp_ref_mat, method = 'kendall', num_threads = 4,
                     print_step = 100, gene_overlap = FALSE, verbose = FALSE){
    #method = "pearson", "kendall", "spearman"
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name

    #################
    library(parallel, verbose = F)
    #Step 1. get overlapped genes
    if (gene_overlap) {
        if (verbose) {
            print('Gene number of exp_sc_mat:')
            print(nrow(exp_sc_mat))
            print('Gene number of exp_ref_mat:')
            print(nrow(exp_ref_mat))
        }
        out.overlap <- get_overlap_genes(exp_sc_mat, exp_ref_mat)
        exp_sc_mat <- out.overlap$exp_sc_mat
        exp_ref_mat <- out.overlap$exp_ref_mat
        if (verbose) {
            print('Number of overlapped genes:')
            print(nrow(exp_sc_mat))
        }
    }

    ###############
    colname_sc <- colnames(exp_sc_mat)
    colname_ref <- colnames(exp_ref_mat)
    exp_sc_mat <- as.matrix(exp_sc_mat)
    exp_ref_mat <- as.matrix(exp_ref_mat)
    #######################################

    #Step 2. calculate corr
    cl <- makeCluster(num_threads, outfile = '')
    clusterEvalQ(cl, library(pcaPP))
    RUN <- parLapply(
        cl = cl,
        dimnames(exp_sc_mat)[[2]],
        .one_get_corr,
        exp_sc_mat = exp_sc_mat,
        exp_ref_mat = exp_ref_mat,
        colname_ref = colname_ref,
        method = method,
        verbose = verbose,
        print_step = print_step
    )
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=num_threads)
    LOG_P_SC_GIVEN_REF <- c()
    for (log_p_sc_given_ref_list in RUN) {
        LOG_P_SC_GIVEN_REF <-
            cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)
    }
    #######################################
    rownames(LOG_P_SC_GIVEN_REF) <- colname_ref
    colnames(LOG_P_SC_GIVEN_REF) <- colname_sc
    ######2019.02.16 start ######
    LOG_P_SC_GIVEN_REF[which(is.na(LOG_P_SC_GIVEN_REF))] <- -1
    ######2019.02.16 end ######

    return(LOG_P_SC_GIVEN_REF)
}


.get_tag_max <- function(P_REF_GIVEN_SC) {
    RN <- rownames(P_REF_GIVEN_SC)
    CN <- colnames(P_REF_GIVEN_SC)
    TAG <- cbind(CN, rep('NA', length(CN)))
    i <- 1
    while (i <= length(CN)) {
        this_rn_index <- which(P_REF_GIVEN_SC[, i] == max(P_REF_GIVEN_SC[, i]))[1]
        TAG[i, 2] <- RN[this_rn_index]
        i <- i + 1
    }
    colnames(TAG) <- c('cell_id', 'tag')
    return(TAG)
}


.combine_tags <- function(df.tags1, df.tags2) {
    # concat reference pval and local pval
    pvalue1 <- df.tags1[, c('scRef.tag', 'pvalue')]
    names(pvalue1) <- c('scRef.tag.1', 'pvalue.1')
    pvalue2 <- df.tags2[, c('scRef.tag', 'pvalue')]
    names(pvalue2) <- c('scRef.tag.2', 'pvalue.2')
    pvalue <- merge(pvalue1, pvalue2, by = 'row.names')
    row.names(pvalue) <- pvalue$Row.names
    pvalue$Row.names <- NULL

    # select more confident tag
    mtx.tag <- as.matrix(pvalue[, c('scRef.tag.1', 'scRef.tag.2')])
    mtx.pval <- as.matrix(pvalue[, c('pvalue.1', 'pvalue.2')])
    mtx.rank <- apply(mtx.pval, 1, rank, ties.method = "first")
    tag.final <-
        apply(as.array(1:dim(mtx.tag)[1]), 1, function(i) {
            mtx.tag[i, mtx.rank[1, i]]
        })
    pval.final <-
        apply(as.array(1:dim(mtx.pval)[1]), 1, function(i) {
            mtx.pval[i, mtx.rank[1, i]]
        })
    tag.final <-
        data.frame(
            scRef.tag = tag.final,
            pvalue = pval.final,
            row.names = dimnames(mtx.tag)[[1]],
            stringsAsFactors = F
        )

    OUT <- list()
    OUT$pvalue <- pvalue
    OUT$tag.final <- tag.final
    return(OUT)

}


generate_ref <- function(exp_sc_mat, TAG, min_cell = 1, M = 'SUM',
                         refnames = FALSE ){
    M <- M
    # print(M)
    min_cell <- min_cell
    refnames <- refnames
    exp_sc_mat <- exp_sc_mat
    TAG <- TAG
    NewRef <- c()
    TAG[, 2] <- as.character(TAG[, 2])
    if (refnames == FALSE) {
        refnames <- names(table(TAG[, 2]))
    }
    else{
        refnames <- refnames
    }
    outnames <- c()
    for (one in refnames) {
        this_col <- which(TAG[, 2] == one)
        if (length(this_col) >= min_cell) {
            outnames <- c(outnames, one)
            if (length(this_col) > 1) {
                if (M == 'SUM') {
                    this_new_ref <- apply(exp_sc_mat[, this_col], 1, sum)
                } else{
                    this_new_ref <- apply(exp_sc_mat[, this_col], 1, mean)
                }
            }
            else{
                this_new_ref <- exp_sc_mat[, this_col]
            }
            NewRef <- cbind(NewRef, this_new_ref)
        }
    }
    if (is.null(dim(NewRef))) {
        return(NULL)
    }
    rownames(NewRef) <- rownames(exp_sc_mat)
    colnames(NewRef) <- outnames
    return(NewRef)
}

getDEgeneF <- function(esetm = NULL, group = NULL, pair = FALSE,
                        block = NULL, p_adj = "fdr", fpkm = T) {
    # limma function
    if (is.null(esetm)) {
        cat(
            "esetm: gene expression matrix",
            "group: factor: \"c\"/\"d\"",
            "pair: TRUE/FALSE*",
            "block: e.g.1 2 2 1 if paired; blank if not",
            "p_adj: p.adjust, fdr* ",
            "fpkm: TRUE/FALSE*",
            sep = "\n"
        )
    } else{
        library(limma, verbose = F)
        if (pair) {
            design <- model.matrix( ~ block + group)
        } else{
            design <- model.matrix( ~ group)
        }
        suppressWarnings(fit <- lmFit(esetm, design))
        if (fpkm) {
            suppressWarnings(fit <- eBayes(fit, trend = T, robust = T))
        } else{
            suppressWarnings(fit <- eBayes(fit))
        }
        x <- topTable(fit, number = nrow(esetm), adjust.method = p_adj,
                      coef = "group2")
        x <- x[!is.na(row.names(x)), ]
        x <- x[!duplicated(row.names(x)), ]
        return(x)
    }
}


.imoprt_outgroup <- function(out.group = 'MCA', normalization = T, use_RUVseq = T) {
    if (class(out.group)[1] %in% c("data.frame", "matrix")) {
        df.out.group <- out.group
    } else {
        if (class(out.group)[1] == "character") {
            if (out.group %in% c("MCA", "HCL")) {
                if (out.group == "MCA") {
                    data(df.MCA)
                    df.out.group <- df.MCA
                } else {
                    data(df.HCL)
                    df.out.group <- df.HCL
                }
            } else {
                if (file.exists(out.group)) {
                    file.out.group <- out.group
                    df.out.group <-
                        read.table(file.out.group, header = T, row.names = 1,
                                   sep = '\t', check.names = F)
                }
            }
        } else {
            stop('Error: incorrect input of outgroup')
        }
    }
    if (normalization) {
        library(Seurat, verbose = F)
        seurat.out.group <-
            CreateSeuratObject(counts = df.out.group, project = "out.group",
                               min.cells = 1, min.features = 5000)
        seurat.out.group <-
            NormalizeData(seurat.out.group, normalization.method = "LogNormalize",
                          scale.factor = 1e6, verbose = F)
        if (use_RUVseq) {
            # get stably expression genes
            seurat.out.group <- FindVariableFeatures(
                seurat.out.group, selection.method = "mvp",
                mean.cutoff = c(3, Inf), dispersion.cutoff = c(-0.05, 0.05), verbose = F)
        }

        return(seurat.out.group)

    } else {
        return(df.out.group)
    }

}


.find_markers_cell_cosg <- function(cell.ref, list_near_cell, exp_ref_mat,
                                    exp_ref_mat.cell, exp_ref_label,
                                    mtx.combat, percent.high.exp,
                                    mtx.combat.use, topN, i) {
    library(Seurat, verbose = F)
    library(COSG)
    cell <- cell.ref[i]
    # print(cell)
    vec.cell <- mtx.combat[, paste0('Ref.', cell)]
    vec.ref <- exp_ref_mat[, cell]
    vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
    # high expression genes
    genes.high <- names(vec.ref.high)
    # genes.high <- rownames(exp_ref_mat.cell)

    # diff in local reference and negative reference
    seurat.Ref.cell <- CreateSeuratObject(counts = exp_ref_mat.cell[genes.high,], project = "Ref")
    seurat.Ref.cell <- NormalizeData(seurat.Ref.cell, verbose = F)
    seurat.Ref.cell@meta.data$original.label <- exp_ref_label
    Idents(seurat.Ref.cell) <- exp_ref_label
    n.neighbor <- min(3, length(list_near_cell[[cell]]))
    markers <-
        cosg(seurat.Ref.cell, groups='all', assay='RNA',
             slot='data', mu=2, n_genes_user=topN*(n.neighbor+1))
    genes.ref <- markers$names[,cell]

    if (n.neighbor > 1) {
        i <- 1
        for (cell_near in rev(list_near_cell[[cell]][1:n.neighbor])) {
            sub_seurat <-
                subset(seurat.Ref.cell, subset = original.label %in% c(cell, cell_near),
                       feature = genes.ref)
            sub_markers <- cosg(sub_seurat, groups='all', assay='RNA',
                                slot='data', mu=1, n_genes_user=topN*(n.neighbor-i+1))
            sub_markers_genes <- sub_markers$names[,cell]
            genes.ref <- intersect(genes.ref, sub_markers_genes)
            coef <- 0.5*(1+n.neighbor-i)
            if (length(genes.ref) < (coef*topN)) {
                genes.ref <- sub_markers_genes
            }
            i <- i + 1
        }
        use.genes <- genes.ref
    } else {
        if (length(genes.ref) > (topN)) {
            df_marker <- data.frame(name = markers$names[,cell], score=markers$scores[,cell])
            colnames(df_marker) <- c('name', 'score')
            df_marker <- df_marker[order(df_marker$score, decreasing = T),]
            genes.ref <- df_marker$name[1:(topN)]
        }
        use.genes <- genes.ref
    }


    # diff in Atlas
    if (!is.null(mtx.combat.use)) {
        if (length(use.genes) > 4*topN) {
            mtx.limma <- cbind(mtx.combat.use, vec.cell)
            bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
            res.limma.MCA <- getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
            res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
            genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
            if (length(genes.diff) < (topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:(topN)]
            }
            if (length(genes.diff) > (3*topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:(3*topN)]
            }
        } else {
            genes.diff <- use.genes
        }
    }else {
        genes.diff <- use.genes
    }

    return(genes.diff)

}


.find_markers_cell <- function(cell.ref, list_near_cell, exp_ref_mat = exp_ref_mat,
                               exp_ref_mat.cell, exp_ref_label,
                               mtx.combat, percent.high.exp,
                               mtx.combat.use, topN, i) {
    library(Seurat, verbose = F)
    cell <- cell.ref[i]
    # print(cell)
    vec.cell <- mtx.combat[, paste0('Ref.', cell)]
    # vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, percent.high.exp)]
    vec.ref <- exp_ref_mat[, cell]
    vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
    # high expression genes
    # genes.high <- intersect(names(vec.cell.high), names(vec.ref.high))
    genes.high <- names(vec.ref.high)
    # genes.high <- rownames(exp_ref_mat.cell)

    # diff in local reference and negative reference
    seurat.Ref.cell <- CreateSeuratObject(counts = exp_ref_mat.cell, project = "Ref")
    seurat.Ref.cell <- NormalizeData(seurat.Ref.cell, verbose = F)
    seurat.Ref.cell@meta.data$original.label <- exp_ref_label
    markers <- FindMarkers(seurat.Ref.cell, ident.1 = cell, group.by = 'original.label',
                           only.pos = T, features = genes.high, min.cells.group = 1,
                           min.pct = 0.1, min.diff.pct = 0,
                           logfc.threshold = 0.4, verbose = F)
    # markers$p_val_fdr <- p.adjust(markers$p_val, method = 'fdr')
    genes.ref.cell <- row.names(markers[(markers$p_val_adj < 0.05),])
    genes.ref <- genes.ref.cell

    # genes.ref = genes.high
    n.neighbor <- min(3, length(list_near_cell[[cell]]) - 1)
    if (n.neighbor > 1) {
        i <- 1
        for (cell_near in rev(list_near_cell[[cell]][1:n.neighbor])) {
            sub_seurat <- subset(seurat.Ref.cell, subset = original.label %in% c(cell, cell_near))
            sub_markers <- FindMarkers(sub_seurat, ident.1 = cell, group.by = 'original.label',
                                       only.pos = T, features = genes.ref.cell, min.cells.group = 1,
                                       min.pct = 0.1, min.diff.pct = 0,
                                       logfc.threshold = 0.25, verbose = F)
            sub_markers <- sub_markers[order(sub_markers$p_val),]
            genes.ref <- intersect(genes.ref, row.names(sub_markers[(sub_markers$p_val < 0.05),]))
            coef <- 0.5*(1+n.neighbor-i)
            if (length(genes.ref) < (coef*topN)) {
                sub_markers <- sub_markers[order(sub_markers$p_val),]
                genes.ref <- row.names(sub_markers)[1:min(coef*topN, nrow(sub_markers))]
            }
            i <- i + 1
        }
        if (length(genes.ref) > (topN)) {
            sub_markers <- sub_markers[order(sub_markers$p_val),]
            genes.ref <- row.names(sub_markers)[1:(topN)]
        }
        use.genes <- genes.ref
    } else {
        if (length(genes.ref) > (topN)) {
            sub_markers <- markers[order(markers$p_val),]
            genes.ref <- row.names(sub_markers)[1:(topN)]
        }
        use.genes <- genes.ref
    }


    # diff in Atlas
    if (!is.null(mtx.combat.use)) {
        if (length(use.genes) > 4*topN) {
            mtx.limma <- cbind(mtx.combat.use, vec.cell)
            bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
            res.limma.MCA <- getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
            res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
            genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
            if (length(genes.diff) < (topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:(topN)]
            }
            if (length(genes.diff) > (3*topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:(3*topN)]
            }
        } else {
            genes.diff <- use.genes
        }
    }else {
        genes.diff <- use.genes
    }

    return(genes.diff)

}

.find_markers_first <- function(exp_ref_mat, exp_ref_mat.cell, exp_ref_label, seurat.out.group,
                                type_ref = 'sum-counts', use_RUVseq = T, method_findmarker = 'COSG',
                                base.topN = 100, percent.high.exp = 0.8, num_threads = 6) {
    library(parallel, verbose = F)
    library(Seurat, verbose = F)
    # check parameters
    topN <- base.topN
    # if (!is.null(base.topN)) {
    #     base.topN <- base.topN
    # } else {
    #     stop('Error in finding markers: provide incorrect parameters')
    # }

    # transform count to fpm
    if (type_ref == 'sum-counts') {
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                    scale.factor = 1e6, verbose = F)
        exp_ref_mat <- as.matrix(seurat.Ref@assays$RNA@data)
    }
    if (type_ref %in% c('fpkm', 'tpm', 'rpkm', 'bulk')) {
        exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
    }
    exp_ref_mat <- na.omit(exp_ref_mat)
    cell.ref <- dimnames(exp_ref_mat)[[2]]

    if (!is.null(seurat.out.group)) {
        ###### regard a outgroup (e.g. MCA/HCA) as reference of DEG
        seurat.MCA <- seurat.out.group
        fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)

        # overlap genes
        fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
        out.overlap <- get_overlap_genes(fpm.MCA, exp_ref_mat)
        fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
        exp_ref_mat <- as.matrix(out.overlap$exp_ref_mat)
        if (use_RUVseq) {
            gene_overlap <- out.overlap$gene_over
            SEG.MCA <- VariableFeatures(seurat.MCA)
            gene.constant <- intersect(gene_overlap, SEG.MCA)
        }
        # print('Number of overlapped genes:')
        # print(nrow(exp_ref_mat))

        cell.MCA <- dimnames(fpm.MCA)[[2]]
        # use RUVseq to remove batch effect
        mtx.in <- cbind(fpm.MCA, exp_ref_mat)
        names.mix <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
        dimnames(mtx.in)[[2]] <- names.mix
        if (use_RUVseq) {
            suppressPackageStartupMessages(library(RUVSeq, verbose = F))
            seqRUVg <- RUVg(as.matrix(mtx.in), gene.constant, k=1, isLog = T)
            mtx.combat <- seqRUVg$normalizedCounts
        } else {
            mtx.combat <- mtx.in
        }
        mtx.combat.use <- mtx.combat[, paste0('MCA.', cell.MCA)]
    } else {
        mtx.combat <- exp_ref_mat
        names.mix <- paste0('Ref.', cell.ref)
        dimnames(mtx.combat)[[2]] <- names.mix
        mtx.combat.use <- NULL
    }

    mat_cor <- cor(exp_ref_mat)
    mat_cor[is.na(mat_cor)] <- 0
    # library(pcaPP)
    # cor.fk(exp_ref_mat)
    list_near_cell <- list()
    for (cell in cell.ref) {
        vec_cor <- mat_cor[cell,][mat_cor[cell,] > 0.7]
        if (length(vec_cor) < 2) {
            list_near_cell[[cell]] <- c()
        } else {
            list_near_cell[[cell]] <- names(sort(vec_cor, decreasing = T))[2:min(length(vec_cor), 4)]
        }
    }

    if (method_findmarker == 'COSG') {
        RUN <- mclapply(
            X = 1:length(cell.ref),
            FUN = .find_markers_cell_cosg,
            cell.ref = cell.ref, list_near_cell = list_near_cell,
            exp_ref_mat = exp_ref_mat,
            exp_ref_mat.cell = exp_ref_mat.cell, exp_ref_label = exp_ref_label,
            mtx.combat = mtx.combat, percent.high.exp = percent.high.exp,
            mtx.combat.use = mtx.combat.use, topN = topN,
            mc.cores = num_threads
        )
    }

    if (method_findmarker == 'Seurat') {
        cl = makeCluster(num_threads, outfile = '', type = 'FORK')
        clusterExport(cl, 'getDEgeneF')
        RUN <- parLapply(
            cl = cl,
            1:length(cell.ref),
            .find_markers_cell,
            cell.ref = cell.ref, list_near_cell = list_near_cell,
            exp_ref_mat = exp_ref_mat,
            exp_ref_mat.cell = exp_ref_mat.cell, exp_ref_label = exp_ref_label,
            mtx.combat = mtx.combat, topN = base.topN, percent.high.exp = percent.high.exp,
            mtx.combat.use = mtx.combat.use
        )
        stopCluster(cl)
    }
    names(RUN) <- cell.ref
    # for (cell in cell.ref) {
    #     if (length(RUN[[cell]] == 0)) {
    #         stop('Error: Failed to find marker genes! Please check whether you correctly install relevant packages.')
    #     }
    # }

    out <- list()
    out[['list.cell.genes']] <- RUN
    out[['list_near_cell']] <- list_near_cell
    return(out)

}


.find_markers_sc_cell <- function(cell.ref, list_near_cell, mtx.combat, LocalRef.sum, percent.high.exp,
                                  select.exp, vec.tag1, list.localNeg, mtx.combat.use, topN, i) {
    library(Seurat, verbose = F)
    cell <- cell.ref[i]
    # print(cell)
    vec.cell <- mtx.combat[, paste0('Ref.', cell)]
    # vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, percent.high.exp)]
    # vec.ref <- LocalRef.sum[, cell]
    # vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
    # high expression genes
    genes.high <- names(select.exp)

    # diff in local reference and negative reference
    seurat.Ref.cell <- CreateSeuratObject(counts = select.exp, project = "Ref")
    seurat.Ref.cell <- NormalizeData(seurat.Ref.cell, verbose = F)
    seurat.Ref.cell@meta.data$original.label <- vec.tag1
    markers <- FindMarkers(seurat.Ref.cell, ident.1 = cell, group.by = 'original.label',
                           only.pos = T, features = genes.high, min.cells.group = 1,
                           min.pct = 0.1, min.diff.pct = 0,
                           logfc.threshold = 0.4, verbose = F)
    # markers$p_val_fdr <- p.adjust(markers$p_val, method = 'fdr')
    genes.ref.cell <- row.names(markers[(markers$p_val_adj < 0.05),])
    genes.ref <- genes.ref.cell

    n.neighbor <- min(3, length(list_near_cell[[cell]]))
    i <- 1
    for (cell_near in rev(list_near_cell[[cell]][1:n.neighbor])) {
        sub_seurat <- subset(seurat.Ref.cell, subset = original.label %in% c(cell, cell_near))
        sub_markers <- FindMarkers(sub_seurat, ident.1 = cell, group.by = 'original.label',
                                   only.pos = T, features = genes.ref.cell, min.cells.group = 1,
                                   min.pct = 0.1, min.diff.pct = 0,
                                   logfc.threshold = 0.25, verbose = F)
        sub_markers <- sub_markers[order(sub_markers$p_val),]
        genes.ref <- intersect(genes.ref, row.names(sub_markers[(sub_markers$p_val < 0.05),]))
        coef <- 0.5*(1+n.neighbor-i)
        if (length(genes.ref) < (coef*topN)) {
            sub_markers <- sub_markers[order(sub_markers$p_val),]
            genes.ref <- row.names(sub_markers)[1:min(coef*topN, nrow(sub_markers))]
        }
        i <- i + 1
    }

    if (length(genes.ref) > (2*topN)) {
        sub_markers <- markers[order(markers$p_val),]
        genes.ref <- row.names(sub_markers)[1:(2*topN)]
    }
    use.genes <- genes.ref

    # diff in Atlas
    if (!is.null(mtx.combat.use)) {
        if (length(use.genes) > 4*topN) {
            mtx.limma <- cbind(mtx.combat.use, vec.cell)
            bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
            res.limma.MCA <- getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
            res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
            genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
            if (length(genes.diff) < (topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:min(topN, nrow(res.limma.MCA))]
            }
            if (length(genes.diff) > (2*topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:(2*topN)]
            }
        } else {
            genes.diff <- use.genes
        }
    } else {
        genes.diff <- use.genes
    }

    return(genes.diff)

}


.find_markers_sc_cell_cosg <- function(cell.ref, list_near_cell, mtx.combat, LocalRef.sum, percent.high.exp,
                                  select.exp, vec.tag1, list.localNeg, mtx.combat.use, topN, i) {
    library(Seurat, verbose = F)
    library(COSG)
    cell <- cell.ref[i]
    # print(cell)
    vec.cell <- mtx.combat[, paste0('Ref.', cell)]
    # vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, percent.high.exp)]
    vec.ref <- LocalRef.sum[, cell]
    vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
    # high expression genes
    # genes.high <- names(select.exp)

    # diff in local reference and negative reference
    seurat.Ref.cell <- CreateSeuratObject(counts = select.exp, project = "Ref")
    seurat.Ref.cell <- NormalizeData(seurat.Ref.cell, verbose = F)
    seurat.Ref.cell@meta.data$original.label <- vec.tag1
    Idents(seurat.Ref.cell) <- vec.tag1
    n.neighbor <- min(3, length(list_near_cell[[cell]]))
    markers <-
        cosg(seurat.Ref.cell, groups='all', assay='RNA',
             slot='data', mu=1, n_genes_user=topN*(n.neighbor+1))
    genes.ref <- markers$names[,cell]

    if (n.neighbor > 1) {
        i <- 1
        for (cell_near in rev(list_near_cell[[cell]][1:n.neighbor])) {
            sub_seurat <-
                subset(seurat.Ref.cell, subset = original.label %in% c(cell, cell_near),
                       feature = genes.ref)
            sub_markers <- cosg(sub_seurat, groups='all', assay='RNA',
                                slot='data', mu=1, n_genes_user=topN*(n.neighbor-i+1))
            sub_markers_genes <- sub_markers$names[,cell]
            genes.ref <- intersect(genes.ref, sub_markers_genes)
            coef <- 0.5*(1+n.neighbor-i)
            if (length(genes.ref) < (coef*topN)) {
                genes.ref <- sub_markers_genes
            }
            i <- i + 1
        }
        use.genes <- genes.ref
    } else {
        if (length(genes.ref) > (topN)) {
            df_marker <- data.frame(name = markers$names[,cell], score=markers$scores[,cell])
            colnames(df_marker) <- c('name', 'score')
            df_marker <- df_marker[order(df_marker$score, decreasing = T),]
            genes.ref <- df_marker$name[1:(topN)]
        }
        use.genes <- genes.ref
    }

    # diff in Atlas
    if (!is.null(mtx.combat.use)) {
        if (length(use.genes) > 4*topN) {
            mtx.limma <- cbind(mtx.combat.use, vec.cell)
            bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
            res.limma.MCA <- getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
            res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
            genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
            if (length(genes.diff) < (topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:min(topN, nrow(res.limma.MCA))]
            }
            if (length(genes.diff) > (2*topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:(2*topN)]
            }
        } else {
            genes.diff <- use.genes
        }
    } else {
        genes.diff <- use.genes
    }

    return(genes.diff)

}


.find_markers_second <- function(select.exp, vec.tag1, LocalRef,
                                 seurat.out.group, list.localNeg, method_findmarker = 'COSG',
                                 use_RUVseq = T, num_threads = 6,
                                 base.topN = 50, percent.high.exp = 0.80) {
    library(parallel, verbose = F)
    library(Seurat, verbose = F)
    # check parameters
    topN <- base.topN
    # if (!is.null(base.topN)) {
    #     base.topN <- base.topN
    # } else {
    #     stop('Error in finding markers: provide incorrect parameters')
    # }

    # transform count to fpm
    seurat.Ref.sum <- CreateSeuratObject(counts = LocalRef, project = "Ref")
    seurat.Ref.sum <- NormalizeData(seurat.Ref.sum, verbose = F)
    LocalRef.sum <- as.matrix(seurat.Ref.sum@assays$RNA@data)
    cell.ref <- dimnames(LocalRef.sum)[[2]]

    if (!is.null(seurat.out.group)) {
        ###### regard a outgroup (e.g. MCA/HCL) as reference of DEG
        seurat.MCA <- seurat.out.group
        fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)

        # overlap genes
        fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
        out.overlap <- get_overlap_genes(fpm.MCA, LocalRef.sum)
        fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
        LocalRef.sum <- as.matrix(out.overlap$exp_ref_mat)
        if (use_RUVseq) {
            gene_overlap <- out.overlap$gene_over
            SEG.MCA <- VariableFeatures(seurat.MCA)
            gene.constant <- intersect(gene_overlap, SEG.MCA)
        }
        # print('Number of overlapped genes:')
        # print(nrow(exp_ref_mat))

        cell.MCA <- dimnames(fpm.MCA)[[2]]
        cell.ref <- dimnames(LocalRef.sum)[[2]]
        # use RUVseq to remove batch effect
        mtx.in <- cbind(fpm.MCA, LocalRef.sum)
        names.mix <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
        dimnames(mtx.in)[[2]] <- names.mix
        if (use_RUVseq) {
            library(RUVSeq, verbose = F)
            seqRUVg <- RUVg(as.matrix(mtx.in), gene.constant, k=3, isLog = T)
            mtx.combat <- seqRUVg$normalizedCounts
        } else {
            mtx.combat <- mtx.in
        }
        mtx.combat <- mtx.in
        mtx.combat.use <- mtx.combat[, paste0('MCA.', cell.MCA)]
    } else {
        mtx.combat <- LocalRef.sum
        names.mix <- paste0('Ref.', cell.ref)
        dimnames(mtx.combat)[[2]] <- names.mix
        mtx.combat.use <- NULL
    }

    mat_cor <- cor(LocalRef.sum)
    mat_cor[is.na(mat_cor)] <- 0
    # library(pcaPP)
    # cor.fk(exp_ref_mat)
    list_near_cell <- list()
    for (cell in cell.ref) {
        vec_cor <- mat_cor[cell,][mat_cor[cell,] > 0.7]
        if (length(vec_cor) < 2) {
            list_near_cell[[cell]] <- c()
        } else {
            list_near_cell[[cell]] <- names(sort(vec_cor, decreasing = T))[2:min(length(vec_cor), 4)]
        }
    }

    if (method_findmarker == 'COSG') {
        RUN <- mclapply(
            X = 1:length(cell.ref),
            FUN = .find_markers_sc_cell_cosg,
            cell.ref = cell.ref, list_near_cell = list_near_cell, mtx.combat = mtx.combat,
            LocalRef.sum = LocalRef.sum, percent.high.exp = percent.high.exp,
            select.exp = select.exp, vec.tag1 = vec.tag1,
            list.localNeg = list.localNeg, mtx.combat.use = mtx.combat.use, topN = topN,
            mc.cores = num_threads
        )
    }

    if (method_findmarker == 'Seurat') {
        cl = makeCluster(num_threads, outfile = '')
        clusterExport(cl, 'getDEgeneF')
        RUN <- parLapply(
            cl = cl,
            1:length(cell.ref),
            .find_markers_sc_cell,
            cell.ref = cell.ref, list_near_cell = list_near_cell, mtx.combat = mtx.combat,
            LocalRef.sum = LocalRef.sum, topN = 50, percent.high.exp = percent.high.exp,
            select.exp = select.exp, vec.tag1 = vec.tag1,
            list.localNeg = list.localNeg, mtx.combat.use = mtx.combat.use
        )
        stopCluster(cl)
    }
    names(RUN) <- cell.ref
    # for (cell in cell.ref) {
    #     if (length(RUN[[cell]] == 0)) {
    #         stop('Error: Failed to find marker genes! Please check whether you correctly install relevant packages.')
    #     }
    # }

    out <- list()
    out[['list.cell.genes']] <- RUN
    out[['exp_ref_mat']] <- LocalRef.sum
    out[['list_near_cell']] <- list_near_cell
    return(out)

}


.cluster_sc <- function(exp_sc_mat, cluster_num_pc = 75, cluster_resolution = 3) {
    # cluster
    library(Seurat, verbose = F)
    # data preparing
    seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
    seurat.unlabeled <-
        NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize",
                      scale.factor = 10000, verbose = F)
    # print(seurat.unlabeled@assays$RNA@data[1:6,1:6])
    seurat.unlabeled <-
        FindVariableFeatures(seurat.unlabeled, selection.method = "vst",
                             nfeatures = 2000, verbose = F)
    seurat.unlabeled <- ScaleData(seurat.unlabeled, verbose = F)

    # PCA
    seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = cluster_num_pc, verbose = F)

    # cluster
    seurat.unlabeled <-
        FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:cluster_num_pc,
                      nn.eps = 0.5, verbose = F)
    seurat.unlabeled <-
        FindClusters(seurat.unlabeled, resolution = cluster_resolution,
                     n.start = 10, n.iter = 100, verbose = F)

    out.cluster <-
        data.frame(
            cluster.id = as.character(seurat.unlabeled@meta.data$seurat_clusters),
            # original.label = seurat.unlabeled@meta.data$original.label,
            row.names = dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
        )
    list_out <- list()
    list_out$out.cluster <- out.cluster
    list_out$seurat.query <- seurat.unlabeled
    return(list_out)

}


.cluster_increase_speed <- function(exp_sc_mat, df.cluster, combine_num_cell = 5, num_threads = 5) {
    library(parallel, verbose = F)
    sc.genes <- row.names(exp_sc_mat)
    df.cluster <- as.matrix(df.cluster)
    cluster.ids <- as.character(unique(df.cluster))
    combine_num_cell <- combine_num_cell

    .merge.one.cluster <- function(cluster.id) {
        library(Seurat)
        # merge cells in one cluster
        cell.ids <- names(df.cluster[df.cluster[, 'cluster.id'] == cluster.id,])
        sub.exp <- exp_sc_mat[, cell.ids]
        # print(dim(sub.exp))

        sub.seurat <- CreateSeuratObject(counts = sub.exp)
        sub.seurat <-
            NormalizeData(sub.seurat, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = F)
        # print(sub.seurat@assays$RNA@data[1:6,1:6])
        sub.seurat <-
            FindVariableFeatures(sub.seurat, selection.method = "disp",
                                 nfeatures = 1000, verbose = F)
        # print('1')
        sub.seurat <- ScaleData(sub.seurat,
                                features = VariableFeatures(sub.seurat),
                                verbose = F)
        # print(head(sub.seurat@assays$RNA@var.features))
        # print('2')

        # PCA
        num.cell <- dim(sub.exp)[2]
        if (num.cell < 5) {
            stop('Error: You should provide a smaller resolution!')
        }
        sub.seurat <- RunPCA(sub.seurat, npcs = min(10, (num.cell-2)), verbose = F)
        sub.pca <- sub.seurat@reductions$pca@cell.embeddings

        num.clusters <- ceiling(num.cell / combine_num_cell)
        if (num.clusters < 2) {
            sub.dict <- data.frame(
                cell.id = dimnames(sub.pca)[[1]],
                cluster.level1 = rep(cluster.id, num.cell),
                cluster.level2 = rep('1', num.cell)
            )
        } else {
            res.cluster <-
                kmeans(sub.pca, centers = num.clusters, nstart = 20, iter.max = 50)
            # names
            sub.dict <- data.frame(
                cell.id = names(res.cluster$cluster),
                cluster.level1 = rep(cluster.id, num.cell),
                cluster.level2 = res.cluster$cluster
            )
        }
        sub.dict$cluster.merge.id <-
            paste(sub.dict$cluster.level1, sub.dict$cluster.level2, sep = '-')
        row.names(sub.dict) <- sub.dict$cell.id

        # function
        .generate_ref <- function(exp_sc_mat, TAG, min_cell = 1, M = 'SUM',
                                  refnames = FALSE ){
            M <- M
            # print(M)
            min_cell <- min_cell
            refnames <- refnames
            exp_sc_mat <- exp_sc_mat
            TAG <- TAG
            NewRef <- c()
            TAG[, 2] <- as.character(TAG[, 2])
            if (refnames == FALSE) {
                refnames <- names(table(TAG[, 2]))
            }
            else{
                refnames <- refnames
            }
            outnames <- c()
            for (one in refnames) {
                this_col <- which(TAG[, 2] == one)
                if (length(this_col) >= min_cell) {
                    outnames <- c(outnames, one)
                    if (length(this_col) > 1) {
                        if (M == 'SUM') {
                            this_new_ref <- apply(exp_sc_mat[, this_col], 1, sum)
                        } else{
                            this_new_ref <- apply(exp_sc_mat[, this_col], 1, mean)
                        }
                    }
                    else{
                        this_new_ref <- exp_sc_mat[, this_col]
                    }
                    NewRef <- cbind(NewRef, this_new_ref)
                }
            }
            if (is.null(dim(NewRef))) {
                return(NULL)
            }
            rownames(NewRef) <- rownames(exp_sc_mat)
            colnames(NewRef) <- outnames
            return(NewRef)
        }

        # merge expression profile
        tag.in <- sub.dict[, c('cell.id', 'cluster.merge.id')]
        sub.exp.merge <- .generate_ref(sub.exp, tag.in)
        sub.exp.merge <- sub.exp.merge[sc.genes, ]
        # print(class(sub.exp.merge))
        if (class(sub.exp.merge)[1] %in% c('numeric', 'integer')) {
            sub.exp.merge <- as.data.frame(sub.exp.merge)
            names(sub.exp.merge) <- unique(sub.dict$cluster.merge.id)
        }

        sub.out <- list()
        sub.out$sub.dict <- sub.dict
        sub.out$sub.exp.merge <- sub.exp.merge
        return(sub.out)

    }


    # split dataset
    cl.input <- list()
    cl <- makeCluster(num_threads, outfile = '')
    # cl <- makeCluster(num_threads, outfile = '', type = 'FORK')
    # clusterExport(cl, '.generate_ref')
    # clusterEvalQ(cl, library(Seurat))
    out.par <- parLapply(
        cl = cl,
        cluster.ids,
        .merge.one.cluster
    )
    stopCluster(cl)

    df.dict <- data.frame(stringsAsFactors = F)
    df.exp.merge <- data.frame(stringsAsFactors = F)
    i = 1
    for (sub.out in out.par) {
        df.dict <- rbind(df.dict, sub.out$sub.dict)
        if (i == 1) {
            df.exp.merge <- sub.out$sub.exp.merge
        } else {
            df.exp.merge <- cbind(df.exp.merge, sub.out$sub.exp.merge)
        }
        i = i + 1
    }
    seurat.exp.merge <-
        CreateSeuratObject(counts = df.exp.merge, project = "exp.merge",
                           min.cells = 1, min.features = 1)
    seurat.exp.merge <-
        NormalizeData(seurat.exp.merge, normalization.method = "LogNormalize", verbose = F)
    df.exp.merge <- seurat.exp.merge@assays$RNA@data

    out.merge <- list()
    out.merge$df.dict <- df.dict
    out.merge$df.exp.merge <- df.exp.merge
    out.merge$out.par <- out.par
    return(out.merge)

}


.confirm_label_auc <- function(exp_sc_mat, list.cell.genes, list_near_cell, scRef.tag, num_threads = 10) {
    library(parallel, verbose = F)
    library(AUCell, verbose = F)
    set.seed(123)
    # confirm label
    exp_sc_mat <- as.matrix(exp_sc_mat)
    df.tag <- as.data.frame(scRef.tag)
    row.names(df.tag) <- df.tag$cell_id
    df.tag$cell_id <- NULL
    names(df.tag) <- 'scRef.tag'

    total_genes <- unique(unlist(list.cell.genes))
    df.auc <- data.frame()
    list.auc.back <- list()
    for (cell in names(list.cell.genes)) {
        if (sum(df.tag$scRef.tag==cell) < 1) {
            next
        }
        sel_cols <- colnames(exp_sc_mat)[(colnames(exp_sc_mat)%in%(colnames(exp_sc_mat)[df.tag$scRef.tag==cell]))]
        total_mtx <- exp_sc_mat[total_genes, sel_cols]
        if (is.null(dim(total_mtx))) {
            total_mtx <- matrix(total_mtx, ncol=1)
            dimnames(total_mtx)[[1]] <- total_genes
            dimnames(total_mtx)[[2]] <- sel_cols
        } else {
            total_mtx <- as.matrix(total_mtx)
        }
        cells_rankings <- AUCell::AUCell_buildRankings(
            total_mtx, verbose = F, plotStats = F)

        cell_markers <- list.cell.genes[[cell]]
        # cells_near <- list_near_cell[[cell]]
        # genes_near <- c()
        # if (length(cells_near) > 0) {
        #     for (i in 1:length(cells_near)) {
        #         genes_near <- c(genes_near, list.cell.genes[[cells_near[i]]])
        #     }
        # } else {
        #     other_ref <- setdiff(names(list.cell.genes), cell)
        #     for (c in sample(other_ref, min(2, length(other_ref)))) {
        #         genes_near <- c(genes_near, list.cell.genes[[c]])
        #     }
        # }

        # genes_near <- list.cell.genes[[cells_near[1]]]
        # genes_near <- unique(c(list.cell.genes[[cells_near[1]]], list.cell.genes[[cells_near[2]]]))
        # num_genes <- min(length(cell_markers), length(genes_near))
        # genes.marker <- sample(cell_markers, num_genes)
        # genes.back <- sample(genes_near, num_genes)
        # genes.marker <- cell_markers
        # genes.back <- genes_near

        genes.marker = cell_markers
        geneSets <- list(geneSet1=genes.marker)
        cells_AUC <- AUCell::AUCell_calcAUC(
            geneSets, cells_rankings, aucMaxRank = ceiling(length(total_genes)*0.05), verbose = F)

        # df.tags1_cd8 <- df.tags1[colnames(exp_sc_mat)%in%
        #                              (colnames(exp_sc_mat)[tag1[,2]==cell]),]
        # df.tags1_cd8$auc <- unlist(cells_AUC@assays@data@listData)
        # View(df.tags1_cd8)
        # background
        # geneSets.back <- list(geneSet1=genes.back)
        # cells_AUC.back <- AUCell_calcAUC(geneSets.back, cells_rankings,
        #                                  aucMaxRank = length(genes.back)/2,
        #                                  verbose = F)
        # list.auc.back[[cell]] <- unlist(cells_AUC.back@assays@data@listData)
        df.auc <- rbind(df.auc, data.frame(AUC = unlist(cells_AUC@assays@data@listData),
                                           # AUC_back = unlist(cells_AUC.back@assays@data@listData),
                                           row.names = colnames(total_mtx)))

    }
    meta.tag <- merge(df.tag, df.auc, by = 'row.names')
    row.names(meta.tag) <- meta.tag$Row.names
    meta.tag$Row.names <- NULL
    # meta.tag$log10Pval <- -log10(meta.tag$pvalue)

    gc()

    list_out <- list()
    list_out$meta.tag <- meta.tag
    list_out$list.auc.back <- list.auc.back
    return(list_out)

}


.one_cutoff_AUC <- function(i, cells, df.tags1, auc_gap, list_tags1_back) {
    library(mclust, verbose = F)
    cell <- cells[i]
    df.sub <- df.tags1[df.tags1$scRef.tag == cell, ]
    sub_AUC <- df.sub[, c('AUC', 'cluster.id')]
    # cell_back <- list_tags1_back[[cell]]
    # back_value <- min(quantile(cell_back,0.9), 0.05)
    df_cluster_median = aggregate(sub_AUC$AUC, by = list(sub_AUC$cluster.id), median)
    colnames(df_cluster_median) <- c('cluster.id', 'median')
    prop_0 <- sum(df.sub$AUC == 0) / dim(df.sub)[1]
    if (sum(df.tags1$scRef.tag == cell) < 6 | prop_0 > 0.95) {
        one_out <- list()
        one_out[[1]] <- cell
        one_out[[2]] <- 0.3 + auc_gap
        one_out[[3]] <- 0.3 + auc_gap
    } else {
        diff_cluster <- max(df_cluster_median$median) - min(df_cluster_median$median)
        diff_cell <- max(df.sub$AUC) - min(df.sub$AUC)
        if (diff_cluster > 0.25 | diff_cell > 0.5) {
            min_G <- 4
            max_G <- 6
        } else {
            min_G <- 2
            max_G <- 5
        }
        model_AUC <- densityMclust(sub_AUC[,1], verbose = F, G = min_G:max_G)
        mean_AUC <- model_AUC$parameters$mean
        len_mean <- length(mean_AUC)
        if (length(table(model_AUC$classification)) != length(mean_AUC)) {
            mean_AUC <- mean_AUC[names(table(model_AUC$classification))]
        }
        diff_AUC <- c()
        for (i in 1:(length(mean_AUC)-1)) {
            diff_AUC <- c(diff_AUC, mean_AUC[i+1] - mean_AUC[i])
        }
        # if (nrow(df.sub) < 200) {
        #     auc_gap <- max(auc_gap, 0.2)
        # }
        auc_gap_tight <- 0.01
        if (sum(diff_AUC > auc_gap_tight, na.rm = T) > 0) {
            cut_class <- max(as.numeric(names(mean_AUC)[1:length(diff_AUC)])[diff_AUC > auc_gap_tight])
            for (i in cut_class:len_mean) {
                if (!(as.character(cut_class+1) %in% names(mean_AUC))) {cut_class = cut_class+1}
            }
            # cut_AUC <- mean_AUC[cut_class]
            df.sub$class_AUC <- model_AUC$classification
            # cut_AUC <- max(df.sub$AUC[df.sub$class_AUC == as.character(cut_class)])
            # cut_AUC <- quantile(df.sub$AUC[df.sub$class_AUC == as.character(cut_class)], 0.99)
            cutoff_auc <- median(df.sub$AUC[df.sub$class_AUC == as.character(cut_class+1)])
        } else {
            cutoff_auc <- 0
        }

        auc_gap_loose <- auc_gap
        if (sum(diff_AUC > auc_gap_loose, na.rm = T) > 0) {
            cut_class <- max(as.numeric(names(mean_AUC)[1:length(diff_AUC)])[diff_AUC > auc_gap_loose])
            # cut_AUC <- mean_AUC[cut_class]
            df.sub$class_AUC <- model_AUC$classification
            # cut_AUC <- max(df.sub$AUC[df.sub$class_AUC == as.character(cut_class)])
            cut_AUC <- quantile(df.sub$AUC[df.sub$class_AUC == as.character(cut_class)], 0.99)
        } else {
            cut_AUC <- 0
        }

        df_cluster_median <- df_cluster_median[order(df_cluster_median$median, decreasing = F),]
        rownames(df_cluster_median) <- 1:nrow(df_cluster_median)
        df_cluster_median$diff <- rep(0, nrow(df_cluster_median))
        for (i in 1:(nrow(df_cluster_median)-1)) {
            df_cluster_median[i, 'diff'] <-
                df_cluster_median[i+1, 'median'] - df_cluster_median[i, 'median']
        }
        auc_gap_cluster <- max(auc_gap-0.03-nrow(df_cluster_median)/1500, 0.05)
        if (sum(df_cluster_median$diff > auc_gap_cluster, na.rm = T) > 0) {
            cut_class <- max(as.numeric(rownames(df_cluster_median))[df_cluster_median$diff > auc_gap_cluster])
            drop_cluster <- df_cluster_median$cluster.id[1:cut_class]
            cut_cluster_AUC <- quantile(df.sub$AUC[df.sub$cluster.id == drop_cluster[cut_class]], 0.9)
        } else {
            cut_cluster_AUC <- 0
        }
        if (nrow(df_cluster_median) >= 20 & sum(df_cluster_median$diff > 0.06) > 0) {
            cut_class <- max(as.numeric(rownames(df_cluster_median))[df_cluster_median$diff > 0.06])
            drop_cluster <- df_cluster_median$cluster.id[1:cut_class]
            cut_class <- min(cut_class, round(nrow(df_cluster_median)/2))
            cut_diff <- df_cluster_median$diff[cut_class]
            if (cut_diff > 2*max(df_cluster_median$diff[(cut_class+1):nrow(df_cluster_median)])) {
                cut_cluster_AUC <- max(quantile(df.sub$AUC[df.sub$cluster.id == drop_cluster[cut_class]], 0.9), cut_cluster_AUC)
            }
        }

        cut_AUC <- max(cut_AUC, cut_cluster_AUC)

        one_out <- list()
        one_out[[1]] <- cell
        one_out[[2]] <- cutoff_auc
        # one_out[[2]] <- c(max(quantile(df.sub_sel$AUC, 0.25), back_value),
        #                   max(quantile(df.sub_sel$diff, 0.25), 0.05))
        one_out[[3]] <- cut_AUC
        # one_out[[3]] <- median(cell_back)
    }
    print(cell)

    return(one_out)
}


.cutoff_AUC <- function(df.tags1, list_tags1_back, exp_sc_mat, threshold,
                        method_findmarker = method_findmarker, num_threads = num_threads) {
    library(parallel, verbose = F)
    cells <- as.character(unique(df.tags1$scRef.tag))
    vec.cutoff <- c()
    vec.neg.cutoff <- c()
    vec.cut <- c()
    if (method_findmarker == 'Seurat') {
        base_thre <- 0.18
    }
    if (method_findmarker == 'COSG') {
        base_thre <- 0.2
    }
    if (threshold <= 5) {
        auc_gap <- (5-threshold)*2 + base_thre
    } else {
        auc_gap <- (5-threshold) + base_thre
    }

    RUN <- mclapply(
        X = 1:length(cells),
        FUN = .one_cutoff_AUC,
        cells = cells,
        df.tags1 = df.tags1,
        auc_gap = auc_gap,
        list_tags1_back = list_tags1_back,
        mc.cores = num_threads
    )

    for (one_out in RUN) {
        cell <- one_out[[1]]
        vec.cutoff <- c(vec.cutoff, one_out[[2]])
        # list.cutoff <- c(quantile(df.sub_sel$AUC, 0.25), quantile(df.sub_sel$diff, 0.25))
        # vec.neg.cutoff <- c(vec.neg.cutoff, one_out[[3]])
        vec.cut <- c(vec.cut, one_out[[3]])
    }
    names(vec.cutoff) <- cells
    # names(vec.neg.cutoff) <- cells
    names(vec.cut) <- cells

    out.cutoff <- list()
    out.cutoff$vec.cutoff <- vec.cutoff
    out.cutoff$vec.neg.cutoff <- vec.neg.cutoff
    out.cutoff$vec.cut <- vec.cut
    return(out.cutoff)
}


scMAGIC <- function(exp_sc_mat, exp_ref_mat, exp_ref_label = NULL,
                    single_round = F, identify_unassigned = T,
                    atlas = NULL, use_RUVseq = T,
                    method_findmarker = 'COSG',
                    percent_high_exp = 0.7, num_marker_gene = 100,
                    cluster_num_pc = 50, cluster_resolution = 3,
                    combine_num_cell = NULL, min_cell = 1,
                    method1 = 'kendall', method2 = NULL,
                    method_HVGene = 'scibet',
                    corr_use_HVGene1 = 2000, corr_use_HVGene2 = 2000,
                    threshold = 5, num_threads = 4, cluster_assign = F,
                    simple_output = T) {

    library(parallel, verbose = F)
    library(Matrix)
    if (method_HVGene == 'scibet') {
        library(scibet)
    }
    if (method_HVGene == 'SciBet_R') {
        library(scibetR)
        library(reticulate)
        np <- import("numpy")
        np.exp2 <- np$exp2
        np.max <- np$max
        np.sum <- np$sum
    }
    if (method_HVGene == 'Seurat') {
        suppressPackageStartupMessages(library(Seurat))
    }
    # check parameters
    # if (!type_ref %in% c('sc-counts', 'sum-counts', 'fpkm', 'tpm', 'rpkm')) {
    #     stop('Error: inexistent input of reference data format')
    # }

    class_query <- class(exp_sc_mat)[[1]]
    if (class_query %in% c("data.frame", "dgCMatrix", "matrix", "dgTMatrix")) {
        if (class_query == "data.frame") {
            exp_sc_mat <- as(as.matrix(exp_sc_mat), 'dgCMatrix')
        }
        if (class_query == "matrix") {
            exp_sc_mat <- as(exp_sc_mat, 'dgCMatrix')
        }
    } else {
        stop("Error: incorrect data format is input to 'exp_sc_mat'! ")
    }
    class_ref <- class(exp_ref_mat)[[1]]
    if (class_ref %in% c("data.frame", "dgCMatrix", "matrix", "dgTMatrix")) {
        if (class_ref == "data.frame") {
            exp_ref_mat <- as(as.matrix(exp_ref_mat), 'dgCMatrix')
        }
        if (class_ref == "matrix") {
            exp_ref_mat <- as(exp_ref_mat, 'dgCMatrix')
        }
    } else {
        stop("Error: incorrect data format is input to 'exp_ref_mat'! ")
    }

    cutoff.1 = 'default'
    cutoff.2 = 'default'
    mod = ''
    # simple_output = T
    type_ref = 'sc-counts'
    out.group = atlas
    opt_speed = F
    num_cell <- ncol(exp_sc_mat)
    if (is.null(combine_num_cell)) {
        if (num_cell > 3000) {
            combine_num_cell = 5
        } else {
            combine_num_cell = 3
        }
    }

    time1 <- Sys.time()
    # get sum-counts format
    if (type_ref == 'sc-counts') {
        print('Sum single cell counts matrix:')
        label.in <- data.frame(cell_id = colnames(exp_ref_mat),
                               tag = as.character(exp_ref_label))
        exp_ref_mat.cell <- exp_ref_mat
        exp_ref_mat.sum <- generate_ref(exp_ref_mat, label.in, M = 'SUM')
        exp_ref_mat <- exp_ref_mat.sum
        type_ref <- 'sum-counts'
    }

    # get overlap genes
    out.overlap <- get_overlap_genes(exp_sc_mat, exp_ref_mat)
    exp_sc_mat <- out.overlap$exp_sc_mat
    exp_ref_mat <- out.overlap$exp_ref_mat
    gene_over <- out.overlap$gene_over
    # print('Number of overlapped genes:')
    # print(nrow(exp_sc_mat))

    if (!is.null(out.group)) {
        seurat.out.group <- .imoprt_outgroup(out.group = out.group, use_RUVseq = use_RUVseq)
        # overlap genes
        gene.overlap <- intersect(gene_over, rownames(seurat.out.group@assays$RNA@counts))
        exp_sc_mat <- exp_sc_mat[gene.overlap, ]
        exp_ref_mat <- exp_ref_mat[gene.overlap, ]
        exp_ref_mat.cell <- exp_ref_mat.cell[gene.overlap,]
    } else {
        seurat.out.group <- NULL
    }

    print('Number of overlapped genes:')
    print(nrow(exp_sc_mat))

    # cluster analysis
    print('Start clustering :')
    list_out.cluster <-
        .cluster_sc(exp_sc_mat,
                    cluster_num_pc = cluster_num_pc,
                    cluster_resolution = cluster_resolution)
    df.cluster <- list_out.cluster$out.cluster
    seurat.query <- list_out.cluster$seurat.query

    # out.merge <-
    #     .cluster_increase_speed(exp_sc_mat, df.cluster,
    #                             combine_num_cell = combine_num_cell, num_threads = num_threads)
    # df.dict <- out.merge$df.dict
    print('Clustering completed!')
    # speed calculation
    if (opt_speed) {
        print('Speed calculation by clustering:')
        # out.merge <-
        #     .cluster_increase_speed(exp_sc_mat, df.cluster,
        #                             combine_num_cell = combine_num_cell, num_threads = num_threads)
        # df.dict <- out.merge$df.dict
        df.exp.merge <- out.merge$df.exp.merge
        df.dict.merge <- unique(df.dict[,c('cluster.level1', 'cluster.merge.id')])
        rownames(df.dict.merge) <- df.dict.merge$cluster.merge.id
        df.dict.merge <- df.dict.merge[colnames(df.exp.merge),]
        min_cell <- ceiling(min_cell / combine_num_cell)
        # query_set <- as.data.frame(t(df.exp.merge))/1.0
        # query_set$label <- df.dict.merge$cluster.level1
    } else {
        df.exp.merge <- exp_sc_mat
    }

    gene_not0 <- rownames(df.exp.merge)[Matrix::rowSums(df.exp.merge)!=0]
    df.exp.merge <- df.exp.merge[gene_not0, ]
    exp_ref_mat <- exp_ref_mat[gene_not0, ]
    exp_ref_mat.cell <- exp_ref_mat.cell[gene_not0, ]

    # find markers of cell types in reference
    topN = num_marker_gene
    percent.high.exp = percent_high_exp
    print('Find marker genes of cell types in reference:')
    out.markers <-
        .find_markers_first(
            exp_ref_mat, exp_ref_mat.cell, exp_ref_label,
            seurat.out.group,
            type_ref = 'sum-counts',
            use_RUVseq = use_RUVseq,
            method_findmarker = method_findmarker, base.topN = topN,
            percent.high.exp = percent.high.exp, num_threads = num_threads
        )
    list.cell.genes <- out.markers[['list.cell.genes']]
    list_near_cell <- out.markers[['list_near_cell']]
    df.exp.merge <- exp_sc_mat

    # rm(exp_sc_mat)
    gc()
    df.exp.merge <- as.matrix(df.exp.merge)
    if (method1 == 'multinomial') {
        query_set <- as.data.frame(t(as.matrix(df.exp.merge)))/1.0
        query_set$label <- df.cluster$cluster.id
    }

    print('First-round annotation:')
    # print(method1)
    if (method1 == 'xxx (based PC)') {
        # pred_tags <- func_sciBet(df.exp.merge, exp_ref_mat.cell, exp_ref_label)
        # tag1 <- data.frame(cell_id = colnames(df.exp.merge), tag = pred_tags)
        # tag1 <- as.matrix(tag1)

        # pred_tags <- func_scClassify(df.exp.merge, exp_ref_mat.cell, exp_ref_label)
        # # pred_tags=func_scPred(df.exp.merge, exp_ref_mat.cell, exp_ref_label)
        # # pred_tags=func_seurat(df.exp.merge, exp_ref_mat.cell, exp_ref_label)
        # table(label_sc, pred_tags)
    } else {
        if (!is.null(corr_use_HVGene1)) {
            if (method_HVGene == 'scibet') {
                num_feature <- corr_use_HVGene1/2
                train_set <- as.data.frame(t(as.matrix(exp_ref_mat.cell)))/1.0
                train_set$label <- exp_ref_label
                HVG <- scibet::SelectGene(train_set, k = num_feature)
            }
            if (method_HVGene == 'SciBet_R') {
                num_feature <- corr_use_HVGene1/2
                train_set <- as.data.frame(t(as.matrix(exp_ref_mat.cell)))/1.0
                train_set$label <- exp_ref_label
                HVG <- scibetR::SelectGene_R(train_set, k = num_feature)
            }
            if (method_HVGene == 'Seurat') {
                HVG <- .get_high_variance_genes(exp_ref_mat, num.genes = corr_use_HVGene1)
            }
            similarity.in <- df.exp.merge[HVG, ]
            ref.in <- exp_ref_mat[HVG, ]
        } else {
            similarity.in <- df.exp.merge
            ref.in <- exp_ref_mat
        }
        if (method1 != 'multinomial') {
            out1 <- .get_cor(ref.in, similarity.in, method = method1, num_threads = num_threads)
            out1 <- t(out1)
            tag1 <- .get_tag_max(out1)
        } else {
            if (!(method_HVGene %in% c('scibet', 'SciBet_R'))) {
                train_set <- as.data.frame(t(as.matrix(exp_ref_mat.cell)))/1.0
            }
            train_set <- train_set[, HVG]
            train_set$label <- exp_ref_label
            test_set <- query_set[, HVG]
            if (method_HVGene == 'SciBet_R') {
                pred_tags <- scibetR::SciBet_R(train_set, test_set, k=num_feature, result = 'list')
                out1 <- scibetR::SciBet_R(train_set, test_set, k=num_feature, result = 'table')
            } else {
                pred_tags <- scibet::SciBet(train_set, test_set, k=num_feature, result = 'list')
                out1 <- scibet::SciBet(train_set, test_set, k=num_feature, result = 'table')
            }
            rownames(out1) <- rownames(test_set)
            out1 <- t(out1)
            tag1 <- data.frame(cell_id = rownames(test_set), tag = pred_tags)
            tag1 <- as.matrix(tag1)
        }
    }
    # out1_scale <- scale(out1)
    # out1_diff <- apply(out1_scale, 2, function(input) {sort(input, decreasing = T)[1] - sort(input, decreasing = T)[2]})

    gc()

    print('Build local reference')
    list_out_1 <- .confirm_label_auc(df.exp.merge, list.cell.genes, list_near_cell,
                                     tag1, num_threads = num_threads)
    df.tags1 <- list_out_1$meta.tag
    # list_tags1_back <- list_out_1$list.auc.back
    cell_ids <- colnames(df.exp.merge)
    df.tags1 <- df.tags1[cell_ids, ]
    # df.tags1$celltype = label_sc
    # df.tags1$diff <- out1_diff
    df.tags1 <- merge(df.tags1, df.cluster, by = 'row.names')
    rownames(df.tags1) <- df.tags1$Row.names
    df.tags1$Row.names <- NULL
    out.cutoff <- .cutoff_AUC(df.tags1, list_tags1_back, exp_sc_mat, threshold,
                              method_findmarker = method_findmarker, num_threads = num_threads)
    df.cutoff.1 <- out.cutoff$vec.cutoff
    neg.cutoff.1 <- out.cutoff$vec.neg.cutoff
    vec.cut_1 <- out.cutoff$vec.cut

    if (single_round) {
        pvalue1 <- data.frame(stringsAsFactors = F)
        for (one_cell in unique(df.tags1$scRef.tag)) {
            df_sub <- df.tags1[df.tags1$scRef.tag == one_cell, c('scRef.tag', 'AUC')]
            df_sub$tag_unassaigned <- df_sub$scRef.tag
            df_sub$tag_unassaigned[df_sub$AUC < vec.cut_1[one_cell]] <- 'Unassigned'
            pvalue1 <- rbind(pvalue1, df_sub)
        }
        names(pvalue1) <- c('scRef.tag.1', 'AUC.1', 'tag_unassaigned.1')
        df_unassigned <- pvalue1[cell_ids, ]
        df_unassigned$tag.final <- pvalue1$scRef.tag.1
        # df_unassigned$tag.final[df_unassigned$tag_unassaigned.1 == 'Unassigned'] <- 'Unassigned'
        if (identify_unassigned) {
            df_unassigned$tag.final[df_unassigned$tag_unassaigned.1 == 'Unassigned'] <- 'Unassigned'
        }

        df.tags <- merge(df_unassigned, df.cluster, by = 'row.names')
        row.names(df.tags) <- df.tags$Row.names
        df.tags$Row.names <- NULL
        df.tags <- df.tags[, c("tag.final", "cluster.id")]
        df.tags$scRef.tag <- df.tags$tag.final
        all_sub_clusters <- unique(df.tags$cluster.id)
        for (sub_cluster in all_sub_clusters) {
            df_sub <- df.tags[df.tags$cluster.id == sub_cluster,]
            sub_table <- sort(table(df_sub$tag.final), decreasing = T)
            if (length(sub_table) == 1) {
                next
            } else {
                if (sub_table[1]/nrow(df_sub) >= 0.8) {
                    df.tags$scRef.tag[df.tags$cluster.id == sub_cluster] <- names(sub_table)[1]
                }
            }
        }

        cell_ids <- colnames(exp_sc_mat)
        df.tags <- df.tags[cell_ids, ]

        if (cluster_assign) {
            all_clusters <- unique(df.tags$cluster.id)
            df.tags$cluster.tags <- df.tags$scRef.tag
            for (cluster in all_clusters) {
                df_sub <- df.tags[df.tags$cluster.id == cluster,]
                sub_table <- sort(table(df_sub$tag.final), decreasing = T)
                df.tags$cluster.tags[df.tags$cluster.id == cluster] <- names(sub_table)[1]
            }
            df.combine <- data.frame(scMAGIC.tag = df.tags$cluster.tags, row.names = rownames(df.tags))
        } else {
            df.combine <- data.frame(scMAGIC.tag = df.tags$scRef.tag, row.names = rownames(df.tags))
        }

        gc()

        time2 <- Sys.time()
        time.scRef <- difftime(time2, time1, units = 'secs')
        output <- list()
        output$tag1 <- tag1
        output$out1 <- out1
        output$combine.out <- df.tags
        output$dict.cluster <- df.cluster
        output$ref.markers <- list.cell.genes
        output$final.out <- df.combine
        output$run.time <- time.scRef

        print('Finish!')

        return(output)

    }

    select.barcode <- c()
    for (cell in names(df.cutoff.1)) {
        sub.cutoff <- df.cutoff.1[cell]
        sub.select <- df.tags1[df.tags1$scRef.tag == cell, ]
        sub.select <- sub.select[sub.select$AUC >= sub.cutoff[1], ]
        select.barcode <- c(select.barcode, row.names(sub.select))
    }
    list.localNeg <- list()

    df.tags1 <- df.tags1[cell_ids, ]
    df.tags1$conf <- rep(0, nrow(df.tags1))
    df.tags1$conf[rownames(df.tags1) %in% select.barcode] <- 1
    df.tags1$expand_conf <- rep(0, nrow(df.tags1))
    table_num <- table(df.tags1$scRef.tag[df.tags1$conf == 1])
    table_prop <- table_num/length(select.barcode)
    not_expend_1 <- names(table_num)[table_num > 300]
    not_expend_2 <- names(table_prop)[table_prop > 0.2]
    not_expend <- union(not_expend_1, not_expend_2)
    for (cell_name in names(table_prop)) {
        if (cell_name %in% not_expend) {
            df.tags1[df.tags1$conf == 1, 'expand_conf'] <- 1
        }
        df_sub <- df.tags1[df.tags1$scRef.tag == cell_name & df.tags1$conf == 1,]
        sub_table <- table(df_sub$cluster.id) / nrow(df_sub)
        expand_clust <- names(sub_table[sub_table > 0.3])
        df.tags1[df.tags1$cluster.id %in% expand_clust & df.tags1$scRef.tag == cell_name, 'expand_conf'] <- 1
        df.tags1[df.tags1$AUC <= vec.cut_1[cell_name] & df.tags1$scRef.tag == cell_name, 'expand_conf'] <- 0
    }
    select.barcode <- rownames(df.tags1)[df.tags1$expand_conf == 1]
    # df_conf=df.tags1[df.tags1$expand_conf==1,]
    select.exp <- df.exp.merge[, cell_ids %in% select.barcode]
    select.tag1 <- tag1[tag1[, 'cell_id'] %in% select.barcode, ]
    LocalRef <- generate_ref(select.exp, select.tag1,  min_cell = min_cell)
    vec.tag1 <- select.tag1[, 'tag']
    print('Cell types in local reference:')
    print(dimnames(LocalRef)[[2]])

    #####
    gc()
    #####
    # find local marker genes
    print('find local marker genes')
    out.markers <-
        .find_markers_second(
            select.exp, vec.tag1, LocalRef,
            seurat.out.group, list.localNeg,
            use_RUVseq = use_RUVseq,
            base.topN = topN, method_findmarker = method_findmarker,
            percent.high.exp = percent.high.exp,
            num_threads = num_threads
        )
    local.cell.genes <- out.markers[['list.cell.genes']]
    list_near_cell <- out.markers[['list_near_cell']]
    local.cell.genes_merge <- local.cell.genes
    # local.cell.genes_merge <- list()
    # for (one_cell in names(local.cell.genes)) {
    #     local.cell.genes_merge[[one_cell]] <-
    #         setdiff(union(list.cell.genes[[one_cell]], local.cell.genes[[one_cell]]), NA)
    # }

    num_local_cell <- length(vec.tag1)
    if (is.null(method2)) {
        if (num_local_cell > 5000) {
            method2 = 'randomforest'
        } else {
            method2 = 'multinomial'
        }
    }
    if (method2 == 'multinomial') {
        query_set <- as.data.frame(t(as.matrix(df.exp.merge)))/1.0
        query_set$label <- df.cluster$cluster.id
    }
    print('Second-round annotation:')
    if (method2 == 'randomforest') {
        pca_query <- as.data.frame(seurat.query@reductions$pca@cell.embeddings)
        train_pca <- pca_query[colnames(select.exp), ]
        train_pca$label <- as.factor(vec.tag1)
        library(randomForest)
        fit.forest <- randomForest(label ~ ., data = train_pca)
        pred_tags <- predict(fit.forest, pca_query)
        tag2 <- data.frame(cell_id = rownames(pca_query), tag = pred_tags)
        tag2 <- as.matrix(tag2)
    } else {
        if (!is.null(corr_use_HVGene2)) {
            if (method_HVGene == 'scibet') {
                num_feature <- corr_use_HVGene2/2
                train_set <- query_set[colnames(select.exp),]
                train_set$label <- vec.tag1
                HVG <- scibet::SelectGene(train_set, k = num_feature)
            }
            if (method_HVGene == 'SciBet_R') {
                num_feature <- corr_use_HVGene2/2
                train_set <- query_set[colnames(select.exp),]
                train_set$label <- vec.tag1
                HVG <- scibetR::SelectGene_R(train_set, k = num_feature)
            }
            if (method_HVGene == 'Seurat') {
                HVG <- .get_high_variance_genes(LocalRef, num.genes = corr_use_HVGene2)
            }
            similarity.in <- df.exp.merge[HVG, ]
            ref.in <- LocalRef[HVG, ]
        } else {
            similarity.in <- df.exp.merge
            ref.in <- exp_ref_mat
        }
        if (method2 != 'multinomial') {
            out2 <- .get_cor(ref.in, similarity.in, method = method2, num_threads = num_threads)
            out2 <- t(out2)
            tag2 <- .get_tag_max(t(out2))
        } else {
            if (!(method_HVGene %in% c('scibet', 'SciBet_R'))) {
                train_set <- query_set[colnames(select.exp),]
            }
            train_set <- train_set[, HVG]
            train_set$label <- vec.tag1
            test_set <- query_set[, HVG]
            if (method_HVGene == 'SciBet_R') {
                pred_tags <- scibetR::SciBet_R(train_set, test_set, k=num_feature, result = 'list')
                # out2 <- scibetR::SciBet_R(train_set, test_set, k=num_feature, result = 'table')
            } else {
                pred_tags <- scibet::SciBet(train_set, test_set, k=num_feature, result = 'list')
                out2 <- scibet::SciBet(train_set, test_set, k=num_feature, result = 'table')
                rownames(out2) <- rownames(test_set)
                out2 <- t(out2)
            }
            tag2 <- data.frame(cell_id = rownames(query_set), tag = pred_tags)
            tag2 <- as.matrix(tag2)
        }
    }

    gc()

    list_out_2 <- .confirm_label_auc(df.exp.merge, local.cell.genes_merge, list_near_cell,
                                     tag2, num_threads = num_threads)
    df.tags2 <- list_out_2$meta.tag
    list_tags2_back <- list_out_2$list.auc.back
    cell_ids <- colnames(df.exp.merge)
    df.tags2 <- df.tags2[cell_ids, ]
    gc()


    pvalue1 <- data.frame(stringsAsFactors = F)
    for (one_cell in unique(df.tags1$scRef.tag)) {
        df_sub <- df.tags1[df.tags1$scRef.tag == one_cell, c('scRef.tag', 'AUC')]
        df_sub$tag_unassaigned <- df_sub$scRef.tag
        df_sub$tag_unassaigned[df_sub$AUC < vec.cut_1[one_cell]] <- 'Unassigned'
        pvalue1 <- rbind(pvalue1, df_sub)
    }
    names(pvalue1) <- c('scRef.tag.1', 'AUC.1', 'tag_unassaigned.1')
    df_unassigned <- pvalue1[cell_ids, ]

    pvalue1 <- data.frame(stringsAsFactors = F)
    for (cell_label in unique(df.tags1$scRef.tag)) {
        df_sub <- df.tags1[df.tags1$scRef.tag == cell_label,]
        # df_sub$AUC_scale <- scale(df_sub$AUC)
        pvalue1 <- rbind(pvalue1, df_sub[, c('scRef.tag', 'AUC')])
    }
    names(pvalue1) <- c('scRef.tag.1', 'AUC.1')
    pvalue2 <- data.frame(stringsAsFactors = F)
    for (cell_label in unique(df.tags2$scRef.tag)) {
        df_sub <- df.tags2[df.tags2$scRef.tag == cell_label,]
        # df_sub$AUC_scale <- scale(df_sub$AUC)
        pvalue2 <- rbind(pvalue2, df_sub[, c('scRef.tag', 'AUC')])
    }
    names(pvalue2) <- c('scRef.tag.2', 'AUC.2')
    pvalue <- merge(pvalue1, pvalue2, by = 'row.names')
    row.names(pvalue) <- pvalue$Row.names
    pvalue$Row.names <- NULL
    pvalue <- pvalue[cell_ids, ]
    mtx.tag <- as.matrix(pvalue[, c('scRef.tag.1', 'scRef.tag.2')])
    mtx.pval <- as.matrix(pvalue[, c('AUC.1', 'AUC.2')])
    mtx.rank <- apply(mtx.pval, 1, rank, ties.method = "first")
    tag.merge <-
        apply(as.array(1:dim(mtx.tag)[1]), 1, function(i) {
            mtx.tag[i, mtx.rank[2, i]]
        })
    pvalue$tag.merge <- tag.merge
    df_auc_unassigned <- merge(pvalue, df_unassigned, by = 'row.names')
    rownames(df_auc_unassigned) <- df_auc_unassigned$Row.names
    df_auc_unassigned$Row.names <- NULL
    df_auc_unassigned <- df_auc_unassigned[cell_ids, ]
    df_auc_unassigned$tag.final <- tag.merge
    if (identify_unassigned) {
        df_auc_unassigned$tag.final[df_auc_unassigned$tag_unassaigned.1 == 'Unassigned'] <- 'Unassigned'
    }

    if (opt_speed) {
        df.cluster <- df.dict[, c("cluster.merge.id", "cluster.level1")]
        df.cluster <- unique(df.cluster)
        df.cluster <- data.frame(cluster.id = df.cluster$cluster.level1,
                                 row.names = df.cluster$cluster.merge.id,
                                 stringsAsFactors = F)
    }

    df.tags <- merge(df_auc_unassigned, df.cluster, by = 'row.names')
    row.names(df.tags) <- df.tags$Row.names
    df.tags$Row.names <- NULL
    df.tags <- df.tags[, c("tag.merge", "tag.final", "cluster.id")]
    df.tags$scRef.tag <- df.tags$tag.final
    all_sub_clusters <- unique(df.tags$cluster.id)
    for (sub_cluster in all_sub_clusters) {
        df_sub <- df.tags[df.tags$cluster.id == sub_cluster,]
        sub_table <- sort(table(df_sub$tag.final), decreasing = T)
        if (length(sub_table) == 1) {
            next
        } else {
            if (sub_table[1]/nrow(df_sub) >= 0.8) {
                df.tags$scRef.tag[df.tags$cluster.id == sub_cluster&df.tags$scRef.tag == 'Unassigned'] <- names(sub_table)[1]
            }
        }
    }

    # if (opt_speed) {
    #     df.tags$cluster.merge.id <- row.names(df.tags)
    #     df.tags.merge <- merge(df.tags, df.dict[, c('cluster.merge.id', 'cell.id')],
    #                            by = 'cluster.merge.id')
    #     df.tags.merge$cluster.merge.id <- NULL
    #     row.names(df.tags.merge) <- df.tags.merge$cell.id
    #     df.tags.merge$cell.id <- NULL
    #     df.tags <- df.tags.merge
    # }
    cell_ids <- colnames(exp_sc_mat)
    df.tags <- df.tags[cell_ids, ]

    if (cluster_assign) {
        all_clusters <- unique(df.tags$cluster.id)
        df.tags$cluster.tags <- df.tags$scRef.tag
        for (cluster in all_clusters) {
            df_sub <- df.tags[df.tags$cluster.id == cluster,]
            sub_table <- sort(table(df_sub$tag.final), decreasing = T)
            df.tags$cluster.tags[df.tags$cluster.id == cluster] <- names(sub_table)[1]
        }
        df.combine <- data.frame(scMAGIC.tag = df.tags$cluster.tags, row.names = rownames(df.tags))
    } else {
        df.combine <- data.frame(scMAGIC.tag = df.tags$scRef.tag, row.names = rownames(df.tags))
    }

    #####
    gc()
    #####
    time2 <- Sys.time()
    time.scMAGIC <- difftime(time2, time1, units = 'secs')
    if (simple_output) {
        output <- df.combine
    } else {
        output <- list()
        output$tag1 <- tag1
        output$out1 <- out1
        output$LocalRef <- LocalRef
        if (identify_unassigned) {
            output$pvalue1 <- df.tags1
            output$pvalue2 <- df.tags
            output$localref.cell <- select.exp
            output$localref.tag <- vec.tag1
            output$dict.cluster <- df.cluster
            # output$cutoff.1 <- df.cutoff.1
            # output$cutoff.neg.1 <- neg.cutoff.1
            # output$cutoff.2 <- df.cutoff.2
            output$ref.markers <- list.cell.genes
            output$local.markers <- local.cell.genes
        }
        output$final.out <- df.combine
        output$run.time <- time.scMAGIC
        if (mod == 'debug') {
            output$df.exp.merge <- df.exp.merge
            output$diff.log10Pval <- diff.log10Pval
        }
    }

    print('Finish!')

    return(output)

}


transformHomoloGene <- function(exp_sc_mat, inTaxID = 9606, outTaxID = 10090) {
    library(homologene)
    genes.in <- rownames(exp_sc_mat)
    res.home <- homologene(genes.in, inTax = inTaxID, outTax = outTaxID)
    res.home <- res.home[!duplicated(res.home[, 1]),]
    res.home <- res.home[!duplicated(res.home[, 2]),]
    genes.out <- res.home[, 1]
    genes.homo <- res.home[, 2]
    exp.out <- exp_sc_mat[genes.out,]
    rownames(exp.out) <- genes.homo
    return(exp.out)
}


scMAGIC_Seurat <- function(seurat.query, seurat.ref, atlas = 'MCA',
                           method_findmarker = 'COSG',
                           percent_high_exp = 0.7, num_marker_gene = 100,
                           cluster_num_pc = 50, cluster_resolution = 3,
                           corr_use_HVGene = 2000, min_cell = 1,
                           method1 = 'kendall', method2 = NULL,
                           threshold = 5, num_threads = 4, cluster_assign = F) {
    exp_sc_mat <- seurat.query@assays$RNA@counts
    exp_ref_mat <- seurat.ref@assays$RNA@counts
    exp_ref_label <- seurat.ref@meta.data$celltype
    output.scMAGIC <- scMAGIC(exp_sc_mat, exp_ref_mat, exp_ref_label,
                              atlas = atlas,
                              corr_use_HVGene1 = corr_use_HVGene,
                              corr_use_HVGene2 = corr_use_HVGene,
                              method_findmarker = method_findmarker,
                              percent_high_exp = percent_high_exp,
                              num_marker_gene = num_marker_gene,
                              cluster_num_pc = cluster_num_pc,
                              cluster_resolution = cluster_resolution,
                              min_cell = min_cell,
                              method1 = method1, method2 = method2,
                              cluster_assign = cluster_assign,
                              threshold = threshold,
                              num_threads = num_threads)
    pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
    seurat.query$prediction_celltype <- pred.scMAGIC
    return(seurat.query)
}

