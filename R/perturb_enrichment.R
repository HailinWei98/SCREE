#' @import Mixscape
#' @import ggplot2
#' @import patchwork
#' @import scales
#' @import dplyr
#' @import reshape2
#' @import plyr
#' @export

CalculatePerturbEnrichment <- function(mtx_dir, sg_dir, NTC = "NTC", prefix = "./", label = "", 
                                       plot.save = TRUE, table.save = TRUE, top = NULL, NTC.cal = TRUE){
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        eccite <- readRDS(mtx_dir)
    } else {
        eccite <- mtx_dir
    }
    
    if (is.character(sg_dir)) {
        message(paste("Reading sgRNA lib file:", sg_dir))
        sg_lib <- read.table(sg_dir, header = T)
    } else {
        sg_lib <- sg_dir
    }
    
    #replicate information
    
    if(!("replicate" %in% colnames(eccite@meta.data))){
        stop("Please add replicate information to the SeuratObject first.")
    }
    
    #perturbation information
    
    if(!("perturbations" %in% colnames(eccite@meta.data))){
        stop("Please add perturbation information to the SeuratObject first.")
    }
    
    if (NTC.cal == FALSE) {
        eccite <- subset(eccite, perturbations != NTC)
    }

    ##Generate propotion matrix
    
    sg_lib <- subset(sg_lib, cell %in% colnames(eccite))
    prtb <- unique(sg_lib$gene)
    pro_mat <- matrix(rep(0, length(prtb) * length(unique(eccite$seurat_clusters))), 
                      nrow = length(prtb))
    rownames(pro_mat) <- prtb
    colnames(pro_mat) <- seq(0, length(unique(eccite$seurat_clusters)) - 1)
    
    #Generate 0-1 matrix
    
    G_mat <- frame2indmatrix(sg_lib, eccite)
    
    G_mat <- as.data.frame(G_mat)
    G_mat$cluster <- eccite$seurat_clusters
    
    for (i in colnames(G_mat)) {
        if (i != "cluster") {
            pro_mat[i, ] <- as.vector(tapply(G_mat[, i], INDEX = G_mat$cluster, FUN = sum))
        }
    }
    
    cluster_pro <- t(pro_mat)/as.vector(table(G_mat$cluster))
    
    bk <- seq(0, 1, by = 0.05)
    #bk <- seq(0, 0.2, by = 0.01)
    lg_bk <- seq(0, 1, 0.5)
    #lg_bk <- seq(0, 0.2, 0.1)
#     if (nrow(pro_mat)/ncol(pro_mat) >= 10) {
#         if (max(cluster_pro) < 0.2) {
#             bk <- seq(0, 0.2, by = 0.01)
#             lg_bk <- seq(0, 0.2, 0.1)
#         }
#     }

    
    cluster_pro_plot <- cluster_pro
    if (is.null(top)) {
        p1 <- pheatmap::pheatmap(cluster_pro_plot, treeheight_col = F, treeheight_row = F, 
                                 main = "Perturbations Enrichment",
                                 colorRampPalette(colors = c("white", "coral1"))(length(bk)),
                                 legend_breaks = lg_bk, breaks = bk, show_rownames = T, fontsize = 15, 
                                 cellwidth = 15, cellheight = 20,
                                 show_colnames = T, silent = T)
    } else {
        top_prtb <- names(head(sort(apply(cluster_pro_plot, 2, max), decreasing = TRUE), top))
        
        if (NTC.cal == TRUE) {
            if (NTC %in% top_prtb) {
                top_prtb <- names(head(sort(apply(cluster_pro_plot, 2, max), decreasing = TRUE), top + 1))
            } else if (NTC %in% colnames(cluster_pro_plot)){
                top_prtb <- c(top_prtb, NTC)
            }
        }

        cluster_pro_plot <- cluster_pro_plot[, top_prtb]
        p1 <- pheatmap::pheatmap(cluster_pro_plot, treeheight_col = F, treeheight_row = F,
                                 main = "Perturbations Enrichment",
                                 colorRampPalette(colors = c("white", "coral1"))(length(bk)),
                                 legend_breaks = lg_bk, breaks = bk, show_rownames = T, fontsize = 15, 
                                 cellwidth = 15, cellheight = 20,
                                 show_colnames = T, silent = T)
    }
    
    if (table.save == TRUE) {
        write.table(cluster_pro, file = file.path(prefix, paste(label, "perturb_ratio.txt", sep = "")), quote = FALSE)
    }
    
    if (plot.save == TRUE) {
        
        dir <- file.path(prefix, "pdf")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }

        dir <- file.path(dir, "perturbation_efficiency")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }

        img_dir <- file.path(prefix, "img")
        if (!(dir.exists(img_dir))) {
            dir.create(img_dir)
        }

        img_dir <- file.path(img_dir, "perturbation_efficiency")
        if (!(dir.exists(img_dir))) {
            dir.create(img_dir)
        }
        
        file <- file.path(dir, paste(label, "perturb_ratio.pdf", sep = ""))
        if (ncol(cluster_pro_plot) > 20) {
            if (nrow(cluster_pro_plot) > 20) {
                cellheight = 15
                cellwidth = 15
            } else {
                cellheight = 20
                cellwidth = 15
            }
            size = 15
        } else {
            cellheight = NA
            cellwidth = NA
            size = 10
        }
        pheatmap(cluster_pro_plot, treeheight_col = F, treeheight_row = F, main = "Perturbations Enrichment",
                 colorRampPalette(colors = c("white", "#E64B35FF"))(length(bk)),
                 legend_breaks = lg_bk, breaks = bk, show_rownames = T, 
                 show_colnames = T, fontsize = size, 
                 cellwidth = cellwidth, cellheight = cellheight, filename = file)
        
        img_file <- file.path(img_dir, paste(label, "perturb_ratio.png", sep = ""))
        pheatmap(cluster_pro_plot, treeheight_col = F, treeheight_row = F, main = "Perturbations Enrichment",
                 colorRampPalette(colors = c("white", "#E64B35FF"))(length(bk)),
                 legend_breaks = lg_bk, breaks = bk, show_rownames = T, 
                 show_colnames = T, fontsize = size, 
                 cellwidth = cellwidth, cellheight = cellheight, filename = img_file)
    }
    
    return(p1)
    
}