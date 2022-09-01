#' @import pheatmap
#' @import psych
#' @export

heatmap <- function(score_dir, pval_dir, score_cut = 0.5, p_val_cut = 0.05, NTC = "NegCtrl", remove_neg = TRUE,
                    num = 0, prefix= "./", label = "", plot.save = TRUE) {
    
    #get score and p_val
    
    if(is.character(score_dir)){
        score <- asread.table(score_dir, header = T, row.names = 1)
    }else{
        score <- score_dir
    }
    if(is.character(pval_dir)){
        p_val <- read.table(pval_dir, header = T, row.names = 1)
    }else{
        p_val <- pval_dir
    }
    
    #filter data frame
    
    if (remove_neg == TRUE) {
        remove_neg <- score[, -which(colnames(score) == NTC)]
    } else {
        remove_neg<-score
    }
    
    p_val_filter <- apply(p_val[, colnames(remove_neg)], 1, function(x) sum(x < p_val_cut))
    score_filter <- apply(remove_neg, 1, function(x) sum(abs(x) > score_cut))
    gene_count <- p_val_filter[p_val_filter > num & score_filter > num]
    filtered <- remove_neg[names(gene_count), ]
                          
    if (nrow(filtered) <= 3) {
        return(warning("Too few genes passed threshold, stop calculating correlation."))
    }
    
    #calculate correlation
                          
    sg_cor <- corr.test(filtered)
    
    bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
                          
    if (ncol(remove_neg) > 20) {
        cell = 15
        size = 15
    } else {
        cell = NA
        size = 10
    }
                          
    p <- pheatmap::pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F,
                  main = "Perturbations Correlation",
                  c(colorRampPalette(colors = c("#4DBBD5FF", "white"))(length(bk)/2),
                    colorRampPalette(colors = c("white", "#E64B35FF"))(length(bk)/2)),
                  legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
                  breaks = bk, show_rownames = T, show_colnames = T, fontsize = size, 
                  filename = NA, silent = TRUE)
                          
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

        file <- file.path(dir, paste(label, "correlation_heatmap.pdf", sep = ""))

        pheatmap::pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F, 
                 main = "Perturbations Correlation",
                 c(colorRampPalette(colors = c("#4DBBD5FF", "white"))(length(bk)/2),
                   colorRampPalette(colors = c("white", "#E64B35FF"))(length(bk)/2)),
                 legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
                 breaks = bk, show_rownames = T, show_colnames = T, fontsize = size, 
                 cellwidth = cell, cellheight = cell, filename = file, silent = FALSE)

        img_file <- file.path(img_dir, paste(label, "correlation_heatmap.png", sep = ""))
        #                  display_numbers = matrix(ifelse(sg_cor$p <= 0.01, "**", ifelse(sg_cor$p <= 0.05 , "*", " ")), nrow(sg_cor$p))
        pheatmap::pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F, 
                 main = "Perturbations Correlation",
                 c(colorRampPalette(colors = c("#4DBBD5FF", "white"))(length(bk)/2),
                   colorRampPalette(colors = c("white", "#E64B35FF"))(length(bk)/2)),
                 legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
                 breaks = bk, show_rownames = T, show_colnames = T, fontsize = size, 
                 cellwidth = cell, cellheight = cell, filename = img_file, silent = TRUE)
        
        return(p)
    } else {
        return(p)
    }
}