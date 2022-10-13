#' Correlation Heatmap
#'
#' Calculate correlation between perturbations based on regulatory score of a union set of potential target genes for all perturbations. The plot function is based on \code{\link[pheatmap]{pheatmap}}
#'
#' @param score Data frame or directory of score from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param pval Data frame or directory of p_value from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param score_cut Score cutoff of \code{improved_scmageck_lr} results. Default is 0.2.
#' @param pval_cut P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.
#' @param NTC The name of negative controls. Default is "NegCtrl".
#' @param remove_neg Logical, remove negative control before calculating correlation or not. Default is \code{TRUE}.
#' @param num Minimum number of perturbation passed score_cut or pval_cut for each gene. Only genes with at least \code{num} numbers of perturbation passed threshold will be retained for correlation calculation. Default is 0. 
#' @param min.filtered.gene Minimum number of genes retained for correlation calculation. Default is 3.
#' @param method Method to calculate correlation, can be one of "spearman", "kendall" or "pearson". Default is "pearson".
#' @param color Color for heatmap. Default is c("#4DBBD5FF", "white", "#E64B35FF"), corresponding to the correlation from -1 to 1. 
#' @param cell Cell width and height settings for heatmap. If set this parameter to a single number or \code{NA}, cell width and height will be the same, and \code{NA} means the cell will be adaptive to the input perturbation numbers. If set this parameter to a vector including two numbers or \code{NA}, the first one will be cell width and the second one will be height. Default is "auto", which means cell width and cell height will be set to 15 when the perturbations numbers > 20 (after removing negative control) while cell width and cell height will be set to NA when the perturbations numbers <= 20 (after removing negative control).
#' @param fontsize Fontsize of the heatmap. Default is "auto", which means the fontsize will be set to 15 when the perturbations numbers > 20 (after removing negative control) while the fontsize will be set to 10 when the perturbations numbers <= 20 (after removing negative control).
#' @param angle Angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315). Default is 90.
#' @param legend.title Include the title of legend in the heatmap or not. Default is FALSE.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is "atuo", according to the cellwidth unless cellwidth is \code{NA}. If cellwidth is \code{NA} or legend.title is \code{TRUE} and width is "auto", the width will be set to 10.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is "atuo", according to the cellheight unless cellheight is \code{NA}. If cellheight is \code{NA} or legend.title is \code{TRUE} and height is "auto", the height will be set to 8.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom utils read.table write.table
#' @import pheatmap
#' @importFrom psych corr.test
#' @import gtable
#' @import grid
#' @import ggplotify
#' @export

heatmap <- function(score, pval, score_cut = 0.2, pval_cut = 0.05, NTC = "NegCtrl", remove_neg = TRUE, num = 0, min.filtered.gene = 3, method = "pearson", color = c("#4DBBD5FF", "white", "#E64B35FF"), cell = "auto", fontsize = "auto", angle = 90, legend.title = FALSE, plot.save = TRUE, prefix= ".", label = "", width = "auto", height = "auto", png_res = 720) {
    
    #get score and p_val
    
    if (is.character(score)) {
        score <- read.table(score, header = T, row.names = 1)
    }
    
    if (is.character(pval)) {
        pval <- read.table(pval, header = T, row.names = 1)
    }
    
    #filter data frame
    
    if (remove_neg == TRUE) {
        remove_neg <- score[, -which(colnames(score) == NTC)]
    } else {
        remove_neg <- score
    }
    
    #A loose strategy
    
    # p_val_filter <- apply(pval[, colnames(remove_neg)], 1, function(x) sum(x < pval_cut))
    # score_filter <- apply(remove_neg, 1, function(x) sum(abs(x) > score_cut))
    # gene_count <- p_val_filter[p_val_filter > num & score_filter > num]
                          
    score_pval <- (pval < pval_cut) * (abs(score) > score_cut)
    gene_count <- rowSums(score_pval)
    gene_count <- gene_count[gene_count > num]
    
    filtered <- remove_neg[names(gene_count), ]
                          
    if (nrow(filtered) <= min.filtered.gene) {
        stop("Too few genes passed threshold, stop calculating correlation.")
    }
    
    #calculate correlation
                          
    sg_cor <- psych::corr.test(filtered, method = method)
    
    bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
                          
    if (length(cell) == 1) {
        if (cell != "auto" & !is.numeric(cell) & !is.na(cell)) {
            warning("Cell must be a numeric parameter or NA or 'auto', using 'auto' instead.")
            cell <- "auto"
        }
        
        if (!is.na(cell) & cell == "auto") {
            if (ncol(remove_neg) > 20) {
                cellwidth = 15
                cellheight = 15
            } else {
                cellwidth = NA
                cellheight = NA
            }
        } else {
            cellwidth = cell
            cellheight = cell
        }
        
    } else if (length(cell) == 2) {
        if ((!is.numeric(cell[1]) & !is.na(cell[1])) | (!is.numeric(cell[2]) & !is.na(cell[2]))) {
            warning("Cell must be a numeric parameter or NA or 'auto', using 'auto' instead.")
            if (ncol(remove_neg) > 20) {
                cellwidth = 15
                cellheight = 15
            } else {
                cellwidth = NA
                cellheight = NA
            }
        } else {
            cellwidth = cell[1]
            cellheight = cell[2]
        }
    }
    
    if (fontsize == "auto") {
        if (ncol(remove_neg) > 20) {
            fontsize = 15
        } else {
            fontsize = 10
        }
    }                          
                          
    p <- pheatmap::pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F,
                            main = "Perturbations Correlation",
                            c(colorRampPalette(colors = c(color[1], color[2]))(length(bk)/2),
                              colorRampPalette(colors = c(color[2], color[3]))(length(bk)/2)),
                            legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
                            cellwidth = cellwidth, cellheight = cellheight,
                            breaks = bk, show_rownames = T, show_colnames = T, fontsize = fontsize, 
                            filename = NA, silent = TRUE, angle_col = angle)
                          
    if (legend.title == TRUE) {
        
        main.grob <- p$gtable$grob[[1]]
        plot.grob <- p$gtable$grob[[2]]
        xlab.grob <- p$gtable$grob[[3]]  
        ylab.grob <- p$gtable$grob[[4]]  
        legend.grob <- p$gtable$grob[[5]]
        #Shift both down by 1 inch
        legend.grob$children[[1]]$y <- legend.grob$children[[1]]$y - unit(0.85, "inches") 
        legend.grob$children[[2]]$y <- legend.grob$children[[2]]$y - unit(0.85, "inches") 
        legend.grob$children[[1]]$x <- legend.grob$children[[1]]$x + unit(0.4, "inches") 
        legend.grob$children[[2]]$x <- legend.grob$children[[2]]$x + unit(0.4, "inches") 
        leg_label <- grid::textGrob("Pearson Correlation", x = 0, y = 0.9, hjust = 0, vjust = 0, 
                              gp = gpar(fontsize = 8, fontface = "bold"))
        legend.grob2 <- grid::addGrob(legend.grob, leg_label)
        
        my_new_gt <- gtable::gtable(widths = p$gtable$widths, height = p$gtable$height)
        my_new_gt$heights[4] <- sum(unit(1, "npc"), unit(-1, "grobheight", plot.grob), 
                                    unit(-10, "bigpts"), unit(-5, "bigpts"), unit(0, "bigpts"), unit(-3, "inches"))
        my_new_gt$widths[3] <- sum(unit(1, "npc"), unit(-1, "grobwidth", plot.grob), unit(-10, "bigpts"), 
                                   -1 * max(unit(1.1, "grobwidth", plot.grob), 
                                            sum(unit(12, "bigpts"), unit(1.32, "grobwidth", plot.grob))), 
                                   unit(-5, "bigpts"), unit(-4, "inches"))
        my_new_gt$widths[5] <- sum(max(unit(1.1, "grobwidth", legend.grob2), 
                                       sum(unit(12, "bigpts"), unit(1.32, "grobwidth", legend.grob2))), unit(3, "inches"))
        gtable <- gtable::gtable_add_grob(my_new_gt,main.grob, 1, 3)
        gtable <- gtable::gtable_add_grob(gtable,plot.grob, 4, 3)
        gtable <- gtable::gtable_add_grob(gtable,xlab.grob,5, 3)
        gtable <- gtable::gtable_add_grob(gtable,ylab.grob,4, 4)
        gtable <- gtable::gtable_add_grob(gtable,legend.grob2, 4, 5)
        p <- ggplotify::as.ggplot(gtable)

    }
                          
    if (plot.save == TRUE) {
        
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }

        dir <- file.path(dir, "perturbation_efficiency")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }

        pdf_dir <- file.path(dir, "pdf")
        if (!(dir.exists(pdf_dir))) {
            dir.create(pdf_dir)
        }
        
        img_dir <- file.path(dir, "img")
        if (!(dir.exists(img_dir))) {
            dir.create(img_dir)
        }

        file <- file.path(pdf_dir, paste(label, "correlation_heatmap.pdf", sep = ""))
        img_file <- file.path(img_dir, paste(label, "correlation_heatmap.png", sep = ""))
        
        if (!is.na(cellwidth) & legend.title == FALSE & width == "auto") {
            width <- cellwidth * ncol(remove_neg) / 72 + 2.5
        } else if ((is.na(cellwidth) | legend.title == TRUE) & width == "auto") {
            width <- 10
        }
        
        if (!is.na(cellheight) & legend.title == FALSE & height == "auto") {
            height <- cellheight * ncol(remove_neg) / 72 + 1.5
        } else if ((is.na(cellheight) | legend.title == TRUE) & height == "auto") {
            height <- 8
        }

        pdf(file, width = width, height = height)
        print(p)
        dev.off()
        png(img_file, width = width, height = height, unit = "in", res = png_res)
        print(p)
        dev.off()
        

        # pheatmap::pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F, 
        #                    main = "Perturbations Correlation",
        #                    c(colorRampPalette(colors = c(color[1], color[2]))(length(bk)/2),
        #                      colorRampPalette(colors = c(color[2], color[3]))(length(bk)/2)),
        #                    legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
        #                    breaks = bk, show_rownames = T, show_colnames = T, fontsize = fontsize, 
        #                    cellwidth = cellwidth, cellheight = cellheight, filename = file, silent = TRUE, 
        #                    angle_col = angle)
        
        #                  display_numbers = matrix(ifelse(sg_cor$p <= 0.01, "**", ifelse(sg_cor$p <= 0.05 , "*", " ")), nrow(sg_cor$p))
        # pheatmap::pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F, 
        #                    main = "Perturbations Correlation",
        #                    c(colorRampPalette(colors = c(color[1], color[2]))(length(bk)/2),
        #                      colorRampPalette(colors = c(color[2], color[3]))(length(bk)/2)),
        #                    legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
        #                    breaks = bk, show_rownames = T, show_colnames = T, fontsize = fontsize, 
        #                    cellwidth = cellwidth, cellheight = cellheight, filename = img_file, silent = TRUE, 
        #                    angle_col = angle)
        
        return(p)
    } else {
        return(p)
    }
}