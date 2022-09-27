#' Calculate Perturbation Enrichment
#'
#' Calculate enrichment ratio for each perturbation in every cluster.
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows. Note that the dataset has to be normalized and scaled, and need a sgRNA information column named "perturbations" in the meta data.
#' @param sg_lib Data frame or directory to a txt file containing 3 columns: cell, barcode, gene. If sgRNA information stored in a matrix-like format or input data frame only has sgRNA frequency of each cell, use \code{\link[SCREE]{sgRNAassign}} to assign sgRNA to each cell.
#' @param NTC The name of negative controls. Default is "NTC".
#' @param NTC.cal Logical, calculate negative control or not. Default is \code{TRUE}, which will remove cells assigned with negative control from the SeuratObject.
#' @param table.save Logical, save perturbation enrichment table or not. Default is \code{TRUE}.
#' @param top Top perturbations with the highest enrichment ratio to visualize. Default is \code{NULL}, using all perturbations to visualize.
#' @param range Range of the ratio. Since for some datasets, the max enrichment ratio among all perturbations is extremely small, the range can be set to a smaller limitation to enhance the enrichment. Default is c(0, 1).
#' @param color Color for heatmap. Default is c("white", "coral1"), corresponding to the score from 0 to 1 or other ranges.
#' @param cell Cell width and height settings for heatmap. If set this parameter to a single number or \code{NA}, cell width and height will be the same, and \code{NA} means the cell will be adaptive to the input perturbation numbers. If set this parameter to a vector including two numbers or \code{NA}, the first one will be cell width and the second one will be height. 
#' @param fontsize Fontsize of the heatmap. Default is 15.
#' @param angle Angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315). Default is 90.
#' @param legend.title Include the title of legend in the heatmap or not. Default is FALSE.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 10.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 8.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom utils read.table write.table
#' @import pheatmap
#' @import ggplotify
#' @import grid
#' @import gtable
#' @export

CalculatePerturbEnrichment <- function(mtx, sg_lib, NTC = "NTC", NTC.cal = TRUE, table.save = TRUE, top = NULL, range = c(0, 1), color = c("white", "coral1"), cell = c(15, 20), fontsize = 15, angle = 90, legend.title = FALSE, plot.save = TRUE, prefix= ".", label = "", width = 10, height = 8, png_res = 720){
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        eccite <- readRDS(mtx)
    } else {
        eccite <- mtx
    }
    
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        sg_lib <- read.table(sg_lib, header = T)
    }
    
    #perturbation information
    
    if (!("perturbations" %in% colnames(eccite@meta.data))) {
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
    
    cluster_pro <- t(pro_mat) / as.vector(table(G_mat$cluster))
    
    bk <- seq(range[1], range[2], by = (range[2] - range[1]) / 20)
    lg_bk <- seq(range[1], range[2], by = (range[2] - range[1]) / 2)
#     if (nrow(pro_mat)/ncol(pro_mat) >= 10) {
#         if (max(cluster_pro) < 0.2) {
#             bk <- seq(0, 0.2, by = 0.01)
#             lg_bk <- seq(0, 0.2, 0.1)
#         }
#     }

    
    cluster_pro_plot <- cluster_pro
    
    if (length(cell) == 1) {
        
        if (!is.numeric(cell) & !is.na(cell)) {
            warning("Cell must be a numeric parameter or NA, set cell to c(15, 20). ")
            cell <- c(15, 20)
        } else {
            cellwidth = cell
            cellheight = cell
        }
  
    } else if (length(cell) == 2) {
        if ((!is.numeric(cell[1]) & !is.na(cell[1])) | (!is.numeric(cell[2]) & !is.na(cell[2]))) {
            warning("Cell must be a numeric parameter or NA, set cell to c(15, 20). ")
            cell <- c(15, 20)
        }
        cellwidth = cell[1]
        cellheight = cell[2]
    }                 
    
    if (is.null(top)) {
        p1 <- pheatmap::pheatmap(cluster_pro_plot, treeheight_col = F, treeheight_row = F, 
                                 main = "Perturbations Enrichment",
                                 colorRampPalette(colors = c(color[1], color[2]))(length(bk)),
                                 legend_breaks = lg_bk, breaks = bk, show_rownames = T, fontsize = fontsize, 
                                 cellwidth = cellwidth, cellheight = cellheight,
                                 show_colnames = T, silent = T, angle_col = angle)
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
                                 colorRampPalette(colors = c(color[1], color[2]))(length(bk)),
                                 legend_breaks = lg_bk, breaks = bk, show_rownames = T, fontsize = fontsize, 
                                 cellwidth = cellwidth, cellheight = cellheight,
                                 show_colnames = T, silent = T, angle_col = angle)
    }
    
    if (table.save == TRUE) {
        write.table(cluster_pro, file = file.path(prefix, paste(label, "perturb_ratio.txt", sep = "")), quote = FALSE)
    }
    
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
        
        img_dir <- file.path(prefix, "img")
        if (!(dir.exists(img_dir))) {
            dir.create(img_dir)
        }
        
        file <- file.path(dir, paste(label, "perturb_ratio.pdf", sep = ""))
        img_file <- file.path(img_dir, paste(label, "perturb_ratio.png", sep = ""))
        # if (ncol(cluster_pro_plot) > 20) {
        #     if (nrow(cluster_pro_plot) > 20) {
        #         cellheight = 15
        #         cellwidth = 15
        #     } else {
        #         cellheight = 20
        #         cellwidth = 15
        #     }
        #     size = 15
        # } else {
        #     cellheight = NA
        #     cellwidth = NA
        #     size = 10
        # }
#         pheatmap(cluster_pro_plot, treeheight_col = F, treeheight_row = F, main = "Perturbations Enrichment",
#                  colorRampPalette(colors = c("white", "#E64B35FF"))(length(bk)),
#                  legend_breaks = lg_bk, breaks = bk, show_rownames = T, 
#                  show_colnames = T, fontsize = size, 
#                  cellwidth = cellwidth, cellheight = cellheight, filename = file)
        
        
        # pheatmap(cluster_pro_plot, treeheight_col = F, treeheight_row = F, main = "Perturbations Enrichment",
        #          colorRampPalette(colors = c("white", "#E64B35FF"))(length(bk)),
        #          legend_breaks = lg_bk, breaks = bk, show_rownames = T, 
        #          show_colnames = T, fontsize = size, 
        #          cellwidth = cellwidth, cellheight = cellheight, filename = img_file)
        
        pdf(file, width = width, height = height)
        print(p)
        dev.off()
        png(img_file, width = width, height = height, unit = "in", res = png_res)
        print(p)
        dev.off()
    }
    
    return(p1)
    
}