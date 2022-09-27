#' Volcano Plot For Potential DEGs
#'
#' Visualize regulatory score and p-value distribution for each perturbation.
#'
#' @param score Data frame or directory of score from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param pval Data frame or directory of p_value from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param selected Perturbation to visualize. Default is \code{NULL}, all perturbations will be chosen.
#' @param score_cut Score cutoff of \code{improved_scmageck_lr} results. Default is 0.2.
#' @param pval_cut P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.
#' @param color Color for volcano plot. Default is c("#4DBBD5FF", "grey", "#E64B35FF"), corresponding to "down-regulated", "no sense" and "up-regulated".
#' @param showCategory Numbers of genes to be labeled for both "down" and "up". Default is 10, and will label 10 "up-regulated" genes and 10 "down-regulated" genes.
#' @param title.size Numeric, title size of the volcano plot. Default is 23.
#' @param legend.text.size Numeric, legend text size of the volcano plot. Default is 16.
#' @param legend.title.size Numeric, legend title size of the volcano plot. Default is 20.
#' @param x.text.size Numeric, x-axis text size of the volcano plot. Default is 16.
#' @param x.title.size Numeric, x-axis title size of the volcano plot. Default is 20.
#' @param y.text.size Numeric, y-axis text size of the volcano plot. Default is 16.
#' @param y.title.size Numeric, y-axis title size of the volcano plot. Default is 20.
#' @param pt.size Point size of the volcano plot. Default is 2.
#' @param raster Logical, convert points to raster format, will be useful to reduce the storage cost of the output figure or pdf. Default is \code{FALSE}.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom utils read.table write.table
#' @import ggplot2
#' @import stringr
#' @import ggrepel
#' @import ggrastr
#' @export

volcano <- function(score, pval, selected = NULL, pval_cut = 0.05, score_cut = 0.2, color = c("#4DBBD5FF","grey","#E64B35FF"), showCategory = 10, title.size = 23, legend.text.size = 16, legend.title.size = 20, x.text.size = 16, x.title.size = 20, y.text.size = 16, y.title.size = 20, pt.size = 2, raster = FALSE, plot.save = TRUE, prefix = ".", label = "", width = 7, height = 7, png_res = 720){
    
    #get score and p_val
    
    if (is.character(score)) {
        score <- read.table(score, header = T, row.names = 1)
    } else {
        score <- score
    }
    
    if (is.character(pval)) {
        pval <- read.table(pval, header = T, row.names = 1)
    }
    
    #create save path
    
    results <- list()
    j = 0
    
    if (plot.save == TRUE) {
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }

        dir <- file.path(dir, "perturbation_efficiency")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }

        dir <- file.path(dir, "volcano")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }
        
        pdf_dir <- file.path(dir, "pdf")
        if (!(dir.exists(pdf_dir))) {
            dir.create(path = pdf_dir)
        }

        img_dir <- file.path(prefix, "img")
        if (!(dir.exists(img_dir))) {
            dir.create(path = img_dir)
        }
        
    }
    
    #prepare color for plot
    
    diff <- color
    names(diff) <- c("down", "non", "up")

    #plot for each perturbation
    
    if (is.null(selected)) {
        selected <- colnames(score)
    }
    
    for (i in selected) {
        
        #classify genes for each perturbation
        
        scmageck <- data.frame(score = score[, i], p_val = pval[, i])
        rownames(scmageck) <- rownames(score)
        scmageck$diff <- "non"
        scmageck$diff[scmageck$score > score_cut & scmageck$p_val < pval_cut] <- "up"
        scmageck$diff[scmageck$score < -score_cut & scmageck$p_val < pval_cut] <- "down"
        a <- plyr::count(scmageck$diff)
        up_gene <- head(scmageck[order(-scmageck$score), ], showCategory)
        up_gene <- subset(up_gene, diff != "non")
        down_gene <- head(scmageck[order(scmageck$score), ], showCategory)
        down_gene <- subset(down_gene, diff != "non")
        top_gene<- unique(c(rownames(up_gene), rownames(down_gene)))
        scmageck$lab <- ""
        scmageck[top_gene, "lab"] <- top_gene
        score_max <- max(abs(scmageck$score))
        
        #plot
        
        if(stringr::str_starts(i, pattern = "chr")) {
            title = paste(i, "Potential Targets", sep = "\n")
        } else {
            title = paste(i, "Potential Targets", sep = " ")
        }
        
        if (raster == TRUE) {
            
            p1 <- ggplot(data = scmageck, 
                         mapping = aes(x = score, y = -log10(p_val + 10e-6), colour = diff, label = lab)) + 
            geom_point_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300)) + 
            theme_test() + 
            labs(title = title) + 
            theme(plot.title = element_text(hjust = 0.5, size = title.size), 
                  text = element_text(hjust = 0.5, face = "bold"),
                  legend.text = element_text(size = legend.text.size),
                  legend.title = element_text(size = legend.title.size),
                  axis.title.x = element_text(size = x.title.size), 
                  axis.title.y = element_text(size = y.title.size),
                  axis.text.y = element_text(hjust = 0.5, size = y.text.size),
                  axis.text.x = element_text(size = x.text.size)) +  
            scale_color_manual(values = diff[sort(unique(scmageck$diff))]) + 
            geom_vline(xintercept = c(-score_cut, score_cut), linetype = "dotted") + 
            geom_label_repel(max.overlaps = 500, show.legend = FALSE, force = T) + 
            guides(colour = guide_legend(title = "Difference", title.hjust = 0.5))

        } else {
            p1 <- ggplot(data = scmageck, 
                         mapping = aes(x = score, y = -log10(p_val + 10e-6), colour = diff, label = lab)) + 
            geom_point(size = pt.size) + 
            theme_test() + 
            labs(title = title) + 
            theme(plot.title = element_text(hjust = 0.5, size = title.size), 
                  text = element_text(hjust = 0.5, face = "bold"),
                  legend.text = element_text(size = legend.text.size),
                  legend.title = element_text(size = legend.title.size),
                  axis.title.x = element_text(size = x.title.size), 
                  axis.title.y = element_text(size = y.title.size),
                  axis.text.y = element_text(hjust = 0.5, size = y.text.size),
                  axis.text.x = element_text(size = x.text.size)) +  
            scale_color_manual(values = diff[sort(unique(scmageck$diff))]) + 
            geom_vline(xintercept = c(-score_cut, score_cut), linetype = "dotted") + 
            geom_label_repel(max.overlaps = 500, show.legend = FALSE, force = T) + 
            guides(colour = guide_legend(title = "Difference", title.hjust = 0.5))
        
        }

        if (score_max < score_cut) {
            p1 <- p1 + xlim(-score_max, score_max)
        }
        #save results and return
        
        j = j + 1
        results[[j]] <- p1
        names(results)[j] <- i
        if (plot.save == TRUE) {
            pdf(file = file.path(dir, paste(i, ".pdf", sep = "")), width = width, height = height)
            print(p1)
            dev.off()

            png(file.path(img_dir, paste(i, ".png", sep = "")), 
                width = width, height = height, unit = "in", res = png_res)
            print(p1)
            dev.off()
        }
    }
    
    return(results) 
}