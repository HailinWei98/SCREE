#' Plot DE Gene Numbers for Each Perturbation
#'
#' Plot potential DE gene numbers for each perturbation based on score and p-value from \code{\link[SCREE]{imporved_scmageck_lr}}.
#'
#' @param score Data frame or directory of score from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param pval Data frame or directory of p_value from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param project Title of the barplot. Default is "perturb".
#' @param top The top perturbations to plot, with the most differential gene numbers. Default is 20, means plot top 20 perturbations.
#' @param select Specific perturbations selected to plot. Default is \code{NULL}, which means to use the top perturbations. If provided a vector of perturbations, \code{top} parameter will be ignored.
#' @param score_cut Score cutoff of \code{improved_scmageck_lr} results. Default is 0.2.
#' @param pval_cut P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.
#' @param y_break Break of the barplot, the first item means the maximum y-axis of the bottom ggplot while the second item means the minimum y-axis of the top ggplot. If you don't need a break of y-axis, set the second item to a number bigger than the max number of the up-regulated and down-regulated genes, among all perturbations to plot. Default is c(50, 200).
#' @param y_height Heights of the two breaked ggplot to be arranged. Default is c(1/5, 4/5).
#' @param title.size Numeric, title size of the barplot. Default is 25.
#' @param legend.text.size Numeric, legend text size of the barplot. Default is 16.
#' @param legend.title.size Numeric, legend title size of the barplot. Default is 20.
#' @param x.text.size Numeric, x-axis text size of the barplot. Default is 16.
#' @param x.title.size Numeric, x-axis title size of the barplot. Default is 20.
#' @param y.text.size Numeric, y-axis text size of the barplot. Default is 16.
#' @param y.title.size Numeric, y-axis title size of the barplot. Default is 20.
#' @param sort_by How to order the perturbations labeled in x-axis. Can be one of "name" (name of the perturbations), "select" (the same order as the vector of the \code{select} parameter) and "number" (total differential gene numbers). 
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @import ggsci
#' @import ggpubr
#' @importFrom plyr count
#' @import reshape2
#' @import ggplot2
#' @importFrom utils read.table write.table
#' 
#' @export

DE_gene_plot<- function(score, pval, project = "perturb", top = 20, select = NULL, score_cut = 0.2, pval_cut = 0.05, y_break = c(50, 200), y_height = c(1/5, 4/5), title.size = 25, legend.text.size = 16, legend.title.size = 20, x.text.size = 16, x.title.size = 20, y.text.size = 16, y.title.size = 20, sort_by = "name", plot.save = TRUE, prefix = ".", label = "", width = 7, height = 7, png_res = 720){
    
    #get score and pval
    
    if (is.character(score)) {
        score <- read.table(score, header = T, row.names = 1)
    }
    
    if (is.character(pval)) {
        pval <- read.table(pval, header = T, row.names = 1)
    }
    
    #get DE gene numbers(up, down, non)
    
    de_genes <- data.frame(non = rep(0, ncol(score)),
                           up = rep(0,ncol(score)), 
                           down = rep(0, ncol(score)))
    rownames(de_genes) <- colnames(score)
    for (i in colnames(score)) {
        scmageck <- data.frame(score = score[, i], pval = pval[, i])
        scmageck$Difference <- "non"
        scmageck$Difference[scmageck$score > score_cut & scmageck$pval < pval_cut] <- "up"
        scmageck$Difference[scmageck$score < -score_cut & scmageck$pval < pval_cut] <- "down"
        a <- plyr::count(scmageck$Difference)
        rownames(a) <- a$x
        for (j in a$x) {
            de_genes[i, j] <- a[j, 2]
        }
    }
    
    new <- data.frame(reshape2::melt(de_genes[, 2:3]))
    new$factor <- rep(rownames(de_genes), 2)
    colnames(new) <- c("Difference", "gene_numbers", "factor")
    new1 <- new
    new1$Difference <- as.character(new1$Difference)
    
    #get total DE gene numbers of all perturbations
    
    for(i in unique(new$factor)){
        a <- subset(new, factor == i)
        all <- sum(a$gene_numbers)
        new1 <- rbind(new1, c("all", all, i))
    }
    
    new2 <- subset(new1, Difference == "all")
    
    if (is.null(select)) {
        if (is.null(top)) {
            top <- ncol(score)
        }
        select <- colnames(score)
        num_level <- head(new2[order(as.numeric(new2$gene_numbers), decreasing = T), ], top)$factor
        new3 <- subset(new, factor %in% num_level)
    } else if (is.character(select)) {
        new2 <- subset(new2, factor %in% select)
        num_level <- new2[order(as.numeric(new2$gene_numbers), decreasing = T), ]$factor
        new3 <- subset(new, factor %in% select)
    } else if (is.numeric(select)) {
        select <- colnames(score[, select])
        new2 <- subset(new2, factor %in% select)
        num_level <- new2[order(as.numeric(new2$gene_numbers), decreasing = T), ]$factor
        new3 <- subset(new, factor %in% select)
    } else {
        stop("Please select in correct format!")
    }
    
    if (sort_by != "name") {
        if (sort_by == "select") {
            new3$factor <- factor(new3$factor, levels = select)
        } else if (sort_by == "number") {
            new3$factor <- factor(new3$factor, levels = num_level)
        } else {
            stop("Please set sort_by from c('name', 'select', 'number').")
        }
    }
#    new3[which(new3$Difference == "down"), ]$gene_numbers <- -new3[which(new3$Difference == "down"), ]$gene_numbers
    
    #get range of y axis
    
    y_max <- max(new3$gene_numbers)
    if (ceiling(y_max/100) == 0) {
        y_max <- y_max
    } else {
        y_max <- ceiling(y_max/100) * 100
    }
    
#     if(ylimit == "auto"){
#         y_max <- max(abs(new3$gene_numbers))
#         if(ceiling(y_max/100) == 0){
#             y_max <- 100
#         }else{
#             y_max <- ceiling(y_max/100)*100
#         }
#         ylimit <- c(-y_max, y_max, y_max/2)
        
#     }else if(is.vector(ylimit)){
#         ylimit <- ylimit
#     }else{
#         warning("Please input correct ylimit, use c(-600, 600, 200) instead.")
#         ylimit = c(-600, 600, 200)
#     }
    
    #plot
    
    if (y_max <= y_break[2]) {
        p1 <- ggplot(new3, aes(x = factor, y = gene_numbers)) +
        geom_bar(stat = 'identity', aes(fill = Difference), position = position_dodge(0.9), colour = "black") +
        theme_test() +
        labs(x = "Perturbation", y = "Gene Numbers", title = project) +
        theme(plot.title = element_text(hjust = 0.5, size = title.size),
              text = element_text(hjust = 0.5, face = "bold"),
              legend.text = element_text(size = legend.text.size),
              legend.title = element_text(size = legend.title.size),
              axis.text.x = element_text(angle = 90, 
                                         hjust = 0.5, 
                                         vjust = 0.5, 
                                         size = x.text.size, 
                                         color = "black"),
              axis.title.x = element_text(size = x.title.size), 
              axis.title.y = element_text(size = y.title.size),
              axis.text.y = element_text(hjust = 0.5, size = y.text.size, color = "black")) + 
        ggsci::scale_fill_npg()
    } else {
        p3 <- ggplot(new3, aes(x = factor, y = gene_numbers)) +
        geom_bar(stat = 'identity', aes(fill = Difference), position = position_dodge(0.9), colour = "black") +
        theme_test() +
        labs(x = "Perturbation", y = "Gene Numbers") +
        theme(text = element_text(hjust = 0.5, face = "bold"),
              legend.text = element_text(size = legend.text.size),
              axis.text.x = element_text(angle = 90, 
                                         hjust = 0.5, 
                                         vjust = 0.5, 
                                         size = x.text.size, 
                                         color = "black"),
              axis.title.x = element_text(size = x.title.size), 
              axis.title.y = element_text(size = y.title.size),
              axis.text.y = element_text(hjust = 0.5, size = y.text.size, color = "black")) + 
        coord_cartesian(ylim = c(0, y_break[1])) + 
        scale_y_continuous(breaks = c(0, y_break[1])) + 
        ggsci::scale_fill_npg()
        p4 <- ggplot(new3, aes(x = factor, y = gene_numbers)) +
        geom_bar(stat = 'identity', aes(fill = Difference), position = position_dodge(0.9), colour = "black") +
        theme_test() +
        labs(x = NULL, y = NULL, title = project) +
        theme(plot.title = element_text(hjust = 0.5, size = title.size),
              text = element_text(hjust = 0.5, face = "bold"),
              legend.text = element_text(size = legend.text.size),
              legend.title = element_text(size = legend.title.size),
              axis.text.y = element_text(hjust = 0.5, size = y.text.size, color = "black"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) + 
        coord_cartesian(ylim = c(y_break[2], y_max)) + 
        scale_y_continuous(breaks = c(y_break[2], y_max)) + 
        ggsci::scale_fill_npg()
        p1 <- ggpubr::ggarrange(p4, p3, heights = y_height, ncol = 1, nrow = 2, 
                                common.legend = TRUE, legend="right", align = "v") 
    }
    
    #save plot
    
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

        pdf(file.path(pdf_dir,
                      paste(label, "DE_gene_cutoff", score_cut, "_p", pval_cut, ".pdf", sep = "")), 
            width = width, height = height)
        print(p1)
        dev.off()

        png(file.path(img_dir, 
                      paste(label, "DE_gene_cutoff", score_cut, "_p", pval_cut, ".png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(p1)
        dev.off()
    }
    
    return(p1)
}
