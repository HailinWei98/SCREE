#' sgRNA Information Visualization
#'
#' Visualize cell numbers for each sgRNA and sgRNA numbers in each cell.
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows. 
#' @param sg_lib Data frame or directory to a txt file containing 3 columns: cell, barcode, gene. If sgRNA information stored in a matrix-like format or input data frame only has sgRNA frequency of each cell, use \code{\link[SCREE]{sgRNAassign}} to assign sgRNA to each cell.
#' @param title.size Numeric, title size of the UMAP. Default is 25.
#' @param legend.text.size Numeric, legend text size of the plot. Default is 12.
#' @param legend.title.size Numeric, legend title size of the plot. Default is 18.
#' @param x.text.size Numeric, x-axis text size of the plot. Default is 16.
#' @param x.title.size Numeric, x-axis title size of the plot. Default is 20.
#' @param y.text.size Numeric, y-axis text size of the plot. Default is 16.
#' @param y.title.size Numeric, y-axis title size of the plot. Default is 20.
#' @param label.size Label size of the the barplot.
#' @param bar_width Bar width. Default is \code{NULL}, set to 90% of the resolution of the data.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 8.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 8.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom utils read.table write.table
#' @import stringr
#' @importFrom plyr count
#' @import ggplot2
#' @import ggsci
#' @export

sgRNA_quality_plot <- function(mtx, sg_lib, title.size = 25, legend.text.size = 12, legend.title.size = 18, x.text.size = 16, x.title.size = 20, y.text.size = 16, y.title.size = 20, label.size = 6, bar_width = NULL, plot.save = TRUE, prefix = ".", label = "", width = 8, height = 8, png_res = 720){
    
    #read files
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        mtx <- readRDS(mtx)
    } 
    
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        sg_lib <- read.table(sg_lib, header = T)
    }

    #remove cells in sgRNA library that are not included in matrix
    
    sg_lib_filtered <- subset(sg_lib, cell %in% intersect(sg_lib$cell, colnames(mtx)))

    #count sgRNA in each cell;count cells of each sgRNA
    
    sg_count <- plyr::count(subset(sg_lib_filtered, cell %in% colnames(mtx))$barcode)
    colnames(sg_count) <- c("sgRNA", "freq")
    sg_count <- sg_count[order(-sg_count$freq), ]
    sg_count$order <- seq(1, nrow(sg_count))

    #prepare for plot each gene
    
    sg_gene <- unique(sg_lib[ ,c("barcode", "gene")])
    rownames(sg_gene) <- sg_gene$barcode
    sg_count$gene <- sg_gene[sg_count$sgRNA, "gene"]

    g1 <- ggplot(data = sg_count, mapping = aes(x = order, y = freq)) +
    geom_line(size = 1) + theme_test() +
    geom_point(data = head(sg_count, 10),
               mapping = aes(x = order, y = freq, color = sgRNA), size = 5) +
    labs(x = "Rank", y = "Cell Numbers", title = "Cell Numbers of sgRNA") +
    theme(plot.title = element_text(hjust = 0.5, size = title.size), 
          legend.text = element_text(size = legend.text.size),
          text = element_text(hjust = 0.5, face = "bold"),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          legend.title = element_text(size = legend.title.size),
          axis.title.x = element_text(size = x.title.size), 
          axis.title.y = element_text(size = y.title.size),
          axis.text.y = element_text(hjust = 0.5, size = y.text.size)) +
    guides(colour = guide_legend(title.hjust = 0.5))
    
    g1 <- suppressMessages(g1 + ggsci::scale_colour_npg(name = "Top 10 sgRNA", breaks = head(sg_count, 10)$sgRNA))
    
    sg_num_count <- plyr::count(mtx[["sgRNA_num"]])
    if (max(mtx$sgRNA_num) > 10) {
        over_10 <- subset(sg_num_count, sgRNA_num > 10)
        sg_num_count <- subset(sg_num_count, sgRNA_num <= 10)
        over_10 <- sum(over_10$freq)
        sg_num_count <- rbind(sg_num_count, c(">10", over_10))
    }
    colnames(sg_num_count) <- c("sgRNA_num", "freq")
    sg_num_count$freq <- as.numeric(sg_num_count$freq)
    sg_num_count$sgRNA_num <- factor(sg_num_count$sgRNA_num, levels = sg_num_count$sgRNA_num, ordered = T)

    g2 <- ggplot(sg_num_count,mapping = aes(x = sgRNA_num, y = freq)) +
    geom_bar(stat = "identity",color = "black", fill = "#8491B4FF", width = bar_width) + theme_test() + 
    labs(x = "sgRNA Numbers", y = "Cell Numbers", title = "sgRNA Numbers in Each Cell") +
    theme(plot.title = element_text(hjust = 0.5, size = title.size),
          text = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = x.text.size),
          axis.title.x = element_text(size = x.title.size), 
          axis.title.y = element_text(size = y.title.size),
          axis.text.y = element_text(hjust = 0.5, size = y.text.size)) +
    geom_text(aes(label = freq), stat = "identity", vjust = -0.5, size = label.size)
    
    #save plot
    
    if (plot.save == TRUE) {
        
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }

        dir <- file.path(dir, "sgRNA")
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

        pdf(file = file.path(pdf_dir, paste(label, "sgRNA_quality.pdf", sep = "")), height = height, width = width)
        print(g1)
        print(g2)
        dev.off()

        png(file.path(img_dir, paste(label, "cell_numbers.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(g1)
        dev.off()

        png(file.path(img_dir, paste(label, "sgRNA_numbers.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(g2)
        dev.off()

        dir2 <- file.path(pdf_dir, paste(label, "sgRNA_quality_of_gene", sep = ""))
        new_dir <- file.path(img_dir, paste(label, "sgRNA_quality_of_gene", sep = ""))
        if (!(dir.exists(new_dir))) {
            dir.create(path = new_dir)
        }
        if (!(dir.exists(dir2))) {
            dir.create(path = dir2)
        }

        for(select_gene in unique(sg_count$gene)){
            pdf(file = file.path(dir2, paste(select_gene, ".pdf", sep = "")), height = height, width = width)
            p1 <- ggplot(data = sg_count, mapping = aes(x = order, y = freq)) +
            geom_line(size = 1) + theme_test() +
            geom_point(data = subset(sg_count, gene == select_gene),
                       mapping = aes(x = order, y = freq, color = sgRNA), size = 5)+
            labs(x = "Rank", y = "Cell Numbers", title = "Cell Numbers of sgRNA", fill = "sgRNA") +
            theme(plot.title = element_text(hjust = 0.5, size = title.size), 
                  legend.text = element_text(size = legend.text.size),
                  text = element_text(hjust = 0.5, face = "bold"),
                  axis.ticks.x = element_blank(), 
                  axis.text.x = element_blank(), 
                  legend.title = element_text(size = legend.title.size),
                  axis.title.x = element_text(size = x.title.size), 
                  axis.title.y = element_text(size = y.title.size),
                  axis.text.y = element_text(hjust = 0.5, size = y.text.size)) + 
            guides(colour = guide_legend(title.hjust = 0.5))
            p1 <- suppressMessages(p1 + ggsci::scale_colour_npg(name = stringr::str_wrap(paste("", select_gene, sep = ""),
                                                                                         width = 30),
                                                                breaks = head(subset(sg_count, 
                                                                                     gene == select_gene), 10)$sgRNA))
            print(p1)
            dev.off()

            png(file.path(new_dir, paste(select_gene, ".png", sep = "")), 
                width = width, height = height, unit = "in", res = png_res)
            print(p1)
            dev.off()
        }
    }

    return(list(g1, g2))
}
