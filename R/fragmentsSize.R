#' Plot Distribution of Fragments Size
#'
#' Calculate fragments size and visualize the distribution.
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param fragments Data frame or directory of fragments file.
#' @param CBCindex Cell barcode index, indicates the column index of cell barcode. Default is 4, in line with a standard fragments file.
#' @param startIndex
#' @param endIndex
#' @param maxSize
#' @param title.size Numeric, title size of the plot. Default is 25.
#' @param x.text.size Numeric, x-axis text size of the plot. Default is 16.
#' @param x.title.size Numeric, x-axis title size of the plot. Default is 20.
#' @param y.text.size Numeric, y-axis text size of the plot. Default is 16.
#' @param y.title.size Numeric, y-axis title size of the plot. Default is 20.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#' 
#' @import ggplot2
#' @import plyr
#' @export

fragmentsSize <- function(mtx, fragments, CBCindex = 4, startIndex = 2, endIndex = 3, maxSize = 1000, title.size = 25, x.text.size = 16, x.title.size = 20, y.text.size = 16, y.title.size = 20, plot.save = TRUE, prefix = ".", label = "", width = 7, height = 7, png_res = 720){
    
    #get fragments information
    
    if (is.character(fragments)) {
        message(paste("Reading fragments file:", fragments))
        fragments <- read.table(gzfile(fragments), header = F)
    } else {
        fragments <- fragments
    }
    
    #get peak matrix
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        peak <- readRDS(mtx)
    } else {
        peak <- mtx
    }

    #filter cells not in peak matrix
    
    fragments <- fragments[(fragments[, CBCindex] %in% colnames(peak)), ]
    
    #calculate fragments size
    
    fragments$size <- fragments[, 3] - fragments[, 2]
    
    #plot
    
    fragments_count <- plyr::count(fragments$size)
    fragments_count <- subset(fragments_count, x <= maxSize)
    g1 <- ggplot(data = fragments_count, mapping = aes(x = x, y = freq)) +
    geom_line(size = 1, color = "#E64B35FF") + 
    theme_test() +
    labs(x = "Fragments Size", y = "Frequency",title = "Fragments Size Distribution") +
    theme(plot.title = element_text(hjust = 0.5, size = title.size), 
          text = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(hjust = 0.5, size = y.text.size),
          axis.text.x = element_text(size = x.text.size), 
          axis.title.y = element_text(size = y.title.size),
          axis.title.x = element_text(size = x.title.size))
    
    if (plot.save == TRUE) {
        
        #save plot

        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }

        dir <- file.path(dir, "ATAC_quality")
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


        pdf(file = file.path(dir, paste(label, "FragmentsSize.pdf", sep = "")), 
            width = width, height = height)
        print(g1)
        dev.off()

        png(file.path(img_dir, paste(label, "FragmentsSize.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(g1)
        dev.off()
    
    }
    
    return(g1)
}