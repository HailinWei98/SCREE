#' function definitions ##### Plot for sgRNA information.
#' @import stringr
#' @import plyr
#' @import ggplot2
#' @export

sgRNA_quality_plot<- function(sg_dir, mtx_dir, LABEL = "", prefix = "./"){
    
    #read files
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        mtx <- readRDS(mtx_dir)
    } else {
        mtx <- mtx_dir
    }
    if (is.character(sg_dir)) {
        message(paste("Reading sgRNA lib file:", sg_dir))
        sg_lib <- read.table(sg_dir, header = T)
    } else {
        sg_lib <- sg_dir
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
    labs(x = "Rank", y = "Cell Numbers", title = "Cell Numbers of ORF") +
    theme(plot.title = element_text(hjust = 0.5, size = 25), 
          legend.text = element_text(size = 12),
          text = element_text(hjust = 0.5, face = "bold"),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) +
    guides(colour = guide_legend(title.hjust = 0.5))
    
    g1 <- suppressMessages(g1 + ggsci::scale_colour_npg(name = "Top 10 ORF", breaks = head(sg_count, 10)$sgRNA))
    
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
    geom_bar(stat = "identity",color = "black", fill = "#8491B4FF") + theme_test() + 
    labs(x = "ORF Numbers", y = "Cell Numbers", title = "ORF Numbers in Each Cell") +
    theme(plot.title = element_text(hjust = 0.5, size = 25),
          text = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 16),
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) +
    geom_text(aes(label = freq), stat = "identity", vjust = -0.5, size = 6)
    
    #save plot
    
    dir <- file.path(prefix, "pdf")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    dir <- file.path(dir, "sgRNA")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "sgRNA")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }

    pdf(file = file.path(dir, paste(LABEL, "sgRNA_quality.pdf", sep = "")), height = 8, width = 8)
    print(g1)
    print(g2)
    dev.off()
    
    png(file.path(img_dir, paste(LABEL, "cell_numbers.png", sep = "")), 
        width = 600, height = 600)
    print(g1)
    dev.off()
    
    png(file.path(img_dir, paste(LABEL, "sgRNA_numbers.png", sep = "")), 
        width = 600, height = 600)
    print(g2)
    dev.off()

    dir2 <- file.path(dir, "sgRNA_quality_of_gene")
    new_dir <- file.path(img_dir, "sgRNA_quality_of_gene")
    if (!(dir.exists(new_dir))) {
        dir.create(path = new_dir)
    }
    if (!(dir.exists(dir2))) {
        dir.create(path = dir2)
    }
    
    for(select_gene in unique(sg_count$gene)){
        pdf(file = file.path(dir2, paste(select_gene, ".pdf", sep = "")), height = 8, width = 8)
        p1 <- ggplot(data = sg_count, mapping = aes(x = order, y = freq)) +
        geom_line(size = 1) + theme_test() +
        geom_point(data = subset(sg_count, gene == select_gene),
                   mapping = aes(x = order, y = freq, color = sgRNA), size = 5)+
        labs(x = "Rank", y = "Cell Numbers", title = "Cell Numbers of sgRNA") +
        theme(plot.title = element_text(hjust = 0.5, size = 25), 
              legend.text = element_text(size = 12),
              text = element_text(hjust = 0.5, face = "bold"),
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), 
              legend.title = element_text(size = 18),
              axis.title.x = element_text(size = 20), 
              axis.title.y = element_text(size = 20),
              axis.text.y = element_text(hjust = 0.5, size = 16)) + 
        guides(colour = guide_legend(title.hjust = 0.5))
        p1 <- suppressMessages(p1 + ggsci::scale_colour_npg(name = stringr::str_wrap(paste("", select_gene, sep = ""),
                                                                                     width = 30),
                                                            breaks = head(subset(sg_count, 
                                                                                 gene == select_gene), 10)$sgRNA))
        print(p1)
        dev.off()
        
        png(file.path(new_dir, paste(select_gene, ".png", sep = "")), 
            width = 600, height = 600)
        print(p1)
        dev.off()
    }
    return(list(g1, g2))
}
