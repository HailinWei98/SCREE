#' function definitions ##### Plot DE gene numbers of each perturbations
#' @export

DE_gene_plot<- function(score_dir, pval_dir, project = "perturb", top = 20, 
                        prefix = "./", label = "", p_val_cut = 0.05, score_cut = 0.5, 
                        y_break = c(50, 200),  height = c(1/5, 4/5)){
    #get score and p_val
    
    if(is.character(score_dir)){
        score <- read.table(score_dir, header = T, row.names = 1)
    }else{
        score <- score_dir
    }
    if(is.character(pval_dir)){
        p_val <- read.table(pval_dir, header = T, row.names = 1)
    }else{
        p_val <- pval_dir
    }
    
    #get DE gene numbers(up, down, non)
    
    de_genes <- data.frame(non = rep(0, ncol(score)),
                           up = rep(0,ncol(score)), down = rep(0, ncol(score)))
    rownames(de_genes) <- colnames(score)
    for(i in colnames(score)){
        scmageck <- data.frame(score = score[, i], p_val = p_val[, i])
        scmageck$Difference <- "non"
        scmageck$Difference[scmageck$score > score_cut & scmageck$p_val < p_val_cut] <- "up"
        scmageck$Difference[scmageck$score < -score_cut & scmageck$p_val < p_val_cut] <- "down"
        a <- plyr::count(scmageck$Difference)
        rownames(a) <- a$x
        for(j in a$x){
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
    new3 <- subset(new, factor %in% head(new2[order(as.numeric(new2$gene_numbers), decreasing = T), ], top)$factor)
#    new3[which(new3$Difference == "down"), ]$gene_numbers <- -new3[which(new3$Difference == "down"), ]$gene_numbers
    
    #get range of y axis
    
    y_max <- max(new3$gene_numbers)
    if(ceiling(y_max/100) == 0){
        y_max <- y_max
    }else{
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
        theme(plot.title = element_text(hjust = 0.5, size = 25),
              text = element_text(hjust = 0.5, face = "bold"),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 20),
              axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 16, color = "black"),
              axis.title.x = element_text(size = 20), 
              axis.title.y = element_text(size = 20),
              axis.text.y = element_text(hjust = 0.5, size = 16, color = "black")) + 
        ggsci::scale_fill_npg()
    } else {
        p3 <- ggplot(new3, aes(x = factor, y = gene_numbers)) +
        geom_bar(stat = 'identity', aes(fill = Difference), position = position_dodge(0.9), colour = "black") +
        theme_test() +
        labs(x = "Perturbation", y = "Gene Numbers") +
        theme(text = element_text(hjust = 0.5, face = "bold"),
              legend.text = element_text(size = 16),
              axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 16, color = "black"),
              axis.title.x = element_text(size = 20), 
              axis.title.y = element_text(size = 20),
              axis.text.y = element_text(hjust = 0.5, size = 16, color = "black")) + 
        coord_cartesian(ylim = c(0, y_break[1])) + 
        scale_y_continuous(breaks = c(0, y_break[1])) + 
        ggsci::scale_fill_npg()
        p4 <- ggplot(new3, aes(x = factor, y = gene_numbers)) +
        geom_bar(stat = 'identity', aes(fill = Difference), position = position_dodge(0.9), colour = "black") +
        theme_test() +
        labs(x = NULL, y = NULL, title = project) +
        theme(plot.title = element_text(hjust = 0.5, size = 25),
              text = element_text(hjust = 0.5, face = "bold"),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 20),
              axis.text.y = element_text(hjust = 0.5, size = 16, color = "black"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) + 
        coord_cartesian(ylim = c(y_break[2], y_max)) + 
        scale_y_continuous(breaks = c(y_break[2], y_max)) + 
        ggsci::scale_fill_npg()
        p1 <- ggpubr::ggarrange(p4, p3, heights = height, ncol = 1, nrow = 2, 
                                common.legend = TRUE, legend="right", align = "v") 
    }
    
    #save plot
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "perturbation_efficiency")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    dir <- file.path(prefix, "pdf")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    dir <- file.path(dir, "perturbation_efficiency")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    pdf(file.path(dir,
                  paste(label, "DE_gene_cutoff", score_cut, "_p", p_val_cut, ".pdf", sep = "")))
    print(p1)
    dev.off()
    
    png(file.path(img_dir, 
                  paste(label, "DE_gene_cutoff", score_cut, "_p", p_val_cut, ".png", sep = "")), 
        width = 600 * 3, height = 600 * 3, res = 72 * 3)
    print(p1)
    dev.off()
    
    return(p1)
}
