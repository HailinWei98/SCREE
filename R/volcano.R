#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @export

volcano <- function(score_dir, pval_dir, p_val_cut = 0.05, score_cut = 0.5, 
                    showCategory = 10, prefix = "./", label = ""){
    
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
    
    #create save path
    
    results <- list()
    j = 0
    
    dir <- file.path(prefix, "pdf")
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
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(path = img_dir)
    }
    
    img_dir <- file.path(img_dir, "perturbation_efficiency")
    if (!(dir.exists(img_dir))) {
        dir.create(path = img_dir)
    }
    
    img_dir <- file.path(img_dir, "volcano")
    if (!(dir.exists(img_dir))) {
        dir.create(path = img_dir)
    }
    
    #prepare color for plot
    
    diff <- c("#4DBBD5FF","grey","#E64B35FF")
    names(diff) <- c("down", "non", "up")

    #plot for each perturbation
    
    for (i in colnames(score)) {
        
        #classify genes for each perturbation
        
        scmageck <- data.frame(score = score[, i], p_val = p_val[, i])
        rownames(scmageck) <- rownames(score)
        scmageck$diff <- "non"
        scmageck$diff[scmageck$score > score_cut & scmageck$p_val < p_val_cut] <- "up"
        scmageck$diff[scmageck$score < -score_cut & scmageck$p_val < p_val_cut] <- "down"
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
        p1 <- ggplot(data = scmageck, 
                     mapping = aes(x = score, y = -log10(p_val + 10e-6), colour = diff, label = lab)) + 
        geom_point(size = 2) + 
        theme_test() + 
        labs(title = title) + 
        theme(plot.title = element_text(hjust = 0.5, size = 23), 
              text = element_text(hjust = 0.5, face = "bold"),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 20),
              axis.title.x = element_text(size = 20), 
              axis.title.y = element_text(size = 20),
              axis.text.y = element_text(hjust = 0.5, size = 16),
              axis.text.x = element_text(size = 16)) +  
        scale_color_manual(values = diff[sort(unique(scmageck$diff))]) + 
        geom_vline(xintercept = c(-score_cut, score_cut), linetype = "dotted") + 
        geom_text_repel(max.overlaps = 500, show.legend = FALSE) + 
        guides(colour = guide_legend(title = "Difference", title.hjust = 0.5))
        
        if (score_max < score_cut) {
            p1 <- p1 + xlim(-score_max, score_max)
        }
        #save results and return
        
        j = j + 1
        results[[j]] <- p1
        names(results)[j] <- i
        pdf(file = file.path(dir, paste(i, ".pdf", sep = "")), useDingbats = T)
        print(p1)
        dev.off()
        
        png(file.path(img_dir, paste(i, ".png", sep = "")), 
            width = 600 * 3, height = 600 * 3, res = 72 * 3)
        print(p1)
        dev.off()
    }
    return(results) 
}