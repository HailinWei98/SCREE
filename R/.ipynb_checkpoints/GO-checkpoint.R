#' Gene Ontology Enrichment Analyse
#'
#' Identify enriched Gene Ontology (GO) terms using potential target genes for each perturbation.
#'
#' @param score Data frame or directory of score from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param pval Data frame or directory of p_value from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param score_cut Score cutoff of \code{improved_scmageck_lr} results. Default is 0.2.
#' @param pval_cut P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.
#' @param DE_gene_to_use Differential gene set to use, can be one of "all" (all potential target genes), "up" (potential target genes with positive score), "down" (potential target genes with negative score). Default is "all".
#' @param database The same as \code{OrgDb} in \code{\link[clusterProfiler]{enrichGO}}. Default is "org.Hs.eg.db".
#' @param gene_type Type of gene names in the score and pval table, can be one of "Symbol", "Ensembl". Default is "Symbol".
#' @param showCategory Numbers of each term to show. Default is 10.
#' @param wrap_width Positive integer giving target line width in characters. A width less than or equal to 1 will put each word on its own line. Default is 150.
#' @param title.size Numeric, title size of the barplot. Default is 25.
#' @param legend.text.size Numeric, legend text size of the barplot. Default is 12.
#' @param legend.title.size Numeric, legend title size of the barplot. Default is 18.
#' @param x.text.size Numeric, x-axis text size of the barplot. Default is 12.
#' @param x.title.size Numeric, x-axis title size of the barplot. Default is 20.
#' @param y.text.size Numeric, y-axis text size of the barplot. Default is 12.
#' @param y.title.size Numeric, y-axis title size of the barplot. Default is 20.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 12.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 8.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom utils read.table write.table
#' @importFrom clusterProfiler enrichGO slice
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import ggplot2
#' @import stringr
#' @export

GOenrichment <- function(score, pval, score_cut = 0.2, pval_cut = 0.05, DE_gene_to_use = "all", database = "org.Hs.eg.db", gene_type = "Symbol", showCategory = 10, wrap_width = 150, title.size = 25, legend.text.size = 12, legend.title.size = 18, x.text.size = 12, x.title.size = 20, y.text.size = 12, y.title.size = 20, plot.save = TRUE, prefix = ".", label = "", width = 12, height = 8, png_res = 720){
    
    #get score and p_val
    
    if (is.character(score)) {
        score <- read.table(score, header = T, row.names = 1)
    }
    
    if (is.character(pval_dir)) {
        pval <- read.table(pval, header = T, row.names = 1)
    }
    
    #get GO results for each perturbation
    
    results <- list()
    j = 0
    
    if (plot.save == TRUE) {
        
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }

        dir <- file.path(dir, "potential_target_gene")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }

        dir <- file.path(dir, "GO")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }
        
        pdf_dir <- file.path(dir, "GO")
        if (!(dir.exists(pdf_dir))) {
            dir.create(path = pdf_dir)
        }

        img_dir <- file.path(dir, "img")
        if (!(dir.exists(img_dir))) {
            dir.create(img_dir)
        }

    }

    #convert gene name
        
    if (gene_type == "Symbol") {
        genes <- bitr(rownames(score), "SYMBOL", "ENTREZID", database)
    } else if (gene_type == "Ensembl") {
        genes <- bitr(rownames(score), "ENSEMBL", "ENTREZID", database)
    } else {
        stop("gene_type must be one of c('Symbol', 'Ensembl')")
    }
    
    for (gene in colnames(score)) {
        
        #get score and p_val for each perturbation
        
        diff_table <- data.frame(score = score[, gene], p_val = pval[, gene])
        rownames(diff_table) <- rownames(score)
        
        #define up-regulated, down-regulated and all differential expressed gene list
        
        up_genes <- rownames(subset(diff_table, score > score_cut & p_val < pval_cut))
        down_genes <- rownames(subset(diff_table, score < -score_cut & p_val < pval_cut))
        DE_genes <- c(up_genes, down_genes)
        
        #prepare gene list for GO enrichment
        
        if (DE_gene_to_use == "all") {
            de <- DE_genes
            title <- "All Potential Targets of"
        } else if (DE_gene_to_use == "up") {
            de <- up_genes
            title <- "Potential Up-regulated Targets of"
        } else if (DE_gene_to_use == "down") {
            de <- down_genes
            title <- "Potential Down-regulated Targets of"
        } else {
            stop("DE_gene_to_use must be one of c('all', 'up', 'down')")
        }
        
        de_genes <- subset(genes, genes[ ,1] %in% de)
        
        if (nrow (de_genes) == 0) {
            next
        }
            
        #GO enrichment
            
        MF_ego <- enrichGO(gene = unique(de_genes[, 2]),
                           OrgDb = database,
                           ont = "MF",
                           keyType = "ENTREZID",
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
        BP_ego <- enrichGO(gene = unique(de_genes[, 2]),
                           OrgDb = database,
                           ont = "BP",
                           keyType = "ENTREZID",
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
        CC_ego <- enrichGO(gene = unique(de_genes[, 2]),
                           OrgDb = database,
                           ont = "CC",
                           keyType = "ENTREZID",
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
        
        if (!is.numeric(showCategory)) {
            stop("showCategory must be numeric")
        }
        MF_results <- na.omit(as.data.frame(MF_ego))
        BP_results <- na.omit(as.data.frame(BP_ego))
        CC_results <- na.omit(as.data.frame(CC_ego))
        
        if (nrow(MF_results) > 0 & 
            nrow(BP_results) > 0 & 
            nrow(CC_results) > 0) {
            
            MF_results <- na.omit(MF_results[1 : showCategory, ])
            BP_results <- na.omit(BP_results[1 : showCategory, ])
            CC_results <- na.omit(CC_results[1 : showCategory, ])
            all_results <- as.data.frame(rbind(MF_results, CC_results, BP_results))
            colnames(all_results) <- colnames(BP_results)
            all_results$ONTOLOGY <- factor(c(rep("MF", nrow(MF_results)), 
                                             rep("CC", nrow(CC_results)), 
                                             rep("BP", nrow(BP_results))), 
                                           levels = c("BP", "CC", "MF"), ordered = T)
            final_color <- c("#8DA1CB", "#FD8D62", "#66C3A5")
            
        } else if (nrow(MF_results) > 0 & 
                   nrow(BP_results) > 0 & 
                   nrow(CC_results) == 0) {
            
            MF_results <- na.omit(MF_results[1 : showCategory, ])
            BP_results <- na.omit(BP_results[1 : showCategory, ])
            all_results <- as.data.frame(rbind(MF_results, BP_results))
            colnames(all_results) <- colnames(BP_results)
                
            all_results$ONTOLOGY <- factor(c(rep("MF", nrow(MF_results)), 
                                             rep("BP", nrow(BP_results))), 
                                           levels = c("BP", "MF"), ordered = T)
            final_color <- c("#8DA1CB", "#66C3A5")
            
        } else if (nrow(MF_results) > 0 & 
                   nrow(BP_results) == 0 & 
                   nrow(CC_results) > 0) {
            
            MF_results <- na.omit(MF_results[1 : showCategory, ])
            CC_results <- na.omit(CC_results[1 : showCategory, ])
            all_results <- as.data.frame(rbind(MF_results, CC_results))
            colnames(all_results) <- colnames(CC_results)
                
            all_results$ONTOLOGY <- factor(c(rep("MF", nrow(MF_results)),
                                             rep("CC", nrow(CC_results))),
                                           levels = c("CC", "MF"), ordered = T)
            final_color <- c("#FD8D62", "#66C3A5")
            
        } else if (nrow(MF_results) == 0 & 
                   nrow(BP_results) > 0 & 
                   nrow(CC_results) > 0) {
            
            BP_results <- na.omit(BP_results[1 : showCategory, ])
            CC_results <- na.omit(CC_results[1 : showCategory, ])
            all_results <- as.data.frame(rbind(CC_results, BP_results))
            colnames(all_results) <- colnames(BP_results)
                
            all_results$ONTOLOGY <- factor(c(rep("CC", nrow(CC_results)),
                                             rep("BP", nrow(BP_results))), 
                                           levels = c("BP", "CC"), ordered = T)
            final_color <- c("#8DA1CB", "#FD8D62")
            
        } else if (nrow(MF_results) == 0 & 
                   nrow(BP_results) == 0 & 
                   nrow(CC_results) > 0) {
             
            all_results <- na.omit(CC_results[1 : showCategory, ])
            all_results$ONTOLOGY <- factor(rep("CC", nrow(all_results)), 
                                           levels = "CC", ordered = T)
            final_color <- "#FD8D62"
            
        } else if (nrow(MF_results) == 0 & 
                   nrow(BP_results) > 0 & 
                   nrow(CC_results) == 0) {
            
            all_results <- na.omit(BP_results[1 : showCategory, ])
            all_results$ONTOLOGY <- factor(rep("BP", nrow(all_results)), 
                                           levels = "BP", ordered = T)
            final_color <- "#8DA1CB"
            
        } else if (nrow(MF_results) > 0 & 
                   nrow(BP_results) == 0 & 
                   nrow(CC_results) == 0) {
            
            all_results <- na.omit(MF_results[1 : showCategory, ])
            all_results$ONTOLOGY <- factor(rep("MF", nrow(all_results)), 
                                           levels = "MF", ordered = T)
            final_color <- "#66C3A5"
            
        } else if (nrow(MF_results) == 0 & 
                   nrow(BP_results) == 0 & 
                   nrow(CC_results) == 0) {
            next
        }
           
        all_results$Description <- stringr::str_wrap(all_results$Description, width = wrap_width)
        
        #plot
        
        g1 <- ggplot(all_results) + 
        geom_bar(aes(x = Description, y = Count, fill = ONTOLOGY), stat = 'identity') +
        labs(x = "GO Terms", y = "Gene Numbers", title = paste(title, gene, sep = " "), fill = "Ontology") + 
        coord_flip() + 
        scale_x_discrete(limits = all_results$Description) + 
        theme_bw() +
        theme(panel.grid = element_blank()) + 
        scale_fill_manual(values = final_color) +
        theme(plot.title = element_text(hjust = 0.5, size = title.size), 
              axis.text.y = element_text(size = y.text.size),
              axis.text.x = element_text(size = x.text.size), 
              axis.title.y = element_text(size = y.title.size),
              axis.title.x = element_text(size = x.title.size),
              legend.title = element_text(size = legend.title.size),
              legend.text = element_text(size = legend.text.size),
              text = element_text(hjust = 0.5, face = "bold"))
        
        #save results and return
        
        j = j + 1
        results[[j]] <- g1
        names(results)[j] <- gene
        
        if (plot.save == TRUE) {
            
            pdf(file = file.path(pdf_dir, paste(label, gene, ".pdf", sep = "")), 
                width = width, height = height)
            print(g1)
            dev.off()

            png(file.path(img_dir, paste(label, gene, ".png", sep = "")), 
                width = width, height = height, unit = "in", res = png_res)
            print(g1)
            dev.off()
            
        }

    }
    return(results)
}