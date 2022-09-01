#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import ggplot2
#' @import stringr
#' @export

GOenrichment <- function(score_dir, pval_dir, p_val_cut = 0.05, score_cut = 0.5, DE_gene_to_use = "all",
                         database = "org.Hs.eg.db", gene_type = "Symbol", showCategory = 10, 
                         prefix = "./", label = ""){
    
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
    
    #get GO results for each perturbation
    
    results <- list()
    j = 0
    
    dir <- file.path(prefix, "pdf")
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
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "potential_target_gene")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "GO")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }

    #convert gene name
        
    if (gene_type == "Symbol") {
        genes <- bitr(rownames(score), "SYMBOL", "ENTREZID", database)
    } else if (gene_type == "Ensembl") {
        genes <- bitr(rownames(score), "ENSEMBL", "ENTREZID", database)
    } else {
        stop("gene_type must be one of c('Symbol', 'Ensembl')")
    }
    
    for(gene in colnames(score)){
        
        #get score and p_val for each perturbation
        
        diff_table <- data.frame(score = score[, gene], p_val = p_val[, gene])
        rownames(diff_table) <- rownames(score)
        
        #define up-regulated, down-regulated and all differential expressed gene list
        
        up_genes <- rownames(subset(diff_table, score > score_cut & p_val < p_val_cut))
        down_genes <- rownames(subset(diff_table, score < -score_cut & p_val < p_val_cut))
        DE_genes <- c(up_genes, down_genes)
        
        #prepare gene list for GO enrichment
        
        if(DE_gene_to_use == "all") {
            de <- DE_genes
            title <- "All Potential Targets of"
        } else if(DE_gene_to_use == "up"){
            de <- up_genes
            title <- "Potential Up-regulated Targets of"
        } else if(DE_gene_to_use == "down"){
            de <- down_genes
            title <- "Potential Down-regulated Targets of"
        } else{
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
            
        } else if(nrow(MF_results) > 0 & 
                  nrow(BP_results) == 0 & 
                  nrow(CC_results) > 0){
            
            MF_results <- na.omit(MF_results[1 : showCategory, ])
            CC_results <- na.omit(CC_results[1 : showCategory, ])
            all_results <- as.data.frame(rbind(MF_results, CC_results))
            colnames(all_results) <- colnames(CC_results)
                
            all_results$ONTOLOGY <- factor(c(rep("MF", nrow(MF_results)),
                                             rep("CC", nrow(CC_results))),
                                           levels = c("CC", "MF"), ordered = T)
            final_color <- c("#FD8D62", "#66C3A5")
            
        } else if(nrow(MF_results) == 0 & 
                  nrow(BP_results) > 0 & 
                  nrow(CC_results) > 0){
            
            BP_results <- na.omit(BP_results[1 : showCategory, ])
            CC_results <- na.omit(CC_results[1 : showCategory, ])
            all_results <- as.data.frame(rbind(CC_results, BP_results))
            colnames(all_results) <- colnames(BP_results)
                
            all_results$ONTOLOGY <- factor(c(rep("CC", nrow(CC_results)),
                                             rep("BP", nrow(BP_results))), 
                                           levels = c("BP", "CC"), ordered = T)
            final_color <- c("#8DA1CB", "#FD8D62")
            
        } else if(nrow(MF_results) == 0 & 
                  nrow(BP_results) == 0 & 
                  nrow(CC_results) > 0){
            
            all_results <- na.omit(CC_results[1 : showCategory, ])
            all_results$ONTOLOGY <- factor(rep("CC", nrow(all_results)), 
                                           levels = "CC", ordered = T)
            final_color <- "#FD8D62"
            
        } else if(nrow(MF_results) == 0 & 
                  nrow(BP_results) > 0 & 
                  nrow(CC_results) == 0){
            
            all_results <- na.omit(BP_results[1 : showCategory, ])
            all_results$ONTOLOGY <- factor(rep("BP", nrow(all_results)), 
                                           levels = "BP", ordered = T)
            final_color <- "#8DA1CB"
            
        } else if(nrow(MF_results) > 0 & 
                  nrow(BP_results) == 0 & 
                  nrow(CC_results) == 0){
            
            all_results <- na.omit(MF_results[1 : showCategory, ])
            all_results$ONTOLOGY <- factor(rep("MF", nrow(all_results)), 
                                           levels = "MF", ordered = T)
            final_color <- "#66C3A5"
            
        } else if(nrow(MF_results) == 0 & 
                  nrow(BP_results) == 0 & 
                  nrow(CC_results) == 0) {
            next
        }
           
        all_results$Description <- stringr::str_wrap(all_results$Description, width = 150)
        
        #plot
        
        g1 <- ggplot(all_results) + 
            geom_bar(aes(x = Description, y = Count, fill = ONTOLOGY), stat = 'identity') +
            labs(x = "GO Terms", y = "Gene Numbers", title = paste(title, gene, sep = " "), fill = "Ontology") + 
            coord_flip() + 
            scale_x_discrete(limits = all_results$Description) + 
            theme_bw() +
            theme(panel.grid = element_blank()) + 
            scale_fill_manual(values = final_color) +
            theme(plot.title = element_text(hjust = 0.5, size = 25), 
                  axis.text.y = element_text(size = 12),
                  axis.text.x = element_text(size = 12), 
                  axis.title.y = element_text(size = 20),
                  axis.title.x = element_text(size = 20),
                  legend.title = element_text(size = 18),
                  legend.text = element_text(size = 12),
                  text = element_text(hjust = 0.5, face = "bold"))
        
        #save results and return
        
        j = j + 1
        results[[j]] <- g1
        names(results)[j] <- gene
        
        pdf(file = file.path(dir, paste(gene, ".pdf", sep = "")), width = 16, height = 8)
        print(g1)
        dev.off()
        
        png(file.path(img_dir, paste(gene, ".png", sep = "")), 
            width = 1200 * 3, height = 600 * 3, res = 72 * 3)
        print(g1)
        dev.off()
    }
    return(results)
}