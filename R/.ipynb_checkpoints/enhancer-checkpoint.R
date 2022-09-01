#' @import Seurat
#' @export

EnhancerGeneExpression <- function(sg_dir, mtx_dir, selected = NULL, prefix = "./", upstream = 2000000, 
                                   downstream = 2000000, gene_annotations = NULL, species = "Hs", version = "v75",
                                   NTC = "NTC", html_config = FALSE) {
    
    # get genome annotations
    
    if(is.null(gene_annotations)){
        if(species == "Hs"){
            if(version == "v75"){
                a <- genes(EnsDb.Hsapiens.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Hsapiens.v79)
            }else if(version == "v86"){
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        }else if(species == "Mm"){
            if(version == "v75"){
                a <- genes(EnsDb.Mmusculus.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Mmusculus.v79)
            }
        }
        gene_anno<- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
    } else {
        gene_anno <- gene_annotations
            }
    
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

    sg_lib <- subset(sg_lib, cell %in% intersect(sg_lib$cell, colnames(mtx)))
    if(is.null(selected)){
        selected <- unique(sg_lib$gene)
        selected <- selected[grep("^chr", selected)]
    }
    
    dir <- file.path(prefix, "pdf")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "enhancer_function")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "enhancer_gene_expression")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }

    img_dir <- file.path(prefix, "img")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
        
    img_dir <- file.path(img_dir, "enhancer_function")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
        
    img_dir <- file.path(img_dir, "enhancer_gene_expression")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
    
    if(is.character(selected)){
        
        if(length(selected) == 1){
                                  
            #get information from selected enhancer
    
            all <- unlist(strsplit(selected, "[.|:|-|_]"))
            chr <- all[1]
            start <- as.numeric(all[2])
            end <- as.numeric(all[3])
        
            #extend region
    
            minbp <- start - upstream
            maxbp <- end + downstream
            
            #get gene model
            
            gene_model <- gene_anno
            gene_model <- gene_model[!is.na(gene_model$chromosome) & 
                                     !is.na(gene_model$start) &
                                     !is.na(gene_model$end) &
                                     !is.na(gene_model$strand) &
                                     !is.na(gene_model$transcript), ]
            gene_model <- gene_model[gene_model$chromosome == chr &
                                     ((gene_model$start > minbp & gene_model$start < maxbp) |
                                      (gene_model$end > minbp & gene_model$end < maxbp) |
                                      (gene_model$start < minbp & gene_model$end > maxbp)), ]
            gene_model <- gene_model[gene_model$transcript %in% rownames(mtx), ]
            
            l <- nrow(gene_model)
            if(l == 0){
                stop("Cannot find genes close to selected region.")
            }
            gene_list <- gene_model$transcript
            
            results <- list()
            select_sg <- subset(sg_lib, gene == selected)
            
            #rename enhancer
            
            new_select <- gsub("[:|-|.]", "_", selected)
            
            mtx_new <- subset(mtx, cells = c(unique(select_sg$cell), colnames(subset(mtx, perturbations == NTC))))
            mtx_new$perturbations <- gsub(selected,
                                          "single enhancer",
                                          mtx_new$perturbations)
            mtx_new$perturbations <- gsub("multiple", "multiple enhancer", mtx_new$perturbations)
            mtx_new$perturbations <- factor(mtx_new$perturbations, 
                                            levels = c("single enhancer", 
                                                       "multiple enhancer", NTC))
            mtx_new@active.ident<- mtx_new$perturbations           
            
            img_dir <- file.path(img_dir, new_select)
            if (!dir.exists(img_dir)) {
                dir.create(path = img_dir)
            }
            
            for(i in 1 : l) {
                gene <- gene_list[i]
                p <- VlnPlot(object = mtx_new,
                             features = gene,
                             pt.size = 0.1) +
                theme(plot.title = element_text(hjust = 0.5, size = 25),
                      axis.text.x = element_text(size = 25),
                      axis.title.x = element_text(size = 22), 
                      axis.title.y = element_text(size = 22),
                      axis.text.y = element_text(hjust = 0.5, size = 25)) + 
                NoLegend()
                results[[i]] <- p
                names(results[i]) <- gene
            }
            
            #save plot

            pdf(file = file.path(dir, paste(new_select, ".pdf", sep = "")), height = 20 * 2/3, width = 18 * 3/4)
            for(i in 1:(ceiling(l/6))){
                if(i == ceiling(l/6)) {
                    up <- l
                } else {
                    up <- i * 6
                }
                p <- plot_grid(plotlist = results[(i * 6 - 5) : up], 
                               ncol = 3, align = "hv") + 
                theme(plot.margin = margin(20, 20, 20, 40))
                p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_select,
                                                                         size = 40, 
                                                                         face = "bold", 
                                                                         rot = 90,
                                                                         vjust = 1)) + 
                theme(plot.margin = margin(20, 20, 20, 30))
                print(p)
            }
            dev.off()
        
            enhancer <- c()
            
            for(i in 1:(ceiling(l/6))){
                if(i == ceiling(l/6)) {
                    up <- l
                } else {
                    up <- i * 6
                }
                p <- plot_grid(plotlist = results[(i * 6 - 5) : up], 
                               ncol = 4, align = "hv") + 
                theme(plot.margin = margin(20, 20, 20, 40))
                p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_select, 
                                                                         size = 40,
                                                                         face = "bold", 
                                                                         rot = 90, 
                                                                         vjust = 1)) + 
                theme(plot.margin = margin(20, 20, 20, 30))
                png(file.path(img_dir, paste(i, ".png", sep = "")), 
                    height = 1400 * 2/3, width = 1300 * 3/4, res = 72)
                print(p)
                dev.off()
                
                enhancer <- c(enhancer, file.path("img/enhancer_function/enhancer_gene_expression", 
                                                  paste(i, ".png", sep = "")))
            }
            
            #generate html config

            names(enhancer) <- seq(1:(ceiling(l/6)))
            enhancer <- paste(names(enhancer), enhancer, collapse = "\" , \"", sep = "\" : \"")
            all_enhancer <- paste("\"enhancer\" : {\"", enhancer, "\"}", sep = "")
            
        } else {
                
            #get results for all perturbations
                
            j <- 0
            results <- list()
            k <- 0
            all_enhancer <- c()
            
            for(perturb in selected){
                
                #get information from selected enhancer
                
                all <- unlist(strsplit(perturb, "[.|:|-|_]"))
                chr <- all[1]
                start <- as.numeric(all[2])
                end <- as.numeric(all[3])
        
                #extend region
    
                minbp <- start - upstream
                maxbp <- end + downstream
                
                #get gene model
            
                gene_model <- gene_anno
                gene_model <- gene_model[!is.na(gene_model$chromosome) & 
                                         !is.na(gene_model$start) &
                                         !is.na(gene_model$end) &
                                         !is.na(gene_model$strand) &
                                         !is.na(gene_model$transcript), ]
                gene_model <- gene_model[gene_model$chromosome == chr &
                                         ((gene_model$start > minbp & gene_model$start < maxbp) |
                                          (gene_model$end > minbp & gene_model$end < maxbp) |
                                          (gene_model$start < minbp & gene_model$end > maxbp)), ]
                gene_model <- gene_model[gene_model$transcript %in% rownames(mtx), ]
            
                l <- nrow(gene_model)
                if(l == 0){
                    next
                }
                gene_list <- gene_model$transcript
            
                #rename enhancer
                
                new_perturb <- gsub("[:|-]", "_", perturb)
                
                result <- list()
                select_sg <- subset(sg_lib, gene == perturb)
                
                mtx_new <- subset(mtx, cells = c(unique(select_sg$cell), colnames(subset(mtx, perturbations == NTC))))
                mtx_new$perturbations <- gsub(perturb,
                                              "single enhancer",
                                              mtx_new$perturbations)
                mtx_new$perturbations <- gsub("multiple", "multiple enhancer", mtx_new$perturbations)
                mtx_new$perturbations <- factor(mtx_new$perturbations, 
                                                levels = c("single enhancer", 
                                                           "multiple enhancer", NTC))
                mtx_new@active.ident<- mtx_new$perturbations
                
                #create save path
                
                dir2 <- file.path(img_dir, new_perturb)
                if (!dir.exists(dir2)) {
                    dir.create(path = dir2)
                }
                
                for(i in 1:l) {
                    gene <- gene_list[i]
                    p <- VlnPlot(object = mtx_new,
                                 features = gene,
                                 pt.size = 0.1) +
                    theme(plot.title = element_text(hjust = 0.5, size = 25),
                          axis.text.x = element_text(size = 25),
                          axis.title.x = element_text(size = 22), 
                          axis.title.y = element_text(size = 22),
                          axis.text.y = element_text(hjust = 0.5, size = 25)) + 
                    NoLegend()
                    
                    result[[i]] <- p
                    names(result[i]) <- gene
                
                }
                
                #save plot

                pdf(file = file.path(dir, paste(new_perturb, ".pdf", sep = "")), height = 20 * 2/3, width = 18 * 3/4)
                for(i in 1:(ceiling(l/6))){
                    if(i == ceiling(l/6)) {
                        up <- l
                    } else {
                        up <- i * 6
                    }
                    p <- plot_grid(plotlist = result[(i * 6 - 5) : up], 
                                   ncol = 4, align = "hv") + 
                    theme(plot.margin = margin(20, 20, 20, 40))
                    p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_perturb,
                                                                             size = 40, 
                                                                             face = "bold", 
                                                                             rot = 90,
                                                                             vjust = 1)) + 
                    theme(plot.margin = margin(20, 20, 20, 30))
                    print(p)
                }
                dev.off()
        
                enhancer <- c()
                
                #generate enhancer directory of html
                
                enhancer_prefix <- file.path("img/enhancer_function/enhancer_gene_expression", new_perturb)
                
                for(i in 1:(ceiling(l/6))){
                    if(i == ceiling(l/6)) {
                        up <- l
                    } else {
                        up <- i * 6
                    }
                    p <- plot_grid(plotlist = result[(i * 6 - 5) : up], 
                                   ncol = 4, align = "hv") + 
                    theme(plot.margin = margin(20, 20, 20, 40))
                    p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_perturb, 
                                                                             size = 40,
                                                                             face = "bold", 
                                                                             rot = 90, 
                                                                             vjust = 1)) + 
                    theme(plot.margin = margin(20, 20, 20, 30))
                    png(file.path(dir2, paste(i, ".png", sep = "")), 
                        height = 1400 * 2/3, width = 1300 * 3/4, res = 72)
                    print(p)
                    dev.off()
                    
                    enhancer <- c(enhancer, paste(file.path(enhancer_prefix, i), ".png", sep = ""))
                    names(enhancer)[i] <- i
                }
                
                j <- j + 1
                results[[j]] <- result
                names(results)[j] <- new_perturb
                
                enhancer <- paste(names(enhancer), enhancer, collapse = "\" , \"", sep = "\" : \"")
                enhancer <- paste("\"", new_perturb, "\" : {\"", enhancer, "\"}", sep = "")
                all_enhancer <- c(all_enhancer, enhancer)
                
            }
            all_enhancer <- paste("", all_enhancer, collapse = ", ", sep = "")
            all_enhancer <- paste("\"enhancer\" : {", all_enhancer, "}", sep = "")
        }
        if (html_config == TRUE) {
            return(list(results, all_enhancer))
        } else {
            return(results)
        }
    }else{
        stop("Please input correct format of selected perturbations")
    }
}
    
#' @import Seurat
#' @export

DirectTargetRatio <- function(score_dir, pval_dir, selected = NULL, score_cut = 0.5, p_val_cut = 0.05, 
                              prefix = "./", upstream = 2000000, downstream = 2000000, gene_annotations = NULL, 
                              species = "Hs", version = "v75") {
        # get genome annotations
    
    if(is.null(gene_annotations)){
        if(species == "Hs"){
            if(version == "v75"){
                a <- genes(EnsDb.Hsapiens.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Hsapiens.v79)
            }else if(version == "v86"){
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        }else if(species == "Mm"){
            if(version == "v75"){
                a <- genes(EnsDb.Mmusculus.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Mmusculus.v79)
            }
        }
        gene_anno<- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
    } else {
        gene_anno <- gene_annotations
            }
    
    #read files
    
    if (is.character(score_dir)) {
        score <- read.table(score_dir, header = T)
    } else {
        score <- score_dir
    }
        
    if (is.character(pval_dir)) {
        p_val <- read.table(pval_dir, header = T)
    } else {
        p_val <- pval_dir
    }

    if(is.null(selected)){
        selected <- colnames(score)
        selected <- selected[grep("^chr", selected)]
    }
    
    dir <- file.path(prefix, "direct_targets")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }
    
    if(is.character(selected)){
                 
        #get results for all perturbations

        i <- 0
        results <- data.frame()

        for(perturb in selected){

            #get information from selected enhancer

            all <- unlist(strsplit(perturb, "[.|:|-|_]"))
            chr <- all[1]
            start <- as.numeric(all[2])
            end <- as.numeric(all[3])

            #extend region

            minbp <- start - upstream
            maxbp <- end + downstream

            #get gene model

            gene_model <- gene_anno
            gene_model <- gene_model[!is.na(gene_model$chromosome) & 
                                     !is.na(gene_model$start) &
                                     !is.na(gene_model$end) &
                                     !is.na(gene_model$strand) &
                                     !is.na(gene_model$transcript), ]
            gene_model <- gene_model[gene_model$chromosome == chr &
                                     ((gene_model$start > minbp & gene_model$start < maxbp) |
                                      (gene_model$end > minbp & gene_model$end < maxbp) |
                                      (gene_model$start < minbp & gene_model$end > maxbp)), ]
            gene_model <- gene_model[gene_model$transcript %in% rownames(score), ]

            gene_list <- gene_model$transcript

            #rename enhancer

            new_perturb <- gsub("[:|-]", "_", perturb)

            select_enhancer <- data.frame(score = score[, perturb], p_val = p_val[, perturb])
            rownames(select_enhancer) <- rownames(score)
            de_enhancer <- subset(select_enhancer, abs(score) > score_cut & p_val < p_val_cut)
            gene_around <- select_enhancer[gene_list, ]
            write.table(gene_around, file = file.path(dir, paste(new_perturb, ".txt", sep = "")), quote = FALSE)
            direct_targets <- subset(gene_around, abs(score) > score_cut & p_val < p_val_cut)
            targets <- data.frame(proximal = nrow(direct_targets), gene_around = nrow(gene_around), 
                                  distal = (nrow(de_enhancer) - nrow(direct_targets)))
            results <- rbind(results, targets)
        }
        
        colnames(results) <- c("proximal", "distal")
        
        rownames(results) <- selected
        
        return(results)
        
    } else {
        stop("Please input correct format of selected perturbations")
    }
}