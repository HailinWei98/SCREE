#' Gene Expression Plot
#'
#' Plot expression of genes near each enhancer.
#'
#' @param sg_lib Data frame or directory to a txt file containing 3 columns: cell, barcode, gene. If sgRNA information stored in a matrix-like format or input data frame only has sgRNA frequency of each cell, use \code{\link[SCREE]{sgRNAassign}} to assign sgRNA to each cell.
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param selected Enhancer regions to visualize gene expression. By default, all enhancers will be chosen.
#' @param species Only support "Hs" and "Mm". Default is "Hs".
#' @param version Version of the reference genome(Ensembl). Default is "v75".
#' @param gene_annotations Gene annotations stored in data frame format, including c("chromosome", "start", "end", "strand", "transcript") as colnames. Default is \code{NULL}, gene annotations are from \code{\link{ensembldb}}.
#' @param upstream The number of nucleotides upstream of the start site of selected region. Default is 2000000.
#' @param downstream The number of nucleotides downstream of the start site of selected region. Default is 2000000.
#' @param NTC The name of negative controls. Default is "NTC".
#' @param title.size Numeric, title size of each gene expression plot. Default is 25.
#' @param x.text.size Numeric, x-axis text size of each gene expression plot. Default is 25.
#' @param x.title.size Numeric, x-axis title size of each gene expression plot. Default is 22.
#' @param y.text.size Numeric, y-axis text size of each gene expression plot. Default is 25.
#' @param y.title.size Numeric, y-axis title size of each gene expression plot. Default is 22.
#' @param annotate.size Numeric, title size of the whole plot including 12 gene expression plot. Default is 40.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 18.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 20.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#' @param html_config Logical, generate and return a list includes the config character string of html. Default is \code{FALSE}.
#'
#' @importFrom utils read.table write.table
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom ensembldb genes
#' @import EnsDb.Hsapiens.v75
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Mmusculus.v75
#' @import EnsDb.Mmusculus.v79
#' @import Seurat
#' @import ggpubr
#' @importFrom cowplot plot_grid
#' @export

EnhancerGeneExpression <- function(sg_lib, mtx, selected = NULL, species = "Hs", version = "v75", gene_annotations = NULL, upstream = 2000000, downstream = 2000000, NTC = "NTC", title.size = 25, x.text.size = 25, x.title.size = 22, y.text.size = 25, y.title.size = 22, annotate.size = 40, plot.save = TRUE, prefix = ".", width = 18, height = 20, png_res = 720, html_config = FALSE) {
    
    # get genome annotations
    
    if (is.null(gene_annotations)) {
        if (species == "Hs") {
            if (version == "v75") {
                a <- genes(EnsDb.Hsapiens.v75)
            } else if (version == "v79") {
                a <- genes(EnsDb.Hsapiens.v79)
            } else if (version == "v86") {
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        } else if (species == "Mm") {
            if (version == "v75") {
                a <- genes(EnsDb.Mmusculus.v75)
            } else if (version == "v79") {
                a <- genes(EnsDb.Mmusculus.v79)
            }
        } else {
            stop("Species must be one of 'Hs' and 'Mm'")
        }

        gene_anno <- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
        
    } else {
        gene_anno <- gene_annotations
    }
    
    #read files
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        mtx <- readRDS(mtx)
    }
        
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        sg_lib <- read.table(sg_lib, header = T)
    }

    sg_lib <- subset(sg_lib, cell %in% intersect(sg_lib$cell, colnames(mtx)))
    if (is.null(selected)) {
        selected <- unique(sg_lib$gene)
        selected <- selected[grep("^chr", selected)]
    }
    
    if (plot.save == TRUE) {
        
        dir <- file.path(prefix, "results")
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

        pdf_dir <- file.path(dir, "pdf")
        if (!dir.exists(pdf_dir)) {
            dir.create(path = pdf_dir)
        }
        
        img_dir <- file.path(dir, "img")
        if (!dir.exists(img_dir)) {
            dir.create(path = img_dir)
        }
        
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
            if (l == 0) {
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
            
            if (plot.save == TRUE) {
                img_dir <- file.path(img_dir, new_select)
                if (!dir.exists(img_dir)) {
                    dir.create(path = img_dir)
                }
            }

            
            for(i in 1 : l) {
                gene <- gene_list[i]
                p <- VlnPlot(object = mtx_new,
                             features = gene,
                             pt.size = 0.1) +
                theme(plot.title = element_text(hjust = 0.5, size = title.size),
                      axis.text.x = element_text(size = x.text.size),
                      axis.title.x = element_text(size = x.title.size), 
                      axis.title.y = element_text(size = y.title.size),
                      axis.text.y = element_text(hjust = 0.5, size = y.text.size)) + 
                NoLegend()
                results[[i]] <- p
                names(results[i]) <- gene
            }
            
            #save plot
            
            if (plot.save == TRUE) {
                pdf(file = file.path(pdf_dir, paste(new_select, ".pdf", sep = "")), 
                    height = height, width = width)
                for (i in 1:(ceiling(l/12))) {
                    if (i == ceiling(l/12)) {
                        up <- l
                    } else {
                        up <- i * 12
                    }
                    p <- plot_grid(plotlist = results[(i * 12 - 11) : up], 
                                   ncol = 4, align = "hv") + 
                    theme(plot.margin = margin(20, 20, 20, 40))
                    p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_select,
                                                                             size = annotate.size, 
                                                                             face = "bold", 
                                                                             rot = 90,
                                                                             vjust = 1)) + 
                    theme(plot.margin = margin(20, 20, 20, 30))
                    print(p)
                }
                dev.off()
            }
            
            enhancer <- c()
            
            for(i in 1:(ceiling(l/12))){
                if(i == ceiling(l/12)) {
                    up <- l
                } else {
                    up <- i * 12
                }
                p <- plot_grid(plotlist = results[(i * 12 - 11) : up], 
                               ncol = 4, align = "hv") + 
                theme(plot.margin = margin(20, 20, 20, 40))
                p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_select, 
                                                                         size = annotate.size,
                                                                         face = "bold", 
                                                                         rot = 90, 
                                                                         vjust = 1)) + 
                theme(plot.margin = margin(20, 20, 20, 30))
                
                if (plot.save == TRUE) {
                    png(file.path(img_dir, paste(i, ".png", sep = "")), 
                        height = height, width = width, unit = "in", res = png_res)
                    print(p)
                    dev.off()
                }
                
                if (html_config == TRUE) {
                    enhancer <- c(enhancer, file.path("results", 
                                                      "enhancer_function", 
                                                      "enhancer_gene_expression", 
                                                      "img", 
                                                      paste(i, ".png", sep = "")))
                }
            }
            
            #generate html config

            if (html_config == TRUE) {
                names(enhancer) <- seq(1:(ceiling(l/12)))
                enhancer <- paste(names(enhancer), enhancer, collapse = "\" , \"", sep = "\" : \"")
                all_enhancer <- paste("\"enhancer\" : {\"", enhancer, "\"}", sep = "")
            }

        } else {
                
            #get results for all perturbations
                
            j <- 0
            results <- list()
            k <- 0
            all_enhancer <- c()
            
            for (perturb in selected) {
                
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
                if (l == 0) {
                    next
                }
                gene_list <- gene_model$transcript
            
                #rename enhancer
                
                new_perturb <- gsub("[:|-]", "_", perturb)
                
                result <- list()
                select_sg <- subset(sg_lib, gene == perturb)
                
                mtx_new <- subset(mtx, cells = c(unique(select_sg$cell), 
                                                 colnames(subset(mtx, perturbations == NTC))))
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
                
                for (i in 1:l) {
                    gene <- gene_list[i]
                    p <- VlnPlot(object = mtx_new,
                                 features = gene,
                                 pt.size = 0.1) +
                    theme(plot.title = element_text(hjust = 0.5, size = title.size),
                          axis.text.x = element_text(size = x.text.size),
                          axis.title.x = element_text(size = x.title.size), 
                          axis.title.y = element_text(size = y.title.size),
                          axis.text.y = element_text(hjust = 0.5, size = y.text.size)) + 
                    NoLegend()
                    
                    result[[i]] <- p
                    names(result[i]) <- gene
                
                }
                
                #save plot

                if (plot.save == TRUE) {
                    pdf(file = file.path(dir, paste(new_perturb, ".pdf", sep = "")), 
                        height = height, width = width)
                    for(i in 1:(ceiling(l/12))){
                        if(i == ceiling(l/12)) {
                            up <- l
                        } else {
                            up <- i * 12
                        }
                        p <- plot_grid(plotlist = result[(i * 12 - 11) : up], 
                                       ncol = 4, align = "hv") + 
                        theme(plot.margin = margin(20, 20, 20, 40))
                        p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_perturb,
                                                                                 size = annotate.size, 
                                                                                 face = "bold", 
                                                                                 rot = 90,
                                                                                 vjust = 1)) + 
                        theme(plot.margin = margin(20, 20, 20, 30))
                        print(p)
                    }
                    dev.off()
                }
                
                if (html_config == TRUE) {
                    
                    enhancer <- c()
                
                    #generate enhancer directory of html

                    enhancer_prefix <- file.path("results", 
                                                 "enhancer_function", 
                                                 "enhancer_gene_expression", 
                                                 "img", 
                                                 new_perturb)

                }
        
                for (i in 1:(ceiling(l/12))) {
                    if (i == ceiling(l/12)) {
                        up <- l
                    } else {
                        up <- i * 12
                    }
                    p <- plot_grid(plotlist = result[(i * 12 - 11) : up], 
                                   ncol = 4, align = "hv") + 
                    theme(plot.margin = margin(20, 20, 20, 40))
                    p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_perturb, 
                                                                             size = annotate.size,
                                                                             face = "bold", 
                                                                             rot = 90, 
                                                                             vjust = 1)) + 
                    theme(plot.margin = margin(20, 20, 20, 30))
                    if (plot.save == TRUE) {
                        png(file.path(dir2, paste(i, ".png", sep = "")), 
                            height = height, width = width, unit = "in", res = png_res)
                        print(p)
                        dev.off()
                    }
  
                    if (html_config == TRUE) {
                        enhancer <- c(enhancer, paste(file.path(enhancer_prefix, i), ".png", sep = ""))
                        names(enhancer)[i] <- i
                    }
                }
                
                j <- j + 1
                results[[j]] <- result
                names(results)[j] <- new_perturb
                
                if (html_config == TRUE) {
                    enhancer <- paste(names(enhancer), enhancer, collapse = "\" , \"", sep = "\" : \"")
                    enhancer <- paste("\"", new_perturb, "\" : {\"", enhancer, "\"}", sep = "")
                    all_enhancer <- c(all_enhancer, enhancer)
                }           
            }
            if (html_config == TRUE) {
                all_enhancer <- paste("", all_enhancer, collapse = ", ", sep = "")
                all_enhancer <- paste("\"enhancer\" : {", all_enhancer, "}", sep = "")
            }
        }
        
        if (html_config == TRUE) {
            return(list(results, all_enhancer))
        } else {
            return(results)
        }
    } else {
        stop("Please input correct format of selected perturbations")
    }
}

#' Calculate Direct Target Numbers
#'
#' Identify potential direct targets for each enhancer.
#'
#' @param score Data frame or directory of score from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param pval Data frame or directory of p_value from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param selected Enhancer regions to calculate direct target ratio. By default, all enhancers will be chosen.
#' @param species Only support "Hs" and "Mm". Default is "Hs".
#' @param version Version of the reference genome(Ensembl). Default is "v75".
#' @param gene_annotations Gene annotations stored in data frame format, including c("chromosome", "start", "end", "strand", "transcript") as colnames. Default is \code{NULL}, gene annotations are from \code{\link{ensembldb}}.
#' @param upstream The number of nucleotides upstream of the start site of selected region. Default is 2000000.
#' @param downstream The number of nucleotides downstream of the start site of selected region. Default is 2000000.
#' @param score_cut Score cutoff of \code{improved_scmageck_lr} results. Default is 0.2.
#' @param pval_cut P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.
#' @param table.save Logical, save the results as table or not. Default is \code{TRUE}. 
#' @param prefix Path to save the table. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#'
#' @importFrom utils read.table write.table
#' @importFrom ensembldb genes
#' @import EnsDb.Hsapiens.v75
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Mmusculus.v75
#' @import EnsDb.Mmusculus.v79
#' @export

DirectTarget <- function(score, pval, selected = NULL, species = "Hs", version = "v75", gene_annotations = NULL, upstream = 2000000, downstream = 2000000, score_cut = 0.2, pval_cut = 0.05, table.save = TRUE, prefix = ".", label = "") {
    
    # get genome annotations
    
    if (is.null(gene_annotations)) {
        if (species == "Hs") {
            if (version == "v75") {
                a <- genes(EnsDb.Hsapiens.v75)
            } else if (version == "v79") {
                a <- genes(EnsDb.Hsapiens.v79)
            } else if (version == "v86") {
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        } else if (species == "Mm") {
            if (version == "v75") {
                a <- genes(EnsDb.Mmusculus.v75)
            } else if (version == "v79") {
                a <- genes(EnsDb.Mmusculus.v79)
            }
        } else {
            stop("Species must be one of 'Hs' and 'Mm'")
        }

        gene_anno <- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
        
    } else {
        gene_anno <- gene_annotations
    }
    
    #read files
    
    if (is.character(score)) {
        score <- read.table(score, header = T)
    }
        
    if (is.character(pval)) {
        pval <- read.table(pval, header = T)
    }

    if (is.null(selected)) {
        selected <- colnames(score)
        selected <- selected[grep("^chr", selected)]
    }
    
    if (table.save == TRUE) {
        
        dir <- file.path(prefix, "results")
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
        
        dir <- file.path(prefix, "direct_targets")
        if (!dir.exists(dir)) {
            dir.create(path = dir)
        }
        
    }
    
    if (is.character(selected)) {
                 
        #get results for all perturbations

        i <- 0
        results <- data.frame()

        for (perturb in selected) {

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

            select_enhancer <- data.frame(score = score[, perturb], p_val = pval[, perturb])
            rownames(select_enhancer) <- rownames(score)
            de_enhancer <- subset(select_enhancer, abs(score) > score_cut & p_val < pval_cut)
            gene_around <- select_enhancer[gene_list, ]
            
            if (table.save == TRUE) {
                write.table(gene_around, 
                            file = file.path(dir, paste(label, new_perturb, ".txt", sep = "")), 
                            quote = FALSE)
            }

            direct_targets <- subset(gene_around, abs(score) > score_cut & p_val < pval_cut)
            targets <- data.frame(proximal = nrow(direct_targets), gene_around = nrow(gene_around), 
                                  distal = (nrow(de_enhancer) - nrow(direct_targets)))
            results <- rbind(results, targets)
        }
        
        colnames(results) <- c("proximal", "gene_around", "distal")
        
        rownames(results) <- selected
        
        return(results)
        
    } else {
        stop("Please input correct format of selected perturbations")
    }
}