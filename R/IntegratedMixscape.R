#' function definitions ##### Integrated some functions in Seurat as mentioned in
#' vignette of Mixscape.Prefer Seurat object after QC.
#' @export

IntegratedMixscape<- function(mtx_dir, sg_dir, sg.to.use = "all", perturb.sig.only = FALSE, mixscape.only = FALSE,
                              nfeature = 2000, dims = 1 : 40, assays = "RNA", algorithm = 1,
                              label.cut = 20, resolution = 0.8,
                              NTC = "NTC", prefix = "./", label = "", plot.save = TRUE, plot.show = TRUE){
    
    #set custom theme of plot
    
    custom_theme <- theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.key.size = unit(0.7, "cm"),
        legend.text = element_text(size = 14))

    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        eccite <- readRDS(mtx_dir)
    } else {
        eccite <- mtx_dir
    }
    
    if (is.character(sg_dir)) {
        message(paste("Reading sgRNA lib file:", sg_dir))
        sg_lib <- read.table(sg_dir, header = T)
    } else {
        sg_lib <- sg_dir
    }

    #replicate information
    
    if(!("replicate" %in% colnames(eccite@meta.data))){
        stop("Please add replicate information to the SeuratObject first.")
    }
    
    #perturbation information
    
    if(!("perturbations" %in% colnames(eccite@meta.data))){
        stop("Please add perturbation information to the SeuratObject first.")
    }
    
    if (mixscape.only == FALSE) {
        
        #add barcode information
    
        if (sg.to.use == "single") {
            sg_in_cell <- data.frame(plyr::count(sg_lib$cell))
            sg_lib <- subset(sg_lib, cell %in% subset(sg_in_cell, freq == 1)[ ,1])
            rownames(sg_lib) <- sg_lib$cell
            eccite <- eccite[ ,colnames(eccite) %in% sg_lib$cell]
        }

        # Prepare RNA assay for dimensionality reduction:

        eccite <- umap(eccite, nfeature, dims, assays, algorithm, reduction.prefix = "", 
                       prefix, label, label.cut, resolution, plot.show = FALSE, plot.save = TRUE, plot.mixscape = TRUE)
        
        p_p <- eccite[[2]]
        p_c <- eccite[[3]]
        p <- p_c / p_p + patchwork::plot_layout(guides = 'auto')
        eccite <- eccite[[1]]
        
        perturb_ratio <- CalculatePerturbEnrichment(eccite, sg_lib, NTC = NTC, prefix = prefix, 
                                                    label = paste(label, "before_", sep = ""), 
                                                    plot.save = TRUE, table.save = TRUE, top = label.cut)

        #calculate cell cycle gene

        s.genes <- cc.genes.updated.2019$s.genes
        g2m.genes <- cc.genes.updated.2019$g2m.genes
        eccite <- CellCycleScoring(eccite, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

        #transform label

        eccite$gene <- eccite$perturbations
        eccite$gene <- as.character(eccite$gene)
        eccite$gene[which(eccite$gene != NTC)] <- "perturbed"

        # Generate plots to check if clustering is driven by biological replicate ID,
        # cell cycle phase or target gene class.

        p1 <- DimPlot(
            object = eccite,
            group.by = 'replicate',
            label = F,
            pt.size = 0.2,
            reduction = "umap", cols = "Dark2", repel = T) +
        scale_color_brewer(palette = "Dark2") +
        ggtitle("Biological Replicate") +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        custom_theme

        p2 <- DimPlot(
            object = eccite,
            group.by = 'Phase',
            label = F, 
            pt.size = 0.2,
            reduction = "umap", repel = T) +
        ggtitle("Cell Cycle Phase") +
        ylab("UMAP 2") +
        xlab("UMAP 1") +
        custom_theme

        p3 <- DimPlot(
            object = eccite,
            group.by = 'gene',
            pt.size = 0.2,
            reduction = "umap",
            split.by = "gene",
            ncol = 1,
            cols = c("grey39","goldenrod3")) +
        ggtitle("Perturbation Status") +
        ylab("UMAP 2") +
        xlab("UMAP 1") +
        custom_theme

        if (plot.show == TRUE) {
            print(((p1 / p2 + patchwork::plot_layout(guides = 'auto')) | p3 ))
        }

        #calculate perturb signature
        eccite <- CalcPerturbSig(
            object = eccite,
            assay = "RNA",
            gd.class ="perturbations",
            nt.cell.class = NTC,
            reduction = "pca",
            ndims = 40,
            num.neighbors = 20,
            split.by = "replicate",
            new.assay.name = "PRTB")

        #### Prepare PRTB assay for dimensionality reduction:
        # Normalize data, find variable features and center data.

        DefaultAssay(object = eccite) <- 'PRTB'

        # Use variable features from RNA assay.

        VariableFeatures(object = eccite) <- VariableFeatures(object = eccite[[assays]])
        eccite <- ScaleData(object = eccite, do.scale = F, do.center = T)

        eccite <- umap(eccite, nfeature, dims, assays = "PRTB", algorithm, reduction.prefix = "prtb",
                       prefix, label, label.cut, resolution, plot.show = FALSE, plot.save, plot.mixscape = TRUE)
        
        q_p <- eccite[[2]]
        q_c <- eccite[[3]]
        q <- q_c / q_p + patchwork::plot_layout(guides = 'auto')
        eccite <- eccite[[1]]

        perturb_ratio <- CalculatePerturbEnrichment(eccite, sg_lib, NTC = NTC, prefix = prefix,
                                                    label = paste(label, "after_", sep = ""), 
                                                    plot.save = TRUE, table.save = TRUE, top = label.cut)
        
        # Generate plots to check if clustering is driven by biological replicate ID,
        # cell cycle phase or target gene class.

        q1 <- DimPlot(
            object = eccite,
            group.by = 'replicate',
            reduction = 'prtbumap',
            pt.size = 0.2, cols = "Dark2", label = F, repel = T) +
        scale_color_brewer(palette = "Dark2") +
        ggtitle("Biological Replicate") +
        ylab("UMAP 2") +
        xlab("UMAP 1") +
        custom_theme

        q2 <- DimPlot(
            object = eccite,
            group.by = 'Phase',
            reduction = 'prtbumap',
            pt.size = 0.2, label = F, repel = T) +
        ggtitle("Cell Cycle Phase") +
        ylab("UMAP 2") +
        xlab("UMAP 1") +
        custom_theme

        q3 <- DimPlot(
            object = eccite,
            group.by = 'gene',
            reduction = 'prtbumap',
            split.by = "gene",
            ncol = 1,
            pt.size = 0.2,
            cols = c("grey39","goldenrod3")) +
        ggtitle("Perturbation Status") +
        ylab("UMAP 2") +
        xlab("UMAP 1") +
        custom_theme

        #save plots.

        if (plot.save == TRUE) {
            dir <- file.path(prefix, "pdf")
            if (!(dir.exists(dir))) {
                dir.create(dir)
            }

            dir <- file.path(dir, "perturbation_efficiency")
            if (!(dir.exists(dir))) {
                dir.create(dir)
            }

            dir <- file.path(dir, "mixscape")
            if (!(dir.exists(dir))) {
                dir.create(dir)
            }

            img_dir <- file.path(prefix, "img")
            if (!(dir.exists(img_dir))) {
                dir.create(img_dir)
            }

            img_dir <- file.path(img_dir, "perturbation_efficiency")
            if (!(dir.exists(img_dir))) {
                dir.create(img_dir)
            }

            img_dir <- file.path(img_dir, "mixscape")
            if (!(dir.exists(img_dir))) {
                dir.create(img_dir)
            }

            pdf(file.path(dir, paste(label, "mixscape_before.pdf", sep = "")))
            print(((p1 / p2 + patchwork::plot_layout(guides = 'auto')) | p3 ))
            dev.off()

            png(file.path(img_dir, paste(label, "mixscape_before.png", sep = "")), 
                width = 600 * 3, height = 600 * 3, res = 72 * 3)
            print(((p1 / p2 + patchwork::plot_layout(guides = 'auto')) | p3 ))
            dev.off()


            pdf(file.path(dir, paste(label, "mixscape_after.pdf", sep = "")))
            print((q1 / q2 + patchwork::plot_layout(guides = 'auto') | q3))
            dev.off()

            png(file.path(img_dir, paste(label, "mixscape_after.png", sep = "")), 
                width = 600 * 3, height = 600 * 3, res = 72 * 3)
            print(((q1 / q2 + patchwork::plot_layout(guides = 'auto')) | q3 ))
            dev.off()
            
            pdf(file.path(dir, paste(label, "cluster.pdf", sep = "")))
            print((p | q))
            dev.off()

            png(file.path(img_dir, paste(label, "cluster.png", sep = "")), 
                width = 600 * 3, height = 600 * 3, res = 72 * 3)
            print((p | q))
            dev.off()
            
            pdf(file = file.path(dir, paste(label, "umap_perturbations.pdf", sep = "")))
            print(q_p)
            dev.off()


            pdf(file = file.path(dir, paste(label, "umap_seurat_clusters.pdf", sep = "")))
            print(q_c)
            dev.off()

            png(file.path(img_dir, paste(label, "umap_perturbations.png", sep = "")), 
                width = 600, height = 600)
            print(q_p)
            dev.off()

            png(file.path(img_dir, paste(label, "umap_seurat_clusters.png", sep = "")), 
                width = 600, height = 600)
            print(q_c)
            dev.off()
        }

        if (plot.show == TRUE) {
            print(((q1 / q2 + patchwork::plot_layout(guides = 'auto')) | q3 ))
            print((p | q))
        }

        if (perturb.sig.only == TRUE) {
            return(eccite)
        } else {
            
            # Run mixscape.
    
            sg_in_cell <- data.frame(plyr::count(sg_lib$cell))
            sg_lib <- subset(sg_lib, cell %in% subset(sg_in_cell, freq == 1)[ ,1])
            rownames(sg_lib) <- sg_lib$cell
            eccite <- eccite[ ,colnames(eccite) %in% sg_lib$cell]        
            
            NT <- data.frame(cell = colnames(eccite), nt = "a")
            NT$nt <- sg_lib[NT$cell, ]$barcode
            eccite <- AddMetaData(eccite, as.factor(NT$nt), col.name = "NT")

            eccite <- suppressWarnings(RunMixscape(
                object = eccite,
                assay = "PRTB",
                slot = "scale.data",
                labels = "perturbations",
                nt.class.name = NTC,
                min.de.genes = 5,
                iter.num = 10,
                de.assay = "RNA",
                verbose = F,
                prtb.type = "KO"))

            if(length(unique(eccite$mixscape_class.global)) == 3){

                # Calculate percentage of KO cells for all target gene classes.

                df <- prop.table(table(eccite$mixscape_class.global, eccite$NT), 2)
                df2 <- reshape2::melt(df)
                df2$Var2 <- as.character(df2$Var2)
                test <- df2[which(df2$Var1 == "KO"), ]
                test <- test[order(test$value, decreasing = T), ]
                new.levels <- test$Var2
                df2$Var2 <- factor(df2$Var2, levels = new.levels )
                df2$Var1 <- factor(df2$Var1, levels = c(NTC, "NP", "KO"))
                df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "_sgRNA")[[1]][1])
                df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "_sgRNA")[[1]][2])
                df3 <- df2[-c(which(df2$gene == NTC)), ]
                df4 <- data.frame(unique(df3[, c(1, 3, 4)]))
                df5 <- subset(df3, gene %in% subset(df4,Var1 == "KO" & value != 0)$gene)

                #only remain genes with non-zero KO

                l <- length(unique(df5$gene))
                if (l < 12) {
                    p <- ggplot(df5, aes(x = guide_number, y = value * 100, fill = Var1)) +
                    geom_bar(stat = "identity") +
                    theme_classic() +
                    scale_fill_manual(values = c("grey49", "grey79", "coral1")) +
                    ylab("% of cells") +
                    xlab("sgRNA")
                    q <- p + theme(axis.text.x = element_text(size = 18, hjust = 1),
                                   axis.text.y = element_text(size = 18),
                                   axis.title = element_text(size = 16),
                                   strip.text = element_text(size = 8, face = "bold")) +
                    facet_wrap(vars(gene), ncol = 3, scales = "free") +
                    labs(fill = "mixscape class") + 
                    theme(legend.title = element_text(size = 14),
                          legend.text = element_text(size = 12))
                    
                    if (plot.show == TRUE) {
                        print(q)
                    }
                    
                    if (plot.save == TRUE) {
                        pdf(file.path(dir, "mixscape_KO_percent.pdf"))
                        print(q)
                        dev.off()

                        ko_dir <- file.path(img_dir, "KO_percent")
                        if (!(dir.exists(ko_dir))) {
                            dir.create(ko_dir)
                        }
                        png(file.path(ko_dir, "1.png"), 
                            width = 600 * 3, height = 600 * 3, res = 72 * 3)
                        print(q)
                        dev.off()
                    }

                } else {
                    plot_list <- list()
                    for(i in 1 : (ceiling(l/12))){
                        df6 <- subset(df5, gene %in% unique(df5$gene)[(i * 12 - 11) : (i * 12)])
                        p <- ggplot(df6, aes(x = guide_number, y = value * 100, fill = Var1)) +
                        geom_bar(stat = "identity") +
                        theme_classic() +
                        scale_fill_manual(values = c("grey49", "grey79", "coral1")) +
                        ylab("% of cells") +
                        xlab("sgRNA")
                        q <- p + theme(axis.text.x = element_text(size = 18, hjust = 1),
                                       axis.text.y = element_text(size = 18),
                                       axis.title = element_text(size = 16),
                                       strip.text = element_text(size = 8, face = "bold")) +
                        facet_wrap(vars(gene), ncol = 3, scales = "free") +
                        labs(fill = "mixscape class") + theme(legend.title = element_text(size = 14),
                                                       legend.text = element_text(size = 12))
                        plot_list[[i]] <- q
                    }
                    
                    if (plot.save == TRUE) {
                        ko_dir <- file.path(img_dir, "KO_percent")
                        if (!(dir.exists(ko_dir))) {
                            dir.create(ko_dir)
                        }
                        pdf(file.path(dir, "mixscape_KO_percent.pdf"))
                        for(i in 1 : (ceiling(l/12))){
                            png(file.path(ko_dir, paste(i, ".png", sep = "")), 
                                width = 600 * 3, height = 600 * 3, res = 72 * 3)
                            print(plot_list[[i]])
                            dev.off()
                        }
                        dev.off()
                        for(i in 1 : (ceiling(l/12))){
                            png(file.path(ko_dir, paste(i, ".png", sep = "")), 
                                width = 600 * 3, height = 600 * 3, res = 72 * 3)
                            print(plot_list[[i]])
                            dev.off()
                        }
                    }
                    
                    if (plot.show == TRUE) {
                        for (i in 1 : (ceiling(l/12))) {
                            print(plot_list[[i]])
                        }
                    }
                }
            } else {
                return(eccite)
            }
        }
    } else {
        
        # Run mixscape.
    
        sg_in_cell <- data.frame(plyr::count(sg_lib$cell))
        sg_lib <- subset(sg_lib, cell %in% subset(sg_in_cell, freq == 1)[ ,1])
        rownames(sg_lib) <- sg_lib$cell
        eccite <- eccite[ ,colnames(eccite) %in% sg_lib$cell]
        
        NT <- data.frame(cell = colnames(eccite), nt = "a")
        NT$nt <- sg_lib[NT$cell, ]$barcode
        eccite <- AddMetaData(eccite, as.factor(NT$nt), col.name = "NT")
        
        eccite <- suppressWarnings(RunMixscape(
            object = eccite,
            assay = "PRTB",
            slot = "scale.data",
            labels = "perturbations",
            nt.class.name = NTC,
            min.de.genes = 5,
            iter.num = 10,
            de.assay = "RNA",
            verbose = F,
            prtb.type = "KO"))

        if(length(unique(eccite$mixscape_class.global)) == 3){

            # Calculate percentage of KO cells for all target gene classes.

            df <- prop.table(table(eccite$mixscape_class.global, eccite$NT), 2)
            df2 <- reshape2::melt(df)
            df2$Var2 <- as.character(df2$Var2)
            test <- df2[which(df2$Var1 == "KO"), ]
            test <- test[order(test$value, decreasing = T), ]
            new.levels <- test$Var2
            df2$Var2 <- factor(df2$Var2, levels = new.levels )
            df2$Var1 <- factor(df2$Var1, levels = c(NTC, "NP", "KO"))
            df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "_sgRNA")[[1]][1])
            df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "_sgRNA")[[1]][2])
            df3 <- df2[-c(which(df2$gene == NTC)), ]
            df4 <- data.frame(unique(df3[, c(1, 3, 4)]))
            df5 <- subset(df3, gene %in% subset(df4,Var1 == "KO" & value != 0)$gene)

            #only remain genes with non-zero KO

            l <- length(unique(df5$gene))
            if (l < 12) {
                p <- ggplot(df5, aes(x = guide_number, y = value * 100, fill = Var1)) +
                geom_bar(stat = "identity") +
                theme_classic() +
                scale_fill_manual(values = c("grey49", "grey79", "coral1")) +
                ylab("% of cells") +
                xlab("sgRNA")
                q <- p + theme(axis.text.x = element_text(size = 18, hjust = 1),
                               axis.text.y = element_text(size = 18),
                               axis.title = element_text(size = 16),
                               strip.text = element_text(size = 8, face = "bold")) +
                facet_wrap(vars(gene), ncol = 3, scales = "free") +
                labs(fill = "mixscape class") + theme(legend.title = element_text(size = 14),
                                                 legend.text = element_text(size = 12))

                if (plot.show == TRUE) {
                    print(q)
                }

                if (plot.save == TRUE) {
                    pdf(file.path(dir, "mixscape_KO_percent.pdf"))
                    print(q)
                    dev.off()

                    ko_dir <- file.path(img_dir, "KO_percent")
                    if (!(dir.exists(ko_dir))) {
                        dir.create(ko_dir)
                    }
                    png(file.path(ko_dir, "1.png"), 
                        width = 600 * 3, height = 600 * 3, res = 72 * 3)
                    print(q)
                    dev.off()
                }

            } else {
                plot_list <- list()
                for(i in 1 : (ceiling(l/12))){
                    df6 <- subset(df5, gene %in% unique(df5$gene)[(i * 12 - 11) : (i * 12)])
                    p <- ggplot(df6, aes(x = guide_number, y = value * 100, fill = Var1)) +
                    geom_bar(stat = "identity") +
                    theme_classic() +
                    scale_fill_manual(values = c("grey49", "grey79", "coral1")) +
                    ylab("% of cells") +
                    xlab("sgRNA")
                    q <- p + theme(axis.text.x = element_text(size = 18, hjust = 1),
                                   axis.text.y = element_text(size = 18),
                                   axis.title = element_text(size = 16),
                                   strip.text = element_text(size = 8, face = "bold")) +
                    facet_wrap(vars(gene), ncol = 3, scales = "free") +
                    labs(fill = "mixscape class") + theme(legend.title = element_text(size = 14),
                                                   legend.text = element_text(size = 12))
                    plot_list[[i]] <- q
                }

                if (plot.save == TRUE) {
                    pdf(file.path(dir, "mixscape_KO_percent.pdf"))
                    for(i in 1 : (ceiling(l/12))){
                        png(file.path(ko_dir, paste(i, ".png", sep = "")), 
                            width = 600 * 3, height = 600 * 3, res = 72 * 3)
                        print(plot_list[[i]])
                        dev.off()
                    }
                    dev.off()
                    ko_dir <- file.path(img_dir, "KO_percent")
                    if (!(dir.exists(ko_dir))) {
                        dir.create(ko_dir)
                    }
                    for(i in 1 : (ceiling(l/12))){
                        png(file.path(ko_dir, paste(i, ".png", sep = "")), 
                            width = 600 * 3, height = 600 * 3, res = 72 * 3)
                        print(plot_list[[i]])
                        dev.off()
                    }
                }

                if (plot.show == TRUE) {
                    for (i in 1 : (ceiling(l/12))) {
                        print(plot_list[[i]])
                    }
                }
            }
        } else {
            return(eccite)
        }
    }
    return(eccite)
}
