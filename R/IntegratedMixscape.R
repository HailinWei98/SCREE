#' Integrated Mixscape Function
#' 
#' Integrated function to evaluate perturbation efficiency based on \code{\link[Seurat]{CalcPerturbSig}} and \code{\link[Seurat]{RunMixscape}}.
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows. Note that the dataset has to be normalized and scaled. 
#' @param sg_lib Data frame or directory to a txt file containing 3 columns: cell, barcode, gene. If sgRNA information stored in a matrix-like format or input data frame only has sgRNA frequency of each cell, use \code{\link[SCREE]{sgRNAassign}} to assign sgRNA to each cell.
#' @param NTC The name of negative controls. Default is "NTC".
#' @param sg.to.use Use all sgRNA in the sg_lib to calculate perturb signature and perturbation enrichment or only use the sgRNA of cells assigned unique sgRNA, can be one of "all", "single". Default is "all".
#' @param sg.split String to split sgRNA name. Default is "_sgRNA", which means sgRNA named in the format like "gene_sgRNA1".
#' @param mixscape.only Logical, only perform the mixscape step or not. Default is \code{FALSE}, if \code{TRUE}, the input SeuratObject should calculate perturb signature first. Users can set this parameter to \code{TRUE} if you don't want to calculate perturb signature again, for the same dataset.
#' @param perturb.sig.only Logical, only calculate perturb signature or not. Default is \code{FALSE}.
#' @param nfeature Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'. Default is 2000.
#' @param selection.method How to choose top variable features. Choose one of : "vst", "mvp", "disp". See more details from the \code{selection.method} in \code{\link[Seurat]{FindVariableFeatures}}.
#' @param npcs Total Number of PCs to compute and store. Default is 50.
#' @param dims Dimensions of reduction to use as input. Default is 1:40.
#' @param assays Assay to used for calculating perturb signature. Default is "RNA".
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Default is 1.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Default is 0.8.
#' @param title.size Numeric, title size of the UMAP. Default is 16.
#' @param legend.key.size Numeric, legend key size of the UMAP. Default is unit(0.7, "cm").
#' @param legend.text.size Numeric, legend text size of the UMAP. Default is 9.
#' @param x.text.size Numeric, x-axis text size of the UMAP. Default is 10.
#' @param x.title.size Numeric, x-axis title size of the UMAP. Default is 12.
#' @param y.text.size Numeric, y-axis text size of the UMAP. Default is 10.
#' @param y.title.size Numeric, y-axis title size of the UMAP. Default is 12.
#' @param pt.size Point size of the UMAP. Default is 0.2.
#' @param raster Logical, convert points to raster format, will be useful to reduce the storage cost of the output figure or pdf. Default is \code{FALSE}.
#' @param label.cut Maximum perturbations to be labeled in the UMAP. Default is 20, which means to label top 20 perturbations ranked by cell numbers.
#' @param NTC.cal Logical, calculate negative control or not. Default is \code{TRUE}, which will remove cells assigned with negative control from the SeuratObject.
#' @param table.save Logical, save perturbation enrichment table or not. Default is \code{TRUE}.
#' @param top Top perturbations with the highest enrichment ratio to visualize. Default is \code{NULL}, using all perturbations to visualize.
#' @param range Range of the ratio. Since for some datasets, the max enrichment ratio among all perturbations is extremely small, the range can be set to a smaller limitation to enhance the enrichment. Default is c(0, 1).
#' @param color Color for heatmap. Default is c("white", "coral1"), corresponding to the score from 0 to 1 or other ranges.
#' @param cell Cell width and height settings for heatmap. If set this parameter to a single number or \code{NA}, cell width and height will be the same, and \code{NA} means the cell will be adaptive to the input perturbation numbers. If set this parameter to a vector including two numbers or \code{NA}, the first one will be cell width and the second one will be height. 
#' @param fontsize Fontsize of the heatmap. Default is 15.
#' @param angle Angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315). Default is 90.
#' @param legend.title Include the title of legend in the heatmap or not. Default is FALSE.
#' @param plot.show Logical, show the UMAP plot or not. Default is \code{TRUE}.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom utils read.table write.table
#' @import Seurat
#' @import ggplot2
#' @importFrom patchwork plot_layout
#' @importFrom plyr count
#' @import reshape2
#' @export

IntegratedMixscape<- function(mtx, sg_lib, NTC = "NTC", sg.to.use = "all", sg.split = "_sgRNA", mixscape.only = FALSE, perturb.sig.only = FALSE, nfeature = 2000, selection.method = "vst", npcs = 50, dims = 1 : 40, assays = "RNA", algorithm = 1, resolution = 0.8, title.size = 16, legend.key.size = unit(0.7, "cm"), legend.text.size = 9, x.text.size = 10, x.title.size = 12, y.text.size = 10, y.title.size = 12, pt.size = 0.2, raster = FALSE, label.cut = 20, NTC.cal = TRUE, table.save = TRUE, top = NULL, range = c(0, 1), color = c("white", "coral1"), cell = c(15, 20), fontsize = 15, angle = 90, legend.title = FALSE, plot.show = TRUE, plot.save = TRUE, prefix = ".", label = "", width = 7, height = 7, png_res = 720){
    
    #set custom theme of plot
    
    custom_theme <- theme(
        plot.title = element_text(size = title.size, hjust = 0.5),
        legend.key.size = legend.key.size,
        legend.text = element_text(size = legend.text.size))

    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        eccite <- readRDS(mtx)
    } else {
        eccite <- mtx
    }
    
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        sg_lib <- read.table(sg_lib, header = T)
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
            sg_lib <- subset(sg_lib, cell %in% subset(sg_in_cell, freq == 1)[ ,"cell"])
            rownames(sg_lib) <- sg_lib$cell
            eccite <- eccite[ ,colnames(eccite) %in% sg_lib$cell]
        }

        # Prepare RNA assay for dimensionality reduction:

        eccite <- umap(mtx = eccite, nfeature = nfeature, selection.method = selection.method, npcs = npcs, dims = dims, assays = assays, algorithm = algorithm, reduction.prefix = "", label.cut = label.cut, resolution = resolution, raster = raster, title.size = title.size, legend.key.size = unit(0.1, "pt"), legend.text.size = 8, x.text.size = x.text.size, x.title.size = x.title.size, y.text.size = y.text.size, y.title.size = y.title.size, pt.size = pt.size, plot.show = FALSE, plot.save = TRUE, prefix = prefix, label = label, plot.return = TRUE, width = width / 2, height = height / 2)
        
        p_p <- eccite[[2]]
        p_c <- eccite[[3]]
        p <- p_c / p_p + patchwork::plot_layout(guides = 'auto')
        eccite <- eccite[[1]]
        
        perturb_ratio <- CalculatePerturbEnrichment(mtx = eccite, sg_lib = sg_lib, NTC = NTC, NTC.cal = NTC.cal, prefix = prefix, 
                                                    label = paste(label, "before_", sep = ""), top = top, range = range, 
                                                    color = color, cell = cell, fontsize = fontsize, 
                                                    angle = angle, legend.title = legend.title, 
                                                    plot.save = TRUE, table.save = TRUE)

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
            pt.size = pt.size,
            reduction = "umap", cols = "Dark2", repel = T, raster = raster) +
        scale_color_brewer(palette = "Dark2") +
        ggtitle("Biological Replicate") +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        custom_theme

        p2 <- DimPlot(
            object = eccite,
            group.by = 'Phase',
            label = F, 
            pt.size = pt.size,
            reduction = "umap", repel = T, raster = raster) +
        ggtitle("Cell Cycle Phase") +
        ylab("UMAP 2") +
        xlab("UMAP 1") +
        custom_theme

        p3 <- DimPlot(
            object = eccite,
            group.by = 'gene',
            pt.size = pt.size,
            reduction = "umap",
            split.by = "gene",
            ncol = 1,
            cols = c("grey39","goldenrod3"), raster = raster) +
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

        eccite <- umap(mtx = eccite, nfeature = nfeature, selection.method = selection.method, npcs = npcs, dims = dims, assays = "PRTB", algorithm = algorithm, reduction.prefix = "prtb", label.cut = label.cut, resolution = resolution, raster = raster, title.size = title.size, legend.key.size = unit(0.1, "pt"), legend.text.size = 8, x.text.size = x.text.size, x.title.size = x.title.size, y.text.size = y.text.size, y.title.size = y.title.size, pt.size = pt.size, plot.show = FALSE, plot.save = FALSE, plot.return = TRUE)
        
        q_p <- eccite[[2]]
        q_c <- eccite[[3]]
        q <- q_c / q_p + patchwork::plot_layout(guides = 'auto')
        eccite <- eccite[[1]]

        perturb_ratio <- CalculatePerturbEnrichment(mtx = eccite, sg_lib = sg_lib, NTC = NTC, NTC.cal = NTC.cal, prefix = prefix, 
                                                    label = paste(label, "after_", sep = ""), top = top, range = range, 
                                                    color = color, cell = cell, fontsize = fontsize, 
                                                    angle = angle, legend.title = legend.title, 
                                                    plot.save = TRUE, table.save = TRUE)
        
        # Generate plots to check if clustering is driven by biological replicate ID,
        # cell cycle phase or target gene class.

        q1 <- DimPlot(
            object = eccite,
            group.by = 'replicate',
            reduction = 'prtbumap',
            pt.size = pt.size, cols = "Dark2", label = F, repel = T, raster = raster) +
        scale_color_brewer(palette = "Dark2") +
        ggtitle("Biological Replicate") +
        ylab("UMAP 2") +
        xlab("UMAP 1") +
        custom_theme

        q2 <- DimPlot(
            object = eccite,
            group.by = 'Phase',
            reduction = 'prtbumap',
            pt.size = pt.size, label = F, repel = T, raster = raster) +
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
            pt.size = pt.size,
            cols = c("grey39", "goldenrod3"), raster = raster) +
        ggtitle("Perturbation Status") +
        ylab("UMAP 2") +
        xlab("UMAP 1") +
        custom_theme

        #save plots.

        if (plot.save == TRUE) {
            
            dir <- file.path(prefix, "results")
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

            pdf_dir <- file.path(dir, "pdf")
            if (!(dir.exists(pdf_dir))) {
                dir.create(pdf_dir)
            }
            
            img_dir <- file.path(dir, "img")
            if (!(dir.exists(img_dir))) {
                dir.create(img_dir)
            }

            pdf(file.path(pdf_dir, paste(label, "mixscape_before.pdf", sep = "")), width = width, height = height)
            print(((p1 / p2 + patchwork::plot_layout(guides = 'auto')) | p3 ))
            dev.off()

            png(file.path(img_dir, paste(label, "mixscape_before.png", sep = "")), 
                width = width, height = height, unit = "in", res = png_res)
            print(((p1 / p2 + patchwork::plot_layout(guides = 'auto')) | p3 ))
            dev.off()


            pdf(file.path(pdf_dir, paste(label, "mixscape_after.pdf", sep = "")), width = width, height = height)
            print((q1 / q2 + patchwork::plot_layout(guides = 'auto') | q3))
            dev.off()

            png(file.path(img_dir, paste(label, "mixscape_after.png", sep = "")), 
                width = width, height = height, unit = "in", res = png_res)
            print(((q1 / q2 + patchwork::plot_layout(guides = 'auto')) | q3 ))
            dev.off()
            
            pdf(file.path(pdf_dir, paste(label, "cluster.pdf", sep = "")), width = width, height = height)
            print((p | q))
            dev.off()

            png(file.path(img_dir, paste(label, "cluster.png", sep = "")), 
                width = width, height = height, unit = "in", res = png_res)
            print((p | q))
            dev.off()
            
#             pdf(file = file.path(pdf_dir, paste(label, "before_umap_perturbations.pdf", sep = "")), 
#                 width = width / 2, height = height / 2)
#             print(p_p)
#             dev.off()

#             pdf(file = file.path(pdf_dir, paste(label, "before_umap_seurat_clusters.pdf", sep = "")), 
#                 width = width / 2, height = height / 2)
#             print(p_c)
#             dev.off()

#             png(file.path(img_dir, paste(label, "before_umap_perturbations.png", sep = "")), 
#                 width = width / 2, height = height / 2, unit = "in", res = png_res)
#             print(p_p)
#             dev.off()

#             png(file.path(img_dir, paste(label, "before_umap_seurat_clusters.png", sep = "")), 
#                 width = width / 2, height = height / 2, unit = "in", res = png_res)
#             print(p_c)
#             dev.off()
            
            pdf(file = file.path(pdf_dir, paste(label, "umap_perturbations.pdf", sep = "")), 
                width = width / 2, height = height / 2)
            print(q_p)
            dev.off()

            pdf(file = file.path(pdf_dir, paste(label, "umap_seurat_clusters.pdf", sep = "")), 
                width = width / 2, height = height / 2)
            print(q_c)
            dev.off()

            png(file.path(img_dir, paste(label, "umap_perturbations.png", sep = "")), 
                width = width / 2, height = height / 2, unit = "in", res = png_res)
            print(q_p)
            dev.off()

            png(file.path(img_dir, paste(label, "umap_seurat_clusters.png", sep = "")), 
                width = width / 2, height = height / 2, unit = "in", res = png_res)
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

            if (length(unique(eccite$mixscape_class.global)) == 3) {

                # Calculate percentage of KO cells for all target gene classes.

                df <- prop.table(table(eccite$mixscape_class.global, eccite$NT), 2)
                df2 <- reshape2::melt(df)
                df2$Var2 <- as.character(df2$Var2)
                test <- df2[which(df2$Var1 == "KO"), ]
                test <- test[order(test$value, decreasing = T), ]
                new.levels <- test$Var2
                df2$Var2 <- factor(df2$Var2, levels = new.levels )
                df2$Var1 <- factor(df2$Var1, levels = c(NTC, "NP", "KO"))
                df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = sg.split)[[1]][1])
                df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = sg.split)[[1]][2])
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
                        pdf(file.path(pdf_dir, paste(label, "mixscape_KO_percent.pdf", sep = "")))
                        print(q)
                        dev.off()

                        ko_dir <- file.path(img_dir, paste(label, "KO_percent", sep = ""))
                        if (!(dir.exists(ko_dir))) {
                            dir.create(ko_dir)
                        }
                        png(file.path(ko_dir, "1.png"), 
                            width = width, height = height, unit = "in", res = png_res)
                        print(q)
                        dev.off()
                    }

                } else {
                    plot_list <- list()
                    for (i in 1 : (ceiling(l/12))) {
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
                        ko_dir <- file.path(img_dir, paste(label, "KO_percent", sep = ""))
                        if (!(dir.exists(ko_dir))) {
                            dir.create(ko_dir)
                        }
                        pdf(file.path(pdf_dir, paste(label, "mixscape_KO_percent.pdf", sep = "")))
                        for (i in 1 : (ceiling(l/12))) {
                            png(file.path(ko_dir, paste(i, ".png", sep = "")), 
                                width = width, height = height, unit = "in", res = png_res)
                            print(plot_list[[i]])
                            dev.off()
                        }
                        dev.off()
                        for (i in 1 : (ceiling(l/12))) {
                            png(file.path(ko_dir, paste(i, ".png", sep = "")), 
                                width = width, height = height, unit = "in", res = png_res)
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

        if (length(unique(eccite$mixscape_class.global)) == 3) {

            # Calculate percentage of KO cells for all target gene classes.

            df <- prop.table(table(eccite$mixscape_class.global, eccite$NT), 2)
            df2 <- reshape2::melt(df)
            df2$Var2 <- as.character(df2$Var2)
            test <- df2[which(df2$Var1 == "KO"), ]
            test <- test[order(test$value, decreasing = T), ]
            new.levels <- test$Var2
            df2$Var2 <- factor(df2$Var2, levels = new.levels )
            df2$Var1 <- factor(df2$Var1, levels = c(NTC, "NP", "KO"))
            df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = sg.split)[[1]][1])
            df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = sg.split)[[1]][2])
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
                    pdf(file.path(pdf_dir, paste(label, "mixscape_KO_percent.pdf", sep = "")))
                    print(q)
                    dev.off()

                    ko_dir <- file.path(img_dir, paste(label, "KO_percent", sep = ""))
                    if (!(dir.exists(ko_dir))) {
                        dir.create(ko_dir)
                    }
                    png(file.path(ko_dir, "1.png"), 
                        width = width, height = height, unit = "in", res = png_res)
                    print(q)
                    dev.off()
                }

            } else {
                plot_list <- list()
                for (i in 1 : (ceiling(l/12))) {
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
                    pdf(file.path(pdf_dir, paste(label, "mixscape_KO_percent.pdf", sep = "")))
                    for (i in 1 : (ceiling(l/12))) {
                        png(file.path(ko_dir, paste(i, ".png", sep = "")), 
                            width =  width, height = height, unit = "in", res = png_res)
                        print(plot_list[[i]])
                        dev.off()
                    }
                    dev.off()
                    ko_dir <- file.path(img_dir, paste(label, "KO_percent", sep = ""))
                    if (!(dir.exists(ko_dir))) {
                        dir.create(ko_dir)
                    }
                    for (i in 1 : (ceiling(l/12))) {
                        png(file.path(ko_dir, paste(i, ".png", sep = "")), 
                            width = width, height = height, unit = "in", res = png_res)
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
