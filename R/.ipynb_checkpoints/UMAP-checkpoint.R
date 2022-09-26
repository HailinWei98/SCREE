#' Run UMAP and Visualization
#'
#' Run UMAP from SeuratObject after normalization and scale, and perform visualization, based on \code{\link[Seurat]}. 
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows. Note that the dataset has to be normalized and scaled, and need a sgRNA information column named "perturbations" in the meta data.
#' @param assays Assay to use. Default is "RNA".
#' @param nfeature Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'. Default is 2000.
#' @param selection.method How to choose top variable features. Choose one of : "vst", "mvp", "disp". See more details from the \code{selection.method} in \code{\link[Seurat]{FindVariableFeatures}}.
#' @param npcs Total Number of PCs to compute and store. Default is 50.
#' @param dims Dimensions of reduction to use as input. Default is 1:40.
#' @param reduction.prefix Prefix of the "pca" and "umap" reduction. Default is "".
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Default is 1.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Default is 0.8.
#' @param title.size Numeric, title size of the UMAP. Default is 25.
#' @param legend.key.size Numeric, legend key size of the UMAP. Default is unit(0.7, "cm").
#' @param legend.text.size Numeric, legend text size of the UMAP. Default is 14.
#' @param x.text.size Numeric, x-axis text size of the UMAP. Default is 16.
#' @param x.title.size Numeric, x-axis title size of the UMAP. Default is 20.
#' @param y.text.size Numeric, y-axis text size of the UMAP. Default is 16.
#' @param y.title.size Numeric, y-axis title size of the UMAP. Default is 20.
#' @param pt.size Point size of the UMAP. Default is 0.2.
#' @param raster Logical, convert points to raster format, will be useful to reduce the storage cost of the output figure or pdf. Default is \code{FALSE}.
#' @param label.cut Maximum perturbations to be labeled in the UMAP. Default is 20, which means to label top 20 perturbations ranked by cell numbers.
#' @param plot.show Logical, show the UMAP plot or not. Default is \code{TRUE}.
#' @param plot.return Logical, return a list include the SeuratObject after running UMAP and two UMAP plots or not. Default is \code{FALSE}.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @import ggplot2
#' @import Seurat
#' @export

umap <- function(mtx, assays = "RNA", nfeature = 2000, selection.method = "vst", npcs = 50, dims = 1:40, reduction.prefix = "", algorithm = 1, resolution = 0.8, title.size = 25, legend.key.size = unit(0.7, "cm"), legend.text.size = 14, x.text.size = 16, x.title.size = 20, y.size.size = 16, y.title.size = 20, pt.size = 0.2, raster = FALSE, label.cut = 20, plot.show = TRUE, plot.return = FALSE, plot.save = TRUE, prefix = ".", label = "", width = 7, height = 7, png_res = 720){
    
    #read file
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        perturb_QC <- readRDS(mtx)
    } else {
        perturb_QC <- mtx
    }

    #perturbation information
    
    if (!("perturbations" %in% colnames(perturb_QC@meta.data))) {
        stop("Please add perturbation information to the SeuratObject first.")
    }
    
    #visualize sgRNA label of each cell in UMAP
    
    perturb_QC <- FindVariableFeatures(perturb_QC, assay = assays, selection.method = selection.method, nfeatures = nfeature)
    perturb_QC <- RunPCA(object = perturb_QC, assay = assays, features = perturb_QC@assays[[assays]]@var.features, npcs = npcs, 
                         reduction.key = paste(reduction.prefix, "pca", sep = ""), 
                         reduction.name = paste(reduction.prefix, "pca", sep = ""))
    perturb_QC <- FindNeighbors(perturb_QC, assay = assays, dims = dims, 
                                reduction = paste(reduction.prefix, "pca", sep = ""))
    perturb_QC <- FindClusters(perturb_QC, resolution = resolution, algorithm = algorithm)
    perturb_QC <- RunUMAP(perturb_QC, dims = dims, assay = assays, 
                          reduction = paste(reduction.prefix, "pca", sep = ""), 
                          reduction.key = paste(reduction.prefix, "umap", sep = ""), 
                          reduction.name = paste(reduction.prefix, "umap", sep = ""))
    
    custom_theme <- theme(
        plot.title = element_text(size = title.size, hjust = 0.5),
        legend.key.size = legend.key.size,
        legend.text = element_text(size = legend.text.size), 
        axis.text.x = element_text(size = x.text.size), 
        axis.text.y = element_text(hjust = 0.5, size = y.text.size), 
        axis.title.x = element_text(size = x.title.size), 
        axis.title.y = element_text(size = y.title.size))
    
    p1 <- DimPlot(perturb_QC, 
                  reduction = paste(reduction.prefix, "umap", sep = ""), 
                  pt.size = pt.size,
                  group.by = "perturbations", 
                  raster = raster) + 
    ggtitle("Perturbations") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
    if(length(unique(perturb_QC$perturbations)) > label.cut) {
        p1 <- p1 + scale_color_discrete(breaks = 
                                        names(head(sort(table(perturb_QC$perturbations), decreasing = T), label.cut)))
    }
    
    p2 <- DimPlot(perturb_QC, 
                  reduction = paste(reduction.prefix, "umap", sep = ""), 
                  pt.size = pt.size, 
                  raster = raster) +
    ggtitle("Seurat clusters") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
    
    if (plot.save == TRUE) {
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }
    
        dir <- file.path(dir, "quality")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }
    
        dir <- file.path(dir, "UMAP")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }
        
        pdf_dir <- file.path(dir, "pdf")
        if (!(dir.exists(pdf_dir))) {
            dir.create(path = pdf_dir)
        }
    
        img_dir <- file.path(dir, "img")
        if (!(dir.exists(img_dir))) {
            dir.create(img_dir)
        }
        
        pdf(file = file.path(pdf_dir, paste(label, "umap_perturbations.pdf", sep = "")))
        print(p1)
        dev.off()
    
    
        pdf(file = file.path(pdf_dir, paste(label, "umap_seurat_clusters.pdf", sep = "")))
        print(p2)
        dev.off()
        
        png(file.path(img_dir, paste(label, "umap_perturbations.png", sep = "")), 
            width = 600, height = 600)
        print(p1)
        dev.off()
        
        png(file.path(img_dir, paste(label, "umap_seurat_clusters.png", sep = "")), 
            width = 600, height = 600)
        print(p2)
        dev.off()
    }
    
    if (plot.show == TRUE) {
        #print(p1 / p2 + patchwork::plot_layout(guides = 'auto'))
        print(list(p1, p2))
    }
    
    if (plot.return == FALSE) {
        return(perturb_QC)
    } else {
        #return(list(perturb_QC, p1 / p2 + patchwork::plot_layout(guides = 'auto')))
        return(list(perturb_QC, p1, p2))
    }

    
}

#' Run UMAP and Visualization
#'
#' Run UMAP from SeuratObject after normalization and scale, and perform visualization for scATAC-seq based input, based on \code{\link[Seurat]}. 
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows. Note that the dataset has to be normalized and scaled, and need a sgRNA information column named "perturbations" in the meta data.
#' @param assays Assay to use. Default is "peak".
#' @param min.cutoff Cutoff for feature to be included in the VariableFeatures for the object. This can be a percentile specified as 'q' followed by the minimum percentile, for example 'q5' to set the top 95% most common features as the VariableFeatures for the object. Alternatively, this can be an integer specifying the minimum number of cells containing the feature for the feature to be included in the set of VariableFeatures. For example, setting to 10 will include features in >10 cells in the set of VariableFeatures. If NULL, include all features in VariableFeatures. If NA, VariableFeatures will not be altered, and only the feature metadata will be updated with the total counts and percentile rank for each feature.
#' @param reduction.key Key for dimension reduction object. Default is "LSI_".
#' @param reduction.name Name for stored dimension reduction object. Default is "lsi".
#' @param nsvs Number of singular values to compute. Default is 50.
#' @param n Number of components to use. Default is 20. If NULL, use all components.
#' @param dims Dimensions of reduction to use as input. Default is \code{NULL}, which means use all components with a correlation < \code{dims.cut} between total counts, from \code{n} components used.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Default is 1.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Default is 0.8.
#' @param title.size Numeric, title size of the UMAP. Default is 25.
#' @param legend.key.size Numeric, legend key size of the UMAP. Default is unit(0.7, "cm").
#' @param legend.text.size Numeric, legend text size of the UMAP. Default is 14.
#' @param x.text.size Numeric, x-axis text size of the UMAP. Default is 16.
#' @param x.title.size Numeric, x-axis title size of the UMAP. Default is 20.
#' @param y.text.size Numeric, y-axis text size of the UMAP. Default is 16.
#' @param y.title.size Numeric, y-axis title size of the UMAP. Default is 20.
#' @param pt.size Point size of the UMAP. Default is 0.2.
#' @param raster Logical, convert points to raster format, will be useful to reduce the storage cost of the output figure or pdf. Default is \code{FALSE}.
#' @param label.cut Maximum perturbations to be labeled in the UMAP. Default is 20, which means to label top 20 perturbations ranked by cell numbers.
#' @param plot.show Logical, show the UMAP plot or not. Default is \code{TRUE}.
#' @param plot.return Logical, return a list include the SeuratObject after running UMAP and two UMAP plots or not. Default is \code{FALSE}.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @import Seurat
#' @import Signac
#' @export

ATACumap <- function(mtx, assays = "peak", min.cutoff = "q5", reduction.key = 'LSI_', reduction.name = 'lsi', nsvs = 50, n = 20, dims = NULL, dims.cut = 0.5, algorithm = 1, resolution = 0.8, title.size = 25, legend.key.size = unit(0.7, "cm"), legend.text.size = 14, x.text.size = 16, x.title.size = 20, y.text.size = 16, y.title.size = 20, pt.size = 0.2, raster = FALSE, label.cut = 20, plot.show = TRUE, plot.return = FALSE, plot.save = TRUE, prefix = ".", label = "", width = 7, height = 7, png_res = 720) {
    
    #read file
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        peak <- readRDS(mtx)
    } else {
        peak <- mtx
    }
    
    peak<- RunTFIDF(peak)
    peak<- FindTopFeatures(peak, assay = assays, min.cutoff = min.cutoff)
    peak <- RunSVD(
        object = peak,
        assay = assays,
        reduction.key = reduction.key,
        reduction.name = reduction.name,
        n = nsvs
    )
    
    a <- DepthCor(peak, assay = assays, reduction = reduction.name, n = n)
    
    if (is.null(dims)) {
        dims <- subset(a$data, abs(counts) < dims.cut)$Component
    }
    
    peak <- FindNeighbors(object = peak, reduction = reduction.name, assay = assays, dims = dims)
    
    peak <- FindClusters(peak, resolution = resolution, algorithm = algorithm)
    
    peak <- RunUMAP(peak, dims = dims, assay = assays, reduction = reduction.name)

    custom_theme <- theme(
        plot.title = element_text(size = title.size, hjust = 0.5),
        legend.key.size = legend.key.size,
        legend.text = element_text(size = legend.text.size), 
        axis.text.x = element_text(size = x.text.size), 
        axis.text.y = element_text(hjust = 0.5, size = y.text.size), 
        axis.title.x = element_text(size = x.title.size), 
        axis.title.y = element_text(size = y.title.size))
    
    p1 <- DimPlot(peak, reduction = "umap", group.by = "perturbations", pt.size = pt.size, raster = raster) + 
    ggtitle("Perturbations") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
    if(length(unique(peak$perturbations)) > label.cut) {
        p1 <- p1 + scale_color_discrete(breaks = names(head(sort(table(peak$perturbations), decreasing = T), label.cut)))
    }
    
    p2 <- DimPlot(peak, reduction = "umap", pt.size = pt.size, raster) +
    ggtitle("Seurat clusters") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
    
    if (plot.save == TRUE) {
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }
    
        dir <- file.path(dir, "ATAC_quality")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }
    
        dir <- file.path(dir, "UMAP")
        if (!(dir.exists(dir))) {
            dir.create(path = dir)
        }
        
        pdf_dir <- file.path(dir, "pdf")
        if (!(dir.exists(pdf_dir))) {
            dir.create(pdf_dir)
        }
    
        img_dir <- file.path(dir, "img")
        if (!(dir.exists(img_dir))) {
            dir.create(img_dir)
        }
        
        pdf(file = file.path(pdf_dir, paste(label, "umap_perturbations.pdf", sep = "")), width = width, height = height)
        print(p1)
        dev.off()
    
    
        pdf(file = file.path(pdf_dir, paste(label, "umap_seurat_clusters.pdf", sep = "")), width = width, height = height)
        print(p2)
        dev.off()
        
        png(file.path(img_dir, paste(label, "umap_perturbations.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(p1)
        dev.off()
        png(file.path(img_dir, paste(label, "umap_seurat_clusters.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(p2)
        dev.off()
    }
    
    if (plot.show == TRUE) {
        #print(p1 / p2  + patchwork::plot_layout(guides = 'auto'))
        print(list(p1, p2))
    }
    
    if (plot.return == FALSE) {
        return(peak)
    } else {
        return(list(peak, p1, p2))
    }
    

}
