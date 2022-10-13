#' Single Cell Quality Control
#'
#' Single cell quality control and visualization.
#' 
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param species Species used to calculate mitochondrial gene percentage, only support "Hs" and "Mm". Default is "Hs".
#' @param gene_frac A paramter for filtering low expressed genes. By default, only genes that have expressions or counts in at least that fractions of cells are kept. Default is 0.01.
#' @param nFeature Limitation of detected feature numbers in each cell, in the format like c(200, 50000). Default is c(200, 50000). If you don't need an upper limit, you can set the upper limit to an extremely large number. 
#' @param nCount Minimal count numbers in each cell. Default is 1000.
#' @param mt Maximum mitochondrial gene percentage of each cell. Default is 10 (means 10%).
#' @param blank_NTC Logical, use blank control as negative control or not. Default is \code{FALSE}.
#' @param title.size Numeric, title size of the Violin plot Default is 20.
#' @param x.text.size Numeric, x-axis text size of the Violin plot Default is 15.
#' @param x.title.size Numeric, x-axis title size of the Violin plot Default is 15.
#' @param y.text.size Numeric, y-axis text size of the Violin plot Default is 20.
#' @param pt.size Point size of the violin plot. Default is 0.1.
#' @param raster Logical, convert points to raster format, will be useful to reduce the storage cost of the output figure or pdf. Default is \code{FALSE}.
#' @param plot.show Logical, show the Violin plot or not. Default is \code{TRUE}.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 8.3.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 8.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom cowplot plot_grid
#' @import ggrastr
#' @import ggplot2
#' @export

scQC <- function(mtx, species = "Hs", gene_frac = 0.01, nFeature = c(200, 50000), nCount = 1000, mt = 10, blank_NTC = FALSE, title.size = 20, x.text.size = 15, x.title.size = 15, y.text.size = 20, pt.size = 0.1, raster = FALSE, plot.show = FALSE, plot.save = TRUE, prefix = ".", label = "", width = 8.3, height = 8, png_res = 720){
    
    #read file
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        perturb <- readRDS(mtx)
    } else {
        perturb <- mtx
    }

    if (!("replicate" %in% colnames(perturb@meta.data) & "perturbations" %in% colnames(perturb@meta.data))) {
        stop("Cannot find metadata named 'replicate' or 'perturbations' in input matrix.")
    }
    
    if (!("percent.mt" %in% colnames(perturb@meta.data))) {
        warning("No metadata named 'percent.mt', we will set the 'percent.mt' of each cell to 0")
        perturb$percent.mt <- 0
    }
    
    perturb$nFeature_RNA <- perturb[[paste("nFeature_", perturb@active.assay, sep = "")]][, 1]
    perturb$nCount_RNA <- perturb[[paste("nCount_", perturb@active.assay, sep = "")]][, 1]

    #filter cells with low quality
        
    if (blank_NTC == TRUE) {
        perturb_QC <- subset(perturb,
                             nFeature_RNA <= nFeature[2] &
                             nFeature_RNA >= nFeature[1] &
                             nCount_RNA >= nCount &
                             percent.mt <= mt)
    } else {
        perturb_QC <- subset(perturb,
                             nFeature_RNA <= nFeature[2] &
                             nFeature_RNA >= nFeature[1] &
                             nCount_RNA >= nCount &
                             percent.mt <= mt &
                             perturbations != 'blank')
    }
    
    perturb_QC <- CreateSeuratObject(counts = GetAssayData(object = perturb_QC, slot = "counts"),
                                     min.cells = gene_frac * ncol(perturb_QC), project = perturb@project.name)
    #calculate percent.mt
    
    if (species == "Hs") {
        perturb_QC[["percent.mt"]] <- PercentageFeatureSet(perturb_QC, pattern = "^MT-")
    } else {
        perturb_QC[["percent.mt"]] <- PercentageFeatureSet(perturb_QC, pattern = "^mt-")
    }
    
    for (meta in colnames(perturb@meta.data)) {
        if (!(meta %in% c("nFeature_RNA", "nCount_RNA", "percent.mt"))) {
            perturb_QC[[meta]] <- perturb[[meta]]
        }
    }
    
    perturb_QC@active.ident <- as.factor(perturb[, colnames(perturb_QC)]$orig.ident)
    perturb@active.ident <- as.factor(perturb$orig.ident)
    #QC plot of the single cell matrix
    
    #plot vlnplot with raster parameter
    
    custom_theme <- theme(plot.title = element_text(hjust = 0.5, 
                                                    size = title.size, 
                                                    face = "bold"),
                          axis.text.x = element_text(size = x.text.size, 
                                                     angle = 45, 
                                                     vjust = 1, 
                                                     hjust = 1, 
                                                     color = "black"),
                          axis.title.x = element_text(size = x.title.size), 
                          axis.text.y = element_text(hjust = 0.5, 
                                                     size = y.text.size, 
                                                     color = "black"))

    a <- perturb@meta.data
    a$Identity <- perturb@active.ident
    p1 <- ggplot(data = a, mapping = aes(x = Identity, y = nFeature_RNA)) + 
    geom_violin(aes(fill = factor(Identity))) + theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "nFeature_RNA") + 
    custom_theme
    
    p2 <- ggplot(data = a, mapping = aes(x = Identity, y = nCount_RNA)) + 
    geom_violin(aes(fill = factor(Identity))) + theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "nCount_RNA") + 
    custom_theme
    
    p3 <- ggplot(data = a, mapping = aes(x = Identity, y = percent.mt)) + 
    geom_violin(aes(fill = factor(Identity))) + theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "percent.mt") + 
    custom_theme
    
    #plot vlnplot with raster parameter, after QC

    a <- perturb_QC@meta.data
    a$Identity <- perturb_QC@active.ident
    q1 <- ggplot(data = a, mapping = aes(x = Identity, y = nFeature_RNA)) + 
    geom_violin(aes(fill = factor(Identity))) + theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "nFeature_RNA") + 
    custom_theme

    q2 <- ggplot(data = a, mapping = aes(x = Identity, y = nCount_RNA)) + 
    geom_violin(aes(fill = factor(Identity))) + theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "nCount_RNA") + 
    custom_theme
    
    q3 <- ggplot(data = a, mapping = aes(x = Identity, y = percent.mt)) + 
    geom_violin(aes(fill = factor(Identity))) + theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "percent.mt") + 
    custom_theme
    
    if (raster == TRUE) {
        p1 <- p1 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        p2 <- p2 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        p3 <- p3 + ggrastr::geom_jitter_rast(size = pt.size, height = 0, raster.dpi = getOption("ggrastr.default.dpi", 300))
        q1 <- q1 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        q2 <- q2 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        q3 <- q3 + ggrastr::geom_jitter_rast(size = pt.size, height = 0, raster.dpi = getOption("ggrastr.default.dpi", 300))
        
    } else {
        p1 <- p1 + geom_jitter(size = pt.size)
        p2 <- p2 + geom_jitter(size = pt.size)
        p3 <- p3 + geom_jitter(size = pt.size, height = 0)
        q1 <- q1 + geom_jitter(size = pt.size)
        q2 <- q2 + geom_jitter(size = pt.size)
        q3 <- q3 + geom_jitter(size = pt.size, height = 0)
    }
    
    p <- cowplot::plot_grid(p1, p2, p3, ncol = 3, align = "h")
    q <- cowplot::plot_grid(q1, q2, q3, ncol = 3, align = "h")
    
    if (plot.show == TRUE) {
        print(p)
        print(q)
    }
    
    if (plot.save == TRUE) {
        
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }

        dir <- file.path(dir, "quality")
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
        
        png(file.path(img_dir, paste(label, "raw_matrix_quality_vlnplot.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(p)
        dev.off()

        png(file.path(img_dir, paste(label, "QC_matrix_quality_vlnplot.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(q)
        dev.off()


        pdf(file = file.path(pdf_dir, paste(label, "QC_matrix_quality_vlnplot.pdf", sep = "")), width = width, height = height)
        print(q)
        dev.off()
        pdf(file = file.path(pdf_dir, paste(label, "raw_matrix_quality_vlnplot.pdf", sep = "")), width = width, height = height)
        print(p)
        dev.off()
    }
    

    return(perturb_QC)
}
