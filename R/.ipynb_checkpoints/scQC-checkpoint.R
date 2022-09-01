#' function definitions ##### Common single-cell quality control.
#' This step will only save the VlnPlot before and after QC but not show.
#' Users can do this step by themselves using Seurat.
#' @import cowplot
#' @export

scQC<- function(mtx_dir, prefix = "./", label = "", species = "Hs", gene_frac = 0.01,
                nFeature = c(200, 5000), nCount = 1000, mt = 10, blank_NTC = FALSE){
    
    dir <- file.path(prefix, "pdf")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    dir <- file.path(dir, "quality")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "quality")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }

    #read file
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        perturb <- readRDS(mtx_dir)
    } else {
        perturb <- mtx_dir
    }

    if(!("replicate" %in% colnames(perturb@meta.data) & "perturbations" %in% colnames(perturb@meta.data))){
        stop("Cannot find meta data named 'replicate' or 'perturbations' in input matrix.")
    }
    if(!("percent.mt" %in% colnames(perturb@meta.data))){
        warning("No meta data named 'percent.mt', we will set the 'percent.mt' of each cell to 0")
        perturb$percent.mt <- 0
    }
    perturb$nFeature_RNA <- perturb[[paste("nFeature_", perturb@active.assay, sep = "")]][, 1]
    perturb$nCount_RNA <- perturb[[paste("nCount_", perturb@active.assay, sep = "")]][, 1]

    #filter cells with low quality
        
    if(blank_NTC == TRUE){
        perturb_QC <- subset(perturb,
                             nFeature_RNA <= nFeature[2] &
                             nFeature_RNA >= nFeature[1] &
                             nCount_RNA >= nCount &
                             percent.mt <= mt)
    }else{
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
    
    if(species == "Hs"){
        perturb_QC[["percent.mt"]] <- PercentageFeatureSet(perturb_QC, pattern = "^MT-")
    }else{
        perturb_QC[["percent.mt"]] <- PercentageFeatureSet(perturb_QC, pattern = "^mt-")
    }
    
    for(meta in colnames(perturb@meta.data)) {
        if (!(meta %in% c("nFeature_RNA", "nCount_RNA", "percent.mt"))) {
            perturb_QC[[meta]] <- perturb[[meta]]
        }
    }
    
    perturb_QC@active.ident <- as.factor(perturb$orig.ident)
    
    #QC plot of the single cell matrix
    
    p1 <- VlnPlot(perturb, features = "nFeature_RNA", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.text.y = element_text(hjust = 0.5, size = 20)) +
    NoLegend()
    p2 <- VlnPlot(perturb, features = "nCount_RNA", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.text.y = element_text(hjust = 0.5, size = 20)) + 
    NoLegend()
    p3 <- VlnPlot(perturb, features = "percent.mt", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(hjust = 0.5, size = 20)) + 
    NoLegend()
    p <- plot_grid(p1, p2, p3, ncol = 3, align = "h")
    

    png(file.path(img_dir, paste(label, "raw_matrix_quality_vlnplot.png", sep = "")), 
        width = 800, height = 600)
    print(p)
    dev.off()
    
    q1 <- VlnPlot(perturb_QC, features = "nFeature_RNA", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) +
    NoLegend()
    q2 <- VlnPlot(perturb_QC, features = "nCount_RNA", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) + 
    NoLegend()
    q3 <- VlnPlot(perturb_QC, features = "percent.mt", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) + 
    NoLegend()
    q <- plot_grid(q1, q2, q3, ncol = 3, align = "h")

    png(file.path(img_dir, paste(label, "QC_matrix_quality_vlnplot.png", sep = "")), 
        width = 800, height = 600)
    print(q)
    dev.off()
    
    
    pdf(file = file.path(dir, paste(label, "QC_matrix_quality_vlnplot.pdf", sep = "")), width = 8.3, height = 8)
    print(q)
    dev.off()
    pdf(file = file.path(dir, paste(label, "raw_matrix_quality_vlnplot.pdf", sep = "")), width = 8.3, height = 8)
    print(p)
    dev.off()
    return(perturb_QC)
}
