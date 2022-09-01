#' @import Signac
#' @export

ATAC_Add_meta_data <- function(sg_dir, mtx_dir, fragments, replicate = 1, cal.FRiP = TRUE){
    
    #Add "sgRNA_num" and "perturbations"
    
    peak <- Add_meta_data(sg_dir, mtx_dir, cal.mt = FALSE, replicate = replicate)

    #calculate FRiP
    
    if(cal.FRiP == TRUE){
        if(is.character(fragments)){
            frag <- CountFragments(fragments = fragments, cells = colnames(peak))
            peak$fragments <- frag$reads_count
            peak <- FRiP(peak, "peaks", total.fragments = "fragments")
        }else{
            stop("Please provide the path of fragments file")
        }
    }else{
        if(!("FRiP" %in% colnames(peak@meta.data))){
            warning("Please make sure that there is meta data named 'FRiP' in input peak matrix while setting 'cal.FRiP = FALSE'. We will set the 'FRiP' of each cell to 1.")
            peak$FRiP <- 1
        }
    }
    
    return(peak)
}

#' @export

ATAC_scQC <- function(mtx_dir, prefix = ".", label = "", peak_frac = 0.01, 
                      nFeature = c(200, 500000), nCount = 1000, FRiP = 0.1, blank_NTC = FALSE){
    
    dir <- file.path(prefix, "results")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    dir <- file.path(dir, "ATAC_quality")
    if (!(dir.exists(dir))){
        dir.create(dir)
    }
    
    img_dir <- file.path(dir, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    dir <- file.path(dir, "pdf")
    if (!(dir.exists(dir))){
        dir.create(dir)
    }
    
    #read file
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        perturb <- readRDS(mtx_dir)
    } else {
        perturb <- mtx_dir
    }
    
    if(!("replicate" %in% colnames(perturb@meta.data) & "perturbations" %in% colnames(perturb@meta.data))){
        stop("Cannot find meta data names 'replicate' or 'perturbations' in input matrix.")
    }
    
    if(!("FRiP" %in% colnames(perturb@meta.data))){
        warning("No meta data named 'FRiP' in input peak matrix. We will set the 'FRiP' of each cell to 1.")   
    }
    perturb$nFeature_peak <- perturb[[paste("nFeature_", perturb@active.assay, sep = "")]][, 1]
    perturb$nCount_peak <- perturb[[paste("nCount_", perturb@active.assay, sep = "")]][, 1]

    #QC plot of the single cell matrix
    
    pdf(file = file.path(dir, paste(label, "raw_matrix_quality_vlnplot.pdf", sep = "")), width = 8, height = 8)
    p1 <- VlnPlot(perturb, features = "nFeature_peak", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) +
    NoLegend()
    p2 <- VlnPlot(perturb, features = "nCount_peak", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) + 
    NoLegend()
    p3 <- VlnPlot(perturb, features = "FRiP", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) + 
    NoLegend()
    p <- plot_grid(p1, p2, p3, ncol = 3, align = "h")
    print(p)
    dev.off()
    
    png(file.path(img_dir, paste(label, "raw_matrix_quality_vlnplot.png", sep = "")), 
        width = 800, height = 600)
    print(p)
    dev.off()

    
    #filter cells with low quality
    
    if(blank_NTC == TRUE){
        perturb_QC <- subset(perturb,
                             nFeature_peak <= nFeature[2] &
                             nFeature_peak >= nFeature[1] &
                             nCount_peak >= nCount &
                             FRiP >= FRiP)
    }else{
        perturb_QC <- subset(perturb,
                             nFeature_peak <= nFeature[2] &
                             nFeature_peak >= nFeature[1] &
                             nCount_peak >= nCount &
                             FRiP >= FRiP &
                             perturbations != 'blank')
    }
    perturb_QC <- CreateSeuratObject(counts = GetAssayData(object = perturb_QC, slot = "counts"), assay = "peak",
                                     min.cells = peak_frac * ncol(perturb_QC), project = perturb@project.name)
    
    for(meta in colnames(perturb@meta.data)) {
        if (!(meta %in% c("nFeature_peak", "nCount_peak"))) {
            perturb_QC[[meta]] <- perturb[[meta]]
        }
    }
    
    pdf(file = file.path(dir, paste(label, "QC_matrix_quality_vlnplot.pdf", sep = "")), width = 8, height = 8)
    q1 <- VlnPlot(perturb_QC, features = "nFeature_peak", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) +
    NoLegend()
    q2 <- VlnPlot(perturb_QC, features = "nCount_peak", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) + 
    NoLegend()
    q3 <- VlnPlot(perturb_QC, features = "FRiP", pt.size = 0.1) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(hjust = 0.5, size = 16)) + 
    NoLegend()
    q <- plot_grid(q1, q2, q3, ncol = 3, align = "h")
    print(q)
    dev.off()
    
    png(file.path(img_dir, paste(label, "QC_matrix_quality_vlnplot.png", sep = "")), 
        width = 800, height = 600)
    print(q)
    dev.off()
    return(perturb_QC)
}

#' @export

CalculateGeneActivity <- function(mtx_dir, fragments, species = "Hs", version = "v75", 
                                  gene_type = "Symbol", protein_coding = TRUE, pro_up = 3000, pro_down = 0){
    
    #get promoter region
    
    pro <- GetPromoter(species, version, gene_type, protein_coding, pro_up, pro_down)
    genebodyandpromoter.coords <- pro[[1]]
    gene.key <- pro[[2]]
    
    #get count matrix
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        perturb <- readRDS(mtx_dir)
    } else {
        perturb <- mtx_dir
    }

    #get fragments
    
    if(is.character(fragments)){
        fragments <- CreateFragmentObject(fragments, cells = colnames(perturb))
    }
    
    #calculate gene activity
    
    gene.activity <- FeatureMatrix(fragments = fragments, 
                                   features = genebodyandpromoter.coords, cells = colnames(perturb))
    
    #generate gene activity matrix
    
    rownames(gene.activity) <- gene.key[rownames(gene.activity)]
    perturb_RNA <- CreateSeuratObject(counts = gene.activity, project = perturb@project.name)
    if("replicate" %in% colnames(perturb@meta.data)){
        replicate <- perturb$replicate
        perturb_RNA$replicate <- replicate
    }
    
    if("perturbations" %in% colnames(perturb@meta.data)){
        perturb_RNA$perturbations <- perturb$perturbations
    }
    return(perturb_RNA)
}

#' @export

GetPromoter <- function(species = "Hs", version = "v75", gene_type = "Symbol", 
                        protein_coding = TRUE, pro_up = 3000, pro_down = 0){

    #get gene ranges from selected reference
    if(species == "Hs"){
        if(version == "v75"){
            gene.ranges <- genes(EnsDb.Hsapiens.v75)
        }else if(version == "v79"){
            gene.ranges <- genes(EnsDb.Hsapiens.v79)
        }else if(version == "v86"){
            gene.ranges <- genes(EnsDb.Hsapiens.v86)
        }  
    }else if(species == "Mm"){
        if(version == "v75"){
            gene.ranges <- genes(EnsDb.Mmusculus.v75)
        }else if(version == "v79"){
            gene.ranges <- genes(EnsDb.Mmusculus.v79)
        }
    }

    seqlevelsStyle(gene.ranges) <- 'UCSC'
    
    if(protein_coding == TRUE){
        gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
    }

    gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

    genebodyandpromoter.coords <- suppressWarnings(trim(Extend(x = gene.ranges, 
                                                               upstream = pro_up, downstream = pro_down)))

    if(gene_type == "Symbol"){
        gene.key <- genebodyandpromoter.coords$gene_name
    }else if(gene_type == "Ensembl"){
        gene.key <- genebodyandpromoter.coords$gene_id
    }else{
        warning("This function only support 'gene Symbol' and 'Ensembl id' as gene names, using gene Symbol instead")
        gene.key <- genebodyandpromoter.coords$gene_name
    }
    

    names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
    
    return(list(genebodyandpromoter.coords, gene.key))
}

