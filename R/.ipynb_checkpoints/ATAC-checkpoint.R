#' Add Meta Data Information onto a SeuratObject
#'
#' Add sgRNA information, fractions of reads in peaks (FRiP) and replicate information onto a SeuratObject for scATAC-seq based input.
#'
#' @param sg_lib Data frame or directory to a txt file containing 3 columns: cell, barcode, gene. If sgRNA information stored in a matrix-like format or input data frame only has sgRNA frequency of each cell, use \code{\link[SCREE]{sgRNAassign}} to assign sgRNA to each cell.
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param fragments Directory of fragments file.
#' @param replicate Directory of a txt file or a vector only containg the replicate information of each cell, in the same order of cells in SeuratObject. If no replicate information, we will consider that all cells are from the same replicate, and this parameter will be set as 1. Default is 1.
#' @param cal.FRiP Logical, calculate FRiP or not. Default is \code{TRUE}.
#'
#' @import Signac
#' @export

ATAC_Add_meta_data <- function(sg_lib, mtx, fragments, replicate = 1, cal.FRiP = TRUE){
    
    #Add "sgRNA_num" and "perturbations"
    
    peak <- Add_meta_data(sg_lib, mtx, cal.mt = FALSE, replicate = replicate)

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

#' Quality Control for scATAC-seq Based Input
#'
#' Perform quality control for scATAC-seq based input based on peaks fraction, nFeature, nCount, FRiP and sgRNA information.
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param chromation.assay Logical, is the assay in the SeuratObject a ChromatinAssay or not. Default is \code{FALSE}
#' @param peak_frac A paramter for filtering low accessibility peaks. By default, only peaks that have counts in at least that fractions of cells are kept. Default is 0.01.
#' @param nFeature Limitation of detected feature numbers in each cell, in the format like c(200, 500000). Default is c(200, 500000). If you don't need an upper limit, you can set the upper limit to an extremely large number. 
#' @param nCount Minimal count numbers in each cell. Default is 1000.
#' @param FRiP Minimal FRiP of each cell. Default is 0.1.
#' @param blank_NTC Logical, use blank control as negative control or not. Default is \code{FALSE}.
#' @param title.size Numeric, title size of the violin plot. Default is 20.
#' @param x.text.size Numeric, x-axis text size of the violin plot. Default is 15.
#' @param x.title.size Numeric, x-axis title size of the violin plot. Default is 15.
#' @param y.text.size Numeric, y-axis text size of the violin plot. Default is 20.
#' @param pt.size Point size of the violin plot. Default is 0.1.
#' @param plot.show Logical, show the Violin plot or not. Default is \code{TRUE}.
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param raster Logical, convert points to raster format, will be useful to reduce the storage cost of the output figure or pdf. Default is \code{FALSE}.
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the output file in inches, for both png and pdf format. Default is 8.3.
#' @param height Height of the graphics region of the output file in inches, for both png and pdf format. Default is 8.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#'
#' @import ggplot2
#' @import Seurat
#' @import ggrastr
#' @import cowplot
#' @export

ATAC_scQC <- function (mtx, chromation.assay = FALSE, peak_frac = 0.01, nFeature = c(200, 500000), nCount = 1000, FRiP = 0.1, blank_NTC = FALSE, title.size = 20, x.text.size = 15, x.title.size = 15, y.text.size = 20, pt.size = 0.1, plot.show = FALSE, plot.save = TRUE, raster = FALSE, prefix = ".", label = "", width = 8.3, height = 8, png_res = 720) {
    
    #read file
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        perturb <- readRDS(mtx)
    } else {
        perturb <- mtx
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
    
    p1 <- ggplot(data = a, mapping = aes(x = Identity, y = nFeature_peak)) + 
    geom_violin(aes(fill = factor(Identity))) + 
    theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "nFeature_peak") + 
    custom_theme

    p2 <- ggplot(data = a, mapping = aes(x = Identity, y = nCount_peak)) + 
    geom_violin(aes(fill=factor(Identity))) + 
    theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "nCount_peak") + 
    custom_theme
    
    p3 <- ggplot(data = a, mapping = aes(x = Identity, y = FRiP)) + 
    geom_violin(aes(fill=factor(Identity))) + 
    theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "FRiP") + 
    custom_theme

    #filter cells with low quality
    
    if (blank_NTC == TRUE) {
        perturb_QC <- subset(perturb,
                             nFeature_peak <= nFeature[2] &
                             nFeature_peak >= nFeature[1] &
                             nCount_peak >= nCount &
                             FRiP >= FRiP)
    } else {
        perturb_QC <- subset(perturb,
                             nFeature_peak <= nFeature[2] &
                             nFeature_peak >= nFeature[1] &
                             nCount_peak >= nCount &
                             FRiP >= FRiP &
                             perturbations != 'blank')
    }
    
    if (chromatin.assay == TRUE) {
        perturb_QC <- CreateSeuratObject(counts = GetAssay(object = perturb_QC), 
                                         assay = "peak", 
                                         min.cells = peak_frac * ncol(perturb_QC), 
                                         project = perturb@project.name)
    } else {
        perturb_QC <- CreateSeuratObject(counts = GetAssayData(object = perturb_QC, slot = "counts"), 
                                         assay = "peak",
                                         min.cells = peak_frac * ncol(perturb_QC), 
                                         project = perturb@project.name)
    }
    
    for(meta in colnames(perturb@meta.data)) {
        if (!(meta %in% c("nFeature_peak", "nCount_peak"))) {
            perturb_QC[[meta]] <- perturb[[meta]]
        }
    }
    
    a <- perturb@meta.data
    a$Identity <- perturb@active.ident
    q1 <- ggplot(data = a, mapping = aes(x = Identity, y = nFeature_peak)) + 
    geom_violin(aes(fill = factor(Identity))) + 
    theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "nFeature_peak") + 
    custom_theme
    
    q2 <- ggplot(data = a, mapping = aes(x = Identity, y = nCount_peak)) + 
    geom_violin(aes(fill=factor(Identity))) + 
    theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "nCount_peak") + 
    custom_theme
    
    q3 <- ggplot(data = a, mapping = aes(x = Identity, y = FRiP)) + 
    geom_violin(aes(fill=factor(Identity))) + 
    theme_classic() + 
    NoLegend() + 
    labs(x = "Identity",y = NULL, title = "FRiP") + 
    custom_theme
    
    if (raster == TRUE) {
        p1 <- p1 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        p2 <- p2 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        p3 <- p3 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        q1 <- q1 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        q2 <- q2 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
        q3 <- q3 + ggrastr::geom_jitter_rast(size = pt.size, raster.dpi = getOption("ggrastr.default.dpi", 300))
    } else {
        p1 <- p1 + geom_jitter(size = pt.size)
        p2 <- p2 + geom_jitter(size = pt.size)
        p3 <- p3 + geom_jitter(size = pt.size)
        q1 <- q1 + geom_jitter(size = pt.size)
        q2 <- q2 + geom_jitter(size = pt.size)
        q3 <- q3 + geom_jitter(size = pt.size)
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
        
        pdf(file = file.path(dir, paste(label, "raw_matrix_quality_vlnplot.pdf", sep = "")), 
            width = width, height = height)
        print(p)
        dev.off()
    
        png(file.path(img_dir, paste(label, "raw_matrix_quality_vlnplot.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(p)
        dev.off()
        
        pdf(file = file.path(dir, paste(label, "QC_matrix_quality_vlnplot.pdf", sep = "")), 
            width = width, height = height)
        print(q)
        dev.off()

        png(file.path(img_dir, paste(label, "QC_matrix_quality_vlnplot.png", sep = "")), 
            width = width, height = height, unit = "in", res = png_res)
        print(q)
        dev.off()
    }
    
    return(perturb_QC)
}

#' Calculate Gene Activity
#'
#' Generate gene activity matrix using peak count matrix.
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param fragments Directory of fragments file or a FragmentObject generated from \code{\link[Signac]{CreateFragmentObject}}.
#' @param species Only support "Hs" and "Mm". Default is "Hs".
#' @param version Version of the reference genome(Ensembl). Default is "v75".
#' @param gene_type Type of gene name, selected from one of c("Symbol", "Ensembl"). Default is "Symbol".
#' @param protein_coding Logical, only use protein coding gene or not. Default is \code{TRUE}.
#' @param pro_up Numeric, the number of nucleotides upstream of the transcription start site that should be included in the promoter region. Default is 3000.
#' @param pro_down Numeric, the number of nucleotides downstream of the transcription start site that should be included in the promoter region. Default is 0.
#' @param sep Vector of separators to use for genomic string. First element is used to separate chromosome and coordinates, second separator is used to separate start and end coordinates.
#'
#' @import Signac
#' @import Seurat
#' @export

CalculateGeneActivity <- function(mtx, fragments, species = "Hs", version = "v75", gene_type = "Symbol", protein_coding = TRUE, pro_up = 3000, pro_down = 0, sep = c("-", "-")){
    
    #get promoter region
    
    pro <- GetPromoter(species, version, gene_type, protein_coding, pro_up, pro_down)
    genebodyandpromoter.coords <- pro[[1]]
    gene.key <- pro[[2]]
    
    #get count matrix
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        perturb <- readRDS(mtx)
    } else {
        perturb <- mtx
    }

    #get fragments
    
    if (is.character(fragments)) {
        fragments <- CreateFragmentObject(fragments, cells = colnames(perturb))
    }
    
    #calculate gene activity
    
    gene.activity <- FeatureMatrix(fragments = fragments, 
                                   features = genebodyandpromoter.coords, 
                                   cells = colnames(perturb), 
                                   sep = sep)
    
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

#' Get Potential Promoter Regions
#'
#' Get potential promoter regions for each gene, using specific reference.
#' 
#' @param species Only support "Hs" and "Mm". Default is "Hs".
#' @param version Version of the reference genome(Ensembl). Default is "v75".
#' @param gene_type Type of gene name, selected from one of c("Symbol", "Ensembl"). Default is "Symbol".
#' @param protein_coding Logical, only use protein coding gene or not. Default is \code{TRUE}.
#' @param pro_up Numeric, the number of nucleotides upstream of the transcription start site that should be included in the promoter region. Default is 3000.
#' @param pro_down Numeric, the number of nucleotides downstream of the transcription start site that should be included in the promoter region. Default is 0.
#'
#' @import EnsDb.Hsapiens.v75
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Mmusculus.v75
#' @import EnsDb.Mmusculus.v79
#' @import ensembldb
#' @import Signac
#' @import GenomeInfoDb
#' @import IRanges
#' @export

GetPromoter <- function(species = "Hs", version = "v75", gene_type = "Symbol", protein_coding = TRUE, pro_up = 3000, pro_down = 0){

    #get gene ranges from selected reference
    if (species == "Hs") {
        if (version == "v75") {
            gene.ranges <- genes(EnsDb.Hsapiens.v75)
        } else if (version == "v79") {
            gene.ranges <- genes(EnsDb.Hsapiens.v79)
        } else if (version == "v86") {
            gene.ranges <- genes(EnsDb.Hsapiens.v86)
        }  
    } else if (species == "Mm") {
        if (version == "v75") {
            gene.ranges <- genes(EnsDb.Mmusculus.v75)
        } else if (version == "v79") {
            gene.ranges <- genes(EnsDb.Mmusculus.v79)
        }
    } else {
        stop("Species must be one of 'Hs' and 'Mm'")
    }

    seqlevelsStyle(gene.ranges) <- 'UCSC'
    
    if (protein_coding == TRUE) {
        gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
    }

    gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

    genebodyandpromoter.coords <- suppressWarnings(trim(Extend(x = gene.ranges, 
                                                               upstream = pro_up, downstream = pro_down)))

    if (gene_type == "Symbol") {
        gene.key <- genebodyandpromoter.coords$gene_name
    } else if (gene_type == "Ensembl") {
        gene.key <- genebodyandpromoter.coords$gene_id
    } else {
        warning("This function only support 'gene Symbol' and 'Ensembl id' as gene names, using gene Symbol instead")
        gene.key <- genebodyandpromoter.coords$gene_name
    }
    

    names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
    
    return(list(genebodyandpromoter.coords, gene.key))
}

