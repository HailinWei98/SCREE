#' Batch Run
#'
#' Get all results using one function easily.
#' 
#' @param sg_dir Data frame or directory to a txt file containing 3 columns: cell, barcode, gene. If sgRNA information stored in a matrix-like format or sinput data frame only has sgRNA frequence of each cell, use \code{\link[SCREEN]{sgRNAassign}} to assign sgRNA to each cell.
#' @param mtx_dir SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param fragments Directory of fragments file used to calculate FRiP for perturb-ATAC input.
#' @param cal.FRiP Logical, calculate FRiP or not. Default is \code{TRUE}.
#' @param species Only support "Hs" and "Mm", if input other species, \code{percent.mt} will be count as "Mm". Default is "Hs".
#' @param version Version of the reference genome(Ensembl), used for perturb-ATAC input and perturb-enhancer input. Default is "v75".
#' @param data_type Type of input data, can be one of c("RNA", "ATAC"). Default is "RNA".
#' @param Mixscape Logical, run \code{IntegratedMixscape} or not. Default is \code{TRUE}.
#' @param prefix Path to save all the results. Default is current directory.
#' @param label The label of the output file.
#' @param gene_type Type of gene name, selected from one of c("Symbol", "Ensembl"). Default is "Symbol".
#' @param protein_coding Logical, only use protein coding gene or not. This parameter is only used for calculating gene activity for perturb-ATAC input. Default is \code{TRUE}.
#' @param frac A paramter for filtering low expressed genes or low accessibility peaks. By default, only genes or peaks that have expressions or counts in at least that fractions of cells are kept. Default is 0.01.
#' @param cal.mt Logical, calculate percentage of mitochondrial gene expression of each cell or not. Default is \code{TRUE}.
#' @param nFeature Limitation of detected feature numbers in each cell, in the format c(200, 5000). Default is c(200, 5000).
#' @param nCount Minimal count numbers in each cell. Default is 1000.
#' @param FRiP Minimal FRiP of each cell. Default is 0.1.
#' @param mt Maximal percentage of mitochondrial gene expression of each cell. Default is 10.
#' @param blank_NTC Logical, use blank control as negative control or not. Default is \code{FALSE}.
#' @param lambda Parameter used in ridge regression of \code{improved_scmageck_lr}. Default is 0.01.
#' @param permutation Permutation times in \code{improved_scmageck_lr}. Default is 10000.
#' @param p_val_cut P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.
#' @param score_cut Score cutoff of \code{improved_scmageck_lr} results. Default is 0.5.
#' @param cicero_p_val_cut P-value cutoff of \code{improved_scmageck_lr} results used for \code{ciceroPlot}. Default is 0.05.
#' @param cicero_socre_cut Score cutoff of \code{improved_scmageck_lr} results used for \code{ciceroPlot}. Default is 0.
#' @param ylimit Limitation of y-axis of \code{DE_gene_plot} in the format c(-600, 600, 200). These numbers mean c(minimum, maximum, interval). Default is "auto", which means that this function will get \code{ylimit} automatically.
#' @param project Title of \code{DE_gene_plot}. Default is "perturb".
#' @param NTC The name of the genes served as negative controls. Default is "NTC".
#' @param replicate Required a vector of replicate information corresponding to each cell with the same order. Default is 1, which means no replicate. 
#' @param select_gene The list of genes for regression in \code{scMAGeCK} step. By default, all genes in the table are subject to regression.
#' @param selected Enhancer regions to visualize for perturb-enhancer or perturbations to chose for perturb-ATAC, in \code{cicero} step. By default, all enhancers or all perturbations will be chosen.
#' @param gene_annotations Gene annotations stored in data frame format, including c("chromosome", "start", "end", "strand", "transcript") as colnames, used for /code{ciceroPlot} step. By default, gene annotations are from \code{ensembldb}.
#' @param pro_annotations Gene annotations stored in data frame format, including c("chromosome", "start", "end", "strand", "transcript") as colnames. By default, gene annotations are from \code{ensembldb}.
#' @param pro_up The number of nucleotides upstream of the transcription start site that should be included in the promoter region, only used for perturb-ATAC data. Default is 3000.
#' @param pro_down The number of nucleotides downstream of the transcription start site that should be included in the promoter region, only used for perturb-ATAC data. Default is 0.
#' @param overlap_cut Maximum overlap nucleotides between peaks and promoters, only used for perturb-ATAC data. Default is 0.
#' @param p_adj_cut Parameter only used for finding DA peaks. Maximum adjust p_value calculated by \code{\link[Seurat]{FindMarkers}}. Default is 0.05.
#' @param logFC_cut Parameter only used for finding DA peaks. Minimum log fold change calculated by \code{\link[Seurat]{FindMarkers}}. Default is 1.
#' @param min.pct Parameter only used for finding DA peaks. Only test genes that are detected in a minimum fraction of min.pct cells in either of the NTC or perturbations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.2. 
#' @param upstream The number of nucleotides upstream of the start site of selected region in \code{ciceroPlot} step. Default is 2000000.
#' @param downstream The number of nucleotides downstream of the start site of selected region in \code{ciceroPlot} step. Default is 2000000.
#' @param test.use Parameter only used for finding DA peaks. Default is "wilcox".For more details, see \code{\link[Seurat]{FindMarkers}}.
#' @param track_size Size of each axis in /code{ciceroPlot} result. Default is c(1,.3,.2,.3). If `include_axis_track=FALSE`, track_size should be a vector with 3 elements.
#' @param include_axis_track Logical, should a genomic axis be plotted? Default is \code{TRUE}.
#' @param connection_color Color for connection lines. A single color, the name of a column containing color values, or the name of a column containing a character or factor to base connection colors on.
#' @param connection_color_legend Logical, should connection color legend be shown?
#' @param connection_width Width of connection lines.
#' @param connection_ymax Connection y-axis height. If NULL, chosen automatically.
#' @param gene_model_color Color for gene annotations.
#' @param alpha_by_coaccess Logical, should the transparency of connection lines be scaled based on co-accessibility score?

#' @import Seurat
#' @import ggplot2
#' @import plyr
#' @import EnsDb.Hsapiens.v75
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Mmusculus.v75
#' @import EnsDb.Mmusculus.v79
#' @import Signac
#' @import Gviz
#' @import ggplotify
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import cowplot
#' @export

SCREEN <- function(sg_dir, mtx_dir, fragments, cal.FRiP = TRUE, species = "Hs", version = "v75", 
                   data_type = "RNA", Mixscape = TRUE, prefix = "./", label = "", gene_type = "Symbol", 
                   protein_coding = TRUE, frac = 0.01, cal.mt = TRUE, nFeature = c(200, 5000), nCount = 1000, 
                   FRiP = 0.1, mt = 10, blank_NTC = FALSE, lambda = 0.01, permutation = NULL,
                   p_val_cut = 0.05, score_cut = 0.5, cicero_p_val_cut = 0.05, cicero_score_cut = 0,
                   ylimit = "auto", project = "perturb", NTC = "NTC", replicate = 1, 
                   select_gene = NULL, selected =  NULL, gene_annotations = NULL, pro_annotations = NULL, 
                   pro_up = 3000, pro_down = 0, overlap_cut = 0, p_adj_cut = 0.05, logFC_cut = 1,
                   min.pct = 0.2, upstream = 2000000, downstream = 2000000, test.use = "wilcox", 
                   track_size = c(1,.3,.2,.3), include_axis_track = TRUE, connection_color = "#7F7CAF",
                   connection_color_legend = TRUE, connection_width = 2, connection_ymax = NULL,
                   gene_model_color = "#81D2C7", alpha_by_coaccess = FALSE,
                   gene_model_shape = c("smallArrow", "box")) {
    
    #get matrix
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        mtx <- readRDS(mtx_dir)
    } else {
        mtx <- mtx_dir
    }
    
    #get sgRNA library
    
    if (is.character(sg_dir)) {
        message(paste("Reading sgRNA lib file:", sg_dir))
        sg_lib <- read.table(sg_dir, header =T )
    } else {
        sg_lib <- sg_dir
    }
    
    #for RNA input
    
    if (data_type == "RNA") {
        mtx <- Add_meta_data(sg_lib, mtx, cal.mt, species, replicate)
        mtx_QC <- scQC(mtx, prefix, label, species, frac, nFeature, nCount, mt, blank_NTC)
        p <- sgRNA_quality_plot(sg_lib, mtx,
                               LABEL = paste(label, "before_QC", sep = ""), prefix)
        q <- sgRNA_quality_plot(sg_lib, mtx_QC,
                                LABEL = paste(label, "after_QC", sep = ""), prefix)
        if(Mixscape == TRUE){
            message("Running Mixscape")
            mixscape <- IntegratedMixscape(sg_lib, mtx_QC, NTC, prefix, label)
        }
        mtx_QC <- normalize_scale(mtx_QC)
        scaled <- GetAssayData(object = mtx_QC, slot = "scale.data")
        
        message("Save RDS file")
        saveRDS(mtx, file = file.path(prefix, "perturb.rds"))
        saveRDS(mtx_QC, file = file.path(prefix, "perturb_QC.rds"))
        
        message("Save sgRNA library to save path")
        write.table(sg_lib, file = file.path(prefix, "sg_lib.txt"),
                    row.names = FALSE, quote = FALSE)
        
        message("Running scMAGeCK")
        result <- improved_scmageck_lr(sg_lib, scaled, NTC, select_gene,
                                        LABEL = paste(label, "improved",sep = ""),
                                        permutation, prefix, lambda)
        score <- result[[1]]
        p_val <- result[[2]]
        
        message("Running DE_gene_plot")
        y <- DE_gene_plot(score, p_val, project,
                          prefix, label, p_val_cut, score_cut,
                          ylimit)
        if(length(grep("^chr", sg_lib$gene)) != 0){
            message("Detect enhancer perturbation in sgRNA library, running ciceroPlot")
            cicero_result <- ciceroPlot(score, p_val, selected, species, versions, gene_annotations, 
                                        cicero_p_val_cut, cicero_score_cut, upstream, downstream, 
                                        track_size, include_axis_track, connection_color,
                                        connection_color_legend, connection_width, connection_ymax, 
                                        gene_model_color, alpha_by_coaccess, gene_model_shape)
        } 
        results <- list()
        results[1] <- result
        results[2] <- y
        results[3] <- cicero_result
        names(results) <- c("scMAGeCK_lr", "DE_gene_plot", "cicero_results")
    } else if(data_type == "ATAC"){
        mtx <- ATAC_Add_meta_data(sg_lib, mtx, fragments, replicate, cal.FRiP)
        mtx_QC <- ATAC_scQC(mtx, prefix, label, frac, nFeature, nCount, FRiP, blank_NTC)
        p <- sgRNA_quality_plot(sg_lib, mtx,
                               LABEL = paste(label, "before_QC", sep = ""), prefix)
        q <- sgRNA_quality_plot(sg_lib, mtx_QC,
                                LABEL = paste(label, "after_QC", sep = ""), prefix)
        
        message("Calculate gene activity")
        RNA <- CalculateGeneActivity(mtx, fragments, species, version, gene_type, protein_coding, pro_up, pro_down)
        
        if(Mixscape == TRUE){
            message("Running Mixscape")
            mixscape <- IntegratedMixscape(sg_lib, RNA, NTC, prefix, label)
        }
        RNA <- normalize_scale(RNA)
        scaled <- GetAssayData(object = RNA, slot = "scale.data")
        
        message("Save RDS file")
        saveRDS(mtx, file = file.path(prefix, "peak.rds"))
        saveRDS(RNA, file = file.path(prefix, "RNA.rds"))
        
        message("Save sgRNA library to save path")
        write.table(sg_lib, file = file.path(prefix, "sg_lib.txt"),
                    row.names = FALSE, quote = FALSE)
        
        message("Running scMAGeCK")
        result <- improved_scmageck_lr(sg_lib, scaled, NTC, select_gene,
                                       LABEL = paste(label, "improved", sep = ""),
                                       permutation, prefix, lambda)
        score <- result[1]
        p_val <- result[2]
        
        message("Running DE_gene_plot of gene activity matrix")
        y <- DE_gene_plot(score, p_val, project,
                          prefix, label, p_val_cut, score_cut,
                          ylimit)
        message("Detect enhancer perturbation in sgRNA library, running ATACciceroPlot")
        cicero_result <- ATACciceroPlot(mtx_QC, score, p_val, selected, species, versions, gene_annotations, 
                                        pro_annotations, pro_up, pro_down, overlap_cut, cicero_p_val_cut, 
                                        cicero_score_cut, p_adj_cut, logFC_cut, NTC, min.pctupstream, downstream, 
                                        test.use, track_size, include_axis_track, connection_color,
                                        connection_color_legend, connection_width, connection_ymax, 
                                        gene_model_color, alpha_by_coaccess, gene_model_shape)
        
        results <- list()
        results[1] <- result
        results[2] <- y
        results[3] <- cicero_result
        names(results) <- c("scMAGeCK_lr", "DE_gene_plot", "cicero_results")
    }
    return(results)
}
