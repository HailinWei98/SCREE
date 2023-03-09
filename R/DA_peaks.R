#' Find Differential Peaks
#' 
#' Find differential accessible peaks between selected perturbation and negative control. Based on \code{\link[Seurat]{FindMarkers}}.
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param selected Perturbation to find DA peaks. 
#' @param NTC The name of negative controls. Default is "NTC".
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the selected perturbation or NTC. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1.
#' @param test.use Denotes which test to use. Available options are: "wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2". See details from \code{\link[Seurat]{FindMarkers}}. Default is "wilcox".
#' @param p_adj_cut Maximum adjust p_value to . Default is 0.05.
#' @param logFC_cut Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' 
#' @import Seurat
#' @export

DApeaks <- function(mtx, selected, NTC = "NTC", min.pct = 0.1, test.use = "wilcox", p_adj_cut = 0.05, logFC_cut = 0.25) {
    
    if (is.character(mtx)) {
        peak <- readRDS(mtx)
    } else {
        peak <- mtx
    }
    
    #rename ident
    
    if ("perturbations" %in% colnames(peak@meta.data)) {
        if(selected %in% peak$perturbations){
            peak@active.ident <- peak$perturbations
        } else {
            stop(paste("Cannot find ", selected, " in perturbations", sep = "'"))
        }
    } else {
        stop("Cannot find 'perturbations' in object, please renamed it or using 'Add_meta_data' function to add it")
    }
    
    #find DA peaks
    
    da_peaks <- FindMarkers(
        object = peak,
        ident.1 = selected,
        ident.2 = NTC,
        min.pct = min.pct,
        test.use = test.use,
        logfc.threshold = logFC_cut
    )
    
    if (nrow(da_peaks) == 0) {
        return(NULL)
    }
    
    da_peak <- subset(da_peaks, p_val_adj <= p_adj_cut)
    if (nrow(da_peak) == 0) {
        return(NULL)
    }
    
    da_peak$chromosome <- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,1]
    da_peak$start <- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,2]
    da_peak$end <- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,3]
    
    return(da_peak)
}

#' Find Potential Enhancer List
#'
#' Identify potential enhancer list from DA peaks.
#'
#' @param da_peak Data frame of DA peaks list, generated from \code{\link[SCREE]{DApeaks}}.
#' @param gene_annotations Gene annotations stored in data frame format, including c("chromosome", "start", "end", "strand", "transcript") as colnames. Default is \code{NULL}, gene annotations are from \code{\link{ensembldb}}.
#' @param overlap_cut Maximum overlap nucleotides between peaks and promoters. Default is 0.
#' @param pro_up The number of nucleotides upstream of the transcription start site that should be included in the promoter region. Default is 3000.
#' @param pro_down The number of nucleotides downstream of the transcription start site that should be included in the promoter region. Default is 0.
#'
#' @export

enhancer <- function(da_peak, gene_anno, overlap_cut, pro_up = 2000, pro_down = 0) {
    
    #get enhancer list
        
    enhancer_list <- data.frame()
    for(i in 1:nrow(da_peak)){
        chr <- da_peak[i, "chromosome"]
        start <- as.numeric(da_peak[i, "start"])
        end <- as.numeric(da_peak[i, "end"])
        # minbp <- start - 2 * (pro_up + pro_down)
        # maxbp <- end + 2 * (pro_up + pro_down)
        minbp <- start
        maxbp <- end
        peak_model <- gene_anno
        peak_model <- peak_model[!is.na(peak_model$chromosome) & 
                                 !is.na(peak_model$start) &
                                 !is.na(peak_model$end) &
                                 !is.na(peak_model$strand), ]
        peak_model$start <- peak_model$start + pro_up
        peak_model$end <- peak_model$start - pro_down
        peak_model <- peak_model[peak_model$chromosome == chr &
                                 ((peak_model$start > minbp & peak_model$start < maxbp) |
                                  (peak_model$end > minbp & peak_model$end < maxbp) |
                                  (peak_model$start < minbp & peak_model$end > maxbp)), ]
        if (nrow(peak_model) == 0) {
            enhancer_list <- rbind(enhancer_list, c(chr, start, end))
        } else {
            pro_seq <- apply(peak_model[,c("start", "end")], 1, function(x){seq(x[1], x[2])})
            overlap <- length(intersect(seq(start, end), pro_seq))
            if (overlap <= overlap_cut) {
                enhancer_list <- rbind(enhancer_list, c(chr, start, end))
            }
        }
    }
    if (nrow(enhancer_list != 0)) {
        colnames(enhancer_list) <- c("chromosome", "start", "end")
        rownames(enhancer_list) <- paste(enhancer_list$chromosome, 
                                         enhancer_list$start, 
                                         enhancer_list$end, 
                                         sep = "-")
    } else {
        enhancer_list <- da_peak[, c("chromosome", "start", "end")]
    }
    
    return(enhancer_list)
}