#' Normalize and Scale Data
#'
#' Normalize and scale the raw UMI count matrix based on \code{\link{Seurat}}.
#'
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.  
#' @param normalization.method Method for normalization, can be one of "LogNormalize", "CLR", "RC". See more details from \code{\link[Seurat]{NormalizeData}}. Default is "LogNormalize".
#' @param scale.factor Sets the scale factor for cell-level normalization. Default is 10000.
#' @param vars.to.regress Variables to regress out (previously latent.vars in RegressOut). Default is c("nCount_RNA", "percent.mt").
#'
#' @import Seurat
#' @export

normalize_scale <- function(mtx, normalization.method = "LogNormalize", scale.factor = 10000, vars.to.regress = c("nCount_RNA", "percent.mt"), features = NULL){
    
    #read file
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx))
        perturb_QC <- readRDS(mtx)
    } else {
        perturb_QC <- mtx
    }

    #normalize and scale on the data
    
    perturb_QC <- NormalizeData(
        object = perturb_QC,
        normalization.method = normalization.method,
        scale.factor = scale.factor)
    if (is.null(features)) {
        features = rownames(perturb_QC)
    }
    
    perturb_QC <- ScaleData(perturb_QC, features = features,
                            vars.to.regress = vars.to.regress)
    return(perturb_QC)
    
}
