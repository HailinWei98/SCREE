#' function definitions ##### Normalize and scale data.
#' @export

normalize_scale <- function(mtx_dir){
    
    #read file
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        perturb_QC <- readRDS(mtx_dir)
    } else {
        perturb_QC <- mtx_dir
    }

    #normalize and scale on the data
    
    perturb_QC <- NormalizeData(
        object = perturb_QC,
        normalization.method = "LogNormalize",
        scale.factor = 10000)
    if("percent.mt" %in% colnames(perturb_QC@meta.data)){
        perturb_QC <- ScaleData(perturb_QC, features = rownames(perturb_QC),
                                vars.to.regress = c("nCount_RNA", "percent.mt"))
    }else{
        perturb_QC <- ScaleData(perturb_QC, features = rownames(perturb_QC),
                                vars.to.regress = c("nCount_RNA"))
    }
    return(perturb_QC)
}
