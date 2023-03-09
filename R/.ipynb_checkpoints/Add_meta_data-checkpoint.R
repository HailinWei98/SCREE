#' Add Metadata Information onto a SeuratObject
#'
#' Add sgRNA information, mitochondrial gene percentage and replicate information onto a SeuratObject.
#'
#' @param sg_lib Data frame or directory to a txt file containing 3 columns: cell, barcode, gene. If sgRNA information stored in a matrix-like format or input data frame only has sgRNA frequency of each cell, use \code{\link[SCREE]{sgRNAassign}} to assign sgRNA to each cell.
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param cal.mt Logical, calculate mitochondrial gene percentage or not. Default is \code{TRUE}.
#' @param species Only support "Hs" and "Mm", if input other species, \code{percent.mt} will be count as "Mm". Default is "Hs".
#' @param replicate Directory of a txt file or a vector only contains the replicate information of each cell, in the same order of cells in SeuratObject. If no replicate information, we will consider that all cells are from the same replicate, and this parameter will be set as 1. Default is NULL.
#' 
#' @importFrom utils read.table write.table
#' @importFrom plyr count
#' @import Seurat
#' 
#' @export

Add_meta_data <- function(sg_lib, mtx, cal.mt = TRUE, species = "Hs", replicate = NULL) {
    
    #read files
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        mtx <- readRDS(mtx)
    }
    
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        sg_lib <- read.table(sg_lib, header = T)
    }

    #remove cells in sgRNA library that are not included in matrix
    
    sg_lib_filtered <- subset(sg_lib, cell %in% intersect(sg_lib$cell, colnames(mtx)))

    #label each cells
    
    lab <- rep("blank", times = ncol(mtx))
    sg_num <- rep(0, times = ncol(mtx))
    names(lab) <- colnames(mtx)
    names(sg_num) <- colnames(mtx)
    sg_in_cell <- data.frame(plyr::count(sg_lib_filtered$cell))

    #find unique label and multiple label
    
    # for (i in 1:nrow(sg_in_cell)) {
    #     x <- sg_in_cell[i, ]
    #     if (x[, 2] == 1) {
    #         lab[x[, 1]] <- sg_lib_filtered[which(sg_lib_filtered$cell == x[, 1]), ]$gene} else {
    #         lab[x[, 1]] <- "multiple"}
    # }
    
    single <- subset(sg_in_cell, freq == 1)$x
    multiple <- subset(sg_in_cell ,freq != 1)$x
    lab[multiple] <- "multiple"
    sg_single <- subset(sg_lib_filtered, cell %in% single)
    rownames(sg_single) <- sg_single$cell
    lab[single] <- lab[single] <- sg_single[intersect(rownames(sg_single), single), "gene"]

    rownames(sg_in_cell) <- sg_in_cell$x
    sg_num[single] <- 1
    sg_num[multiple] <- sg_in_cell[intersect(rownames(sg_in_cell), multiple), "freq"]
    

#     for (i in 1:nrow(sg_in_cell)) {
#         if (sg_in_cell[i, 1] %in% names(sg_num)) {
#             sg_num[sg_in_cell[i, 1]] <- sg_in_cell[i, 2]
#         }
#     }
        
    #add the label information to Seurat object
    
    lab <- as.factor(lab)
    mtx <- AddMetaData(mtx, lab, col.name = "perturbations")
    mtx <- AddMetaData(mtx, sg_num, col.name = "sgRNA_num")

    #Add replicate information to Seurat object
    
    if (!("replicate" %in% colnames(mtx@meta.data))) {
        if (!is.null(replicate)) {
            if (is.vector(replicate) == FALSE) {
                stop("Replicate information must be a vector.")
            } else {
                if (is.null(names(replicate))) {
                    names(replicate) <- colnames(mtx)
                }
            }
            
            if (is.character(replicate) & !is.vector(replicate)) {
                replicate <- read.table(replicate, header = F)
                replicate <- as.vector(replicate[, 1])
            }

        } else {
            replicate <- 1
        }
        mtx <- AddMetaData(mtx, as.factor(replicate), col.name = "replicate")
    }

    
    #calculate percent.mt
    
    if (cal.mt == TRUE) {
        if (species == "Hs") {
            mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^MT-")
        } else {
            mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^mt-")
        }
    } else {
        if (!("percent.mt" %in% colnames(mtx@meta.data))) {
            message("No 'percent.mt' in meta.data and 'cal.mt == FALSE', set all 'percnet.mt' to 0")
            mtx[["percent.mt"]] <- 0
        }
    }
    
    return(mtx)
}
