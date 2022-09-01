#' function definitions ##### Create a 0-1 index matrix to store the infromation of
#' perturbations.0 means cells without the perturbation while 1 means cells with
#' the perturbation

frame2indmatrix <- function(bc_d, targetobj) {
    #rnm = unique(bc_d$cell)
    rnm = colnames(targetobj)
    cnm = unique(bc_d$gene)
    scalef = targetobj
    message(paste(length(unique(bc_d$cell)), "..."))
    message(paste(ncol(scalef), "..."))
    #rnm = rnm[!is.na(rnm)] #remove NA
    #rnm = rnm[rnm %in% colnames(scalef)]
    if (length(rnm) == 0) {
        stop("Cell names do not match in expression matrix and barcode.")
    }
    cnm = cnm[!is.na(cnm)]#remove NA
    ind_matrix = matrix(rep(0, length(rnm) * length(cnm)), nrow = length(rnm))
    rownames(ind_matrix) = rnm
    colnames(ind_matrix) = cnm
    row <- bc_d[, 'cell']
    col <- bc_d[, 'gene']
    test <- (row %in% rnm) & (col %in% cnm)
    idx <- as.matrix(data.frame(row[test], col[test]))
    #idx <- cbind(row[test], col[test])
    ind_matrix[idx] <- 1
    return(ind_matrix)
}
