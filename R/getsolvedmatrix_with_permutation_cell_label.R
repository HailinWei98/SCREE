#' function definitions ##### Return list of 2 matrix:
#' liner regression results of perturbation score and
#' p-value calculated by permutation

getsolvedmatrix_with_permutation_cell_label <- function(Xm, Ym, lambda = 0.01, npermutation = 1000) {
    Amat_ret = getsolvedmatrix(Xm, Ym, lambda = lambda)
    Amat_ret_higher = Amat_ret * 0
    
    # permute N times randomly shuffle cell labels
    
    for(npm in 1:npermutation){
        if (npm%%100 == 0) {
            message(paste("Permutation:", npm, "/", npermutation, "..."))
    }
    cells_shu = sample(colnames(Ym), ncol(Ym))
    Xm_s = Xm[cells_shu, ]
    Amat_random = getsolvedmatrix(Xm_s, Ym, lambda = lambda)

    Amat_ret_higher = Amat_ret_higher + (abs(Amat_random) > abs(Amat_ret)) * 1
    }
    Amat_ret_higher = Amat_ret_higher/npermutation
    return(list(Amat_ret, Amat_ret_higher))
}
