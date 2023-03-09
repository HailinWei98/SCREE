#' function definitions ##### Return a list with 2 matrix:
#' index matrix, but values of NegCtrl columns are all 1;
#' genes * cells matrix, similar to scale_data matrix,
#' but all values over quantile with
#' setting of outlier_threshold will be replace with the outlier_threshold quantile

single_gene_matrix_regression <- function(targetobj, ngctrlgene = c("NonTargetingControlGuideForHuman"),
                                          indmatrix = NULL, selected_genes_list = NULL, NTC_baseline = TRUE) {
    
    # return X matrix and Y matrix for regression note that all the ngctrlgene are merged into one
    # column, 'NegCtrl' if indmatrix is provided, the Xmat will be constructed from indmatrix
    
    outlier_threshold = 0.95

    #because select cells and genes before input,in this step, there is no need to select
    
    YmatT = targetobj
    select_genes <- rownames(YmatT)
    if (!is.null(selected_genes_list)) {
        select_genes = select_genes[select_genes %in% selected_genes_list]
        if (length(select_genes) == 0) {
            stop("No genes left for regression. Check your selected gene list.")
        }
    }
    message(paste("Selected genes:", length(select_genes)))

    if (NTC_baseline == TRUE) {
        if (length(ngctrlgene) == 1) {
            colnames(indmatrix)[colnames(indmatrix) == ngctrlgene]<- "NegCtrl"
            indmatrix[, "NegCtrl"] <- 1
        } else if (length(ngctrlgene) > 1) {
            indmatrix <- indmatrix[, -which(colnames(indmatrix) %in% ngctrlgene[-1])] #only remain one columns of NTC
            colnames(indmatrix)[colnames(indmatrix) == ngctrlgene[1]] <- "NegCtrl"
            indmatrix[, "NegCtrl"] <- 1
        }
    }

    YmatT = YmatT[select_genes, ]
    
    # remove outliers
        
    Ymat_outlier = apply(YmatT, 1, function(X) {
        return(quantile(X, probs = outlier_threshold))
    })
    YmatT <- YmatT - Ymat_outlier
    YmatT[YmatT > 0] <- 0
    YmatT <- YmatT + Ymat_outlier

    return(list(indmatrix, YmatT))
}

