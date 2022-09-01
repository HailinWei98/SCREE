#' function definitions ##### Improved scMAGeCK_LR, save and return two data.frame:
#' rownames indicates selected_gene while colnames indicates perturbations.
#' @export

improved_scmageck_lr <- function(BARCODE, RDS, NEGCTRL, SELECT_GENE = NULL, LABEL = NULL,
                        PERMUTATION = NULL, SAVEPATH = "./", LAMBDA = 0.01, NTC_baseline = TRUE) {
    if (!is.null(LABEL)) {
        data_label = LABEL
    } else {
        data_label = "sample1"
    }

    if (!is.null(PERMUTATION)) {
        n_permutation = as.integer(PERMUTATION)
    } else {
        n_permutation = 10000
    }

    # read cell assignment and libray file ####
    
    if (is.character(BARCODE)) {
        bc_dox = read.table(BARCODE, header = TRUE, as.is = TRUE)
    } else {
        bc_dox = BARCODE
    }
    
    

    if (sum(colnames(bc_dox) %in% c("cell", "barcode", "gene")) != 3) {
        stop("cell, barcode, or gene column names not found in barcode file.")
    }

    message(paste("Total barcode records:", nrow(bc_dox)))

    # load neg control guides ####
    
    ngctrlgenelist = strsplit(NEGCTRL, ",")[[1]] #NTC need to be a string split by ","
    message(paste("Neg Ctrl guide:", paste(ngctrlgenelist, collapse = ";")))

    # read Seurat RDS file ####
    
    if (is.character(RDS)) {
        message(paste("Reading RDS file:", RDS))
        targetobj = readRDS(RDS)
    } else {
        targetobj = RDS
    }
    
    if (isS4(targetobj)) {
        targetobj <- SeuratObject::GetAssayData(object = targetobj, slot = "scale.data")
    }
    
    # check if names are consistent
    
    nmatch = sum(bc_dox[, "cell"] %in% colnames(x = targetobj))
    if (nmatch == 0) {
        stop("No cell names match in expression matrix and barcode file.")
    }
    bc_dox <- subset(bc_dox, cell %in% colnames(x = targetobj))

    # convert to ind_matrix ####
    
    ind_matrix <- frame2indmatrix(bc_dox, targetobj) #return TRUE and FALSE matrix
    message(paste("Index matrix dimension:", nrow(ind_matrix), ",", ncol(ind_matrix)))

    # try to perform matrix regresson on single genes ####
    
    mat_for_single_reg = single_gene_matrix_regression(targetobj, selected_genes_list = SELECT_GENE,
                                                       ngctrlgene = ngctrlgenelist, indmatrix = ind_matrix, 
                                                       NTC_baseline = NTC_baseline)
    Xmat = mat_for_single_reg[[1]]

    # Xmat[,which(colnames(Xmat)%in%ngctrlgenelist)[1]]=1 # already integrated into function
    
    Ymat = mat_for_single_reg[[2]]
    
    Ymat = Ymat[, rownames(Xmat)]
    rm(list = c("mat_for_single_reg", "targetobj", "ind_matrix"))

    # remove values in Y mat
    
    Amat_pm_lst = getsolvedmatrix_with_permutation_cell_label(Xmat, Ymat, lambda = LAMBDA, 
                                                              npermutation = n_permutation)
    Amat = Amat_pm_lst[[1]]
    Amat_pval = Amat_pm_lst[[2]]
   
    if (!is.null(SAVEPATH)) {
        write.table(data.frame(Perturbedgene = rownames(Amat), Amat), 
                    file = file.path(SAVEPATH, paste(data_label, "_score.txt", sep = "")), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
        write.table(data.frame(Perturbedgene = rownames(Amat), Amat_pval), 
                    file = file.path(SAVEPATH, paste(data_label, "_score_pval.txt", sep = "")), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
    return(list(data.frame(Perturbedgene = rownames(Amat), Amat), 
                data.frame(Perturbedgene = rownames(Amat), Amat_pval)))
}
