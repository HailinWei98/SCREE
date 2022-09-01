#' @export

DApeaks <- function(object, selected, NTC = "NTC", min.pct = 0.2, 
                    test.use = "wilcox", p_adj_cut = 0.05, logFC_cut = 1){
    
    if(is.character(object)){
        peak <- readRDS(object)
    }else{
        peak <- object
    }
    
    #rename ident
    
    if("perturbations" %in% colnames(peak@meta.data)){
        if(selected %in% peak$perturbations){
            peak@active.ident <- peak$perturbations
        }else{
            stop(paste("Cannot find ", selected, " in perturbations", sep = "'"))
        }
    }else{
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
    if(nrow(da_peaks) == 0){
        return(NULL)
    }
    da_peak <- subset(da_peaks, p_val_adj <= p_adj_cut)
    if(nrow(da_peak) == 0){
        return(NULL)
    }
    
    da_peak$chromosome <- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,1]
    da_peak$start <- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,2]
    da_peak$end <- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,3]
    
    return(da_peak)
}

#' @export

enhancer <- function(da_peak, pro_anno, overlap_cut, pro_up = 3000, pro_down = 0){
    
    #get enhancer list
        
    enhancer_list <- data.frame()
    for(i in 1:nrow(da_peak)){
        chr <- da_peak[i, "chromosome"]
        start <- as.numeric(da_peak[i, "start"])
        end <- as.numeric(da_peak[i, "end"])
        minbp <- start - 2 * (pro_up + pro_down)
        maxbp <- end + 2 * (pro_up + pro_down)
        peak_model <- pro_anno
        peak_model <- peak_model[!is.na(peak_model$chromosome) & 
                                 !is.na(peak_model$start) &
                                 !is.na(peak_model$end) &
                                 !is.na(peak_model$strand), ]
        peak_model <- peak_model[peak_model$chromosome == chr &
                                 ((peak_model$start > minbp & peak_model$start < maxbp) |
                                  (peak_model$end > minbp & peak_model$end < maxbp) |
                                  (peak_model$start < minbp & peak_model$end > maxbp)), ]
        if(nrow(peak_model) == 0){
            enhancer_list <- rbind(enhancer_list, c(chr, start, end))
        }else{
            pro_seq <- apply(peak_model[,c("start", "end")], 1, function(x){seq(x[1], x[2])})
            overlap <- length(intersect(seq(start, end), pro_seq))
            if(overlap <= overlap_cut){
                enhancer_list <- rbind(enhancer_list, c(chr, start, end))
            }
        }
    }
    if(nrow(enhancer_list != 0)){
        colnames(enhancer_list) <- c("chromosome", "start", "end")
        rownames(enhancer_list) <- paste(enhancer_list$chromosome, enhancer_list$start, enhancer_list$end, sep = "-")
    }else{
        enhancer_list <- da_peak[, c("chromosome", "start", "end")]
    }
    
    return(enhancer_list)
}