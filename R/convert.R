#' @export

sgRNAassign <- function(sg_lib, type = "CountMatrix", row_names = FALSE, freq_cut = 20, freq_percent = 0.8, 
                        freq = "freq", cell = "cell", barcode = "barcode"){
    
    #get sg_lib
    
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        if(row_names == FALSE){
            sg_lib <- read.table(sg_lib, header = T)
        }else{
            sg_lib <- read.table(sg_lib, header = T, row.names = 1)
        }
    }
    
    #assign sgRNA to each cell
    
    if(type == "CountMatrix"){
        
        #filter results
        
        sg_lib <- as.data.frame(sg_lib)
        sg_lib <- sg_lib * (sg_lib > freq_cut)
        sg_lib2 <- as.data.frame(t(t(sg_lib)/(colSums(sg_lib) + (colSums(sg_lib) == 0))))
        
        #convert matrix
        
        a <- reshape2::melt(sg_lib, variable.name = "cell", value.name = "count")
        a2 <- reshape2::melt(sg_lib2, variable.name = "cell", value.name = "count")
        
        #get sgRNA information
        
        a$sgRNA <- rep(rownames(sg_lib), ncol(sg_lib))
        a2$sgRNA <- rep(rownames(sg_lib2), ncol(sg_lib2))
        
        #filter freq_percent
        
        b <- subset(a2, count > freq_percent)
        c <- subset(a2, !(cell %in% b$cell))
        c <- subset(c, count > 0)
        sg_lib <- rbind(b[, c("cell", "sgRNA")], c[, c("cell", "sgRNA")])

        sg_lib$gene <- gsub("_.*\\b", "", sg_lib$sgRNA)
        colnames(sg_lib) <- c("cell", "barcode", "gene")
        
    } else if (type == "DataFrame"){
        if (sum(colnames(sg_lib) %in% c(freq, cell, barcode)) != 3) {
            stop("freq, barcode, or cell column names not found in input sgRNA library.")
        }
        new_lib <- data.frame(cell = sg_lib[, cell], barcode = sg_lib[, barcode], freq = sg_lib[, freq])
        a <- subset(new_lib, freq > freq_cut)
        
        b <- a$cellBC[duplicated(a$cell)]

        c <- subset(a, cell %in% b)

        d <- subset(a, !(cell %in% b))

        e <- data.frame()
        
        for(i in unique(c$cell)){
            test <- subset(c, cell == i)
            test$freq <- test$freq/sum(test$freq)
            if(max(test$freq) >= freq_percent){
                test <- subset(test, freq >= freq_percent)
            }
            e <- rbind(test[, c("cell", "barcode")], e)
        }
        sg_lib <- rbind(e, d[, c("cell", "barcode")])
        sg_lib$gene <- gsub("_.*\\b", "", sg_lib$barcode)
        colnames(sg_lib) <- c("cell", "barcode", "gene")
    }else{
        stop("Please select from c('CountMatrix', 'DataFrame')")
    }
    return(sg_lib)  
}