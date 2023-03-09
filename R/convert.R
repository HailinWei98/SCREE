#' Assign sgRNAs to Cells
#'
#' Assign sgRNAs to cells based on sgRNA count matrix or data frame includes sgRNA counts in each cell.
#'
#' @param sg_lib Data frame or count matrix or directory to a txt file containing information of sgRNA count in each cell. 
#' @param type Type of input data frame, can be one of "CountMatrix" or "DataFrame". "CountMatrix" means a sgRNA count matrix similar to gene count matrix, whose rows are sgRNAs and columns are cells. "DataFrame" means a data frame includes three columns: "freq", "cell", "barcode". "freq" is sgRNA count in each cell, "barcode" is sgRNA name. Default is \code{"CountMatrix"}.
#' @param row_names Logical, is rownames included in the txt file or not. Default is \code{FALSE}.
#' @param freq_cut Cutoff of sgRNA count numbers. For each cell, only sgRNA with counts more than freq_cut will be retained. Default is \code{20}.
#' @param freq_percent Cutoff of sgRNA frequency. For each cell, only sgRNA with frequency more than freq_percent will be retained. Default is \code{0.8}.
#' @param freq Name of the "freq" column, only used for "DataFrame" type. Default is \code{"freq"}.
#' @param cell Name of the "cell" column, only used for "DataFrame" type. Default is \code{"cell"}.
#' @param barcode Name of the "barcode" column, only used for "DataFrame" type. Default is \code{"barcode"}.
#' @param unique Logical, only retain cells with the top target that passed the \code{freq_percent} or not. Default is \code{FALSE}.
#' 
#' @import reshape2
#' @importFrom utils read.table write.table
#' @export

sgRNAassign <- function (sg_lib, type = "CountMatrix", row_names = FALSE, freq_cut = 20, freq_percent = 0.8, freq = "freq", cell = "cell", barcode = "barcode", unique = FALSE) {
    
    #get sg_lib
    
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        if (row_names == FALSE) {
            sg_lib <- read.table(sg_lib, header = T)
        } else {
            sg_lib <- read.table(sg_lib, header = T, row.names = 1)
        }
    }
    
    #assign sgRNA to each cell
    
    if (type == "CountMatrix") {
        
        #filter results
        
        sg_lib <- as.data.frame(sg_lib)
        sg_lib <- sg_lib * (sg_lib >= freq_cut)
        sg_lib2 <- as.data.frame(t(t(sg_lib) / (colSums(sg_lib) + (colSums(sg_lib) == 0))))
        
        #convert matrix
        
        a <- reshape2::melt(sg_lib, variable.name = "cell", value.name = "count")
        a2 <- reshape2::melt(sg_lib2, variable.name = "cell", value.name = "count")
        
        #get sgRNA information
        
        a$sgRNA <- rep(rownames(sg_lib), ncol(sg_lib))
        a2$sgRNA <- rep(rownames(sg_lib2), ncol(sg_lib2))
        
        #filter freq_percent
        
        b <- subset(a2, count > 0)
        b <- subset(b, count >= freq_percent)
        c <- subset(a2, !(cell %in% b$cell))
        c <- subset(c, count > 0)
        if (unique == TRUE) {
            sg_lib <- b[, c("cell", "sgRNA")]
        } else {
            sg_lib <- rbind(b[, c("cell", "sgRNA")], c[, c("cell", "sgRNA")])
        }
        

        sg_lib$gene <- gsub("_.*\\b", "", sg_lib$sgRNA)
        colnames(sg_lib) <- c("cell", "barcode", "gene")
        
    } else if (type == "DataFrame") {
        if (sum(colnames(sg_lib) %in% c(freq, cell, barcode)) != 3) {
            stop("freq, barcode, or cell column names not found in input sgRNA library.")
        }
        new_lib <- data.frame(cell = sg_lib[, cell], barcode = sg_lib[, barcode], freq = sg_lib[, freq])
        a <- subset(new_lib, freq >= freq_cut)
        
        b <- a$cell[duplicated(a$cell)]

        c <- subset(a, cell %in% b)

        d <- subset(a, !(cell %in% b))

        e <- data.frame()
        
        for (i in unique(c$cell)) {
            test <- subset(c, cell == i)
            test$freq <- test$freq / sum(test$freq)
            if (max(test$freq) >= freq_percent) {
                test <- subset(test, freq >= freq_percent)
            }
            e <- rbind(test, e)
            
        }
        if (unique == TRUE) {
            e <- subset(e, freq >= freq_percent)
        }
        sg_lib <- rbind(e[, c("cell", "barcode")], d[, c("cell", "barcode")])
        sg_lib$gene <- gsub("_.*\\b", "", sg_lib$barcode)
        colnames(sg_lib) <- c("cell", "barcode", "gene")
    } else {
        stop("Please select from c('CountMatrix', 'DataFrame')")
    }
    
    return(sg_lib)  
}