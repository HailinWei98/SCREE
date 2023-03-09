#' Regulatory Potential Plot
#'
#' Visualize the relationship between potential enhancers and genes, using results from \code{\link[SCREE]{improved_scmageck_lr}}. Modifies from \code{\link[cicero]{plot_connections}}.
#'
#' @param score Data frame or directory of score from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param pval Data frame or directory of p_value from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.
#' @param selected Enhancer regions to visualize. Default is \code{NULL}, all enhancers will be chosen.
#' @param species Only support "Hs" and "Mm". Default is "Hs".
#' @param version Version of the reference genome(Ensembl). Default is "v75".
#' @param gene_annotations Gene annotations stored in data frame format, including c("chromosome", "start", "end", "strand", "transcript") as colnames. Default is \code{NULL}, gene annotations are from \code{\link{ensembldb}}.
#' @param score_cut Score cutoff of \code{improved_scmageck_lr} results. Default is 0.
#' @param pval_cut P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.
#' @param upstream The number of nucleotides upstream of the start site of selected region. Default is 2000000.
#' @param downstream The number of nucleotides downstream of the start site of selected region. Default is 2000000.
#' @param track_size Size of each axis. Default is c(1,.3,.2,.3). If `include_axis_track=FALSE`, track_size should be a vector with 3 elements.
#' @param include_axis_track Logical, should a genomic axis be plotted? Default is \code{TRUE}.
#' @param connection_color Color for connection lines. A single color, the name of a column containing color values, or the name of a column containing a character or factor to base connection colors on. Default is "#7F7CAF".
#' @param connection_color_legend Logical, should connection color legend be shown? Default is \code{TRUE}.
#' @param connection_width Width of connection lines. Default is 2
#' @param connection_ymax Connection y-axis height. Default is NULL, chosen automatically.
#' @param gene_model_color Color for gene annotations. Default is "#81D2C7".
#' @param alpha_by_coaccess Logical, should the transparency of connection lines be scaled based on co-accessibility score? Default is \code{FALSE}.
#' @param gene_model_shap Character scalar. The shape in which to display the track items. Currently only box, arrow, ellipse, and smallArrow are implemented. Default is c("smallArrow", "box"). 
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#' @param html_config Logical, generate and return a list includes the config character string of html. Default is \code{FALSE}.
#' 
#' @import ggplot2
#' @importFrom utils read.table write.table
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @import Gviz
#' @importFrom ensembldb genes
#' @import EnsDb.Hsapiens.v75
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Mmusculus.v75
#' @import EnsDb.Mmusculus.v79
#' @import GenomeInfoDb 
#' @import ggplotify
#' @importFrom IRanges IRanges
#' @export

ciceroPlot <- function(score, pval, selected = NULL, species = "Hs", version = "v75", gene_annotations = NULL, score_cut = 0, pval_cut = 0.05, upstream = 2000000, downstream = 2000000, track_size = c(1, .3, .2, .3), include_axis_track = TRUE, connection_color = "#7F7CAF", connection_color_legend = TRUE, connection_width = 2, connection_ymax = NULL, gene_model_color = "#81D2C7", alpha_by_coaccess = FALSE, gene_model_shape = c("smallArrow", "box"), plot.save = TRUE, prefix = ".", label = "", width = 7, height = 7, png_res = 720, html_config = FALSE){
    
    color_names = NULL
    
    # get genome annotations
    
    if (is.null(gene_annotations)) {
        if (species == "Hs") {
            if (version == "v75") {
                a <- genes(EnsDb.Hsapiens.v75)
            } else if (version == "v79") {
                a <- genes(EnsDb.Hsapiens.v79)
            } else if (version == "v86") {
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        } else if (species == "Mm") {
            if (version == "v75") {
                a <- genes(EnsDb.Mmusculus.v75)
            } else if (version == "v79") {
                a <- genes(EnsDb.Mmusculus.v79)
            }
        } else {
            stop("Species must be one of 'Hs' and 'Mm'")
        }

        gene_anno <- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
        
    } else {
        gene_anno <- gene_annotations
    }
        
    #get score and p-value
    
    if (is.character(score)) {
        score <- read.table(score, header = T, row.names = 1)
    }
    
    if (is.character(pval)) {
        pval <- read.table(pval, header = T, row.names = 1)
    }
    
    if (is.null(selected)) {
        selected <- colnames(score)[which(colnames(score) != "NegCtrl")]
        selected <- selected[grep("^chr", selected)]
    }
    
    if (plot.save == TRUE) {
        
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }
        
        dir <- file.path(dir, "enhancer_function")
        if (!dir.exists(dir)) {
            dir.create(path = dir)
        }

        dir <- file.path(dir, "cicero")
        if (!dir.exists(dir)) {
            dir.create(path = dir)
        }

        img_dir <- file.path(dir, "img")
        if (!dir.exists(img_dir)) {
            dir.create(path = img_dir)
        }
        
        dir <- file.path(dir, "pdf")
        if (!dir.exists(dir)) {
            dir.create(path = dir)
        }
    }

        
    if (is.character(selected)) {
        
        if (length(selected) == 1) {
                                  
            #get score and p-value of selected enhancer
                
            conn_input <- data.frame(score = score[, selected], pval = pval[, selected])
            rownames(conn_input) <- rownames(score)
            conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < pval_cut)
            if (nrow(conn_input) == 0) {
                stop(paste("Cannot find genes passed threshold in ", selected, " scMAGeCK results", sep = "'"))
            }
            #get information from selected enhancer
    
            all <- unlist(strsplit(selected, "[.|:|-|_]"))
            chr <- all[1]
            start <- as.numeric(all[2])
            end <- as.numeric(all[3])
        
            #extend region
    
            minbp <- start - upstream
            maxbp <- end + downstream
            
            #get results
            
            gg <- get_results(chr, start, end, minbp, maxbp, gene_anno, track_size, gene_model_shape,
                              conn_input, connection_color, include_axis_track, score, 
                              score_cut, connection_width, alpha_by_coaccess, color_names)
            if (is.null(gg)) {
                stop("Cannot find gene model close to selected region.")
            }
            
            selected <- paste(chr, start, end, sep = "_")
            
            gg <- gg + 
            labs(y = "Regulatory Potential", 
                 title = selected) +
            theme(plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 20, face = "bold"),
                  axis.title.y = element_text(size = 18, angle = 90, hjust = 0.75))
            
            #save plot
            
            if (plot.save == TRUE) {
                pdf(file = file.path(dir, paste(selected, ".pdf", sep = "")), 
                    width = width, height = height)
                print(gg)
                dev.off()

                png(file.path(img_dir, paste(selected, ".png", sep = "")), 
                    width = width, height = height, unit = "in", res = png_res)
                print(gg)
                dev.off()
            }
            
            if (html_config == TRUE) {
                #generate html config
            
                cicero <- paste(selected, 
                                file.path("results/enhancer_function/cicero/img", 
                                          paste(selected, ".png", sep = "")), 
                                sep = "\" : \"")
                cicero_base64 <- paste(selected, 
                                       knitr::image_uri(file.path(prefix, "results/enhancer_function/cicero/img", 
                                                                  paste(selected, ".png", sep = ""))), 
                                       sep = "\" : \"")
                cicero <- paste("\"cicero\" : {\"", cicero, "\"}", sep = "")
                cicero_base64 <- paste("\"cicero\" : {\"", cicero_base64, "\"}", sep = "")
                return(list(gg, cicero, cicero_base64))  
                
            } else {
                return(gg)  
            }
            
        } else {
                
            #get results for all perturbations
            
            results <- list()
            cicero <- c()
            cicero_base64 <- c()
            j <- 0
            for (perturb in selected) {
                
                #get score and p-value of selected enhancer
                
                conn_input <- data.frame(score = score[, perturb], pval = pval[, perturb])
                rownames(conn_input) <- rownames(score)
                conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < pval_cut)
                if (nrow(conn_input) == 0) {
                    warning(paste("Cannot find genes passed threshold in ", 
                                  perturb, 
                                  " scMAGeCK results", 
                                  sep = "'"))
                    next
                }
                
                #get information from selected enhancer
                
                all <- unlist(strsplit(perturb, "[.|:|-|_]"))
                chr <- all[1]
                start <- as.numeric(all[2])
                end <- as.numeric(all[3])
        
                #extend region
    
                minbp <- start - upstream
                maxbp <- end + downstream
            
                #get results
            
                gg <- get_results(chr, start, end, minbp, maxbp, gene_anno, track_size, gene_model_shape,
                                  conn_input, connection_color, include_axis_track, score, 
                                  score_cut, connection_width, alpha_by_coaccess, color_names)
                if (is.null(gg)) {
                    next
                }
                j <- j + 1
                
                perturb <- paste(chr, start, end, sep = "_")
                
                gg <- gg + 
                labs(y = "Regulatory Potential", 
                     title = perturb) +
                theme(plot.title = element_text(hjust = 0.5), 
                      text = element_text(size = 20, face = "bold"),
                      axis.title.y = element_text(size = 18, angle = 90, hjust = 0.75))
                results[[j]] <- gg
                names(results)[j] <- perturb
                
                #save plot
                
                if (plot.save == TRUE) {
                    pdf(file = file.path(dir, paste(perturb, ".pdf", sep = "")), 
                        width = width, height = height)
                    print(gg)
                    dev.off()

                    png(file.path(img_dir, paste(perturb, ".png", sep = "")), 
                        width = width, height = height, res = png_res)
                    print(gg)
                    dev.off()
                }
                
                #get html config
                
                if (html_config == TRUE) {
                    cicero <- c(cicero, 
                                file.path("results/enhancer_function/cicero/img", 
                                          paste(perturb, ".png", sep = "")))
                    cicero_base64 <- c(cicero, 
                                       knitr::image_uri(file.path(prefix, "results/enhancer_function/cicero/img", 
                                                                  paste(perturb, ".png", sep = ""))))
                    names(cicero)[j] <- perturb
                    names(cicero_base64)[j] <- perturb
                }
                
            }
            
            if (html_config == TRUE) {
                cicero <- paste(names(cicero), cicero, collapse = "\" , \"", sep = "\" : \"")
                cicero <- paste("\"cicero\" : {\"", cicero, "\"}", sep = "")
                cicero_base64 <- paste(names(cicero_base64), cicero_base64, collapse = "\" , \"", sep = "\" : \"")
                cicero_base64 <- paste("\"cicero\" : {\"", cicero_base64, "\"}", sep = "")
                return(list(results, cicero, cicero_base64))
            } else {
                return(results)
            }
        }
    } else {
        stop("Please input correct format of selected perturbations")
    }
}
    
generate_plotting_subset <- function(connections, chr, minbp, maxbp) {
    
    connections$Peak1 <- as.character(connections$Peak1)
    connections$Peak2 <- as.character(connections$Peak2)

    pcolor_map <- data.frame(Peak1 = connections$Peak1,
                             peak_color = connections$peak_color)
    pcolor_map <- pcolor_map[!duplicated(pcolor_map), ]
    connections$peak_color <- NULL

    if (sum(!c("chr_1", "chr_2", "bp1_1", "bp2_1", "bp2_1", "bp2_2") %in% names(connections)) != 0 ) {
        suppressWarnings(connections$chr <- NULL)
        suppressWarnings(connections$bp1 <- NULL)
        suppressWarnings(connections$bp2 <- NULL)
        suppressWarnings(connections$chr_2 <- NULL)
        suppressWarnings(connections$bp1_2 <- NULL)
        suppressWarnings(connections$bp2_2 <- NULL)
        
        connections <- cbind(connections, df_for_coords(connections$Peak1)[ ,c(1, 2, 3)])
        cons2 <- df_for_coords(connections$Peak2) #slow
        cons2$Peak <- NULL
        names(cons2) <- c("chr_2", "bp1_2", "bp2_2")
        connections <- cbind(connections, cons2) #slow
    } else {
        if (!grepl("chr", connections$chr_1[1])) {
            connections$chr_1 <- paste0("chr", connections$chr_1)
        }
        if (!grepl("chr", connections$chr_2[1])) {
            connections$chr_2 <- paste0("chr", connections$chr_2)
        }
        names(connections)[names(connections) == "chr_1"] <- "chr"
        names(connections)[names(connections) == "bp1_1"] <- "bp1"
        names(connections)[names(connections) == "bp2_1"] <- "bp2"
    }

    sub <- connections[connections$chr_2 == chr & connections$bp1 <= maxbp &
                       connections$bp2 <= maxbp & connections$bp1 >= minbp &
                       connections$bp2 >= minbp & connections$bp1_2 <= maxbp &
                       connections$bp2_2 <= maxbp & connections$bp1_2 >= minbp &
                       connections$bp2_2 >= minbp, ]


    sub <- sub[!duplicated(sub), ]

    sub <- merge(sub, pcolor_map, all.x = TRUE)
    sub$peak_color <- as.character(sub$peak_color)
    sub$peak_color[is.na(sub$peak_color)] <- "black"

    return(sub)
}
                            
df_for_coords <- function(coord_strings) {
    coord_strings <- gsub(",", "", coord_strings)
    coord_cols <- as.data.frame(split_peak_names(coord_strings),
                                stringsAsFactors = FALSE )
    names(coord_cols) <- c("chr", "bp1", "bp2")
    coord_cols$Peak <- coord_strings
    coord_cols$bp1 <- as.numeric(coord_cols$bp1)
    coord_cols$bp2 <- as.numeric(coord_cols$bp2)
    coord_cols
}
                        
split_peak_names <- function(inp) {
    out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                    ":|-|_", 3)
    out[ ,1] <- stringi::stri_reverse(out[, 1])
    out[ ,2] <- stringi::stri_reverse(out[, 2])
    out[ ,3] <- stringi::stri_reverse(out[, 3])
    out[ ,c(3, 2, 1), drop = FALSE]
}

make_peak_track <- function(df) {
    df <- df[!duplicated(df[ ,c("chr", "bp1", "bp2", "peak_color")]), ]

    if (sum(duplicated(df[,c("chr", "bp1", "bp2")])) > 0)
        stop(paste("Multiple peak colors correspond to a single peak. Be sure that",
                   "your peak_color column name assigns colors for Peak1 only",
                   collapse = " "))
    
    df2 <- df[!duplicated(df[ ,"Peak1"]), ]
    gr <-  GenomicRanges::GRanges(as.character(df2$chr),
                                  IRanges::IRanges(as.numeric(as.character(df2$bp1)),
                                                   as.numeric(as.character(df2$bp2))), score = df2$coaccess)
    return(gr)
}
    
plotBedpe <- function(bedpedata,
                      chrom,
                      chromstart,
                      chromend,
                      ymax,
                      score_cut,
                      width,
                      alpha_by_coaccess,
                      color_names = NULL)
{ ###### All borrowed and modified from Sushi package.

    if (nrow(bedpedata) == 0) {
        warning("Nothing to plot")
        return()
    }

    bedpedata  <- bedpedata[ ,c("chrom1", "start1", "stop1", "chrom2", "start2",
                                "stop2", "height", "width", "color")]

    # normalize height
    
    maxheight <- ymax

    bedpedata$alpha <- .6
    if(alpha_by_coaccess) {
        bedpedata$alpha <- (bedpedata$height - score_cut)/maxheight
    }
    bedpedata$height <- bedpedata$height/maxheight
    
    # remove any rows with 0 height
    
    bedpedata <- bedpedata[abs(bedpedata$height) > 0, ]

    # reclass data
    
    if (any(class(bedpedata) == "data.table")) {
        for(i in c("start1", "stop1", "start2", "stop2")) {
            bedpedata[[i]] <- as.numeric(as.character((bedpedata[[i]])))
        }
    } else {
        for(i in c("start1", "stop1", "start2", "stop2")) {
            bedpedata[ ,i] <- as.numeric(as.character((bedpedata[ ,i])))
        }
    }

    # add position columns
    
    bedpedata$pos1 = apply(bedpedata[, c("start1", "stop1")], 1, mean)
    bedpedata$pos2 = apply(bedpedata[, c("start2", "stop2")], 1, mean)

    totalrange <- as.numeric(as.character(chromend)) -
    as.numeric(as.character(chromstart))
    if (nrow(bedpedata) == 0) warning("Nothing to plot")

    if (!is.null(color_names)) {
        boxSize <- .3
        spacing <- 0.2
        vspace <- .05
        for (i in seq_len(length(color_names))) {
            grid::grid.lines(unit(c(spacing,spacing + boxSize), "inches"),
                             c(1 - vspace * i, 1 - vspace * i),
                             gp = grid::gpar(col = color_names[i], lwd = width))
            grid::grid.text(x = unit(.1 + (boxSize + spacing), "inches"),
                            y = 1 - vspace * i, just = c(0, 0.5),
                            label = names(color_names)[i])
        }
    }

    # plot the data
    
    grid::grid.function(function(x) list(x = x, y = (score_cut/(ymax))), 
                        gp = grid::gpar(col = "black", lty = "dashed", lwd = width))
    for (row in (seq_len(nrow(bedpedata)))) {
        x1     = bedpedata$pos1[row]
        x2     = bedpedata$pos2[row]
        height = bedpedata$height[row]
        width  = bedpedata$width[row]
        color  = bedpedata$color[row]
        alpha  = bedpedata$alpha[row]
        plotpair(x1, x2, height, totalrange, width, color, chromstart, alpha)
    }
}

plotpair <- function(start, end, height, totalrange,
                     width, color, chromstart, alpha) {
    
    #scale values for plotting
    
    x1 = (min(start, end) - as.numeric(as.character(chromstart)))/totalrange
    x2 = (max(start,end) - as.numeric(as.character(chromstart)))/totalrange
    hx1 <- (x1 + x2)/2
    hy1 <- height/.725

    grid::grid.bezier(x = c(x1, hx1, hx1, x2), y = c(0, hy1, hy1, 0),
                      default.units = "npc",
                      gp = grid::gpar(col = color, lwd = width,
                                    alpha = (alpha * .9 + .1),fontsizecex = 10))
}

get_results <- function(chr, start, end, minbp, maxbp, gene_anno, track_size, gene_model_shape,
                        conn_input, connection_color, include_axis_track, score, 
                        score_cut, connection_width, alpha_by_coaccess, color_names){
    #get gene model

    gene_model <- gene_anno
    gene_model <- gene_model[!is.na(gene_model$chromosome) & 
                             !is.na(gene_model$start) &
                             !is.na(gene_model$end) &
                             !is.na(gene_model$strand) &
                             !is.na(gene_model$transcript), ]
    gene_model <- gene_model[gene_model$chromosome == chr &
                             ((gene_model$start > minbp & gene_model$start < maxbp) |
                              (gene_model$end > minbp & gene_model$end < maxbp) |
                              (gene_model$start < minbp & gene_model$end > maxbp)), ]
    gene_model <- gene_model[gene_model$transcript %in% rownames(score), ]
    

    gene_model <- subset(gene_model, transcript %in% rownames(conn_input))
    grtrack <- Gviz::GeneRegionTrack(gene_model, chromosome = chr, geneSymbols = TRUE,
                                     name = "", fill = "#81D2C7",
                                     col = "#81D2C7",  fontcolor = "black",
                                     fontcolor.group = "black", fontsize.group = 12,
                                     fontsize = 6, shape = gene_model_shape, 
                                     collapseTranscripts = "longest", cex.group = 1)
        
    #generate peak connection data frame
        
    peak <- paste(chr, start, end, sep = "_")
    connection_df <- data.frame()
    if(nrow(gene_model) == 0){
        return(NULL)
    }
    for(j in 1:nrow(gene_model)){
        gene <- gene_model[j, ]
        grange <- paste(chr, gene$start, gene$end, sep = "_")
        coaccess <- conn_input[gene$transcript, "score"]
        connection_df <- rbind(connection_df, c(grange, peak, as.numeric(coaccess)))
    }
    colnames(connection_df) <- c("Peak1", "Peak2", "coaccess")
    connection_df$peak_color <- "black"
    sub <- generate_plotting_subset(connection_df, chr, minbp, maxbp)
    sub$coaccess <- as.numeric(sub$coaccess)
    
    if(nrow(sub) == 0){
        return(NULL)
    }
    
    #get grang of peaks and generate dataTrack
    
    gr <- make_peak_track(sub)
    bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
    dk <- c(colorRampPalette(colors = c("#4DBBD5FF", "white"))(length(bk)/2),
            colorRampPalette(colors = c("white", "#E64B35FF"))(length(bk)/2))
    dtrack <- Gviz::DataTrack(gr, type = c("heatmap"), chromosome = chr, gradient = dk, 
                              ylim = c(-1 , 1), yTicksAt = c(-1, -0.5, 0, 0.5), 
                              name = "scMAGeCK score", fontsize = 18)
    
    #rename sub
    
    if (!nrow(sub) == 0) {
        if (connection_color %in% names(sub)) {
            color_levs <- levels(as.factor(sub[, connection_color]))
            color_names <- rep("temp", length(color_levs))
            names(color_names) <- color_levs
            new_connection_color <- get_colors(sub[, connection_color])
            for (n in color_levs) {
                color_names[n] <- new_connection_color[which(sub[,connection_color] == n)[1]]
            }
            connection_color <- new_connection_color
        
        }
        sub$color <- connection_color

        sub$width <- connection_width

        sub <- sub[ ,c("chr", "bp1", "bp2", "chr_2", "bp1_2", "bp2_2", "coaccess",
                      "width", "color")]
        names(sub) <- c("chrom1", "start1", "stop1", "chrom2", "start2", "stop2",
                        "height", "width", "color")
    } else {
        warning("No connections above score_cutoff")
    
    }
    sub$height <- abs(sub$height)
    ctrack <- CustomTrack(plottingFunction = function(GdObject, prepare) {
        Gviz::displayPars(GdObject) <- list(ylim = c(0, max(abs(as.numeric(connection_df$coaccess)))))
        if (!prepare) {
            plotBedpe(sub, chrom = chr, chromstart = minbp, chromend = maxbp,
                      max(abs(as.numeric(connection_df$coaccess))), score_cut,
                      connection_width, alpha_by_coaccess, color_names)
        }
        return(invisible(GdObject))}, name = "regulatory potential", 
                          fontsize.group = 6, fontsize = 18, cex.title = 1.3)
        
    #in order to show all the gene names
    
    minbp <- minbp - 500000
    
    #get atrack
    
    atrack <- Gviz::GenomeAxisTrack(fontsize = 12, name = "")
        
    #plot
    
    if (include_axis_track == TRUE) {
        gg <- as.ggplot(function()plotTracks(list(ctrack, dtrack, atrack, grtrack), margin = 10,
                                             showTitle = F, from = minbp, to = maxbp, 
                                             chromosome = chr, sizes = track_size, 
                                             transcriptAnnotation = "symbol", background.title = "transparent",
                                             col.border.title = "transparent", lwd.border.title = "transparent",
                                             col.axis = "black", col.title = "black",
                                             fontcolor.legend = "black"))
    } else {
        
        gg <- as.ggplot(function()plotTracks(list(ctrack, dtrack, grtrack), margin = 10,
                                             showTitle = F, from = minbp, to = maxbp, 
                                             chromosome = chr, sizes = track_size, 
                                             transcriptAnnotation = "symbol", background.title = "transparent",
                                             col.border.title = "transparent", lwd.border.title = "transparent",
                                             col.axis = "black", col.title = "black",
                                             fontcolor.legend = "black"))
    }
    return(gg)
}
                        
                        
#' Regulatory Potential Plot
#'
#' Visualize the relationship between potential enhancers and genes, using results from \code{\link[SCREE]{improved_scmageck_lr}}. Modifies from \code{\link[cicero]{plot_connections}}. For scATAC-seq based input, potential enhancers are defined as DA peaks without overlap with promoter regions.
#'                     
#' @param mtx SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.
#' @param score Data frame or directory of score from \code{improved_scmageck_lr}.
#' @param pval Data frame or directory of p_value from \code{improved_scmageck_lr}.
#' @param selected Perturbation to find DA peaks. Default is \code{NULL}, all perturbations will be chosen.
#' @param DApeaks Data frame or directory of DA peaks table generated before. Need a additional column of perturbation information, named "TF". Default is \code{NULL}, will find DA peaks for selected perturbations.
#' @param species Only support "Hs" and "Mm". Default is "Hs".
#' @param version Version of the reference genome(Ensembl). Default is "v75".
#' @param gene_annotations Gene annotations stored in data frame format, including c("chromosome", "start", "end", "strand", "transcript") as colnames. Default is \code{NULL}, gene annotations are from \code{\link{ensembldb}}.
#' @param pro_up The number of nucleotides upstream of the transcription start site that should be included in the promoter region. Default is 3000.
#' @param pro_down The number of nucleotides downstream of the transcription start site that should be included in the promoter region. Default is 0.
#' @param overlap_cut Maximum overlap nucleotides between peaks and promoters. Default is 0.
#' @param NTC The name of negative controls. Default is "NTC".
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the selected perturbation or NTC. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1.
#' @param test.use Denotes which test to use. Available options are: "wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2". See details from \code{\link[Seurat]{FindMarkers}}. Default is "wilcox".
#' @param p_adj_cut Maximum adjust p_value. Default is 0.05.
#' @param logFC_cut Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25. Increasing logfc.threshold speeds up the function, but can miss weaker.
#' @param score_cut Score cutoff of \code{improved_scmageck_lr} results. Default is 0.
#' @param pval_cut P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.
#' @param upstream The number of nucleotides upstream of the start site of selected region. Default is 2000000.
#' @param downstream The number of nucleotides downstream of the start site of selected region. Default is 2000000.
#' @param track_size Size of each axis. Default is c(1,.3,.2,.3). If `include_axis_track=FALSE`, track_size should be a vector with 3 elements.
#' @param include_axis_track Logical, should a genomic axis be plotted? Default is \code{TRUE}.
#' @param connection_color Color for connection lines. A single color, the name of a column containing color values, or the name of a column containing a character or factor to base connection colors on. Default is "#7F7CAF".
#' @param connection_color_legend Logical, should connection color legend be shown? Default is \code{TRUE}.
#' @param connection_width Width of connection lines. Default is 2
#' @param connection_ymax Connection y-axis height. Default is NULL, chosen automatically.
#' @param gene_model_color Color for gene annotations. Default is "#81D2C7".
#' @param alpha_by_coaccess Logical, should the transparency of connection lines be scaled based on co-accessibility score? Default is \code{FALSE}.
#' @param gene_model_shap Character scalar. The shape in which to display the track items. Currently only box, arrow, ellipse, and smallArrow are implemented. Default is c("smallArrow", "box"). 
#' @param plot.save Logical, save plots or not. Default is \code{TRUE}. 
#' @param table.save Logical, save DA peaks table or not. Default is \code{TRUE}.                         
#' @param prefix Path to save the plots. Default is current directory.
#' @param label The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".
#' @param width Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param height Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.
#' @param png_res The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.
#' @param html_config Logical, generate and return a list includes the config character string of html. Default is \code{FALSE}.
#' 
#' @import ggplot2
#' @import Gviz
#' @importFrom utils read.table write.table
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom ensembldb genes
#' @import GenomeInfoDb
#' @importFrom IRanges IRanges
#' @import ggplotify
#' @export 
    
ATACciceroPlot<- function(mtx, score, pval, selected =  NULL, DA_peaks = NULL, species = "Hs", version = "v75", gene_annotations = NULL, pro_up = 2000, pro_down = 0, overlap_cut = 0, NTC = "NTC", min.pct = 0.1, test.use = "wilcox", p_adj_cut = 0.05, logFC_cut = 0.25, score_cut = 0, pval_cut = 0.05, upstream = 2000000, downstream = 2000000, track_size = c(1, .3, .2, .3), include_axis_track = TRUE, connection_color = "#7F7CAF", connection_color_legend = TRUE, connection_width = 2, connection_ymax = NULL, gene_model_color = "#81D2C7", alpha_by_coaccess = FALSE, gene_model_shape = c("smallArrow", "box"), plot.save = TRUE, table.save = TRUE, prefix = ".", label = "", width = 7, height = 7, png_res = 720, html_config = FALSE) {
    
    color_names = NULL
    
    # get genome annotations
    
    if (is.null(gene_annotations)) {
        if (species == "Hs") {
            if (version == "v75") {
                a <- genes(EnsDb.Hsapiens.v75)
            } else if (version == "v79") {
                a <- genes(EnsDb.Hsapiens.v79)
            } else if (version == "v86") {
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        } else if (species == "Mm") {
            if (version == "v75") {
                a <- genes(EnsDb.Mmusculus.v75)
            } else if (version == "v79") {
                a <- genes(EnsDb.Mmusculus.v79)
            }
        } else {
            stop("Species must be one of 'Hs' and 'Mm'")
        }

        gene_anno <- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
        
    } else {
        gene_anno <- gene_annotations
    }
        
    #get score and p-value
    
    if (is.character(score)) {
        score <- read.table(score, header = T, row.names = 1)
    }
    
    if (is.character(pval)) {
        pval <- read.table(pval, header = T, row.names = 1)
    }
    
    #create store list
    
    results <- list()
        
    #get selected TFs list
    
    if (is.null(selected)) {
        selected <- colnames(score)[which(colnames(score) != "NegCtrl")]
    }
    
    if (plot.save == TRUE | table.save == TRUE) {
        dir <- file.path(prefix, "results")
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }
        
        dir <- file.path(dir, "enhancer_function")
        if (!dir.exists(dir)) {
            dir.create(path = dir)
        }

        dir <- file.path(dir, "cicero")
        if (!dir.exists(dir)) {
            dir.create(path = dir)
        }

        img_dir <- file.path(dir, "img")
        if (!dir.exists(img_dir)) {
            dir.create(path = img_dir)
        }
        
        pdf_dir <- file.path(dir, "pdf")
        if (!dir.exists(pdf_dir)) {
            dir.create(path = pdf_dir)
        }
    }
    
    if (is.character(selected)) {
        
        if (length(selected) == 1) {
            
            #get DA peaks and enhancer list
            
            if (is.null(DA_peaks)) {
                da_peak <- DApeaks(mtx, selected, NTC, min.pct, test.use, p_adj_cut, logFC_cut)
            } else if (is.character(DA_peaks)) {
                DApeaks <- read.table(DA_peaks, header = T, row.names = 1)
                if (!(selected %in% unique(DApeaks$TF))) {
                    stop("Please select correct name of perturbations.")
                }
                da_peak <- subset(DApeaks, TF == selected)
            } else {
                da_peak <- subset(DA_peaks, TF == selected)
            }
            
            if(is.null(da_peak)){
                stop(paste("Cannot find DA peaks pass threshold in ", 
                           selected, 
                           " FindMarkers results", 
                           sep = "'"))
            } else {
                
                if (table.save == TRUE & is.null(DA_peaks)){
                    write.table(da_peak, file = file.path(dir, "DApeaks.txt"), quote = FALSE)
                }
                
                enhancer_list <- enhancer(da_peak, gene_anno, overlap_cut, pro_up = pro_up, pro_down = pro_down)
                
                #generate config character of html output
                
                if (html_config == TRUE) {
                    new_da <- data.frame(Perturbations = rep(gene, nrow(da_peak)),
                                         DApeaks = paste(da_peak$chromosome, 
                                                         da_peak$start, 
                                                         da_peak$end, 
                                                         sep = "_"),
                                         log2FC = da_peak$avg_log2FC,
                                         p_val_adj = da_peak$p_val_adj)
                    new_da <- paste("<tr><td>", new_da$Perturbations, 
                                    "</td><td>", new_da$DApeaks, 
                                    "</td><td>", round(new_da$log2FC, 3), 
                                    "</td><td>", new_da$p_val_adj, 
                                    "</td></tr>", sep = "")
                }

                all_enhancer <- paste(enhancer_list$chromosome, enhancer_list$start, enhancer_list$end, sep = "_")
                
            }
        
            if (nrow(da_peak) == nrow(enhancer_list)) {
                warning("All DA peaks has overlap with promoters, using all DA peaks passed threshold as input list")
            }

                
            #get score and p-value of selected enhancer
                
            conn_input <- data.frame(score = score[, selected], pval = pval[, selected])
            rownames(conn_input) <- rownames(score)
            conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < pval_cut)
            if (nrow(conn_input) == 0) {
                stop(paste("Cannot find genes passed threshold in ", selected, " scMAGeCK results", sep = "'"))
            }
            
            #get results for all candidate enhancer
                
            j <- 0
            
            if (plot.save == TRUE) {
                dir1 <- file.path(pdf_dir, selected)
                if (!dir.exists(dir1)) {
                    dir.create(path = dir1)
                }

                dir2 <- file.path(img_dir, selected)
                if (!dir.exists(dir2)) {
                    dir.create(path = dir2)
                }
            }

            if (html_config == TRUE) {
                
                cicero <- c()
                cicero_base64 <- c()
                cicero_prefix <- file.path("results/enhancer_function/cicero/img", selected)
            
            }

            for (i in 1:nrow(enhancer_list)) {
                chr <- enhancer_list[i, "chromosome"]
                start <- as.numeric(enhancer_list[i, "start"])
                end <- as.numeric(enhancer_list[i, "end"])
                    
                #extend region
                    
                minbp <- start - upstream
                maxbp <- end + downstream
                    
                #get results
                    
                gg <- get_results(chr, start, end, minbp, maxbp, gene_anno, track_size, gene_model_shape,
                                  conn_input, connection_color, include_axis_track, score, 
                                  score_cut, connection_width, alpha_by_coaccess, color_names)
                if (is.null(gg)) {
                    next
                }
                j <- j + 1
                
                peak <- paste(chr, start, end, sep = "_")
                
                #generate enhancer directory of html

                gg <- gg + 
                labs(y = "Regulatory Potential", 
                     title = paste(selected, peak, sep = ":")) +
                theme(plot.title = element_text(hjust = 0.5), 
                      text = element_text(size = 20, face = "bold"),
                      axis.title.y = element_text(size = 18, angle = 90, hjust = 0.75))
                results[[j]] <- gg
                names(results)[j] <- peak
                
                #save plot
                
                if (plot.save == TRUE) {
                    pdf(file = file.path(dir1, paste(peak, ".pdf", sep = "")), 
                        width = width, height = height)
                    print(gg)
                    dev.off()

                    png(file.path(dir2, paste(peak, ".png", sep = "")), 
                        width = width, height = height, unit = "in", res = png_res)
                    print(gg)
                    dev.off()
                }
                
                if (html_config == TRUE) {
                    
                    cicero <- c(cicero, paste(file.path(cicero_prefix, peak), ".png", sep = ""))
                    cicero_base64 <- c(cicero_base64, paste(file.path(prefix, cicero_prefix, peak), ".png", sep = ""))
                    names(cicero)[j] <- peak
                    names(cicero_base64)[j] <- peak
                 
                }
                
            }
            
            #generate html config
            
            if (html_config == TRUE) {
                cicero <- paste(names(cicero), cicero, collapse = "\" , \"", sep = "\" : \"")
                cicero <- paste("\"", selected, "\" : {\"", cicero, "\"}", sep = "")
                cicero <- paste("\"cicero\" : {\"", cicero, "\"}", sep = "")
                cicero_base64 <- paste(names(cicero_base64), cicero_base64, collapse = "\" , \"", sep = "\" : \"")
                cicero_base64 <- paste("\"", selected, "\" : {\"", cicero_base64, "\"}", sep = "")
                cicero_base64 <- paste("\"cicero\" : {\"", cicero_base64, "\"}", sep = "")
                DA <- paste("", new_da, collapse = "", sep = "")
                DA <- paste("\"DApeaks\" : {\"", DA, "\"}", sep = "")
            }
            
        } else {
                
            #get results for all perturbations
                
            k <- 0
            results <- list()
            DA <- c()
            da_index <- 0
            html_results <- data.frame()
            all_cicero <- c()
            all_cicero_base64 <- c()
            if (table.save == TRUE) {
                all_peak <- data.frame()
            }
            
            for(TF in selected){
                
                #get DA peaks and enhancer list
                if (is.null(DA_peaks)) {
                    da_peak <- DApeaks(mtx, selected = TF, NTC, min.pct, test.use, p_adj_cut, logFC_cut)
                } else if (is.character(DA_peaks)) {
                    DApeaks <- read.table(DA_peaks, header = T, row.names = 1)
                    if (!(TF %in% unique(DApeaks$TF))) {
                        warning(paste(TF, " is not in the column of 'TF'."))
                        next
                    }
                    da_peak <- subset(DApeaks, TF == TF)
                } else {
                    da_peak <- subset(DA_peaks, TF == TF)
                }

                if(is.null(da_peak)){
                    warning(paste("Cannot find DA peaks pass threshold in ", 
                                  TF, " FindMarkers results", sep = "'"))
                    next
                } else {
                    enhancer_list <- enhancer(da_peak, 
                                              gene_anno, 
                                              overlap_cut, 
                                              pro_up = pro_up, 
                                              pro_down = pro_down)
                    if (table.save == TRUE & is.null(DA_peaks)) {
                        da_peak$TF <- TF
                        all_peak <- rbind(all_peak, da_peak)
                        colnames(all_peak) <- colnames(da_peak)
                    }
                    
                    #generate enhancer and DApeaks information of html
                    
                    if (html_config == TRUE) {
                        new_da <- data.frame(Perturbations = rep(TF, nrow(da_peak)), 
                                             DApeaks = paste(da_peak$chromosome, 
                                                             da_peak$start, 
                                                             da_peak$end, 
                                                             sep = "_"), 
                                             log2FC = da_peak$avg_log2FC,
                                             p_val_adj = da_peak$p_val_adj)
                        html_results <- rbind(html_results, new_da)
                        names(html_results) <- c("Perturbations", "DApeaks", "log2FC", "p_val_adj")
                        new_da <- paste("<tr><td>", new_da$Perturbations, 
                                        "</td><td>", new_da$DApeaks, 
                                        "</td><td>", round(new_da$log2FC, 3), 
                                        "</td><td>", new_da$p_val_adj, 
                                        "</td></tr>", sep = "")
                        new_da <- paste("", new_da, collapse = "", sep = "")
                        DA <- c(DA, new_da)
                        da_index <- da_index + 1
                        names(DA)[da_index] <- TF
                    }
    
                }
                
                if (nrow(da_peak) == nrow(enhancer_list)) {
                    warning(paste("All DA peaks has overlap with promoters in ", 
                                  TF, ", using all DA peaks passed threshold as input list", sep = "'"))
                }

                #get score and p-value of selected enhancer
                
                conn_input <- data.frame(score = score[ ,TF], pval = pval[ ,TF])
                rownames(conn_input) <- rownames(score)
                conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < pval_cut)
                if (nrow(conn_input) == 0) {
                    warning(paste("Cannot find genes passed threshold in ", TF, " scMAGeCK results", sep = "'"))
                        next
                }
                
                #get results for all candidate enhancer
                    
                sub_results <- list()
                
                j <- 0
                
                cicero <- c()
                cicero_base64 <- c()
                
                if (plot.save == TRUE) {

                    dir1 <- file.path(pdf_dir, TF)
                    if (!dir.exists(dir1)) {
                        dir.create(path = dir1)
                    }

                    dir2 <- file.path(img_dir, TF)
                    if (!dir.exists(dir2)) {
                        dir.create(path = dir2)
                    }

                }
                
                for (i in 1:nrow(enhancer_list)) {
                    chr <- enhancer_list[i, "chromosome"]
                    start <- as.numeric(enhancer_list[i, "start"])
                    end <- as.numeric(enhancer_list[i, "end"])
                    
                    #extend region
                    
                    minbp <- start - upstream
                    maxbp <- end + downstream
                    
                    #get results
                    
                    gg <- get_results(chr, start, end, minbp, maxbp, gene_anno, track_size, gene_model_shape,
                                      conn_input, connection_color, include_axis_track, score, 
                                      score_cut, connection_width, alpha_by_coaccess, color_names)
                    if (is.null(gg)) {
                        next
                    }
                    j <- j + 1
                    
                    peak <- paste(chr, start, end, sep = "_")
                                      
                    gg <- gg + 
                    labs(y = "Regulatory Potential", 
                         title = paste(TF, peak, sep = ":")) +
                    theme(plot.title = element_text(hjust = 0.5), 
                          text = element_text(size = 20, face = "bold"),
                          axis.title.y = element_text(size = 18, angle = 90, hjust = 0.75))
                    sub_results[[j]] <- gg
                    names(sub_results)[j] <- peak
                    
                                    
                    #save plot
                    
                    if (plot.save == TRUE) {
                        pdf(file = file.path(dir1, paste(peak, ".pdf", sep = "")), 
                            width = width, height = height)
                        print(gg)
                        dev.off()

                        png(file.path(dir2, paste(peak, ".png", sep = "")), 
                            width = width, height = height, unit = "in", res = png_res)
                        print(gg)
                        dev.off()
                    }
                    
                    #generate enhancer directory of html
                    
                    if (html_config == TRUE) {
                        cicero_prefix <- file.path("results/enhancer_function/cicero/img", TF)
                        cicero <- c(cicero, paste(file.path(cicero_prefix, peak), ".png", sep = ""))
                        cicero_base64 <- c(cicero_base64, 
                                           knitr::image_uri(paste(file.path(prefix, cicero_prefix, peak), ".png", sep = "")))
                        names(cicero)[j] <- peak
                        names(cicero_base64)[j] <- peak
                    }

                }
                
                if (html_config == TRUE) {
                    if (length(cicero) != 0) {
                        cicero <- paste(names(cicero), cicero, collapse = "\" , \"", sep = "\" : \"")
                        cicero <- paste("\"", TF, "\" : {\"", cicero, "\"}", sep = "")
                        cicero_base64 <- paste(names(cicero_base64), cicero_base64, collapse = "\" , \"", sep = "\" : \"")
                        cicero_base64 <- paste("\"", TF, "\" : {\"", cicero_base64, "\"}", sep = "")
                        all_cicero <- c(all_cicero, cicero)
                        all_cicero_base64 <- c(all_cicero_base64, cicero_base64)
                    }
                }

                if(length(sub_results) != 0){
                    k <- k + 1
                    results[[k]] <- sub_results
                    names(results)[k] <- TF
                }
            }
            
            #save DA peaks table
            
            if (table.save == TRUE & is.null(DA_peaks)) {
                write.table(all_peak, file = file.path(dir, "DApeaks.txt"), quote = FALSE)
            }
            
            #generate results of html
            
            if (html_config == TRUE) {
                if (nrow(html_results) != 0) {
                    all_da <- paste("<tr><td>", html_results$Perturbations, 
                                    "</td><td>", html_results$DApeaks, 
                                    "</td><td>", round(html_results$log2FC, 3), 
                                    "</td><td>", html_results$p_val_adj, 
                                    "</td></tr>", sep = "")
                    all_da <- paste(all_da, collapse = "", sep = "")
                    DA <- c(all_da, DA)
                    names(DA)[1] <- "All_Results"
                } else {
                    DA <- paste("\"All_Results\" : ", "\"\"", sep = "")
                }
                
                DA <- paste(names(DA), DA, collapse = "\" , \"", sep = "\" : \"")
                DA <- paste("\"DApeaks\" : {\"", DA, "\"}", sep = "")
                all_cicero <- paste("", all_cicero, collapse = ", ", sep = "")
                all_cicero <- paste("\"cicero\" : {", all_cicero, "}", sep = "")
                all_cicero_base64 <- paste("", all_cicero_base64, collapse = ", ", sep = "")
                all_cicero_base64 <- paste("\"cicero\" : {", all_cicero_base64, "}", sep = "")
            }
        }
        if (html_config == FALSE) {
            return(results)
        } else {
            return(list(results, DA, all_cicero, all_cicero_base64))
        }
    } else {
        stop("Please input correct format of selected perturbations")
    }
}