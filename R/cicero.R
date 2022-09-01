#' @export
#' @import Gviz
#' @import ensembldb
#' @import ggplotify

ciceroPlot<- function(score_dir, pval_dir, selected = NULL, species = "Hs", version = "v75", gene_annotations = NULL, 
                      p_val_cut = 0.05, score_cut = 0, prefix = "./", label = "", 
                      upstream = 2000000, downstream = 2000000,
                      track_size = c(1, .3, .2, .3), include_axis_track = TRUE, connection_color = "#7F7CAF",
                      connection_color_legend = TRUE, connection_width = 2, connection_ymax = NULL, 
                      gene_model_color = "#81D2C7", alpha_by_coaccess = FALSE, 
                      gene_model_shape = c("smallArrow", "box"), html_config = FALSE){
    
    color_names = NULL
    
    # get genome annotations
    
    if(is.null(gene_annotations)){
        if(species == "Hs"){
            if(version == "v75"){
                a <- genes(EnsDb.Hsapiens.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Hsapiens.v79)
            }else if(version == "v86"){
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        }else if(species == "Mm"){
            if(version == "v75"){
                a <- genes(EnsDb.Mmusculus.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Mmusculus.v79)
            }
        }
        gene_anno<- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
    }else{
        gene_anno <- gene_annotations
            }
        
    #get score and p-value
    
    if(is.character(score_dir)){
        score <- read.table(score_dir, header = T, row.names = 1)
    }else{
        score <- score_dir
    }
    if(is.character(pval_dir)){
        pval <- read.table(pval_dir, header = T, row.names = 1)
    }else{
        pval <- pval_dir
    }
    
    if(is.null(selected)){
        selected <- colnames(score)[which(colnames(score) != "NegCtrl")]
        selected <- selected[grep("^chr", selected)]
    }
        
    dir <- file.path(prefix, "pdf")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "enhancer_function")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "cicero")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }

    img_dir <- file.path(prefix, "img")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
        
    img_dir <- file.path(img_dir, "enhancer_function")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
        
    img_dir <- file.path(img_dir, "cicero")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
        
    if(is.character(selected)){
        if(length(selected) == 1){
                                  
            #get score and p-value of selected enhancer
                
            conn_input <- data.frame(score = score[, selected], pval = pval[, selected])
            rownames(conn_input) <- rownames(score)
            conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < p_val_cut)
            if(nrow(conn_input) == 0){
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
            if(is.null(gg)){
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
            
            pdf(file = file.path(dir, paste(selected, ".pdf", sep = "")))
            print(gg)
            dev.off()
        
            png(file.path(img_dir, paste(selected, ".png", sep = "")), 
                width = 600 * 3, height = 600 * 3, res = 72 * 3)
            print(gg)
            dev.off()
            
            #generate html config
            
            cicero <- paste(selected, file.path("img/enhancer_function/cicero", paste(selected, ".png", sep = "")), 
                            sep = "\" : \"")
            cicero <- paste("\"cicero\" : {\"", cicero, "\"}", sep = "")
            
            if (html_config == TRUE) {
              return(list(gg, cicero))  
            } else {
              return(gg)  
            }
            
        } else {
                
            #get results for all perturbations
            
            results <- list()
            cicero <- c()
            j <- 0
            for(perturb in selected){
                
                #get score and p-value of selected enhancer
                
                conn_input <- data.frame(score = score[, perturb], pval = pval[, perturb])
                rownames(conn_input) <- rownames(score)
                conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < p_val_cut)
                if(nrow(conn_input) == 0){
                    warning(paste("Cannot find genes passed threshold in ", perturb, " scMAGeCK results", sep = "'"))
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
                if(is.null(gg)){
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
            
                pdf(file = file.path(dir, paste(perturb, ".pdf", sep = "")))
                print(gg)
                dev.off()
        
                png(file.path(img_dir, paste(perturb, ".png", sep = "")), 
                    width = 600 * 3, height = 600 * 3, res = 72 * 3)
                print(gg)
                dev.off()
                
                #get html config
                cicero <- c(cicero, file.path("img/enhancer_function/cicero", paste(perturb, ".png", sep = "")))
                names(cicero)[j] <- perturb

                
            }
            
            cicero <- paste(names(cicero), cicero, collapse = "\" , \"", sep = "\" : \"")
            cicero <- paste("\"cicero\" : {\"", cicero, "\"}", sep = "")
            if (html_config == TRUE) {
                return(list(results, cicero))
            } else {
                return(results)
            }
        }
    }else{
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

    if(sum(!c("chr_1", "chr_2", "bp1_1", "bp2_1", "bp2_1", "bp2_2") %in% names(connections)) != 0 ) {
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
        if(!grepl("chr", connections$chr_1[1])) {
            connections$chr_1 <- paste0("chr", connections$chr_1)
        }
        if(!grepl("chr", connections$chr_2[1])) {
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
        
    peak <- paste(chr, start, end, sep = "-")
    connection_df <- data.frame()
    if(nrow(gene_model) == 0){
        return(NULL)
    }
    for(j in 1:nrow(gene_model)){
        gene <- gene_model[j, ]
        grange <- paste(chr, gene$start, gene$end, sep = "-")
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
    dtrack <- DataTrack(gr, type = c("heatmap"), chromosome = chr, gradient = dk, 
                        ylim = c(-1 , 1), yTicksAt = c(-1, -0.5, 0, 0.5), name = "scMAGeCK score", ,fontsize = 18)
    
    #rename sub
    
    if (!nrow(sub) == 0) {
        if (connection_color %in% names(sub)) {
            color_levs <- levels(as.factor(sub[, connection_color]))
            color_names <- rep("temp", length(color_levs))
            names(color_names) <- color_levs
            new_connection_color <- get_colors(sub[, connection_color])
            for(n in color_levs) {
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
        }else {
        warning("No connections above score_cutoff")
        }
    sub$height <- abs(sub$height)
    ctrack <- CustomTrack(plottingFunction = function(GdObject, prepare) {
        Gviz::displayPars(GdObject) <- list(ylim = c(0, max(abs(as.numeric(connection_df$coaccess)))))
        if(!prepare) {
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
    
    if(include_axis_track == TRUE){
        gg <- as.ggplot(function()plotTracks(list(ctrack, dtrack, atrack, grtrack), margin = 10,
                                             showTitle = F, from = minbp, to = maxbp, 
                                             chromosome = chr, sizes = track_size, 
                                             transcriptAnnotation = "symbol", background.title = "transparent",
                                             col.border.title = "transparent", lwd.border.title = "transparent",
                                             col.axis = "black", col.title = "black",
                                             fontcolor.legend = "black"))
    }else{
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
                        
#' @export    
    
ATACciceroPlot<- function(object, score_dir, pval_dir, selected =  NULL, species = "Hs", version = "v75",
                          gene_annotations = NULL, pro_annotations = NULL, pro_up = 3000, pro_down = 0, 
                          overlap_cut = 0, p_val_cut = 0.05, score_cut = 0, p_adj_cut = 0.05, logFC_cut = 1, 
                          NTC = "NTC", min.pct = 0.2, upstream = 2000000, downstream = 2000000, test.use = "wilcox", 
                          track_size = c(1, .3, .2, .3), include_axis_track = TRUE, connection_color = "#7F7CAF", 
                          connection_color_legend = TRUE, connection_width = 2, connection_ymax = NULL, 
                          gene_model_color = "#81D2C7", alpha_by_coaccess = FALSE, 
                          gene_model_shape = c("smallArrow", "box"), prefix = "./", label = "", html_config = FALSE){
    
    color_names = NULL
    
    # get genome annotations
    
    if(is.null(gene_annotations)){
        if(species == "Hs"){
            if(version == "v75"){
                a <- genes(EnsDb.Hsapiens.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Hsapiens.v79)
            }else if(version == "v86"){
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        }else if(species == "Mm"){
            if(version == "v75"){
                a <- genes(EnsDb.Mmusculus.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Mmusculus.v79)
            }
        }
        gene_anno <- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
    }else{
        gene_anno <- gene_annotations
            }
    
    #get promoters annotations
        
    if(is.null(pro_annotations)){
        if(species == "Hs"){
            if(version == "v75"){
                b <- promoters(EnsDb.Hsapiens.v75)
            }else if(version == "v79"){
                b <- promoters(EnsDb.Hsapiens.v79)
            }else if(version == "v86"){
                b <- promoters(EnsDb.Hsapiens.v86)
            }  
        } else if (species == "Mm"){
            if(version == "v75"){
                b <- promoters(EnsDb.Mmusculus.v75)
            }else if(version == "v79"){
                b <- promoters(EnsDb.Mmusculus.v79)
            }
        }
        pro_anno <- data.frame(b)
        pro_anno$chromosome <- paste0("chr", pro_anno$seqnames)
    }else{
        pro_anno <- pro_annotations
            }
        
    #get score and p-value
    
    if(is.character(score_dir)){
        score <- read.table(score_dir, header = T, row.names = 1)
    }else{
        score <- score_dir
    }
    if(is.character(pval_dir)){
        pval <- read.table(pval_dir, header = T, row.names = 1)
    }else{
        pval <- pval_dir
    }    
    
    #create store list
    
    results <- list()
        
    #get selected TFs list
    
    if(is.null(selected)){
        selected <- colnames(score)[which(colnames(score) != "NegCtrl")]
    }
    
    dir <- file.path(prefix, "pdf")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "enhancer_function")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "cicero")
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }

    img_dir <- file.path(prefix, "img")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
        
    img_dir <- file.path(img_dir, "enhancer_function")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
        
    img_dir <- file.path(img_dir, "cicero")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
    
    if(is.character(selected)){
        
        if(length(selected) == 1){
            #get DA peaks and enhancer list
                
            da_peak <- DApeaks(object, selected, NTC, min.pct, test.use, p_adj_cut, logFC_cut)
            if(is.null(da_peak)){
                stop(paste("Cannot find DA peaks pass threshold in ", selected, " FindMarkers results", sep = "'"))
            }else{
                enhancer_list <- enhancer(da_peak, pro_anno, overlap_cut, pro_up = pro_up, pro_down = pro_down)
                
                #generate config character of html output
                
                new_da <- data.frame(Perturbations = rep(gene, nrow(da_peak)),
                                     DApeaks = paste(da_peak$chromosome, da_peak$start, da_peak$end, sep = "_"),
                                     log2FC = da_peak$avg_log2FC,
                                     p_val_adj = da_peak$p_val_adj)
                new_da <- paste("<tr><td>", new_da$Perturbations, 
                                "</td><td>", new_da$DApeaks, 
                                "</td><td>", new_da$log2FC, 
                                "</td><td>", new_da$p_val_adj, 
                                "</td></tr>", sep = "")
                all_enhancer <- paste(enhancer_list$chromosome, enhancer_list$start, enhancer_list$end, sep = "_")
                
            }
            if(nrow(da_peak) == nrow(enhancer_list)){
                warning("All DA peaks has overlap with promoters, using all DA peaks passed threshold as input list")
            }

                
            #get score and p-value of selected enhancer
                
            conn_input <- data.frame(score = score[, selected], pval = pval[, selected])
            rownames(conn_input) <- rownames(score)
            conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < p_val_cut)
            if(nrow(conn_input) == 0){
                stop(paste("Cannot find genes passed threshold in ", selected, " scMAGeCK results", sep = "'"))
            }
            
            #get results for all candidate enhancer
                
            j <- 0
            
            dir1 <- file.path(dir, selected)
            if (!dir.exists(dir1)) {
                dir.create(path = dir1)
            }
            
            dir2 <- file.path(img_dir, selected)
            if (!dir.exists(dir2)) {
                dir.create(path = dir2)
            }
            
            cicero <- c()
            cicero_prefix <- file.path("img/enhancer_function/cicero", selected)
            
            for(i in 1:nrow(enhancer_list)){
                chr <- enhancer_list[i, "chromosome"]
                start <- as.numeric(enhancer_list[i, "start"])
                end <- as.numeric(enhancer_list[i, "end"])
                    
                #extend region
                    
                minbp<- start - upstream
                maxbp<- end + downstream
                    
                #get results
                    
                gg <- get_results(chr, start, end, minbp, maxbp, gene_anno, track_size, gene_model_shape,
                                  conn_input, connection_color, include_axis_track, score, 
                                  score_cut, connection_width, alpha_by_coaccess, color_names)
                if(is.null(gg)){
                    next
                }
                j <- j + 1
                
                peak <- paste(chr, start, end, sep = "_")
                
                #generate enhancer directory of html
                
                cicero <- c(cicero, paste(file.path(cicero_prefix, peak), ".png", sep = ""))
                names(cicero)[j] <- peak
                
                gg <- gg + 
                labs(y = "Regulatory Potential", 
                     title = paste(selected, peak, sep = ":")) +
                theme(plot.title = element_text(hjust = 0.5), 
                      text = element_text(size = 20, face = "bold"),
                      axis.title.y = element_text(size = 18, angle = 90, hjust = 0.75))
                results[[j]] <- gg
                names(results)[j] <- peak
                
                #save plot
                
                pdf(file = file.path(dir1, paste(peak, ".pdf", sep = "")))
                print(gg)
                dev.off()
        
                png(file.path(dir2, paste(peak, ".png", sep = "")), 
                    width = 600 * 3, height = 600 * 3, res = 72 * 3)
                print(gg)
                dev.off()
            }
            
            #generate html config
            
            cicero <- paste(names(cicero), cicero, collapse = "\" , \"", sep = "\" : \"")
            cicero <- paste("\"", selected, "\" : {\"", cicero, "\"}", sep = "")
            cicero <- paste("\"cicero\" : {\"", cicero, "\"}", sep = "")
            DA <- paste("", new_da, collapse = "", sep = "")
            DA <- paste("\"DApeaks\" : {\"", DA, "\"}", sep = "")
            
        } else {
                
            #get results for all perturbations
                
            k <- 0
            results <- list()
            DA <- c()
            da_index <- 0
            html_results <- data.frame()
            all_cicero <- c()
            
            for(TF in selected){
                
                #get DA peaks and enhancer list
                    
                da_peak <- DApeaks(object, selected = TF, NTC, min.pct, test.use, p_adj_cut, logFC_cut)
                if(is.null(da_peak)){
                    warning(paste("Cannot find DA peaks pass threshold in ", 
                                  TF, " FindMarkers results", sep = "'"))
                    next
                } else {
                    enhancer_list <- enhancer(da_peak, pro_anno, overlap_cut, pro_up = pro_up, pro_down = pro_down)
                    
                    #generate enhancer and DApeaks information of html
                    
                    new_da <- data.frame(Perturbations = rep(TF, nrow(da_peak)), 
                                         DApeaks = paste(da_peak$chromosome, da_peak$start, da_peak$end, sep = "."), 
                                         log2FC = da_peak$avg_log2FC,
                                         p_val_adj = da_peak$p_val_adj)
                    html_results <- rbind(html_results, new_da)
                    names(html_results) <- c("Perturbations", "DApeaks", "log2FC", "p_val_adj")
                    new_da <- paste("<tr><td>", new_da$Perturbations, 
                                    "</td><td>", new_da$DApeaks, 
                                    "</td><td>", new_da$log2FC, 
                                    "</td><td>", new_da$p_val_adj, 
                                    "</td></tr>", sep = "")
                    new_da <- paste("", new_da, collapse = "", sep = "")
                    DA <- c(DA, new_da)
                    da_index <- da_index + 1
                    names(DA)[da_index] <- TF
                    
                }
                if(nrow(da_peak) == nrow(enhancer_list)){
                    warning(paste("All DA peaks has overlap with promoters in ", 
                                  TF, ", using all DA peaks passed threshold as input list", sep = "'"))
                }

                #get score and p-value of selected enhancer
                
                conn_input <- data.frame(score = score[ ,TF], pval = pval[ ,TF])
                rownames(conn_input) <- rownames(score)
                conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < p_val_cut)
                if(nrow(conn_input) == 0){
                    warning(paste("Cannot find genes passed threshold in ", TF, " scMAGeCK results", sep = "'"))
                        next
                }
                
                #get results for all candidate enhancer
                    
                sub_results <- list()
                
                j <- 0
                
                cicero <- c()
                
                for(i in 1:nrow(enhancer_list)){
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
                    if(is.null(gg)){
                        next
                    }
                    j <- j + 1
                    
                    peak <- paste(chr, start, end, sep = "_")
                    
                    #generate enhancer directory of html
                    
                    dir1 <- file.path(dir, TF)
                    if (!dir.exists(dir1)) {
                        dir.create(path = dir1)
                    }

                    dir2 <- file.path(img_dir, TF)
                    if (!dir.exists(dir2)) {
                        dir.create(path = dir2)
                    }
                
                    cicero_prefix <- file.path("img/enhancer_function/cicero", TF)
                    cicero <- c(cicero, paste(file.path(cicero_prefix, peak), ".png", sep = ""))
                    names(cicero)[j] <- peak
                    
                    gg <- gg + 
                    labs(y = "Regulatory Potential", 
                         title = paste(TF, peak, sep = ":")) +
                    theme(plot.title = element_text(hjust = 0.5), 
                          text = element_text(size = 20, face = "bold"),
                          axis.title.y = element_text(size = 18, angle = 90, hjust = 0.75))
                    sub_results[[j]] <- gg
                    names(sub_results)[j] <- peak
                    
                                    
                    #save plot
                
                    pdf(file = file.path(dir1, paste(peak, ".pdf", sep = "")))
                    print(gg)
                    dev.off()
        
                    png(file.path(dir2, paste(peak, ".png", sep = "")), 
                        width = 600 * 3, height = 600 * 3, res = 72 * 3)
                    print(gg)
                    dev.off()
                    
                }
                if (length(cicero) != 0) {
                    cicero <- paste(names(cicero), cicero, collapse = "\" , \"", sep = "\" : \"")
                    cicero <- paste("\"", TF, "\" : {\"", cicero, "\"}", sep = "")
                    all_cicero <- c(all_cicero, cicero)
                }

                
                if(length(sub_results) != 0){
                    k <- k + 1
                    results[[k]] <- sub_results
                    names(results)[k] <- TF
                }
            }
            
            #generate results of html
            
            if (nrow(html_results) != 0) {
                all_da <- paste("<tr><td>", html_results$Perturbations, 
                                "</td><td>", html_results$DApeaks, 
                                "</td><td>", html_results$log2FC, 
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
        }
        if (html_config == FALSE) {
            return(results)
        } else {
            return(list(results, DA, all_cicero))
        }
    }else{
        stop("Please input correct format of selected perturbations")
    }
}