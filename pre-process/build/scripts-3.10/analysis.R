library(SCREE)

args <- commandArgs(T)
type <- args[1]
replicate <- args[2]
mtx <- args[3]
sg_lib <- args[4]
sg_format <- args[5]
freq_cut <- as.numeric(args[6])
freq_percent <- as.numeric(args[7])
unique <- args[8]
project <- args[9]
species <- args[10]
fragments <- args[11]
prefix <- args[12]
label <- args[13]
feature_frac <- as.numeric(args[14])
nFeature <- as.numeric(args[15])
nCount <- as.numeric(args[16])
percent.mt <- as.numeric(args[17])
blank <- args[18]
NTC <- args[19]
sg.split <- args[20]
raster <- args[21]
ref_version <- args[22]
gene_type <- args[23]
pval_cut <- as.numeric(args[24])
score_cut <- as.numeric(args[25])
article <- args[26]
article_name <- args[27]
data <- args[28]
data_name <- args[29]
padj_cut <- as.numeric(args[30])
logFC_cut <- as.numeric(args[31])

if (is.character(replicate)) {
    if (replicate == "1") {
        replicate <- 1
    } else {
        replicate <- read.table(replicate, header = F)
    }
}

if (is.character(mtx)) {
    if (stringr::str_ends(mtx, ".rds")) {
        mtx <- readRDS(mtx)
    } else {
        mtx <- Read10X(mtx)
    }
}

if (unique == "True") {
    unique <- TRUE
} else {
    unique <- FALSE
}

if (blank == "True") {
    blank <- TRUE
} else {
    blank <- FALSE
}

if (raster == "True") {
    raster <- TRUE
} else {
    raster <- FALSE
}

if (sg_format == "matrix") {
    sg_lib <- read.table(sg_lib, header = T, row.names = 1)
    sg_lib <- sgRNAassign(sg_lib = sg_lib, freq_cut = freq_cut, freq_percent = freq_percent, unique = unique)
} else if (sg_format == "10X") {
    sg_lib <- mtx$`CRISPR Guide Capture`
    mtx <- mtx$`Gene Expression`
    sg_lib <- sgRNAassign(sg_lib = sg_lib, freq_cut = freq_cut, freq_percent = freq_percent, unique = unique)
} else if (sg_format == "dataframe") {
    sg_lib <- read.table(sg_lib, header = T)
    sg_lib <- sgRNAassign(sg_lib = sg_lib, type = "DataFrame", freq_cut = freq_cut, freq_percent = freq_percent, unique = unique)
} else if (sg_format == "table") {
    sg_lib <- read.table(sg_lib, header = T)
}

if (type == "RNA") {
    
    mtx <- CreateSeuratObject(counts = mtx, project = project)
    mtx <- Add_meta_data(sg_lib = sg_lib, mtx = mtx, species = species, replicate = replicate)
    saveRDS(mtx, file = file.path(prefix, paste(label, "RNA.rds")))
    sgRNA_quality_plot(sg_lib = sg_lib, mtx = mtx, bar_width = bar_width, prefix = prefix, label = label)
    mtx_QC <- scQC(mtx = mtx, species = species, prefix = prefix, label = label, gene_frac = feature_frac, nFeature = c(nFeature, 500000), nCount = nCount, mt = percent.mt, blank_NTC = blank)
    mtx_QC <- normalize_scale(mtx = mtx_QC)
    saveRDS(mtx_QC, file = file.path(prefix, paste(label, "RNA_QC.rds")))
    mixscape <- IntegratedMixscape(sg_lib = sg_lib, mtx = mtx_QC, NTC = NTC, sg.split = sg.split, prefix = prefix, label = label)
    saveRDS(mixscape, file = file.path(prefix, paste(label, "mixscape.rds")))
    
    results <- improved_scmageck_lr(BARCODE = sg_lib, RDS = mtx_QC, NEGCTRL = NTC, SELECT_GENE = NULL, LABEL = "improved", PERMUTATION = 10000, SAVEPATH = prefix, LAMBDA = 0.01)
    score <- results[[1]][, -1]
    pval <- results[[2]][, -1]
    
    DE_gene_plot(score = score, pval = pval, project = project, prefix = prefix, label = label, pval_cut = pval_cut, score_cut = score_cut, sort_by = "number", y_break = c(50, 150), width = 8, height = 6)
    
    volcano(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, height = 6, width = 6, showCategory = 5)
    
    if (ncol(score) <= 50) {
        heatmap(score = score, pval = pval, prefix = prefix, label = label)
    }
    
    if (species == "Hs") {
        
        GOenrichment(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, DE_gene_to_use = "all", database = "org.Hs.eg.db", gene_type = gene_type, showCategory = 10)

    } else {
        
        GOenrichment(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, DE_gene_to_use = "all", database = "org.Mm.eg.db", gene_type = gene_type, showCategory = 10)

    }
    
    config <- config_generation(mtx = mtx, mtx_QC = mtx_QC, sg_lib = sg_lib, score = score, pval = pval, project = project, prefix = prefix, label = label, species = species, version = ref_version, type = "RNA", NTC = NTC, article = article, data = data, article_name = article_name, data_name = data_name, gene_type = gene_type, score_cut = score_cut, pval_cut = pval_cut, DA = NULL, cicero = NULL, enhancer = NULL, base64 = FALSE)
    
    html_output(html_dir = "RNA.html", config, prefix = prefix, label = label)
    
} else if (type == "enhancer") {
    
    mtx <- CreateSeuratObject(counts = mtx, project = project)
    mtx <- Add_meta_data(sg_lib = sg_lib, mtx = mtx, species = species, replicate = replicate)
    saveRDS(mtx, file = file.path(prefix, paste(label, "RNA.rds")))
    sgRNA_quality_plot(sg_lib = sg_lib, mtx = mtx, bar_width = bar_width, prefix = prefix, label = label)
    mtx_QC <- scQC(mtx = mtx, species = species, prefix = prefix, label = label, gene_frac = feature_frac, nFeature = c(nFeature, 500000), nCount = nCount, mt = percent.mt, blank_NTC = blank)
    mtx_QC <- normalize_scale(mtx = mtx_QC)
    saveRDS(mtx_QC, file = file.path(prefix, paste(label, "RNA_QC.rds")))
    if (length(unique(sg_lib$barcode)) > 100) {
        mixscape <- IntegratedMixscape(sg_lib = sg_lib, mtx = mtx_QC, NTC = NTC, sg.split = sg.split, prefix = prefix, label = label, perturb.sig.only = TRUE)
    } else {
        mixscape <- IntegratedMixscape(sg_lib = sg_lib, mtx = mtx_QC, NTC = NTC, sg.split = sg.split, prefix = prefix, label = label)
    }
    
    saveRDS(mixscape, file = file.path(prefix, paste(label, "mixscape.rds")))
    
    results <- improved_scmageck_lr(BARCODE = sg_lib, RDS = mtx_QC, NEGCTRL = NTC, SELECT_GENE = NULL, LABEL = "improved", PERMUTATION = 10000, SAVEPATH = prefix, LAMBDA = 0.01)
    score <- results[[1]][, -1]
    pval <- results[[2]][, -1]
    
    DE_gene_plot(score = score, pval = pval, project = project, prefix = prefix, label = label, pval_cut = pval_cut, score_cut = score_cut, sort_by = "number", y_break = c(50, 150), width = 8, height = 6)
    
    volcano(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, height = 6, width = 6, showCategory = 5)
    
    if (ncol(score) <= 50) {
        heatmap(score = score, pval = pval, prefix = prefix, label = label)
    }
    
    if (species == "Hs") {
        
        GOenrichment(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, DE_gene_to_use = "all", database = "org.Hs.eg.db", gene_type = gene_type, showCategory = 10)

    } else {
        
        GOenrichment(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, DE_gene_to_use = "all", database = "org.Mm.eg.db", gene_type = gene_type, showCategory = 10)

    }
    

    cicero <- ciceroPlot(score = score, pval = pval, selected = NULL, species = species, version = ref_version, gene_annotations = NULL, pval_cut = pval_cut, score_cut = score_cut, upstream = 2000000, downstream = 2000000, track_size = c(1,.3,.2,.3), include_axis_track = TRUE, prefix = prefix, label = label, html_config = TRUE)
    
    config <- config_generation(mtx = mtx, mtx_QC = mtx_QC, sg_lib = sg_lib, score = score, pval = pval, project = project, prefix = prefix, label = label, species = species, version = ref_version, type = "enhancer", NTC = NTC, article = article, data = data, article_name = article_name, data_name = data_name, gene_type = gene_type, score_cut = score_cut, pval_cut = pval_cut, DA = NULL, cicero = cicero[[2]], enhancer = NULL, base64 = FALSE)
    
    html_output(html_dir = "enhancer.html", config, prefix = prefix, label = label)

} else if (type == "ATAC") {
    
    peak <- CreateChromatinAssay(counts = mtx, fragments = fragments, sep = c("_", "_"))
    peak <- CreateSeuratObject(counts = peak, project = project, assay = "peak")
    peak <- ATAC_Add_meta_data(sg_lib = sg_lib, mtx = peak, species = species, replicate = replicate, fragments = fragments)
    saveRDS(peak, file = file.path(prefix, paste(label, "peak.rds")))
    fragmentsSize(mtx = peak, fragments = fragments, CBCindex = 4, startIndex = 2, endIndex = 3, maxSize = 1000, prefix = prefix, label = label)
    sgRNA_quality_plot(sg_lib = sg_lib, mtx = peak, bar_width = bar_width, prefix = prefix, label = label)
    peak_QC <- ATAC_scQC(mtx = peak, chromatin.assay = TRUE, species = species, prefix = prefix, label = label, peak_frac = feature_frac, nFeature = c(nFeature, 500000), nCount = nCount, FRiP = percent.mt / 100, blank_NTC = blank)
    saveRDS(peak, file = file.path(prefix, paste(label, "peak_QC.rds")))
    mtx <- CalculateGeneActivity(mtx = peak, fragments = fragments, species = species, version = , gene_type = gene_type, protein_coding = TRUE, pro_up = 2000, pro_down = 0, sep = c("_", "_"))
    if (species == "Hs") {
        mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^MT-")
    } else {
        mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^mt-")
    }
    saveRDS(mtx, file = file.path(prefix, paste(label, "RNA.rds")))
    mtx_QC <- scQC(mtx = mtx, species = species, prefix = prefix, label = label, gene_frac = feature_frac, nFeature = c(nFeature, 500000), nCount = nCount, mt = percent.mt, blank_NTC = blank)
    mtx_QC <- normalize_scale(mtx = mtx_QC)
    saveRDS(mtx_QC, file = file.path(prefix, paste(label, "RNA_QC.rds")))
    mixscape <- IntegratedMixscape(sg_lib = sg_lib, mtx = mtx_QC, NTC = NTC, sg.split = sg.split, prefix = prefix, label = label)
    saveRDS(mixscape, file = file.path(prefix, paste(label, "mixscape.rds")))
    
    results <- improved_scmageck_lr(BARCODE = sg_lib, RDS = mtx_QC, NEGCTRL = NTC, SELECT_GENE = NULL, LABEL = "improved", PERMUTATION = 10000, SAVEPATH = prefix, LAMBDA = 0.01)
    score <- results[[1]][, -1]
    pval <- results[[2]][, -1]
    
    DE_gene_plot(score = score, pval = pval, project = project, prefix = prefix, label = label, pval_cut = pval_cut, score_cut = score_cut, sort_by = "number", y_break = c(50, 150), width = 8, height = 6)
    
    volcano(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, height = 6, width = 6, showCategory = 5)
    
    if (ncol(score) <= 50) {
        heatmap(score = score, pval = pval, prefix = prefix, label = label)
    }
    
    if (species == "Hs") {
        
        GOenrichment(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, DE_gene_to_use = "all", database = "org.Hs.eg.db", gene_type = gene_type, showCategory = 10)

    } else {
        
        GOenrichment(score = score, pval = pval, selected = NULL, prefix = prefix, label = label, score_cut = score_cut, pval_cut = pval_cut, DE_gene_to_use = "all", database = "org.Mm.eg.db", gene_type = gene_type, showCategory = 10)

    }
    
    cicero <- ATACciceroPlot(mtx = peak_QC, score = score, pval = pval, selected = NULL, species = species, version = ref_version, gene_annotations = NULL, pro_up = 2000, pro_down = 0, overlap_cut = 0, pval_cut = pval_cut, score_cut = score_cut, p_adj_cut = padj_cut, logFC_cut = logFC_cut, NTC = NTC, min.pct = 0.1, upstream = 2000000, downstream = 2000000, test.use = "wilcox", track_size = c(1,.3,.2,.3), include_axis_track = TRUE, prefix = prefix, label = label, html_config = TRUE)
    
    config <- config_generation(mtx = mtx, mtx_QC = mtx_QC, sg_lib = sg_lib, score = score, pval = pval, project = project, prefix = prefix, label = label, species = species, version = ref_version, type = "ATAC", NTC = NTC, article = article, data = data, article_name = article_name, data_name = data_name, gene_type = gene_type, score_cut = score_cut, pval_cut = pval_cut, DA = cicero[[2]], cicero = cicero[[3]], enhancer = NULL, base64 = FALSE)
    
    html_output(html_dir = "ATAC.html", config, prefix = prefix, label = label)
    
}